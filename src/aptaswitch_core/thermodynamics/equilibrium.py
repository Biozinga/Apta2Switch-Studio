"""Ideal dilute-solution equilibrium with concentrations expressed only in mol/L.

This module consumes formation constants; it does not predict molecular energies.
Every complex column of ``stoichiometry`` gives the number of each strand type.
One monomer column per strand is required, with formation free energy / log K zero.
Distinct complexes with the same composition remain distinct chemical species.

Mass action is parameterized by logarithms of free monomer concentrations. Solving
the logarithmic conservation equations avoids exponentiating large association
constants. No clipping, solvent-density guess, or additional symmetry factor is
applied to the equilibrium concentrations.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
from numbers import Real

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import least_squares
from scipy.special import logsumexp


GAS_CONSTANT_KCAL_MOL_K = 0.00198720425864083


class EquilibriumConvergenceError(RuntimeError):
    """The requested conservation accuracy could not be established."""


@dataclass(frozen=True)
class EquilibriumResult:
    """Concentrations and independently checked conservation diagnostics.

    Log concentrations retain very dilute species that underflow float64 in M.
    Impossible complexes (including any absent strand) have concentration zero
    and log concentration -inf. Arrays are read-only. ``iterations`` counts
    optimizer function evaluations across all continuation stages.
    """

    concentrations_molar: NDArray[np.float64]
    log_concentrations_molar: NDArray[np.float64]
    mass_balance_residual_molar: NDArray[np.float64]
    max_relative_mass_balance_error: float
    max_log_mass_action_error: float
    iterations: int


def formation_energies_to_log_constants(
    stoichiometry: ArrayLike,
    formation_free_energies_kcal_mol: ArrayLike,
    *,
    temperature_kelvin: float = 310.15,
    standard_concentration_molar: float = 1.0,
) -> NDArray[np.float64]:
    """Convert formation energies to natural log association constants in M.

    ``delta_G`` is the standard energy for sum(a_i monomer_i) -> complex,
    including association and symmetry corrections exactly once, upstream.
    With n strands: log(K_M) = -delta_G/(R*T) + (1-n)*log(c_standard).
    A 1 M standard uses c_standard=1. A mole-fraction standard at a specified
    water concentration can use, for example, c_standard=55.1 explicitly.
    This conversion adds no other association or symmetry correction.
    """
    matrix = _stoichiometry(stoichiometry)
    energy = _vector(formation_free_energies_kcal_mol, matrix.shape[1], "formation energies")
    if np.any(np.isnan(energy)) or np.any(np.isneginf(energy)):
        raise ValueError("Formation energies must be finite or +inf for an impossible complex.")
    temperature = _real_scalar(temperature_kelvin, "temperature_kelvin")
    standard = _real_scalar(standard_concentration_molar, "standard_concentration_molar")
    if not np.isfinite(temperature) or temperature <= 0:
        raise ValueError("temperature_kelvin must be finite and positive.")
    if not np.isfinite(standard) or standard <= 0:
        raise ValueError("standard_concentration_molar must be finite and positive.")
    with np.errstate(over="ignore"):
        constants = -energy / (GAS_CONSTANT_KCAL_MOL_K * temperature)
        constants += (1.0 - matrix.sum(axis=0)) * np.log(standard)
    if np.any(np.isposinf(constants)) or np.any(np.isnan(constants)):
        raise ValueError("Formation energy conversion exceeds the supported numerical range.")
    if np.any(np.isneginf(constants) & np.isfinite(energy)):
        raise ValueError("Formation energy conversion exceeds the supported numerical range.")
    return constants


def solve_equilibrium(
    stoichiometry: ArrayLike,
    total_concentrations_molar: ArrayLike,
    *,
    log_equilibrium_constants: ArrayLike | None = None,
    formation_free_energies_kcal_mol: ArrayLike | None = None,
    temperature_kelvin: float = 310.15,
    standard_concentration_molar: float = 1.0,
    relative_tolerance: float = 1e-10,
    max_iterations: int = 1000,
) -> EquilibriumResult:
    """Solve mass action and all strand conservation equations.

    Matrix shape is (strand types, complexes), totals are strand concentrations
    in M, and exactly one of the two energy/constant arguments is required.
    Log constants are natural logarithms of the numerical value of K expressed
    in M**(1-n); they have already undergone any standard-state conversion.
    Thermodynamic conversion options apply only to formation-energy input.

    Conservation is checked against each nonzero strand total separately, so a
    dilute strand cannot hide behind the tolerance of a more abundant strand.
    Failure raises EquilibriumConvergenceError, never a partial or zero result.
    """
    matrix = _stoichiometry(stoichiometry)
    totals = _vector(total_concentrations_molar, matrix.shape[0], "strand totals")
    if not np.all(np.isfinite(totals)) or np.any(totals < 0):
        raise ValueError("Strand totals must be finite nonnegative concentrations in mol/L.")
    if (log_equilibrium_constants is None) == (formation_free_energies_kcal_mol is None):
        raise ValueError("Provide exactly one of log_equilibrium_constants or formation_free_energies_kcal_mol.")
    if formation_free_energies_kcal_mol is not None:
        constants = formation_energies_to_log_constants(
            matrix, formation_free_energies_kcal_mol,
            temperature_kelvin=temperature_kelvin,
            standard_concentration_molar=standard_concentration_molar,
        )
    else:
        constants = _vector(log_equilibrium_constants, matrix.shape[1], "log equilibrium constants")
        if np.any(np.isnan(constants)) or np.any(np.isposinf(constants)):
            raise ValueError("Log equilibrium constants must be finite or -inf for an impossible complex.")
        temperature = _real_scalar(temperature_kelvin, "temperature_kelvin")
        standard = _real_scalar(standard_concentration_molar, "standard_concentration_molar")
        if temperature != 310.15 or standard != 1.0:
            raise ValueError("Temperature and standard concentration only apply to formation-energy input; log constants are already converted.")
    monomers = _monomer_indices(matrix)
    if np.any(constants[monomers] != 0):
        raise ValueError("Monomer formation energies / log equilibrium constants must be zero.")
    tolerance = _real_scalar(relative_tolerance, "relative_tolerance")
    if not np.isfinite(tolerance) or not 100 * np.finfo(float).eps <= tolerance < 1:
        raise ValueError("relative_tolerance must be finite, at least 100 machine epsilons, and less than 1.")
    if isinstance(max_iterations, (bool, np.bool_)) or not isinstance(max_iterations, (int, np.integer)) or max_iterations < 1:
        raise ValueError("max_iterations must be a positive integer.")

    present = totals > 0
    allowed = np.isfinite(constants) & ~np.any(matrix[~present] > 0, axis=0)
    log_concentrations = np.full(matrix.shape[1], -np.inf)
    evaluations = 0
    if np.any(present):
        active_matrix = matrix[np.ix_(present, allowed)]
        active_constants = constants[allowed]
        log_totals = np.log(totals[present])
        if np.all(active_matrix.sum(axis=0) == 1):
            log_concentrations[monomers[present]] = log_totals
        else:
            with np.errstate(divide="ignore"):
                log_stoichiometry = np.log(active_matrix)

            # Start with weak associations and continue to the exact supplied
            # constants. Starting directly at extreme K can erase free-monomer
            # contributions from the Jacobian and hide a limiting strand. This
            # homotopy changes only the numerical path, not the final model.
            trial_log = active_constants + active_matrix.T @ log_totals
            permitted_log = np.min(
                np.where(active_matrix > 0, log_totals[:, None] - log_stoichiometry, np.inf),
                axis=0,
            )
            association_order = active_matrix.sum(axis=0) - 1
            associating = association_order > 0
            shift = max(0.0, float(np.max(
                (trial_log[associating] - permitted_log[associating] + np.log(active_matrix.shape[1]))
                / association_order[associating]
            )))
            stages = int(np.ceil(shift / 10))
            initial = log_totals.copy()
            optimizer_tolerance = max(1e-14, min(1e-10, tolerance * 0.01))
            for stage in range(stages + 1):
                remaining = int(max_iterations) - evaluations
                if remaining < 1:
                    raise EquilibriumConvergenceError("Equilibrium solver exhausted its iteration budget before reaching the supplied constants.")
                stage_shift = shift * (1 - stage / stages) if stages else 0.0
                stage_constants = active_constants - stage_shift * association_order
                residual_and_jacobian = _balanced_residual(
                    active_matrix, totals[present], stage_constants, initial,
                )
                result = least_squares(
                    lambda value: residual_and_jacobian(value)[0], initial,
                    jac=lambda value: residual_and_jacobian(value)[1],
                    x_scale="jac", ftol=optimizer_tolerance,
                    xtol=optimizer_tolerance, gtol=optimizer_tolerance,
                    max_nfev=remaining,
                )
                evaluations += int(result.nfev)
                if not result.success or not np.all(np.isfinite(result.x)):
                    raise EquilibriumConvergenceError(f"Equilibrium solver did not converge: {result.message}")
                if np.max(np.abs(residual_and_jacobian(result.x)[0])) > tolerance:
                    raise EquilibriumConvergenceError("Equilibrium continuation failed logarithmic conservation validation.")
                initial = result.x
            log_concentrations[allowed] = active_constants + active_matrix.T @ initial

    with np.errstate(under="ignore", over="ignore"):
        concentrations = np.exp(log_concentrations)
    mass_residual = matrix @ concentrations - totals
    mass_error = float(np.max(np.abs(mass_residual[present]) / totals[present], initial=0.0))
    if not np.all(np.isfinite(concentrations)) or mass_error > tolerance or np.any(mass_residual[~present] != 0):
        raise EquilibriumConvergenceError(
            f"Equilibrium failed conservation validation: maximum relative error {mass_error:.3g}, requested {tolerance:.3g}."
        )
    if np.any(allowed):
        # Omit absent rows: 0 * (-inf) would otherwise contaminate this identity.
        mass_action_residual = (
            log_concentrations[allowed]
            - matrix[np.ix_(present, allowed)].T @ log_concentrations[monomers[present]]
            - constants[allowed]
        )
        action_error = float(np.max(np.abs(mass_action_residual), initial=0.0))
        rounding_allowance = 32 * np.finfo(float).eps * max(1.0, float(np.max(np.abs(constants[allowed]))))
        if not np.isfinite(action_error) or action_error > max(tolerance, rounding_allowance):
            raise EquilibriumConvergenceError("Equilibrium failed logarithmic mass-action validation.")
    else:
        action_error = 0.0
    for value in (concentrations, log_concentrations, mass_residual):
        value.flags.writeable = False
    return EquilibriumResult(concentrations, log_concentrations, mass_residual, mass_error, action_error, evaluations)


def _vector(value: ArrayLike, length: int, label: str) -> NDArray[np.float64]:
    array = _real_array(value, label)
    if array.shape != (length,):
        raise ValueError(f"{label} must have shape ({length},).")
    return array


def _stoichiometry(value: ArrayLike) -> NDArray[np.float64]:
    matrix = _real_array(value, "Stoichiometry")
    if matrix.ndim != 2 or min(matrix.shape) < 1:
        raise ValueError("Stoichiometry must be a nonempty (strand types, complexes) matrix.")
    if not np.all(np.isfinite(matrix)) or np.any(matrix < 0) or np.any(matrix != np.floor(matrix)):
        raise ValueError("Stoichiometry must contain finite nonnegative integers.")
    if np.any(matrix.sum(axis=0) == 0):
        raise ValueError("Each complex must contain at least one strand.")
    return matrix


def _real_array(value: ArrayLike, label: str) -> NDArray[np.float64]:
    # Preserve original scalar types: coercing a mixed list to float first
    # would turn True into 1, parse strings, or discard complex components.
    original = np.asarray(value, dtype=object)
    if any(isinstance(item, (bool, np.bool_)) or not isinstance(item, Real)
           for item in original.flat):
        raise ValueError(f"{label} must contain real numbers, without booleans, strings, or complex values.")
    try:
        return np.asarray(original, dtype=float)
    except (OverflowError, TypeError, ValueError) as error:
        raise ValueError(f"{label} values must be representable as real float64 numbers.") from error


def _real_scalar(value: Real, label: str) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{label} must be a real number.")
    try:
        return float(value)
    except (OverflowError, ValueError) as error:
        raise ValueError(f"{label} must be representable as a real float64 number.") from error


def _monomer_indices(matrix: NDArray[np.float64]) -> NDArray[np.int64]:
    indices = []
    for strand in range(matrix.shape[0]):
        columns = np.flatnonzero((matrix[strand] == 1) & (matrix.sum(axis=0) == 1))
        if len(columns) != 1:
            raise ValueError("Exactly one monomer complex is required for each strand type.")
        indices.append(int(columns[0]))
    return np.asarray(indices, dtype=np.int64)


def _balanced_residual(matrix, totals, constants, initial):
    """Conservation in a basis that cancels dominant complexes exactly.

    For a strongly bound AB at equal A/B totals, separate ordinary mass sums
    lose all information about the tiny free monomers. Subtracting their two
    conservation equations *symbolically* instead gives A_free = B_free.
    The same change of basis works for arbitrary multistrand compositions.
    Fractions are used only for this small linear transformation; concentration
    optimization remains float64. Positive and negative contributions are
    subsequently compared in log space, never by subtracting near-equal sums.
    """
    strand_count, complex_count = matrix.shape
    estimated_log_c = constants + matrix.T @ initial
    order = np.argsort(-estimated_log_c, kind="stable")
    basis = []
    for column in order:
        candidate = basis + [int(column)]
        if np.linalg.matrix_rank(matrix[:, candidate]) == len(candidate):
            basis = candidate
        if len(basis) == strand_count:
            break
    if len(basis) != strand_count:
        raise EquilibriumConvergenceError("Unable to construct an independent conservation basis.")

    # Exact inverse of the selected integer stoichiometry columns. Inexact
    # elimination could leave spurious 1e-16 * [dominant complex] terms that
    # overwhelm free monomers many orders of magnitude below that level.
    augmented = [
        [Fraction(int(matrix[row, col])) for col in basis]
        + [Fraction(int(row == col)) for col in range(strand_count)]
        for row in range(strand_count)
    ]
    for col in range(strand_count):
        pivot = next(row for row in range(col, strand_count) if augmented[row][col])
        augmented[col], augmented[pivot] = augmented[pivot], augmented[col]
        divisor = augmented[col][col]
        augmented[col] = [value / divisor for value in augmented[col]]
        for row in range(strand_count):
            if row != col:
                factor = augmented[row][col]
                augmented[row] = [value - factor * reference for value, reference in zip(augmented[row], augmented[col])]
    inverse = [row[strand_count:] for row in augmented]
    log_positive = np.full((strand_count, complex_count + 1), -np.inf)
    log_negative = np.full_like(log_positive, -np.inf)
    for row in range(strand_count):
        coefficients = [
            sum(inverse[row][strand] * int(matrix[strand, col]) for strand in range(strand_count))
            for col in range(complex_count)
        ]
        # Place the transformed total on the opposite side of the equality.
        coefficients.append(-sum(inverse[row][strand] * Fraction(float(totals[strand])) for strand in range(strand_count)))
        for col, value in enumerate(coefficients):
            if value:
                target = log_positive if value > 0 else log_negative
                target[row, col] = math.log(abs(value.numerator)) - math.log(value.denominator)

    def residual_and_jacobian(log_free):
        log_c = np.append(constants + matrix.T @ log_free, 0.0)
        positive = log_positive + log_c
        negative = log_negative + log_c
        log_plus = logsumexp(positive, axis=1)
        log_minus = logsumexp(negative, axis=1)
        derivative = (
            np.exp(positive[:, :-1] - log_plus[:, None])
            - np.exp(negative[:, :-1] - log_minus[:, None])
        ) @ matrix.T
        return log_plus - log_minus, derivative

    return residual_and_jacobian
