"""Independent ViennaRNA reference for native RNA regression tests."""
from collections import Counter
from dataclasses import asdict
import math

from aptaswitch_core.thermodynamics.structures import format_structure


def direct_vienna(spec, conditions, count, seed):
    """Public Vienna API, with exact PF salt factors and explicit integer MFE."""
    import RNA
    from aptaswitch_core.thermodynamics.runtime_identity import identify_vienna_runtime
    if RNA.__version__ != "2.7.2":
        raise RuntimeError("The reference requires ViennaRNA 2.7.2")
    RNA.params_load_RNA_Turner2004()
    md = RNA.md(temperature=conditions.celsius, salt=1.021, betaScale=1., pf_smooth=1,
                dangles=0, special_hp=1, noLP=0, noGU=0, noGUclosure=0, logML=0,
                circ=0, circ_penalty=1, gquad=0, uniq_ML=1, energy_set=0,
                backtrack=1, backtrack_type="F", compute_bpp=1, max_bp_span=-1,
                min_loop_size=3, window_size=-1, oldAliEn=0, ribo=0,
                cv_fact=1., nc_fact=1., sfact=1.07, saltMLLower=6,
                saltMLUpper=24, saltDPXInit=99999)
    fc = RNA.fold_compound(spec.sequence, md)
    fc.exp_params_rescale(0.)
    thermal = fc.exp_params.kT / 1000
    salt = -.114 * math.log(conditions.sodium) * conditions.kelvin / 310.15
    factor = math.exp(-salt / thermal)
    closes_pair = (RNA.DECOMP_PAIR_HP, RNA.DECOMP_PAIR_IL, RNA.DECOMP_PAIR_ML)
    fc.sc_add_exp_f(lambda i, j, k, l, decomposition, data: factor if decomposition in closes_pair else 1.)
    # Official MFE is a separate, rounded objective; never substitute it for
    # the native minimum of continuous PF-state energies.
    integer_salt = math.floor(100 * salt + .5)
    fc.sc_add_f(lambda i, j, k, l, decomposition, data: integer_salt if decomposition in closes_pair else 0)
    mfe_structure, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)
    _, free_energy = fc.pf()
    upper, n = fc.bpp(), spec.n_nucleotides
    matrix = [[0.] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            matrix[i][j] = matrix[j][i] = upper[i + 1][j + 1]
    for i in range(n):
        matrix[i][i] = 1 - math.fsum(matrix[i])
    RNA.init_rand(seed)
    samples = list(fc.pbacktrack(count)) if count else []
    result = dict(log_partition=-free_energy / thermal, free_energy_kcal_mol=free_energy,
                thermal_energy_kcal_mol=thermal, matrix=matrix, sample_counts=dict(Counter(samples)),
                minimum_energy_kcal_mol=mfe, minimum_structure=mfe_structure,
                minimum_objective="official integer-centikcal MFE; per-pair salt rounded to centikcal",
                rna_runtime=asdict(identify_vienna_runtime(RNA)), pair_salt_energy_kcal_mol=salt)
    from aptaswitch_core.thermodynamics.hybrid_motifs import HybridMotifModel
    model = HybridMotifModel(sodium_correction="published-pair")
    result["continuous_pf_energy_of_official_minimum_kcal_mol"] = model.energy(spec, conditions, mfe_structure)
    if n <= 18:
        # Independent traversal of all structures, evaluated with Vienna's
        # public integer energy routine and explicit integer salt contribution.
        # A second compound has no callbacks, preventing double-counting salt.
        from aptaswitch_core.thermodynamics.enumeration import enumerate_structures
        RNA.params_load_RNA_Turner2004()
        evaluator = RNA.fold_compound(spec.sequence, md)
        energies = {}
        for state in enumerate_structures(spec, pair_allowed=model.pair_allowed, min_hairpin_size=3):
            structure = format_structure(state, (n,))
            pairs = sum(j > i for i, j in enumerate(state.partners))
            energy = round(100 * evaluator.eval_structure(structure)) + pairs * integer_salt
            energies[structure] = energy
        minimum = min(energies.values())
        result["integer_exhaustive_minimum"] = {
            "energy_kcal_mol": minimum / 100., "number_of_minimum_structures": sum(e == minimum for e in energies.values()),
            "official_witness_is_minimum": energies[mfe_structure] == minimum,
            "absolute_error_kcal_mol": abs(mfe - minimum / 100.), "number_of_structures": len(energies)}
        result["integer_structure_energies_centikcal"] = energies
    return result
