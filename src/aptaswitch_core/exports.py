"""Run export helpers."""

from __future__ import annotations

from .localization import tr

import csv
import json
from dataclasses import asdict
from pathlib import Path
from html import escape

from .models import DesignCandidate, RunResult, FunctionalDefectWeights
from .scoring import candidate_score_100


CSV_HEADERS = [
    "Rank",
    "Trial",
    "Score_selection",
    "Defect",
    "ON_yield_pct",
    "Leak_pct",
    "dG_OFF_kcal_mol",
    "dG_ON_kcal_mol",
    "ddG_activation_kcal_mol",
    "dG_RBS_linker_kcal_mol",
    "Bad_RBS_linker",
    "STOP_codon",
    "First_STOP",
    "Reporter",
    "Trigger_DNA_5to3",
    "Switch_RNA_5to3",
    "Status",
    "Exclusion_reason",
    "Defect_structural",
    "Defect_concentration",
    "Defect_regional_balanced",
    "Defect_OFF_stem",
    "Defect_OFF_toehold",
    "Defect_ON_RBS_linker",
    "Defect_OFF",
    "Defect_ON",
    "Defect_metric_for_score",
    "Defect_ON_RBS_AUG",
    "Functional_weight_OFF_toehold",
    "Functional_weight_OFF_stem",
    "Functional_weight_ON_RBS_AUG",
    "Functional_weight_ON_RBS_linker",
    "Defect_metric_scheme",
    "Defect_OFF_stem_before_AUG",
    "Defect_OFF_stem_after_AUG",
    "Score_raw",
    "Score_scale",
    "Trigger_linearity_defect",
    "Trigger_unpaired_fraction",
    "Weight_trigger_linearity",
    "Thermodynamics_engine",
    "Thermodynamics_model",
    "Thermodynamics_parameter_digest",
    "ON_RBS_AUG_joint_unpaired_probability",
    "ON_RBS_AUG_opening_kcal_mol",
    "ON_RBS_AUG_opening_impossible",
    "OSTIR_OFF_reference_rate",
    "OSTIR_ON_constrained_rate",
    "OSTIR_log10_ON_OFF",
    "OSTIR_exploratory_merit",
    "Weight_rbs_calculator",
]


def _export_thermodynamics(result: RunResult) -> dict:
    """Identify historical rows without inventing a parameter fingerprint."""
    return {"engine": result.request.engine, "model": result.request.thermo.material,
            **result.thermodynamics}


def candidate_row(rank: int, candidate: DesignCandidate, *,
                  functional_defect_weights: FunctionalDefectWeights | None = None,
                  score_scale: str = "bounded_v1", trigger_linearity_weight: float | None = None,
                  rbs_calculator_weight: float | None = None,
                  thermodynamics: dict | None = None) -> list:
    regional = candidate.regional_defects
    functional_weights = functional_defect_weights or FunctionalDefectWeights(**regional.get("weights", {}))
    stem_parts = regional.get("off_stem_parts", {})
    thermodynamics = thermodynamics or {}
    accessibility = getattr(candidate, "rbs_accessibility", {})
    rbs = getattr(candidate, "rbs_calculator", {})
    return [
        "" if candidate.excluded else rank,
        candidate.trial,
        candidate_score_100(candidate, score_scale=score_scale),
        candidate.defect,
        candidate.on_yield_pct,
        candidate.leak_pct,
        candidate.delta_g_off,
        candidate.delta_g_on,
        candidate.ddg_activation,
        candidate.delta_g_rbs_linker,
        "" if candidate.has_bad_rbs_linker is None else ("yes" if candidate.has_bad_rbs_linker else "no"),
        "yes" if candidate.has_stop_codon else "no",
        candidate.first_stop or "",
        candidate.reporter_label,
        candidate.trigger_seq,
        candidate.switch_seq,
        "excluded" if candidate.excluded else "evaluated",
        tr(candidate.exclusion_reason),
        candidate.defect_structural,
        candidate.defect_concentration,
        regional.get("balanced"),
        regional.get("off_stem"),
        regional.get("off_toehold"),
        regional.get("on_rbs_linker"),
        regional.get("off"),
        regional.get("on"),
        "functional",
        regional.get("on_rbs_aug"),
        functional_weights.off_toehold,
        functional_weights.off_stem,
        functional_weights.on_rbs_aug,
        functional_weights.on_rbs_linker,
        regional.get("scheme"),
        stem_parts.get("before_aug", {}).get("normalized_defect"),
        stem_parts.get("after_aug", {}).get("normalized_defect"),
        candidate.score_raw,
        "bounded_v1",
        candidate.trigger_linearity.get("normalized_defect"),
        candidate.trigger_linearity.get("unpaired_fraction"),
        trigger_linearity_weight,
        thermodynamics.get("engine"), thermodynamics.get("model"),
        thermodynamics.get("parameter_digest"),
        accessibility.get("probability"),
        accessibility.get("opening_free_energy_kcal_mol"),
        accessibility.get("constrained_ensemble_empty"),
        rbs.get("off_rate"),
        rbs.get("on_rate"),
        rbs.get("log10_on_off"),
        rbs.get("merit"),
        rbs_calculator_weight,
    ]


def inverse_comparison_table(result: RunResult) -> tuple[list[str], list[list]]:
    from .inverse_extension import ranked_inverse_candidates
    columns = ["Trial", "Original_switch_rank", "Switch_score", "Conservation_impact",
               "Structure_rank", "Compromise_rank", "Compromise_score", "Performance_weight",
               "Pareto", "Trigger_5to3", "Extension_5prime", "Extension_3prime",
               "Original_aptamer_trigger_2D", "Extended_aptamer_trigger_2D"]
    rows = []
    for rank, row in enumerate(ranked_inverse_candidates(result), 1):
        metrics = row["metrics"]
        rows.append([row["trial"], row["original_switch_rank"], row["switch_score"], row["conservation_score"],
                     row["structure_rank"], rank, row["compromise_score"], row["performance_weight"], row["pareto"],
                     metrics["extended_trigger"], metrics["extension_5prime"], metrics["extension_3prime"],
                     metrics["reference_mfe_structure"], metrics["mfe_structure"]])
    return columns, rows


def write_json(result: RunResult, output_dir: Path) -> Path:
    path = output_dir / f"{result.run_id}_results.json"
    path.write_text(json.dumps(result.to_dict(), indent=2, allow_nan=False), encoding="utf-8")
    return path


def write_csv(result: RunResult, output_dir: Path) -> Path:
    path = output_dir / f"{result.run_id}_ranked.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(CSV_HEADERS)
        for rank, candidate in enumerate(result.candidates, start=1):
            writer.writerow(candidate_row(rank, candidate,
                                          functional_defect_weights=result.request.functional_defect_weights,
                                          score_scale=result.request.score_scale,
                                          trigger_linearity_weight=result.request.scoring.trigger_linearity,
                                          rbs_calculator_weight=getattr(result.request.scoring, "rbs_calculator", None),
                                          thermodynamics=_export_thermodynamics(result)))
    return path


def write_xlsx(result: RunResult, output_dir: Path) -> Path | None:
    try:
        import openpyxl
        from openpyxl.styles import Alignment, Font, PatternFill
        from openpyxl.utils import get_column_letter
    except ImportError:
        return None

    path = output_dir / f"{result.run_id}_ranked.xlsx"
    workbook = openpyxl.Workbook()
    sheet = workbook.active
    sheet.title = tr('Candidats AptaSwitch')
    sheet.append(CSV_HEADERS)
    for rank, candidate in enumerate(result.candidates, start=1):
        sheet.append(candidate_row(rank, candidate,
                                   functional_defect_weights=result.request.functional_defect_weights,
                                   score_scale=result.request.score_scale,
                                   trigger_linearity_weight=result.request.scoring.trigger_linearity,
                                   rbs_calculator_weight=getattr(result.request.scoring, "rbs_calculator", None),
                                   thermodynamics=_export_thermodynamics(result)))

    header_fill = PatternFill(start_color="1F2937", end_color="1F2937", fill_type="solid")
    header_font = Font(color="FFFFFF", bold=True)
    for cell in sheet[1]:
        cell.fill = header_fill
        cell.font = header_font
        cell.alignment = Alignment(horizontal="center")

    for idx, header in enumerate(CSV_HEADERS, start=1):
        width = 100 if header == "Switch_RNA_5to3" else 36 if header == "Trigger_DNA_5to3" else max(12, min(36, len(header) + 3))
        sheet.column_dimensions[get_column_letter(idx)].width = width
    sheet.freeze_panes = "A2"
    sheet.auto_filter.ref = sheet.dimensions
    definitions = workbook.create_sheet("Defect metrics")
    definitions.append(["Field", "Definition", "Unit / interpretation"])
    definitions.append(["Defect_metric_for_score", result.request.defect_metric,
                        "functional: weighted defects of the four functional switch regions"])
    definitions.append(["Defect", "Recorded global normalized tube ensemble defect under the identified model", "fraction, 0–1; retained in its original column"])
    definitions.append(["Defect_structural", "Structural contribution to global tube defect, weighted by available target-complex concentration", "fraction, 0–1; not the defect of an isolated complex"])
    definitions.append(["Defect_concentration", "Target-complex concentration shortfall contribution to global tube defect", "fraction, 0–1"])
    definitions.append(["Defect_regional_balanced", "sum(functional_weight * functional_defect) / sum(functional_weights)", "fraction, 0–1; default 1/1/1/0.75 gives 26.67% to OFF toehold, OFF stem and ON RBS–AUG, and 20% to ON RBS–linker"])
    definitions.append(["Defect_OFF_stem", "Mean target-pair mismatch probability over both arms of the stems immediately before and after AUG", "fraction, 0–1; 6 bp before and 9 bp after AUG in the classic 30-nt architecture"])
    definitions.append(["Defect_OFF_toehold", "Mean target-state mismatch probability over OFF toehold bases", "fraction, 0–1; ensemble probabilities"])
    definitions.append(["Defect_ON_RBS_linker", "Mean target-state mismatch probability over ON RBS-to-linker bases", "fraction, 0–1; ensemble probabilities"])
    definitions.append(["Defect_ON_RBS_AUG", "Mean probability of being paired from the RBS-loop prefix through the end of AUG in ON", "fraction, 0–1; expected state is unpaired; functional_v2 only"])
    definitions.append(["Defect_OFF_stem_before_AUG", "Mean target-pair mismatch probability over both arms of the stem immediately before AUG in OFF", "fraction, 0–1; informative component of Defect_OFF_stem"])
    definitions.append(["Defect_OFF_stem_after_AUG", "Mean target-pair mismatch probability over both arms of the stem immediately after AUG in OFF", "fraction, 0–1; informative component of Defect_OFF_stem"])
    definitions.append(["Defect_metric_scheme", "Recorded definition of measured functional regions", "functional_v2 identifies the four functional regions"])
    definitions.append(["Defect_OFF", "Weighted mean of off_stem and off_toehold using their saved weights", "fraction, 0–1; blank if both OFF weights are zero"])
    definitions.append(["Defect_ON", "Weighted mean of on_rbs_aug and on_rbs_linker using their saved weights", "fraction, 0–1; blank if both ON weights are zero"])
    zone_weights = result.request.functional_defect_weights
    weight_label = "Functional"
    for key, fraction in zone_weights.fractions().items():
        definitions.append([f"{weight_label} weight: {key}", getattr(zone_weights, key),
                            f"relative weight; {fraction:.6%} of the selected compromise"])
    definitions.append(["Overlapping ON regions", "on_rbs_aug is contained in on_rbs_linker; bases in this shared interval contribute to both criteria", "Intentional additional weight on the RBS-to-AUG initiation region"])
    definitions.append(["OFF sequestration", "Stem integrity around AUG and its opposing bulge is used as a structural proxy", "AUG and the bulge are not required to pair in the target structure; this metric does not prove translational repression"])
    definitions.append(["Score_selection", "50 + (100/pi)*atan(Score_raw/100)", "bounded comparison score, 0–100; a monotonic scale, not a probability of experimental success"])
    definitions.append(["Score_raw", "w_on*ON - w_leak*Leak - w_defect*100*selected_defect - w_ddg*ddG + w_rbs*dG_RBS_linker + 100*w_rbs_calculator*OSTIR_exploratory_merit - w_trigger_linearity*100*Trigger_linearity_defect", "The trigger-linearity term applies only to switch-first inverse extension; selected_defect follows Defect_metric_for_score. A missing OSTIR measurement must not be reconstructed or silently added to a historical score."])
    definitions.append(["Score_scale", "bounded_v1", "Bounded display scale; Score_raw retains the comparison sum for exact ranking"])
    definitions.append(["Trigger_linearity_defect", "Mean probability of being paired in the isolated trigger ensemble: 1 - mean(P_unpaired)", "fraction, 0–1; lower is more accessible; calculated only for inverse extension under the saved switch-design conditions"])
    definitions.append(["Trigger_unpaired_fraction", "Mean probability of being unpaired in the isolated trigger ensemble", "fraction, 0–1; not the probability of an entirely linear conformation"])
    definitions.append(["Weight_trigger_linearity", result.request.scoring.trigger_linearity, "relative penalty weight applied to Score_raw only during switch-first inverse extension"])
    definitions.append(["Leak_pct", "100 * [complex containing exactly 1 trigger + 1 switch] / initial [switch] in the aptamer + trigger + switch tube", f"Association percentage, not OFF translation leakage; aptamer concentration equals trigger concentration. Competing complexes up to {result.request.thermo.validation_max_complex_size} strands contribute to equilibrium, but only trigger-switch dimers count in the numerator."])
    definitions.append(["OSTIR_OFF_reference_rate", "OSTIR prediction for the actual RNA switch around its intended AUG, without structural constraints", "OSTIR arbitrary units; RNA reference prediction, not measured translation or protein production"])
    definitions.append(["OSTIR_ON_constrained_rate", "OSTIR prediction with only the intended trigger-occupied sites unavailable for intramolecular RNA pairing", "Ideal occupied RNA projection; no DNA binding energy, tube concentration or measured hybrid initiation rate"])
    definitions.append(["OSTIR_log10_ON_OFF", "log10(OSTIR_ON_constrained_rate / OSTIR_OFF_reference_rate)", "Exploratory RNA-model contrast; 1 means a predicted factor of 10, not an experimentally validated induction ratio"])
    definitions.append(["OSTIR_exploratory_merit", "on/(on+K) * max(0, (on-off)/(on+off)); K=exp(7.279194329)", "Bounded ranking heuristic, 0–1. K is a normalization pivot, not a validated ON target; larger is preferred."])
    definitions.append(["Weight_rbs_calculator", getattr(result.request.scoring, "rbs_calculator", None), "Relative bonus weight: +100 * weight * OSTIR_exploratory_merit; not an experimentally established optimal weight"])
    definitions.append(["OSTIR reference model", "OSTIR 1.1.3 / ViennaRNA 2.7.2; RNA-Turner2004; E. coli anti-Shine-Dalgarno; 37 C; 1.021 M monovalent salt; 0 M magnesium", "Actual candidate RNA in the native +/-35 nt window around its intended AUG, clipped at sequence ends; distinct from user-selected hybrid tube conditions. Recorded per-candidate parameters are in RBS Calculator model and JSON."])
    definitions.append(["ON_RBS_AUG_joint_unpaired_probability", "Z(window unpaired) / Z(unrestricted connected ON)", "All bases from the RBS-loop prefix through AUG simultaneously unpaired; conditional on the DNA–RNA complex already formed, not a product of marginals"])
    definitions.append(["ON_RBS_AUG_opening_kcal_mol", "-RT ln(P_joint_unpaired)", "Informative opening cost; not a kinetic barrier and not an additional score term"])
    definitions.append(["ON_RBS_AUG_opening_impossible", "True only for an empty constrained ensemble", "Probability=0 and blank opening energy mean +infinity when true; numerical underflow retains a finite opening cost when false"])
    definitions.append(["RBS stability during optimization", "global_defect + w_rbs * max(0, -dG_RBS_linker) / 100 + w_rbs_calculator * (1 - OSTIR_exploratory_merit)", "Refinement only within fixed local-defect and tube guards; isolated RNA MFE, no hard -8 kcal/mol exclusion. OSTIR adds an exploratory preference; check recorded optimization settings for historical runs."])
    definitions.append(["Missing metrics", "Blank cells mean unavailable, not calculated or a state with zero weight; they never mean zero", "Excluded candidates have no computed score or defect"])
    definitions.append(["Scope", "Ensemble metrics differ from the single displayed MFE structure", "Neither metric guarantees translation or experimental performance"])
    model_sheet = workbook.create_sheet("Thermodynamic model")
    model_sheet.append(["Field", "Recorded value"])
    for key, value in _export_thermodynamics(result).items():
        model_sheet.append([key, json.dumps(value, allow_nan=False) if isinstance(value, (dict, list)) else value])
    model_sheet.column_dimensions["A"].width = 28
    model_sheet.column_dimensions["B"].width = 110
    rbs_candidates = [(candidate.trial, getattr(candidate, "rbs_calculator", {}))
                      for candidate in result.candidates if getattr(candidate, "rbs_calculator", {})]
    if rbs_candidates:
        rbs_sheet = workbook.create_sheet("RBS Calculator model")
        rbs_sheet.append(["Trial", "Field", "Recorded value"])
        for trial, measurement in rbs_candidates:
            for key in ("method", "reference", "window_start", "window_end", "window_sequence",
                        "window_aug_start", "window_occupied_indices", "on_constraints", "scope", "normalization"):
                if key in measurement:
                    value = measurement[key]
                    rbs_sheet.append([trial, key, json.dumps(value, allow_nan=False)
                                      if isinstance(value, (dict, list)) else value])
        rbs_sheet.column_dimensions["A"].width = 10
        rbs_sheet.column_dimensions["B"].width = 30
        rbs_sheet.column_dimensions["C"].width = 110
        rbs_sheet.freeze_panes = "A2"
    if result.optimization:
        optimization_sheet = workbook.create_sheet("Optimization")
        optimization_sheet.append(["Trial", "Stop reason", "Budget kind", "Best defect", "Evaluations", "Proposals", "Seed", "Options"])
        for outcome in result.optimization.get("trials", []):
            optimization_sheet.append([outcome.get("trial"), outcome.get("stop_reason"), outcome.get("budget_kind"),
                                       outcome.get("best_loss"), outcome.get("evaluations"), outcome.get("proposals"),
                                       outcome.get("options", {}).get("seed"), json.dumps(outcome.get("options", {}), allow_nan=False)])
        optimization_sheet.freeze_panes = "A2"
        for column in "ABCDEFG":
            optimization_sheet.column_dimensions[column].width = 20
        optimization_sheet.column_dimensions["H"].width = 95
    definitions.column_dimensions["A"].width = 30
    definitions.column_dimensions["B"].width = 95
    definitions.column_dimensions["C"].width = 72
    definitions.freeze_panes = "A2"
    for row in definitions:
        for cell in row:
            cell.alignment = Alignment(vertical="top", wrap_text=True)
    for cell in definitions[1]:
        cell.fill, cell.font = header_fill, header_font
    if result.request.inverse_extension is not None:
        comparison = workbook.create_sheet(tr("Comparaison des extensions"))
        columns, rows = inverse_comparison_table(result)
        comparison.append(columns)
        for row in rows:
            comparison.append(row)
        comparison.freeze_panes = "A2"
        comparison.auto_filter.ref = comparison.dimensions
        metadata = workbook.create_sheet(tr("Conditions de comparaison"))
        metadata.append(["Field", "Value"])
        for key, value in result.inverse_extension_analysis.get("model", {}).items():
            metadata.append([key, value])
        metadata.append(["status", result.inverse_extension_analysis.get("status")])
        metadata.append(["ranking", "weighted normalized ranks; default performance weight 0.5; lower is better"])
    workbook.save(path)
    return path


def write_top_svg(result: RunResult, output_dir: Path) -> Path | None:
    top = next((candidate for candidate in result.candidates if not candidate.excluded), None)
    if top is None:
        return None
    path = output_dir / f"{result.run_id}_top1.svg"
    trigger_width = max(180, min(900, len(top.trigger_seq) * 14))
    switch_width = max(260, min(1200, len(top.switch_seq) * 7))
    width = max(900, switch_width + 120)
    height = 380
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
  <rect width="100%" height="100%" fill="#f8fafc"/>
  <text x="48" y="48" font-family="Arial" font-size="24" font-weight="700" fill="#111827">{escape(tr('Conception Apta2Switch-Studio la mieux classée'))}</text>
  <text x="48" y="78" font-family="Arial" font-size="14" fill="#374151">{escape(tr('Score {score:.4f}/100 | ON {on:.2f}% | fuite {leak:.2f}% | ΔG RBS–linker {dg:.2f} kcal/mol', score=candidate_score_100(top, score_scale=result.request.score_scale), on=top.on_yield_pct, leak=top.leak_pct, dg=top.delta_g_rbs_linker))}</text>
  <rect x="48" y="118" width="{trigger_width}" height="34" rx="4" fill="#2563eb"/>
  <text x="60" y="141" font-family="Arial" font-size="13" font-weight="700" fill="#ffffff">{escape(tr('Trigger ADN 5′ vers 3′'))}</text>
  <text x="48" y="176" font-family="Menlo, Consolas, monospace" font-size="12" fill="#111827">{top.trigger_seq}</text>
  <rect x="48" y="218" width="{switch_width}" height="34" rx="4" fill="#059669"/>
  <text x="60" y="241" font-family="Arial" font-size="13" font-weight="700" fill="#ffffff">{escape(tr('Toehold switch ARN 5′ vers 3′'))}</text>
  <text x="48" y="276" font-family="Menlo, Consolas, monospace" font-size="11" fill="#111827">{top.switch_seq}</text>
  <text x="48" y="330" font-family="Arial" font-size="12" fill="#4b5563">{escape(tr('Ce SVG est un aperçu compact du résultat, pas une simulation de structure.'))}</text>
</svg>
"""
    path.write_text(svg, encoding="utf-8")
    return path


def write_run_log(result: RunResult, output_dir: Path) -> Path:
    path = output_dir / f"{result.run_id}.log"
    lines = [
        tr('Identifiant du run : {v0}', v0=result.run_id),
        tr('Statut : {v0}', v0=tr(result.status)),
        tr('Début : {v0}', v0=result.started_at),
        tr('Fin : {v0}', v0=result.completed_at),
        tr('Moteur : {v0}', v0=result.request.engine),
        tr('Essais demandés : {v0}', v0=result.request.trials),
        tr('Candidats : {v0}', v0=len(result.candidates)),
        tr('Exclure les candidats avec un STOP prématuré : {v0}', v0=result.request.exclude_stop_candidates),
        tr('Exclure les candidats au RBS–linker non linéaire : {v0}', v0=result.request.exclude_non_linear_candidates),
        tr('Candidats exclus (analyses restantes ignorées) : {v0}', v0=sum((candidate.excluded for candidate in result.candidates))),
        tr('Concentration du trigger (M) : {v0}', v0=result.request.thermo.resolved_trigger_concentration_m),
        tr('Concentration du toehold switch (M) : {v0}', v0=result.request.thermo.resolved_switch_concentration_m),
        f"Defect_metric_for_score: {result.request.defect_metric}",
        f"Score_scale_recorded: {result.request.score_scale}",
        "Score_selection_export: bounded_v1; 50 + (100/pi)*atan(Score_raw/100); 0–100 comparison score, not a probability",
        "Score_raw: w_on*ON - w_leak*Leak - w_defect*100*selected_defect - w_ddg*ddG + w_rbs*dG_RBS_linker + 100*w_rbs_calculator*OSTIR_exploratory_merit - w_trigger_linearity*100*Trigger_linearity_defect (last term only for inverse extension; missing historical OSTIR measurements are not rescored)",
        f"Weight_trigger_linearity: {result.request.scoring.trigger_linearity}",
        f"Weight_rbs_calculator: {getattr(result.request.scoring, 'rbs_calculator', None)}",
        f"Trigger_linearity_applies: {result.request.inverse_extension is not None}",
    ]
    lines.extend([
        f"Thermodynamics: {json.dumps(_export_thermodynamics(result), allow_nan=False)}",
        f"Optimization: {json.dumps(result.optimization, allow_nan=False)}",
        "Defect_metric_scheme: functional_v2",
        f"Functional_defect_weights: {json.dumps(asdict(result.request.functional_defect_weights))}",
    ])
    if result.warnings:
        lines.append(tr('Avertissements :'))
        lines.extend(f"- {warning}" for warning in result.warnings)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_project_file(result: RunResult, output_dir: Path) -> Path:
    path = output_dir / ".aptaswitch.json"
    data = {
        "format": "aptaswitch-studio-project-v1",
        "run_id": result.run_id,
        "request": result.request.to_dict(),
        "thermodynamics": result.thermodynamics,
        "optimization": result.optimization,
        "results_file": f"{result.run_id}_results.json",
        "ranked_csv": f"{result.run_id}_ranked.csv",
        "ranked_xlsx": f"{result.run_id}_ranked.xlsx",
        "top_svg": f"{result.run_id}_top1.svg" if any(not candidate.excluded for candidate in result.candidates) else None,
    }
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")
    return path


def export_run(result: RunResult) -> dict[str, str]:
    output_dir = Path(result.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "json": str(write_json(result, output_dir)),
        "csv": str(write_csv(result, output_dir)),
        "log": str(write_run_log(result, output_dir)),
        "project": str(write_project_file(result, output_dir)),
    }
    if result.request.inverse_extension is not None:
        comparison = output_dir / f"{result.run_id}_inverse_extension.csv"
        with comparison.open("w", newline="", encoding="utf-8-sig") as handle:
            writer = csv.writer(handle)
            columns, rows = inverse_comparison_table(result)
            writer.writerow(columns)
            writer.writerows(rows)
        paths["inverse_csv"] = str(comparison)
    xlsx = write_xlsx(result, output_dir)
    if xlsx:
        paths["xlsx"] = str(xlsx)
    svg = write_top_svg(result, output_dir)
    if svg:
        paths["svg"] = str(svg)
    return paths
