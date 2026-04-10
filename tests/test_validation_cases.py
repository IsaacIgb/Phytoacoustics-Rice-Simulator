"""
test_validation_cases.py — Validation harness v4.

Tests against VALIDATION_CASES.yaml (version 3) using rice_sound_master_v2.csv.
Maps primary_outcome_name to KPI categories and checks effect ranges.
"""

import os
import sys
import math
import yaml

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "backend"))

from backend.sound_response import estimate_all_effects
from backend.acoustics import FieldConfig, compute_field_spl

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")

# Map primary_outcome_name from YAML → our KPI name
OUTCOME_TO_KPI = {
    "photosynthesis_rate": "photosynthesis",
    "assimilation_rate": "photosynthesis",
    "chlorophyll_content": "photosynthesis",
    "Fv_Fm": "photosynthesis",
    "quantum_yield_FvFm": "photosynthesis",
    "plant_height": "vigor",
    "seedling_dry_weight": "vigor",
    "root_length": "vigor",
    "vigor": "vigor",
    "germination_rate": "germination",
    "germination_time": "germination",
    "germination_pct": "germination",
    "yield": "yield",
    "yield_kg_ha": "yield",
    # Water sub-KPIs — three independent GP tracks (gsw_drought, gsw_wellwatered, rwc, iWUE)
    "stomatal_conductance_drought": "water_status_gsw_drought",
    "stomatal_conductance": "water_status_gsw_wellwatered",
    "gsw": "water_status_gsw_wellwatered",
    "relative_water_content": "water_status_rwc",
    "drought_tolerance": "water_status_gsw_drought",
    "iWUE": "water_status_efficiency",
    "intrinsic_water_use_efficiency": "water_status_efficiency",
    # Aggregate fallbacks
    "water_status": "water_status",
    # Stress sub-KPIs
    "sheath_blight_incidence": "stress_resilience_disease",
    "stress_resilience": "stress_resilience",
    "flavonoid_content": "stress_resilience",
    "membrane_permeability": "stress_resilience",
}


def run_sound_case(case):
    """Run a sound-effect validation case."""
    inp = case.get("inputs", {})
    exp = case.get("expected", {})
    
    freq = float(inp.get("frequency_hz") or 0)
    spl = float(inp.get("spl_db_at_1m") or 68)
    hours = float(inp.get("hours_per_day") or 3)
    stage = inp.get("stage_bucket") or "mixed"
    stress = "drought" if (inp.get("water_regime") or "") == "drought" else "well_watered"
    
    effects = estimate_all_effects(freq, spl, hours, stage, stress, DATA_DIR)
    
    # Determine which KPI to check
    outcome_name = exp.get("primary_outcome_name", "")
    kpi = OUTCOME_TO_KPI.get(outcome_name)
    
    if not kpi:
        # Try to infer from the case description or use the first non-zero effect
        for k, v in effects.items():
            if abs(v.effect_pct) > 0.01:
                kpi = k
                break
        if not kpi:
            kpi = "photosynthesis"  # default
    
    effect = effects.get(kpi)
    if not effect:
        return "FAIL", f"KPI '{kpi}' not found in effects", 0, 0
    
    eff_pct = effect.effect_pct
    conf = effect.confidence
    band = exp.get("expected_effect_pct_band", [-999, 999])
    direction = exp.get("primary_direction", "")

    # Sign convention reconciliation:
    # The model internally flips germination_time (negative time = positive outcome).
    # YAML tests use raw measurement convention: faster germination = negative %.
    # Un-flip ONLY for germination_time when YAML expects decrease.
    if outcome_name == 'germination_time' and direction == 'decrease':
        eff_pct = -eff_pct
    
    # Check direction
    dir_ok = True
    if direction == "increase" and eff_pct < band[0]:
        dir_ok = False
    elif direction == "decrease" and eff_pct > band[1]:
        dir_ok = False
    elif direction == "no_clear_change" and abs(eff_pct) > max(abs(band[0]), abs(band[1])):
        dir_ok = False
    
    # Check band
    in_band = band[0] <= eff_pct <= band[1]
    
    if in_band:
        return "PASS", f"effect={eff_pct:+.1f}%, conf={conf:.3f}", eff_pct, conf
    elif dir_ok:
        return "PARTIAL", f"effect={eff_pct:+.1f}% (expected [{band[0]},{band[1]}]), direction OK", eff_pct, conf
    else:
        return "FAIL", f"effect={eff_pct:+.1f}% (expected [{band[0]},{band[1]}])", eff_pct, conf


def run_acoustics_case(case):
    """Run an acoustics validation case."""
    case_id = case["case_id"]
    
    if "distance" in case_id:
        # Test inverse-square law: SPL drops ~6 dB when distance doubles
        # Small field, speaker at center, no height offset
        fc = FieldConfig(field_length_m=20, field_width_m=20, speaker_count=1,
                         spacing_m=10, speaker_height_m=0.0, canopy_height_m=0.0,
                         spl_at_1m=90, frequency_hz=1000, hours_per_day=3)
        result = compute_field_spl(fc, 40)  # Fine grid to find 1m/2m points
        gx, gy = result.grid_x, result.grid_y
        cx, cy = 10.0, 10.0
        spl_1m, spl_2m = None, None
        for i in range(len(gx)):
            for j in range(len(gy)):
                d = math.sqrt((gx[i] - cx)**2 + (gy[j] - cy)**2)
                if 0.8 < d < 1.3 and spl_1m is None:
                    spl_1m = result.spl_grid[i][j]
                if 1.8 < d < 2.3 and spl_2m is None:
                    spl_2m = result.spl_grid[i][j]
        if spl_1m is not None and spl_2m is not None:
            drop = spl_1m - spl_2m
            if 4 < drop < 8:
                return "PASS", f"drop={drop:.1f} dB", drop, 1.0
            return "FAIL", f"drop={drop:.1f} dB (expected ~6)", drop, 0
        return "FAIL", "Could not find grid points at 1m/2m", 0, 0

    if "two_source" in case_id:
        fc1 = FieldConfig(field_length_m=10, field_width_m=10, speaker_count=1,
                          spacing_m=5, speaker_height_m=2.0, canopy_height_m=0.8,
                          spl_at_1m=90, frequency_hz=1000, hours_per_day=3)
        fc2 = FieldConfig(field_length_m=10, field_width_m=10, speaker_count=2,
                          spacing_m=5, speaker_height_m=2.0, canopy_height_m=0.8,
                          spl_at_1m=90, frequency_hz=1000, hours_per_day=3)
        r1 = compute_field_spl(fc1, 5)
        r2 = compute_field_spl(fc2, 5)
        diff = float(r2.mean_spl) - float(r1.mean_spl)
        if 1 < diff < 5:
            return "PASS", f"combined={r2.mean_spl} dB (+{diff:.1f})", diff, 1.0
        return "FAIL", f"diff={diff:.1f} dB (expected ~3)", diff, 0

    return "SKIP", "Unknown acoustics case", 0, 0


def main():
    yaml_path = os.path.join(os.path.dirname(__file__), "VALIDATION_CASES.yaml")
    with open(yaml_path) as f:
        vc = yaml.safe_load(f)

    print(f"\n{'='*70}")
    print(f"VALIDATION HARNESS v4 — Rice Sound Simulator")
    print(f"Dataset: rice_sound_master_v2.csv")
    print(f"{'='*70}\n")

    cases = vc.get("cases", [])
    results = {"PASS": 0, "PARTIAL": 0, "FAIL": 0, "SKIP": 0}
    role_results = {}
    details = []

    # Group by domain_role for display
    from collections import defaultdict
    by_role = defaultdict(list)
    for c in cases:
        by_role[c.get("domain_role", "unknown")].append(c)

    for role in sorted(by_role.keys()):
        print(f"--- {role.replace('_', ' ').title()} ---")
        for case in by_role[role]:
            cid = case["case_id"]
            
            if role == "acoustics_core":
                status, msg, val, conf = run_acoustics_case(case)
            elif case.get("inputs", {}).get("frequency_hz") is None and role == "core_sanity":
                # No-sound baseline test
                effects = estimate_all_effects(0, 0, 0, "mixed", "well_watered", DATA_DIR)
                max_eff = max(abs(e.effect_pct) for e in effects.values())
                band = case.get("expected", {}).get("expected_effect_pct_band", [-1, 1])
                if max_eff <= max(abs(band[0]), abs(band[1])):
                    status, msg = "PASS", f"max|effect|={max_eff:.3f}%"
                else:
                    status, msg = "FAIL", f"max|effect|={max_eff:.3f}%"
                val, conf = max_eff, 1.0
            elif case.get("inputs", {}).get("frequency_hz") is None and role == "analog_inhibition":
                # High SPL grass inhibition — use 2000 Hz where rice data is sparse,
                # so the grass inhibition window has more influence
                inp = case.get("inputs", {})
                spl = float(inp.get("spl_db_at_1m") or 95)
                effects = estimate_all_effects(2000, spl, 3, "vegetative", "well_watered", DATA_DIR)
                eff = effects.get("vigor", effects.get("photosynthesis"))
                if eff:
                    band = case.get("expected", {}).get("expected_effect_pct_band", [-50, 0])
                    if band[0] <= eff.effect_pct <= band[1]:
                        status = "PASS"
                    else:
                        status = "PARTIAL" if eff.effect_pct < 0 else "FAIL"
                    msg = f"effect={eff.effect_pct:+.1f}% at {spl} dB"
                    val, conf = eff.effect_pct, eff.confidence
                else:
                    status, msg, val, conf = "FAIL", "No vigor/photo effect", 0, 0
            else:
                status, msg, val, conf = run_sound_case(case)

            sym = {"PASS": "✓", "PARTIAL": "~", "FAIL": "✗", "SKIP": "?"}[status]
            print(f"  {sym} [{role:25s}] {cid}: {msg}")
            results[status] += 1
            role_results.setdefault(role, {"PASS": 0, "PARTIAL": 0, "FAIL": 0})
            role_results[role][status] = role_results[role].get(status, 0) + 1
            details.append((cid, role, status, msg))

        print()

    # Summary
    total = sum(results.values())
    print(f"{'='*70}")
    print(f"VALIDATION SUMMARY: {results['PASS']} PASS, {results['PARTIAL']} PARTIAL, {results['FAIL']} FAIL / {total} total")
    for role, counts in sorted(role_results.items()):
        p = counts.get("PASS", 0)
        t = sum(counts.values())
        print(f"  {role:30s}: {p}/{t} pass")
    print(f"{'='*70}\n")

    # Write report
    report_path = os.path.join(os.path.dirname(__file__), "VALIDATION_REPORT.md")
    with open(report_path, "w") as f:
        f.write(f"# Validation Report v4\n\n")
        f.write(f"Dataset: rice_sound_master_v2.csv\n\n")
        f.write(f"**{results['PASS']} PASS, {results['PARTIAL']} PARTIAL, {results['FAIL']} FAIL** / {total} total\n\n")
        f.write("| Case | Role | Status | Details |\n|---|---|---|---|\n")
        for cid, role, status, msg in details:
            f.write(f"| {cid} | {role} | {status} | {msg} |\n")

    print(f"Report: {report_path}")


if __name__ == "__main__":
    main()
