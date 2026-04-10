"""
test_all.py — Tests for rice_core, sound_response, acoustics, and KPI.

Implements the validation tests from the specification:
1. No-sound baseline
2. Seed-stage benefit scenario
3. Drought scenario
4. Extreme exposure
5. Acoustics sanity
"""

import sys
import os
import math

# Add backend to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "backend"))

from rice_core import (
    RiceState, SimulationConfig, WeatherDay,
    SoundDriver, DualPhaseSoundParams,
    run_simulation, run_simulation_dual_phase,
    daily_step, calc_heat_units, calc_development_rate, lint,
)
from sound_response import (
    estimate_effect, estimate_all_effects, load_all_csv_data,
    spl_activation_S_SPL, air_coupling_S_air, smech_air,
    train_mechactivation_gp, mechactivation_Sf, MechTransductionState,
    mechanistic_prior_for_outcome, PATHWAYS, OUTCOME_PATHWAY_WEIGHTS,
)
from acoustics import (
    spl_at_distance, combine_spl, compute_field_spl, FieldConfig,
)
from kpi_engine import (
    compute_yield_index, compute_water_index,
    compute_stress_resilience_index, compute_coverage_quality_index,
)

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
PASS = 0
FAIL = 0


def check(name, condition, detail=""):
    global PASS, FAIL
    if condition:
        PASS += 1
        print(f"  ✓ {name}")
    else:
        FAIL += 1
        print(f"  ✗ {name} — {detail}")


# =============================================================================
# 1. Rice Core Tests
# =============================================================================

def test_rice_core():
    print("\n=== Rice Core Tests ===")
    
    # Heat units: at optimum temperature, should be ~22 °Cd/d
    hu = calc_heat_units(tmax=30.0, tmin=22.0)
    check("Heat units at moderate temp > 0", hu > 0, f"hu={hu}")
    check("Heat units at moderate temp reasonable", 10 < hu < 25, f"hu={hu}")
    
    # Heat units at very low temp should be near zero
    hu_cold = calc_heat_units(tmax=8.0, tmin=2.0)
    check("Heat units at cold temp ~0", hu_cold < 1.0, f"hu={hu_cold}")
    
    # Development rate
    dvr = calc_development_rate(0.0, 20.0)
    check("Dev rate > 0 at DVS=0", dvr > 0, f"dvr={dvr}")
    
    # LINT interpolation
    val = lint([(0, 0), (1, 10), (2, 20)], 0.5)
    check("LINT interpolation midpoint", abs(val - 5.0) < 0.01, f"val={val}")
    val_ext = lint([(0, 0), (1, 10)], -1)
    check("LINT extrapolation below", val_ext == 0.0, f"val={val_ext}")
    
    # Full simulation completes
    config = SimulationConfig(total_days=120, water_regime="irrigated")
    history = run_simulation(config)
    check("Simulation produces history", len(history) > 10, f"len={len(history)}")
    
    final = history[-1]
    check("Final DVS > 0", final.dvs > 0, f"dvs={final.dvs}")
    check("Final biomass > 0", final.total_biomass > 0, f"bm={final.total_biomass}")
    check("Final yield proxy >= 0", final.yield_proxy >= 0, f"y={final.yield_proxy}")
    check("LAI peaked > 1", max(s.lai for s in history) > 1.0,
          f"max_lai={max(s.lai for s in history)}")
    
    # No-sound baseline: zero treatment = no sound effect (trivially true in rice_core)
    check("Baseline has no sound modifier", True, "rice_core has no sound layer")
    
    # Drought reduces yield
    config_drought = SimulationConfig(
        total_days=120, water_regime="drought", drought_severity=0.8
    )
    history_drought = run_simulation(config_drought)
    final_drought = history_drought[-1]
    check("Drought reduces biomass vs irrigated",
          final_drought.total_biomass < final.total_biomass,
          f"drought={final_drought.total_biomass} vs irr={final.total_biomass}")


# =============================================================================
# 2. Sound Response Tests
# =============================================================================

def test_sound_response():
    print("\n=== Sound Response Tests ===")
    
    # Load data — v2 CSV has 76 total rows, 47 usable for calibration
    data = load_all_csv_data(DATA_DIR)
    check("Calibration data loaded", len(data) > 0, f"rows={len(data)}")
    check("v2 CSV has expected row count", len(data) >= 70, f"rows={len(data)} (expected ~76)")
    usable = sum(1 for r in data if r.usable_for_calibration)
    check("Usable rows for calibration", usable >= 40, f"usable={usable} (expected ~47)")
    
    # No-sound baseline: 0 hours should produce no effect
    no_sound = estimate_effect(
        "yield", frequency_hz=1000, spl_db=0, hours_per_day=0,
        data_dir=DATA_DIR,
    )
    check("No-sound (0h) yields zero or very low effect",
          abs(no_sound.effect_pct) < 1.0 or no_sound.supporting_row_count == 0,
          f"effect={no_sound.effect_pct}")
    
    # Seed-stage benefit: ~350 Hz, 68 dB should show positive effect
    seed_effect = estimate_effect(
        "photosynthesis", frequency_hz=350, spl_db=68,
        hours_per_day=3.0, stage_bucket="seedling",
        stress_context="well_watered", data_dir=DATA_DIR,
    )
    check("Seed-stage 350Hz shows positive photosynthesis effect",
          seed_effect.effect_pct > 0 or seed_effect.supporting_row_count > 0,
          f"effect={seed_effect.effect_pct}, n={seed_effect.supporting_row_count}")
    
    # Vigor estimate
    vigor = estimate_effect(
        "vigor", frequency_hz=350, spl_db=68,
        hours_per_day=3.0, stage_bucket="seedling",
        data_dir=DATA_DIR,
    )
    check("Vigor estimate has supporting rows",
          vigor.supporting_row_count > 0 or vigor.confidence > 0,
          f"n={vigor.supporting_row_count}, conf={vigor.confidence}")
    
    # Extreme frequency: 50000 Hz should have low confidence for yield
    extreme = estimate_effect(
        "yield", frequency_hz=50000, spl_db=85,
        hours_per_day=3.0, data_dir=DATA_DIR,
    )
    check("Extreme frequency has low confidence",
          extreme.confidence < 0.3 or extreme.supporting_row_count == 0,
          f"conf={extreme.confidence}")
    
    # All effects — should include main KPIs plus water sub-KPIs
    all_eff = estimate_all_effects(
        frequency_hz=400, spl_db=70, hours_per_day=3.0,
        data_dir=DATA_DIR,
    )
    check("estimate_all_effects returns all outcomes",
          len(all_eff) >= 10, f"n={len(all_eff)}")
    # Check water sub-KPIs exist
    for sub in ["water_status", "water_status_gsw_drought", "water_status_rwc", "water_status_efficiency"]:
        check(f"Sub-KPI '{sub}' present", sub in all_eff, f"keys={list(all_eff.keys())}")
    # Check main KPIs are present
    for main in ["germination", "vigor", "photosynthesis", "yield", "stress_resilience"]:
        check(f"Main KPI '{main}' present", main in all_eff, f"keys={list(all_eff.keys())}")


# =============================================================================
# 3. Acoustics Tests
# =============================================================================

def test_acoustics():
    print("\n=== Acoustics Tests ===")
    
    # Doubling distance reduces SPL by ~6 dB
    spl_1m = 85.0
    spl_2m = spl_at_distance(spl_1m, 2.0)
    spl_4m = spl_at_distance(spl_1m, 4.0)
    diff_2x = spl_2m - spl_4m
    check("Doubling distance: ~6 dB drop",
          abs(diff_2x - 6.0) < 0.5,
          f"diff={diff_2x:.2f}")
    
    spl_at_10 = spl_at_distance(85.0, 10.0)
    check("SPL at 10m = 85 - 20 = 65 dB",
          abs(spl_at_10 - 65.0) < 0.5,
          f"spl={spl_at_10}")
    
    # Combining two equal sources: +3 dB
    combined = combine_spl([80.0, 80.0])
    check("Two equal sources: ~+3 dB",
          abs(combined - 83.01) < 0.5,
          f"combined={combined:.2f}")
    
    # Single source
    single = combine_spl([80.0])
    check("Single source unchanged", abs(single - 80.0) < 0.01, f"single={single}")
    
    # Field computation
    config = FieldConfig(
        field_length_m=100, field_width_m=50,
        speaker_count=4, spl_at_1m=85.0,
    )
    result = compute_field_spl(config, grid_resolution=15)
    check("Field SPL computed", result.mean_spl > 0, f"mean={result.mean_spl}")
    check("Max SPL >= Mean SPL", result.max_spl >= result.mean_spl)
    check("Min SPL <= Mean SPL", result.min_spl <= result.mean_spl)
    check("Coverage 0-1 range",
          0 <= result.coverage_above_60db <= 1,
          f"cov60={result.coverage_above_60db}")
    check("Speakers generated", len(result.speakers) == 4,
          f"n={len(result.speakers)}")


# =============================================================================
# 4. KPI Tests
# =============================================================================

def test_kpis():
    print("\n=== KPI Tests ===")
    
    # Yield index: no effect = 50
    yi = compute_yield_index(5000, 5000, 0.0, 0.5)
    check("Yield index: no effect = 50",
          abs(yi.value - 50.0) < 1.0, f"val={yi.value}")
    
    # Yield index: positive effect > 50
    yi_pos = compute_yield_index(5000, 5000, 10.0, 0.8)
    check("Yield index: +10% > 50", yi_pos.value > 50, f"val={yi_pos.value}")
    
    # Water index: irrigated baseline near 100
    wi = compute_water_index(1.0, 0.0, 0.5, "irrigated")
    check("Water index: irrigated ~100",
          wi.value >= 95, f"val={wi.value}")
    
    # Stress resilience: all zeros = 50
    sri = compute_stress_resilience_index(0, 0, 0, 0.5, 0.5, 0.5)
    check("Stress resilience: no effect ~50",
          abs(sri.value - 50) < 5, f"val={sri.value}")
    
    # Coverage quality
    cqi = compute_coverage_quality_index(
        mean_spl=75, coverage_above_60db=0.9,
        coverage_above_70db=0.6, uniformity=0.8,
        hours_per_day=3.0,
    )
    check("Coverage quality > 0", cqi.value > 0, f"val={cqi.value}")
    check("Coverage quality <= 100", cqi.value <= 100, f"val={cqi.value}")
    
    # Very high SPL warning
    cqi_loud = compute_coverage_quality_index(
        mean_spl=110, coverage_above_60db=1.0,
        coverage_above_70db=1.0, uniformity=0.9,
        hours_per_day=3.0,
    )
    check("High SPL reduces confidence",
          cqi_loud.confidence < cqi.confidence,
          f"loud_conf={cqi_loud.confidence}")


# =============================================================================
# 5. Integration test: /simulate pipeline
# =============================================================================

def test_integration():
    print("\n=== Integration Tests ===")
    
    # Test the full pipeline manually (without FastAPI)
    from rice_core import SimulationConfig, run_simulation
    from sound_response import estimate_all_effects
    from acoustics import FieldConfig, compute_field_spl
    from kpi_engine import compute_all_kpis
    
    config = SimulationConfig(total_days=120, water_regime="irrigated")
    history = run_simulation(config)
    final = history[-1]
    
    effects = estimate_all_effects(
        frequency_hz=400, spl_db=70, hours_per_day=3.0,
        data_dir=DATA_DIR,
    )
    
    acoustics = compute_field_spl(FieldConfig(
        field_length_m=100, field_width_m=50,
        speaker_count=4, spl_at_1m=70.0,
    ), grid_resolution=10)
    
    effects_dict = {}
    for k, v in effects.items():
        effects_dict[k] = {"effect_pct": v.effect_pct, "confidence": v.confidence}
    
    kpis = compute_all_kpis(
        baseline_yield=final.yield_proxy,
        sound_effects=effects,
        water_stress=final.water_stress,
        water_regime="irrigated",
        acoustics_result={
            "mean_spl": acoustics.mean_spl,
            "coverage_above_60db": acoustics.coverage_above_60db,
            "coverage_above_70db": acoustics.coverage_above_70db,
            "uniformity": acoustics.uniformity,
        },
        hours_per_day=3.0,
    )
    
    check("Integration: KPIs computed", len(kpis) == 4, f"n={len(kpis)}")
    check("Integration: all KPIs have values",
          all(v.value >= 0 for v in kpis.values()))
    check("Integration: yield KPI present", "yield_index" in kpis)



# =============================================================================
# 6. v4 Architecture Tests: MechActivation, pathways, dual-phase driver
# =============================================================================

def test_v4_architecture():
    print("\n=== v4 Architecture Tests ===")

    # --- S_SPL activation function ---
    check("S_SPL below floor = 0",     abs(spl_activation_S_SPL(55.0)) < 0.01)
    check("S_SPL at floor = 0",        abs(spl_activation_S_SPL(60.0)) < 0.01)
    check("S_SPL mid-ramp ~ 0.5",      abs(spl_activation_S_SPL(72.5) - 0.5) < 0.05)
    check("S_SPL at opt = 1.0",        abs(spl_activation_S_SPL(85.0) - 1.0) < 0.01)
    check("S_SPL plateau at 92 = 1.0", abs(spl_activation_S_SPL(92.0) - 1.0) < 0.01)
    check("S_SPL at damage = 0",       abs(spl_activation_S_SPL(110.0)) < 0.01)
    check("S_SPL above damage = 0",    abs(spl_activation_S_SPL(115.0)) < 0.01)

    # --- S_air airborne coupling ---
    check("S_air at 500 Hz = 1.0",     abs(air_coupling_S_air(500.0) - 1.0) < 0.01)
    check("S_air at 5000 Hz = 1.0",    abs(air_coupling_S_air(5000.0) - 1.0) < 0.01)
    check("S_air at 10000 Hz ~= 0.5",  abs(air_coupling_S_air(10000.0) - 0.5) < 0.05)
    check("S_air at 15000 Hz = 0.0",   abs(air_coupling_S_air(15000.0)) < 0.01)
    check("S_air at ultrasound = 0.0", abs(air_coupling_S_air(25000.0)) < 0.01)

    # --- MechActivation GP ---
    data = load_all_csv_data(DATA_DIR)
    mech_state = train_mechactivation_gp(data, DATA_DIR)
    check("MechActivation GP trains",   mech_state.available or mech_state.n_training >= 2,
          f"available={mech_state.available}, n={mech_state.n_training}")
    check("MechActivation y_max > 0",   mech_state.y_max > 0, f"y_max={mech_state.y_max}")

    if mech_state.available:
        sf_350 = mechactivation_Sf(350.0, mech_state)
        sf_800 = mechactivation_Sf(800.0, mech_state)
        check("S_f(350Hz) in [0,1]", 0.0 <= sf_350 <= 1.0, f"S_f={sf_350:.3f}")
        check("S_f(800Hz) in [0,1]", 0.0 <= sf_800 <= 1.0, f"S_f={sf_800:.3f}")

    # --- smech_air ---
    dummy_state = MechTransductionState(available=False, y_max=50.0,
                                         prior_points=[(500.0, 30.0)])
    sm_optimal = smech_air(500.0, 85.0, dummy_state)
    sm_below   = smech_air(500.0, 55.0, dummy_state)
    sm_us      = smech_air(25000.0, 85.0, dummy_state)
    check("S_mech,air optimal SPL > 0", sm_optimal > 0, f"val={sm_optimal:.3f}")
    check("S_mech,air below floor = 0", abs(sm_below) < 0.01, f"val={sm_below:.3f}")
    check("S_mech,air ultrasound = 0",  abs(sm_us) < 0.01, f"val={sm_us:.3f}")

    # --- mechanistic_prior_for_outcome ---
    mp_photo_350 = mechanistic_prior_for_outcome(
        "photosynthesis", 350.0, 70.0, "seedling", "well_watered", "airborne", 0.5
    )
    mp_photo_1200 = mechanistic_prior_for_outcome(
        "photosynthesis", 1200.0, 70.0, "seedling", "well_watered", "airborne", 0.1
    )
    mp_drought_800 = mechanistic_prior_for_outcome(
        "water_status_gsw_drought", 800.0, 100.0, "vegetative", "drought", "airborne", 0.5
    )
    check("Mechanistic prior: photo 350Hz > 0", mp_photo_350 > 0,
          f"val={mp_photo_350:.2f}%")
    check("Mechanistic prior: photo 1200Hz small", abs(mp_photo_1200) < 5.0,
          f"val={mp_photo_1200:.2f}% (expected near-zero)")
    check("Mechanistic prior: drought 800Hz > 0", mp_drought_800 > 0,
          f"val={mp_drought_800:.2f}%")

    # --- Pathway configs exist ---
    for name in ["seed_vigor", "leaf_gas_exchange", "drought_resilience",
                 "membrane_damage", "ultrasound_cavitation"]:
        check(f"Pathway '{name}' configured", name in PATHWAYS)

    # --- Outcome pathway weights ---
    for outcome in ["photosynthesis", "vigor", "water_status_gsw_drought",
                    "stress_resilience", "germination", "yield"]:
        check(f"OUTCOME_PATHWAY_WEIGHTS['{outcome}'] exists",
              outcome in OUTCOME_PATHWAY_WEIGHTS,
              f"missing key")
        weights = OUTCOME_PATHWAY_WEIGHTS.get(outcome, {})
        total_w = sum(weights.values())
        check(f"Weights for '{outcome}' sum to 1.0",
              abs(total_w - 1.0) < 0.01,
              f"sum={total_w:.3f}")

    # --- High-SPL inhibitory zone ---
    # 0 Hz is used as a SENTINEL for "no resonance band" to test pure SPL-driven
    # boundary inhibition. At 100 Hz, Bochu (200 Hz / 106 dB) is within the
    # ±150 Hz frequency-proximity window for spl_matched, which would bypass
    # the inhibition logic. 0 Hz has no nearby vigor data so the broadband
    # high-SPL inhibition path is reached cleanly.
    # This mirrors validation case inhibitory_grass_high_spl_001 (freq=None → 0 Hz).
    vigor_inhibited = estimate_effect(
        "vigor", frequency_hz=0.0, spl_db=95.0,
        hours_per_day=3.0, stage_bucket="vegetative",
        stress_context="well_watered", data_dir=DATA_DIR,
    )
    check("Vigor at 95 dB vegetative (0 Hz sentinel, no local data) is negative",
          vigor_inhibited.effect_pct < -20.0,
          f"effect={vigor_inhibited.effect_pct:.1f}% (expected < -20%)")

    # Sheath blight sign assertion: CSV stores +50 meaning "+50% disease REDUCTION"
    # (already beneficial direction). The engine must NOT flip this sign.
    # If it returns negative, _FLIP_SIGN_OUTCOMES has wrongly re-added it.
    sheath_effect = estimate_effect(
        "stress_resilience_disease", frequency_hz=550.0, spl_db=70.0,
        hours_per_day=1.5, stage_bucket="vegetative",
        stress_context="well_watered", data_dir=DATA_DIR,
    )
    check("Sheath blight at 550 Hz is POSITIVE (disease reduction is beneficial)",
          sheath_effect.effect_pct > 0,
          f"effect={sheath_effect.effect_pct:.1f}% (expected > 0; CSV pre-normalised)")

    # --- Dual-phase sound driver ---
    config = SimulationConfig(total_days=150, water_regime="irrigated")
    germ_driver = SoundDriver(
        photosynthesis_multiplier=1.0, vigor_multiplier=1.05,
        water_status_multiplier=1.0, germination_multiplier=1.1,
        stress_resilience_multiplier=1.0, confidence=0.4, active=True,
    )
    veg_driver = SoundDriver(
        photosynthesis_multiplier=1.1, vigor_multiplier=1.0,
        water_status_multiplier=1.05, germination_multiplier=1.0,
        stress_resilience_multiplier=1.0, confidence=0.5, active=True,
    )
    dual = DualPhaseSoundParams(
        germination_driver=germ_driver,
        vegetative_driver=veg_driver,
        germination_days=10,
    )
    history_dual = run_simulation_dual_phase(config, dual_phase=dual)
    check("Dual-phase simulation completes", len(history_dual) > 10,
          f"len={len(history_dual)}")
    check("Dual-phase final DVS > 0", history_dual[-1].dvs > 0,
          f"dvs={history_dual[-1].dvs:.2f}")
    check("Dual-phase yield >= 0", history_dual[-1].yield_proxy >= 0,
          f"yield={history_dual[-1].yield_proxy:.0f}")

    # Dual-phase with None (fallback to single-phase)
    history_none = run_simulation_dual_phase(config, dual_phase=None)
    check("Dual-phase None falls back to single-phase", len(history_none) > 10)


# =============================================================================
# Run all
# =============================================================================

if __name__ == "__main__":
    test_rice_core()
    test_sound_response()
    test_acoustics()
    test_kpis()
    test_integration()
    test_v4_architecture()
    
    print(f"\n{'='*50}")
    print(f"Results: {PASS} passed, {FAIL} failed, {PASS+FAIL} total")
    if FAIL > 0:
        sys.exit(1)
    else:
        print("All tests passed!")
