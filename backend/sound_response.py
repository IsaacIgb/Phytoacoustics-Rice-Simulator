"""
sound_response.py — Empirical-window sound-response engine (v2).

Architecture:
  This engine does NOT attempt to fit a continuous dose-response curve
  or build a mechanistic biophysical model. The phytoacoustics literature
  lacks the quantitative coupling constants needed for either approach.

  Instead, it uses the "empirical windowing" approach documented in the
  PAFT (Plant Acoustic Frequency Technology) literature:
    1. BIOPHYSICAL BOUNDARIES — hard physics constraints from cross-species
       membrane mechanosensory research (activation floor, damage ceiling).
    2. EMPIRICAL WINDOWS — discrete frequency×SPL×stage bands with
       associated effect directions and magnitude ranges, derived from
       published rice and cereal studies.
    3. CSV EVIDENCE LOOKUP — granular row-level data from attached CSVs
       for refinement within windows.
    4. CROSS-SPECIES TIERS — monocot/grass data used at reduced confidence
       for directional support only, never as primary calibration.

  The result is a literature-guided hypothesis generator, not a
  predictive model. It can say "the literature suggests X in this
  frequency range" but cannot predict from first principles.

Sources:
  - Bochu et al. 2003 (rice seed, 400 Hz / 4000 Hz)
  - Jusoh & Ramlee 2023 (rice seedling, 350-380 Hz)
  - Jeong et al. 2008 (rice drought, 50-250 Hz gene expression, 800 Hz RWC)
  - Hassan et al. 2014 (rice drought, 800-1500 Hz)
  - Qi et al. 2010 (rice chlorophyll, 550 Hz)
  - Hou et al. 2009 (rice yield, PAFT 550 Hz)
  - Munasinghe et al. 2023 (rice stomatal, 350 Hz)
  - Sri Lankan rice 2024 (3000-5000 Hz vegetative)
  - Cross-species: wheat, corn, oat, barley, sorghum, general grasses
  See MODEL_NOTES.md and VALIDATION_CASES.yaml for full source list.
"""

import csv
import math
import os
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Tuple

# =============================================================================
# Layer 1: BIOPHYSICAL BOUNDARIES
# Hard constraints from conserved plant mechanosensory physics.
# These apply regardless of species because they are about membrane
# mechanics, not species-specific resonance.
# =============================================================================

# Minimum SPL to trigger measurable mechanosensory response.
# Source: cross-species consensus, 60 dB observed as minimum threshold
# for measurable vibrations in plant biological structures.
ACTIVATION_FLOOR_DB = 60.0

# SPL above which cellular membrane damage begins.
# Source: cross-species consensus, >110 dB damages cell membranes.
DAMAGE_CEILING_DB = 110.0

# SPL above which continuous broadband exposure stunts growth in grasses.
# Source: Grass_General_90dB, >90 dB decreases vegetative growth by >40%.
GRASS_INHIBITION_DB = 90.0

# Frequency band where most biological activity occurs across species.
# Source: cross-species consensus, 100-2000 Hz most active.
ACTIVE_BAND_LOW_HZ = 100.0
ACTIVE_BAND_HIGH_HZ = 2000.0

# Frequency above which confidence decays for airborne sound effects.
# Sri Lankan rice shows effects at 3-5 kHz, wheat at 5 kHz.
# Wheat 12 kHz shows no effect. Rather than a hard ceiling,
# confidence decays linearly from 5 kHz to 15 kHz.
AIRBORNE_CONFIDENCE_DECAY_START_HZ = 5000.0
AIRBORNE_CONFIDENCE_DECAY_END_HZ = 15000.0

# Ultrasound boundary — different mechanism (liquid cavitation)
ULTRASOUND_FLOOR_HZ = 20000.0


# =============================================================================
# Layer 2: EMPIRICAL WINDOWS
# Discrete frequency×SPL×stage bands from published studies.
# Each window has: frequency range, SPL range, stage, crop tier,
# expected direction, magnitude band, and confidence.
# =============================================================================

@dataclass
class EmpiricalWindow:
    """A discrete empirical evidence window."""
    window_id: str
    freq_low_hz: float
    freq_high_hz: float
    spl_low_db: Optional[float]  # None = any SPL above activation floor
    spl_high_db: Optional[float]
    stages: List[str]            # e.g. ["seed", "seedling"]
    crop_tier: str               # "rice", "monocot_grass", "cross_species"
    direction: str               # "increase", "decrease", "none"
    effect_pct_low: float
    effect_pct_high: float
    confidence: float            # 0-1
    outcomes: List[str]          # which outcomes this window applies to
    source: str                  # citation
    notes: str = ""


# Rice-specific windows (highest priority)
RICE_WINDOWS = [
    EmpiricalWindow(
        window_id="rice_photo_350_400",
        freq_low_hz=300, freq_high_hz=450,
        spl_low_db=60, spl_high_db=85,
        stages=["seedling", "vegetative"],
        crop_tier="rice", direction="increase",
        effect_pct_low=15, effect_pct_high=45,
        confidence=0.75,
        outcomes=["photosynthesis", "vigor"],
        source="Jusoh 2023 (350-380 Hz assimilation/height)",
    ),
    EmpiricalWindow(
        window_id="rice_wue_350_weak",
        freq_low_hz=300, freq_high_hz=400,
        spl_low_db=None, spl_high_db=None,
        stages=["vegetative"],
        crop_tier="rice", direction="increase",
        effect_pct_low=0, effect_pct_high=15,
        confidence=0.20,
        outcomes=["water_status"],
        source="Munasinghe 2023 (350 Hz iWUE, no SPL/no effect size)",
        notes="Weak evidence: no SPL, no numeric effect size reported",
    ),
    EmpiricalWindow(
        window_id="rice_yield_paft_400_700",
        freq_low_hz=400, freq_high_hz=700,
        spl_low_db=65, spl_high_db=90,
        stages=["mixed", "vegetative", "reproductive"],
        crop_tier="rice", direction="increase",
        effect_pct_low=3, effect_pct_high=15,
        confidence=0.55,
        outcomes=["yield"],
        source="Hou 2009 (550 Hz PAFT yield +5.7%); Qi 2010 (chlorophyll +10%)",
    ),
    EmpiricalWindow(
        window_id="rice_chlorophyll_500_600",
        freq_low_hz=450, freq_high_hz=650,
        spl_low_db=70, spl_high_db=90,
        stages=["vegetative"],
        crop_tier="rice", direction="increase",
        effect_pct_low=5, effect_pct_high=15,
        confidence=0.50,
        outcomes=["photosynthesis"],
        source="Qi 2010 (550 Hz chlorophyll +10%)",
    ),
    EmpiricalWindow(
        window_id="rice_drought_800_1500",
        freq_low_hz=800, freq_high_hz=1500,
        spl_low_db=60, spl_high_db=85,
        stages=["vegetative"],
        crop_tier="rice", direction="increase",
        effect_pct_low=0, effect_pct_high=20,
        confidence=0.35,
        outcomes=["water_status", "stress_resilience"],
        source="Jeong 2008 (800 Hz RWC, qualitative); Hassan 2014 (800-1500 Hz drought)",
        notes="Qualitative evidence only; magnitude uncertain",
    ),
    EmpiricalWindow(
        window_id="rice_germ_400_seed",
        freq_low_hz=350, freq_high_hz=450,
        spl_low_db=95, spl_high_db=110,
        stages=["seed", "germination"],
        crop_tier="rice", direction="increase",
        effect_pct_low=5, effect_pct_high=25,
        confidence=0.40,
        outcomes=["germination"],
        source="Bochu 2003 (400 Hz 106 dB germination)",
        notes="High SPL needed for seed penetration",
    ),
    EmpiricalWindow(
        window_id="rice_gene_expression_50_250",
        freq_low_hz=50, freq_high_hz=300,
        spl_low_db=60, spl_high_db=80,
        stages=["seedling", "vegetative"],
        crop_tier="rice", direction="increase",
        effect_pct_low=0, effect_pct_high=10,
        confidence=0.25,
        outcomes=["stress_resilience"],
        source="Jeong 2008 (50 Hz gene expression; 250 Hz drought prep)",
        notes="Gene-level only, no direct phenotypic magnitude",
    ),
    EmpiricalWindow(
        window_id="rice_sri_lankan_3k_5k",
        freq_low_hz=3000, freq_high_hz=5000,
        spl_low_db=None, spl_high_db=None,
        stages=["vegetative"],
        crop_tier="rice", direction="increase",
        effect_pct_low=0, effect_pct_high=15,
        confidence=0.20,
        outcomes=["vigor"],
        source="Sri Lankan rice 2024 (3-5 kHz vegetative volume)",
        notes="Sri Lankan varieties; may not generalize to all rice",
    ),
    # INHIBITORY windows
    EmpiricalWindow(
        window_id="rice_injury_4000_high_spl",
        freq_low_hz=3500, freq_high_hz=5000,
        spl_low_db=105, spl_high_db=130,
        stages=["seed", "germination", "seedling", "vegetative"],
        crop_tier="rice", direction="decrease",
        effect_pct_low=-30, effect_pct_high=-5,
        confidence=0.50,
        outcomes=["germination", "vigor", "photosynthesis"],
        source="Bochu 2003 (4000 Hz 111 dB injury)",
    ),
]

# Cross-species grass/monocot windows (lower confidence tier)
CROSS_SPECIES_WINDOWS = [
    EmpiricalWindow(
        window_id="grass_high_spl_inhibition",
        freq_low_hz=0, freq_high_hz=20000,
        spl_low_db=90, spl_high_db=130,
        stages=["vegetative", "seedling", "reproductive"],
        crop_tier="cross_species", direction="decrease",
        effect_pct_low=-50, effect_pct_high=-10,
        confidence=0.80,  # High confidence — well-documented biophysical boundary
        outcomes=["vigor", "photosynthesis", "yield"],
        source="Grass_General_90dB (>40% growth reduction)",
        notes="Continuous exposure >90 dB stunts grasses — conserved membrane physics",
    ),
    EmpiricalWindow(
        window_id="grass_1000hz_defense_shift",
        freq_low_hz=800, freq_high_hz=1500,
        spl_low_db=None, spl_high_db=None,
        stages=["vegetative"],
        crop_tier="cross_species", direction="increase",
        effect_pct_low=0, effect_pct_high=10,
        confidence=0.20,
        outcomes=["stress_resilience"],
        source="Grass_General_1000Hz (defense metabolite shift)",
        notes="Growth-to-defense transition boundary in grasses",
    ),
]

ALL_WINDOWS = RICE_WINDOWS + CROSS_SPECIES_WINDOWS

# Confidence discount for cross-species extrapolation
# Confidence discount for cross-species empirical WINDOWS only.
# CSV data uses concordance-based scoring instead (see _concordance_discount).
# Windows represent structural biophysical boundaries (>90 dB inhibition, etc.)
# which are about conserved membrane physics, so a moderate flat discount applies.
CROSS_SPECIES_WINDOW_DISCOUNT = 0.6


# =============================================================================
# Layer 3: CSV EVIDENCE LOOKUP
# Granular row-level data for refinement within windows.
# =============================================================================

@dataclass
class SoundStudyRow:
    study_id: str
    frequency_hz: float
    spl_db: float
    spl_db_known: bool
    hours_per_day: float
    total_days: float
    sound_type: str
    outcome_name: str
    reported_effect_size: float
    calculated_effect_size_pct: float
    effect_direction: str
    crop_type: str           # from crop_tier column
    stage_bucket: str
    stress_context: str
    confidence_weight: str
    completeness_class: str
    usable_for_calibration: bool
    kpi_category: str        # VIGOR, PHOTOSYNTHESIS, WATER_STATUS, etc.
    notes: str


def _safe_float(val, default=0.0):
    try:
        return float(val)
    except (ValueError, TypeError):
        return default


def _infer_stress(water_regime):
    wr = (water_regime or '').lower()
    if 'drought' in wr:
        return 'drought'
    return 'well_watered'


def load_all_csv_data(data_dir=None):
    """Load rice_sound_master_v2.csv (canonical dataset)."""
    if data_dir is None:
        data_dir = os.path.join(os.path.dirname(__file__), "..", "data")

    rows = []
    path = os.path.join(data_dir, "rice_sound_master_v2.csv")
    if not os.path.exists(path):
        return rows

    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            stage = r.get("growth_stage", "mixed")
            stress = _infer_stress(r.get("water_regime", ""))

            eff_size = _safe_float(r.get("reported_effect_size", ""))
            calc_size = _safe_float(r.get("calculated_effect_size_pct", ""))
            direction = ""
            if eff_size > 0 or calc_size > 0:
                direction = "increase"
            elif eff_size < 0 or calc_size < 0:
                direction = "decrease"
            else:
                direction = "unknown"

            spl_known = r.get("spl_db_known", "no").lower() == "yes"

            row = SoundStudyRow(
                study_id=r.get("study_id", ""),
                frequency_hz=_safe_float(r.get("frequency_hz", "")),
                spl_db=_safe_float(r.get("spl_db", "")),
                spl_db_known=spl_known,
                hours_per_day=_safe_float(r.get("hours_per_day", "")),
                total_days=_safe_float(r.get("total_days", "")),
                sound_type=r.get("sound_type", ""),
                outcome_name=r.get("outcome_name", ""),
                reported_effect_size=eff_size,
                calculated_effect_size_pct=calc_size,
                effect_direction=direction,
                crop_type=r.get("crop_tier", "rice"),
                stage_bucket=stage,
                stress_context=stress,
                confidence_weight=r.get("confidence_weight", "low"),
                completeness_class=r.get("completeness_class", "insufficient"),
                usable_for_calibration=r.get("usable_for_calibration", "no").lower() == "yes",
                kpi_category=r.get("kpi_category", ""),
                notes=r.get("notes", ""),
            )
            rows.append(row)

    return rows


_ALL_DATA: Optional[List[SoundStudyRow]] = None

def get_all_data(data_dir=None):
    global _ALL_DATA
    if _ALL_DATA is None:
        _ALL_DATA = load_all_csv_data(data_dir)
    return _ALL_DATA


# =============================================================================
# KPI category mapping
# The v2 CSV has an explicit kpi_category column. We use it directly
# rather than mapping outcome_name → KPI. This is cleaner and avoids
# the old problem of unmapped outcome names being silently dropped.
#
# For backward compatibility, OUTCOME_MAPPING is kept but now only used
# as fallback when kpi_category is empty.
# =============================================================================

KPI_CATEGORIES = {
    "germination":                  ["GERMINATION", "VIGOR"],
    "vigor":                        ["VIGOR"],
    "photosynthesis":               ["PHOTOSYNTHESIS"],
    "water_status":                 ["WATER_STATUS", "DROUGHT_TOLERANCE"],
    "water_status_openness":        ["WATER_STATUS", "DROUGHT_TOLERANCE"],
    "water_status_gsw_drought":     ["DROUGHT_TOLERANCE"],
    "water_status_gsw_wellwatered": ["WATER_STATUS"],
    "water_status_rwc":             ["DROUGHT_TOLERANCE"],
    "water_status_efficiency":      ["WATER_STATUS"],
    "stress_resilience":            ["STRESS_RESILIENCE"],
    "stress_resilience_disease":    ["STRESS_RESILIENCE"],
    "yield":                        ["YIELD"],
}

# Sub-outcome filters control which CSV rows feed each sub-KPI's Lorentzian.
#
# Key separation for water status (follows CWSI theory, Idso et al. 1981):
#   gsw_drought  = active stomatal regulation under water deficit (most sensitive)
#   gsw_wellwatered = stomatal response under normal irrigation (different physiology)
#   RWC = passive hydration state (slower, smaller magnitude response)
#
# These MUST NOT be pooled into one Lorentzian because:
#   1. gsw responds 2.4× more than RWC (Jeong 2014: +89% vs +37% at 1500 Hz)
#   2. gsw and RWC operate on different timescales (minutes vs hours-days)
#   3. Drought gsw and well-watered gsw have different frequency peaks
#      (Jeong peaks at ~1000-1500 Hz, Jusoh oscillates at 350-380 Hz)
#
# Reference: Nature meta-analysis (2016, 164 studies): "the gs response ratio
# was greater than those of A, Tr, RWC, and LWP under different drought intensities"
SUB_OUTCOME_FILTERS = {
    "water_status_gsw_drought": {
        "include": ["stomatal_conductance_drought"],
        "exclude": [],
    },
    "water_status_gsw_wellwatered": {
        "include": ["stomatal_conductance"],
        "exclude": ["stomatal_conductance_drought"],
    },
    "water_status_rwc": {
        "include": ["relative_water_content"],
        "exclude": [],
    },
    # Legacy aggregate — still used if queried directly
    "water_status_openness": {
        "include": ["stomatal_conductance", "stomatal_conductance_drought",
                    "gsw", "relative_water_content"],
        "exclude": ["intrinsic_water_use_efficiency", "iWUE"],
    },
    "water_status_efficiency": {
        "include": ["intrinsic_water_use_efficiency", "iWUE"],
        "exclude": [],
    },
    "stress_resilience_disease": {
        "include": ["sheath_blight_incidence"],
        "exclude": [],
    },
    "vigor": {
        "include": [],
        "exclude": ["fresh_weight_increase_rate"],
    },
}

# Legacy fallback only
OUTCOME_MAPPING = {
    "germination": ["germination_pct", "germination_index", "germination_time",
                    "germination_time_h", "speed_of_germination", "water_uptake_rate",
                    "germination_potential", "germination_energy", "germination_speed"],
    "vigor": ["seedling_dry_weight", "seedling_dry_weight_g", "plant_height",
              "plant_height_stem", "root_length", "root_total_length",
              "fresh_weight", "fresh_weight_increase_rate",
              "vegetative_volume", "total_dry_weight", "shoot_dry_weight",
              "early_vegetative_growth", "vegetative_growth_rate"],
    "photosynthesis": ["photosynthesis_rate", "assimilation_rate", "Asat",
                       "chlorophyll_content", "quantum_yield_FvFm"],
    "water_status": ["relative_water_content", "intrinsic_water_use_efficiency",
                     "iWUE", "stomatal_conductance", "stomatal_conductance_drought",
                     "gsw", "stomatal_conductance_gsw"],
    "water_status_openness": ["relative_water_content", "stomatal_conductance",
                              "stomatal_conductance_drought", "gsw"],
    "water_status_efficiency": ["intrinsic_water_use_efficiency", "iWUE"],
    "stress_resilience": ["ROS_scavenging_enzymes", "flavonoid_content",
                          "total_flavonoid_content", "SOD_CAT_POD_enzyme_activity",
                          "drought_stress_genes", "gene_expression",
                          "secondary_metabolites", "membrane_permeability",
                          "sheath_blight_incidence"],
    "stress_resilience_disease": ["sheath_blight_incidence"],
    "yield": ["yield_kg_ha"],
}


# =============================================================================
# Core estimator
# =============================================================================

@dataclass
class EffectEstimate:
    """Result from the sound effect estimator."""
    effect_pct: float
    effect_direction: str
    confidence: float
    supporting_row_count: int
    matched_windows: List[str]
    notes: str


def _check_biophysical_boundaries(frequency_hz, spl_db, stage_bucket="mixed",
                                    has_spl_matched_data=False):
    """
    Layer 1: Apply hard biophysical constraints.
    
    When has_spl_matched_data=True, the caution zone (90-110 dB) is
    bypassed because the Lorentzian already has data measured at this SPL.
    Only hard limits (activation floor <60 dB, damage ceiling >110 dB)
    remain as absolute constraints.
    """
    notes = []

    if spl_db > 0 and spl_db < ACTIVATION_FLOOR_DB:
        notes.append(f"SPL {spl_db:.0f} dB below {ACTIVATION_FLOOR_DB:.0f} dB activation floor")
        return 0.0, notes

    if spl_db > DAMAGE_CEILING_DB:
        notes.append(f"SPL {spl_db:.0f} dB above {DAMAGE_CEILING_DB:.0f} dB damage ceiling — membrane damage expected")
        return -1.0, notes

    # Caution zone (90-110 dB): only when extrapolating beyond studied SPLs.
    if spl_db > GRASS_INHIBITION_DB and stage_bucket not in ("seed", "germination"):
        if not has_spl_matched_data:
            reduction = (spl_db - GRASS_INHIBITION_DB) / (DAMAGE_CEILING_DB - GRASS_INHIBITION_DB)
            modifier = max(0.3, 1.0 - reduction * 0.5)
            notes.append(f"SPL {spl_db:.0f} dB caution zone, no matched data — modifier={modifier:.2f}")
            return modifier, notes

    if AIRBORNE_CONFIDENCE_DECAY_START_HZ < frequency_hz < ULTRASOUND_FLOOR_HZ:
        if frequency_hz >= AIRBORNE_CONFIDENCE_DECAY_END_HZ:
            notes.append(f"Frequency {frequency_hz:.0f} Hz: no documented airborne effect above {AIRBORNE_CONFIDENCE_DECAY_END_HZ:.0f} Hz")
            return 0.0, notes
        elif frequency_hz > AIRBORNE_CONFIDENCE_DECAY_START_HZ:
            decay = 1.0 - (frequency_hz - AIRBORNE_CONFIDENCE_DECAY_START_HZ) / (AIRBORNE_CONFIDENCE_DECAY_END_HZ - AIRBORNE_CONFIDENCE_DECAY_START_HZ)
            notes.append(f"High frequency {frequency_hz:.0f} Hz: confidence scaled to {decay:.0%}")
            return decay, notes

    return 1.0, notes


def _match_windows(frequency_hz, spl_db, stage_bucket, outcome_name,
                   stress_context="well_watered"):
    """
    Layer 2: Find matching empirical windows.
    Returns list of (window, match_quality) tuples.
    """
    target_outcomes = OUTCOME_MAPPING.get(outcome_name, [outcome_name])
    matches = []

    for w in ALL_WINDOWS:
        # Frequency match
        if frequency_hz < w.freq_low_hz or frequency_hz > w.freq_high_hz:
            continue

        # SPL match (if window specifies SPL bounds)
        if w.spl_low_db is not None and spl_db > 0:
            if spl_db < w.spl_low_db or spl_db > w.spl_high_db:
                continue

        # Stage match
        stage_match = stage_bucket in w.stages or "mixed" in w.stages
        if not stage_match and stage_bucket not in ["mixed", "any"]:
            continue

        # Outcome match
        outcome_match = any(o in target_outcomes for o in w.outcomes)
        if not outcome_match and outcome_name not in w.outcomes:
            continue

        # Calculate match quality
        quality = w.confidence

        # Cross-species window discount.
        # Unlike CSV data (which gets concordance-checked against the rice
        # Lorentzian), windows represent structural biophysical boundaries
        # (e.g., >90 dB grass inhibition, 1000 Hz defense shift).
        # These are about conserved membrane physics, not species-specific
        # resonance, so a moderate flat discount is appropriate.
        if w.crop_tier != "rice":
            quality *= 0.6  # Moderate, not the old 0.4 blanket

        # Stage precision bonus
        if stage_bucket in w.stages:
            quality *= 1.0
        else:
            quality *= 0.6

        # Frequency centrality — prefer center of window
        window_center = (w.freq_low_hz + w.freq_high_hz) / 2
        window_width = w.freq_high_hz - w.freq_low_hz
        if window_width > 0:
            dist_from_center = abs(frequency_hz - window_center) / window_width
            centrality = max(0.3, 1.0 - dist_from_center)
            quality *= centrality

        # Drought context bonus for drought windows
        if stress_context == "drought" and "drought" in (w.notes + w.source).lower():
            quality *= 1.3

        matches.append((w, quality))

    return matches


def _fit_lorentzian(points):
    """
    Fit a Lorentzian R(f) = A / (1 + ((f-f0)/γ)²) to data points.
    points: list of (freq, effect, weight).
    Returns (f0, A, gamma, fit_confidence).
    
    The grid search covers the full frequency range of the data,
    not just ±200 Hz from center. This is critical for KPIs like
    water_status where data spans 250-1500 Hz.
    """
    n = len(points)
    if n == 0:
        return 0, 0, 150, 0

    if n == 1:
        return points[0][0], points[0][1], 150.0, points[0][2] * 0.6

    freqs = [p[0] for p in points]
    effects = [abs(p[1]) for p in points]
    f_min, f_max = min(freqs), max(freqs)
    f_center = sum(f * w for f, e, w in points) / max(0.001, sum(w for _, _, w in points))
    max_eff = max(effects) if effects else 30

    # Search range covers the full data extent with padding
    f0_lo = max(20, int(f_min - 100))
    f0_hi = int(f_max + 100)
    f0_step = max(5, (f0_hi - f0_lo) // 80)  # ~80 steps across range

    # Amplitude range covers observed effects
    A_hi = min(200, int(max_eff * 1.5) + 10)
    A_step = max(1, A_hi // 40)

    # Gamma range: narrow peaks (40 Hz) to very broad (1000 Hz)
    gamma_hi = max(400, int((f_max - f_min) * 0.8))

    best_err = float('inf')
    best_params = (f_center, 30.0, 150.0)

    for f0_try in range(f0_lo, f0_hi, f0_step):
        for A_try in range(1, A_hi, A_step):
            for gamma_try in range(40, gamma_hi, 20):
                err = 0
                for f, eff, w in points:
                    predicted = A_try / (1 + ((f - f0_try) / gamma_try) ** 2)
                    err += w * (predicted - eff) ** 2
                if err < best_err:
                    best_err = err
                    best_params = (f0_try, A_try, gamma_try)

    f0, A, gamma = best_params
    rmse = math.sqrt(best_err / max(1, sum(w for _, _, w in points)))
    fit_confidence = min(1.0, sum(w for _, _, w in points) / n) * max(0.3, 1.0 - rmse / 20.0)
    return f0, A, gamma, fit_confidence


def _concordance_discount(row_freq, row_effect, rice_f0, rice_A, rice_gamma):
    """
    Compute a concordance-based confidence discount for a cross-species
    data point, based on how well it fits the rice Lorentzian trend.

    If the cross-species data AGREES with the rice resonance curve
    (same direction, similar magnitude), it gets near-full confidence
    because it reinforces the established trend.

    If it DISAGREES (opposite direction, or wildly different magnitude),
    it gets heavily discounted because it may reflect species-specific
    morphological differences rather than conserved mechanisms.

    Returns a multiplier between 0.15 (strong disagreement) and 1.0 (perfect agreement).
    """
    if rice_A == 0 or rice_gamma == 0:
        # No rice trend established — fall back to moderate discount
        return 0.5

    # What the rice curve predicts at this frequency
    rice_prediction = rice_A / (1 + ((row_freq - rice_f0) / rice_gamma) ** 2)

    # Direction concordance: do they agree on sign?
    same_direction = (row_effect > 0 and rice_prediction > 0.5) or \
                     (row_effect < 0 and rice_prediction < -0.5) or \
                     (abs(row_effect) < 1 and abs(rice_prediction) < 1)

    if not same_direction:
        # Opposite directions — this is a species-specific divergence
        # (e.g., oat inhibition at 300 Hz vs rice benefit)
        return 0.15

    # Magnitude concordance: how close is the magnitude?
    if rice_prediction > 0.5:
        ratio = min(row_effect, rice_prediction) / max(row_effect, rice_prediction, 0.1)
        # ratio near 1.0 = similar magnitude, near 0 = very different
        return 0.5 + 0.5 * max(0, ratio)
    else:
        # Both near zero or rice predicts nothing here — moderate trust
        # Cross-species data in a region where rice has no strong prediction
        # gets moderate credit for providing directional information
        return 0.6


def _lookup_csv_evidence(frequency_hz, spl_db, outcome_name, stage_bucket,
                         stress_context, data_dir):
    """
    Layer 3: Lorentzian resonance model + CSV evidence.

    Two-pass approach:
      Pass 1: Fit Lorentzian to RICE-ONLY data points.
      Pass 2: Score cross-species data by concordance with the rice
              curve, then re-fit with all data weighted accordingly.

    This replaces the blanket 40% cross-species discount with a
    trend-concordance check: cross-species data that reinforces the
    rice trend gets near-full confidence, data that contradicts it
    gets heavily discounted.

    Physics basis: the mechanosensory machinery (MS ion channels, Ca²⁺
    signalling) is conserved across monocot grasses. The optimal
    frequency may differ due to morphology, but the response direction
    at a given frequency should be similar if the mechanism is the same.

    Returns (effect, confidence, row_count, notes).
    """
    data = get_all_data(data_dir)

    # Use kpi_category from CSV directly, with OUTCOME_MAPPING as fallback
    target_kpi_cats = KPI_CATEGORIES.get(outcome_name, [outcome_name.upper()])
    target_outcomes = OUTCOME_MAPPING.get(outcome_name, [outcome_name])

    CONF_MAP = {"high": 1.0, "medium": 0.6, "low": 0.3}
    COMP_MAP = {"full": 1.0, "partial": 0.5, "insufficient": 0.2}

    # === Pass 1: Collect rice-only and cross-species points separately ===
    rice_points = []
    cross_points = []

    for row in data:
        if not row.usable_for_calibration:
            continue
        if row.frequency_hz <= 0:
            continue

        # Match by kpi_category (primary) or outcome_name (fallback)
        matched = False
        if row.kpi_category and row.kpi_category in target_kpi_cats:
            matched = True
        elif row.outcome_name in target_outcomes:
            matched = True
        if not matched:
            continue

        # Apply sub-outcome filter if this is a sub-KPI query
        sub_filter = SUB_OUTCOME_FILTERS.get(outcome_name)
        if sub_filter:
            if sub_filter["include"] and row.outcome_name not in sub_filter["include"]:
                continue
            if sub_filter["exclude"] and row.outcome_name in sub_filter["exclude"]:
                continue

        # Use calculated_effect_size_pct (primary) or reported_effect_size
        eff = row.calculated_effect_size_pct
        if eff == 0:
            eff = row.reported_effect_size
        if eff == 0:
            continue

        # Sign normalisation: some outcomes are measured such that
        # negative = beneficial (e.g., germination time -50% = faster,
        # disease incidence -50% = healthier, membrane permeability -13% = less damage).
        # Flip these so the Lorentzian sees all beneficial effects as positive.
        FLIP_SIGN_OUTCOMES = {
            'germination_time', 'germination_time_h', 'speed_of_germination',
            'sheath_blight_incidence', 'membrane_permeability',
        }
        if row.outcome_name in FLIP_SIGN_OUTCOMES:
            eff = -eff

        conf = CONF_MAP.get(row.confidence_weight, 0.3) * COMP_MAP.get(row.completeness_class, 0.2)

        # SPL proximity weighting — only when both sides have known SPL
        spl_w = 1.0
        if row.spl_db_known and spl_db > 0 and row.spl_db > 0:
            spl_diff = abs(spl_db - row.spl_db)
            spl_w = max(0.1, 1.0 - spl_diff / 30.0)

        point = (row.frequency_hz, eff, conf * spl_w)

        if row.crop_type == "rice":
            rice_points.append(point)
        else:
            cross_points.append((row.frequency_hz, eff, conf * spl_w, row.crop_type))

    # === Separate audible from ultrasound data ===
    # Ultrasound (≥20 kHz) operates via liquid cavitation, not airborne resonance.
    # The Lorentzian is physics-correct for airborne resonance only.
    # Ultrasound data gets direct distance-weighted lookup instead.
    #
    # Ultrasound rows are often marked usable_for_calibration=no because they
    # shouldn't feed the Lorentzian — but they ARE valid for ultrasound-specific
    # calibration. We collect them separately from ALL data, not just usable rows.
    audible_rice = [p for p in rice_points if p[0] < ULTRASOUND_FLOOR_HZ]
    audible_cross = [(f, e, c, cr) for f, e, c, cr in cross_points if f < ULTRASOUND_FLOOR_HZ]

    # Collect ultrasound points from ALL rows (bypass usability filter)
    # Apply outcome-specific routing so germination tests use germination
    # data and vigor tests use vigor data, not pooled together.
    ULTRASOUND_GERMINATION_OUTCOMES = {
        'germination_time', 'germination_time_h', 'germination_pct',
        'germination_speed', 'germination_potential', 'germination_energy',
        'germination_index', 'speed_of_germination', 'water_uptake_rate',
    }
    ULTRASOUND_VIGOR_OUTCOMES = {
        'seedling_dry_weight', 'seedling_dry_weight_g', 'plant_height',
        'root_length', 'fresh_weight', 'shoot_dry_weight',
    }

    ultrasound_points = []
    for row in data:
        if row.frequency_hz < ULTRASOUND_FLOOR_HZ:
            continue
        # Match by kpi_category or outcome
        matched = False
        if row.kpi_category and row.kpi_category in target_kpi_cats:
            matched = True
        elif row.outcome_name in target_outcomes:
            matched = True
        if not matched:
            continue

        # Outcome-specific routing for ultrasound
        if outcome_name == "germination" and row.outcome_name not in ULTRASOUND_GERMINATION_OUTCOMES:
            continue
        if outcome_name == "vigor" and row.outcome_name not in ULTRASOUND_VIGOR_OUTCOMES:
            continue

        eff = row.calculated_effect_size_pct
        if eff == 0:
            eff = row.reported_effect_size
        if eff == 0:
            continue

        # Sign normalisation for ultrasound outcomes
        FLIP_SIGN_OUTCOMES = {
            'germination_time', 'germination_time_h', 'speed_of_germination',
        }
        if row.outcome_name in FLIP_SIGN_OUTCOMES:
            eff = -eff

        conf = CONF_MAP.get(row.confidence_weight, 0.3) * COMP_MAP.get(row.completeness_class, 0.2)
        if row.crop_type != "rice":
            conf *= 0.6
        ultrasound_points.append((row.frequency_hz, eff, conf))

    # === If query is ultrasound, use direct lookup ===
    if frequency_hz >= ULTRASOUND_FLOOR_HZ:
        if not ultrasound_points:
            return 0.0, 0.0, 0, ["No ultrasound data for this outcome"]
        # Distance-weighted average of ultrasound points
        total_w = 0
        total_e = 0
        for f, eff, w in ultrasound_points:
            freq_dist = abs(frequency_hz - f)
            prox = max(0.05, 1.0 - freq_dist / 30000.0)
            total_w += w * prox
            total_e += eff * w * prox
        if total_w > 0:
            effect = total_e / total_w
            conf = min(1.0, total_w / len(ultrasound_points))
        else:
            effect, conf = 0, 0
        return effect, conf, len(ultrasound_points), [f"Ultrasound direct lookup ({len(ultrasound_points)} points)"]

    # === Nearest-neighbour bypass for oscillating data ===
    # iWUE data oscillates wildly across 30 Hz (±60%) due to standing-wave
    # artefacts. A Lorentzian cannot fit this pattern. Use direct nearest-
    # neighbour lookup instead: return the effect of the closest data point.
    NEAREST_NEIGHBOR_OUTCOMES = {'water_status_efficiency'}
    if outcome_name in NEAREST_NEIGHBOR_OUTCOMES and audible_rice:
        best_dist = float('inf')
        best_eff = 0
        best_conf = 0
        for f, eff, w in audible_rice:
            d = abs(frequency_hz - f)
            if d < best_dist:
                best_dist = d
                best_eff = eff
                best_conf = w
        # Confidence decays with distance from nearest point
        decay = max(0.1, 1.0 - best_dist / 50.0)  # 50 Hz range
        return best_eff, best_conf * decay, len(audible_rice), \
            [f"Nearest-neighbour: {best_eff:+.1f}% at {best_dist:.0f} Hz distance"]

    # === Pass 1b: Fit Lorentzian to rice-only AUDIBLE data ===
    if audible_rice:
        rice_f0, rice_A, rice_gamma, _ = _fit_lorentzian(audible_rice)
    else:
        rice_f0, rice_A, rice_gamma = 0, 0, 150

    # === Pass 2: Score cross-species points by concordance ===
    all_points = list(audible_rice)

    for freq, eff, base_conf, crop in audible_cross:
        discount = _concordance_discount(freq, eff, rice_f0, rice_A, rice_gamma)
        all_points.append((freq, eff, base_conf * discount))

    if not all_points:
        return 0.0, 0.0, 0, []

    # === Fit final Lorentzian to all concordance-weighted data ===
    f0, A, gamma, fit_confidence = _fit_lorentzian(all_points)

    # Evaluate at query frequency
    effect = A / (1 + ((frequency_hz - f0) / gamma) ** 2)

    # Confidence decays with distance from data
    all_freqs = [p[0] for p in all_points]
    min_data_dist = min(abs(frequency_hz - f) for f in all_freqs)
    distance_penalty = max(0.1, 1.0 - min_data_dist / (gamma * 3))
    final_confidence = fit_confidence * distance_penalty

    n_rice = len(audible_rice)
    n_cross = len(audible_cross)
    concordant = sum(1 for f, e, c, cr in audible_cross
                     if _concordance_discount(f, e, rice_f0, rice_A, rice_gamma) > 0.5)

    notes = [f"Lorentzian: f0={f0:.0f}Hz, A={A:.0f}%, γ={gamma:.0f}Hz"]
    if n_cross > 0:
        notes.append(f"Cross-species: {concordant}/{n_cross} concordant with rice trend")

    return effect, final_confidence, n_rice + n_cross, notes


def estimate_effect(
    outcome_name: str,
    frequency_hz: float,
    spl_db: float,
    hours_per_day: float,
    stage_bucket: str = "mixed",
    stress_context: str = "well_watered",
    data_dir: str = None,
) -> EffectEstimate:
    """
    Estimate the effect of sound treatment on a specific rice outcome.

    Architecture (3-layer):
      1. Biophysical boundaries — hard constraints
      2. Empirical windows — discrete literature-based bands
      3. CSV evidence — granular refinement

    Returns EffectEstimate with effect_pct, direction, confidence, and notes.
    """
    # Handle no-treatment case
    if hours_per_day <= 0 or (spl_db <= 0 and frequency_hz < ULTRASOUND_FLOOR_HZ):
        return EffectEstimate(0.0, "none", 0.0, 0, [], "No treatment applied.")

    # === Layer 3 first: CSV evidence (need to know if SPL-matched data exists) ===
    csv_effect, csv_weight, csv_count, csv_notes = _lookup_csv_evidence(
        frequency_hz, spl_db, outcome_name, stage_bucket, stress_context, data_dir)

    # Check if we have data near this SPL in the CSV
    has_spl_matched = csv_count > 0 and csv_weight > 0.05

    # === Layer 1: Biophysical boundaries ===
    bio_modifier, bio_notes = _check_biophysical_boundaries(
        frequency_hz, spl_db, stage_bucket, has_spl_matched_data=has_spl_matched)

    if bio_modifier == 0.0:
        return EffectEstimate(0.0, "none", 0.0, 0, [],
                              "; ".join(bio_notes) or "Below activation threshold.")

    if bio_modifier < 0:
        return EffectEstimate(
            effect_pct=-15.0,
            effect_direction="decrease",
            confidence=0.45,
            supporting_row_count=0,
            matched_windows=["biophysical_damage_ceiling"],
            notes="; ".join(bio_notes) + " — Hard biophysical override: cellular damage expected.",
        )

    # === Layer 2: Empirical windows ===
    window_matches = _match_windows(frequency_hz, spl_db, stage_bucket,
                                     outcome_name, stress_context)

    # === Combine layers ===
    all_notes = list(bio_notes)
    matched_window_ids = []

    if not window_matches and csv_count == 0:
        return EffectEstimate(0.0, "none", 0.0, 0, [],
                              "No matching empirical windows or resonance data.")

    # Check if any inhibitory window with explicit SPL criteria matched.
    # Inhibitory windows with specific SPL bounds take priority over
    # positive windows without SPL bounds (e.g., injury at 4000 Hz/111 dB
    # overrides Sri Lankan positive window at 3-5 kHz with no SPL check).
    has_inhibitory_spl_match = any(
        w.direction == "decrease" and w.spl_low_db is not None
        for w, q in window_matches
    )

    if has_inhibitory_spl_match:
        # Keep only inhibitory windows
        window_matches = [(w, q) for w, q in window_matches if w.direction == "decrease"]

    # Window-based estimate: weighted average of window midpoints
    window_effect = 0.0
    window_confidence = 0.0
    if window_matches:
        total_q = sum(q for _, q in window_matches)
        for w, q in window_matches:
            midpoint = (w.effect_pct_low + w.effect_pct_high) / 2
            window_effect += midpoint * q
            matched_window_ids.append(w.window_id)
        if total_q > 0:
            window_effect /= total_q
            window_confidence = min(1.0, total_q / len(window_matches))

    # Blend window and resonance estimates
    # When no windows match, trust CSV 100% — no reason to dilute.
    # When windows match, blend based on data density.
    if csv_count > 0 and csv_weight > 0:
        csv_confidence = csv_weight
        if not window_matches:
            # No window guidance — CSV drives everything
            final_effect = csv_effect
            final_confidence = csv_confidence
        elif csv_count >= 3:
            # Strong Lorentzian + windows for directional check
            blend_csv = 0.85
            blend_win = 0.15
            final_effect = csv_effect * blend_csv + window_effect * blend_win
            final_confidence = csv_confidence * blend_csv + window_confidence * blend_win
        else:
            blend_csv = 0.65
            blend_win = 0.35
            final_effect = csv_effect * blend_csv + window_effect * blend_win
            final_confidence = csv_confidence * blend_csv + window_confidence * blend_win
    elif window_matches:
        final_effect = window_effect
        final_confidence = window_confidence
    else:
        final_effect = 0.0
        final_confidence = 0.0

    # Apply biophysical modifier (can flip sign for damage zone)
    final_effect *= bio_modifier

    # Dosage scaling — only for AIRBORNE sound, not ultrasound.
    # Ultrasound seed treatments have completely different duration profiles
    # (minutes of liquid soaking, not hours of chronic airborne exposure).
    # The ultrasound data already reflects the study's treatment duration.
    if hours_per_day > 0 and frequency_hz < ULTRASOUND_FLOOR_HZ:
        if hours_per_day < 1.0:
            dose_factor = math.sqrt(hours_per_day)
        elif hours_per_day <= 4.0:
            dose_factor = 1.0
        else:
            dose_factor = 1.0 + 0.1 * math.log(hours_per_day / 4.0)
            dose_factor = min(dose_factor, 1.3)
        final_effect *= dose_factor

    # Cap extreme effects with confidence penalty
    # With v2 CSV, measured effects reach +89% (Jeong gsw) and -50% (Wang germ time).
    # Cap at ±100% as safety limit — anything beyond is extrapolation artefact.
    if abs(final_effect) > 100:
        final_confidence *= 0.5
        final_effect = max(-100, min(100, final_effect))
        all_notes.append("Effect capped at ±100%")

    # Direction
    if final_effect > 0.5:
        direction = "increase"
    elif final_effect < -0.5:
        direction = "decrease"
    else:
        direction = "none"

    final_confidence = round(max(0.0, min(1.0, final_confidence)), 3)

    return EffectEstimate(
        effect_pct=round(final_effect, 2),
        effect_direction=direction,
        confidence=final_confidence,
        supporting_row_count=csv_count,
        matched_windows=matched_window_ids,
        notes="; ".join(all_notes + [f"Windows: {', '.join(matched_window_ids)}" if matched_window_ids else "No window match"]),
    )


def estimate_all_effects(frequency_hz, spl_db, hours_per_day,
                         stage_bucket="mixed", stress_context="well_watered",
                         data_dir=None):
    """Estimate effects for all supported outcomes including sub-KPIs."""
    results = {}

    # Primary KPIs
    for outcome in ["germination", "vigor", "photosynthesis", "yield"]:
        results[outcome] = estimate_effect(
            outcome, frequency_hz, spl_db, hours_per_day,
            stage_bucket, stress_context, data_dir)

    # Water status: three independent Lorentzian tracks per CWSI theory.
    # gsw_drought  — active stomatal regulation under deficit (Jeong 2014)
    # gsw_wellwatered — stomatal response under irrigation (Jusoh 2023)
    # rwc — passive hydration state (Jeong 2014, slower/smaller response)
    # iWUE — water use efficiency (Jusoh 2023, nearest-neighbour)
    ws_gsw_d = estimate_effect("water_status_gsw_drought", frequency_hz, spl_db,
                                hours_per_day, stage_bucket, stress_context, data_dir)
    ws_gsw_w = estimate_effect("water_status_gsw_wellwatered", frequency_hz, spl_db,
                                hours_per_day, stage_bucket, stress_context, data_dir)
    ws_rwc = estimate_effect("water_status_rwc", frequency_hz, spl_db,
                              hours_per_day, stage_bucket, stress_context, data_dir)
    ws_eff = estimate_effect("water_status_efficiency", frequency_hz, spl_db,
                             hours_per_day, stage_bucket, stress_context, data_dir)

    # Store all sub-KPIs for validation routing
    results["water_status_gsw_drought"] = ws_gsw_d
    results["water_status_gsw_wellwatered"] = ws_gsw_w
    results["water_status_rwc"] = ws_rwc
    results["water_status_efficiency"] = ws_eff
    # Legacy aggregate for backward compatibility
    results["water_status_openness"] = ws_gsw_d if stress_context == "drought" else ws_gsw_w

    # Aggregate water_status for the UI: context-dependent weighting.
    # Under drought: gsw_drought is the primary indicator (CWSI-aligned),
    #   RWC is confirmatory, iWUE is secondary.
    # Under well-watered: iWUE and gsw_wellwatered are the relevant measures.
    if stress_context == "drought":
        components = [
            (ws_gsw_d, 0.6),   # Primary: stomatal regulation under deficit
            (ws_rwc, 0.25),    # Confirmatory: tissue hydration state
            (ws_eff, 0.15),    # Secondary: efficiency
        ]
    else:
        components = [
            (ws_gsw_w, 0.4),   # Stomatal openness
            (ws_eff, 0.4),     # Water use efficiency
            (ws_rwc, 0.2),     # Hydration state (less relevant when well-watered)
        ]

    agg_pct = sum(c.effect_pct * w for c, w in components)
    agg_conf = sum(c.confidence * w for c, w in components)
    agg_count = sum(c.supporting_row_count for c, _ in components)
    if any(c.effect_pct != 0 for c, _ in components):
        agg_dir = "increase" if agg_pct > 0.5 else "decrease" if agg_pct < -0.5 else "none"
        results["water_status"] = EffectEstimate(
            round(agg_pct, 2), agg_dir, round(agg_conf, 3), agg_count, [],
            f"CWSI-weighted: gsw primary under {'drought' if stress_context=='drought' else 'irrigation'}")
    else:
        results["water_status"] = EffectEstimate(0, "none", 0, 0, [], "No water status data")

    # Stress resilience: compute sub-KPIs
    sr_disease = estimate_effect("stress_resilience_disease", frequency_hz, spl_db,
                                 hours_per_day, stage_bucket, stress_context, data_dir)
    sr_general = estimate_effect("stress_resilience", frequency_hz, spl_db,
                                 hours_per_day, stage_bucket, stress_context, data_dir)

    results["stress_resilience_disease"] = sr_disease
    results["stress_resilience"] = sr_general

    return results


if __name__ == "__main__":
    print("=== Sound Response Engine v2 ===")
    for freq in [350, 550, 800, 2000, 5000, 15000]:
        effects = estimate_all_effects(freq, 70, 3.0, "seedling", "well_watered")
        print(f"\n{freq} Hz, 70 dB, seedling:")
        for name, est in effects.items():
            if est.effect_pct != 0 or est.confidence > 0:
                print(f"  {name}: {est.effect_pct:+.1f}% (conf={est.confidence:.2f}, windows={est.matched_windows})")
