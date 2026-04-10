"""
sound_response.py — Mechanistic GP sound-response engine (v4).

v4: mechanistic prior with mechanotransduction latent and pathway activations.

Architecture (four layers + two special pathways):

  1. BIOPHYSICAL BOUNDARIES
     Hard physics constraints: activation floor, damage ceiling, SPL caution zone.

  2. MECHANOTRANSDUCTION LATENT  S_mech,air(f, SPL)
     Captures mechanosensory activation as a function of frequency and SPL.
       S_mech,air = S_f(f) × S_SPL(SPL) × S_air(f)
     where:
       S_f(f)    = normalised MechActivation GP mean (frequency selectivity)
       S_SPL     = piecewise-linear SPL activation window [60, 85, 100, 110] dB
       S_air(f)  = high-frequency airborne coupling decay [0, 1]
     Ultrasound (>=20 kHz) bypasses this entirely.

  3. PATHWAY ACTIVATIONS + MECHANISTIC PRIOR (replaces empirical-window prior)
     Five biological pathways (seed_vigor, leaf_gas_exchange, drought_resilience,
     membrane_damage, ultrasound_cavitation) each have a log-Gaussian activation
     over frequency gated by stage, regime, and SPL. Pathway activations are
     combined with OUTCOME_PATHWAY_WEIGHTS to produce a per-outcome prior mean.
     This prior is biologically motivated: it encodes WHERE in the spectrum each
     outcome is expected to respond, not just WHICH discrete windows exist.

  4. GAUSSIAN PROCESS REGRESSION (Matérn ν=5/2, heteroscedastic, local density)
     GP fitted to residuals (observed - mechanistic_prior) with Matérn ν=5/2
     kernel, heteroscedastic per-study noise, local-density confidence.
     Frequency aggregation prevents pseudo-replication.

  Special pathways:
     - Ultrasound (>=20 kHz): cavitation mechanism, separate lookup
     - iWUE (water_status_efficiency): nearest-neighbour (standing-wave data)

Sources:
  - Bochu et al. 2003 (rice seed, 400 Hz / 4000 Hz)
  - Jusoh & Ramlee 2023 (rice seedling, 350-380 Hz)
  - Jeong et al. 2014 (rice drought, 250-1500 Hz, stomata + RWC)
  - Hassan et al. 2014 (rice drought, 800-1500 Hz)
  - Qi et al. 2010 (rice chlorophyll, 550 Hz)
  - Hou et al. 2009 (rice yield, PAFT 550 Hz)
  - Munasinghe et al. 2023 (rice stomatal, 350 Hz)
  - Sri Lankan rice 2024 (3000-5000 Hz vegetative)
  - Kim et al. 2019 (flavonoid content, 250/800/1000 Hz)
  - Cross-species: wheat, corn, oat, barley, sorghum, general grasses
  See MODEL_NOTES.md and VALIDATION_CASES.yaml for full source list.
"""

import csv
import math
import os
import warnings
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Tuple

import numpy as np

try:
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import Matern, ConstantKernel
    _SKLEARN_AVAILABLE = True
except ImportError:
    _SKLEARN_AVAILABLE = False
    warnings.warn(
        "scikit-learn not found. sound_response will use mechanistic-prior-only "
        "predictions with reduced accuracy. Install with: pip install scikit-learn",
        ImportWarning,
        stacklevel=2,
    )


# =============================================================================
# Layer 1: BIOPHYSICAL BOUNDARIES
# =============================================================================

ACTIVATION_FLOOR_DB = 60.0
DAMAGE_CEILING_DB = 110.0
GRASS_INHIBITION_DB = 90.0
ACTIVE_BAND_LOW_HZ = 100.0
ACTIVE_BAND_HIGH_HZ = 2000.0
AIRBORNE_CONFIDENCE_DECAY_START_HZ = 5000.0
AIRBORNE_CONFIDENCE_DECAY_END_HZ = 15000.0
ULTRASOUND_FLOOR_HZ = 20000.0

# SPL window for S_SPL (mechanotransduction activation)
_SPL_L_MIN = 60.0     # activation floor
_SPL_L_OPT = 85.0     # plateau start
_SPL_L_MAX = 100.0    # plateau end / damage ramp start
_SPL_L_DAM = 110.0    # damage ceiling (full inhibition)


# =============================================================================
# Layer 2: MECHANOTRANSDUCTION LATENT
# =============================================================================

def spl_activation_S_SPL(spl_db: float) -> float:
    """
    S_SPL(SPL) — piecewise linear SPL activation gate ∈ [0, 1].

    Models the mechanosensory activation window:
      - 0 below 60 dB (below membrane deflection threshold)
      - linear ramp [0→1] from 60 to 85 dB
      - 1.0 plateau from 85 to 100 dB (optimal mechanosensory range)
      - linear ramp [1→0] from 100 to 110 dB (damage zone onset)
      - 0 above 110 dB (cellular damage — handled by biophysical boundary)

    Note: the hard damage-ceiling behaviour (returning -15% effect) is
    handled separately in _check_biophysical_boundaries. This function
    only provides the continuous activation gate used by S_mech,air.
    """
    if spl_db <= _SPL_L_MIN:
        return 0.0
    if spl_db <= _SPL_L_OPT:
        return (spl_db - _SPL_L_MIN) / (_SPL_L_OPT - _SPL_L_MIN)
    if spl_db <= _SPL_L_MAX:
        return 1.0
    if spl_db < _SPL_L_DAM:
        return 1.0 - (spl_db - _SPL_L_MAX) / (_SPL_L_DAM - _SPL_L_MAX)
    return 0.0


def air_coupling_S_air(frequency_hz: float) -> float:
    """
    S_air(f) — high-frequency airborne coupling decay ∈ [0, 1].

    Above 5 kHz, airborne mechanical coupling to plant tissues weakens
    as wavelengths become shorter than leaf dimensions. This decay is
    applied to airborne frequencies only; ultrasound uses the cavitation
    pathway which bypasses this.

      - 1.0 for f <= 5000 Hz
      - linear decay to 0.0 between 5000 and 15000 Hz
      - 0.0 for f >= 15000 Hz
    """
    if frequency_hz >= ULTRASOUND_FLOOR_HZ:
        return 0.0  # ultrasound bypasses airborne pathway
    if frequency_hz <= AIRBORNE_CONFIDENCE_DECAY_START_HZ:
        return 1.0
    if frequency_hz >= AIRBORNE_CONFIDENCE_DECAY_END_HZ:
        return 0.0
    return 1.0 - ((frequency_hz - AIRBORNE_CONFIDENCE_DECAY_START_HZ) /
                   (AIRBORNE_CONFIDENCE_DECAY_END_HZ - AIRBORNE_CONFIDENCE_DECAY_START_HZ))


@dataclass
class MechTransductionState:
    """
    State object for the MechActivation GP.
    Produced by train_mechactivation_gp; passed into smech_air.
    """
    gp: Optional[object] = None         # fitted GaussianProcessRegressor or None
    y_max: float = 100.0                # normalisation: 95th-pct GP mean on 50-5000 Hz
    prior_points: list = field(default_factory=list)  # (freq, effect) for IDW fallback
    n_training: int = 0
    available: bool = False


# MechActivation outcomes pool: outcomes that reflect cellular stress/activation
# status rather than phenotypic growth. These are pooled into one GP over
# frequency to capture the mechanosensory responsiveness curve.
_MECH_KPI_CATS = {"STRESS_RESILIENCE", "DROUGHT_TOLERANCE", "PHOTOSYNTHESIS"}
_MECH_OUTCOMES = {
    # Fv/Fm: photosystem efficiency — direct mechanosensory indicator
    "quantum_yield_FvFm",
    # Stress metabolites (sign: beneficial = positive)
    "total_flavonoid_content",
    "SOD_CAT_POD_enzyme_activity",
    "secondary_metabolites",
    # Disease resilience
    "sheath_blight_incidence",
    # Drought stress response
    "relative_water_content",
    "stomatal_conductance_drought",
    "drought_stress_genes",
    # Membrane stress (note: stored negative = some damage at high SPL)
    "cell_membrane_permeability",
    # Photosynthesis as a proxy for chloroplast mechanosensory response
    "photosynthesis_rate",
    "chlorophyll_content",
}
# Outcomes that are beneficial when negative — flip so GP sees all as positive
_MECH_FLIP = {"cell_membrane_permeability", "germination_time",
               "germination_time_h", "speed_of_germination"}


def _collect_mech_points(data) -> List[Tuple[float, float, float]]:
    """
    Collect (freq, mech_activation_pct, noise_alpha) for MechActivation GP.

    Pools sign-normalised stress/activation outcomes from the CSV.
    Only audible frequencies (<20 kHz); ultrasound cavitation is separate.
    """
    _CONF_MAP = {"high": 1.0, "medium": 0.6, "low": 0.3}
    _COMP_MAP = {"full": 1.0, "partial": 0.5, "insufficient": 0.2}
    _BASE_NOISE = 40.0

    # Allow all rows (usable or not) for MechActivation, but down-weight non-usable
    points = []
    for row in data:
        if row.frequency_hz <= 0 or row.frequency_hz >= ULTRASOUND_FLOOR_HZ:
            continue
        if row.outcome_name not in _MECH_OUTCOMES:
            continue

        eff = row.calculated_effect_size_pct or row.reported_effect_size
        if not isinstance(eff, (int, float)) or eff == 0.0:
            continue
        if row.outcome_name in _MECH_FLIP:
            eff = -eff

        conf_w = _CONF_MAP.get(row.confidence_weight, 0.3)
        comp_w = _COMP_MAP.get(row.completeness_class, 0.2)
        usable_mult = 1.0 if row.usable_for_calibration else 2.0  # higher noise if not primary
        species_mult = 1.0 if row.crop_type == "rice" else 2.78
        reliability = max(0.01, conf_w * comp_w)
        alpha = min(500.0, _BASE_NOISE * (1.0 / reliability) * usable_mult * species_mult)

        points.append((row.frequency_hz, float(eff), alpha))

    return points


def train_mechactivation_gp(data, data_dir=None) -> MechTransductionState:
    """
    Train a Matérn ν=5/2 GP over frequency for the MechActivation KPI.

    The MechActivation GP captures the frequency-selectivity of the plant's
    mechanosensory response by pooling all stress/activation outcomes. Its
    output S_f(f) ∈ [0,1] reflects how responsive the plant's stress pathways
    are at a given frequency, regardless of which specific outcome is queried.

    Returns MechTransductionState (available=False if insufficient data).
    """
    state = MechTransductionState()
    points = _collect_mech_points(data)
    if len(points) < 2:
        return state  # not enough data for a GP

    # Aggregate same-frequency points (precision-weighted) to prevent pseudo-replication
    from collections import defaultdict
    freq_groups = defaultdict(list)
    for f, e, a in points:
        freq_groups[round(f, 0)].append((f, e, a))

    agg = []
    for pts in freq_groups.values():
        precs = [1.0 / max(p[2], 1e-6) for p in pts]
        tot = sum(precs)
        mf = sum(p[0] * pr for p, pr in zip(pts, precs)) / tot
        me = sum(p[1] * pr for p, pr in zip(pts, precs)) / tot
        ma = max(1.0, 1.0 / tot)
        agg.append((mf, me, ma))

    state.prior_points = [(f, e) for f, e, _ in agg]
    state.n_training = len(agg)

    if not _SKLEARN_AVAILABLE or len(agg) < 2:
        state.available = False
        return state

    freqs = np.array([p[0] for p in agg])
    effects = np.array([p[1] for p in agg])
    alphas = np.array([p[2] for p in agg])

    amp_init = max(float(np.mean(effects ** 2)), 1.0)
    if len(agg) < 3:
        kernel = (ConstantKernel(amp_init, constant_value_bounds="fixed")
                  * Matern(200.0, length_scale_bounds="fixed", nu=2.5))
        n_restarts = 0
    else:
        kernel = (ConstantKernel(amp_init, constant_value_bounds=(amp_init * 0.01, amp_init * 100))
                  * Matern(200.0, length_scale_bounds=(100.0, 1200.0), nu=2.5))
        n_restarts = 3

    gp = GaussianProcessRegressor(kernel=kernel, alpha=alphas,
                                   normalize_y=False, n_restarts_optimizer=n_restarts,
                                   random_state=42)
    try:
        gp.fit(freqs.reshape(-1, 1), effects)
    except Exception:
        state.available = False
        return state

    # Compute Y_max: 95th-pct of GP mean over 50-5000 Hz (audible range)
    grid = np.linspace(50, 5000, 200).reshape(-1, 1)
    means = gp.predict(grid)
    pos_means = means[means > 0]
    if len(pos_means) > 0:
        state.y_max = float(np.percentile(pos_means, 95))
    else:
        state.y_max = 30.0  # fallback

    state.gp = gp
    state.available = True
    return state


def mechactivation_mean(frequency_hz: float, state: MechTransductionState) -> float:
    """
    Return the MechActivation GP mean at frequency_hz (in percent).
    Falls back to inverse-distance-weighted average if GP unavailable.
    """
    if state.available and state.gp is not None:
        mean_val = float(state.gp.predict([[frequency_hz]])[0])
        return mean_val
    # IDW fallback
    if not state.prior_points:
        return 0.0
    weights = [1.0 / max(abs(frequency_hz - f), 50.0) for f, _ in state.prior_points]
    tot = sum(weights)
    return sum(w * e for (_, e), w in zip(state.prior_points, weights)) / tot


def mechactivation_Sf(frequency_hz: float, state: MechTransductionState) -> float:
    """
    S_f(f) ∈ [0, 1] — normalised mechanosensory frequency selectivity.

    Clips y_mech(f) / Y_max into [0, 1]. Frequencies where the plant shows
    little stress-pathway activation return values near 0; frequencies where
    measured stress/activation effects are large return values near 1.
    """
    y = mechactivation_mean(frequency_hz, state)
    if state.y_max <= 0:
        return 0.0
    return max(0.0, min(1.0, y / state.y_max))


def smech_air(frequency_hz: float, spl_db: float,
               state: MechTransductionState) -> float:
    """
    S_mech,air(f, SPL) = S_f(f) × S_SPL(SPL) × S_air(f) ∈ [0, 1].

    Combines:
      - Frequency selectivity from the MechActivation GP
      - SPL activation gate (optimal window 85-100 dB)
      - Airborne high-frequency coupling decay

    Returns 0 for ultrasound frequencies (handled by cavitation pathway).
    """
    if frequency_hz >= ULTRASOUND_FLOOR_HZ:
        return 0.0
    sf = mechactivation_Sf(frequency_hz, state)
    s_spl = spl_activation_S_SPL(spl_db)
    s_air = air_coupling_S_air(frequency_hz)
    return sf * s_spl * s_air


# Module-level MechActivation state (lazy-initialised per data_dir)
_mech_state_cache: Dict[str, MechTransductionState] = {}


def _get_mech_state(data, data_dir: str) -> MechTransductionState:
    """Return the MechActivation GP state, training it once per data_dir."""
    key = data_dir or "__default__"
    if key not in _mech_state_cache:
        _mech_state_cache[key] = train_mechactivation_gp(data, data_dir)
    return _mech_state_cache[key]


# =============================================================================
# Layer 3: PATHWAY ACTIVATIONS + MECHANISTIC PRIOR
# =============================================================================

@dataclass
class PathwayConfig:
    """
    Configuration for one biological sound-response pathway.

    Each pathway has a log-Gaussian activation shape over frequency,
    gated by stage, water regime, medium, and minimum SPL.
    """
    name: str
    mu_log_hz: float        # log(center_freq) — natural log
    sigma_log_hz: float     # width of log-Gaussian
    applicable_stages: List[str]   # e.g. ["seedling", "vegetative"]
    applicable_regimes: List[str]  # ["well_watered", "drought", "any"]
    medium: str             # "airborne" or "ultrasound"
    min_spl_db: float = 0.0
    peak_effect_pct: float = 30.0  # amplitude of the pathway prior (%)


# Five pathways derived from the phytoacoustics literature
PATHWAYS: Dict[str, PathwayConfig] = {
    # Seed/vigor: Bochu 2003 at 400 Hz/106 dB; peak in germination and seedling
    "seed_vigor": PathwayConfig(
        name="seed_vigor",
        mu_log_hz=math.log(380),   # centre ~380 Hz (Jusoh/Bochu cluster)
        sigma_log_hz=0.30,         # broad: 200-700 Hz at 1σ
        applicable_stages=["seed", "germination", "seedling", "mixed"],
        applicable_regimes=["any", "well_watered", "drought"],
        medium="airborne",
        min_spl_db=0.0,
        peak_effect_pct=25.0,
    ),
    # Leaf gas exchange: Jusoh 2023 at 350-380 Hz, 60-80 dB, seedling/vegetative
    "leaf_gas_exchange": PathwayConfig(
        name="leaf_gas_exchange",
        mu_log_hz=math.log(357),   # centre ~357 Hz (Jusoh peak)
        sigma_log_hz=0.12,         # narrow: very selective around 300-430 Hz
        applicable_stages=["seedling", "vegetative", "mixed"],
        applicable_regimes=["any", "well_watered"],
        medium="airborne",
        min_spl_db=0.0,
        peak_effect_pct=40.0,      # matches Jusoh assimilation +43%
    ),
    # Drought resilience: Jeong 2014 at 250-1500 Hz; peaks near 800-1000 Hz
    "drought_resilience": PathwayConfig(
        name="drought_resilience",
        mu_log_hz=math.log(800),   # centre ~800 Hz (Jeong RWC/gsw cluster)
        sigma_log_hz=0.55,         # wide: 250-2500 Hz at 1σ (Jeong spans 250-1500)
        applicable_stages=["vegetative", "reproductive", "mixed"],
        applicable_regimes=["any", "drought", "well_watered"],
        medium="airborne",
        min_spl_db=0.0,
        peak_effect_pct=45.0,      # matches Jeong gsw +71% at 800 Hz
    ),
    # Membrane damage: Bochu 2003 at 4000 Hz/111 dB — inhibitory at high SPL
    "membrane_damage": PathwayConfig(
        name="membrane_damage",
        mu_log_hz=math.log(4000),  # centre ~4000 Hz (Bochu injury)
        sigma_log_hz=0.25,
        applicable_stages=["seed", "germination", "seedling", "vegetative", "mixed"],
        applicable_regimes=["any", "well_watered", "drought"],
        medium="airborne",
        min_spl_db=100.0,          # only active above 100 dB (high-SPL injury)
        peak_effect_pct=-20.0,     # negative = inhibition
    ),
    # Ultrasound cavitation: Wang 2022/2020 at 35-40 kHz seed treatments
    "ultrasound_cavitation": PathwayConfig(
        name="ultrasound_cavitation",
        mu_log_hz=math.log(38000), # centre ~38 kHz
        sigma_log_hz=0.18,
        applicable_stages=["seed", "germination"],
        applicable_regimes=["any", "well_watered", "drought"],
        medium="ultrasound",
        min_spl_db=0.0,
        peak_effect_pct=20.0,
    ),
}


def pathway_activation(
    pathway: PathwayConfig,
    frequency_hz: float,
    spl_db: float,
    stage_bucket: str,
    stress_context: str,
    medium: str,
    smech_air_value: float,
) -> float:
    """
    Compute activation p(x) for one pathway at the given conditions.

    For airborne pathways:
      p = log_gaussian(f) × stage_gate × regime_gate × spl_gate × smech_scale
    For ultrasound:
      p = log_gaussian(f) × stage_gate × regime_gate   (ignores S_mech)

    Returns a value in approximately [0, 1] before scaling by peak_effect_pct.
    """
    if frequency_hz <= 0:
        return 0.0

    # Stage gate
    if ("any" not in pathway.applicable_stages and
            stage_bucket not in pathway.applicable_stages):
        # Allow "mixed" as a wildcard query to pass all pathways
        if stage_bucket != "mixed":
            return 0.0

    # Regime gate
    if "any" not in pathway.applicable_regimes:
        if stress_context not in pathway.applicable_regimes:
            return 0.0

    # Medium gate
    if pathway.medium == "ultrasound" and frequency_hz < ULTRASOUND_FLOOR_HZ:
        return 0.0
    if pathway.medium == "airborne" and frequency_hz >= ULTRASOUND_FLOOR_HZ:
        return 0.0

    # Minimum SPL gate (for membrane_damage, only at very high SPL)
    if spl_db < pathway.min_spl_db:
        return 0.0

    # Log-Gaussian frequency activation
    log_f = math.log(max(frequency_hz, 1.0))
    z = (log_f - pathway.mu_log_hz) / pathway.sigma_log_hz
    log_gauss = math.exp(-0.5 * z * z)  # in [0, 1]

    if pathway.medium == "ultrasound":
        return log_gauss

    # Airborne: scale by S_mech,air (captures SPL + airborne coupling)
    # Use a soft minimum so the pathway is still detectable even when
    # S_mech is low (e.g. at very low SPL), but the magnitude is suppressed.
    smech_scale = max(0.15, smech_air_value) if spl_db >= _SPL_L_MIN else 0.0
    return log_gauss * smech_scale


# Outcome-pathway weights.
# Rules (from spec):
#   - Primary objective: preserve direction + approximate magnitude at anchor cases
#   - Do not allow cross-species windows to dominate (weight >0.3 without rice support)
#   - Adjust at most 0.1 at a time per weight; document reason
OUTCOME_PATHWAY_WEIGHTS: Dict[str, Dict[str, float]] = {
    # Photosynthesis: dominated by leaf gas exchange (Jusoh 2023 seedling +40%)
    # Secondary seed_vigor allows moderate prior at germination stages.
    "photosynthesis": {
        "leaf_gas_exchange": 0.8,
        "seed_vigor": 0.2,
    },
    # Vigor (height, root length, dry weight): seed_vigor is primary
    # leaf_gas_exchange provides the secondary overlap seen in Jusoh data.
    # drought_resilience at 0.2 avoids creating false peaks in 800-1500 Hz.
    "vigor": {
        "seed_vigor": 0.6,
        "leaf_gas_exchange": 0.2,
        "drought_resilience": 0.2,
    },
    # Water status sub-KPIs:
    "water_status_gsw_drought": {
        "drought_resilience": 1.0,
    },
    "water_status_gsw_wellwatered": {
        "leaf_gas_exchange": 1.0,
    },
    "water_status_rwc": {
        "drought_resilience": 0.8,
        "leaf_gas_exchange": 0.2,
    },
    "water_status_efficiency": {
        # iWUE uses nearest-neighbour — mechanistic prior not used.
        # Kept here for completeness; _gp_predict will bypass.
        "leaf_gas_exchange": 1.0,
    },
    # Stress resilience: drought_resilience primary (Jeong flavonoid, Hassanien sheath blight)
    # membrane_damage at 0.3 captures the Bochu 400 Hz membrane permeability data.
    "stress_resilience": {
        "drought_resilience": 0.7,
        "membrane_damage": 0.3,
    },
    "stress_resilience_disease": {
        "drought_resilience": 0.8,
        "leaf_gas_exchange": 0.2,
    },
    # Germination: seed_vigor is the dominant pathway at 200-500 Hz
    # drought_resilience has a small role (drought stress genes aid germination)
    "germination": {
        "seed_vigor": 0.9,
        "drought_resilience": 0.1,
    },
    # Yield: PAFT resonance at 400-700 Hz (Hou 2009, Qi 2010)
    # leaf_gas_exchange proxy captures the photosynthesis-yield link
    "yield": {
        "leaf_gas_exchange": 0.6,
        "seed_vigor": 0.4,
    },
}

# Fallback weights for outcomes not in the table above
_DEFAULT_PATHWAY_WEIGHTS = {
    "leaf_gas_exchange": 0.5,
    "seed_vigor": 0.3,
    "drought_resilience": 0.2,
}


def mechanistic_prior_for_outcome(
    outcome_name: str,
    frequency_hz: float,
    spl_db: float,
    stage_bucket: str,
    stress_context: str,
    medium: str,
    smech_air_value: float,
) -> float:
    """
    Compute the mechanistic prior mean k_prior(x) for a given outcome.

    Combines pathway activations weighted by OUTCOME_PATHWAY_WEIGHTS.
    Returns an effect percentage reflecting biological expectation at this
    frequency/SPL/stage/regime before any data correction.

    Falls back to 0.0 if no applicable pathways have activation.
    """
    weights = OUTCOME_PATHWAY_WEIGHTS.get(outcome_name, _DEFAULT_PATHWAY_WEIGHTS)
    total_activation = 0.0
    total_weight = 0.0

    for pathway_name, w in weights.items():
        if pathway_name not in PATHWAYS:
            continue
        pw = PATHWAYS[pathway_name]
        act = pathway_activation(
            pw, frequency_hz, spl_db, stage_bucket,
            stress_context, medium, smech_air_value,
        )
        total_activation += act * w * pw.peak_effect_pct
        total_weight += w

    if total_weight <= 0:
        return 0.0
    return total_activation / total_weight


# =============================================================================
# Empirical windows (retained for _match_windows traceability only)
# The mechanistic prior replaces windows as the GP prior mean.
# Windows are still reported in matched_windows for debugging.
# =============================================================================

@dataclass
class EmpiricalWindow:
    """Empirical evidence window — used only for traceability annotation in v4."""
    window_id: str
    freq_low_hz: float
    freq_high_hz: float
    spl_low_db: Optional[float]
    spl_high_db: Optional[float]
    stages: List[str]
    crop_tier: str
    direction: str
    effect_pct_low: float
    effect_pct_high: float
    confidence: float
    outcomes: List[str]
    source: str
    notes: str = ""


RICE_WINDOWS = [
    EmpiricalWindow(
        window_id="rice_photo_350_400",
        freq_low_hz=300, freq_high_hz=450, spl_low_db=60, spl_high_db=85,
        stages=["seedling", "vegetative"], crop_tier="rice", direction="increase",
        effect_pct_low=15, effect_pct_high=45, confidence=0.75,
        outcomes=["photosynthesis", "vigor"],
        source="Jusoh 2023 (350-380 Hz assimilation/height)",
    ),
    EmpiricalWindow(
        window_id="rice_drought_800_1500",
        freq_low_hz=800, freq_high_hz=1500, spl_low_db=60, spl_high_db=85,
        stages=["vegetative"], crop_tier="rice", direction="increase",
        effect_pct_low=5, effect_pct_high=50, confidence=0.65,
        outcomes=["water_status", "stress_resilience"],
        source="Jeong 2014 (250-1500 Hz stomata + RWC)",
        notes="10 full/high Jeong 2014 rows",
    ),
    EmpiricalWindow(
        window_id="rice_yield_paft_400_700",
        freq_low_hz=400, freq_high_hz=700, spl_low_db=65, spl_high_db=90,
        stages=["mixed", "vegetative", "reproductive"], crop_tier="rice",
        direction="increase", effect_pct_low=3, effect_pct_high=15, confidence=0.55,
        outcomes=["yield"],
        source="Hou 2009 (550 Hz PAFT yield +5.7%)",
    ),
    EmpiricalWindow(
        window_id="rice_germ_400_seed",
        freq_low_hz=350, freq_high_hz=450, spl_low_db=95, spl_high_db=110,
        stages=["seed", "germination"], crop_tier="rice", direction="increase",
        effect_pct_low=5, effect_pct_high=25, confidence=0.40,
        outcomes=["germination"],
        source="Bochu 2003 (400 Hz 106 dB germination)",
    ),
    EmpiricalWindow(
        window_id="rice_injury_4000_high_spl",
        freq_low_hz=3500, freq_high_hz=5000, spl_low_db=105, spl_high_db=130,
        stages=["seed", "germination", "seedling", "vegetative"],
        crop_tier="rice", direction="decrease",
        effect_pct_low=-30, effect_pct_high=-5, confidence=0.50,
        outcomes=["germination", "vigor", "photosynthesis"],
        source="Bochu 2003 (4000 Hz 111 dB injury)",
    ),
]

CROSS_SPECIES_WINDOWS = [
    EmpiricalWindow(
        window_id="grass_high_spl_inhibition",
        freq_low_hz=0, freq_high_hz=20000, spl_low_db=90, spl_high_db=130,
        stages=["vegetative", "seedling", "reproductive"],
        crop_tier="cross_species", direction="decrease",
        effect_pct_low=-50, effect_pct_high=-10, confidence=0.80,
        outcomes=["vigor", "photosynthesis", "yield"],
        source="Grass_General_90dB (>40% growth reduction at >90 dB)",
    ),
]

ALL_WINDOWS = RICE_WINDOWS + CROSS_SPECIES_WINDOWS


def _match_windows(frequency_hz, spl_db, stage_bucket, outcome_name,
                   stress_context="well_watered"):
    """Match empirical windows — for traceability annotation only in v4."""
    target_outcomes = OUTCOME_MAPPING.get(outcome_name, [outcome_name])
    matches = []
    for w in ALL_WINDOWS:
        if frequency_hz < w.freq_low_hz or frequency_hz > w.freq_high_hz:
            continue
        if w.spl_low_db is not None and spl_db > 0:
            if spl_db < w.spl_low_db or spl_db > w.spl_high_db:
                continue
        stage_match = stage_bucket in w.stages or "mixed" in w.stages
        if not stage_match and stage_bucket not in ["mixed", "any"]:
            continue
        outcome_match = (any(o in target_outcomes for o in w.outcomes) or
                         outcome_name in w.outcomes)
        if not outcome_match:
            continue
        quality = w.confidence
        if w.crop_tier != "rice":
            quality *= 0.6
        matches.append((w, quality))
    return matches


def _window_prior_mean(frequency_hz, spl_db, outcome_name, stage_bucket,
                        stress_context="well_watered"):
    """
    Compute window-based prior mean at a given frequency (fallback floor for v4).

    Used when the mechanistic prior is near-zero but direct empirical windows
    cover the frequency — e.g. rice_yield_paft_400_700 at 550 Hz for yield.
    This preserves direct literature evidence that the mechanistic pathway
    model may underweight due to pathway geometry.
    """
    matches = _match_windows(
        frequency_hz, spl_db, stage_bucket, outcome_name, stress_context
    )
    if not matches:
        return 0.0
    total_q = sum(q for _, q in matches)
    if total_q <= 0.0:
        return 0.0
    return sum(
        (w.effect_pct_low + w.effect_pct_high) / 2.0 * q
        for w, q in matches
    ) / total_q



# =============================================================================
# CSV data loading
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
    crop_type: str
    stage_bucket: str
    stress_context: str
    confidence_weight: str
    completeness_class: str
    usable_for_calibration: bool
    kpi_category: str
    notes: str


def _safe_float(val, default=0.0):
    try:
        return float(val)
    except (ValueError, TypeError):
        return default


def _infer_stress(water_regime):
    wr = (water_regime or "").lower()
    if "drought" in wr:
        return "drought"
    return "well_watered"


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
            if eff_size > 0 or calc_size > 0:
                direction = "increase"
            elif eff_size < 0 or calc_size < 0:
                direction = "decrease"
            else:
                direction = "unknown"
            spl_known = r.get("spl_db_known", "no").lower() == "yes"
            rows.append(SoundStudyRow(
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
            ))
    return rows


_ALL_DATA: Optional[List[SoundStudyRow]] = None


def get_all_data(data_dir=None):
    global _ALL_DATA
    if _ALL_DATA is None:
        _ALL_DATA = load_all_csv_data(data_dir)
    return _ALL_DATA


# =============================================================================
# KPI routing
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
                          "sheath_blight_incidence", "cell_membrane_permeability",
                          "H2O2_content"],
    "stress_resilience_disease": ["sheath_blight_incidence"],
    "yield": ["yield_kg_ha"],
}


# =============================================================================
# EffectEstimate public type (unchanged)
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


# =============================================================================
# Layer 1: Biophysical boundary check
# =============================================================================

def _check_biophysical_boundaries(frequency_hz, spl_db, stage_bucket="mixed",
                                   has_spl_matched_data=False):
    notes = []
    if spl_db > 0 and spl_db < ACTIVATION_FLOOR_DB:
        notes.append(f"SPL {spl_db:.0f} dB below {ACTIVATION_FLOOR_DB:.0f} dB activation floor")
        return 0.0, notes
    if spl_db > DAMAGE_CEILING_DB:
        notes.append(f"SPL {spl_db:.0f} dB above {DAMAGE_CEILING_DB:.0f} dB damage ceiling")
        return -1.0, notes
    if spl_db > GRASS_INHIBITION_DB and stage_bucket not in ("seed", "germination"):
        if not has_spl_matched_data:
            # Broadband high-SPL inhibition (grass membranes, Grass_General_90dB: -40%).
            # Spec formula: inhibition = -0.3 - 0.2 * frac (frac in [0,1] over 90-110 dB).
            # Interpretation: inhibition IS the fractional effect size (×100 = %).
            # At 90 dB: -0.30 → -30%; at 95 dB: -0.35 → -35%; at 110 dB: -0.50 → -50%.
            frac = (spl_db - GRASS_INHIBITION_DB) / (DAMAGE_CEILING_DB - GRASS_INHIBITION_DB)
            modifier = -0.3 - 0.2 * frac
            notes.append(
                f"Broadband high-SPL inhibition (grass membranes) applied. "
                f"SPL {spl_db:.0f} dB, modifier={modifier:.3f}"
            )
            return modifier, notes
    if AIRBORNE_CONFIDENCE_DECAY_START_HZ < frequency_hz < ULTRASOUND_FLOOR_HZ:
        if frequency_hz >= AIRBORNE_CONFIDENCE_DECAY_END_HZ:
            notes.append(f"Frequency {frequency_hz:.0f} Hz: no documented airborne effect")
            return 0.0, notes
        decay = 1.0 - ((frequency_hz - AIRBORNE_CONFIDENCE_DECAY_START_HZ) /
                        (AIRBORNE_CONFIDENCE_DECAY_END_HZ -
                         AIRBORNE_CONFIDENCE_DECAY_START_HZ))
        notes.append(f"High frequency {frequency_hz:.0f} Hz: confidence scaled {decay:.0%}")
        return decay, notes
    return 1.0, notes


# =============================================================================
# GP data collection
# =============================================================================

_CONF_MAP = {"high": 1.0, "medium": 0.6, "low": 0.3}
_COMP_MAP = {"full": 1.0, "partial": 0.5, "insufficient": 0.2}
_BASE_NOISE = 30.0
_MAX_NOISE = 500.0
_CROSS_SPECIES_NOISE_MULT = (1.0 / 0.6) ** 2

# _FLIP_SIGN_OUTCOMES: outcomes where a MORE NEGATIVE raw value is BENEFICIAL.
# Flipping these makes the GP see all beneficial effects as positive uniformly.
#
# Included:
#   germination_time, germination_time_h:  shorter time = faster = better  → flip
#   speed_of_germination: stored negative when faster in some studies      → flip
#
# EXPLICITLY EXCLUDED (and why):
#   sheath_blight_incidence: The CSV column calculated_effect_size_pct stores
#     +50 meaning "+50% disease REDUCTION" (already beneficial direction).
#     The reported_effect_size also says "+50% disease reduction".
#     Flipping would yield -50 (= harmful), which is wrong.
#     DESIGN DECISION: sheath_blight_incidence is pre-normalised in
#     rice_sound_master_v2.csv so that positive = beneficial. Do NOT add it
#     here. If the CSV changes convention (raw % incidence), this must be
#     revisited. Unit test: test_v4_architecture asserts +50% at 550 Hz.
#
#   cell_membrane_permeability: stored as -13.1 (mild damage at high SPL).
#     Negative = some membrane damage = slightly harmful. Correct as-is.
_FLIP_SIGN_OUTCOMES = {
    "germination_time", "germination_time_h", "speed_of_germination",
}

# _BOUNDARY_LAYER_STUDIES: Studies that encode biophysical BOUNDARY conditions
# rather than frequency-specific dose-response data. Excluded from GP training
# because they already shape the _check_biophysical_boundaries() rules:
#   GrassGeneral_90dB → justifies GRASS_INHIBITION_DB = 90 dB rule.
#     Including it as a GP point would distort the vigor kernel length scale
#     (its -40% at 95 Hz drags the kernel to narrow around the 350-400 Hz
#     cluster, underweighting evidence at 800-1500 Hz).
# See MODEL_NOTES.md "Data Sources" for the full rationale.
_BOUNDARY_LAYER_STUDIES = {"GrassGeneral_90dB"}


def _collect_gp_points(data, outcome_name, spl_db_query, stage_bucket,
                        stress_context, audible_only=True):
    """Collect (frequency_hz, effect_pct, noise_alpha) for GP fitting."""
    target_kpi_cats = KPI_CATEGORIES.get(outcome_name, [outcome_name.upper()])
    target_outcomes = OUTCOME_MAPPING.get(outcome_name, [outcome_name])
    points = []
    for row in data:
        if not row.usable_for_calibration:
            continue
        if row.frequency_hz <= 0:
            continue
        if audible_only and row.frequency_hz >= ULTRASOUND_FLOOR_HZ:
            continue
        matched = (
            (row.kpi_category and row.kpi_category in target_kpi_cats) or
            row.outcome_name in target_outcomes
        )
        if not matched:
            continue
        if row.study_id in _BOUNDARY_LAYER_STUDIES:
            continue
        sub_filter = SUB_OUTCOME_FILTERS.get(outcome_name)
        if sub_filter:
            if sub_filter["include"] and row.outcome_name not in sub_filter["include"]:
                continue
            if sub_filter["exclude"] and row.outcome_name in sub_filter["exclude"]:
                continue
        eff = row.calculated_effect_size_pct or row.reported_effect_size
        if eff == 0.0:
            continue
        if row.outcome_name in _FLIP_SIGN_OUTCOMES:
            eff = -eff
        conf_w = _CONF_MAP.get(row.confidence_weight, 0.3)
        comp_w = _COMP_MAP.get(row.completeness_class, 0.2)
        reliability = conf_w * comp_w
        spl_w = 1.0
        if row.spl_db_known and spl_db_query > 0 and row.spl_db > 0:
            spl_diff = abs(spl_db_query - row.spl_db)
            spl_w = max(0.2, 1.0 - spl_diff / 40.0)
        effective_rel = max(0.01, reliability * spl_w)
        species_mult = _CROSS_SPECIES_NOISE_MULT if row.crop_type != "rice" else 1.0
        alpha = min(_MAX_NOISE, _BASE_NOISE * (1.0 / effective_rel) * species_mult)
        points.append((row.frequency_hz, eff, alpha))
    return points


# =============================================================================
# Frequency aggregation (v3.1+)
# =============================================================================

def _aggregate_same_frequency(points, tol_hz=1.0):
    """Collapse pseudo-replication: precision-weighted merge of same-freq points."""
    if not points:
        return []
    sorted_pts = sorted(points, key=lambda p: p[0])
    aggregated = []
    group = [sorted_pts[0]]
    for pt in sorted_pts[1:]:
        if abs(pt[0] - group[0][0]) <= tol_hz:
            group.append(pt)
        else:
            aggregated.append(_merge_freq_group(group))
            group = [pt]
    aggregated.append(_merge_freq_group(group))
    return aggregated


def _merge_freq_group(group):
    if len(group) == 1:
        return group[0]
    precisions = [1.0 / max(p[2], 1e-6) for p in group]
    total_prec = sum(precisions)
    mean_freq = sum(p[0] * pr for p, pr in zip(group, precisions)) / total_prec
    mean_eff = sum(p[1] * pr for p, pr in zip(group, precisions)) / total_prec
    combined_alpha = max(1.0, 1.0 / total_prec)
    return (mean_freq, mean_eff, combined_alpha)


# =============================================================================
# Ultrasound lookup
# =============================================================================

def _ultrasound_lookup(frequency_hz, outcome_name, data):
    """Distance-weighted lookup for ultrasound (>=20 kHz)."""
    target_kpi_cats = KPI_CATEGORIES.get(outcome_name, [outcome_name.upper()])
    target_outcomes = OUTCOME_MAPPING.get(outcome_name, [outcome_name])
    _US_GERM = {"germination_time", "germination_time_h", "germination_pct",
                "germination_speed", "germination_potential", "germination_energy",
                "germination_index", "speed_of_germination", "water_uptake_rate"}
    _US_VIGOR = {"seedling_dry_weight", "seedling_dry_weight_g", "plant_height",
                 "root_length", "fresh_weight", "shoot_dry_weight"}
    us_pts = []
    for row in data:
        if row.frequency_hz < ULTRASOUND_FLOOR_HZ:
            continue
        matched = (
            (row.kpi_category and row.kpi_category in target_kpi_cats) or
            row.outcome_name in target_outcomes
        )
        if not matched:
            continue
        if outcome_name == "germination" and row.outcome_name not in _US_GERM:
            continue
        if outcome_name == "vigor" and row.outcome_name not in _US_VIGOR:
            continue
        eff = row.calculated_effect_size_pct or row.reported_effect_size
        if eff == 0.0:
            continue
        if row.outcome_name in _FLIP_SIGN_OUTCOMES:
            eff = -eff
        conf_w = _CONF_MAP.get(row.confidence_weight, 0.3)
        comp_w = _COMP_MAP.get(row.completeness_class, 0.2)
        base_conf = conf_w * comp_w
        if row.crop_type != "rice":
            base_conf *= 0.6
        us_pts.append((row.frequency_hz, eff, base_conf))
    if not us_pts:
        return 0.0, 0.0, 0, ["No ultrasound data for this outcome"]
    total_w, total_e = 0.0, 0.0
    for f, eff, w in us_pts:
        prox = max(0.05, 1.0 - abs(frequency_hz - f) / 30000.0)
        total_w += w * prox
        total_e += eff * w * prox
    effect = total_e / total_w if total_w > 0 else 0.0
    conf = min(1.0, total_w / len(us_pts)) if total_w > 0 else 0.0
    return effect, conf, len(us_pts), [f"Ultrasound nearest-neighbour ({len(us_pts)} pts)"]


# =============================================================================
# Layer 4: Gaussian Process estimator (v4 — uses mechanistic prior)
# =============================================================================

def _gp_predict(frequency_hz, spl_db, outcome_name, stage_bucket,
                 stress_context, data_dir):
    """
    GP prediction using mechanistic prior mean (v4).

    For ultrasound (>=20 kHz): delegates to ultrasound_lookup.
    For iWUE: uses nearest-neighbour (standing-wave oscillation).
    For all other audible outcomes:
      1. Compute S_mech,air(f, SPL) via the MechActivation GP.
      2. Compute mechanistic prior for the outcome.
      3. Fit GP to residuals (observed - mechanistic_prior).
      4. Return posterior mean + uncertainty-based confidence.
    """
    data = get_all_data(data_dir)

    # --- Ultrasound ---
    if frequency_hz >= ULTRASOUND_FLOOR_HZ:
        ef, co, np_, no = _ultrasound_lookup(frequency_hz, outcome_name, data)
        return ef, co, np_, no, True  # ultrasound: SPL matching is N/A

    raw_points = _collect_gp_points(
        data, outcome_name, spl_db, stage_bucket, stress_context, audible_only=True
    )
    if outcome_name == "water_status_efficiency":
        # Do NOT aggregate iWUE frequencies. The Jusoh 2023 iWUE data
        # oscillates ±60% over just 30 Hz (350→353→357→359→380 Hz) due to
        # standing-wave interference artefacts in the growth chamber. This
        # oscillation is REAL SIGNAL for the nearest-neighbour lookup; merging
        # 357 and 359 Hz would collapse the peak and lose the ±62.7% at 359 Hz.
        points = raw_points
    else:
        points = _aggregate_same_frequency(raw_points, tol_hz=1.0)
    n = len(points)
    n_raw = len(raw_points)

    # SPL-specific match: True only when CSV data exists NEAR THIS FREQUENCY
    # AND THIS SPL (within ±14 dB). Frequency proximity ±150 Hz ensures that
    # studies at different frequency bands don't bypass the inhibitory zone for
    # unrelated frequencies (e.g. Bochu at 200 Hz/106 dB does not count for a
    # query at 0 Hz/95 dB; but Jeong at 800 Hz/100 dB DOES count at 800 Hz).
    _spl_threshold_w = 0.65   # spl_w > 0.65 ↔ SPL within ~14 dB of query
    _freq_prox_hz   = 150.0   # frequency proximity window for SPL match
    spl_matched = (
        n > 0 and spl_db > 0 and
        any(
            row.spl_db_known and row.spl_db > 0 and
            (1.0 - abs(spl_db - row.spl_db) / 40.0) > _spl_threshold_w and
            abs(row.frequency_hz - frequency_hz) <= _freq_prox_hz
            for row in data
            if row.usable_for_calibration and row.frequency_hz < ULTRASOUND_FLOOR_HZ
            and (row.kpi_category in KPI_CATEGORIES.get(outcome_name, [outcome_name.upper()])
                 or row.outcome_name in OUTCOME_MAPPING.get(outcome_name, [outcome_name]))
        )
    )

    # --- Mechanistic prior ---
    mech_state = _get_mech_state(data, data_dir or "")
    smech = smech_air(frequency_hz, spl_db, mech_state)
    medium = "airborne"
    mech_prior = mechanistic_prior_for_outcome(
        outcome_name, frequency_hz, spl_db, stage_bucket, stress_context,
        medium, smech,
    )
    # Window floor: if mechanistic prior is near-zero but direct empirical
    # windows cover this frequency, use the window midpoint as a prior floor.
    # This preserves literature evidence that pathway geometry may underweight
    # (e.g. yield at 550 Hz from Hou 2009 / rice_yield_paft_400_700 window).
    window_prior = _window_prior_mean(
        frequency_hz, spl_db, outcome_name, stage_bucket, stress_context
    )
    if abs(mech_prior) < 2.0 and abs(window_prior) >= 2.0:
        prior_at_query = window_prior  # direct evidence floor
    else:
        prior_at_query = mech_prior

    # --- iWUE nearest-neighbour ---
    if outcome_name == "water_status_efficiency":
        if not points:
            return prior_at_query, 0.1, 0, ["iWUE: no data, returning mechanistic prior"], spl_matched
        best_dist, best_eff, best_alpha = float("inf"), 0.0, _MAX_NOISE
        for f, eff, alpha in points:
            d = abs(frequency_hz - f)
            if d < best_dist:
                best_dist, best_eff, best_alpha = d, eff, alpha
        decay = max(0.05, 1.0 - best_dist / 50.0)
        base_conf = min(1.0, _BASE_NOISE / best_alpha)
        conf = base_conf * decay * (0.3 + 0.7 * min(1.0, n / 4.0))
        return (best_eff, conf, n,
                [f"iWUE nearest-neighbour: {best_eff:+.1f}% at dist={best_dist:.0f} Hz"],
                spl_matched)

    if n == 0:
        window_matches = _match_windows(
            frequency_hz, spl_db, stage_bucket, outcome_name, stress_context
        )
        matched_ids = [w.window_id for w, _ in window_matches]
        conf = 0.15 if window_matches or abs(prior_at_query) > 1.0 else 0.0
        return (
            prior_at_query, conf, 0,
            [f"Mechanistic prior only (no CSV data); S_mech={smech:.3f}; "
             f"windows: {', '.join(matched_ids) or 'none'}"],
            spl_matched,
        )

    freqs = np.array([p[0] for p in points])
    effects = np.array([p[1] for p in points])
    alphas = np.array([p[2] for p in points])

    def _hybrid_prior(f_hz):
        """Hybrid prior: mechanistic if strong, else window floor."""
        m = mechanistic_prior_for_outcome(
            outcome_name, f_hz, spl_db, stage_bucket, stress_context, medium,
            smech_air(f_hz, spl_db, mech_state),
        )
        if abs(m) < 2.0:
            w = _window_prior_mean(f_hz, spl_db, outcome_name, stage_bucket, stress_context)
            if abs(w) >= 2.0:
                return w
        return m

    # Residuals relative to HYBRID prior (same function used for query and training)
    prior_at_train = np.array([_hybrid_prior(f) for f in freqs])
    residuals = effects - prior_at_train

    if not _SKLEARN_AVAILABLE:
        dists = np.abs(freqs - frequency_hz)
        weights = 1.0 / (dists + 50.0)
        weights /= weights.sum()
        effect = float(prior_at_query + np.dot(weights, residuals))
        conf = min(0.4, n / 8.0)
        return effect, conf, n, ["sklearn unavailable -- IDW fallback used"], spl_matched

    X_train = freqs.reshape(-1, 1)
    amplitude_init = max(float(np.mean(residuals ** 2)), 1.0)

    if n < 3:
        kernel = (ConstantKernel(amplitude_init, constant_value_bounds="fixed")
                  * Matern(200.0, length_scale_bounds="fixed", nu=2.5))
        n_restarts = 0
    else:
        kernel = (
            ConstantKernel(amplitude_init,
                           constant_value_bounds=(amplitude_init * 0.01,
                                                  amplitude_init * 100.0))
            * Matern(200.0, length_scale_bounds=(100.0, 1200.0), nu=2.5)
        )
        n_restarts = 5

    gp = GaussianProcessRegressor(
        kernel=kernel, alpha=alphas, normalize_y=False,
        n_restarts_optimizer=n_restarts, random_state=42,
    )
    try:
        gp.fit(X_train, residuals)
    except Exception as exc:
        return (prior_at_query, 0.1, n,
                [f"GP fit failed ({exc}); returning mechanistic prior"],
                spl_matched)

    residual_pred, residual_std = gp.predict(np.array([[frequency_hz]]), return_std=True)
    effect = prior_at_query + float(residual_pred[0])
    std = float(residual_std[0])

    effect_scale = max(abs(effect), 15.0)
    proximity_conf = float(np.exp(-0.5 * (std / effect_scale) ** 2))

    if n >= 3:
        l_scale = float(gp.kernel_.k2.length_scale)
    else:
        l_scale = 200.0
    n_local = sum(1 for f in freqs if abs(f - frequency_hz) <= l_scale)
    n_factor = min(1.0, n_local / 2.0)

    confidence = proximity_conf * (0.3 + 0.7 * n_factor)
    confidence = max(0.0, min(1.0, confidence))

    kernel_desc = str(gp.kernel_) if n >= 3 else "Matern(l=200 Hz) [fixed, n<3]"
    notes = [
        f"GP v4: mean={effect:.1f}%, std={std:.1f}%, n_clusters={n} (raw={n_raw})",
        f"S_mech={smech:.3f}, prior={prior_at_query:.1f}%, l={l_scale:.0f} Hz, n_local={n_local}",
        f"Kernel: {kernel_desc}",
    ]
    return effect, confidence, n, notes, spl_matched


# =============================================================================
# Public API: estimate_effect (unchanged interface)
# =============================================================================

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
    Estimate sound treatment effect on a rice outcome.

    v4 architecture:
      1. Biophysical boundaries (hard physics)
      2. S_mech,air + mechanistic pathway prior
      3. GP residual correction on CSV data
      4. Dosage scaling
    """
    if hours_per_day <= 0 or (spl_db <= 0 and frequency_hz < ULTRASOUND_FLOOR_HZ):
        return EffectEstimate(0.0, "none", 0.0, 0, [], "No treatment applied.")

    gp_effect, gp_confidence, n_points, gp_notes, spl_matched = _gp_predict(
        frequency_hz, spl_db, outcome_name, stage_bucket, stress_context, data_dir
    )

    # has_spl_matched_data passed from _gp_predict (returned as 4th element)
    has_spl_matched_data = n_points > 0 and spl_matched
    bio_modifier, bio_notes = _check_biophysical_boundaries(
        frequency_hz, spl_db, stage_bucket,
        has_spl_matched_data=has_spl_matched_data,
    )

    if bio_modifier == 0.0:
        return EffectEstimate(0.0, "none", 0.0, 0, [],
                              "; ".join(bio_notes) or "Below activation threshold.")

    if bio_modifier < 0:
        if spl_db > DAMAGE_CEILING_DB:
            # Hard damage ceiling: membrane rupture, fixed -15% override
            return EffectEstimate(
                effect_pct=-15.0, effect_direction="decrease", confidence=0.45,
                supporting_row_count=0, matched_windows=["biophysical_damage_ceiling"],
                notes="; ".join(bio_notes) + " — Hard biophysical override: cellular damage.",
            )
        else:
            # High-SPL inhibitory zone (90-110 dB, no matched data).
            # bio_modifier is the spec fractional value (e.g. -0.35 at 95 dB).
            # Multiply by 100 to get percentage: -0.35 → -35%.
            # Ref: Grass_General_90dB (-40% vegetative growth at >90 dB continuous).
            inhibition_pct = bio_modifier * 100.0  # -0.35 → -35%
            inhibition_pct = max(-50.0, inhibition_pct)  # floor at -50%
            confidence = 0.35  # moderate: cross-species grass study
            return EffectEstimate(
                effect_pct=round(inhibition_pct, 2), effect_direction="decrease",
                confidence=confidence, supporting_row_count=0,
                matched_windows=["biophysical_high_spl_inhibition"],
                notes="; ".join(bio_notes),
            )

    if n_points == 0 and gp_confidence == 0.0:
        return EffectEstimate(0.0, "none", 0.0, 0, [],
                              "No matching pathways, windows, or CSV data.")

    final_effect = gp_effect * bio_modifier

    if hours_per_day > 0 and frequency_hz < ULTRASOUND_FLOOR_HZ:
        if hours_per_day < 1.0:
            dose_factor = math.sqrt(hours_per_day)
        elif hours_per_day <= 4.0:
            dose_factor = 1.0
        else:
            dose_factor = min(1.3, 1.0 + 0.1 * math.log(hours_per_day / 4.0))
        final_effect *= dose_factor

    if abs(final_effect) > 100:
        gp_confidence *= 0.5
        final_effect = max(-100.0, min(100.0, final_effect))
        gp_notes.append("Effect capped at +-100%")

    direction = ("increase" if final_effect > 0.5
                 else "decrease" if final_effect < -0.5 else "none")

    window_matches = _match_windows(
        frequency_hz, spl_db, stage_bucket, outcome_name, stress_context
    )
    matched_window_ids = [w.window_id for w, _ in window_matches]
    all_notes = bio_notes + gp_notes
    if matched_window_ids:
        all_notes.append(f"Windows(ref): {', '.join(matched_window_ids)}")

    return EffectEstimate(
        effect_pct=round(final_effect, 2),
        effect_direction=direction,
        confidence=round(max(0.0, min(1.0, gp_confidence)), 3),
        supporting_row_count=n_points,
        matched_windows=matched_window_ids,
        notes="; ".join(all_notes),
    )


# =============================================================================
# Public API: estimate_all_effects
# =============================================================================

def estimate_all_effects(frequency_hz, spl_db, hours_per_day,
                          stage_bucket="mixed", stress_context="well_watered",
                          data_dir=None):
    """Estimate effects for all supported outcomes including water sub-KPIs."""
    results = {}
    for outcome in ["germination", "vigor", "photosynthesis", "yield"]:
        results[outcome] = estimate_effect(
            outcome, frequency_hz, spl_db, hours_per_day,
            stage_bucket, stress_context, data_dir,
        )

    ws_gsw_d = estimate_effect("water_status_gsw_drought", frequency_hz, spl_db,
                                hours_per_day, stage_bucket, stress_context, data_dir)
    ws_gsw_w = estimate_effect("water_status_gsw_wellwatered", frequency_hz, spl_db,
                                hours_per_day, stage_bucket, stress_context, data_dir)
    ws_rwc = estimate_effect("water_status_rwc", frequency_hz, spl_db,
                              hours_per_day, stage_bucket, stress_context, data_dir)
    ws_eff = estimate_effect("water_status_efficiency", frequency_hz, spl_db,
                              hours_per_day, stage_bucket, stress_context, data_dir)

    results["water_status_gsw_drought"] = ws_gsw_d
    results["water_status_gsw_wellwatered"] = ws_gsw_w
    results["water_status_rwc"] = ws_rwc
    results["water_status_efficiency"] = ws_eff
    results["water_status_openness"] = (
        ws_gsw_d if stress_context == "drought" else ws_gsw_w
    )

    if stress_context == "drought":
        components = [(ws_gsw_d, 0.6), (ws_rwc, 0.25), (ws_eff, 0.15)]
    else:
        components = [(ws_gsw_w, 0.4), (ws_eff, 0.4), (ws_rwc, 0.2)]

    agg_pct = sum(c.effect_pct * w for c, w in components)
    agg_conf = sum(c.confidence * w for c, w in components)
    agg_count = sum(c.supporting_row_count for c, _ in components)

    if any(c.effect_pct != 0 for c, _ in components):
        agg_dir = ("increase" if agg_pct > 0.5 else "decrease" if agg_pct < -0.5 else "none")
        results["water_status"] = EffectEstimate(
            round(agg_pct, 2), agg_dir, round(agg_conf, 3), agg_count, [],
            f"CWSI-weighted: {'gsw_drought primary' if stress_context == 'drought' else 'gsw_wellwatered+iWUE primary'}",
        )
    else:
        results["water_status"] = EffectEstimate(0, "none", 0, 0, [], "No water status data")

    sr_disease = estimate_effect("stress_resilience_disease", frequency_hz, spl_db,
                                  hours_per_day, stage_bucket, stress_context, data_dir)
    sr_general = estimate_effect("stress_resilience", frequency_hz, spl_db,
                                  hours_per_day, stage_bucket, stress_context, data_dir)
    results["stress_resilience_disease"] = sr_disease
    results["stress_resilience"] = sr_general

    return results


# =============================================================================
# Smoke test
# =============================================================================

if __name__ == "__main__":
    data_dir = os.path.join(os.path.dirname(__file__), "..", "data")
    print("=== Sound Response Engine v4 (Mechanistic GP) ===")
    print(f"scikit-learn: {_SKLEARN_AVAILABLE}")
    data = get_all_data(data_dir)
    mech = train_mechactivation_gp(data, data_dir)
    print(f"MechActivation GP: available={mech.available}, n={mech.n_training}, y_max={mech.y_max:.1f}%")
    for f in [300, 350, 550, 800, 1000, 1500]:
        sf = mechactivation_Sf(f, mech)
        sm = smech_air(f, 75.0, mech)
        print(f"  {f} Hz: S_f={sf:.3f}  S_mech(75dB)={sm:.3f}")
    print()
    for freq in [350, 550, 800, 1000]:
        effects = estimate_all_effects(freq, 75, 8.0, "mixed", "well_watered", data_dir)
        print(f"{freq} Hz, 75 dB:")
        for name in ["photosynthesis", "vigor", "water_status", "stress_resilience"]:
            e = effects[name]
            print(f"  {name}: {e.effect_pct:+.1f}% conf={e.confidence:.3f}")
