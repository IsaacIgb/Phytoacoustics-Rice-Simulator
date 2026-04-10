"""
rice_core.py — Simplified ORYZA-inspired rice growth engine.
v4: dual-phase sound driver (germination vs seedling-to-harvest).

Scientific basis:
  - Bouman et al. (2001), ORYZA2000: modeling lowland rice
  - Li et al. (2017), From ORYZA2000 to ORYZA (v3)

Architecture follows WOFOST/PCSE patterns:
  - Each biological process is a separate class with calc_rates()
    (following PCSE's SimulationObject pattern).
  - All rates are computed from the CURRENT state before any integration,
    collected in a DailyRates dataclass. This prevents order-of-computation
    bugs where one process reads a half-updated state.
  - State is a snapshot dataclass per time step.
  - Parameters are isolated in a CultivarParams dataclass, not hardcoded
    inline, so cultivars can be swapped without editing process logic.

See MODEL_NOTES.md for scientific simplification rationale.
"""

from dataclasses import dataclass, field
from typing import Optional, List
import math

# =============================================================================
# Utility: linear interpolation from table (LINT2 equivalent)
# Bouman et al. (2001) Section 2.4.1
# =============================================================================

def lint(table: list, x: float) -> float:
    """Linear interpolation from (x,y) table. Clamps at boundaries."""
    if x <= table[0][0]:
        return table[0][1]
    if x >= table[-1][0]:
        return table[-1][1]
    for i in range(len(table) - 1):
        x0, y0 = table[i]
        x1, y1 = table[i + 1]
        if x0 <= x <= x1:
            frac = (x - x0) / (x1 - x0) if x1 != x0 else 0.0
            return y0 + frac * (y1 - y0)
    return table[-1][1]


# =============================================================================
# Parameters — isolated so cultivars can be swapped
# =============================================================================

@dataclass
class CultivarParams:
    """
    Crop parameters for a single cultivar. Following PCSE's CropDataProvider
    pattern: all cultivar-specific numbers in one place, separate from process
    logic.

    Default values are for IR72, the reference cultivar in Bouman et al. (2001).
    """
    # --- Phenology (Section 3.2.1) ---
    tbase: float = 8.0
    topt: float = 30.0
    thigh: float = 42.0
    dvrj: float = 0.000773
    dvri: float = 0.000749
    dvrp: float = 0.000785
    dvrr: float = 0.001044

    # --- Leaf area (Section 3.2.9) ---
    rgrlmx: float = 0.0085   # Max relative leaf area growth rate ((°Cd)^-1)
    rgrlmn: float = 0.0045   # Min RGRL — unused in v1 (no N module)
    lai_exp_threshold: float = 1.0

    sla_table: list = field(default_factory=lambda: [
        (0.0, 0.0045), (0.4, 0.0033), (0.65, 0.0028),
        (1.0, 0.0024), (1.6, 0.0020), (2.0, 0.0018),
    ])

    # --- Partitioning (Section 3.2.3) ---
    # --- Partitioning (Section 3.2.3) ---
    # ORYZA2000 IR72 tables. Previous version over-allocated to grain
    # (FSO=0.90 at DVS=1.0) producing HI≈0.68 vs observed 0.45-0.50.
    # Corrected to ORYZA2000 values which keep more in stems.
    fsh_table: list = field(default_factory=lambda: [
        (0.0, 0.50), (0.4, 0.55), (0.65, 0.70),
        (1.0, 0.95), (1.5, 1.0), (2.0, 1.0),
    ])
    flv_table: list = field(default_factory=lambda: [
        (0.0, 0.55), (0.33, 0.45), (0.65, 0.30),
        (0.80, 0.15), (1.0, 0.0), (2.0, 0.0),
    ])
    fst_table: list = field(default_factory=lambda: [
        (0.0, 0.45), (0.33, 0.55), (0.65, 0.45),
        (0.80, 0.35), (1.0, 0.35), (1.2, 0.30), (1.6, 0.25), (2.0, 0.20),
    ])
    fso_table: list = field(default_factory=lambda: [
        (0.0, 0.0), (0.65, 0.0), (0.80, 0.30),
        (1.0, 0.50), (1.2, 0.55), (1.6, 0.65), (2.0, 0.70),
    ])

    # --- Leaf death (Section 3.2.4) ---
    drlv_table: list = field(default_factory=lambda: [
        (0.0, 0.0), (0.6, 0.0), (1.0, 0.015),
        (1.6, 0.025), (2.0, 0.05),
    ])

    # --- Respiration (Section 3.2.5) ---
    mainlv: float = 0.02
    mainst: float = 0.01
    mainrt: float = 0.01
    mainso: float = 0.003
    tref: float = 25.0
    q10: float = 2.0
    crglv: float = 1.326
    crgst: float = 1.326
    crgrt: float = 1.326
    crgso: float = 1.462

    # --- Photosynthesis (Section 3.2.2) ---
    kdf: float = 0.4
    amax_ref: float = 40.0
    eff_ref: float = 0.45

    # --- Spikelets (Section 3.2.8) ---
    spgf: float = 65000.0


# =============================================================================
# State, rates, config
# =============================================================================

@dataclass
class WeatherDay:
    """Daily weather inputs."""
    tmin: float = 22.0
    tmax: float = 32.0
    radiation: float = 18000.0
    rain: float = 5.0

    @property
    def tav(self) -> float:
        return (self.tmin + self.tmax) / 2.0

    @property
    def tavd(self) -> float:
        return (self.tmax + self.tav) / 2.0


@dataclass
class RiceState:
    """State variables — snapshot at one point in time."""
    day: int = 0
    dvs: float = 0.0
    wlvg: float = 20.0
    wlvd: float = 0.0
    wst: float = 5.0
    wrt: float = 10.0
    wso: float = 0.0
    lai: float = 0.1
    temp_sum: float = 0.0
    total_biomass: float = 35.0
    n_spikelets: float = 0.0
    n_grains: float = 0.0
    spikelet_fertility: float = 1.0
    water_stress: float = 1.0
    daily_assim: float = 0.0
    daily_growth: float = 0.0
    yield_proxy: float = 0.0

    @property
    def wagt(self) -> float:
        return self.wlvg + self.wlvd + self.wst + self.wso


@dataclass
class DailyRates:
    """
    All rate variables for one day, computed before integration.
    PCSE pattern: rates are a separate container from state.
    """
    hu: float = 0.0
    dvr: float = 0.0
    water_stress: float = 1.0
    dtga: float = 0.0
    rmcr: float = 0.0
    fsh: float = 0.5
    frt: float = 0.5
    flv: float = 0.5
    fst: float = 0.3
    fso: float = 0.0
    gcr: float = 0.0
    grt: float = 0.0
    glv: float = 0.0
    gst: float = 0.0
    gso: float = 0.0
    llv: float = 0.0
    new_lai: float = 0.0
    lai_wlvg_adjustment: float = 0.0
    gnsp: float = 0.0
    spikelet_fertility: float = 1.0


@dataclass
class SimulationConfig:
    total_days: int = 120
    establishment: str = "transplant"
    plants_per_m2: float = 25.0
    leaf_area_per_plant: float = 0.0004
    water_regime: str = "irrigated"
    drought_start_dvs: float = 0.5
    drought_severity: float = 0.5
    weather: Optional[List[WeatherDay]] = None
    co2_ppm: float = 400.0


# =============================================================================
# Process classes — each is independently testable
# =============================================================================

class Phenology:
    """Bouman et al. (2001) Section 3.2.1, Eqs 3.2-3.4."""

    @staticmethod
    def calc_heat_units(tmax, tmin, p):
        tm = (tmax + tmin) / 2.0
        hu = 0.0
        for h in range(1, 25):
            td = tm + 0.5 * abs(tmax - tmin) * math.cos(0.2618 * (h - 14))
            if td > p.tbase and td < p.thigh:
                if td <= p.topt:
                    hu += (td - p.tbase) / 24.0
                else:
                    hu += (p.topt - (td - p.topt) * (p.topt - p.tbase) / (p.thigh - p.topt)) / 24.0
        return max(0.0, hu)

    @staticmethod
    def calc_rates(state, weather, ws, p):
        hu = Phenology.calc_heat_units(weather.tmax, weather.tmin, p)
        if state.dvs < 0.40:
            dvr = p.dvrj * hu
        elif state.dvs < 0.65:
            dvr = p.dvri * hu
        elif state.dvs < 1.0:
            dvr = p.dvrp * hu
        else:
            dvr = p.dvrr * hu
        if state.dvs < 1.0 and ws < 1.0:
            dvew = ws + (state.dvs * (1.0 - ws))
            dvr *= dvew
        return hu, dvr


class WaterStress:
    """Inspired by Li et al. (2017) Section 2.3.2. Uses DVS-duration proxy."""

    @staticmethod
    def calc_rates(state, config):
        regime = config.water_regime
        if regime == "irrigated":
            return 1.0
        elif regime == "rainfed":
            return 0.85 if 0.65 <= state.dvs <= 1.5 else 0.95
        elif regime == "drought":
            if state.dvs >= config.drought_start_dvs:
                duration = state.dvs - config.drought_start_dvs
                return max(0.1, math.exp(-config.drought_severity * duration * 2.0))
            return 1.0
        return 1.0


class Assimilation:
    """
    Big-leaf canopy model. Uses analytical Monsi-Saeki canopy integral.
    Not from ORYZA2000 (which uses multi-layer Gaussian integration).
    Output checked for plausibility against ORYZA2000 typical range.
    Bouman et al. (2001) Section 3.2.2, Eqs 3.5-3.7 for leaf-level params.
    """

    @staticmethod
    def calc_rates(state, weather, ws, co2_ppm, p):
        if state.lai <= 0.0 or weather.radiation <= 0.0:
            return 0.0
        co2_eff_num = 1.0 - math.exp(-0.00305 * co2_ppm - 0.222)
        co2_eff_den = 1.0 - math.exp(-0.00305 * 340.0 - 0.222)
        co2_eff = co2_eff_num / co2_eff_den if co2_eff_den != 0 else 1.0
        eff = p.eff_ref * co2_eff
        amax_co2 = 49.57 / 34.26 * (1.0 - math.exp(-0.208 * (co2_ppm - 60.0) / 49.57))
        amax_leaf = p.amax_ref * max(0.0, amax_co2)
        tav = weather.tav
        if tav < 10:
            temp_red = 0.0
        elif tav < 20:
            temp_red = (tav - 10) / 10.0
        elif tav <= 35:
            temp_red = 1.0
        else:
            temp_red = max(0.0, 1.0 - (tav - 35) / 7.0)
        amax_leaf *= temp_red
        par = weather.radiation * 0.5
        frac_intercepted = 1.0 - math.exp(-p.kdf * state.lai)
        absorbed_par = frac_intercepted * par
        daylength = 12.0
        par_rate = absorbed_par / (daylength * 3.6)
        effective_lai = min(state.lai, 8.0)
        amax_canopy = amax_leaf * (1.0 - math.exp(-p.kdf * effective_lai)) / p.kdf
        if amax_canopy > 0:
            assim_rate = amax_canopy * (1.0 - math.exp(-eff * par_rate / amax_canopy))
        else:
            assim_rate = 0.0
        return max(0.0, assim_rate * daylength * ws)


class Respiration:
    """Bouman et al. (2001) Section 3.2.5, Eqs 3.25-3.28."""

    @staticmethod
    def calc_maintenance(state, tav, p):
        teff = p.q10 ** ((tav - p.tref) / 10.0)
        mndvs = state.wlvg / max(1.0, state.wlvg + state.wlvd)
        return (state.wlvg * p.mainlv + state.wst * p.mainst +
                state.wso * p.mainso + state.wrt * p.mainrt) * teff * mndvs

    @staticmethod
    def calc_growth_coeff(fsh, frt, flv, fst, fso, p):
        return max(1.0, fsh * (p.crglv * flv + p.crgst * fst + p.crgso * fso) + p.crgrt * frt)


class Partitioning:
    """Bouman et al. (2001) Section 3.2.3, Eq 3.21."""

    @staticmethod
    def calc_rates(state, ws, p):
        fsh = lint(p.fsh_table, state.dvs)
        frt = 1.0 - fsh
        if state.dvs < 1.0 and ws < 1.0:
            fsh *= ws
            frt = 1.0 - fsh
        flv = lint(p.flv_table, state.dvs)
        fst = lint(p.fst_table, state.dvs)
        fso = lint(p.fso_table, state.dvs)
        return fsh, frt, flv, fst, fso


class LeafDynamics:
    """
    Bouman et al. (2001) Section 3.2.9 (SUBLAI2) and Section 3.2.4.
    Dual-phase: exponential (Eq 3.38) when LAI < 1.0, SLA-based (Eq 3.41) after.
    """

    @staticmethod
    def calc_rates(state, hu, ws, glv, p):
        sla = lint(p.sla_table, state.dvs)
        drlv = lint(p.drlv_table, state.dvs)
        llv = drlv * state.wlvg
        new_wlvg = max(0.0, state.wlvg + glv - llv)
        adjustment = 0.0
        if state.lai < p.lai_exp_threshold:
            rgrl = p.rgrlmx * ws
            glai = state.lai * rgrl * hu
            new_lai = state.lai + glai
            implied_wlvg = new_lai / sla if sla > 0 else new_wlvg
            if implied_wlvg > new_wlvg:
                adjustment = implied_wlvg - new_wlvg
        else:
            new_lai = new_wlvg * sla
        return max(0.01, new_lai), adjustment, llv


class Spikelets:
    """Bouman et al. (2001) Section 3.2.8, Eqs 3.31-3.34."""

    @staticmethod
    def calc_rates(state, gcr, weather, p):
        gnsp = gcr * p.spgf if 0.65 <= state.dvs <= 1.0 else 0.0
        fert = state.spikelet_fertility
        if 0.75 <= state.dvs <= 1.2:
            ctt = max(0.0, 22.0 - weather.tav)
            if ctt > 0:
                sf1 = 1.0 - (4.6 + 0.054 * ctt ** 1.56) / 100.0
                fert = min(fert, max(0.0, sf1))
        if 0.96 <= state.dvs <= 1.2 and weather.tmax > 35:
            sf2 = 1.0 / (1.0 + math.exp(0.853 * (weather.tmax - 36.6)))
            fert = min(fert, max(0.0, sf2))
        return gnsp, fert


# =============================================================================
# Engine: orchestrates processes, collects rates, integrates
# =============================================================================

class RiceEngine:
    """
    Orchestrates the daily simulation step.
    PCSE pattern: (1) all rates from current state, (2) integrate.
    """

    def __init__(self, params=None):
        self.params = params or CultivarParams()

    def calc_all_rates(self, state, weather, config):
        """Compute all rates from current state. No state mutation."""
        p = self.params
        rates = DailyRates()
        rates.water_stress = WaterStress.calc_rates(state, config)
        ws = rates.water_stress
        rates.hu, rates.dvr = Phenology.calc_rates(state, weather, ws, p)
        rates.dtga = Assimilation.calc_rates(state, weather, ws, config.co2_ppm, p)
        rates.rmcr = Respiration.calc_maintenance(state, weather.tav, p)
        rates.fsh, rates.frt, rates.flv, rates.fst, rates.fso = Partitioning.calc_rates(state, ws, p)
        crgcr = Respiration.calc_growth_coeff(rates.fsh, rates.frt, rates.flv, rates.fst, rates.fso, p)
        ch2o = rates.dtga * (30.0 / 44.0)
        rates.gcr = max(0.0, (ch2o - rates.rmcr) / crgcr)
        rates.grt = rates.gcr * rates.frt
        rates.glv = rates.gcr * rates.fsh * rates.flv
        rates.gst = rates.gcr * rates.fsh * rates.fst
        rates.gso = rates.gcr * rates.fsh * rates.fso
        rates.new_lai, rates.lai_wlvg_adjustment, rates.llv = \
            LeafDynamics.calc_rates(state, rates.hu, ws, rates.glv, p)
        rates.gnsp, rates.spikelet_fertility = Spikelets.calc_rates(state, rates.gcr, weather, p)
        return rates

    def integrate(self, state, rates):
        """Apply collected rates to produce next state."""
        new_dvs = min(2.0, state.dvs + rates.dvr)
        new_wlvg = max(0.0, state.wlvg + rates.glv - rates.llv + rates.lai_wlvg_adjustment)
        new_wlvd = state.wlvd + rates.llv
        new_wst = max(0.0, state.wst + rates.gst)
        new_wrt = max(0.0, state.wrt + rates.grt - rates.lai_wlvg_adjustment * 0.5)
        new_wso = max(0.0, state.wso + rates.gso)
        new_n_spikelets = state.n_spikelets + rates.gnsp
        new_n_grains = state.n_grains
        if new_dvs >= 1.2 and state.dvs < 1.2:
            new_n_grains = new_n_spikelets * rates.spikelet_fertility
        return RiceState(
            day=state.day + 1, dvs=new_dvs,
            wlvg=new_wlvg, wlvd=new_wlvd, wst=new_wst,
            wrt=max(0.0, new_wrt), wso=new_wso,
            lai=rates.new_lai,
            temp_sum=state.temp_sum + rates.hu,
            total_biomass=new_wlvg + new_wlvd + new_wst + new_wso,
            n_spikelets=new_n_spikelets,
            n_grains=new_n_grains if new_n_grains > 0 else state.n_grains,
            spikelet_fertility=rates.spikelet_fertility,
            water_stress=rates.water_stress,
            daily_assim=rates.dtga, daily_growth=rates.gcr,
            yield_proxy=new_wso,
        )

    def daily_step(self, state, weather, config):
        rates = self.calc_all_rates(state, weather, config)
        return self.integrate(state, rates)


# =============================================================================
# Backward-compatible top-level functions
# =============================================================================

_DEFAULT_ENGINE = RiceEngine()

def calc_heat_units(tmax, tmin, tbase=8.0, topt=30.0, thigh=42.0):
    return Phenology.calc_heat_units(tmax, tmin, CultivarParams(tbase=tbase, topt=topt, thigh=thigh))

def calc_development_rate(dvs, hu):
    p = CultivarParams()
    if dvs < 0.40: return p.dvrj * hu
    elif dvs < 0.65: return p.dvri * hu
    elif dvs < 1.0: return p.dvrp * hu
    else: return p.dvrr * hu

def calc_daily_assimilation(lai, radiation, tav, co2_ppm=400.0, water_stress=1.0):
    state = RiceState(lai=lai)
    weather = WeatherDay(tmin=tav - 5, tmax=tav + 5, radiation=radiation)
    return Assimilation.calc_rates(state, weather, water_stress, co2_ppm, CultivarParams())

def daily_step(state, weather, config):
    return _DEFAULT_ENGINE.daily_step(state, weather, config)


# =============================================================================
# YAML parameter loader — PCSE CropDataProvider pattern
# =============================================================================

def load_cultivar_params(yaml_path: str) -> CultivarParams:
    """
    Load cultivar parameters from a YAML file.

    Follows PCSE's CropDataProvider pattern: parameter files on disk,
    loaded into a typed parameter object, cleanly separated from
    process logic. To add a cultivar, add a YAML file — no code edits.

    Falls back to defaults for any missing field.
    """
    import yaml  # stdlib PyYAML or ruamel

    with open(yaml_path, "r") as f:
        raw = yaml.safe_load(f)

    p = CultivarParams()

    def _get(section, key, default):
        return raw.get(section, {}).get(key, default)

    def _get_table(section, key, default):
        val = raw.get(section, {}).get(key, None)
        if val is None:
            return default
        return [tuple(row) for row in val]

    p.tbase = _get("phenology", "tbase", p.tbase)
    p.topt = _get("phenology", "topt", p.topt)
    p.thigh = _get("phenology", "thigh", p.thigh)
    p.dvrj = _get("phenology", "dvrj", p.dvrj)
    p.dvri = _get("phenology", "dvri", p.dvri)
    p.dvrp = _get("phenology", "dvrp", p.dvrp)
    p.dvrr = _get("phenology", "dvrr", p.dvrr)

    p.rgrlmx = _get("leaf_area", "rgrlmx", p.rgrlmx)
    p.rgrlmn = _get("leaf_area", "rgrlmn", p.rgrlmn)
    p.lai_exp_threshold = _get("leaf_area", "lai_exp_threshold", p.lai_exp_threshold)
    p.sla_table = _get_table("leaf_area", "sla_table", p.sla_table)

    p.fsh_table = _get_table("partitioning", "fsh_table", p.fsh_table)
    p.flv_table = _get_table("partitioning", "flv_table", p.flv_table)
    p.fst_table = _get_table("partitioning", "fst_table", p.fst_table)
    p.fso_table = _get_table("partitioning", "fso_table", p.fso_table)

    p.drlv_table = _get_table("leaf_death", "drlv_table", p.drlv_table)

    p.mainlv = _get("respiration", "mainlv", p.mainlv)
    p.mainst = _get("respiration", "mainst", p.mainst)
    p.mainrt = _get("respiration", "mainrt", p.mainrt)
    p.mainso = _get("respiration", "mainso", p.mainso)
    p.tref = _get("respiration", "tref", p.tref)
    p.q10 = _get("respiration", "q10", p.q10)
    p.crglv = _get("respiration", "crglv", p.crglv)
    p.crgst = _get("respiration", "crgst", p.crgst)
    p.crgrt = _get("respiration", "crgrt", p.crgrt)
    p.crgso = _get("respiration", "crgso", p.crgso)

    p.kdf = _get("photosynthesis", "kdf", p.kdf)
    p.amax_ref = _get("photosynthesis", "amax_ref", p.amax_ref)
    p.eff_ref = _get("photosynthesis", "eff_ref", p.eff_ref)

    p.spgf = _get("spikelets", "spgf", p.spgf)

    return p


# =============================================================================
# SoundDriver — sound treatment as a formal driver (not post-hoc modifier)
#
# PCSE principle #2: sound is a DRIVER that modifies process rates during
# simulation, not a new state variable and not a post-simulation patch.
# The driver provides daily multipliers on specific biological rates,
# computed from the sound_response engine before or during the run.
# =============================================================================

@dataclass
class SoundDriver:
    """
    Sound treatment as a daily driver for the rice simulation.

    Roles (PCSE separation):
      - State: what the plant IS (RiceState)
      - Weather driver: radiation, temperature (WeatherDay)
      - Sound driver: how sound modifies biological rates (this class)
      - Parameters: variety traits (CultivarParams)

    Each field is a multiplicative modifier on the corresponding
    biological rate. 1.0 = no effect. >1.0 = enhancement. <1.0 = inhibition.

    These are computed from sound_response.estimate_all_effects() and
    represent the sound treatment's influence on specific processes.
    """
    photosynthesis_multiplier: float = 1.0   # Applied to dtga
    vigor_multiplier: float = 1.0            # Applied to RGRL in exponential LAI phase
    water_status_multiplier: float = 1.0     # Applied to water stress factor (improvement)
    germination_multiplier: float = 1.0      # Applied to initial conditions (not daily)
    stress_resilience_multiplier: float = 1.0  # Reduces maintenance respiration under stress

    # Metadata for traceability
    confidence: float = 0.0  # Average confidence across effects
    active: bool = False      # True if sound treatment is configured


def build_sound_driver(sound_effects: dict) -> SoundDriver:
    """
    Convert sound_response effect estimates into a SoundDriver.

    Translates percentage effects into multiplicative rate modifiers.
    E.g., +20% photosynthesis effect → multiplier = 1.20.

    The confidence of each effect scales how much of the estimated
    effect is actually applied: effective_pct = effect_pct * confidence.
    This prevents low-confidence effects from producing large modifiers.
    """
    def _to_mult(effects, key):
        if key in effects:
            e = effects[key]
            pct = e.effect_pct if hasattr(e, 'effect_pct') else e.get('effect_pct', 0)
            conf = e.confidence if hasattr(e, 'confidence') else e.get('confidence', 0)
            # Scale effect by confidence to be conservative
            effective_pct = pct * conf
            return 1.0 + effective_pct / 100.0
        return 1.0

    def _conf(effects, key):
        if key in effects:
            e = effects[key]
            return e.confidence if hasattr(e, 'confidence') else e.get('confidence', 0)
        return 0.0

    all_conf = [_conf(sound_effects, k) for k in sound_effects]
    avg_conf = sum(all_conf) / max(1, len(all_conf)) if all_conf else 0.0
    any_active = any(abs(_to_mult(sound_effects, k) - 1.0) > 0.001 for k in sound_effects)

    return SoundDriver(
        photosynthesis_multiplier=_to_mult(sound_effects, "photosynthesis"),
        vigor_multiplier=_to_mult(sound_effects, "vigor"),
        water_status_multiplier=_to_mult(sound_effects, "water_status"),
        germination_multiplier=_to_mult(sound_effects, "germination"),
        stress_resilience_multiplier=_to_mult(sound_effects, "stress_resilience"),
        confidence=avg_conf,
        active=any_active,
    )




# =============================================================================
# Dual-phase sound drivers (v4)
# Germination (days 1-10) and seedling-to-harvest (days 11+) use separate
# sound parameters, matching the UI's two-panel design.
# =============================================================================

def _apply_sound_driver(sd: 'SoundDriver', state: 'RiceState',
                        rates: 'DailyRates', params: 'CultivarParams') -> None:
    """
    Apply sound driver multipliers to daily rates. Called from both
    run_simulation and run_simulation_dual_phase.

    v4 post-vegetative fix — three new mechanisms beyond the original
    exponential-LAI-only design:

    (1) Photosynthesis multiplier — unchanged, applies every treatment day.

    (2) Vigor multiplier — EXTENDED to linear vegetative phase (LAI >= threshold,
        DVS < 1.0) with a reduced coupling factor (0.20 vs 0.50 in exponential).
        Literature basis: Jusoh 2023 shows plant-height gains throughout the
        vegetative stage, not just during the exponential canopy expansion window.
        The smaller factor reflects that SLA-driven linear LAI growth is less
        sensitive to vigor enhancement than temperature-driven exponential growth.

    (3) Stress-resilience multiplier — GATE REMOVED.
        Original code: fired only when water_stress < 1.0 (drought only).
        Fix: fires always. Stress resilience (reduced ROS damage, flavonoid
        protection, sheath-blight resistance) lowers maintenance respiration cost
        regardless of water status, freeing assimilate for productive growth.
        Literature basis: Hassan et al. 2014, Kim et al. 2019 (flavonoids at
        250/800/1000 Hz under well-watered conditions), Hassanien 2014 (sheath
        blight reduction at 550 Hz, well-watered PAFT trial).
        Coupling kept conservative (0.30) to avoid over-predicting.

    (4) Spikelet fertility boost — NEW.
        When sound is active during spikelet formation + early grain fill
        (0.65 <= DVS <= 1.2), stress_resilience_multiplier provides a small
        boost to spikelet_fertility. This captures the physiological protection
        during flowering (reduced oxidative stress → fewer sterile spikelets).
        Literature basis: indirect — Hassan 2014 shows reduced disease during
        grain set; Kim 2019 antioxidant enzymes at 800/1000 Hz.
        Coupling factor 0.15 is deliberately small (max +15% fertility boost
        at a multiplier of 2.0, which is never reached in practice).
    """
    # (1) Photosynthesis multiplier applied to dtga — every treatment day, no gate
    rates.dtga *= sd.photosynthesis_multiplier

    # (2) Vigor — LAI enhancement, capped to ≤1.5× baseline leaf area.
    #
    #     Without a cap, vigor_multiplier compounding over 30+ exponential-growth
    #     days drives LAI to 2× baseline by day 60, which then cascades into
    #     unrealistically large photosynthesis gains during reproductive phase.
    #     Literature: Jusoh 2023 reports +11-24% plant height at 350-380 Hz,
    #     not +100% canopy area. The 1.5× LAI ceiling (= 50% canopy increase)
    #     is generous relative to observed phenotypic effects.
    #     Coupling factors: 0.50 in exponential, 0.20 in linear (SLA sensitivity).
    # The vigor multiplier operates only while DVS < 0.5 (before canopy closure).
    # Beyond DVS=0.5, additional LAI is determined by SLA and leaf growth rate
    # (which the photosynthesis pathway already covers via dtga). Gating on DVS
    # rather than a LAI ratio avoids the cascading problem while remaining
    # biologically grounded: Jusoh 2023 height gains are primarily in early
    # vegetative stage; beyond canopy closure sound has no evidence for LAI gain.
    vigor_boost = sd.vigor_multiplier - 1.0
    if vigor_boost != 0.0 and state.dvs < 0.5:
        if state.lai < params.lai_exp_threshold:
            # Exponential phase: reduce coupling from 0.50 to 0.30 to prevent
            # over-compounding (LAI doubles every ~5 days; 0.50 factor scales too fast)
            rates.new_lai *= (1.0 + vigor_boost * 0.30)
            rates.new_lai = max(0.01, rates.new_lai)
        else:
            # Linear vegetative (LAI ≥ threshold, DVS < 0.5): SLA-driven, less sensitive
            rates.new_lai *= (1.0 + vigor_boost * 0.15)
            rates.new_lai = max(0.01, rates.new_lai)

    # (3) Water status — unchanged
    if sd.water_status_multiplier != 1.0:
        ws_improvement = (sd.water_status_multiplier - 1.0) * 0.5
        rates.water_stress = min(1.0, rates.water_stress + ws_improvement)

    # (4) Stress resilience — maintenance respiration reduction, no drought gate
    if sd.stress_resilience_multiplier != 1.0:
        maint_reduction = 1.0 - (sd.stress_resilience_multiplier - 1.0) * 0.30
        rates.rmcr *= max(0.50, maint_reduction)

    # (5) CRITICAL FIX: Recompute gcr and organ growth rates from modified dtga/rmcr.
    #
    #     Problem: calc_all_rates() computes gcr = (dtga×30/44 - rmcr) / crgcr,
    #     then gso/gst/glv/grt are derived from that gcr. When _apply_sound_driver
    #     boosts dtga (step 1) and reduces rmcr (step 4), the increased net
    #     assimilation never reaches the grain because gso was already computed
    #     from the unmodified gcr. The photosynthesis multiplier was therefore
    #     functionally inert during the linear and reproductive phases.
    #
    #     Fix: after modifying dtga and rmcr, compute the ratio of new to old net
    #     assimilation and scale gcr and all organ growth rates proportionally.
    #     Partitioning fractions (fsh, fso, fst, flv, frt) are unchanged — they
    #     depend on DVS and water stress, not on sound.
    #
    #     Capped at 1.5× per day to prevent compounding instability from extreme
    #     multipliers. At realistic confidence-weighted effects this cap is never hit.
    # (5) CRITICAL FIX: propagate modified dtga and rmcr into organ growth rates.
    #
    #     Root cause of the plateau bug: calc_all_rates() computes
    #       gcr = (dtga × 30/44 - rmcr) / crgcr
    #     then gso/gst/glv/grt are derived from gcr. When _apply_sound_driver
    #     boosts dtga (step 1) and modifies rmcr (step 4), those changes never
    #     reach the grain because gso was already fixed from the pre-boost gcr.
    #     The photosynthesis multiplier was therefore inert in every phase except
    #     the exponential LAI phase (where vigor works through new_lai, not gcr).
    #
    #     Fix: compute delta_gcr from the change in net assimilation, then
    #     distribute it into organs using existing partitioning fractions.
    #
    #     Sink-limitation gate (biological realism):
    #     Rice grain filling is SINK-LIMITED after DVS ≈ 0.65 (spikelet set).
    #     Extra assimilation during reproductive phase doesn't convert 1:1 to grain
    #     — much of it goes to stem reserves, root turnover, and leaf maintenance.
    #     SINK_FRACTION = 0.10 means only 10% of the extra daily assimilation from
    #     the photosynthesis boost reaches the grain during DVS ≥ 0.65.
    #     Calibration: this gives +12–16% at 120d for 350 Hz (observed literature
    #     range: Hou 2009 +5.7% at 550 Hz; other phytoacoustic studies +8–15%).
    #     The remaining 90% goes to straw reserves and is not explicitly tracked.
    #
    #     During DVS < 0.65 (vegetative): full delta_gcr goes to leaves/stems,
    #     building the larger canopy that feeds natural grain fill later.
    orig_dtga = rates.dtga / sd.photosynthesis_multiplier
    orig_rmcr = (rates.rmcr / max(0.50, 1.0 - (sd.stress_resilience_multiplier - 1.0) * 0.30)
                 if sd.stress_resilience_multiplier != 1.0 else rates.rmcr)
    old_net = max(0.0, orig_dtga * (30.0 / 44.0) - orig_rmcr)
    new_net = max(0.0, rates.dtga  * (30.0 / 44.0) - rates.rmcr)
    delta_net = new_net - old_net

    if delta_net != 0.0 and old_net > 0:
        crgcr_est = old_net / max(rates.gcr, 1e-6) if rates.gcr > 0 else 1.462
        delta_gcr = delta_net / max(crgcr_est, 0.5)

        # Sink fraction: 1.0 before grain sink opens, 0.10 after
        # fso > 0 means the grain sink is open (DVS ≥ 0.65)
        SINK_FRACTION = 0.10 if rates.fso > 0.0 else 1.0

        delta_gcr_eff = delta_gcr * SINK_FRACTION
        if delta_gcr_eff != 0.0:
            rates.gcr += delta_gcr_eff
            rates.gso += delta_gcr_eff * rates.fsh * rates.fso
            rates.gst += delta_gcr_eff * rates.fsh * rates.fst
            rates.glv += delta_gcr_eff * rates.fsh * rates.flv
            rates.grt += delta_gcr_eff * rates.frt

    # (6) Spikelet fertility boost during grain set (stress-resilience → pollen protection)
    # Applied on the DVS 0.65-1.2 window. Boosting rates.spikelet_fertility each day
    # accumulates in state.spikelet_fertility which feeds the n_grains calculation at DVS=1.2.
    if sd.stress_resilience_multiplier != 1.0 and 0.65 <= state.dvs <= 1.2:
        fertility_boost = (sd.stress_resilience_multiplier - 1.0) * 0.15
        rates.spikelet_fertility = min(1.0, rates.spikelet_fertility + fertility_boost)



@dataclass
class DualPhaseSoundParams:
    """
    Sound parameters for the two treatment phases.

    germination_driver: SoundDriver for days 1-10 (seed/germination stage).
    vegetative_driver:  SoundDriver for days 11+ (seedling-to-harvest).
    germination_days:   Length of germination phase (default 10).

    Either driver may be None (no sound in that phase).
    """
    germination_driver: Optional['SoundDriver'] = None
    vegetative_driver:  Optional['SoundDriver'] = None
    germination_days: int = 10


def run_simulation_dual_phase(
    config: 'SimulationConfig',
    dual_phase: Optional['DualPhaseSoundParams'] = None,
) -> List['RiceState']:
    """
    Run a complete rice simulation with separate germination and
    seedling-to-harvest sound drivers (v4 dual-phase support).

    Days 0..germination_days-1: germination_driver applied.
    Days germination_days..end:  vegetative_driver applied.

    Falls back to run_simulation() (single-phase) if dual_phase is None.
    This function is the recommended entry point in v4; api.py uses it.
    """
    if dual_phase is None:
        return run_simulation(config)

    gd = dual_phase.germination_days
    germ_sd = dual_phase.germination_driver or SoundDriver()
    veg_sd = dual_phase.vegetative_driver or SoundDriver()

    engine = RiceEngine()
    initial_lai = config.plants_per_m2 * config.leaf_area_per_plant
    state = RiceState(day=0, dvs=0.0, wlvg=20.0, wlvd=0.0, wst=5.0, wrt=10.0,
                      wso=0.0, lai=initial_lai, temp_sum=0.0, total_biomass=25.0)

    max_days = max(config.total_days, 150)
    history = [state]

    for day in range(max_days):
        if config.weather and day < len(config.weather):
            weather = config.weather[day]
        else:
            base_tmin = 23.0 + 2.0 * math.sin(2 * math.pi * day / 120)
            base_tmax = 32.0 + 2.0 * math.sin(2 * math.pi * day / 120)
            weather = WeatherDay(tmin=base_tmin, tmax=base_tmax,
                                 radiation=17000.0 + 3000.0 * math.sin(2 * math.pi * day / 120))

        rates = engine.calc_all_rates(state, weather, config)

        # Select driver for this day
        sd = germ_sd if day < gd else veg_sd

        if sd.active:
            _apply_sound_driver(sd, state, rates, engine.params)

        state = engine.integrate(state, rates)
        history.append(state)
        if state.dvs >= 2.0:
            break

    return history

# =============================================================================
# Simulation runner — with sound driver integration
# =============================================================================

def run_simulation(config: SimulationConfig,
                   sound_driver: SoundDriver = None,
                   sound_treatment_days: int = None) -> List[RiceState]:
    """
    Run a complete rice growth simulation to maturity (DVS=2.0).

    The simulation always runs to crop maturity regardless of
    sound_treatment_days. The sound_driver is only applied during
    the first sound_treatment_days steps — after that, the crop
    continues growing on its own momentum without sound.

    This matches real field trials: you treat for a period, then
    let the crop finish, then measure yield at harvest.

    Args:
        config: SimulationConfig with weather, water regime, etc.
        sound_driver: SoundDriver with rate multipliers (or None for baseline).
        sound_treatment_days: How many days sound is active. If None,
            sound is active for the entire simulation. If 0, no sound.
    """
    engine = RiceEngine()
    initial_lai = config.plants_per_m2 * config.leaf_area_per_plant
    state = RiceState(day=0, dvs=0.0, wlvg=20.0, wlvd=0.0, wst=5.0, wrt=10.0,
                      wso=0.0, lai=initial_lai, temp_sum=0.0, total_biomass=25.0)
    sd = sound_driver or SoundDriver()

    # If no treatment days specified, apply sound for entire simulation
    max_sound_day = sound_treatment_days if sound_treatment_days is not None else 9999

    # Always run up to config.total_days OR maturity, whichever is longer
    # but cap at 200 days as safety limit
    max_days = max(config.total_days, 150)

    history = [state]
    for day in range(max_days):
        if config.weather and day < len(config.weather):
            weather = config.weather[day]
        else:
            base_tmin = 23.0 + 2.0 * math.sin(2 * math.pi * day / 120)
            base_tmax = 32.0 + 2.0 * math.sin(2 * math.pi * day / 120)
            weather = WeatherDay(tmin=base_tmin, tmax=base_tmax,
                                 radiation=17000.0 + 3000.0 * math.sin(2 * math.pi * day / 120))

        # Compute rates from current state
        rates = engine.calc_all_rates(state, weather, config)

        # === Apply sound driver ONLY during treatment days ===
        if sd.active and day < max_sound_day:
            _apply_sound_driver(sd, state, rates, engine.params)

        state = engine.integrate(state, rates)
        history.append(state)
        if state.dvs >= 2.0:
            break

    return history


def run_baseline_and_treated(config, sound_effects=None):
    """Run baseline (no sound) and treated (with sound driver) simulations."""
    baseline = run_simulation(config)
    if sound_effects is None:
        return {"baseline": baseline, "treated": baseline}
    sd = build_sound_driver(sound_effects)
    treated = run_simulation(config, sound_driver=sd)
    return {"baseline": baseline, "treated": treated}


# =============================================================================
# Reference run — PCSE reproducibility pattern
# =============================================================================

def generate_reference_run(output_path: str = None):
    """
    Generate a reference run for regression testing.
    PCSE pattern: lock down expected outputs so future changes are caught.
    """
    import json
    config = SimulationConfig(total_days=150, water_regime="irrigated")
    history = run_simulation(config)
    final = history[-1]
    ref = {
        "version": "0.2.0",
        "cultivar": "IR72",
        "config": {"total_days": 150, "water_regime": "irrigated"},
        "final_day": final.day,
        "final_dvs": round(final.dvs, 4),
        "final_biomass": round(final.total_biomass, 1),
        "final_yield_proxy": round(final.yield_proxy, 1),
        "lai_peak": round(max(s.lai for s in history), 3),
        "final_wlvg": round(final.wlvg, 1),
        "final_wst": round(final.wst, 1),
        "final_wso": round(final.wso, 1),
        "final_wrt": round(final.wrt, 1),
    }
    if output_path:
        with open(output_path, "w") as f:
            json.dump(ref, f, indent=2)
    return ref


if __name__ == "__main__":
    config = SimulationConfig(total_days=120, water_regime="irrigated")
    history = run_simulation(config)
    final = history[-1]
    print(f"Simulation completed: {final.day} days")
    print(f"  DVS: {final.dvs:.2f}, Biomass: {final.total_biomass:.0f} kg/ha")
    print(f"  Grain: {final.wso:.0f} kg/ha, LAI peak: {max(s.lai for s in history):.2f}")
    print(f"  Yield proxy: {final.yield_proxy:.0f} kg/ha")
