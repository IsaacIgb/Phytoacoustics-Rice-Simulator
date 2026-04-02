"""
kpi_engine.py — KPI computation layer.

Every KPI is based on explicit variables and formulas.
Each KPI returns value, confidence, and explanation.

KPIs:
  1. Yield Index — relative yield change from sound treatment
  2. Water Index — relative improvement in water-related metrics
  3. Stress Resilience Index — improvement in stress physiology
  4. Coverage Quality Index — field acoustic exposure quality
"""

from dataclasses import dataclass
from typing import List, Optional
import math


@dataclass
class KPIResult:
    """A single KPI result with value, confidence, and explanation."""
    name: str
    value: float           # Normalized 0-100 or relative % change
    confidence: float      # 0-1
    explanation: str
    components: dict        # Underlying variable values


def compute_yield_index(
    baseline_yield: float,
    treated_yield: float,
    sound_yield_effect_pct: float,
    sound_confidence: float,
) -> KPIResult:
    """
    Yield Index: relative yield improvement from sound treatment.
    
    Formula:
      yield_index = ((treated_yield - baseline_yield) / baseline_yield) * 100
      
    The treated_yield is baseline_yield * (1 + sound_yield_effect_pct/100).
    Confidence comes from the sound evidence and simulation quality.
    """
    if baseline_yield <= 0:
        return KPIResult(
            name="Yield Index",
            value=0.0,
            confidence=0.0,
            explanation="Baseline yield is zero; cannot compute index.",
            components={"baseline": baseline_yield, "treated": treated_yield},
        )
    
    # Apply sound effect to baseline
    effective_treated = baseline_yield * (1.0 + sound_yield_effect_pct / 100.0)
    
    pct_change = ((effective_treated - baseline_yield) / baseline_yield) * 100.0
    
    # Confidence is the sound evidence confidence, reduced if effect is very large
    conf = sound_confidence
    if abs(pct_change) > 30:
        conf *= 0.5  # Large effects get confidence penalty
    
    # Normalize to 0-100 scale: 50 = no change, >50 = improvement
    normalized = 50.0 + pct_change
    normalized = max(0.0, min(100.0, normalized))
    
    return KPIResult(
        name="Yield Index",
        value=round(normalized, 1),
        confidence=round(conf, 3),
        explanation=f"Sound treatment estimated to change yield by {pct_change:+.1f}%. "
                    f"Based on {sound_confidence:.0%} confidence evidence.",
        components={
            "baseline_yield_kg_ha": round(baseline_yield, 0),
            "estimated_treated_yield_kg_ha": round(effective_treated, 0),
            "pct_change": round(pct_change, 2),
        },
    )


def compute_water_index(
    water_stress_baseline: float,
    sound_water_effect_pct: float,
    sound_confidence: float,
    regime: str = "irrigated",
) -> KPIResult:
    """
    Water Index: improvement in water-related metrics.
    
    Formula:
      water_improvement = sound_water_effect_pct (from evidence)
      base_water_score = water_stress_baseline * 100
      water_index = base_water_score + water_improvement
      
    For irrigated: baseline stress is ~1.0 (no stress)
    For drought: baseline stress < 1.0, sound may improve RWC
    """
    base_score = water_stress_baseline * 100.0
    improvement = sound_water_effect_pct * (sound_confidence)
    
    water_index = base_score + improvement
    water_index = max(0.0, min(100.0, water_index))
    
    # Confidence higher under drought (more evidence for drought context)
    conf = sound_confidence
    if regime == "irrigated" and abs(sound_water_effect_pct) < 2:
        conf *= 0.5  # Less meaningful under no-stress conditions
    
    return KPIResult(
        name="Water Index",
        value=round(water_index, 1),
        confidence=round(conf, 3),
        explanation=f"Water status score: {base_score:.0f}/100 baseline. "
                    f"Sound estimated to improve by {improvement:+.1f} points.",
        components={
            "baseline_water_stress": round(water_stress_baseline, 3),
            "sound_water_effect_pct": round(sound_water_effect_pct, 2),
            "regime": regime,
        },
    )


def compute_stress_resilience_index(
    sound_photosynthesis_effect: float,
    sound_stress_effect: float,
    sound_vigor_effect: float,
    photo_confidence: float,
    stress_confidence: float,
    vigor_confidence: float,
    water_stress: float = 1.0,
) -> KPIResult:
    """
    Stress Resilience Index: composite of physiological improvements.
    
    Formula:
      resilience = weighted_mean(
        photosynthesis_effect * 0.4,
        stress_resilience_effect * 0.3,
        vigor_effect * 0.3
      )
      
    Normalized to 0-100 scale (50 = no change).
    """
    # Weight by both magnitude and confidence
    w_photo = 0.4
    w_stress = 0.3
    w_vigor = 0.3
    
    weighted_effect = (
        sound_photosynthesis_effect * w_photo * photo_confidence +
        sound_stress_effect * w_stress * stress_confidence +
        sound_vigor_effect * w_vigor * vigor_confidence
    )
    
    total_conf_weight = (
        w_photo * photo_confidence +
        w_stress * stress_confidence +
        w_vigor * vigor_confidence
    )
    
    if total_conf_weight > 0:
        avg_effect = weighted_effect / total_conf_weight
    else:
        avg_effect = 0.0
    
    avg_confidence = total_conf_weight / (w_photo + w_stress + w_vigor)
    
    # Under stress, resilience matters more
    if water_stress < 0.8:
        stress_relevance = 1.5  # Amplify under drought
    else:
        stress_relevance = 1.0
    
    normalized = 50.0 + avg_effect * stress_relevance
    normalized = max(0.0, min(100.0, normalized))
    
    return KPIResult(
        name="Stress Resilience Index",
        value=round(normalized, 1),
        confidence=round(avg_confidence, 3),
        explanation=f"Composite stress resilience: photosynthesis ({sound_photosynthesis_effect:+.1f}%), "
                    f"oxidative defense ({sound_stress_effect:+.1f}%), "
                    f"vigor ({sound_vigor_effect:+.1f}%).",
        components={
            "photosynthesis_effect_pct": round(sound_photosynthesis_effect, 2),
            "stress_resilience_effect_pct": round(sound_stress_effect, 2),
            "vigor_effect_pct": round(sound_vigor_effect, 2),
            "water_stress": round(water_stress, 3),
        },
    )


def compute_coverage_quality_index(
    mean_spl: float,
    coverage_above_60db: float,
    coverage_above_70db: float,
    uniformity: float,
    hours_per_day: float,
) -> KPIResult:
    """
    Coverage Quality Index: how well the field is acoustically exposed.
    
    Formula:
      coverage_score = coverage_above_60db * 40 + coverage_above_70db * 20
      uniformity_score = uniformity * 20
      dosage_score = min(hours_per_day / 3.0, 1.0) * 20
      
      quality_index = coverage_score + uniformity_score + dosage_score
      
    Scale: 0-100. Higher = better exposure.
    """
    # Coverage component (0-60 points)
    cov_score = coverage_above_60db * 40.0 + coverage_above_70db * 20.0
    
    # Uniformity component (0-20 points)
    uni_score = uniformity * 20.0
    
    # Dosage component (0-20 points)
    # Evidence suggests 2-4 h/day is typical. Scale linearly to 3h.
    dose_ratio = min(hours_per_day / 3.0, 1.0) if hours_per_day > 0 else 0.0
    dose_score = dose_ratio * 20.0
    
    total = cov_score + uni_score + dose_score
    total = max(0.0, min(100.0, total))
    
    # Confidence is high for coverage (physics-based) but varies with uniformity
    conf = 0.8 if uniformity > 0.5 else 0.6
    
    # Warn if mean SPL is very high (potential harm)
    notes = ""
    if mean_spl > 100:
        notes = " Warning: mean SPL >100 dB may cause tissue damage."
        conf *= 0.5
    
    return KPIResult(
        name="Coverage Quality Index",
        value=round(total, 1),
        confidence=round(conf, 3),
        explanation=f"Coverage >60dB: {coverage_above_60db*100:.0f}%, "
                    f">70dB: {coverage_above_70db*100:.0f}%, "
                    f"uniformity: {uniformity:.2f}, "
                    f"dosage: {hours_per_day}h/day.{notes}",
        components={
            "coverage_60db_pct": round(coverage_above_60db * 100, 1),
            "coverage_70db_pct": round(coverage_above_70db * 100, 1),
            "uniformity": round(uniformity, 3),
            "mean_spl_db": round(mean_spl, 1),
            "hours_per_day": hours_per_day,
        },
    )


def compute_all_kpis(
    baseline_yield: float,
    sound_effects: dict,
    water_stress: float,
    water_regime: str,
    acoustics_result: dict,
    hours_per_day: float,
) -> dict:
    """
    Compute all four KPIs from simulation results and sound effects.
    
    sound_effects: dict of outcome -> EffectEstimate (or dict with effect_pct, confidence)
    acoustics_result: dict with mean_spl, coverage_above_60db, etc.
    """
    def _get(effects, key, attr="effect_pct"):
        if key in effects:
            e = effects[key]
            if hasattr(e, attr):
                return getattr(e, attr)
            elif isinstance(e, dict):
                return e.get(attr, 0.0)
        return 0.0
    
    def _conf(effects, key):
        return _get(effects, key, "confidence")
    
    yield_kpi = compute_yield_index(
        baseline_yield=baseline_yield,
        treated_yield=baseline_yield,
        sound_yield_effect_pct=_get(sound_effects, "yield"),
        sound_confidence=_conf(sound_effects, "yield"),
    )
    
    water_kpi = compute_water_index(
        water_stress_baseline=water_stress,
        sound_water_effect_pct=_get(sound_effects, "water_status"),
        sound_confidence=_conf(sound_effects, "water_status"),
        regime=water_regime,
    )
    
    stress_kpi = compute_stress_resilience_index(
        sound_photosynthesis_effect=_get(sound_effects, "photosynthesis"),
        sound_stress_effect=_get(sound_effects, "stress_resilience"),
        sound_vigor_effect=_get(sound_effects, "vigor"),
        photo_confidence=_conf(sound_effects, "photosynthesis"),
        stress_confidence=_conf(sound_effects, "stress_resilience"),
        vigor_confidence=_conf(sound_effects, "vigor"),
        water_stress=water_stress,
    )
    
    coverage_kpi = compute_coverage_quality_index(
        mean_spl=acoustics_result.get("mean_spl", 0),
        coverage_above_60db=acoustics_result.get("coverage_above_60db", 0),
        coverage_above_70db=acoustics_result.get("coverage_above_70db", 0),
        uniformity=acoustics_result.get("uniformity", 0),
        hours_per_day=hours_per_day,
    )
    
    return {
        "yield_index": yield_kpi,
        "water_index": water_kpi,
        "stress_resilience_index": stress_kpi,
        "coverage_quality_index": coverage_kpi,
    }
