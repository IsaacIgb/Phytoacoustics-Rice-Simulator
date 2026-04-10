# Mathematical Specification — Phytoacoustics Rice Simulator (v4)

This document specifies every equation in the simulator with rationale. All section references to ORYZA2000 refer to Bouman et al. (2001).

---

## 1. Phenology

**Eq P1** — Heat unit accumulation (bilinear temperature response):
```
HU = max(0, min(Tav, Topt) - Tbase) × correction_for_Thigh
Tav = (Tmax + Tmin) / 2
```
Clamped at Tbase=8°C (no growth below), Topt=30°C (full rate), Thigh=42°C (decline).

**Eq P2** — Development rate (DVR):
```
DVS < 1.0 (vegetative):  DVR = DVRJ × HU
DVS ≥ 1.0 (reproductive): DVR = DVRR × HU
DVS accumulates from 0 (emergence) to 2.0 (maturity).
```

---

## 2. Photosynthesis

**Eq A1** — Big-leaf canopy assimilation (analytical Monsi-Saeki integral):
```
DTGA = Amax_canopy × (1 - exp(-K × LAI)) / K × f(Rad) × f(Tav) × water_stress
```
Amax_canopy = Amax_leaf × temperature_scaling; K = 0.6 extinction coefficient.

**Eq A2** — CO2 elevation modifier:
```
CO2_factor = 1 + 0.0008 × (CO2_ppm - 350)   [approx +0.08% per ppm above 350]
```

---

## 3. Respiration

**Eq R1** — Maintenance respiration (Q10 temperature response):
```
RMCR = Q10^((Tav - Tref) / 10) × Σ_organs(maintenance_coeff × organ_weight)
```

**Eq R2** — Growth respiration coefficient (assimilate cost per unit DM):
```
CRG = Σ_organs(fraction × growth_coeff)
GCR = max(0, (CH2O - RMCR) / CRG)
```

---

## 4. Leaf area dynamics

**Eq L1** — Exponential phase (LAI < 1.0):
```
dLAI/dt = LAI × RGRL × f(Tav) × water_stress
RGRL = RGRLMX (nitrogen not modelled)
```

**Eq L2** — Linear phase (LAI ≥ 1.0):
```
dLAI/dt = GLV × SLA(DVS)
GLV = GCR × FSH × FLV
```

**Eq L3** — Leaf death:
```
LLV = DRLV(DVS) × WLVG
```

---

## 5. Sound-Response Engine (v4)

### Architecture

Four layers + two special pathways:

1. **Biophysical boundaries** — hard physics limits
2. **Mechanotransduction latent S_mech,air** — SPL × frequency × airborne coupling
3. **Pathway activations + mechanistic prior** — biological frequency selectivity
4. **Gaussian Process on residuals** — per-outcome data correction

Special: ultrasound cavitation (≥20 kHz); iWUE nearest-neighbour.

### Layer 1: Biophysical boundaries

**Eq SE1** — Activation floor:
```
if SPL > 0 and SPL < 60 dB → effect = 0   (no mechanosensory response)
```

**Eq SE2** — Damage ceiling:
```
if SPL > 110 dB → effect = -15%, conf = 0.45   (membrane rupture, hard override)
```

**Eq SE2a** — Broadband high-SPL inhibition (90–110 dB, no SPL-matched data):
```
frac = (SPL - 90) / (110 - 90)                   # 0 at 90 dB, 1 at 110 dB
modifier = -0.3 - 0.2 × frac                      # -0.30 at 90 → -0.50 at 110
effect_pct = modifier × 100                        # -30% at 90 → -50% at 110 dB
```
Applied when `spl_matched = False`: no CSV data within ±14 dB SPL **and** ±150 Hz of the query.
Justified by GrassGeneral_90dB (-40% vegetative growth). Excluded from GP (boundary-layer study).

**Eq SE2b** — S_SPL(SPL) — piecewise linear activation gate ∈ [0, 1]:
```
SPL < 60 dB:           S_SPL = 0
60 ≤ SPL < 85 dB:      S_SPL = (SPL - 60) / 25
85 ≤ SPL ≤ 100 dB:     S_SPL = 1.0
100 < SPL < 110 dB:    S_SPL = 1 - (SPL - 100) / 10
SPL ≥ 110 dB:          S_SPL = 0
```

**Eq SE2c** — S_air(f) — airborne coupling decay ∈ [0, 1]:
```
f ≤ 5000 Hz:            S_air = 1.0
5000 < f < 15000 Hz:    S_air = 1 - (f - 5000) / 10000
f ≥ 15000 Hz:           S_air = 0.0
f ≥ 20000 Hz:           ultrasound pathway (bypasses airborne)
```

**Eq SE2d** — S_f(f) — MechActivation GP frequency selectivity ∈ [0, 1]:
```
MechActivation GP (Matérn ν=5/2) trained on pooled stress/activation outcomes:
  Fv/Fm, total_flavonoid_content, sheath_blight_incidence*,
  stomatal_conductance_drought, relative_water_content,
  photosynthesis_rate, chlorophyll_content, SOD/CAT/POD enzymes
  (* sheath_blight_incidence is pre-normalised in CSV: +50 = beneficial)

Y_max = 95th percentile of GP mean on 50–5000 Hz grid.
S_f(f) = clip(y_mech(f) / Y_max, 0, 1)
```

**Eq SE2e** — S_mech,air(f, SPL) ∈ [0, 1]:
```
S_mech,air = S_f(f) × S_SPL(SPL) × S_air(f)
Returns 0 for f ≥ 20000 Hz (ultrasound pathway applies).
```

### Layer 3: Pathway activations

**Eq SE3a** — Log-Gaussian pathway activation p(pathway, f, SPL, stage, regime):
```
z = (ln(f) - mu_log_hz) / sigma_log_hz
log_gauss = exp(-0.5 × z²)

Stage gate:  0 if stage_bucket not in pathway.applicable_stages (mixed passes all)
Regime gate: 0 if stress_context not in pathway.applicable_regimes
SPL gate:    0 if spl_db < pathway.min_spl_db

For airborne:
  smech_scale = max(0.15, S_mech,air) if spl_db ≥ 60 dB else 0
  p = log_gauss × smech_scale

For ultrasound pathway:
  p = log_gauss   (S_mech,air not used)
```

Five pathway configurations:

| Pathway | mu_log_hz | sigma | min_spl | peak_pct |
|---------|-----------|-------|---------|----------|
| seed_vigor | ln(380) | 0.30 | 0 dB | +25% |
| leaf_gas_exchange | ln(357) | 0.12 | 0 dB | +40% |
| drought_resilience | ln(800) | 0.55 | 0 dB | +45% |
| membrane_damage | ln(4000) | 0.25 | 100 dB | -20% |
| ultrasound_cavitation | ln(38000) | 0.18 | 0 dB | +20% |

**Eq SE3b** — Mechanistic prior for outcome:
```
prior_mech(f) = Σ_pathway [ p(pathway, f, ...) × weight × peak_pct ] / Σ weight
```

**Eq SE3c** — Hybrid prior (window floor):
```
prior_window(f) = quality-weighted midpoint of matching empirical windows
If |prior_mech(f)| < 1% AND |prior_window(f)| ≥ 1%:
    prior(f) = prior_window(f)      # direct-evidence floor
Else:
    prior(f) = prior_mech(f)
```
Preserves direct literature evidence when pathway geometry is too narrow (e.g. Hou 2009 yield at 550 Hz, 3.6σ from leaf_gas_exchange centre).

OUTCOME_PATHWAY_WEIGHTS (initial values, comment-documented in code):

| Outcome | Pathways (weight) |
|---------|------------------|
| photosynthesis | leaf_gas_exchange(0.8), seed_vigor(0.2) |
| vigor | seed_vigor(0.6), leaf_gas_exchange(0.2), drought_resilience(0.2) |
| water_status_gsw_drought | drought_resilience(1.0) |
| water_status_gsw_wellwatered | leaf_gas_exchange(1.0) |
| water_status_rwc | drought_resilience(0.8), leaf_gas_exchange(0.2) |
| stress_resilience | drought_resilience(0.7), membrane_damage(0.3) |
| germination | seed_vigor(0.9), drought_resilience(0.1) |
| yield | leaf_gas_exchange(0.6), seed_vigor(0.4) |

Weight adjustment rules: change at most 0.1 per weight per edit; document reason in code comment.

### Layer 4: Gaussian Process

**Eq SE4** — Frequency aggregation (pre-GP, prevents pseudo-replication):
```
Group points within ±1 Hz. For each cluster:
  precision_i = 1 / α_i
  mean_freq   = Σ(f_i × precision_i) / Σ precision_i
  mean_effect = Σ(eff_i × precision_i) / Σ precision_i
  combined_α  = 1 / Σ precision_i

Exception: water_status_efficiency (iWUE) — aggregation skipped.
  The Jusoh 2023 iWUE data oscillates ±60% over 30 Hz (standing-wave artefact).
  This oscillation is real signal for nearest-neighbour lookup; merging
  357/359 Hz would collapse the ±62.7% at 359 Hz.
```

**Eq SE5** — Per-study noise variance (heteroscedastic):
```
conf_w = {high: 1.0, medium: 0.6, low: 0.3}
comp_w = {full: 1.0, partial: 0.5, insufficient: 0.2}
reliability = conf_w × comp_w
spl_w = max(0.2, 1 - |SPL_query - SPL_row| / 40)   (when both known)
effective_rel = max(0.01, reliability × spl_w)
species_mult = (1/0.6)² ≈ 2.78 if crop_type ≠ rice
α_i = min(500, BASE_NOISE × (1 / effective_rel) × species_mult)   BASE_NOISE=30
```

**Eq SE6** — GP kernel (Matérn ν=5/2):
```
k(f, f') = C × (1 + √5 r/l + 5r²/3l²) × exp(-√5 r/l),  r = |f - f'|
C initialised from mean(residuals²)   [amplitude scaled to actual residual variance]
l initial = 200 Hz, bounds = [100, 1200] Hz  (n ≥ 3: optimised; n < 3: fixed)
5-restart L-BFGS-B optimisation of log-marginal-likelihood (n ≥ 3)
```

**Eq SE7** — GP residuals and prediction:
```
residual_i = effect_i - prior(f_i)      [hybrid prior, same function as query]

GP.fit(freqs_aggregated, residuals, alpha=α_i)

GP.predict(f_query) → (μ_residual, σ_residual)
effect = prior(f_query) + μ_residual
std    = σ_residual
```
Consistent prior: the same `_hybrid_prior()` helper computes prior at both training points and query, preventing residual/prediction mismatch.

**Eq SE8** — GP confidence with local density:
```
effect_scale   = max(|effect|, 15)
proximity_conf = exp(-0.5 × (std / effect_scale)²)

l_scale  = fitted length scale (or 200 Hz if n < 3)
n_local  = count of aggregated training clusters within l_scale of f_query
n_factor = min(1.0, n_local / 2.0)

confidence = proximity_conf × (0.3 + 0.7 × n_factor)
```
Local density: a study at 800 Hz earns confidence only for queries near 800 Hz. This means a high-quality study at 800 Hz earns the same local credit as one at 350 Hz, eliminating the low-band cluster bias.

**Eq SE9** — Dosage scaling (airborne only):
```
hours < 1.0:         dose_factor = √hours
1.0 ≤ hours ≤ 4.0:  dose_factor = 1.0
hours > 4.0:         dose_factor = min(1.3, 1.0 + 0.1 × ln(hours/4))
Ultrasound:          dose_factor = 1.0 always
```

**Eq SE10** — Effect cap:
```
|effect| > 100% → clamp, confidence × 0.5
```

**Eq SE11** — High-frequency confidence decay (5–15 kHz):
```
f < 5 kHz:      modifier = 1.0
5 < f < 15 kHz: modifier = 1 - (f - 5000) / 10000
f ≥ 15 kHz:     modifier = 0
f ≥ 20 kHz:     ultrasound pathway
```

### Water status sub-KPIs (CWSI theory)

**Eq WS1** — Aggregate water_status:
```
Drought:      WS = 0.60 × gsw_drought + 0.25 × rwc + 0.15 × iWUE
Well-watered: WS = 0.40 × gsw_wellwatered + 0.40 × iWUE + 0.20 × rwc
```

### Sign normalisation

Outcomes where more-negative raw value = beneficial: germination_time, germination_time_h, speed_of_germination → sign flipped before GP.

`sheath_blight_incidence` is NOT in this list. The CSV stores +50 meaning "+50% disease reduction" (already beneficial). Flipping would yield -50 = harmful. Documented in `_FLIP_SIGN_OUTCOMES` comment in `sound_response.py`. Unit test `test_v4_architecture` asserts positive effect at 550 Hz.

### Ultrasound pathway (≥20 kHz)

**Eq US1** — Distance-weighted average:
```
For each matching row at f_row:
  prox = max(0.05, 1 - |f_query - f_row| / 30000)
  weighted sum of effect × conf × prox
effect = total_effect / total_weight
```

### iWUE nearest-neighbour

**Eq NN1** — iWUE lookup:
```
Select row with minimum |f_query - f_point|
effect = effect at that point
base_conf = min(1.0, BASE_NOISE / α_nearest)
confidence = base_conf × max(0.05, 1 - distance/50) × (0.3 + 0.7 × min(1, n/4))
```

---

## 6. Acoustics

**Eq AC1** — Inverse-square SPL at distance d:
```
SPL(d) = SPL_1m - 20 × log10(d)
```

**Eq AC2** — Coherent source combination:
```
SPL_combined = 10 × log10(Σ 10^(SPL_i / 10))
```

---

## 7. KPI computation

**Yield index:** Baseline 50; increases with positive yield effect × confidence.
**Water index:** Near 100 under irrigation; reflects sound-induced stomatal improvement.
**Stress resilience index:** Weighted combination of photosynthesis, stress, and vigor effects.
**Coverage quality:** SPL coverage above 60/70 dB × uniformity × dosage.

---

## 8. Sound-to-crop modifiers (daily loop)

```
Photosynthesis multiplier → rates.dtga  ×= multiplier
Vigor multiplier          → rates.new_lai ×= (1 + boost × 0.5)  during exponential phase
Water status multiplier   → rates.water_stress += (mult - 1) × 0.5
Stress resilience         → rates.rmcr ×= max(0.5, 1 - (mult-1) × 0.3)  under stress
```

Dual-phase: germination driver for days < germ_days; vegetative driver for days ≥ germ_days.

---

## Cross-reference: ORYZA2000 vs simulator

| Simulator section | ORYZA2000 | Notes |
|-------------------|-----------|-------|
| P1, P2 | Section 3.2.1 | Bilinear HU instead of full phenology model |
| A1, A2 | Section 3.2.2 | Big-leaf approximation |
| R1, R2 | Section 3.2.5 | Tropical rice coefficients |
| L1–L3 | Section 3.2.9 | Exponential phase is PCSE adaptation |
| SE1–SE11, WS1, US1, NN1 | — | Sound engine, no ORYZA equivalent |
| AC1, AC2 | — | Acoustics, no ORYZA equivalent |
