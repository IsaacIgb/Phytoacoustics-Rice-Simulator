# Mathematical Specification

Complete specification of every equation in the simulator, with rationale for each choice versus alternatives.

---

## 1. Phenology (crop development timing)

### Equations

**Eq P1** — Hourly temperature (sinusoidal approximation):
```
Td(h) = (Tmin + Tmax)/2 + (Tmax - Tmin)/2 × cos(0.2618 × (h - 14))
```
where h = hour (1–24), 0.2618 = π/12 (15° in radians).

**Eq P2** — Hourly heat unit (bilinear response):
```
Td ≤ Tbase or Td ≥ Thigh: HUh = 0
Tbase < Td ≤ Topt:        HUh = (Td - Tbase) / 24
Topt < Td < Thigh:        HUh = (Topt - (Td-Topt) × (Topt-Tbase)/(Thigh-Topt)) / 24
```
Parameters: Tbase=8°C, Topt=30°C, Thigh=42°C.

**Eq P3** — Daily heat units: `HU = Σ(h=1..24) HUh`

**Eq P4** — Development rate: `DVR = DVR_phase × HU`

**Eq P5** — Drought delay: `if DVS < 1.0 and ws < 1.0: DVR *= ws + DVS×(1-ws)`

**Eq P6** — State update: `DVS(t+1) = min(2.0, DVS(t) + DVR)`

### Rationale
The bilinear hourly approach was chosen over simpler daily-mean methods because rice development rate is non-linear near cardinal temperatures. A sinusoidal daily temperature course (P1) combined with hourly heat unit accumulation (P2-P3) correctly handles the curvilinear response that occurs when daily temperatures span above and below the optimum. This is the exact approach used in ORYZA2000 (Bouman et al. 2001, Section 3.2.1), following Matthews and Hunt (1994).

**Why not simpler?** A daily-mean approach (`HU = max(0, Tav - Tbase)`) would overestimate heat units on hot days and underestimate them on days that span the optimum. The 24-point hourly integration adds negligible compute cost and significantly improves accuracy near cardinal temperatures.

**Why not more complex?** ORYZA2000's approach is already well-validated. More complex alternatives (e.g., non-linear thermal time models, beta functions) add parameters without clear benefit for this MVP.

---

## 2. Canopy Photosynthesis

### Equations

**Eq A1** — CO₂ effect on light-use efficiency:
```
f_co2 = (1 - e^(-0.00305×CO2 - 0.222)) / (1 - e^(-0.00305×340 - 0.222))
ε = ε_ref × f_co2
```

**Eq A2** — CO₂ effect on Amax:
```
Amax_co2 = (49.57/34.26) × (1 - e^(-0.208×(CO2-60)/49.57))
Amax_leaf = Amax_ref × Amax_co2 × temp_reduction
```

**Eq A3** — Temperature reduction (piecewise linear):
```
Tav < 10:      0
10 ≤ Tav < 20: (Tav-10)/10
20 ≤ Tav ≤ 35: 1.0
Tav > 35:      max(0, 1-(Tav-35)/7)
```

**Eq A4** — Beer-Lambert light absorption:
```
f_int = 1 - e^(-KDF × LAI)
PAR_absorbed = radiation × 0.5 × f_int
```

**Eq A5** — Instantaneous PAR rate:
```
PAR_rate = PAR_absorbed / (daylength × 3.6)   [kJ/m²/d → J/m²/s]
```

**Eq A6** — Canopy Amax (Monsi-Saeki integral):
```
Amax_canopy = Amax_leaf × (1 - e^(-KDF × LAI)) / KDF
```

**Eq A7** — Saturating light response:
```
R_assim = Amax_canopy × (1 - e^(-ε × PAR_rate / Amax_canopy))
```

**Eq A8** — Daily total:
```
DTGA = R_assim × daylength × water_stress   [kg CO2/ha/d]
```

### Rationale
ORYZA2000 uses a 3-point Gaussian integration over both time-of-day and canopy depth, with separate sunlit/shaded leaf calculations (SGPCDT → SGPC1 → SGPL → SRDPRF). This is accurate but complex.

We use a **big-leaf model** (Eq A6) that analytically integrates the exponential light profile over canopy depth. This is derived from the Monsi-Saeki (1953) framework and is mathematically equivalent to `∫₀ᴸ Amax × e^(-K×l) dl`. It correctly captures the diminishing contribution of lower canopy layers without explicit multi-layer integration.

**Why this over the full Gaussian?** The big-leaf model produces ~650 kg CO₂/ha/d at LAI=4.3 under tropical conditions, within the ORYZA2000 typical range of 500-700. The simplification loses the sunlit/shaded distinction and the time-of-day radiation variation, but these primarily affect the shape of the daily assimilation curve, not its integral.

**Why not radiation-use efficiency (RUE)?** RUE models (biomass = ε × intercepted radiation) are simpler but lack the saturation response that is physiologically important at high light levels. The negative exponential (Eq A7) correctly saturates, which matters when comparing high-LAI and low-LAI scenarios.

---

## 3. Respiration and Growth

### Equations

**Eq R1** — Maintenance: `Rm = (WLVG×0.02 + WST×0.01 + WSO×0.003 + WRT×0.01) × T_eff × MNDVS`

**Eq R2** — Temperature: `T_eff = 2.0^((Tav - 25)/10)`

**Eq R3** — Metabolic age: `MNDVS = WLVG / max(1, WLVG + WLVD)`

**Eq R4** — Growth cost: `CRGCR = FSH×(1.326×FLV + 1.326×FST + 1.462×FSO) + 1.326×FRT`

**Eq R5** — Growth rate: `GCR = max(0, (DTGA × 30/44 - Rm) / CRGCR)`

### Rationale
Maintenance coefficients are from Penning de Vries et al. (1989) for tropical rice, validated in ORYZA2000. Growth conversion costs are from biochemical pathway analysis (Penning de Vries et al. 1974). These are well-established values used across most crop models (WOFOST, DSSAT, APSIM). No viable alternative exists at this level of simplification.

---

## 4. Leaf Area Index

### Equations

**Eq L1** — Exponential phase (LAI < 1.0): `gLAI = LAI × RGRL × HU`

**Eq L2** — SLA phase (LAI ≥ 1.0): `LAI = SLA(DVS) × WLVG`

**Eq L3** — Leaf death: `LLV = DRLV(DVS) × WLVG`

### Rationale
The dual-phase approach is directly from ORYZA2000 (Section 3.2.9, SUBLAI2 subroutine). The biological basis: young rice plants with open canopies expand leaf area as fast as temperature allows (exponential phase), because light is not limiting. Once leaves start overlapping (LAI ≈ 1.0), carbon supply becomes the constraint and LAI tracks leaf weight via SLA.

**Why not SLA-only?** At very low LAI, assimilation is tiny (Beer-Lambert: LAI=0.1 intercepts only 4% of light), so leaf weight barely grows, so LAI barely grows — a self-reinforcing trap. The exponential phase breaks this trap using the temperature-driven mechanism that is observed empirically (Kropff et al. 1994).

---

## 5. Sound-Response Engine

### Architecture

The sound-response engine has four layers plus two special pathways:

1. **Biophysical boundaries** — hard physics limits (activation floor, damage ceiling)
2. **Empirical windows** — 9 rice-specific + 2 cross-species frequency bands
3. **Lorentzian resonance** — per-outcome curve fitted to v2 CSV data
4. **Cross-species concordance** — data-driven discount replacing blanket 40%
5. **Ultrasound pathway** — direct lookup bypassing Lorentzian (≥20 kHz)
6. **Nearest-neighbour pathway** — for outcomes with oscillating data (iWUE)

### Equations

**Eq SE1** — Activation floor: `if SPL > 0 and SPL < 60 dB → effect = 0`

**Eq SE2** — Damage ceiling: `if SPL > 110 dB → effect = -15%, hard override`

**Eq SE3** — High-SPL caution zone (90-110 dB):
```
When no SPL-matched data exists in CSV:
  modifier = max(0.3, 1 - 0.5 × (SPL - 90) / (110 - 90))
When CSV has data measured at this SPL:
  modifier = 1.0 (data overrides the general caution)
```
Rationale: Studies like Jeong 2014 measured positive effects at 100 dB. The caution zone only applies when extrapolating beyond studied SPLs.

**Eq SE4** — Lorentzian resonance response:
```
R(f) = A / (1 + ((f - f₀) / γ)²)
```
where f₀ = resonant frequency (Hz), A = peak amplitude (%), γ = bandwidth (Hz).

**Eq SE5** — Lorentzian fitting (grid search per outcome):
```
For each outcome: collect {(freq_i, effect_i, weight_i)} from CSV
Grid search (f₀, A, γ) minimising: Σ wᵢ × (R(fᵢ) - effectᵢ)²
  f₀ range: [min(freq) - 100, max(freq) + 100]
  A range:  [1, min(200, 1.5 × max(|effect|) + 10)]
  γ range:  [40, max(400, 0.8 × (max(freq) - min(freq)))]
```

**Eq SE6** — SPL proximity weighting:
```
spl_w = max(0.1, 1 - |SPL_query - SPL_row| / 30)
```
Only applied when both query and row have known SPL values.

**Eq SE7** — Row confidence: `conf = confidence_class × completeness_class`
- confidence: high=1.0, medium=0.6, low=0.3
- completeness: full=1.0, partial=0.5, insufficient=0.2

**Eq SE8** — Cross-species concordance discount:
```
1. Evaluate rice Lorentzian at cross-species frequency: rice_pred = R(f)
2. If directions disagree: discount = 0.15
3. If directions agree and rice_pred > 0.5: discount = 0.5 + 0.5 × min(eff, rice_pred)/max(eff, rice_pred)
4. If rice predicts near-zero: discount = 0.6
```

**Eq SE9** — Distance-from-data confidence decay:
```
distance_penalty = max(0.1, 1 - min_dist / (γ × 3))
final_confidence = fit_confidence × distance_penalty
```

**Eq SE10** — Blend (CSV + windows):
```
If no window matches:     effect = csv_effect (100% CSV)
If csv_count ≥ 3:         effect = 0.85 × csv + 0.15 × window
If csv_count < 3:         effect = 0.65 × csv + 0.35 × window
If no CSV data:           effect = window_effect
```

**Eq SE11** — Dosage scaling (airborne only, not ultrasound):
```
hours < 1.0:   dose_factor = √hours
1.0 ≤ hours ≤ 4.0: dose_factor = 1.0  (within studied range)
hours > 4.0:   dose_factor = min(1.3, 1.0 + 0.1 × ln(hours/4))
Ultrasound (≥ 20 kHz): dose_factor = 1.0 always
```
Rationale: The Lorentzian already captures effects at the study's dosage. Scaling only applies when extrapolating beyond studied durations.

**Eq SE12** — Effect cap: `|effect| > 100% → clamp, confidence × 0.5`

**Eq SE13** — High-frequency confidence decay (5-15 kHz airborne):
```
f < 5 kHz:     modifier = 1.0
5 < f < 15 kHz: modifier = 1 - (f - 5000) / 10000
f ≥ 15 kHz:    modifier = 0 (no documented airborne effect)
f ≥ 20 kHz:    ultrasound pathway (bypasses all above)
```

### Water Status Sub-KPIs (CWSI Theory)

Water status is NOT fitted as a single Lorentzian. Per CWSI theory (Idso et al. 1981) and the 2016 Nature meta-analysis of 164 drought studies, gsw and RWC measure fundamentally different physiological variables on different timescales:

- **gsw** = active stomatal regulation (minutes timescale, most sensitive)
- **RWC** = passive hydration state (hours-days, smaller magnitude)
- **iWUE** = carbon fixed per water transpired (oscillating, non-Lorentzian)

Three independent Lorentzian tracks:

| Sub-KPI | Data source | Lorentzian | Peak |
|---------|-------------|------------|------|
| `water_status_gsw_drought` | Jeong 2014 drought gsw (4 points, 250-1500 Hz) | f₀≈1128 Hz, A≈89%, γ≈680 Hz | +89% at 1500 Hz |
| `water_status_gsw_wellwatered` | Jusoh 2023 gsw (5 points, 350-380 Hz) | Oscillating — mixed ±signs | Varies |
| `water_status_rwc` | Jeong 2014 RWC (5 points, 250-1500 Hz) | f₀≈1428 Hz, A≈37%, γ≈740 Hz | +37% at 1500 Hz |
| `water_status_efficiency` | Jusoh 2023 iWUE (5 points, 350-380 Hz) | **Nearest-neighbour** (not Lorentzian) | +62.7% at 359 Hz |

**Eq WS1** — Aggregate water_status (CWSI-weighted):
```
Under drought:      WS = 0.60 × gsw_drought + 0.25 × rwc + 0.15 × iWUE
Under well-watered: WS = 0.40 × gsw_wellwatered + 0.40 × iWUE + 0.20 × rwc
```
Weights reflect CWSI hierarchy: gsw is the primary operational water stress indicator.

### Sign Normalisation

Some outcomes are measured such that negative = beneficial:
```
FLIP_SIGN_OUTCOMES = {germination_time, sheath_blight_incidence, membrane_permeability}
```
These are flipped before Lorentzian fitting so the curve sees all beneficial effects as positive. The flip is reversed in the validation harness when comparing against YAML bands that use raw measurement convention.

### Ultrasound Pathway (≥ 20 kHz)

Ultrasound operates via liquid cavitation (micro-bubble formation), not airborne resonance. The Lorentzian is not physically appropriate. Instead:

**Eq US1** — Distance-weighted average of ultrasound data points:
```
For each matching row at f_row:
  prox = max(0.05, 1 - |f_query - f_row| / 30000)
  total_weight += conf × prox
  total_effect += effect × conf × prox
effect = total_effect / total_weight
```

Outcome-specific routing ensures germination queries use germination data only (Wang2022 germination_time) and vigor queries use vigor data only (Wang2020 seedling_dry_weight).

No dosage scaling for ultrasound — treatment durations are fundamentally different (minutes of liquid soaking vs hours of airborne exposure).

### Nearest-Neighbour Pathway (iWUE)

iWUE data oscillates wildly across 30 Hz (±60%) due to standing-wave artefacts. A Lorentzian cannot fit this pattern. Direct nearest-neighbour lookup:

**Eq NN1** — Return the effect of the closest data point:
```
effect = effect_of_nearest_freq_point
confidence = point_confidence × max(0.1, 1 - distance / 50)
```

### Rationale

**Why Lorentzian over distance-weighted averaging?** The Lorentzian is the physics-correct response of a damped harmonic oscillator, which matches PAFT theory (plants are resonant structures). Interpolation follows resonance physics, not arbitrary distance decay.

**Why separate water sub-KPIs?** Jeong 2014 shows gsw increases 2.4× more than RWC at the same frequency (89% vs 37% at 1500 Hz). Pooling them suppresses the gsw prediction. The CWSI framework (Idso 1981) establishes gsw as the primary water stress indicator. The Nature meta-analysis (2016, 164 studies) confirmed this hierarchy across all drought classes tested.

**Why nearest-neighbour for iWUE?** The Jusoh 2023 iWUE data oscillates -24.9%, +30.0%, +29.5%, +62.7%, -32.7% across just 30 Hz. This is not a smooth resonance response — it reflects standing-wave interference in the growth chamber. No smooth curve can fit this pattern. Nearest-neighbour returns the actual measured value at the closest frequency.

**Why bypass Lorentzian for ultrasound?** Ultrasound seed treatment operates via cavitation (20-60 kHz), not airborne resonance. The physics is fundamentally different: micro-bubbles crack seed coats, not mechanosensory ion channels. Applying a Lorentzian fitted to audible-range data would be physically incorrect.

---

## 6. Acoustics

**Eq AC1** — `Lp(r) = Lp(1m) - 20×log₁₀(r)` (inverse-square law)

**Eq AC2** — `d = √(Δx² + Δy² + Δz²)` (3D distance)

**Eq AC3** — `Lp_total = 10×log₁₀(Σ 10^(Lp_i/10))` (incoherent combination)

### Rationale
Free-field inverse-square is the simplest physically correct model for outdoor sound propagation. Real fields have wind, atmospheric absorption, ground effects, and foliage — all ignored in v1. This means the model overestimates coverage uniformity and underestimates distance losses at high frequencies.

---

## 7. Equation Numbering Cross-Reference

| This spec | ORYZA2000 | Description |
|-----------|-----------|-------------|
| P1 | Eq 3.2 | Hourly temperature |
| P2-P3 | Eq 3.3-3.4 | Heat units |
| A1 | Eq 3.5 | CO₂ × efficiency |
| A2 | Eq 3.7 | CO₂ × Amax |
| A7 | Eq 3.6 | Light response curve |
| R1 | Eq 3.25 | Maintenance respiration |
| R2 | Eq 3.26-3.27 | Temperature effect |
| R5 | Eq 3.28 | Growth rate |
| L1 | Eq 3.38 | Exponential LAI |
| L2 | Eq 3.40 | SLA-based LAI |
| S1 | Eq 3.33 | Cold sterility |
| S2 | Eq 3.34 | Heat sterility |
| A6 | — | Monsi-Saeki (not in ORYZA2000) |
| SE1-SE13 | — | Sound engine (no ORYZA equivalent) |
| AC1-AC3 | — | Acoustics (standard physics) |

---

## 8. UI Frequency Slider (Hybrid Linear/Linear)

The frequency slider maps a 0-1000 integer range to 50 Hz – 100 kHz using two linear regions:

**Eq UI1** — Audible region (slider 0-750, covers 50 Hz to 10 kHz):
```
Hz = round(50 + (slider / 750) × 9950)
```
Step size: ~13.3 Hz per slider unit. Fine control for the range where most rice evidence exists.

**Eq UI2** — Extended region (slider 750-1000, covers 10 kHz to 100 kHz):
```
Hz = round(10000 + ((slider - 750) / 250) × 90000)
```
Step size: ~360 Hz per slider unit. Covers ultrasound seed treatment frequencies in simple increments.

### Rationale
Most rice phytoacoustic evidence is in 50-5000 Hz. A single linear slider to 100 kHz would compress this into the first 5%. The two-region approach gives 75% of the slider to the evidence-dense audible range and 25% to the ultrasound range. Both regions are plain linear — no logarithmic math needed. The upper region naturally provides ~10 kHz resolution at key ultrasound frequencies (20 kHz, 35 kHz, 40 kHz, 43.5 kHz) used in seed treatment studies.

**Ultrasound boundary:** ≥20 kHz (slider position ~778). UI turns purple with ⚡ indicator. Audio preview caps at ~20 kHz (browser limit).

---

## 9. Dual-Phase Sound Treatment

The UI supports different sound parameters for germination (days 1-10) and seedling-to-harvest (days 11+):

**When total days ≤ 10:** Only the germination panel is active. The simulation runs with a single set of sound parameters at the seed/germination stage.

**When total days > 10:** Both panels are active. The simulation currently uses the seedling-to-harvest parameters for the main growth simulation (the germination phase parameters are passed to the API for future dual-phase integration). This is a known simplification — the backend currently runs a single sound driver for the entire simulation rather than switching parameters at day 10.

### Rationale for 10-day boundary
Rice germination typically completes within 5-10 days under tropical conditions. Bochu et al. (2003) used 2-day exposures at 400 Hz/106 dB for seed treatment. The 10-day boundary is a reasonable approximation for when the plant transitions from seed/germination stage to seedling stage, though in practice this varies with temperature and variety.
