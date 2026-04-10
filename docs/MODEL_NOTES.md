# Model Notes

## How This Simulator Works

This simulator combines two independent layers:

1. A **rice growth engine** (simplified from ORYZA2000) that models baseline crop development.
2. A **sound-response engine (v4)** that estimates how pure-tone sound treatment modifies rice physiology.

Neither layer is a fully predictive biophysical model. The rice growth engine is a documented simplification of ORYZA2000. The sound-response engine is a **literature-guided hypothesis generator** using a mechanistic Gaussian Process architecture, not a dose-response curve or mechanistic model.

---

## Rice Growth Engine

### Scientific basis
- **Primary source:** Bouman et al. (2001), *ORYZA2000: modeling lowland rice*
- **Drought extensions:** Li et al. (2017), *From ORYZA2000 to ORYZA (v3)*
- **Architecture pattern:** WOFOST/PCSE (modular process classes, rate-then-integrate, separate parameters)

### What was implemented from ORYZA2000

| Process | ORYZA2000 reference | Implementation |
|---------|---------------------|----------------|
| Phenology | Section 3.2.1, Eqs 3.2-3.4 | Hourly heat units with bilinear temperature response |
| Leaf area growth | Section 3.2.9, Eqs 3.37-3.38 | Dual-phase: exponential when LAI<1.0, SLA-based after |
| Canopy photosynthesis | Section 3.2.2, Eqs 3.5-3.7 | Big-leaf model (analytical Monsi-Saeki integral) |
| Partitioning | Section 3.2.3, Eq 3.21 | DVS-dependent tables, drought shifts to roots |
| Respiration | Section 3.2.5, Eqs 3.25-3.28 | Tropical rice maintenance + growth coefficients |
| Spikelet sterility | Section 3.2.8, Eqs 3.31-3.34 | Cold/heat sterility |
| Drought stress | Li et al. (2017) Section 2.3.2 | DVS-duration proxy (not soil water balance) |

### Key simplifications
- **Canopy photosynthesis:** Big-leaf approximation. Output ~650 kg CO2/ha/d at LAI=4.3 (ORYZA2000: 500–700).
- **Leaf area transition:** Immediate switch at LAI=1.0 instead of ORYZA2000's smooth blend.
- **Drought stress:** Regime-based factor, not soil water balance.
- **Nitrogen/photoperiod:** Not modelled.

### Dual-phase sound drivers (v4)
`DualPhaseSoundParams` holds two independent `SoundDriver` instances:
- **GerminationSoundDriver** — applied for days 0 to `germination_days` (default 10). Stage bucket `"germination"`. Typically higher SPL for seed coat penetration (Bochu 2003: 400 Hz / 106 dB).
- **VegetativeSoundDriver** — applied from `germination_days` onward. Stage bucket `"seedling"` or `"vegetative"`. Moderate SPL for stomatal and photosynthetic effects (Jusoh 2023: 350 Hz / 68 dB).

`run_simulation_dual_phase()` is the v4 entry point in `rice_core.py`. It selects the correct driver each day. `api.py` reads separate germination parameters from the request and builds both drivers before calling this function.

---

## Sound-Response Engine (v4 — Mechanistic GP)

v4 adds mechanotransduction latent variables and biological pathway activations as the GP prior mean, replacing the discrete empirical windows used in v3.

### Architecture (four layers + two special pathways)

#### Layer 1: Biophysical Boundaries (unchanged)
- **Activation floor:** 60 dB — minimum SPL to trigger mechanosensory response
- **Damage ceiling:** 110 dB — cell membrane damage, hard override to -15%
- **High-SPL broadband inhibition (90–110 dB):** When no SPL-matched positive rice data exist within ±150 Hz of the query, a broadband inhibitory modifier is applied: `modifier = -0.3 - 0.2 × frac` where `frac = (SPL - 90) / 20`. At 95 dB: -35%. This encodes the Grass_General_90dB finding (-40% vegetative growth at >90 dB continuous broadband exposure) as a biophysical boundary rather than a GP training point.
- **High-frequency decay:** 5–15 kHz confidence decay for airborne sound

#### Layer 2: Mechanotransduction Latent S_mech,air(f, SPL)

  S_mech,air(f, SPL) = S_f(f) × S_SPL(SPL) × S_air(f)

- **S_f(f):** Normalised output of the **MechActivation GP** — a Matérn ν=5/2 GP trained on pooled stress/activation outcomes (Fv/Fm, flavonoid content, sheath blight incidence, stomatal conductance, photosynthesis rate, chlorophyll content, SOD/CAT/POD enzymes). Normalised to [0,1] via 95th-pct Y_max over 50–5000 Hz. Captures the frequency selectivity of the plant's mechanosensory response.
- **S_SPL(SPL):** Piecewise-linear gate: 0 below 60 dB, linear ramp to 1 at 85 dB, plateau 85–100 dB, linear decay to 0 at 110 dB.
- **S_air(f):** Airborne coupling decay: 1 below 5 kHz, linear decay to 0 at 15 kHz.

S_mech,air = 0 for ultrasound (≥20 kHz); cavitation pathway applies instead.

#### Layer 3: Pathway Activations + Mechanistic Prior

Five biological pathways replace the discrete empirical windows as the GP prior mean:

| Pathway | Centre (Hz) | σ (log-Hz) | Stages | Regime | min SPL | Peak effect |
|---------|-------------|------------|--------|--------|---------|-------------|
| seed_vigor | ~380 | 0.30 | seed, germination, seedling | any | 0 dB | +25% |
| leaf_gas_exchange | ~357 | 0.12 | seedling, vegetative | any/well-watered | 0 dB | +40% |
| drought_resilience | ~800 | 0.55 | vegetative | any/drought | 0 dB | +45% |
| membrane_damage | ~4000 | 0.25 | all | any | 100 dB | -20% |
| ultrasound_cavitation | ~38000 | 0.18 | seed | any | 0 dB | +20% |

Each pathway activation uses a log-Gaussian frequency shape gated by stage, regime, SPL, and scaled by S_mech,air. Pathway activations are combined via `OUTCOME_PATHWAY_WEIGHTS` (e.g. photosynthesis: 80% leaf_gas_exchange + 20% seed_vigor) to form the per-outcome mechanistic prior.

**Hybrid prior floor:** When the mechanistic prior is near-zero (< 1%) but a direct empirical window exists at the query frequency, the window midpoint is used as the prior. This preserves direct literature evidence that pathway geometry may underweight — e.g. Hou 2009 yield at 550 Hz (window `rice_yield_paft_400_700`), which is 3.6σ from the leaf_gas_exchange pathway centre.

#### Layer 4: Gaussian Process Regression (unchanged from v3.1)

Each study in `rice_sound_master_v2.csv` contributes a training point `(frequency_hz, effect_pct)` with per-study noise variance α derived from `confidence_weight`, `completeness_class`, SPL proximity, and crop species. The GP:

- Treats every study fairly regardless of frequency
- Returns honest high uncertainty in sparse regions
- Is fitted to **residuals** (observed − hybrid prior), so far from data it returns to the prior
- Handles cross-species data via higher noise variance (α × 2.78)

**Kernel:** Matérn ν=5/2, length scale [100, 1200] Hz, amplitude initialised from residual RMS.
**Frequency aggregation:** Multiple studies at the same frequency are precision-weighted into one cluster before fitting, preventing pseudo-replication from narrowing the kernel.
**Local density confidence:** Confidence counts only training clusters within one fitted length scale of the query. A study at 800 Hz earns credit near 800 Hz only.

#### Special pathways (unchanged from v3.1)

- **Ultrasound (≥20 kHz):** Liquid cavitation mechanism. Distance-weighted direct lookup, bypasses GP and S_mech entirely.
- **iWUE (water_status_efficiency):** Nearest-neighbour with 50 Hz decay. The Jusoh 2023 iWUE data oscillates ±60% over 30 Hz due to standing-wave chamber artefacts. No smooth kernel can represent this. Aggregation is also skipped for iWUE — the oscillation is real signal that nearest-neighbour relies on.

### Sign normalisation

Outcomes where a MORE NEGATIVE raw measured value is BENEFICIAL are flipped before GP fitting (flip: germination_time, germination_time_h, speed_of_germination). `sheath_blight_incidence` is **NOT flipped** because `rice_sound_master_v2.csv` already stores it as +50 meaning "+50% disease reduction" (beneficial direction). Flipping would incorrectly make 550 Hz look harmful. This is documented in the CSV `notes` column and in the code comment for `_FLIP_SIGN_OUTCOMES`. If the CSV convention changes to raw incidence values (higher = worse), sheath_blight_incidence must be added back to `_FLIP_SIGN_OUTCOMES`.

### Why the Lorentzian was replaced (v3.0), what was fixed in v3.1, and the v4 upgrade

**Lorentzian (v2):** Imposed a bell-curve shape. With 39/46 usable rows in 250–499 Hz, the curve peaked at 350–400 Hz and produced artefactually confident tails at untested adjacent frequencies.

**GP v3.0:** Replaced Lorentzian with a GP using empirical windows as prior mean. Fixed honest uncertainty. Retained bugs: double sign-flip on sheath blight; `GrassGeneral_90dB` distorting the vigor kernel; global (not local) density confidence.

**GP v3.1:** Fixed pseudo-replication (frequency aggregation), local density confidence, sheath blight sign, and `GrassGeneral_90dB` exclusion from GP training (it belongs to the biophysical boundary layer — see below).

**v4:** Replaced empirical-window prior with mechanotransduction latent + five biological pathways. Added dual-phase sound drivers. Added MechActivation GP. Added hybrid prior floor for direct-evidence window preservation.

### Water status: three independent GPs (CWSI theory, unchanged)

Per Idso et al. (1981) and the 2016 Nature meta-analysis:

| Sub-KPI | Data | Method |
|---------|------|--------|
| gsw_drought | Jeong 2014 (250–1500 Hz) | GP |
| gsw_wellwatered | Jusoh 2023 (350–380 Hz) | GP |
| water_status_rwc | Jeong 2014 (250–1500 Hz) | GP |
| water_status_efficiency (iWUE) | Jusoh 2023 (350–380 Hz) | Nearest-neighbour |

Aggregate `water_status` CWSI-weighted: drought = 60% gsw_drought + 25% rwc + 15% iWUE; well-watered = 40% gsw_wellwatered + 40% iWUE + 20% rwc.

### Honest limitations

1. **46 usable data points** across all KPIs. GP is well-constrained only near data; uncertainty grows rapidly elsewhere — now reported honestly via local density confidence.
2. **Pathway parameters** (mu_log_hz, sigma_log_hz, peak_pct) are set from literature anchors, not optimised. They are prior beliefs, not fitted parameters.
3. **No published equation** predicts how a given frequency affects rice from first principles. This engine is principled interpolation — more than a lookup table, less than a predictive model.
4. **iWUE** results are highly local (±50 Hz). Do not extrapolate.

---

## Data sources

### Rice-specific CSV data
- Jusoh & Ramlee 2023 (350–380 Hz seedling assimilation/height/iWUE/gsw)
- Qi et al. 2010 (550 Hz chlorophyll)
- Hou et al. 2009 (550 Hz PAFT yield)
- Bochu et al. 2003 (400 Hz/106 dB germination, 4000 Hz/111 dB injury)
- Jeong et al. 2014 (250–1500 Hz drought gsw + RWC — 10 high-quality rice points)
- Jeong et al. 2008 (50–250 Hz gene expression)
- Hassan et al. 2014 (800–1500 Hz drought)
- Munasinghe et al. 2023 (350 Hz iWUE, weak evidence)
- Kim et al. 2019 (250/800/1000 Hz flavonoid content)
- Sri Lankan rice 2024 (3–5 kHz vegetative)
- Ultrasound: Wang 2022, Wang 2020, UltrasonicSoak 2025, UltrasoundGerm 2024

### Cross-species data (higher noise variance in GP)
- Wheat, corn, oat, barley, sorghum, general grasses

### Boundary-layer evidence (excluded from GP training)
- **GrassGeneral_90dB** (-40% vegetative growth at >90 dB continuous broadband exposure, grass). This study justifies `GRASS_INHIBITION_DB = 90 dB` and the broadband high-SPL inhibitory modifier in `_check_biophysical_boundaries()`. Including it as a GP training point would shrink the vigor kernel length scale around 90 dB/low Hz, underweighting evidence at 800–1500 Hz from Jeong 2014. It is stored in `_BOUNDARY_LAYER_STUDIES` in `sound_response.py`.

---

## WOFOST/PCSE Architecture

- **Process classes:** Phenology, WaterStress, Assimilation, Respiration, Partitioning, LeafDynamics, Spikelets
- **Rate-then-integrate:** DailyRates computed from current state before integration
- **CultivarParams:** Separate dataclass, loadable from YAML (data/IR72.yaml)
- **SoundDriver / DualPhaseSoundParams:** Sound treatment as formal drivers modifying biological rates
- **RiceEngine:** Orchestrator class

---

## Theoretical Framework

Three phytoacoustic theories inform the model:

**PAFT Resonance (350–800 Hz):** Plants are damped mechanical resonators. Optimal frequency induces cellular resonance. The leaf_gas_exchange pathway (357 Hz centre) and seed_vigor pathway (380 Hz) encode this. The GP's Matérn kernel learns the effective frequency selectivity from data without imposing a bell curve.

**Mechanotransduction Ca²⁺ Pathway:** Acoustic waves deform plasma membranes → mechanosensitive ion channels open → Ca²⁺ influx → CDPKs → H⁺-ATPases → stomatal changes → enhanced photosynthesis. The drought_resilience pathway (800 Hz, Jeong 2014) captures the drought-stress gene expression and stomatal response driven by this pathway at higher frequencies.

**Ultrasound Cavitation (≥20 kHz):** Micro-bubbles formed by cavitation crack seed coats, increase hydration diffusion, enhance alpha-amylase. Entirely separate from airborne resonance. Handled by dedicated nearest-neighbour lookup in `_ultrasound_lookup()`.
