# Model Notes

## How This Simulator Works

This simulator combines two independent layers:
1. A **rice growth engine** (simplified from ORYZA2000) that models baseline crop development.
2. A **sound-response engine** that estimates how pure-tone sound treatment modifies rice physiology.

Neither layer is a fully predictive biophysical model. The rice growth engine is a documented simplification of ORYZA2000. The sound-response engine is a **literature-guided hypothesis generator** using empirical windowing, not a dose-response curve or mechanistic model.

---

## Rice Growth Engine

### Scientific basis
- **Primary source:** Bouman et al. (2001), *ORYZA2000: modeling lowland rice*
- **Drought extensions:** Li et al. (2017), *From ORYZA2000 to ORYZA (v3)*
- **Architecture pattern:** WOFOST/PCSE (modular process classes, rate-then-integrate, separate parameters)

### What was implemented from ORYZA2000

| Process | ORYZA2000 reference | Implementation |
|---------|---------------------|----------------|
| Phenology | Section 3.2.1, Eqs 3.2-3.4 (SUBDD, PHENOL) | Hourly heat units with bilinear temperature response |
| Leaf area growth | Section 3.2.9, Eqs 3.37-3.38 (SUBLAI2) | Dual-phase: exponential when LAI<1.0, SLA-based after |
| Canopy photosynthesis | Section 3.2.2, Eqs 3.5-3.7 | Big-leaf model (analytical Monsi-Saeki integral) |
| Partitioning | Section 3.2.3, Eq 3.21 | DVS-dependent tables, drought shifts to roots |
| Respiration | Section 3.2.5, Eqs 3.25-3.28 | Tropical rice maintenance + growth coefficients |
| Spikelet sterility | Section 3.2.8, Eqs 3.31-3.34 | Cold/heat sterility |
| Drought stress | Li et al. (2017) Section 2.3.2 | DVS-duration proxy (not soil water balance) |

### Key simplifications
- **Canopy photosynthesis:** Big-leaf approximation instead of ORYZA2000's multi-layer Gaussian integration. Uses `Amax_canopy = Amax_leaf × (1-exp(-K×LAI))/K`. Output checked: ~650 kg CO2/ha/d at LAI=4.3 (ORYZA2000 typical: 500-700).
- **Leaf area transition:** Immediate switch at LAI=1.0 instead of ORYZA2000's smooth weighted blend.
- **Drought stress:** Regime-based factor using exponential decay over DVS-duration, not soil water balance. The factor of 2.0 in the exponent is arbitrary.
- **Nitrogen:** Not modelled. RGRL uses RGRLMX (0.0085 (°Cd)⁻¹) always.
- **Photoperiod:** Not modelled. Assumes optimal.
- **Leaf-weight reconciliation during exponential phase:** When exponential LAI growth implies more leaf weight than assimilate-driven growth produced, leaf weight is adjusted upward. This is an ad hoc fix, not from ORYZA2000.

---

## Sound-Response Engine

### What it is NOT

This engine does **not**:
- Fit a continuous dose-response curve across frequencies
- Build a mechanistic biophysical model of membrane mechanotransduction
- Predict from first principles what happens at untested frequencies

### What it IS: Lorentzian Resonance + CWSI Water Sub-KPIs

The engine uses a 4-layer architecture with two special pathways:

#### Layer 1: Biophysical Boundaries (hard physics)
- **Activation floor:** 60 dB — minimum SPL to trigger mechanosensory response
- **Damage ceiling:** 110 dB — cell membrane damage, hard override to -15%
- **High-SPL caution:** 90-110 dB — soft modifier when NO SPL-matched data exists in CSV. When studies measured effects at this SPL (e.g., Jeong 2014 at 100 dB), the data overrides the caution zone.
- **High-frequency decay:** 5-15 kHz confidence decay for airborne sound

#### Layer 2: Empirical Windows (9 rice-specific + 2 cross-species bands)
Discrete frequency×SPL×stage bands providing directional guidance. These are blended with the Lorentzian but the Lorentzian takes priority when it has ≥3 data points.

#### Layer 3: Lorentzian Resonance (per-outcome curve fitted to v2 CSV)
Each KPI has its own Lorentzian `R(f) = A / (1 + ((f - f₀) / γ)²)` fitted by grid search to the `rice_sound_master_v2.csv` dataset (76 rows, 47 usable). The grid search covers the full data frequency range with adaptive amplitude and bandwidth limits.

#### Layer 4: Cross-Species Concordance
Non-rice data is scored by how well it agrees with the rice Lorentzian trend. Concordant data gets 0.7-0.96 confidence; discordant data gets 0.15.

#### Special: Ultrasound Pathway (≥ 20 kHz)
Ultrasound seed treatment operates via liquid cavitation, not airborne resonance. The Lorentzian is bypassed. Distance-weighted direct lookup from matching ultrasound rows. No dosage scaling (different treatment duration profiles). Outcome-specific routing: germination data for germination queries, vigor data for vigor queries.

#### Special: Nearest-Neighbour for iWUE
Intrinsic water use efficiency oscillates wildly across 30 Hz (±60%) due to standing-wave artefacts. A Lorentzian cannot fit this. Direct nearest-neighbour returns the actual measured value.

### Water Status: Three Independent Tracks (CWSI Theory)

Per Crop Water Stress Index theory (Idso et al. 1981) and the 2016 Nature meta-analysis of 164 drought studies:

| Track | Physiology | Data | Timescale |
|-------|-----------|------|-----------|
| `gsw_drought` | Active stomatal regulation | Jeong 2014 (250-1500 Hz) | Minutes |
| `gsw_wellwatered` | Stomatal opening | Jusoh 2023 (350-380 Hz) | Minutes |
| `rwc` | Passive hydration state | Jeong 2014 (250-1500 Hz) | Hours-days |
| `iWUE` | Carbon per water | Jusoh 2023 (350-380 Hz) | Nearest-neighbour |

These are NOT pooled. Each has its own Lorentzian. The 2.4× magnitude difference between gsw (+89%) and RWC (+37%) at 1500 Hz (Jeong 2014) confirms they must be separate.

Aggregate `water_status` is derived with CWSI-weighted combination:
- Drought: gsw_drought 60%, RWC 25%, iWUE 15%
- Well-watered: gsw_wellwatered 40%, iWUE 40%, RWC 20%

### Sign Normalisation

Outcomes where negative = beneficial (germination_time, sheath_blight_incidence, membrane_permeability) are sign-flipped before Lorentzian fitting so the curve sees all beneficial effects as positive.

### Honest Limitations

1. **The Lorentzian is a smoothing approximation.** It cannot capture local spikes within 30 Hz ranges (e.g., 21.1% height at 357 Hz surrounded by 8-11% at adjacent frequencies).
2. **47 usable data points** across all KPIs. Most outcomes have 3-5 points. The Lorentzian is well-constrained only near data; extrapolation decays rapidly.
3. **No published equation exists** that can predict how a given sound frequency will affect rice. This simulator is a principled interpolation engine grounded in resonance physics — more than a lookup table, less than a predictive model.

---

## Data Sources

### Rice-specific CSV data
- Jusoh & Ramlee 2023 (350-380 Hz seedling assimilation/height)
- Qi et al. 2010 (550 Hz chlorophyll)
- Hou et al. 2009 (550 Hz PAFT yield)
- Bochu et al. 2003 (400 Hz/106 dB germination, 4000 Hz/111 dB injury)
- Jeong et al. 2008 (50-250 Hz gene expression, 800 Hz drought RWC)
- Hassan et al. 2014 (800-1500 Hz drought)
- Munasinghe et al. 2023 (350 Hz iWUE, weak evidence)
- Sri Lankan rice 2024 (3-5 kHz vegetative)
- Ultrasound studies: Wang 2022, Wang 2020, UltrasonicSoak 2025, UltrasoundGerm 2024

### Cross-species CSV data (40% confidence discount)
- Wheat: 300 Hz, 1250 Hz (no effect), 3-5 kHz germination, 5 kHz tillering, 12 kHz (no effect)
- Corn: 100-300 Hz phonotropism, 300 Hz germination, 4500 Hz circadian
- Oat: 300 Hz, 600 Hz (root inhibition)
- Barley: 43.5 kHz ultrasound germination
- Sorghum: >20 kHz ultrasound seed coat
- General grasses: 70 dB height optimum, >90 dB growth inhibition, 100 Hz/92 dB vegetative, 1000 Hz defense shift

---

## WOFOST/PCSE Architecture

The rice growth engine follows PCSE patterns:
- **Process classes:** Phenology, WaterStress, Assimilation, Respiration, Partitioning, LeafDynamics, Spikelets
- **Rate-then-integrate:** DailyRates computed from current state before any integration
- **CultivarParams:** Separate dataclass, loadable from YAML (data/IR72.yaml)
- **SoundDriver:** Sound treatment as a formal driver modifying biological rates during simulation
- **RiceEngine:** Orchestrator class coordinating all processes

---

## User Interface Design Decisions

### Dual-phase sound input
The UI separates sound treatment into germination (days 1-10) and seedling-to-harvest (days 11+) because the literature shows different optimal parameters for each phase. Seed-stage treatments use higher SPL (95-110 dB for seed coat penetration, Bochu 2003) while vegetative-stage treatments use lower SPL (60-85 dB for stomatal/photosynthetic effects, Jusoh 2023). The dual-phase UI lets users configure both independently.

**Known limitation:** The backend currently uses only one set of sound parameters per simulation run. Dual-phase backend support (switching sound driver at day 10) is planned but not yet implemented.

### Hybrid frequency slider
The slider uses two linear regions: 50 Hz–10 kHz (75% of range, ~13 Hz steps) and 10 kHz–100 kHz (25%, ~360 Hz steps). This gives fine control over the evidence-dense audible range while covering ultrasound seed treatment frequencies. The purple colour and ⚡ indicator above 20 kHz visually distinguish ultrasound from audible frequencies. Audio preview caps at ~20 kHz (browser limit).

### Days-first input ordering
Total days is the first input because it determines the UI state: ≤10 days shows only the germination panel, >10 days shows both panels. This prevents users from configuring seedling-to-harvest parameters when only germination days are selected.

---

## Theoretical Framework

The sound-response engine is informed by three distinct theories from the phytoacoustics literature. The Lorentzian resonance model is the mathematical unification of the first two; the third operates through a separate physical pathway.

### Theory 1: PAFT Resonance (explains WHERE effects peak)

Plant Acoustic Frequency Technology (Hou Tianzhen) posits that plants emit their own acoustic signals — from internal fluid dynamics and xylem cavitation — and that applying external sound matching these internal frequencies induces cellular resonance. The plant is a damped mechanical resonator with natural frequencies determined by its physical dimensions, water content, and tissue stiffness.

The Lorentzian curve `R(f) = A / (1 + ((f - f₀) / γ)²)` is the physics-correct response function for a damped harmonic oscillator. When we fit it to rice data, the peak at f₀ ≈ 345 Hz may correspond to the natural resonant frequency of rice seedling leaf structures. The bandwidth γ ≈ 110 Hz reflects how damped the system is — a narrow γ means the plant is a selective resonator, a wide γ means it responds broadly.

**Supporting studies:** Jusoh 2023 (peak assimilation response at 350-380 Hz), Qi 2010 (chlorophyll at 550 Hz — on the falling edge of the resonance), Hou 2009 (PAFT field trials at 550 Hz midpoint), Sri Lankan rice (3-5 kHz — possibly a secondary resonance mode).

### Theory 2: Mechanotransduction / Ca²⁺ Pathway (explains HOW effects work)

This is the cellular mechanism that explains why resonance produces biological changes:

1. Acoustic waves create mechanical stress → deform the plasma membrane
2. Membrane tension opens mechanosensitive (MS) ion channels
3. Rapid influx of cytosolic calcium (Ca²⁺) occurs
4. Ca²⁺ activates calcium-dependent protein kinases (CDPKs)
5. CDPKs activate plasma membrane H⁺-ATPases
6. Altered osmotic balance in guard cells → changes in stomatal conductance
7. Downstream: enhanced water use efficiency, increased carbon assimilation, stress gene expression

Blocking Ca²⁺ channels in sound-treated cells neutralises beneficial effects (confirmed experimentally), proving this pathway is necessary.

**Why the Lorentzian is compatible:** Mechanical membrane stress is maximised when the forcing frequency matches the tissue's natural frequency — at resonance. The Lorentzian IS the membrane stress response curve. The peak of the Lorentzian corresponds to maximum Ca²⁺ influx; the tails correspond to frequencies where membrane deflection is insufficient to open MS channels.

**Supporting studies:** Jeong 2008 (gene expression at 50-250 Hz, drought stress prep), Arabidopsis CML38 calcium-binding gene upregulation at 1000 Hz, duckweed photosynthesis-antenna protein genes at 100-500 Hz.

### Theory 3: Ultrasound Cavitation (separate mechanism)

Completely distinct from airborne resonance. Ultrasound (≥20 kHz) applied to seeds in liquid creates micro-cavitation bubbles that:
- Physically micro-fracture the seed coat without damaging the embryo
- Increase water hydration diffusion through the seed coat
- Enhance alpha-amylase enzymatic activity (starch breakdown)
- Weaken cell wall rigidity

This is brute-force mechanical disruption, not frequency-selective resonance. The optimal parameters (frequency, duration, liquid medium) are determined by cavitation physics, not plant resonance. Our model handles ultrasound through a separate pathway that bypasses the Lorentzian.

**Supporting studies:** Wang 2022 (40 kHz rice germination time -50%), Wang 2020 (35 kHz seedling dry weight +12%), barley 43.5 kHz (100% germination), sorghum >20 kHz (seed coat fracture).

### How the model unifies these theories

| Theory | Frequency range | Mechanism | Model component |
|--------|----------------|-----------|-----------------|
| PAFT Resonance | 50-5000 Hz | Structural resonance | Lorentzian curve (f₀, A, γ) |
| Mechanotransduction | 50-5000 Hz | Ca²⁺ ion channel cascade | Empirical windows (direction, magnitude bounds) |
| Ultrasound cavitation | ≥20 kHz | Liquid micro-cavitation | Separate ultrasound pathway |
| Biophysical limits | All | Membrane physics | Hard boundaries (60 dB floor, 110 dB ceiling) |

The Lorentzian resonance model is the mathematical shape of theories 1+2 combined. Theory 1 predicts WHERE the effect peaks (the resonant frequency f₀). Theory 2 predicts WHAT happens at the cellular level (the Ca²⁺ cascade). The Lorentzian amplitude A and bandwidth γ are empirically fitted to published rice data, not derived from first-principles membrane mechanics (those coupling constants don't exist yet in the literature).

---

# Assumptions Used in Implementation

## Rice Biology
- Phenological development uses ORYZA2000 heat unit approach with cardinal temperatures 8/30/42°C (Bouman et al. 2001, Section 3.2.1). Photoperiod sensitivity omitted in v1.
- Canopy photosynthesis uses a big-leaf model: `Amax_canopy = Amax_leaf × (1 - exp(-K×LAI)) / K`. This is standard canopy radiation theory (Monsi-Saeki framework), not an ORYZA2000 formula. ORYZA2000 uses multi-layer Gaussian integration. Output checked for plausibility (~650 kg CO2/ha/d at LAI=4.3).
- Leaf area growth follows ORYZA2000's dual-phase approach (Section 3.2.9): exponential temperature-driven growth when LAI < 1.0 (Eq 3.37-3.38, RGRL=0.0085 (°Cd)⁻¹), switching to SLA-based when LAI ≥ 1.0 (Eq 3.40).
- Drought stress uses a DVS-duration proxy inspired by Li et al. (2017) DTF concept, not a soil water balance. The factor of 2.0 in the exponent is arbitrary.

## Sound Response
- The engine uses **empirical windowing**, not continuous dose-response fitting. This matches the PAFT literature approach.
- Biophysical boundaries (60 dB floor, 110 dB ceiling, 90 dB inhibition) are hard constraints from cross-species membrane physics, not species-specific tuning.
- Empirical windows are discrete frequency×SPL×stage bands, not smooth curves. There are gap regions where the engine returns zero.
- Cross-species data (wheat, corn, oat, barley, sorghum, grasses) is used at 40% confidence discount for directional support only.
- The confidence score measures proximity to evidence, not prediction accuracy.
- Interpolation between data points is a mathematical artefact of distance weighting, not a biological prediction.
- When SPL exceeds the damage ceiling (110 dB), a hard inhibitory result is returned regardless of what windows or CSV data suggest.
- Seed-stage high-SPL germination treatments (95-110 dB) are exempted from the general grass inhibition rule (>90 dB) because they use a different mechanism (seed coat penetration vs chronic exposure).

## Acoustics
- Free-field inverse-distance law: Lp(r) = Lp(1m) - 20×log10(r).
- Ignores: wind, atmospheric absorption, ground reflection, foliage scattering, speaker directivity.

## KPIs
- Yield Index: 0-100 scale, 50 = no change.
- Water Index: baseline stress score + sound improvement.
- Stress Resilience Index: weighted composite of photosynthesis (40%), oxidative defense (30%), vigor (30%).
- Coverage Quality Index: coverage >60dB + >70dB + uniformity + dosage.

---

## RREC Farm Baseline Integration

### Data source
USDA-ARS Dale Bumpers National Rice Research Center, Stuttgart AR. 2022 growing season. Published as Farag et al. (2024), Dryad doi:10.5061/dryad.v41ns1s4z.

### What the data provides
- **152 plots** with measured yield (metric tons/hectare) across 21 rice cultivars
- **Phenology dates**: emergence (DOY 132-133), heading (68-98 days from emergence)
- **Daily weather**: 136 days of Tmin, Tmax (°F→°C converted), solar radiation (MJ/m²→kJ/m²)
- **Five experiments** grouped into three user-facing categories:
  - Hybrid (RT7521FP): median 13.7 t/ha — commercial high-yield hybrid
  - Inbred (Jeff, Santee Gold): median 7.95 t/ha — standard breeding lines
  - Traditional (Tiara, Mini Core): median 6.3 t/ha — diverse/heritage varieties

### How it's used
1. **Real weather replaces synthetic weather.** The simulation is driven by actual 2022 Stuttgart AR temperatures and solar radiation instead of sinusoidal approximations. This improves realism of the growth trajectory.

2. **RREC yields serve as the benchmark for comparison.** The UI shows three horizontal bars: RREC farm yield (real), simulator baseline (ORYZA with real weather, no sound), and sound-treated (ORYZA with sound driver). The user sees the sound effect as a delta against real farm data.

3. **Toggle views by cultivar type.** The user can switch between comparing against hybrid, inbred, or traditional baselines. The simulation itself doesn't change — only which RREC reference yield is shown.

### What it does NOT do
- The RREC data does not replace the ORYZA growth engine. ORYZA still runs the daily simulation and applies sound effects to specific biological processes.
- The RREC cultivars are US varieties (not IR72). The simulator's cultivar parameters remain IR72 from ORYZA2000. This means the simulator baseline (~7.3 t/ha) doesn't perfectly match any single RREC experiment, but falls between the inbred and traditional medians.
- No sound treatment was applied in the RREC experiments. The farm data is a no-sound baseline only.

### Day alignment
Day 0 of the simulation = emergence date in the RREC data (DOY 132, May 12 2022). When the user selects N treatment days, they see only those N days of growth trajectory, not the full season.

---

## Cross-Species Concordance Discount (replaces blanket 40%)

### Problem with the blanket discount
The original implementation applied a flat 0.4 multiplier to all non-rice confidence scores. This penalised two very different situations equally:
- Wheat at 300 Hz showing a positive effect (concordant with rice — same mechanism, similar frequency range)
- Oat at 300 Hz showing root inhibition (discordant — opposite response at the same frequency)

Both received the same 60% penalty, which is scientifically incorrect.

### New approach: trend-concordance scoring
The cross-species discount is now computed by checking how well each non-rice data point fits the Lorentzian curve established by rice-only data:

**Two-pass process:**
1. Fit the Lorentzian to rice-only CSV data points
2. For each cross-species data point, compute what the rice curve predicts at that frequency, then score concordance

**Concordance rules:**
- **Same direction, similar magnitude** (e.g., corn +10% at 200 Hz, rice predicts +18%): discount 0.7-0.95 (near-full trust)
- **Same direction, different magnitude** (e.g., wheat +15% at 300 Hz, rice predicts +39%): discount 0.5-0.7 (moderate trust)
- **Near-zero region** (e.g., wheat 1250 Hz no effect, rice predicts ~0%): discount 0.6-0.85 (the agreement on "nothing happens here" is informative)
- **Opposite direction** (e.g., oat -20% at 300 Hz, rice predicts +39%): discount 0.15 (species-specific divergence, heavily penalised)

**Biological rationale:** The mechanosensory machinery (MS ion channels, Ca²⁺ signalling) is conserved across monocot grasses. If a cross-species study shows the same response direction as rice at the same frequency, this reinforces the hypothesis that the conserved mechanism is driving the effect. If it shows the opposite, this suggests species-specific morphological differences are dominant, and the data should not be used as a rice proxy.

**For empirical windows** (biophysical boundaries like >90 dB inhibition): a moderate flat 0.6 discount is retained because these windows represent conserved membrane physics, not frequency-specific resonance.

---

## Honest Positioning

No published equation exists that predicts how a given sound frequency will affect rice crop physiology. The phytoacoustics literature reports results at specific frequencies but does not generalise them into a computable function. The PAFT resonance framework is a theory about why certain frequencies work, not a predictive equation. The mechanotransduction pathway is described qualitatively but lacks the quantitative coupling constants needed for first-principles modelling.

This simulator is a first attempt at bridging that gap. It uses a Lorentzian resonance curve — the physics-correct response of a damped harmonic oscillator — fitted to published rice data, with cross-species evidence scored by concordance with the rice trend. The result is a principled interpolation engine: more than a lookup table, less than a predictive model. Between studied frequencies, it gives physics-grounded estimates. Beyond the evidence range, it gives physics-grounded guesses that decay toward zero.

Whether those estimates are correct at unstudied frequencies is an empirical question that only new experiments can answer. The confidence score quantifies proximity to evidence, not prediction accuracy.
