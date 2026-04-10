# Psychoacoustic Rice Growth Simulator

No published model can yet predict, from first principles, how a given sound frequency will reshape rice crop physiology. This project is a research tool that gets as close as current evidence allows: a mechanistically-informed interpolation engine that couples sound treatments into an ORYZA-style rice growth model. It is **more than a lookup table, less than a predictive model**.

Built for researchers probing phytoacoustic treatments — not as a farm recommendation engine.

---

## What it does

The simulator combines four systems:

1. **Rice growth engine**
   Daily-step crop simulation simplified from ORYZA2000 (Bouman et al. 2001), driven by real RREC 2022 weather from Stuttgart, AR. Baseline IR72 yield is ~7.3 t/ha and matches USDA-ARS field observations.

2. **Sound-response engine (v4)**
   Estimates how pure-tone sound modifies rice physiology using:
   - A **mechanotransduction latent** S_mech,air(f, SPL) that encodes mechanosensory activation from frequency and local SPL.
   - A small set of **pathway activations** (seed-vigor, leaf gas exchange, drought resilience, membrane damage, ultrasound cavitation) that act as a mechanistic prior over outcomes.
   - Outcome-specific **Gaussian Processes on residuals** (Matérn 5/2, heteroscedastic noise) that correct this prior where data exist.

   It produces effect percentages and confidence scores for photosynthesis, vigor, water status sub-KPIs, stress resilience and yield. iWUE uses nearest-neighbour (standing-wave artefacts), and ultrasound (≥20 kHz) uses a separate cavitation lookup.

3. **Acoustics model**
   Free-field inverse-distance SPL distribution across a configurable speaker layout, with canopy-level SPL estimates and simple coverage metrics.

4. **RREC farm baseline**
   Plot-level yield and phenology for 152 plots and 21 cultivars at the USDA-ARS Dale Bumpers National Rice Research Center (Farag et al. 2024). Results are shown as deltas against typical hybrid, inbred and traditional yields.

---

## How to interpret it

- The sound engine is a **literature-guided hypothesis generator**, not a dose-response oracle.
- Mechanistic structure appears as:
  - SPL floor/ceiling and high-SPL caution band.
  - Mechanotransduction latent and pathway prior that encode which bands and stages should respond.
- The GPs then interpolate **residuals** around that prior, with uncertainty that grows rapidly away from data.

Confidence scores reflect **proximity to available evidence and model uncertainty**, not real-world accuracy at a given farm.

For full equations and assumptions, see:

- `docs/MATH_SPEC.md` — complete mathematical specification
- `docs/MODEL_NOTES.md` — architecture, theoretical framework, limitations

---

## Quick start

```bash
git clone https://github.com/your-user/rice-sound-sim.git
cd rice-sound-sim
```

Create and activate a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate        # macOS / Linux
# or
.\.venv\Scripts\activate         # Windows
```

Install dependencies:

```bash
pip install --upgrade pip
pip install -r requirements.txt  # numpy, scikit-learn, pyyaml, fastapi, uvicorn, pydantic
```

Run the backend API:

```bash
python backend/api.py
```

Then open `http://localhost:8000` in your browser.

---

## Sound algorithm (v4) in one paragraph

For each day and sound regime, the engine:

1. Computes local SPL at the plant from speaker layout and inverse-distance acoustics.
2. Computes S_mech,air(f, SPL) from:
   - an SPL activation curve (60–100 dB band, damage ceiling at 110 dB),
   - a MechActivation GP over frequency trained on stress-pathway KPIs,
   - a high-frequency airborne decay (5–15 kHz).
3. Uses S_mech,air to drive pathway activations and an outcome-specific mechanistic prior.
4. Fits outcome GPs to **residuals** (observed − prior) with heteroscedastic noise from study reliability and species.
5. Maps predicted effects into daily modifiers on photosynthesis, LAI growth, water stress and stress penalties in the ORYZA-style crop engine.

Ultrasound seed soaks (20–40 kHz) and iWUE remain special pathways that bypass the GP and the mechanotransduction latent.

---

## Dual-phase sound drivers

The UI separates sound treatment into:
- **Germination phase (days 1–10):** `GerminationSoundDriver` — stage bucket `"germination"`, typically higher SPL for seed coat penetration (Bochu 2003: 400 Hz / 106 dB).
- **Seedling-to-harvest (days 11+):** `VegetativeSoundDriver` — stage bucket `"seedling"` or `"vegetative"`, moderate SPL for stomatal and photosynthetic effects (Jusoh 2023: 350 Hz / 68 dB).

Both phases are configured independently in the API. The backend uses `DualPhaseSoundParams` and `run_simulation_dual_phase()` in `rice_core.py`.

---

## Validation

Run:

```bash
python tests/test_all.py               # 94 unit / integration tests
python tests/test_validation_cases.py  # 25 benchmark cases
```

Current results: **21 PASS, 4 PARTIAL, 0 FAIL / 25 total**

The validation harness checks rice anchor cases (Jusoh 350–380 Hz, Jeong 250–1500 Hz, Bochu 400/4000 Hz, Wang ultrasound) plus cross-species and "no effect" controls.

The 4 PARTIALs are all directionally correct:

| Case | Status | Note |
|------|--------|------|
| jusoh_standing_wave_353hz | PARTIAL | Standing-wave artefact cannot be captured by smooth GP |
| bochu_seed_inhibitory_4000hz | PARTIAL | Biophysical ceiling returns -15%, band expects [-13,-2] |
| ultrasound_seed_germination | PARTIAL | Rice GP overpredicts for barley at 43.5 kHz |
| barley_ultrasound_germination | PARTIAL | Same cross-species boundary |

See `tests/VALIDATION_REPORT.md` and `tests/VALIDATION_CASES.yaml` for full case definitions.

---

## Project structure

```
backend/
  sound_response.py     Mechanistic GP + pathway priors + CWSI water sub-KPIs (v4)
  rice_core.py          ORYZA-inspired rice growth engine + dual-phase sound drivers (v4)
  acoustics.py          Free-field canopy SPL model
  kpi_engine.py         Explicit-formula KPI computation
  api.py                HTTP backend — dual-phase pipeline
frontend/
  index.html            Single-page web UI
tests/
  test_all.py           94 unit/integration tests
  test_validation_cases.py  25 validation benchmarks
  VALIDATION_CASES.yaml     Benchmark case definitions
  VALIDATION_REPORT.md      Latest test results
data/
  IR72.yaml             Cultivar parameters
  rice_sound_master_v2.csv  Primary evidence dataset (76 rows, 46 usable)
  rrec_baseline_2022.json   RREC farm yields (152 plots, 21 cultivars)
  RREC_DAILY_2022.xlsx      Real weather data
docs/
  MATH_SPEC.md          Complete mathematical specification
  MODEL_NOTES.md        Architecture, theory, assumptions, limitations
```

---

## Scientific basis

**Rice biology:** Bouman et al. (2001) ORYZA2000; Li et al. (2017) ORYZA v3; Penning de Vries et al. (1989) maintenance coefficients; Monsi & Saeki (1953) canopy radiation.

**Sound evidence:** Jusoh & Ramlee 2023, Qi et al. 2010, Hou et al. 2009, Bochu et al. 2003, Jeong et al. 2008 & 2014, Hassan et al. 2014, Munasinghe et al. 2023, Wang 2022, Wang 2020, Kim et al. 2019, plus wheat/corn/oat/barley/sorghum/grass studies. Full source list in `MODEL_NOTES.md`.

**Architecture:** WOFOST/PCSE modular process design (Wageningen). GP regression: scikit-learn GaussianProcessRegressor, Matérn ν=5/2. v4 mechanotransduction latent and pathway activations.

---

For feedback or queries: [IsaacIgb](https://github.com/isaacigb)
