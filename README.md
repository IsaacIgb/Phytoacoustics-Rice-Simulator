# Psychoacoustic Rice Growth Simulator

No published model or digital tool can yet demonstrate, from first principles, how a given sound frequency will influence crop physiology. This project is a research tool that gets as close as current evidence allows: a mechanistically-informed interpolation engine that couples sound treatments into an ORYZA-style rice growth model. It is **more than a lookup table, but less than a predictive model (due to scarcity of input data and not lack of will)**.

Built for researchers probing phytoacoustic treatments.

---

## What it does

For each day and sound regime, it:

1. Estimates plant-level SPL from the speaker layout using inverse-distance acoustics.
2. Computes the mechanotransduction latent `S_mech,air(f, SPL)` using:
   - an SPL activation curve,
   - a frequency-trained mechanosensory response model,
   - a high-frequency airborne decay term.
3. Uses this latent to drive pathway activations and form an outcome-specific mechanistic prior.
4. Fits Gaussian Processes to the residual effects with heteroscedastic noise based on study reliability and species.
5. Maps predicted effect sizes and confidence scores into daily modifiers on photosynthesis, LAI growth, water stress, stress penalties, and yield in the ORYZA-style crop engine.

Ultrasound seed soaks (`20–40 kHz`) and iWUE remain special-case pathways that bypass the main GP-mechanotransduction route.

## How to interpret it

The sound engine is a **literature-guided hypothesis generator**, not a dose-response oracle.

Confidence scores reflect **proximity to available evidence and model uncertainty**, not real-world accuracy at a given farm.

At the end of each simulation, you will see your sound exposure's ('treatment') results benchmarked against the RREC and non-treatment baselines for plot-level yield and phenology data. You will also see and using also see a coverage map of canopy-level SPL and simple coverage metrics. 

For full equations and assumptions, see:

- `docs/MATH_SPEC.md` — complete mathematical specification
- `docs/MODEL_NOTES.md` — architecture, theoretical framework, limitations

---

## Dual-phase UI

The UI separates sound treatment into:
- **Germination phase (days 1–10):** `GerminationSoundDriver` — stage bucket `"germination"`, typically higher SPL for seed coat penetration (Bochu 2003: 400 Hz / 106 dB).
- **Seedling-to-harvest (days 11+):** `VegetativeSoundDriver` — stage bucket `"seedling"` or `"vegetative"`, moderate SPL for stomatal and photosynthetic effects (Jusoh 2023: 350 Hz / 68 dB).

Both phases are configured independently in the API. The backend uses `DualPhaseSoundParams` and `run_simulation_dual_phase()` in `rice_core.py`.

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
