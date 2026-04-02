# Rice Sound Simulator

No published equation exists that can predict how a given sound frequency will affect rice crop physiology. This simulator is a first attempt: a principled interpolation engine grounded in resonance physics and calibrated against published evidence — more than a lookup table, less than a predictive model. Whether it's right at unstudied frequencies is an empirical question that only new experiments can answer.

---

Research-grade MVP simulating the relative effect of sound treatments on poacea crops' (namely rice) growth, water status, stress physiology, and field acoustic exposure. Built for researchers exploring phytoacoustic treatment design, not as a farm recommendation engine.

## What It Does

Combines four independent systems:

1. **Rice growth engine** — daily-step crop simulation simplified from ORYZA2000 (Bouman et al. 2001), driven by real RREC 2022 weather data from Stuttgart AR. Harvest index 0.47, baseline yield ~7.3 t/ha with real weather — validated against RREC farm observations.

2. **Sound-response engine** — estimates how sound treatment modifies rice physiology using Lorentzian resonance fitting with CWSI-based water sub-KPIs, calibrated against `rice_sound_master_v2.csv` (76 rows, 47 usable from 11 rice studies + cross-species cereal/grass data). Water status uses three independent Lorentzian tracks per CWSI theory (Idso 1981): gsw_drought, gsw_wellwatered, and RWC. Ultrasound (≥20 kHz) uses a separate cavitation pathway. Produces effect percentages with confidence scores for 6 biological outcomes.

3. **Acoustics model** — free-field inverse-distance SPL distribution across a configurable speaker layout, producing coverage heatmaps and uniformity metrics.

4. **RREC farm baseline** — real plot-level yield and phenology data from the USDA-ARS Dale Bumpers National Rice Research Center (Farag et al. 2024, 152 plots, 21 cultivars). Results are displayed as comparisons against real farm yields, with toggle views for hybrid (13.7 t/ha), inbred (7.95 t/ha), and traditional (6.3 t/ha) cultivar categories.

## What It Is NOT

This simulator does not predict absolute outcomes from first principles. The sound-response engine is a **literature-guided hypothesis generator**: it can say "published studies suggest +35% photosynthesis at 350 Hz for rice seedlings" but it cannot predict what happens at an untested frequency from membrane biophysics alone. The Lorentzian resonance model provides physics-grounded interpolation between studied frequencies, but the resonance parameters (f₀, A, γ) are empirically fitted, not derived from cellular mechanics.

Confidence scores reflect **proximity to available evidence**, not prediction accuracy.

For full transparency on every equation, simplification, and design choice, see:
- `docs/MATH_SPEC.md` — complete mathematical specification (48 numbered equations with rationale)
- `docs/MODEL_NOTES.md` — algorithm architecture, theoretical framework, assumptions, honest limitations

## Setup

### 1. Clone or download the project

```bash
git clone https://github.com/your-user/rice-sound-sim.git
cd rice-sound-sim
```

---

### 2. Install Python 3

**macOS**

```bash
brew install python
```

**Windows**

1. Download Python 3 from https://www.python.org/downloads/
2. During install, tick **“Add Python to PATH”**.

You should then have `python3` and `pip3` available in a terminal (on Windows you may need to open a new terminal after install).

In the unlikely event you're running an older version of python (pre 2020), bash your commands using python and pip, not python3 or pip3.

---

### 3. Create and activate a virtual environment

**macOS (Terminal / zsh)**

```bash
python3 -m venv .venv
source .venv/bin/activate
```

**Windows (Command Prompt or PowerShell)**

```bat
python3 -m venv .venv
.\.venv\Scripts\activate
```

Your prompt should now start with `(.venv)`.

---

### 4. Install Python dependencies (inside the venv)

```bash
pip3 install --upgrade pip
pip3 install numpy pyyaml
```

---

### 5. Run the backend API

From the project root, with the virtual environment active:

```bash
cd /path/to/dowloaded/rice-sound-sim
source .venv/bin/activate
python3 backend/api.py    
```
Then go to you default browser and paste http://localhost:8000 into the searchbar to open application

For Each new terminal session:

```bash
cd /path/to/downloaded/rice-sound-sim
source .venv/bin/activate      

# or on Windows:
# .\.venv\Scripts\activate
python3 backend/api.py
```

## How the Sound Algorithm Works

The sound-response engine uses a 4-layer architecture plus two special pathways (detailed in `MODEL_NOTES.md`):

**Layer 1 — Biophysical boundaries.** Hard physics constraints: 60 dB activation floor, 110 dB damage ceiling. The 90-110 dB caution zone applies only when extrapolating beyond studied SPLs — when the CSV has data measured at the query SPL, the data overrides the caution. See `docs/MATH_SPEC.md` Eqs SE1-SE3.

**Layer 2 — Empirical windows.** Nine rice-specific and two cross-species frequency×SPL×stage bands providing directional guidance. Blended with the Lorentzian at 15% weight when the Lorentzian has ≥3 data points (85/15 CSV/window blend).

**Layer 3 — Lorentzian resonance.** Per-outcome curve `R(f) = A / (1 + ((f - f₀) / γ)²)` fitted by grid search to v2 CSV data. The grid covers the full data frequency range with adaptive amplitude and bandwidth limits. See `docs/MATH_SPEC.md` Eqs SE4-SE5.

**Layer 4 — Cross-species concordance.** Non-rice data scored by agreement with the rice Lorentzian trend. Concordant data gets 0.7-0.96 confidence; discordant gets 0.15. Replaces the old blanket 40% discount.

**Water sub-KPIs (CWSI theory).** Water status uses three independent Lorentzians per crop water stress index theory (Idso et al. 1981): `gsw_drought` (stomatal conductance under deficit, Jeong 2014), `gsw_wellwatered` (Jusoh 2023), and `rwc` (relative water content, Jeong 2014). iWUE uses nearest-neighbour lookup (data oscillates too wildly for Lorentzian). Aggregate uses CWSI-weighted combination: drought = 60% gsw + 25% RWC + 15% iWUE.

**Ultrasound pathway (≥ 20 kHz).** Bypasses the Lorentzian entirely — cavitation mechanism, not resonance. Direct distance-weighted lookup from ultrasound data. No dosage scaling. Outcome-specific routing (germination vs vigor tracks).

**Theoretical basis:** The model unifies three phytoacoustic theories — PAFT resonance, mechanotransduction/Ca²⁺ pathway, and ultrasound cavitation. See `docs/MODEL_NOTES.md` for the full theoretical framework.

## User Interface

Single-page web interface with a **days-first** input flow:

**Treatment duration** — total days (1-120) determines the UI state:
- ≤10 days: germination-only mode (one sound panel)
- &gt;10 days: dual-phase mode — "Germination Sound (Days 1-10)" and "Seedling→Harvest Sound (Days 11+)" side by side, each with independent frequency, SPL, hours/day, and audio preview

This separation exists because the science separates these domains: germination studies use ultrasound/high-SPL seed treatments (Bochu 400 Hz/106 dB, barley 43.5 kHz) while post-germination studies use audible-range airborne sound at moderate SPL (Jusoh 350 Hz/68 dB).

**Frequency slider** — three-region linear mapping optimised for sub-1 kHz resolution (where most evidence is):
- First 50% (0-500): 50 Hz to 1 kHz, ~2 Hz steps. Fine control over the evidence-dense audible range.
- Next 30% (500-800): 1 kHz to 10 kHz, ~30 Hz steps. Coverage of cross-species mid-range data.
- Last 20% (800-1000): 10 kHz to 60 kHz, ~250 Hz steps. Purple ⚡ indicator above 20 kHz marks ultrasound.

See `docs/MATH_SPEC.md` Section 8 for the slider equations.

**Field layout and water regime** are uniform inputs (the physical field doesn't change between growth phases). Paddy field visualiser shows a live preview with red pulsing speakers on a translucent green grid.

## Rice Growth Engine

Simplified from ORYZA2000 following WOFOST/PCSE architecture patterns: modular process classes (Phenology, Assimilation, Respiration, Partitioning, LeafDynamics, Spikelets, WaterStress), rate-then-integrate separation, cultivar parameters loadable from YAML.

Key processes: hourly heat unit accumulation with bilinear temperature response, big-leaf canopy photosynthesis (Monsi-Saeki analytical integral), DVS-dependent dry matter partitioning, dual-phase LAI growth (exponential below LAI=1.0, SLA-based above).

Sound enters the simulation as a **SoundDriver** — a formal driver that modifies biological rates (photosynthesis, vigor, water status, stress resilience) during each daily step, not post-hoc. Effect percentages are scaled by confidence before becoming multipliers.

See `docs/MATH_SPEC.md` Sections 1-4 and 6 for all rice engine equations. See `docs/MODEL_NOTES.md` for documented simplifications vs ORYZA2000.

## Validation

```bash
python tests/test_all.py               # 49 unit/integration tests
python tests/test_validation_cases.py  # 25 validation benchmarks (v3)
```

Results (v3): **19 PASS, 6 PARTIAL, 0 FAIL / 25 total**

| Domain role | Pass/Total | Key cases |
|---|---|---|
| rice_primary | 11/15 | Jusoh assimilation/height/iWUE, Jeong RWC/gsw/FvFm, Hou yield, Hassan sheath blight, Wang germination/dry weight |
| rice_secondary | 1/1 | Munasinghe weak evidence |
| analog_negative_control | 2/2 | Wheat 1250 Hz, wheat 12 kHz (no false positives) |
| analog_inhibition | 1/2 | Oat concordance discount (grass high-SPL is PARTIAL) |
| analog_ultrasound_positive | 0/1 | Barley 43.5 kHz (PARTIAL — rice data overpredicts for barley) |
| core_sanity | 1/1 | No sound = no effect |
| core_boundary | 1/1 | Extreme SPL damage ceiling |
| acoustics_core | 2/2 | Inverse-square 6 dB drop, source combination +3 dB |

The 6 PARTIALs are all directionally correct — magnitude outside tight bands due to Lorentzian smoothing, composite KPI averaging, or cross-species transfer boundaries.

See `tests/VALIDATION_REPORT.md` and `tests/VALIDATION_CASES.yaml` for full case definitions and results.

## Project Structure

```
backend/
  rice_core.py          ORYZA-inspired rice growth engine (PCSE architecture)
  sound_response.py     Lorentzian resonance + CWSI water sub-KPIs + ultrasound pathway
  acoustics.py          Free-field canopy SPL model
  kpi_engine.py         Explicit-formula KPI computation
  api.py                HTTP backend (stdlib, no frameworks)
frontend/
  index.html            Single-page web UI with export/reset
tests/
  test_all.py           49 unit/integration tests
  test_validation_cases.py  25 validation benchmarks (v3)
  VALIDATION_CASES.yaml     Benchmark case definitions
  VALIDATION_REPORT.md      Latest test results
data/
  IR72.yaml             Cultivar parameters (PCSE CropDataProvider pattern)
  rice_sound_master_v2.csv  Primary evidence dataset (76 rows, 47 usable)
  rrec_baseline_2022.json   RREC farm yields (152 plots, 21 cultivars)
  RREC_DAILY_2022.xlsx      Real weather data driving the simulation
  reference_run_IR72_irrigated.json  Locked baseline for regression testing
docs/
  MATH_SPEC.md          Complete mathematical specification
  MODEL_NOTES.md        Architecture, theory, assumptions, limitations
```

## Documentation

| File | Contents |
|---|---|
| `docs/MATH_SPEC.md` | Every equation with rationale. Sections: Phenology, Photosynthesis, Respiration, LAI, Sound-Response (Lorentzian + CWSI water sub-KPIs + ultrasound + nearest-neighbour), Acoustics, UI |
| `docs/MODEL_NOTES.md` | Algorithm architecture, CWSI theory (Idso 1981), water sub-KPI rationale (Nature 2016 meta-analysis), data sources, honest limitations |
| `tests/VALIDATION_CASES.yaml` | 25 benchmark cases with domain roles and expected ranges (v3) |
| `tests/VALIDATION_REPORT.md` | Latest results: 19 PASS, 6 PARTIAL, 0 FAIL |

## Scientific Basis

**Rice biology:** Bouman et al. (2001) ORYZA2000; Li et al. (2017) ORYZA v3; Penning de Vries et al. (1989) maintenance coefficients; Monsi & Saeki (1953) canopy radiation.

**Sound evidence:** Jusoh & Ramlee 2023, Qi et al. 2010, Hou et al. 2009, Bochu et al. 2003, Jeong et al. 2008, Hassan et al. 2014, Munasinghe et al. 2023, Wang 2022, Wang 2020, plus wheat/corn/oat/barley/sorghum/grass studies. Full source list in `MODEL_NOTES.md`.

**Architecture:** WOFOST/PCSE modular process design (Wageningen).

---

For feedback or queries: [IsaacIgb](https://github.com/isaacigb)
