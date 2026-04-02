"""
api.py — HTTP backend for the rice sound simulator.

Uses real RREC 2022 weather data from Stuttgart AR when available.
Runs BOTH baseline and sound-treated simulations using the SoundDriver.
Returns RREC farm baseline data for comparison in the UI.
"""

import os, sys, json, math
from http.server import HTTPServer, SimpleHTTPRequestHandler
from urllib.parse import urlparse

sys.path.insert(0, os.path.dirname(__file__))

from rice_core import (
    RiceState, SimulationConfig, WeatherDay, SoundDriver,
    run_simulation, build_sound_driver,
)
from sound_response import estimate_all_effects
from acoustics import FieldConfig, compute_field_spl
from kpi_engine import compute_all_kpis

FRONTEND_DIR = os.path.join(os.path.dirname(__file__), "..", "frontend")
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")

# Load RREC baseline once at startup
_RREC = None
def get_rrec_baseline():
    global _RREC
    if _RREC is None:
        path = os.path.join(DATA_DIR, "rrec_baseline_2022.json")
        if os.path.exists(path):
            with open(path) as f:
                _RREC = json.load(f)
        else:
            _RREC = {}
    return _RREC


def build_rrec_weather(total_days):
    """Build WeatherDay list from RREC daily weather data."""
    rrec = get_rrec_baseline()
    weather_data = rrec.get("weather_daily", [])
    weather_list = []
    for i in range(total_days):
        if i < len(weather_data):
            w = weather_data[i]
            weather_list.append(WeatherDay(
                tmin=w["tmin_c"],
                tmax=w["tmax_c"],
                radiation=w["radiation_kj"],
            ))
        else:
            # Fallback to synthetic if RREC data runs out
            base_tmin = 21.0 + 2.0 * math.sin(2 * math.pi * i / 120)
            base_tmax = 32.0 + 2.0 * math.sin(2 * math.pi * i / 120)
            weather_list.append(WeatherDay(
                tmin=base_tmin, tmax=base_tmax,
                radiation=17000.0 + 3000.0 * math.sin(2 * math.pi * i / 120),
            ))
    return weather_list


def run_simulation_pipeline(params: dict) -> dict:
    """Run the full simulation pipeline and return results dict."""

    total_days = int(params.get("total_days", 120))
    sound_treatment_days = total_days  # Days slider = how long sound plays

    # Build weather for a full season (150 days) regardless of treatment duration
    weather = build_rrec_weather(150)

    # 1. Run BASELINE rice simulation (no sound) — always to maturity
    rice_config = SimulationConfig(
        total_days=150,  # Run full season
        water_regime=params.get("water_regime", "irrigated"),
        drought_severity=float(params.get("drought_severity", 0.5)),
        weather=weather,
    )
    baseline_history = run_simulation(rice_config)
    baseline_final = baseline_history[-1]

    # 2. Compute sound effects
    freq = float(params.get("frequency_hz", 1000))
    spl = float(params.get("spl_db_at_1m", 85))
    hours = float(params.get("hours_per_day", 3))
    stage = params.get("stage_bucket", "mixed")
    stress_ctx = "drought" if rice_config.water_regime == "drought" else "well_watered"

    sound_effects = estimate_all_effects(
        frequency_hz=freq, spl_db=spl, hours_per_day=hours,
        stage_bucket=stage, stress_context=stress_ctx, data_dir=DATA_DIR,
    )

    # 3. Build SoundDriver and run TREATED simulation
    # Sound is active for sound_treatment_days, then crop finishes to maturity
    sd = build_sound_driver(sound_effects)
    treated_history = run_simulation(rice_config, sound_driver=sd,
                                     sound_treatment_days=sound_treatment_days)
    treated_final = treated_history[-1]

    # 4. Derive canopy height from simulated stem weight (WST)
    # Rice plant height follows a power-law with stem dry weight:
    #   h = a × WST^b
    # Calibrated for lowland indica rice (IR72-type):
    #   a = 0.0387, b = 0.368
    # Fitted to: WST=5 → 0.07m (emergence), WST=3800 → 0.85m (heading),
    #            WST=6000 → 0.95m (maturity). Matches IR72 field observations.
    #
    # Uses TREATED simulation WST because the acoustics should reflect
    # the actual canopy the sound is hitting. Sound-treated crops are
    # taller (more photosynthesis → more stem weight → taller canopy).
    midpoint_day = min(sound_treatment_days // 2, len(treated_history) - 1)
    midpoint_wst = treated_history[midpoint_day].wst
    canopy_height = 0.0387 * max(1.0, midpoint_wst) ** 0.368
    canopy_height = round(max(0.05, min(canopy_height, 1.1)), 3)

    # 5. Compute field acoustics
    field_config = FieldConfig(
        field_length_m=float(params.get("field_length_m", 100)),
        field_width_m=float(params.get("field_width_m", 50)),
        speaker_count=int(params.get("speaker_count", 4)),
        spacing_m=float(params.get("spacing_m", 25)),
        speaker_height_m=float(params.get("speaker_height_m", 2.0)),
        canopy_height_m=canopy_height,
        spl_at_1m=spl, frequency_hz=freq, hours_per_day=hours,
    )
    grid_res = int(params.get("grid_resolution", 40))
    acoustics = compute_field_spl(field_config, grid_res)

    # 5. Sound effects dict
    sound_effects_dict = {}
    for k, v in sound_effects.items():
        sound_effects_dict[k] = {
            "effect_pct": v.effect_pct, "confidence": v.confidence,
            "effect_direction": v.effect_direction,
            "supporting_row_count": v.supporting_row_count, "notes": v.notes,
        }

    # 6. Acoustics dict
    acoustics_dict = {
        "mean_spl": acoustics.mean_spl, "min_spl": acoustics.min_spl,
        "max_spl": acoustics.max_spl, "median_spl": acoustics.median_spl,
        "coverage_above_60db": acoustics.coverage_above_60db,
        "coverage_above_70db": acoustics.coverage_above_70db,
        "uniformity": acoustics.uniformity,
    }

    # 7. KPIs
    kpis_raw = compute_all_kpis(
        baseline_yield=baseline_final.yield_proxy, sound_effects=sound_effects,
        water_stress=baseline_final.water_stress,
        water_regime=rice_config.water_regime,
        acoustics_result=acoustics_dict, hours_per_day=hours,
    )
    kpis_response = {}
    for k, v in kpis_raw.items():
        kpis_response[k] = {
            "name": v.name, "value": v.value, "confidence": v.confidence,
            "explanation": v.explanation, "components": v.components,
        }

    # 8. Time series for BOTH baseline and treated
    def ts_from(history):
        return {
            "days": [s.day for s in history],
            "dvs": [round(s.dvs, 3) for s in history],
            "biomass": [round(s.total_biomass, 1) for s in history],
            "lai": [round(s.lai, 3) for s in history],
            "wso": [round(s.wso, 1) for s in history],
        }

    # 9. Heatmap
    heatmap = {
        "grid_x": acoustics.grid_x.tolist(),
        "grid_y": acoustics.grid_y.tolist(),
        "spl_values": acoustics.spl_grid.tolist(),
        "speakers": [{"x": round(s.x, 1), "y": round(s.y, 1)} for s in acoustics.speakers],
    }

    # 10. RREC baseline data for comparison
    rrec = get_rrec_baseline()
    rrec_comparison = {}
    if rrec:
        rrec_comparison = {
            "source": rrec.get("source", ""),
            "citation": rrec.get("citation", ""),
            "year": rrec.get("year", 2022),
            "season": rrec.get("season", {}),
            "overall": rrec.get("overall", {}),
            "categories": rrec.get("categories", {}),
        }

    # 11. Yield deltas
    baseline_yield_tha = baseline_final.yield_proxy / 1000  # kg/ha → t/ha
    treated_yield_tha = treated_final.yield_proxy / 1000
    yield_change_pct = ((treated_yield_tha - baseline_yield_tha) / max(0.1, baseline_yield_tha)) * 100

    # 12. Confidence notes
    confidence_notes = []
    if all(v.confidence < 0.1 for v in sound_effects.values()):
        confidence_notes.append("Very low evidence support for this frequency/SPL combination.")
    if freq < 100 or freq > 5000:
        confidence_notes.append("Frequency outside primary evidence range (100-5000 Hz).")
    if spl > 100:
        confidence_notes.append("Very high SPL may cause tissue damage; effects uncertain.")
    if baseline_final.dvs < 1.5:
        confidence_notes.append(f"Simulation ended before grain filling (DVS={baseline_final.dvs:.2f}).")

    return {
        "inputs": params,
        "sound_treatment_days": sound_treatment_days,
        "baseline_final_day": baseline_final.day,
        "baseline_dvs": round(baseline_final.dvs, 3),
        "baseline_biomass": round(baseline_final.total_biomass, 0),
        "baseline_yield_tha": round(baseline_yield_tha, 2),
        "baseline_yield": round(baseline_final.yield_proxy, 0),
        "baseline_lai_max": round(max(s.lai for s in baseline_history), 2),
        "treated_yield_tha": round(treated_yield_tha, 2),
        "treated_biomass": round(treated_final.total_biomass, 0),
        "yield_change_pct": round(yield_change_pct, 1),
        "baseline_timeseries": ts_from(baseline_history),
        "treated_timeseries": ts_from(treated_history),
        "kpis": kpis_response,
        "sound_effects": sound_effects_dict,
        "acoustics_summary": acoustics_dict,
        "heatmap_data": heatmap,
        "rrec_baseline": rrec_comparison,
        "confidence_notes": confidence_notes,
    }


class SimHandler(SimpleHTTPRequestHandler):
    def do_GET(self):
        path = urlparse(self.path).path
        if path == "/" or path == "/index.html":
            fpath = os.path.join(FRONTEND_DIR, "index.html")
            if os.path.exists(fpath):
                self.send_response(200)
                self.send_header("Content-Type", "text/html; charset=utf-8")
                self.end_headers()
                with open(fpath, "rb") as f:
                    self.wfile.write(f.read())
                return
        if path == "/health":
            self._json_response(200, {"status": "ok"})
            return
        self.send_response(404)
        self.end_headers()

    def do_POST(self):
        if urlparse(self.path).path == "/simulate":
            content_len = int(self.headers.get("Content-Length", 0))
            body = self.rfile.read(content_len)
            try:
                params = json.loads(body) if body else {}
                result = run_simulation_pipeline(params)
                self._json_response(200, result)
            except Exception as e:
                import traceback; traceback.print_exc()
                self._json_response(500, {"error": str(e)})
            return
        self.send_response(404)
        self.end_headers()

    def do_OPTIONS(self):
        self.send_response(200)
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "POST, GET, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.end_headers()

    def _json_response(self, code, data):
        body = json.dumps(data).encode()
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, fmt, *args):
        pass


def main():
    port = int(os.environ.get("PORT", 8000))
    server = HTTPServer(("0.0.0.0", port), SimHandler)
    print(f"Rice Sound Simulator running at http://localhost:{port}")
    print("Press Ctrl+C to stop.")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down.")
        server.server_close()


if __name__ == "__main__":
    main()
