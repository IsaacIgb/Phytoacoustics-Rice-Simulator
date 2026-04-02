# Validation Report v3

Dataset: rice_sound_master_v2.csv

**19 PASS, 6 PARTIAL, 0 FAIL** / 25 total

| Case | Role | Status | Details |
|---|---|---|---|
| acoustics_distance_sanity_001 | acoustics_core | PASS | drop=6.1 dB |
| acoustics_two_source_sanity_001 | acoustics_core | PASS | combined=80.2 dB (+3.4) |
| inhibitory_grass_high_spl_001 | analog_inhibition | PARTIAL | effect=-3.9% at 95.0 dB |
| inhibitory_oat_root_300hz_001 | analog_inhibition | PASS | effect=+5.1%, conf=0.102 |
| no_effect_wheat_1250hz_001 | analog_negative_control | PASS | effect=+0.2%, conf=0.018 |
| no_effect_wheat_12000hz_001 | analog_negative_control | PASS | effect=+0.1%, conf=0.006 |
| barley_ultrasound_germination_43500hz_001 | analog_ultrasound_positive | PARTIAL | effect=+46.7% (expected [2,18]), direction OK |
| out_of_range_extreme_pure_tone_001 | core_boundary | PASS | effect=-15.0%, conf=0.450 |
| baseline_no_sound_001 | core_sanity | PASS | max|effect|=0.000% |
| jusoh_seedling_assimilation_350hz_001 | rice_primary | PASS | effect=+40.5%, conf=0.443 |
| jusoh_seedling_height_357hz_001 | rice_primary | PASS | effect=+15.1%, conf=0.314 |
| jusoh_iwue_peak_359hz_001 | rice_primary | PASS | effect=+62.7%, conf=1.000 |
| jusoh_standing_wave_353hz_001 | rice_primary | PARTIAL | effect=+39.2% (expected [-1,21]), direction OK |
| jeong_drought_rwc_800hz_001 | rice_primary | PASS | effect=+21.5%, conf=0.687 |
| jeong_drought_rwc_250hz_001 | rice_primary | PASS | effect=+10.5%, conf=0.687 |
| jeong_drought_gsw_1500hz_001 | rice_primary | PASS | effect=+87.8%, conf=0.369 |
| jeong_drought_fvfm_800hz_001 | rice_primary | PASS | effect=+7.3%, conf=0.263 |
| bochu_seedling_height_400hz_001 | rice_primary | PARTIAL | effect=+37.4% (expected [12,18]), direction OK |
| bochu_seed_inhibitory_4000hz_002 | rice_primary | PARTIAL | effect=-15.0% (expected [-13,-2]), direction OK |
| hou_yield_550hz_001 | rice_primary | PASS | effect=+6.9%, conf=0.583 |
| hassan_sheath_blight_550hz_001 | rice_primary | PASS | effect=-50.0%, conf=0.600 |
| ultrasound_seed_germination_20000hz_001 | rice_primary | PARTIAL | effect=+46.7% (expected [27,43]), direction OK |
| wang_seed_germination_time_40000hz_001 | rice_primary | PASS | effect=-46.7%, conf=0.650 |
| wang_seedling_dry_weight_35000hz_001 | rice_primary | PASS | effect=+12.0%, conf=1.000 |
| weak_evidence_munasinghe_350hz_001 | rice_secondary | PASS | effect=-24.9%, conf=0.733 |
