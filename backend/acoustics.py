"""
acoustics.py — Simple canopy exposure model.

Uses free-field inverse-distance SPL approximation:
  Lp(r) = Lp_1m - 20*log10(r)

Combines multiple sources logarithmically.
Ignores reflections, wind, and air absorption in v1.

See ASSUMPTIONS_USED.md for documented simplifications.
"""

import math
import numpy as np
from dataclasses import dataclass
from typing import List, Tuple, Optional


@dataclass
class Speaker:
    """A single speaker in the field."""
    x: float          # Position x (m)
    y: float          # Position y (m)
    height: float     # Speaker height above ground (m)
    spl_at_1m: float  # SPL at 1 m reference distance (dB)


@dataclass
class FieldConfig:
    """Field and speaker layout configuration."""
    field_length_m: float = 100.0
    field_width_m: float = 50.0
    speaker_count: int = 4
    spacing_m: float = 25.0
    speaker_height_m: float = 2.0
    canopy_height_m: float = 0.8
    spl_at_1m: float = 85.0
    frequency_hz: float = 1000.0
    hours_per_day: float = 3.0


@dataclass
class AcousticsResult:
    """Results from the acoustics computation."""
    spl_grid: np.ndarray       # 2D grid of SPL values (dB)
    grid_x: np.ndarray         # X coordinates of grid points
    grid_y: np.ndarray         # Y coordinates of grid points
    mean_spl: float            # Area-weighted mean SPL across field
    min_spl: float             # Minimum SPL in field
    max_spl: float             # Maximum SPL in field
    median_spl: float          # Median SPL
    coverage_above_60db: float # Fraction of field above 60 dB
    coverage_above_70db: float # Fraction of field above 70 dB
    uniformity: float          # 1 - (std/mean), higher = more uniform
    speakers: List[Speaker]    # Speaker positions used


def spl_at_distance(spl_1m: float, distance_m: float) -> float:
    """
    Free-field SPL at distance r from a point source.
    
    Lp(r) = Lp(1m) - 20*log10(r)
    
    This is the inverse-square law for sound pressure in free field.
    Valid assumption for outdoor open-field conditions.
    
    Documented simplification: ignores atmospheric absorption,
    ground reflections, foliage scattering, and speaker directivity.
    """
    if distance_m <= 0.01:
        distance_m = 0.01  # Avoid log(0); cap at 1 cm
    return spl_1m - 20.0 * math.log10(distance_m)


def combine_spl(spl_values: List[float]) -> float:
    """
    Combine multiple incoherent sound sources logarithmically.
    
    Lp_total = 10 * log10( sum(10^(Lp_i/10)) )
    
    This assumes incoherent (uncorrelated) sources, which is
    appropriate for speakers playing the same pure tone only if
    they are not phase-locked. For v1 this is acceptable.
    """
    if not spl_values:
        return 0.0
    total = sum(10.0 ** (spl / 10.0) for spl in spl_values)
    if total <= 0:
        return 0.0
    return 10.0 * math.log10(total)


def generate_speaker_positions(config: FieldConfig) -> List[Speaker]:
    """
    Generate speaker positions in a grid pattern within the field.
    
    Speakers are placed in a regular grid with the specified spacing,
    centered within the field boundaries.
    """
    speakers = []
    
    if config.speaker_count <= 0:
        return speakers
    
    if config.speaker_count == 1:
        # Single speaker at field center
        speakers.append(Speaker(
            x=config.field_length_m / 2,
            y=config.field_width_m / 2,
            height=config.speaker_height_m,
            spl_at_1m=config.spl_at_1m,
        ))
        return speakers
    
    # Calculate grid dimensions
    aspect = config.field_length_m / max(config.field_width_m, 0.1)
    n_along = max(1, int(math.sqrt(config.speaker_count * aspect)))
    n_across = max(1, int(config.speaker_count / n_along))
    
    # Adjust to not exceed speaker_count
    while n_along * n_across > config.speaker_count:
        if n_along > n_across:
            n_along -= 1
        else:
            n_across -= 1
    
    # Add remaining speakers along the longer axis
    remaining = config.speaker_count - (n_along * n_across)
    
    # Spacing within field
    margin_x = config.spacing_m / 2
    margin_y = config.spacing_m / 2
    
    usable_length = config.field_length_m - 2 * margin_x
    usable_width = config.field_width_m - 2 * margin_y
    
    dx = usable_length / max(1, n_along - 1) if n_along > 1 else 0
    dy = usable_width / max(1, n_across - 1) if n_across > 1 else 0
    
    for i in range(n_along):
        for j in range(n_across):
            x = margin_x + i * dx if n_along > 1 else config.field_length_m / 2
            y = margin_y + j * dy if n_across > 1 else config.field_width_m / 2
            speakers.append(Speaker(
                x=x, y=y,
                height=config.speaker_height_m,
                spl_at_1m=config.spl_at_1m,
            ))
    
    # Place remaining speakers along center line
    for k in range(remaining):
        x = margin_x + (k + 1) * usable_length / (remaining + 1)
        speakers.append(Speaker(
            x=x, y=config.field_width_m / 2,
            height=config.speaker_height_m,
            spl_at_1m=config.spl_at_1m,
        ))
    
    return speakers


def compute_field_spl(config: FieldConfig,
                      grid_resolution: int = 50) -> AcousticsResult:
    """
    Compute SPL distribution across the field canopy.
    
    For each grid point at canopy height, calculates the 3D distance
    to each speaker and computes the combined SPL from all sources.
    
    Grid resolution determines the number of points along each axis.
    """
    speakers = generate_speaker_positions(config)
    
    if not speakers:
        empty = np.zeros((grid_resolution, grid_resolution))
        gx = np.linspace(0, config.field_length_m, grid_resolution)
        gy = np.linspace(0, config.field_width_m, grid_resolution)
        return AcousticsResult(
            spl_grid=empty, grid_x=gx, grid_y=gy,
            mean_spl=0, min_spl=0, max_spl=0, median_spl=0,
            coverage_above_60db=0, coverage_above_70db=0,
            uniformity=0, speakers=speakers,
        )
    
    # Create grid
    gx = np.linspace(0, config.field_length_m, grid_resolution)
    gy = np.linspace(0, config.field_width_m, grid_resolution)
    spl_grid = np.zeros((grid_resolution, grid_resolution))
    
    canopy_z = config.canopy_height_m
    
    for ix, x in enumerate(gx):
        for iy, y in enumerate(gy):
            spl_contributions = []
            for spk in speakers:
                # 3D distance from speaker to canopy point
                dx = x - spk.x
                dy = y - spk.y
                dz = spk.height - canopy_z
                distance = math.sqrt(dx**2 + dy**2 + dz**2)
                
                spl_at_point = spl_at_distance(spk.spl_at_1m, distance)
                spl_contributions.append(spl_at_point)
            
            spl_grid[ix, iy] = combine_spl(spl_contributions)
    
    # Summary statistics
    flat = spl_grid.flatten()
    mean_spl = float(np.mean(flat))
    min_spl = float(np.min(flat))
    max_spl = float(np.max(flat))
    median_spl = float(np.median(flat))
    
    coverage_60 = float(np.mean(flat >= 60.0))
    coverage_70 = float(np.mean(flat >= 70.0))
    
    std_spl = float(np.std(flat))
    uniformity = 1.0 - (std_spl / max(abs(mean_spl), 1.0))
    uniformity = max(0.0, min(1.0, uniformity))
    
    return AcousticsResult(
        spl_grid=spl_grid,
        grid_x=gx,
        grid_y=gy,
        mean_spl=round(mean_spl, 1),
        min_spl=round(min_spl, 1),
        max_spl=round(max_spl, 1),
        median_spl=round(median_spl, 1),
        coverage_above_60db=round(coverage_60, 3),
        coverage_above_70db=round(coverage_70, 3),
        uniformity=round(uniformity, 3),
        speakers=speakers,
    )


if __name__ == "__main__":
    config = FieldConfig(
        field_length_m=100, field_width_m=50,
        speaker_count=4, spacing_m=25,
        speaker_height_m=2.0, canopy_height_m=0.8,
        spl_at_1m=85.0, frequency_hz=1000, hours_per_day=3.0,
    )
    result = compute_field_spl(config, grid_resolution=20)
    print(f"Mean SPL: {result.mean_spl} dB")
    print(f"Min: {result.min_spl}, Max: {result.max_spl}")
    print(f"Coverage >60dB: {result.coverage_above_60db*100:.0f}%")
    print(f"Coverage >70dB: {result.coverage_above_70db*100:.0f}%")
    print(f"Uniformity: {result.uniformity:.2f}")
    print(f"Speakers: {len(result.speakers)}")
