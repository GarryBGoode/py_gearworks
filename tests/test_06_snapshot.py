# Copyright 2026 Gergely Bencsik
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#     http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Geometry snapshot tests over the representative gear set (REWORK_PLAN M0).

The current code is the oracle: snapshots record its output, and later code must
reproduce it. Per gear case, three things are recorded:

1. Exact geometry: the analytic 2D tooth profile (before spline fitting) at a few z
   levels, as sample points plus curve length. Compared by distance from the stored
   points to the current curve, so a change of curve parametrization doesn't break
   it, only a change of shape does.
2. Solid: validity, volume, area, face and edge count, bounding box.
3. Fit deviation: distance from exact flank points to the built solid's faces. This
   is regression-checked against the snapshot (no worse than recorded, with margin).
   Detailed fitting metrics are separate (REWORK_PLAN 5.2).

Regenerate with ``pytest tests/test_06_snapshot.py --update-snapshots``. Snapshots
store the commit they were made from; regenerate only when a geometry change is
intended, and say why in the commit.

Only the "exact geometry adapter" section below touches package internals. When the
internals change (M2, M3), update that section, not the snapshots.
"""

import json
import platform
import subprocess
from pathlib import Path

import numpy as np
import pytest
from build123d import Compound, Vector

from gear_cases import case_params_xfail_known

SNAPSHOT_DIR = Path(__file__).parent / "snapshots"
SNAPSHOT_VERSION = 1

N_Z_LEVELS = 3
N_PROFILE_POINTS = 101  # stored points per profile
N_PROFILE_DENSE = 20001  # current-curve samples to measure distance against
N_FLANK_POINTS = 15  # flank points per z level for the solid deviation check

# Tolerances, for module 1 gears.
PROFILE_ATOL = 1e-5  # exact curves; well above polyline sagitta (~1e-7)
LENGTH_RTOL = 1e-6
SOLID_RTOL = 1e-4  # volume, area: approximations, platforms may differ slightly
BBOX_ATOL = 1e-3
FLANK_DEV_FACTOR = 1.5  # fit deviation may not grow more than this...
FLANK_DEV_ATOL = 1e-6  # ...plus this


# --- exact geometry adapter -------------------------------------------------------


def _z_levels(gear):
    z_vals = gear.gearcore.z_vals
    return np.linspace(z_vals[0], z_vals[-1], N_Z_LEVELS)


def _world(gear, profile, points):
    return gear.gearcore.transform(profile.transform(points))


def exact_profile(gear, z):
    """Full single-tooth profile (dedendum arc, flank, addendum arc, mirrored
    flank) at height z, as a callable of t in [0, 1] returning world points, and
    its length."""
    profile = gear.gearcore.curve_gen_at_z(z)
    chain = profile.profile

    def curve(t):
        return _world(gear, profile, chain(t))

    return curve, float(chain.length)


def exact_flank(gear, z):
    """One tooth flank at height z, as a callable of t in [0, 1] returning world
    points."""
    profile = gear.gearcore.curve_gen_at_z(z)

    def curve(t):
        return _world(gear, profile, profile.tooth_curve(t))

    return curve


# --- measurements -----------------------------------------------------------------


def point_to_polyline_distance(points, polyline):
    """Distance of each point to a polyline.

    Parameters
    ----------
    points : np.ndarray
        Shape (m, 3).
    polyline : np.ndarray
        Shape (n, 3), n >= 2.

    Returns
    -------
    np.ndarray
        Shape (m,).
    """
    a = polyline[:-1]
    ab = polyline[1:] - a
    ab_sq = np.maximum(np.sum(ab * ab, axis=1), 1e-300)
    dist = np.empty(len(points))
    for i, p in enumerate(points):
        u = np.clip(np.sum((p - a) * ab, axis=1) / ab_sq, 0.0, 1.0)
        closest = a + u[:, None] * ab
        dist[i] = np.sqrt(np.min(np.sum((p - closest) ** 2, axis=1)))
    return dist


def measure_profiles(gear):
    profiles = []
    t = np.linspace(0, 1, N_PROFILE_POINTS)
    for z in _z_levels(gear):
        curve, length = exact_profile(gear, z)
        profiles.append(
            {"z": float(z), "length": length, "points": curve(t).tolist()}
        )
    return profiles


def measure_solid(part):
    bbox = part.bounding_box()
    return {
        "is_valid": bool(part.is_valid),
        "volume": float(part.volume),
        "area": float(part.area),
        "n_faces": len(part.faces()),
        "n_edges": len(part.edges()),
        "bbox_min": list(bbox.min.to_tuple()),
        "bbox_max": list(bbox.max.to_tuple()),
    }


def measure_flank_deviation(gear, part):
    shell = Compound(part.faces())
    t = np.linspace(0.02, 0.98, N_FLANK_POINTS)
    dist = []
    for z in _z_levels(gear):
        points = exact_flank(gear, z)(t)
        dist.extend(shell.distance_to(Vector(*p)) for p in points)
    return {"max": float(np.max(dist)), "mean": float(np.mean(dist))}


def measure(case):
    gear = case.build()
    part = gear.build_part()
    return {
        "profiles": measure_profiles(gear),
        "solid": measure_solid(part),
        "flank_deviation": measure_flank_deviation(gear, part),
    }


# --- snapshot I/O -----------------------------------------------------------------


def _git(*args):
    try:
        out = subprocess.run(
            ["git", *args],
            cwd=Path(__file__).parent,
            capture_output=True,
            text=True,
            check=True,
        )
        return out.stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def _meta():
    status = _git("status", "--porcelain", "--untracked-files=no")
    return {
        "snapshot_version": SNAPSHOT_VERSION,
        "commit": _git("rev-parse", "HEAD"),
        "dirty": None if status is None else bool(status),
        "python": platform.python_version(),
        "platform": platform.platform(),
    }


def snapshot_path(case):
    return SNAPSHOT_DIR / f"{case.name}.json"


def write_snapshot(case, data):
    SNAPSHOT_DIR.mkdir(exist_ok=True)
    content = {"meta": _meta(), "case": case.name, "purpose": case.purpose, **data}
    snapshot_path(case).write_text(json.dumps(content, indent=1) + "\n")


def read_snapshot(case):
    path = snapshot_path(case)
    if not path.exists():
        pytest.fail(
            f"No snapshot for {case.name}. Generate it with "
            "`pytest tests/test_06_snapshot.py --update-snapshots`."
        )
    return json.loads(path.read_text())


# --- comparison -------------------------------------------------------------------


def compare_profiles(gear, ref_profiles):
    errors = []
    t_dense = np.linspace(0, 1, N_PROFILE_DENSE)
    for ref in ref_profiles:
        curve, length = exact_profile(gear, ref["z"])
        dev = point_to_polyline_distance(np.array(ref["points"]), curve(t_dense))
        if dev.max() > PROFILE_ATOL:
            errors.append(
                f"profile z={ref['z']:.4g}: max distance from snapshot "
                f"{dev.max():.3e} > {PROFILE_ATOL:.0e} (at point {int(dev.argmax())})"
            )
        if length != pytest.approx(ref["length"], rel=LENGTH_RTOL):
            errors.append(
                f"profile z={ref['z']:.4g}: length {length:.9g} != snapshot "
                f"{ref['length']:.9g}"
            )
    return errors


def compare_solid(solid, ref):
    errors = []
    if not solid["is_valid"]:
        errors.append("solid is not valid")
    for key in ("volume", "area"):
        if solid[key] != pytest.approx(ref[key], rel=SOLID_RTOL):
            errors.append(f"{key} {solid[key]:.9g} != snapshot {ref[key]:.9g}")
    for key in ("n_faces", "n_edges"):
        if solid[key] != ref[key]:
            errors.append(f"{key} {solid[key]} != snapshot {ref[key]}")
    for key in ("bbox_min", "bbox_max"):
        if solid[key] != pytest.approx(ref[key], abs=BBOX_ATOL):
            errors.append(f"{key} {solid[key]} != snapshot {ref[key]}")
    return errors


def compare_flank_deviation(dev, ref):
    limit = ref["max"] * FLANK_DEV_FACTOR + FLANK_DEV_ATOL
    if dev["max"] > limit:
        return [
            f"flank deviation max {dev['max']:.3e} > limit {limit:.3e} "
            f"(snapshot {ref['max']:.3e})"
        ]
    return []


@pytest.mark.parametrize("case", case_params_xfail_known())
def test_snapshot(case, update_snapshots):
    if update_snapshots:
        write_snapshot(case, measure(case))
        pytest.skip(f"snapshot written: {snapshot_path(case).name}")

    ref = read_snapshot(case)
    gear = case.build()
    part = gear.build_part()

    errors = compare_profiles(gear, ref["profiles"])
    errors += compare_solid(measure_solid(part), ref["solid"])
    errors += compare_flank_deviation(
        measure_flank_deviation(gear, part), ref["flank_deviation"]
    )
    assert not errors, f"{case.name} differs from snapshot:\n  " + "\n  ".join(errors)
