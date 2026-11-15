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

"""Representative gear set for the M0 safety net (REWORK_PLAN section 3).

Each case is a named factory that builds a fresh gear through the user-level API
only, so the set survives internal refactors. Factories are used instead of gear
instances because gears are mutable (``mesh_to()``, ``build_part()``) and tests must
not leak state into each other.

Geometry snapshots, fitting metrics and curated matrices all draw from this list.
Add cases here, don't define representative gears inline in individual tests.

Conventions (as in the wrapper): ``cone_angle`` is the full cone angle, i.e. twice
the pitch cone angle. Module is 1 unless the case says otherwise.
"""

import dataclasses
from typing import Any, Callable

import numpy as np
import pytest

import py_gearworks as pgw
from py_gearworks.defs import PI, RIGHT, UP


@dataclasses.dataclass(frozen=True)
class GearCase:
    """A named representative gear.

    Attributes
    ----------
    name : str
        Stable identifier, used as pytest id and as snapshot key. Don't rename
        without regenerating snapshots.
    factory : Callable
        Returns a new gear object on each call.
    purpose : str
        What the case is meant to stress.
    known_issue : str
        Non-empty if the case currently produces wrong output (see TODO.md). Tests
        should xfail on it rather than snapshot a broken baseline.
    """

    name: str
    factory: Callable[[], Any]
    purpose: str
    known_issue: str = ""

    def build(self):
        """Return a fresh gear object."""
        return self.factory()


# Rack tooth becomes sharp at h = (pi/4) / tan(alpha) = 2.158 module for alpha=20deg.
_SHARP_RACK_HEIGHT = (PI / 4) / np.tan(20 * PI / 180)
_SHARP_GEAR_PROFILE_SHIFT = 0.5

# Bevel pinion / gear pair, 90deg shaft angle.
_BEVEL_PAIR_CONE_ANGLES = pgw.cone_angle_from_teeth(12, 60)

# Roots-blower style lobes: pitch circle divided into 2*z arcs, each rolled out by one
# full turn of the rolling circle -> r_roll = r_pitch / (2*z). Addendum and dedendum
# are the rolling circle diameter: 2 * r_roll = m/2.
_ROOTS_Z = 3
_ROOTS_ROLL_COEFF = 1 / (2 * _ROOTS_Z)
_ROOTS_HEIGHT_COEFF = 0.5


GEAR_CASES = [
    GearCase(
        name="spur_z6_sharp",
        factory=lambda: pgw.SpurGear(
            number_of_teeth=6,
            addendum_coefficient=_SHARP_RACK_HEIGHT,
            dedendum_coefficient=_SHARP_RACK_HEIGHT,
            profile_shift=_SHARP_GEAR_PROFILE_SHIFT,
            tip_truncation=0,
            enable_undercut=True,
            height=2,
        ),
        purpose=(
            "Extremely small pinion, sharp-rack tooth height, sharp tips. Deep "
            "undercut, involute and trochoid fitting, disappearing addendum and "
            "dedendum arcs."
        ),
    ),
    GearCase(
        name="spur_z10_shift_truncated",
        factory=lambda: pgw.SpurGear(
            number_of_teeth=10,
            profile_shift=0.8,
            tip_truncation=0.1,
            height=2,
        ),
        purpose=(
            "Small pinion with large profile shift, the tip would be pointed "
            "(pointed around x=0.7), so tip truncation is active."
        ),
    ),
    GearCase(
        name="spur_z35",
        factory=lambda: pgw.SpurGear(number_of_teeth=35, height=4),
        purpose=(
            "Vanilla spur gear. Root is just below the base circle, which tests "
            "finding the undercut-involute crossover when the undercut is tiny or "
            "absent."
        ),
    ),
    GearCase(
        name="helical_z3_worm_approx",
        factory=lambda: pgw.HelicalGear(
            number_of_teeth=3,
            helix_angle=PI / 2 * 0.85,
            height=20,
            z_anchor=0.5,
        ),
        purpose=(
            "Large-twist helical (worm-approximation screw gear from "
            "examples.worm_approx). Stresses vertical surface fitting."
        ),
    ),
    GearCase(
        name="ring_z120_fillets",
        factory=lambda: pgw.SpurRingGear(
            number_of_teeth=120,
            root_fillet=0.2,
            tip_fillet=0.1,
            height=4,
        ),
        purpose="Large ring gear with root and tip fillets.",
    ),
    GearCase(
        name="bevel_z16_undercut",
        factory=lambda: pgw.BevelGear(
            number_of_teeth=16,
            cone_angle=PI / 2,
            enable_undercut=True,
            root_fillet=0,
            height=2,
        ),
        purpose="45deg bevel with undercut.",
    ),
    GearCase(
        name="bevel_z16_sharp_root",
        factory=lambda: pgw.BevelGear(
            number_of_teeth=16,
            cone_angle=PI / 2,
            enable_undercut=False,
            root_fillet=0,
            height=2,
        ),
        purpose="45deg bevel without undercut or root fillet: sharp root corner.",
    ),
    GearCase(
        name="bevel_z16_fillets",
        factory=lambda: pgw.BevelGear(
            number_of_teeth=16,
            cone_angle=PI / 2,
            enable_undercut=False,
            root_fillet=0.2,
            tip_fillet=0.1,
            height=2,
        ),
        purpose="45deg bevel with root and tip fillets.",
    ),
    GearCase(
        name="bevel_z60_flat",
        factory=lambda: pgw.BevelGear(
            number_of_teeth=60,
            cone_angle=_BEVEL_PAIR_CONE_ANGLES[1],
            height=4,
        ),
        purpose="Large, close-to-flat bevel gear. Pair of bevel_z12_pinion.",
    ),
    GearCase(
        name="bevel_z12_pinion",
        factory=lambda: pgw.BevelGear(
            number_of_teeth=12,
            cone_angle=_BEVEL_PAIR_CONE_ANGLES[0],
            height=4,
        ),
        purpose="Small, steep bevel pinion. Pair of bevel_z60_flat.",
    ),
    GearCase(
        name="cycloid_z18",
        factory=lambda: pgw.CycloidGear(number_of_teeth=18, height=2),
        purpose="Regular cycloid gear (rolling circles at half pitch radius).",
    ),
    GearCase(
        name="cycloid_z3_roots",
        factory=lambda: pgw.CycloidGear(
            number_of_teeth=_ROOTS_Z,
            inside_cycloid_coefficient=_ROOTS_ROLL_COEFF,
            outside_cycloid_coefficient=_ROOTS_ROLL_COEFF,
            addendum_coefficient=_ROOTS_HEIGHT_COEFF,
            dedendum_coefficient=_ROOTS_HEIGHT_COEFF,
            tip_truncation=0,
            root_fillet=0,
            tip_fillet=0,
            height=4,
        ),
        purpose=(
            "Roots-blower style 3-lobe cycloid. Full epi- and hypocycloid arcs, "
            "addendum and dedendum arcs vanish."
        ),
    ),
]

GEAR_CASES_BY_NAME = {case.name: case for case in GEAR_CASES}

# Paired cases, for mesh-level checks (contact ratio, intersection).
GEAR_CASE_PAIRS = [
    ("bevel_z12_pinion", "bevel_z60_flat"),
    ("spur_z10_shift_truncated", "spur_z35"),
]


def case_params(cases=GEAR_CASES):
    """Return the cases as ``pytest.param`` objects with the case name as id.

    Usage: ``@pytest.mark.parametrize("case", case_params())``.
    """
    return [pytest.param(case, id=case.name) for case in cases]


def case_params_xfail_known(cases=GEAR_CASES):
    """Like ``case_params()``, but cases with a ``known_issue`` are marked xfail.

    Use for checks the known issue breaks (e.g. solid validity).
    """
    return [
        pytest.param(
            case,
            id=case.name,
            marks=(
                [pytest.mark.xfail(reason=case.known_issue, strict=True)]
                if case.known_issue
                else []
            ),
        )
        for case in cases
    ]


def build_pair(pair):
    """Build the two gears of a ``GEAR_CASE_PAIRS`` entry, the first meshed to the
    second.

    Parameters
    ----------
    pair : tuple[str, str]
        Case names.

    Returns
    -------
    tuple
        The two gear objects.
    """
    gear1 = GEAR_CASES_BY_NAME[pair[0]].build()
    gear2 = GEAR_CASES_BY_NAME[pair[1]].build()
    gear1.mesh_to(gear2, target_dir=RIGHT)
    return gear1, gear2


def pair_id(pair):
    """Pytest id of a ``GEAR_CASE_PAIRS`` entry."""
    return f"{pair[0]}--{pair[1]}"


if __name__ == "__main__":
    from ocp_vscode import show

    GAP = 5

    def _view_row(items, y):
        """Place (label, parts) groups next to each other along x, each group's
        bounding box starting where the previous one ended, centered on y."""
        shown = []
        x_cursor = 0.0
        for label, parts in items:
            bbox = parts[0].bounding_box()
            for p in parts[1:]:
                bbox = bbox.add(p.bounding_box())
            offset = (x_cursor - bbox.min.X) * RIGHT + (y - bbox.center().Y) * UP
            for k, p in enumerate(parts):
                p = p.translate(offset)
                p.label = f"{label}_{k}" if len(parts) > 1 else label
                shown.append(p)
            x_cursor += bbox.size.X + GAP
        return shown

    singles = []
    for case in GEAR_CASES:
        gear = case.build()
        singles.append((case.name, [gear.build_part(), gear.build_boundary_wire(1)]))
    row_singles = _view_row(singles, y=0)

    pairs = []
    for pair in GEAR_CASE_PAIRS:
        gear1, gear2 = build_pair(pair)
        pairs.append(
            (
                pair_id(pair),
                [gear1.build_part(), gear2.build_part()],
            )
        )
    singles_half_height = max(
        p.bounding_box().size.Y / 2 for _, parts in singles for p in parts
    )
    pairs_half_height = max(
        p.bounding_box().size.Y / 2 for _, parts in pairs for p in parts
    )
    row_pairs = _view_row(
        pairs, y=-(singles_half_height + pairs_half_height + 2 * GAP)
    )

    show(*row_singles, *row_pairs)
