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

# This is the example shown in README.md.
# Its content between the SNIPPET START/END markers is kept in sync with the
# README automatically - see scripts/sync_readme_examples.py.

# --8<-- [start:snippet]
from py_gearworks import *
from ocp_vscode import show
from build123d import *

# create 2 spur gears
gear1 = SpurGear(
    number_of_teeth=12,
    module=2,
    height=4,
    profile_shift=0.3,
)
gear2 = SpurGear(
    number_of_teeth=23,
    module=2,
    height=4,
)

# move and align gear 1 next to gear 2 in the Y direction
# backlash can be optionally specified
# angle_bias conrtols location within backlash range (-1 to 1)
# backlash is a coefficient of module
# there will be 0.2 mm distance between inactive tooth sides in this example
gear1.mesh_to(gear2, target_dir=UP, backlash=0.1, angle_bias=1)

# generate build123d Part objects
gear_part_1 = gear1.build_part()
gear_part_2 = gear2.build_part()

# center-bores are recommended to be added separately via build123d workflow
# center_location_top can be used as a build123d location object
# location * Hole means placement of Hole at that location (build123d syntax)
hole_obj_1 = gear1.center_location_top * Hole(radius=2, depth=4)
gear_part_1 = gear_part_1.cut(hole_obj_1)
hole_obj_2 = gear2.center_location_top * Hole(radius=2, depth=4)
gear_part_2 = gear_part_2.cut(hole_obj_2)

# export to STEP files (build123d export function)
# note: export retains positioning, gear1 will not be at origin
export_step(gear_part_1, "gear1.step")
export_step(gear_part_2, "gear2.step")

# export to DXF files (build123d export function)
gear_wire_1 = gear1.build_boundary_wire()
gear_wire_2 = gear2.build_boundary_wire()

exporter = ExportDXF(unit=Unit.MM, line_weight=0.5)
exporter.add_layer("Layer 1", line_type=LineType.CONTINUOUS)
exporter.add_shape(gear_wire_1, layer="Layer 1")
exporter.add_shape(gear_wire_2, layer="Layer 1")
exporter.write("gears.dxf")

# visualize parts
show(gear_part_1, gear_part_2)
# --8<-- [end:snippet]
