# Copyright 2024 Gergely Bencsik
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#     http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from build123d import *
from py_gearworks import *
from ocp_vscode import show, set_port

# set_port(3939)


def write_svg(part, example_number, viewpos=(-100, -100, 70)):
    """Save an image of the BuildPart object as SVG"""

    # builder: BuildPart = BuildPart._get_context()

    visible, hidden = part.project_to_viewport(viewpos)
    max_dimension = max(*Compound(children=visible + hidden).bounding_box().size)
    exporter = ExportSVG(scale=100 / max_dimension)
    exporter.add_layer("Visible")
    # exporter.add_layer("Hidden", line_color=(99, 99, 99), line_type=LineType.ISO_DOT)
    exporter.add_shape(visible, layer="Visible")
    # exporter.add_shape(hidden, layer="Hidden")
    exporter.write(f"assets/general_ex{example_number}.svg")


def example_1():
    # [Ex. 1]
    gear1 = SpurGear(number_of_teeth=12, module=1.2, height=5)
    gear2 = SpurGear(number_of_teeth=24, module=1.2, height=5)
    gear1.mesh_to(gear2, target_dir=RIGHT)
    gear_part_1 = gear1.build_part()
    gear_part_2 = gear2.build_part()
    # [Ex. 1]
    write_svg(
        Compound([gear_part_1, gear_part_2]), example_number=1, viewpos=(0, -40, 70)
    )


def example_2():
    # [Ex. 2]
    gear1 = SpurGear(number_of_teeth=12)
    gear2 = SpurRingGear(number_of_teeth=24)
    gear1.mesh_to(gear2, target_dir=RIGHT + UP)
    gear_part_1 = gear1.build_part()
    gear_part_2 = gear2.build_part()
    # [Ex. 2]
    write_svg(
        Compound([gear_part_1, gear_part_2]), example_number=2, viewpos=(0, -40, 70)
    )


def example_3():
    # [Ex. 3]
    gear1 = SpurGear(number_of_teeth=8)
    gear2 = SpurGear(number_of_teeth=8, profile_shift=0.7, tip_truncation=0)
    gear3 = SpurGear(number_of_teeth=8, profile_shift=0.7, tip_truncation=0.2)
    gear1.mesh_to(gear2, target_dir=LEFT)
    gear3.mesh_to(gear2, target_dir=RIGHT)
    gear_part_1 = gear1.build_part()
    gear_part_2 = gear2.build_part()
    gear_part_3 = gear3.build_part()

    # [Ex. 3]
    write_svg(
        Compound([gear_part_1, gear_part_2, gear_part_3]),
        example_number=3,
        viewpos=(0, -40, 70),
    )


def example_4():
    # [Ex. 4]
    gear1 = HelicalGear(number_of_teeth=12, height=5, helix_angle=PI / 6)
    gear2 = HelicalGear(number_of_teeth=24, height=5, helix_angle=-PI / 6)
    gear1.mesh_to(gear2, target_dir=RIGHT)
    gear_part_1 = gear1.build_part()
    gear_part_2 = gear2.build_part()
    # [Ex. 4]
    write_svg(
        Compound([gear_part_1, gear_part_2]), example_number=4, viewpos=(0, -40, 20)
    )


def example_45():
    # [Ex. 45]
    gear1 = HelicalGear(number_of_teeth=12, height=15, helix_angle=PI / 4, z_anchor=0.5)
    gear2 = HelicalGear(number_of_teeth=24, height=15, helix_angle=PI / 4, z_anchor=0.5)
    gear1.mesh_to(gear2, target_dir=RIGHT)
    gear_part_1 = gear1.build_part()
    gear_part_2 = gear2.build_part()
    # [Ex. 45]
    write_svg(
        Compound([gear_part_1, gear_part_2]), example_number=45, viewpos=(40, -40, 100)
    )


def example_5():
    # [Ex. 5]
    gear1 = SpurGear(number_of_teeth=24, height=10, crowning=200)
    gear_part_1 = gear1.build_part()
    # [Ex. 5]
    write_svg(gear_part_1, example_number=5, viewpos=(0, -40, 20))


def example_6():
    # [Ex. 6]
    n1 = 13
    n2 = 32
    cone_angle1, cone_angle2 = cone_angle_from_teeth(n1, n2, axis_angle=PI / 2)
    gear1 = BevelGear(number_of_teeth=n1, cone_angle=cone_angle1, height=5)
    gear2 = BevelGear(number_of_teeth=n2, cone_angle=cone_angle2, height=5)
    gear1.mesh_to(gear2, target_dir=RIGHT)
    gear_part_1 = gear1.build_part()
    gear_part_2 = gear2.build_part()
    # [Ex. 6]
    write_svg(
        Compound([gear_part_1, gear_part_2]), example_number=6, viewpos=(-50, -100, 30)
    )


def example_7():
    # [Ex. 7]
    n1 = 12
    n2 = 31
    cone_angle1, cone_angle2 = cone_angle_from_teeth(n1, n2)
    gear1 = BevelGear(
        number_of_teeth=n1, cone_angle=cone_angle1, height=5, helix_angle=0.5
    )
    gear2 = BevelGear(
        number_of_teeth=n2, cone_angle=cone_angle2, height=5, helix_angle=-0.5
    )
    gear1.mesh_to(gear2, target_dir=RIGHT + UP)
    gear_part_1 = gear1.build_part()
    gear_part_2 = gear2.build_part()
    # [Ex. 7]
    write_svg(
        Compound([gear_part_1, gear_part_2]), example_number=7, viewpos=(0, -100, 70)
    )


def example_8():
    # [Ex. 8]
    n1 = 12
    n2 = 31
    gear1 = CycloidGear(number_of_teeth=n1, height=5, inside_cycloid_coefficient=0.25)
    gear2 = CycloidGear(number_of_teeth=n2, height=5, inside_cycloid_coefficient=0.5)
    gear1.adapt_cycloid_radii(gear2)
    gear1.mesh_to(gear2, target_dir=RIGHT)
    gear_part_1 = gear1.build_part()
    gear_part_2 = gear2.build_part()
    # [Ex. 8]
    write_svg(
        Compound([gear_part_1, gear_part_2]), example_number=8, viewpos=(0, -20, 70)
    )


def example_9():
    # [Ex. 9]
    n_gear = 12
    n_rack = 40
    gear = SpurGear(number_of_teeth=n_gear, module=2, height=5, profile_shift=0.25)
    rack = InvoluteRack(number_of_teeth=n_rack, module=2, height=5)
    rack.mesh_to(gear, target_dir=LEFT, offset=19)
    gear_part = gear.build_part()
    rack_part = rack.build_part()
    # [Ex. 9]
    write_svg(Compound([gear_part, rack_part]), example_number=9, viewpos=(0, -20, 70))


def example_10():
    # [Ex. 10]
    n_gear = 12
    n_rack = 40
    gear = HelicalGear(
        number_of_teeth=n_gear,
        module=2,
        height=25,
        helix_angle=PI / 3,
        herringbone=True,
    )
    rack = HelicalRack(
        number_of_teeth=n_rack,
        module=2,
        height=25,
        helix_angle=PI / 3,
        herringbone=True,
    )
    rack.mesh_to(gear, target_dir=LEFT, offset=0)
    gear_part = gear.build_part()
    rack_part = rack.build_part()
    # [Ex. 10]
    write_svg(Compound([gear_part, rack_part]), example_number=10, viewpos=(0, -45, 70))


def example_11():
    # [Ex. 11]
    gear1 = SpurGear(12, module=2, profile_shift=0.7, backlash=0.1)
    gear2 = SpurGear(13, module=2, profile_shift=0.6, backlash=0.1)
    gear3 = SpurGear(
        14,
        module=2,
        profile_shift=0.5,
    )
    gear1.mesh_to(gear2, target_dir=RIGHT, angle_bias=1)
    gear3.mesh_to(gear2, target_dir=LEFT, backlash=0, angle_bias=0)
    gear_part_1 = gear1.build_part()
    gear_part_2 = gear2.build_part()
    gear_part_3 = gear3.build_part()
    # [Ex. 11]
    write_svg(
        Compound([gear_part_1, gear_part_2, gear_part_3]),
        example_number=11,
        viewpos=(0, -30, 100),
    )


if __name__ == "__main__":
    example_1()
    example_2()
    example_3()
    example_4()
    example_45()
    example_5()
    example_6()
    example_7()
    example_8()
    example_9()
    example_10()
    example_11()
