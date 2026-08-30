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


from py_gearworks import *
from ocp_vscode import show, Animation
from copy import deepcopy
from build123d import *
from scipy.spatial.transform import Rotation as scp_Rotation
from scipy.optimize import minimize

ring_dia = 90


n_sun_1 = 33
n_sun_2 = 36

n_planet = 8
n_ring = n_sun_1 + 2 * n_planet + 3
n_ring = n_ring - n_ring % 3 + 3  # make n_ring multiple of 3

module = ring_dia / n_ring
num_planets = 3

height = 10
drill_hole_dia = 3
hollow_center_dia = (n_sun_1 - 5) * module

print(f"n_ring: {n_ring}")
print(f"n_sun_1: {n_sun_1}")
print(f"n_sun_2: {n_sun_2}")
print(f"n_planet: {n_planet}")
print(f"module: {module}")

ratio = (n_sun_2 - n_sun_1) / n_sun_2 / (1 + n_sun_1 / n_ring)
print(f"ratio: {ratio}")
print(f"inverse ratio: {1/ratio}")

gear_planet = SpurGear(
    number_of_teeth=n_planet,
    module=module,
    height=height,
    profile_shift=0.25,
    backlash=0.0,
    z_anchor=0,
)
gear_ring = SpurRingGear(
    number_of_teeth=n_ring,
    module=module,
    profile_shift=-0.0,
    height=height,
    backlash=0.0,
    z_anchor=0,
)
gear_ring.angle = (gear_ring.pitch_angle / 2) * ((n_planet + 1) % 2)


gear_sun_1 = SpurGear(
    number_of_teeth=n_sun_1,
    module=module,
    height=height / 2,
    profile_shift=1.2 + 0.05,
    addendum_coefficient=0.8,
    backlash=0.00,
    z_anchor=0,
)
gear_sun_2 = SpurGear(
    number_of_teeth=n_sun_2,
    module=module,
    height=height / 2,
    profile_shift=-0.35 - 0.1,
    backlash=0.00,
    z_anchor=0,
)

sun_ref_gear = SpurGear(
    number_of_teeth=n_sun_1,
    module=module,
    height=height,
    profile_shift=1.2 + 0.05,
    backlash=0.00,
    z_anchor=0.0,
)
gear_sun_2.center = OUT * (height / 2)
gear_sun_2.angle = gear_sun_1.gearcore.shape_recipe(gear_sun_1.p2z(1)).transform.angle


def distance_cost_func(x):
    shift_1 = x[0]
    shift_2 = x[1]
    shift_planet = x[2]

    gear_sun_1.inputparam.profile_shift = shift_1
    gear_sun_1.calc_params()
    gear_sun_2.inputparam.profile_shift = shift_2
    gear_sun_2.calc_params()
    gear_planet.inputparam.profile_shift = shift_planet
    gear_planet.calc_params()

    l1 = calc_involute_mesh_distance(
        gear_sun_1.r_base,
        gear_planet.r_base,
        -gear_sun_1.angle_base,
        -gear_planet.angle_base,
        gear_planet.pitch_angle,
        inside_ring=False,
        backlash=0.0,
    )
    l2 = calc_involute_mesh_distance(
        gear_sun_2.r_base,
        gear_planet.r_base,
        -gear_sun_2.angle_base,
        -gear_planet.angle_base,
        gear_planet.pitch_angle,
        inside_ring=False,
        backlash=0.0,
    )
    l3 = calc_involute_mesh_distance(
        gear_planet.r_base,
        gear_ring.r_base,
        -gear_planet.angle_base,
        -gear_ring.angle_base,
        gear_ring.pitch_angle,
        inside_ring=True,
        backlash=0.0,
    )

    cost_val = (l1 - l2) ** 2 + (l1 - l3) ** 2 + (l2 - l3) ** 2
    return cost_val


sol1 = minimize(
    distance_cost_func,
    x0=[
        gear_sun_1.inputparam.profile_shift,
        gear_sun_2.inputparam.profile_shift,
        gear_planet.inputparam.profile_shift,
    ],
    bounds=[(-1, 1.5), (-1, 1.5), (-0.2, 1.5)],
)

print(f"Optimization result: {sol1}")

gear_sun_2.center = OUT * (height / 2)

placement_unit_angle = 2 * PI / n_ring / (1 + n_sun_1 / n_ring)

placement_angles = [
    np.round((i * 2 * PI / num_planets) / placement_unit_angle) * placement_unit_angle
    for i in range(num_planets)
]

placement_directions = [
    scp_Rotation.from_euler("z", placement_angles[i]).apply([1, 0, 0])
    for i in range(num_planets)
]

a_gear_sun_1 = Compound(
    children=[
        gear_sun_1.build_part().solid()
        - gear_sun_1.face_location_top * Cylinder(hollow_center_dia / 2, height)
    ],
    label="gear_sun_1",
)

a_gear_sun_2 = Compound(
    children=[
        gear_sun_2.build_part().solid()
        - gear_sun_2.face_location_top * Cylinder(hollow_center_dia / 2, height)
    ],
    label="gear_sun_2",
)

a_gear_ring = Compound(
    children=[
        gear_ring.build_part().solid()
        - gear_ring.face_location_top * Cylinder(drill_hole_dia / 2, height)
    ],
    label="gear_ring",
)

a_planet_gears = []
for i in range(num_planets):
    dir = placement_directions[i]
    planet_gear = Compound(
        children=[
            gear_planet.build_part().solid()
            - gear_planet.face_location_top
            * Cylinder(drill_hole_dia / 2, height=height * 2)
        ],
        label=f"gear_planet_{i+1}",
    )
    gear_planet.mesh_to(gear_ring, target_dir=dir)
    planet_gear.location = gear_planet.center_location_bottom
    gear_planet.center = ORIGIN
    gear_planet.angle = 0
    a_planet_gears.append(planet_gear)


gear_drive = Compound(
    children=[a_gear_sun_1, a_gear_sun_2, a_gear_ring, *a_planet_gears],
    label="gear_drive",
)

duration = n_sun_1 // 3 * 2
n = duration * 30
time_track = np.linspace(0, duration, n + 1)
gear_sun_track_1 = (
    np.linspace(0, gear_sun_1.pitch_angle * 180 / np.pi, n + 1) * duration
)
gear_sun_track_2 = (
    np.linspace(0, gear_sun_2.pitch_angle * 180 / np.pi, n + 1) * duration
)
gear_ring_track = np.linspace(0, -gear_ring.pitch_angle * 180 / np.pi, n + 1) * duration
gear_planet_track = (
    np.linspace(0, -gear_planet.pitch_angle * 180 / np.pi, n + 1) * duration
)

show(gear_drive)
# assign the tracks to the gears
animation = Animation(gear_drive)
animation.add_track("/gear_drive/gear_sun_1", "rz", time_track, gear_sun_track_1)
animation.add_track("/gear_drive/gear_sun_2", "rz", time_track, gear_sun_track_2)
animation.add_track("/gear_drive/gear_ring", "rz", time_track, gear_ring_track)
for i in range(num_planets):
    animation.add_track(
        f"/gear_drive/gear_planet_{i+1}", "rz", time_track, gear_planet_track
    )
animation.add_track(
    "/gear_drive", "rz", time_track, -gear_sun_track_1 * n_sun_2 / n_sun_2
)


animation.animate(speed=1)
