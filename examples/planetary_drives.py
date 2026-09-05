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

# This example is a bit of a playground for testing planetary gears.
# The goal is to create a set with any number of teeth on the sun and ring
# and insert appropriate number of planets.
# Also there's profile shifts.


n_ring = 115
n_sun = 11

shift_tooth = np.mod(n_ring - n_sun, 2) * 0.5 - 1
shift_sun = 0.15
shift_ring = -0.0
n_planet = int(np.floor((n_ring - n_sun) / 2 + shift_tooth))
shift_planet = (-shift_tooth - shift_sun + shift_ring) / 2
height = 5

backlash = 0.00

drill_hole_dia = 3


gearset = PlanetaryGearset(n_sun, n_ring, n_planet, 1)
gearset.num_planets = gearset.max_num_planets(addendum_ratio=1.1 + shift_planet)

gear_sun = SpurGear(
    number_of_teeth=n_sun, height=height, profile_shift=shift_sun, backlash=backlash
)
gear_planet = SpurGear(
    number_of_teeth=n_planet,
    height=height,
    profile_shift=shift_planet,
    backlash=backlash,
)
gear_ring = SpurRingGear(
    number_of_teeth=n_ring, profile_shift=shift_ring, height=height, backlash=backlash
)
gear_ring.angle = (gear_ring.pitch_angle / 2) * ((n_planet + 1) % 2)

placement_angles = gearset.get_planet_angles_distributed()

placement_directions = [
    scp_Rotation.from_euler("z", placement_angles[i]).apply(RIGHT)
    for i in range(gearset.num_planets)
]

a_gear_sun = Compound(
    children=[
        gear_sun.build_part().solid()
        - gear_sun.face_location_top * Cylinder(drill_hole_dia / 2, 10)
    ],
    label="gear_sun",
)
a_gear_ring = Compound(
    children=[
        gear_ring.build_part().solid()
        - gear_ring.face_location_top * Cylinder(drill_hole_dia / 2, 10)
    ],
    label="gear_ring",
)

a_planet_gears = []
for i in range(gearset.num_planets):
    dir = placement_directions[i]
    planet_gear = Compound(
        children=[
            gear_planet.build_part().solid()
            - gear_planet.face_location_top * Cylinder(drill_hole_dia / 2, 10)
        ],
        label=f"gear_planet_{i+1}",
    )
    gear_planet.mesh_to(gear_sun, target_dir=dir)
    planet_gear.location = gear_planet.center_location_bottom
    gear_planet.center = ORIGIN
    gear_planet.angle = 0
    a_planet_gears.append(planet_gear)


gear_drive = Compound(
    children=[a_gear_sun, a_gear_ring, *a_planet_gears],
    label="gear_drive",
)

duration = n_ring
n = duration * 30
time_track = np.linspace(0, duration, n + 1)
gear_sun_track = np.linspace(0, gear_sun.pitch_angle * 180 / np.pi, n + 1) * duration
gear_ring_track = np.linspace(0, -gear_ring.pitch_angle * 180 / np.pi, n + 1) * duration
gear_planet_track = (
    np.linspace(0, -gear_planet.pitch_angle * 180 / np.pi, n + 1) * duration
)
show(gear_drive, deviation=0.01)
# assign the tracks to the gears
animation = Animation(gear_drive)
animation.add_track("/gear_drive/gear_sun", "rz", time_track, gear_sun_track)
animation.add_track("/gear_drive/gear_ring", "rz", time_track, gear_ring_track)
for i in range(gearset.num_planets):
    animation.add_track(
        f"/gear_drive/gear_planet_{i+1}", "rz", time_track, gear_planet_track
    )
animation.add_track("/gear_drive", "rz", time_track, -gear_ring_track)


animation.animate(speed=1)
