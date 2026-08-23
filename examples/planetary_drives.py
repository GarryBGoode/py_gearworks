from py_gearworks import *
from ocp_vscode import show, Animation
from copy import deepcopy
from build123d import *
from scipy.spatial.transform import Rotation as scp_Rotation
from scipy.optimize import minimize

n_sun = 8
n_ring = 83

n_planet = (n_ring - n_sun) // 2 - 0
num_planets = 3

height = 5
drill_hole_dia = 3

gear_sun = SpurGear(
    number_of_teeth=n_sun, height=height, profile_shift=0.5, backlash=0.025
)
gear_planet = SpurGear(
    number_of_teeth=n_planet, height=height, profile_shift=0.0, backlash=0.025
)
gear_ring = SpurRingGear(number_of_teeth=n_ring, height=height, backlash=0.025)
gear_ring.angle = (gear_ring.pitch_angle / 2) * ((n_planet + 1) % 2)

placement_unit_angle = 2 * PI / n_ring / (1 + n_sun / n_ring)

placement_angles = [
    np.round((i * 2 * PI / num_planets) / placement_unit_angle) * placement_unit_angle
    for i in range(num_planets)
]

placement_directions = [
    scp_Rotation.from_euler("z", placement_angles[i]).apply([1, 0, 0])
    for i in range(num_planets)
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
for i in range(num_planets):
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
# assign the tracks to the gears
animation = Animation(gear_drive)
animation.add_track("/gear_drive/gear_sun", "rz", time_track, gear_sun_track)
animation.add_track("/gear_drive/gear_ring", "rz", time_track, gear_ring_track)
for i in range(num_planets):
    animation.add_track(
        f"/gear_drive/gear_planet_{i+1}", "rz", time_track, gear_planet_track
    )
# animation.add_track("/gear_drive", "rz", time_track, -gear_ring_track)

show(gear_drive)
animation.animate(speed=1)
