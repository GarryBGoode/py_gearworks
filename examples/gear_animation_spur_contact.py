import py_gearworks as pgw
from build123d import *
from ocp_vscode import *
import numpy as np
from scipy.optimize import minimize

n1 = 7
n2 = 71

n_spline = 8

gamma1 = np.arctan2(n1, n2)
gamma2 = np.pi / 2 - gamma1

gear1 = pgw.SpurGear(
    number_of_teeth=n1,
    height=5,
    profile_shift=0.45,
    module=2.5,
    tip_fillet=0.2,
    root_fillet=0.1,
)
gear2 = pgw.SpurGear(
    number_of_teeth=n2,
    height=5,
    profile_shift=-0.2,
    module=2.5,
    tip_fillet=0.2,
    root_fillet=0.1,
)


a_gear1 = Compound(
    children=[
        gear1.build_part(n_points_hz=n_spline).solid()
        - gear1.face_location_top * Cylinder(2, 10)
    ],
    label="gear1",
)
a_gear2 = Compound(
    children=[
        gear2.build_part(n_points_hz=n_spline).solid()
        - gear2.face_location_top * Cylinder(2, 10)
    ],
    label="gear2",
)
gear2.center = gear2.pitch_radius * pgw.UP
gear1.mesh_to(gear2, target_dir=pgw.DOWN, backlash=0.1, angle_bias=-1)

line_of_contact = pgw.generate_line_of_contact(gear1, gear2, z_level=1)
edge_loc = pgw.line_to_b123d(line_of_contact[0])

gears = Compound(children=[a_gear1, a_gear2, edge_loc], label="gears")


a_gear1.location = gear1.center_location_bottom
a_gear2.location = gear2.center_location_bottom


duration = 4
n = duration * 300
time_track = np.linspace(0, duration, n + 1)
gear1_track = np.linspace(0, -gear1.pitch_angle * 180 / np.pi, n + 1)
gear2_track = np.linspace(0, gear2.pitch_angle * 180 / np.pi, n + 1)
# assign the tracks to the gears
animation = Animation(gears)
animation.add_track("/gears/gear1", "rz", time_track, gear1_track)
animation.add_track("/gears/gear2", "rz", time_track, gear2_track)


show(gears, deviation=0.02, angular_tolerance=0.055)

# Start animation
animation.animate(speed=1)
