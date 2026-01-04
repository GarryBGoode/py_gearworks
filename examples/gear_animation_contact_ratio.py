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
    tip_fillet=0.0,
    root_fillet=0.0,
)
gear2 = pgw.SpurGear(
    number_of_teeth=n2,
    height=5,
    profile_shift=-0.2,
    module=2.5,
    tip_fillet=0.0,
    root_fillet=0.0,
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

gear1.mesh_to(gear2, target_dir=pgw.LEFT, backlash=0.2, angle_bias=-1)


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


show(gears, deviation=0.005, angular_tolerance=0.015)

# Start animation
animation.animate(speed=1)


# Used for generating screenshots for documentation via tcv_screenshots
def main():
    from tcv_screenshots import save_model, get_saved_models

    N_saver = 30
    angle_1 = gear1.pitch_angle * 180 / np.pi / N_saver
    angle_2 = gear2.pitch_angle * 180 / np.pi / N_saver

    for k in range(N_saver):
        b_gear1 = (
            gear1.center_location_bottom
            * Rotation((0, 0, angle_1 * k))
            * gear1.center_location_bottom.inverse()
            * a_gear1
        )
        b_gear2 = (
            gear2.center_location_bottom
            * Rotation((0, 0, -angle_2 * k))
            * gear2.center_location_bottom.inverse()
            * a_gear2
        )
        gears = Compound(children=[b_gear1, b_gear2], label="gears")

        save_model(
            gears,
            f"bevelgears_{k:02d}",
            {
                "cadWidth": 800,
                "height": 600,
                "position": (20, -20, 15),
                "target": (0, 0, 5),
                "zoom": 1.5,
            },
        )
    return get_saved_models()
