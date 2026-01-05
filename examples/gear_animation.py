import py_gearworks as pgw
from build123d import *
from ocp_vscode import *
import numpy as np

n1 = 11
n2 = 33

cone_angle_1, cone_angle_2 = pgw.cone_angle_from_teeth(n1, n2)

gear1 = pgw.BevelGear(
    number_of_teeth=n1,
    cone_angle=cone_angle_1,
    helix_angle=np.pi / 6,
    height=5,
    profile_shift=0.3,
    backlash=0.025,
)
gear2 = pgw.BevelGear(
    number_of_teeth=n2,
    cone_angle=cone_angle_2,
    helix_angle=-np.pi / 6,
    height=5,
    profile_shift=-0.3,
    backlash=0.025,
)

a_gear1 = Compound(
    children=[gear1.build_part().solid() - gear1.face_location_top * Cylinder(2, 10)],
    label="gear1",
)
a_gear2 = Compound(
    children=[gear2.build_part().solid() - gear2.face_location_top * Cylinder(2, 10)],
    label="gear2",
)

gears = Compound(children=[a_gear1, a_gear2], label="gears")


duration = 2
n = duration * 30
time_track = np.linspace(0, duration, n + 1)
gear1_track = np.linspace(0, -gear1.pitch_angle * 180 / np.pi, n + 1) * duration
gear2_track = np.linspace(0, gear2.pitch_angle * 180 / np.pi, n + 1) * duration
# assign the tracks to the gears
animation = Animation(gears)
animation.add_track("/gears/gear1", "rz", time_track, gear1_track)
animation.add_track("/gears/gear2", "rz", time_track, gear2_track)

# Animations have a limitation that the movement can only be translation or rotation
# along a principal axis - see "rz" above.
# However, if a build123d location is applied to the compound, the animation will
# respect that location - which may include complex translation and re-orientation.
# That is why mesh_to() method needs to be called after the animation tracks are defined,
# and applied to the Compound after creation.
gear1.mesh_to(gear2, target_dir=pgw.LEFT, angle_bias=-1)
a_gear1.location = gear1.center_location_bottom
a_gear2.location = gear2.center_location_bottom


show(gears)

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
