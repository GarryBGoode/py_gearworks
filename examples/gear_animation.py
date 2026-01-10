import py_gearworks as pgw
from build123d import *
from ocp_vscode import *
import numpy as np

n1 = 13
n2 = 53

cone_angle_1, cone_angle_2 = pgw.cone_angle_from_teeth(n1, n2)

gear1 = pgw.BevelGear(
    number_of_teeth=n1,
    cone_angle=cone_angle_1,
    helix_angle=np.pi / 6,
    height=5,
    profile_shift=0.21,
    backlash=0.05,
)
gear2 = pgw.BevelGear(
    number_of_teeth=n2,
    cone_angle=cone_angle_2,
    helix_angle=-np.pi / 6,
    height=5,
    profile_shift=-0.21,
    backlash=0.05,
)


a_gear1 = Compound(
    children=[gear1.build_part().solid() - gear1.face_location_top * Cylinder(2, 10)],
    label="gear1",
)
a_gear2 = Compound(
    children=[gear2.build_part().solid() - gear2.face_location_top * Cylinder(2, 10)],
    label="gear2",
)

# adding 1% above top level for rendering reasons
# otherwise it partially collides with gear surfaces and renders incorrectly
z_value = 0


def get_minmax_point(gear):
    invo_curve, trf = gear.involute_curve_at_z(z_value)
    invo_min_point = trf(invo_curve(0))
    invo_max_point = trf(invo_curve(1))
    return invo_min_point, invo_max_point


gear2.center = gear2.pitch_radius * pgw.RIGHT
gear1.mesh_to(gear2, target_dir=pgw.LEFT, angle_bias=1)
a_gear1.location = gear1.center_location_bottom
a_gear2.location = gear2.center_location_bottom

min_point_1, max_point_1 = get_minmax_point(gear1)
min_point_2, max_point_2 = get_minmax_point(gear2)

addendum_1 = pgw.curve_to_edges(gear1.circle_at_point(max_point_1))
dedendum_1 = pgw.curve_to_edges(gear1.circle_at_point(min_point_1))
addendum_2 = pgw.curve_to_edges(gear2.circle_at_point(max_point_2))
dedendum_2 = pgw.curve_to_edges(gear2.circle_at_point(min_point_2))


loc_gears = pgw.generate_line_of_contact(gear2, gear1, z_value)
edge_loc1 = pgw.curve_to_edges(loc_gears[0])
edge_loc2 = pgw.curve_to_edges(loc_gears[1])

edge_loc1.color = Color("red")
edge_loc2.color = Color("blue")

rb = (
    np.sin(gear1.gamma)
    * np.cos(gear1.inputparam.pressure_angle)
    * gear1.radius_spherical
)


print(f"Contact ratio: {contact_ratio:.3f}")

gears = Compound(
    children=[
        a_gear1,
        a_gear2,
        edge_loc1,
        edge_loc2,
        addendum_1[0],
        dedendum_1[0],
        addendum_2[0],
        dedendum_2[0],
    ],
    label="gears",
)


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
