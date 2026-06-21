from py_gearworks import *
from ocp_vscode import show, Animation
from copy import deepcopy
from build123d import *

# This is an attempt to create hypoid gears using modifications
# on regular bevel gears.
# Hasn't yet worked out well, but may be a starting point for future attempts.

n1, n2 = 20, 60
cone1, cone2 = cone_angle_from_teeth(n1, n2)

height = 5
beta = PI / 12

h_a = 0.8
h_d = 1.0


def conic_angle_convert(t, angle, cone_data: ConicData):

    R = cone_data.spherical_radius_untransformed
    angle2 = angle * R / (R - t)
    return angle2


gear1 = InvoluteGear(number_of_teeth=n1, height=height, cone_angle=cone1, z_anchor=0.5)
cone_1 = deepcopy(gear1.cone_data)

gear1.gearcore.shape_recipe.transform.angle = lambda t: conic_angle_convert(
    t, np.tan(beta) / (n1 / 2) * t, cone_1
)
gear1.gearcore.shape_recipe.tooth_generator.pitch_intersect_angle = (
    lambda t: conic_angle_convert(t, gear1.pitch_angle / 4, cone_1)
)
gear1.gearcore.shape_recipe.limits.h_a = lambda t: conic_angle_convert(t, h_a, cone_1)
gear1.gearcore.shape_recipe.limits.h_d = lambda t: conic_angle_convert(t, h_d, cone_1)
# gear1.calc_params()


gear2 = InvoluteGear(number_of_teeth=n2, height=height, cone_angle=cone2, z_anchor=0.5)
cone_2 = deepcopy(gear2.cone_data)
gear2.gearcore.shape_recipe.transform.angle = lambda t: conic_angle_convert(
    t, np.tan(beta) / (n2 / 2) * t, cone_2
)
gear2.gearcore.shape_recipe.tooth_generator.pitch_intersect_angle = (
    lambda t: gear2.pitch_angle / 2
    - conic_angle_convert(t, gear2.pitch_angle / 4, cone_2)
)
gear2.gearcore.shape_recipe.limits.h_a = lambda t: conic_angle_convert(t, h_a, cone_2)
gear2.gearcore.shape_recipe.limits.h_d = lambda t: conic_angle_convert(t, h_d, cone_2)
# gear2.calc_params()


a_gear1 = Compound(
    children=[gear1.build_part().solid() - gear1.face_location_top * Cylinder(2, 10)],
    label="gear1",
)
a_gear2 = Compound(
    children=[gear2.build_part().solid() - gear2.face_location_top * Cylinder(2, 10)],
    label="gear2",
)

gears = Compound(
    children=[
        a_gear1,
        a_gear2,
    ],
    label="gears",
)
gear1.mesh_to(gear2, target_dir=RIGHT)


# norm_vec = scp_Rotation.from_euler("y", gear1.gamma).apply(OUT)
# rot_beta = scp_Rotation.from_rotvec(-beta * 2 * norm_vec).as_matrix()
rot_beta = scp_Rotation.from_euler("z", -2 * beta / np.cos(gear1.gamma)).as_matrix()
gear1.gearcore.transform.orientation = rot_beta @ gear1.gearcore.transform.orientation

a_gear1.location = gear1.center_location_at_z(0)
a_gear2.location = gear2.center_location_at_z(0)


duration = 2
n = duration * 30
time_track = np.linspace(0, duration, n + 1)
gear1_track = np.linspace(0, -gear1.pitch_angle * 180 / np.pi, n + 1) * duration
gear2_track = np.linspace(0, gear2.pitch_angle * 180 / np.pi, n + 1) * duration
# assign the tracks to the gears
animation = Animation(gears)
animation.add_track("/gears/gear1", "rz", time_track, gear1_track)
animation.add_track("/gears/gear2", "rz", time_track, gear2_track)
show(gears)
animation.animate(speed=1)
