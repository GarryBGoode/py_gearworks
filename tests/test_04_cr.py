import py_gearworks as pgw
import pytest as pytest
import numpy as np


@pytest.mark.parametrize("num_teeth_1", [8, 55])
@pytest.mark.parametrize("num_teeth_2", [21, 33])
@pytest.mark.parametrize("module", [0.5, 2])  # test if module is used correctly
@pytest.mark.parametrize("beta", [0, np.pi / 6])  # spiral angle
# negative value for undercut, only for this test though
@pytest.mark.parametrize("root_fillet", [-1, 0, 0.3])
@pytest.mark.parametrize("tip_fillet", [0, 0.25])
@pytest.mark.parametrize("conic", [False, True])
def test_contact_ratio(
    num_teeth_1,
    num_teeth_2,
    module,
    beta,
    height,
    root_fillet,
    tip_fillet,
    conic,
):

    if conic:
        cone_angle_1, cone_angle_2 = pgw.cone_angle_from_teeth(num_teeth_1, num_teeth_2)
    else:
        cone_angle_1 = cone_angle_2 = 0.0

    if root_fillet < 0:
        undercut = True
        root_fillet = 0.0
    else:
        undercut = False

    gear1 = pgw.InvoluteGear(
        number_of_teeth=num_teeth_1,
        module=module,
        height=height,
        helix_angle=beta,
        root_fillet=root_fillet,
        tip_fillet=tip_fillet,
        enable_undercut=undercut,
        cone_angle=cone_angle_1,
    )
    gear2 = pgw.InvoluteGear(
        number_of_teeth=num_teeth_2,
        module=module,
        height=height,
        helix_angle=-beta,
        root_fillet=root_fillet,
        tip_fillet=tip_fillet,
        cone_angle=cone_angle_2,
        enable_undercut=undercut,
    )

    gear1.mesh_to(gear2, target_dir=pgw.LEFT, backlash=0.05, angle_bias=1)
    # for debug
    # gear1_part = gear1.build_part()
    # gear2_part = gear2.build_part()

    # loa_gears = pgw.generate_line_of_action(gear1, gear2, z_level=1)
    # edge_loa1 = pgw.curve_to_edges(loa_gears[0])
    # edge_loa2 = pgw.curve_to_edges(loa_gears[1])

    # loc_gears = pgw.generate_line_of_contact(gear2, gear1, 1)
    # edge_loc1 = pgw.curve_to_edges(loc_gears[0])
    # edge_loc2 = pgw.curve_to_edges(loc_gears[1])

    contact_ratio = pgw.get_contact_ratio_2D(gear1, gear2, z_ratio=1.0)

    # mainly looking for failure of calculation here, not value
    assert contact_ratio > 0.1


if __name__ == "__main__":
    test_contact_ratio(13, 53, 2, 0, 5, 0.2, 0.2, True)
