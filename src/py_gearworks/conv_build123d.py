# Copyright 2026 Gergely Bencsik
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


from py_gearworks.core import *
from py_gearworks.function_generators import *
from py_gearworks.curve import *
import build123d as bd
from py_gearworks.conv_spline import *
from py_gearworks.base_classes import *
from scipy.spatial.transform import Rotation as scp_Rotation
from OCP.BRepAlgoAPI import BRepAlgoAPI_BuilderAlgo
from OCP.BRepBuilderAPI import (
    BRepBuilderAPI_MakeFace,
    BRepBuilderAPI_MakeSolid,
    BRepBuilderAPI_MakeWire,
    BRepBuilderAPI_Sewing,
)
from OCP.BRepPrimAPI import BRepPrimAPI_MakePrism
from OCP.ShapeFix import ShapeFix_Face, ShapeFix_Solid
from OCP.TopAbs import TopAbs_SHELL
from OCP.TopoDS import TopoDS
from OCP.TopTools import TopTools_ListOfShape
from OCP.gp import gp_Ax3, gp_Dir, gp_Pln, gp_Pnt, gp_Sphere, gp_Vec
from typing import Callable
import dataclasses
import numpy as np
import time
import logging
import warnings
import copy


class GearBuilder(GearToNurbs):
    """A class for building build123d Part objects from gear profiles.

    The class inherits from GearToNurbs, which is responsible for generating the NURBS
    surface points and weights, this class is responsible for converting to build123d.
    Conversion happens in a reference space, with scaling of 1 (module of 1) and on the
    XY plane, default orientation. A transformation is applied after conversion to
    represent the final part.

    Parameters
    ----------
    gear : pgw.Gear
        The gear object to build.
    n_points_hz : int, optional
        Number of points used for spline approximation for each segment of the 2D gear
        profile that is not a line or an arc. Lines and arcs use exact NURB
        representation with 2 and 3 points, respectively. The default is 4.
    n_points_vert : int, optional
        Number of 2D profile slices used for generating 3D surfaces. The default is 4.
    oversampling_ratio : float, optional
        Ratio of the number of evaluations of analytical functions to the number of
        unknown points in spline approximation. Affects both horizontal points and
        vertical slices. For spline approximation, the endpoints are fixed, so the
        unkown points are the mid-points. Minimum value is 2, the default is 3. When
        fractional, the number of evaluations is rounded up.
        Example: for a 3-point spline and oversampling of 3, the unkown point is the
        middle one, the number of evaluations are the 2 end points + 3 in the middle,
        so 5 in total.
    side_surface_extension_ratio : float, optional
        Bevel gears only. Side surfaces are extended beyond the top and bottom by this
        ratio of the gear height, to be trimmed by the spherical covers. Doubled on
        each failed attempt. The default is 0.01.
    cover_extension_ratio : float, optional
        Bevel gears only. Spherical cover rings span from the outer ring (r_o) to
        beyond the tooth profile, extended by this ratio of the distance between r_o
        and the farther of the addendum and dedendum circles. They are trimmed by the
        side surfaces. The default is 0.1.
    """

    def __init__(
        self,
        gear: pgw.Gear,
        n_points_hz: int = 4,
        n_points_vert: int = 4,
        oversampling_ratio: float = 3,
        side_surface_extension_ratio: float = 0.01,
        cover_extension_ratio: float = 0.1,
    ):
        start_builder = time.time()
        if gear.cone.cone_angle == 0:
            # gear construction by creating all outside surfaces
            # and then using them to define a Solid
            super().__init__(
                gear=gear,
                n_points_hz=n_points_hz,
                n_points_vert=n_points_vert,
                oversampling_ratio=oversampling_ratio,
            )
            logging.info(
                f"Spline generation time: {time.time()-start_builder:.5f} seconds"
            )
            bot_cover = self.generate_cover(
                self.nurb_profile_stacks[0][0], self.gear_stacks[0][0]
            )
            z_bot = self.gear_stacks[0][0].transform.center[2]
            z_top = self.gear_stacks[-1][-1].transform.center[2]

            if self.is_prismatic():
                # profile does not change along z (no twist, crowning etc.),
                # the gear is a simple extrusion of the bottom cover
                time_extrude = time.time()
                prism = BRepPrimAPI_MakePrism(
                    bot_cover.wrapped, gp_Vec(0, 0, z_top - z_bot)
                )
                self.solid = bd.Solid(TopoDS.Solid(prism.Shape()))
                logging.info(f"Extrusion time: {time.time()-time_extrude:.5f} seconds")
            else:
                top_cover = self.generate_cover(
                    self.nurb_profile_stacks[-1][-1], self.gear_stacks[-1][-1]
                )
                surfaces = self.gen_side_surfaces()
                if gear.tooth_param.inside_teeth:
                    surfaces.append(self.gen_outside_ring())
                surfaces.append(bot_cover)
                surfaces.append(top_cover)

                time_solid_stitch = time.time()
                self.solid = solid_from_faces(
                    surfaces,
                    ref_face=bot_cover,
                    outward=lambda p: bd.Vector(0, 0, z_bot - z_top),
                )
                logging.info(
                    f"Solid stitching time: {time.time()-time_solid_stitch:.5f} seconds"
                )
        else:
            # gear construction by trimming over-extended side surfaces and spherical
            # cover rings against each other, then sewing them together with the
            # center faces (flat discs, or cone ring for inside teeth)
            z_vals_save = copy.deepcopy(gear.z_vals)
            zdiff = gear.z_vals[-1] - gear.z_vals[0]
            current_ratio = side_surface_extension_ratio

            # The trimming is performed by a general fuse of faces, which might fail
            # on edge cases. A few different extension ratios are attempted.
            for attempt in range(4):  # initial attempt + up to 3 retries
                # extend z_vals to ensure side surfaces cross the covers.
                gear.z_vals = copy.deepcopy(z_vals_save)
                gear.z_vals[-1] += current_ratio * zdiff
                gear.z_vals[0] -= current_ratio * zdiff
                super().__init__(
                    gear=gear,
                    n_points_hz=n_points_hz,
                    n_points_vert=n_points_vert,
                    oversampling_ratio=oversampling_ratio,
                )
                # restore original z_vals
                self.gear.z_vals = copy.deepcopy(z_vals_save)
                side_surfaces = self.gen_side_surfaces()

                start_trim = time.time()
                try:
                    # cover parameters depend on accurate (original) z_vals
                    self.solid = self.gen_bevel_solid(
                        side_surfaces, cover_extension_ratio
                    )
                    logging.info(
                        f"Trimming and stitching time: "
                        f"{time.time()-start_trim:.5f} seconds"
                    )
                    break
                except RuntimeError as err:
                    if attempt < 3:
                        warnings.warn(
                            f"Bevel gear construction failed (attempt {attempt + 1}): "
                            f"{err} Retrying with extension ratio "
                            f"{current_ratio * 2:.4f}.",
                            RuntimeWarning,
                            stacklevel=2,
                        )
                        current_ratio *= 2
            else:
                raise RuntimeError(
                    "Trimming gear surfaces by spherical covers failed "
                    "after 4 attempts."
                )

        self.part = bd.Part() + self.solid
        self.part_transformed = apply_transform_part(self.part, self.gear.transform)

        logging.info(
            f"Total time for generation: {time.time()-start_builder:.5f} seconds"
        )

    def gen_bevel_solid(
        self, side_surfaces: list[bd.Face], cover_extension_ratio: float
    ) -> bd.Solid:
        """Builds the Solid of a bevel gear from the over-extended side surfaces.

        Spherical cover rings (bottom and top) and side surfaces are trimmed against
        each other via a general fuse, which only intersects faces instead of solids.
        The trimmed faces are sewn together with the center faces: flat discs for
        outside teeth, a cone ring for inside teeth."""
        gearcopy = copy.deepcopy(self.gear)
        gearcopy.transform = GearTransform()
        covers = [
            SphereCover.from_gear(gearcopy, z, cover_extension_ratio)
            for z in (gearcopy.z_vals[0], gearcopy.z_vals[-1])
        ]
        cover_bot, cover_top = covers
        cover_faces = [cover.make_face() for cover in covers]

        if self.gear.tooth_param.inside_teeth:
            height = cover_top.z_o - cover_bot.z_o
            cone = bd.Solid.make_cone(
                cover_bot.r_o,
                cover_top.r_o,
                np.abs(height),
                plane=bd.Plane(
                    origin=(0, 0, cover_bot.z_o),
                    x_dir=cover_bot.x_dir,
                    z_dir=(0, 0, np.sign(height)),
                ),
            )
            center_faces = [cone.faces().sort_by(bd.Axis.Z)[1]]
            # the body is between the cone and the teeth, outward is away from teeth
            radial_sign = np.sign(cover_bot.r_o - cover_bot.r_end)

            def outward(p: bd.Vector):
                return bd.Vector(p.X, p.Y, 0) * radial_sign

        else:
            center_faces = [cover.make_disc() for cover in covers]

            def outward(p: bd.Vector):
                return bd.Vector(0, 0, cover_bot.z_o - cover_top.z_o)

        fuse = BRepAlgoAPI_BuilderAlgo()
        arguments = TopTools_ListOfShape()
        # side surfaces go in a single compound, so they are not intersected with
        # each other, only with the covers
        for shape in [bd.Compound(side_surfaces), *cover_faces, *center_faces]:
            arguments.Append(shape.wrapped)
        fuse.SetArguments(arguments)
        fuse.SetRunParallel(True)
        fuse.Build()
        if not fuse.IsDone():
            raise RuntimeError("General fuse of gear faces failed.")

        tol = 1e-5 * max(cover.radius for cover in covers)
        faces = []
        for face in side_surfaces:
            # drop the over-extended parts of side surfaces
            faces.extend(
                fragment
                for fragment in fuse_fragments(fuse, face)
                if not any(
                    is_beyond_covers(vertex, covers, tol)
                    for vertex in fragment.vertices()
                )
            )
        for cover, cover_face in zip(covers, cover_faces):
            # the ring is split into the part bounded by the r_o circle and the tooth
            # profile, and the over-extended part beyond the tooth profile
            ring_fragments = [
                fragment
                for fragment in fuse_fragments(fuse, cover_face)
                if any(cover.is_on_r_o(vertex, tol) for vertex in fragment.vertices())
            ]
            if len(ring_fragments) != 1:
                raise RuntimeError("Side surfaces did not cleanly split the cover.")
            faces.extend(ring_fragments)
        center_fragments = [
            fragment for face in center_faces for fragment in fuse_fragments(fuse, face)
        ]
        faces.extend(center_fragments)

        return solid_from_faces(faces, ref_face=center_fragments[0], outward=outward)

    def gen_side_surfaces(self):
        surface_gen_start = time.time()
        n_teeth = self.gear.tooth_param.num_teeth_act
        surfaces = []

        patches_z = [
            [*surfdata_z.get_patches()][:-3] for surfdata_z in self.side_surf_data
        ][: len(self.gear.z_vals) - 1]

        for j in range(n_teeth):
            # rotating the control points is much cheaper than rotating the Face
            rot = rot_z(self.gear.tooth_param.pitch_angle * j)
            for patches in patches_z:
                for patch in patches:
                    # shape: vert x horiz x xyz
                    points = patch["points"] @ rot.T
                    weights = patch["weights"]
                    vpoints = [
                        nppoint2Vector(points[i]) for i in range(points.shape[0])
                    ]
                    face = bd.Face.make_bezier_surface(vpoints, weights.tolist())
                    surfaces.append(face)

        logging.info(
            f"Surface generation time: {time.time()-surface_gen_start:.5f} seconds"
        )
        return surfaces

    def gen_outside_ring(self):
        r_o = -self.gear.shape_recipe.limits.h_o + self.gear.tooth_param.num_teeth / 2
        ring_base = bd.Edge.make_circle(radius=r_o, plane=bd.Plane.XY)

        edge_ring = bd.Edge.make_line(
            bd.Vector((r_o, 0, self.gear.z_vals[0])),
            bd.Vector((r_o, 0, self.gear.z_vals[-1])),
        )

        ring_surf = bd.Face.sweep(profile=edge_ring, path=ring_base)
        return ring_surf

    def generate_cover(
        self, nurb_stack: GearRefProfileExtended, gear_stack: GearRefProfileExtended
    ):

        if self.gear.cone.cone_angle != 0:

            if not self.gear.tooth_param.inside_teeth:
                curve = crv.NURBSCurve.from_curve_chain(nurb_stack.tooth_profile_closed)
                splines = gen_splines(curve)

                face_tooth = bd.Face.make_surface(bd.Wire(splines))

                cover_edge_curve = crv.NURBSCurve(
                    nurb_stack.rd_curve, nurb_stack.rd_connector
                )
                cover_edge = bd.Edge() + gen_splines(cover_edge_curve)

                num_teeth = self.gear.tooth_param.num_teeth_act
                cover_edge = cover_edge + [
                    cover_edge.rotate(
                        axis=bd.Axis.Z,
                        angle=(j + 1) * nurb_stack.pitch_angle * 180 / PI,
                    )
                    for j in range(num_teeth - 1)
                ]
                cover_face = bd.Face(bd.Wire(cover_edge))

                out_face = cover_face + [
                    face_tooth.rotate(
                        axis=bd.Axis.Z, angle=j * nurb_stack.pitch_angle * 180 / PI
                    )
                    for j in range(num_teeth)
                ]
                return out_face
            else:
                curve = crv.NURBSCurve.from_curve_chain(nurb_stack.profile_closed)
                splines = gen_splines(curve)
                face_tooth = bd.Face.make_surface(bd.Wire(splines))
                num_teeth = self.gear.tooth_param.num_teeth_act
                face_fuse_time = time.time()
                out_face = bd.Face.fuse(
                    *[
                        face_tooth.rotate(
                            axis=bd.Axis.Z, angle=j * nurb_stack.pitch_angle * 180 / PI
                        )
                        for j in range(num_teeth)
                    ]
                )
                logging.info(
                    f"Face fuse time: {time.time()-face_fuse_time:.5f} seconds"
                )
                return out_face
        else:

            num_teeth = self.gear.tooth_param.num_teeth_act
            curve = cover_profile_curve(nurb_stack)
            profile_edge_time = time.time()
            # rotating control points and joining all edges once is much cheaper
            # than rotating and fusing Wires of each tooth
            edges = []
            for j in range(num_teeth):
                rot = rot_z(nurb_stack.pitch_angle * j)
                for sub_curve in curve.get_curves():
                    edges.append(
                        bd.Edge.make_bezier(
                            *nppoint2Vector(sub_curve.points @ rot.T),
                            weights=sub_curve.weights.tolist(),
                        )
                    )
            profile_wire = make_wire_ordered(edges)
            logging.info(
                f"Profile edge generation time: {time.time()-profile_edge_time:.5f} seconds"
            )

            face_fuse_time = time.time()
            z_val = gear_stack.transform.center[2]

            if self.gear.tooth_param.inside_teeth:
                r_o = (
                    -self.gear.shape_recipe.limits.h_o
                    + self.gear.tooth_param.num_teeth / 2
                )
                ring = bd.Edge.make_circle(radius=r_o, plane=bd.Plane.XY.offset(z_val))
                face = make_planar_face(
                    z_val, bd.Wire(ring), inner_wires=[profile_wire]
                )
            else:
                face = make_planar_face(z_val, profile_wire)
            logging.info(f"Face fuse time: {time.time()-face_fuse_time:.5f} seconds")
            return face

    def is_prismatic(self, tol=1e-7):
        """Returns True if the bottom and top 2D profiles are identical, meaning the
        gear can be generated by extruding the bottom cover."""
        curves_bot = cover_profile_curve(self.nurb_profile_stacks[0][0]).get_curves()
        curves_top = cover_profile_curve(self.nurb_profile_stacks[-1][-1]).get_curves()
        if len(curves_bot) != len(curves_top):
            return False
        for c_bot, c_top in zip(curves_bot, curves_top):
            if c_bot.points.shape != c_top.points.shape:
                return False
            if not np.allclose(c_bot.points[:, :2], c_top.points[:, :2], atol=tol):
                return False
            if not np.allclose(c_bot.weights, c_top.weights, atol=tol):
                return False
        return True


class GearBuilder_old(GearToNurbs):
    """A class for building Part objects from gear profiles."""

    def __init__(
        self,
        gear: pgw.Gear,
        n_points_hz=4,
        n_points_vert=4,
        oversampling_ratio=2.5,
        add_plug=False,
    ):
        super().__init__(
            gear=gear,
            n_points_hz=n_points_hz,
            n_points_vert=n_points_vert,
            oversampling_ratio=oversampling_ratio,
        )
        surfaces = []
        ro_surfaces = []

        start = time.time()
        for k in range(len(self.gear.z_vals) - 1):
            surfdata_z = self.side_surf_data[k]

            for patch in surfdata_z.get_patches():
                points = patch["points"]
                weights = patch["weights"]
                vpoints = [nppoint2Vector(points[k]) for k in range(points.shape[0])]
                surfaces.append(bd.Face.make_bezier_surface(vpoints, weights.tolist()))
            ro_surfaces.append(surfaces[-2])
        self.surfaces = surfaces
        top_points, top_weights = (
            self.side_surf_data[-1].points[-1, :, :],
            self.side_surf_data[-1].weights[-1, :],
        )
        top_curve = crv.NURBSCurve.from_points(
            top_points, knots=self.side_surf_data[-1].knots, weights=top_weights
        )
        splines = [self.gen_splines(curve) for curve in top_curve.get_curves()]
        top_surface = bd.Face.make_surface(bd.Wire(splines))

        bot_points, bot_weights = (
            self.side_surf_data[0].points[0, :, :],
            self.side_surf_data[0].weights[0, :],
        )
        bot_curve = crv.NURBSCurve.from_points(
            bot_points, knots=self.side_surf_data[0].knots, weights=bot_weights
        )
        splines = [self.gen_splines(curve) for curve in bot_curve.get_curves()]
        bot_surface = bd.Face.make_surface(bd.Wire(splines))

        if len(ro_surfaces) > 1:
            ro_surface = bd.Face.fuse(*ro_surfaces)
        else:
            ro_surface = ro_surfaces[0]
        ro_spline_top = self.gen_splines(top_curve.get_curves()[-2])
        ro_spline_bot = self.gen_splines(bot_curve.get_curves()[-2])
        surfaces.insert(0, bot_surface)
        surfaces.append(top_surface)
        shell = bd.Shell(surfaces)
        solid1 = bd.Solid(shell)
        solid1 = fix_attempt(solid1)

        logging.log(
            logging.INFO, f"Gear 1-tooth solid build time: {time.time()-start:.5f}"
        )
        start = time.time()

        self.profile_solid = solid1

        n_teeth = self.gear.tooth_param.num_teeth_act
        bin_n_teeth = bin(n_teeth)[2:]
        shape_dict = []
        solid2_to_fuse = []
        angle_construct = 0.0
        angle_idx = 0
        tol = 1e-4

        axis1 = bd.Axis.Z

        for k in range(len(bin_n_teeth)):

            if k == 0:
                shape_dict.append(solid1)
                angle = 0
            else:
                angle = self.gear.tooth_param.pitch_angle * RAD2DEG * (2 ** (k - 1))
                rotshape = shape_dict[k - 1].rotate(axis1, angle)
                fuse_shape = (
                    shape_dict[k - 1].fuse(rotshape, glue=False, tol=tol).clean()
                )
                fuse_shape = fix_attempt(fuse_shape)
                shape_dict.append(fuse_shape)

            if bin_n_teeth[-(k + 1)] == "1":

                angle_construct = (
                    angle_idx * self.gear.tooth_param.pitch_angle * RAD2DEG
                )
                rotshape = shape_dict[k].rotate(axis1, angle_construct)

                solid2_to_fuse.append(rotshape)
                angle_idx = angle_idx + 2**k

        if len(solid2_to_fuse) > 1:
            self.solid = bd.Solid.fuse(*solid2_to_fuse, glue=False, tol=tol).clean()
        else:
            self.solid = solid2_to_fuse[0].clean()

        self.solid = fix_attempt(self.solid)

        plug_surfaces = []
        plug_splines_top = []
        plug_splines_bot = []
        if add_plug:
            for k in range(n_teeth):
                plug_surfaces.append(
                    ro_surface.rotate(
                        axis1, self.gear.tooth_param.pitch_angle * RAD2DEG * k
                    )
                )
                plug_splines_bot.append(
                    ro_spline_bot.rotate(
                        axis1, self.gear.tooth_param.pitch_angle * RAD2DEG * k
                    )
                )
                plug_splines_top.append(
                    ro_spline_top.rotate(
                        axis1, self.gear.tooth_param.pitch_angle * RAD2DEG * k
                    )
                )
            plug_top = bd.Face.make_surface(bd.Wire(plug_splines_top))
            plug_bot = bd.Face.make_surface(bd.Wire(plug_splines_bot))
            plug_surfaces.insert(0, plug_bot)
            plug_surfaces.append(plug_top)
            plug = bd.Solid(bd.Shell(plug_surfaces))
            plug = fix_attempt(plug)
            self.solid = self.solid.fuse(plug).clean()
            self.solid = fix_attempt(self.solid)

        logging.log(
            logging.INFO, f"Gear solid fuse time: {time.time()-start:.5f} seconds"
        )
        self.solid = bd.BasePartObject(self.solid).fix()
        self.solid_transformed = apply_transform_part(self.solid, self.gear.transform)
        self.part_transformed = self.solid_transformed

    def gen_splines(self, curve_bezier: Curve):
        vectors = nppoint2Vector(curve_bezier.points)
        weights = curve_bezier.weights.tolist()
        return bd.Edge.make_bezier(*vectors, weights=weights)


def apply_transform_part(part: bd.Part, transform: GearTransform):
    location1 = transform2Location(transform)
    part = part.scale(transform.scale)
    part = location1 * part
    return part


def apply_animation(gear: pgw.Gear, part: bd.Part, time: float = 1):
    pass


def fix_attempt(solid):
    if not solid.is_valid():
        warnings.warn("Invalid solid found", RuntimeWarning, stacklevel=2)
        solid = solid.fix()
    return solid


def nppoint2Vector(p: np.ndarray):
    if p.size == 3:
        return bd.Vector((p[0], p[1], p[2]))
    else:
        return [bd.Vector((p[k, 0], p[k, 1], p[k, 2])) for k in range(p.shape[0])]


def np2v(p: np.ndarray):
    # shorthand for npppoint2Vector
    return nppoint2Vector(p)


def gen_splines(curve_bezier: Curve):
    if isinstance(curve_bezier, NURBSCurve) or isinstance(curve_bezier, CurveChain):
        splines = []
        for curve in curve_bezier.get_curves():
            if curve.active:
                vectors = nppoint2Vector(curve.points)
                weights = curve.weights.tolist()
                splines.append(bd.Edge.make_bezier(*vectors, weights=weights))
        return splines
    else:
        vectors = nppoint2Vector(curve_bezier.points)
        weights = curve_bezier.weights.tolist()
        return bd.Edge.make_bezier(*vectors, weights=weights)


def rot_z(angle: float):
    """Rotation matrix around the Z axis."""
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def cover_profile_curve(nurb_stack: GearRefProfileExtended):
    """NURBS curve of one tooth of the 2D profile, used for flat covers."""
    curve = crv.NURBSCurve.from_curve_chain(nurb_stack.profile)
    curve.del_inactive_curves()
    curve.enforce_continuity()
    return curve


def make_wire_ordered(edges: list[bd.Edge]):
    """Joins edges that are already ordered head-to-tail into a Wire.
    Skips the edge sorting and cleaning of bd.Wire(), falls back to it on failure."""
    wire_builder = BRepBuilderAPI_MakeWire()
    for edge in edges:
        wire_builder.Add(edge.wrapped)
    if wire_builder.IsDone():
        return bd.Wire(wire_builder.Wire())
    return bd.Wire(edges)


def make_planar_face(z: float, outer_wire: bd.Wire, inner_wires=()):
    """Makes a face on the XY plane offset to z. Skips the planarity search and
    wire fixing of bd.Face(), since the plane is known."""
    plane = gp_Pln(gp_Pnt(0, 0, z), gp_Dir(0, 0, 1))
    face_builder = BRepBuilderAPI_MakeFace(plane, outer_wire.wrapped, True)
    if not inner_wires:
        return bd.Face(face_builder.Face())
    for inner_wire in inner_wires:
        face_builder.Add(inner_wire.wrapped)
    # inner wires need to be oriented opposite to outer wire
    face_fix = ShapeFix_Face(face_builder.Face())
    face_fix.FixOrientation()
    face_fix.Perform()
    return bd.Face(TopoDS.Face(face_fix.Result()))


@dataclasses.dataclass
class SphereCover:
    """Spherical ring of a bevel gear's top or bottom cover. It spans from the r_o
    circle (inner or outer ring of the tooth profile) to beyond the tooth profile. The
    sphere is centered on the Z axis, latitudes are measured from the sphere center's XY
    plane.
    """

    center: np.ndarray
    radius: float
    angle: float
    lat_o: float
    lat_end: float

    @classmethod
    def from_gear(cls, gear: pgw.Gear, z: float, extension_ratio: float):
        """Cover of the gear at z. The gear should have default (identity)
        transform."""
        center, radius = gear.sphere_data_at_z(z)
        profile = gear.curve_gen_at_z(z)

        def latitude(p):
            d = p - center
            return np.arctan2(d[2], np.hypot(d[0], d[1]))

        lat_o = latitude(profile.transform(profile.ro_curve(0)))
        # addendum or dedendum circle, whichever is farther from r_o
        lat_far = max(
            latitude(profile.transform(profile.ra_curve(0.5))),
            latitude(profile.transform(profile.rd_curve(0.5))),
            key=lambda lat: np.abs(lat - lat_o),
        )
        lat_end = np.clip(
            lat_far + extension_ratio * (lat_far - lat_o), -PI / 2, PI / 2
        )
        return cls(
            center=center,
            radius=np.abs(radius),
            angle=gear.shape_recipe(z).transform.angle,
            lat_o=lat_o,
            lat_end=lat_end,
        )

    @property
    def r_o(self):
        """Radius of the r_o circle."""
        return self.radius * np.cos(self.lat_o)

    @property
    def r_end(self):
        """Radius of the ring's end circle, beyond the tooth profile."""
        return self.radius * np.cos(self.lat_end)

    @property
    def z_o(self):
        """Z coordinate of the r_o circle."""
        return self.center[2] + self.radius * np.sin(self.lat_o)

    @property
    def x_dir(self):
        """Direction of the seam, rotated along with the tooth profile."""
        return (np.cos(self.angle), np.sin(self.angle), 0)

    def make_face(self) -> bd.Face:
        axis = gp_Ax3(gp_Pnt(*self.center), gp_Dir(0, 0, 1), gp_Dir(*self.x_dir))
        face_builder = BRepBuilderAPI_MakeFace(
            gp_Sphere(axis, self.radius),
            0,
            2 * PI,
            min(self.lat_o, self.lat_end),
            max(self.lat_o, self.lat_end),
        )
        return bd.Face(face_builder.Face())

    def make_disc(self) -> bd.Face:
        """Flat disc bounded by the r_o circle."""
        plane = bd.Plane(origin=(0, 0, self.z_o), x_dir=self.x_dir, z_dir=(0, 0, 1))
        circle = bd.Edge.make_circle(self.r_o, plane=plane)
        return make_planar_face(self.z_o, bd.Wire(circle))

    def signed_distance(self, p: bd.Vertex):
        """Positive outside of the sphere, negative inside."""
        return np.linalg.norm(np.array([p.X, p.Y, p.Z]) - self.center) - self.radius

    def is_on_r_o(self, p: bd.Vertex, tol: float):
        return (
            np.abs(np.hypot(p.X, p.Y) - self.r_o) < tol and np.abs(p.Z - self.z_o) < tol
        )


def is_beyond_covers(p: bd.Vertex, covers: list[SphereCover], tol: float):
    """True if p is clearly outside the region between the 2 (concentric) covers.
    Inside the region p is outside of one sphere and inside the other."""
    d0 = covers[0].signed_distance(p)
    d1 = covers[1].signed_distance(p)
    return min(np.abs(d0), np.abs(d1)) > tol and np.sign(d0) == np.sign(d1)


def fuse_fragments(fuse: BRepAlgoAPI_BuilderAlgo, face: bd.Face) -> list[bd.Face]:
    """Fragments of an input face after a general fuse operation."""
    if fuse.IsDeleted(face.wrapped):
        return []
    modified = fuse.Modified(face.wrapped)
    if modified.IsEmpty():
        return [face]
    return [bd.Face(TopoDS.Face(fragment)) for fragment in modified]


def solid_from_faces(
    faces: list[bd.Face],
    ref_face: bd.Face,
    outward: Callable[[bd.Vector], bd.Vector],
):
    """Sews faces into a closed Solid. Orientation of the Solid is set via ref_face,
    whose outward direction at a point p is known to be outward(p).
    This avoids the costly point classification of ShapeFix_Solid.
    Raises RuntimeError if the faces don't form a single closed shell."""
    sewing = BRepBuilderAPI_Sewing()
    for face in faces:
        sewing.Add(face.wrapped)
    sewing.Perform()
    sewed_shape = sewing.SewedShape()
    if sewed_shape.ShapeType() != TopAbs_SHELL or sewing.NbFreeEdges() > 0:
        raise RuntimeError(
            f"Sewing faces did not result in a closed shell, "
            f"found {sewing.NbFreeEdges()} free edges."
        )
    shell = TopoDS.Shell(sewed_shape)
    solid = BRepBuilderAPI_MakeSolid(shell).Solid()

    ref_sewn = sewing.Modified(ref_face.wrapped)
    for face in bd.Solid(solid).faces():
        if face.wrapped.IsSame(ref_sewn):
            normal = face.normal_at(0.5, 0.5)
            if normal.dot(outward(face.position_at(0.5, 0.5))) < 0:
                solid.Reverse()
            return bd.Solid(solid)

    # reference face not found, fall back to the slow but robust method
    return bd.Solid(ShapeFix_Solid().SolidFromShell(shell))


def transform2Location(transform: GearTransform):
    rot1 = scp_Rotation.from_matrix(transform.orientation)
    degrees = rot1.as_euler("zyx", degrees=True)
    loc = bd.Location(
        transform.center,
        [degrees[0] + transform.angle * 180 / PI, degrees[1], degrees[2]],
        bd.Extrinsic.ZYX,
    )

    return loc


def generate_boundary_edges(
    nurbprofile: GearRefProfile,
    transform: GearTransform = None,
    angle_range: float = 2 * PI,
):
    if transform is None:
        # identity transform by default
        transform = GearTransform()
    nurb_profile = gearprofile_to_nurb(nurbprofile)
    # don't want to get more inputs about num of teeth, but without rounding this can
    # lose 1 tooth
    N = int(np.round(angle_range / nurbprofile.pitch_angle))

    curves = []
    for i in range(N):
        # angle = i * profile.pitch_angle
        curves.extend(
            [
                nurb.apply_transform(transform)
                for nurb in nurb_profile.profile.copy().get_curves()
            ]
        )
        transform.angle += nurbprofile.pitch_angle

    nurbs_curve = crv.NURBSCurve(*curves)
    nurbs_curve.enforce_continuity()

    return gen_splines(nurbs_curve)


def arc_to_b123d(arc: crv.ArcCurve) -> bd.Edge:
    """Converts a py_gearworks ArcCurve to a build123d Edge object."""
    if (arc.t_0, arc.t_1) != (0, 1):
        # arc can be extended, better make a new one
        arc2 = crv.ArcCurve.from_2_point_center(arc.p0, arc.p1, arc.center)
    else:
        arc2 = arc

    loc = bd.Location(
        arc2.center,
        [arc2.roll * 180 / PI, arc2.pitch * 180 / PI, arc2.yaw * 180 / PI],
        bd.Intrinsic.XYZ,
    )

    if arc2.angle < 0:
        start = arc2.angle * 180 / PI
        end = 0
    else:
        start = 0
        end = arc2.angle * 180 / PI

    bd_arc = bd.Edge.make_circle(
        radius=arc2.radius,
        plane=bd.Plane(loc),
        start_angle=start,
        end_angle=end,
    )
    return bd_arc


def line_to_b123d(line: crv.LineCurve) -> bd.Edge:
    """Converts a py_gearworks LineCurve to a build123d Edge object."""
    return bd.Edge.make_line(np2v(line.p0), np2v(line.p1))


def curve_to_edges(curve: crv.Curve):
    if isinstance(curve, crv.CurveChain):
        return [curve_to_edges(curve) for curve in curve.get_curves()]
    elif isinstance(curve, crv.NURBSCurve) | isinstance(curve, crv.NurbCurve):
        return gen_splines(curve)
    elif isinstance(curve, crv.ArcCurve):
        return [arc_to_b123d(curve)]
    elif isinstance(curve, crv.LineCurve):
        return [line_to_b123d(curve)]
    elif isinstance(curve, crv.TransformedCurve):
        if curve.t_0 != 0 or curve.t_1 != 1:
            nurb = crv.convert_curve_nurbezier(curve)
            return gen_splines(nurb)
        else:
            nurb = crv.convert_curve_nurbezier(curve.target_curve)
            nurb.apply_transform(curve.transform_method)
            return gen_splines(nurb)
    else:
        nurb = crv.convert_curve_nurbezier(curve)
        return gen_splines(nurb)
