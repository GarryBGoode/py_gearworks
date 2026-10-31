from py_gearworks import *
import numpy as np
from ocp_vscode import show
from scipy.spatial.transform import Rotation
from scipy.optimize import brentq, root, minimize_scalar
from matplotlib import pyplot as plt
from functools import partial
import shapely
from shapely.geometry import Polygon

module = 2
pitch = module * np.pi
addendum_coefficient = 1.0
dedendum_coefficient = 1.25
pressure_angle = np.radians(15)

n_teeth = 35
r_pitch = module * n_teeth / 2


def make_cut_profile(
    addendum_coefficient, dedendum_coefficient, pitch, pressure_angle, module
):
    add = addendum_coefficient * module
    ded = dedendum_coefficient * module
    ref_point_a_0 = RIGHT * add
    ref_point_a = RIGHT * add + DOWN * (pitch / 4 - np.tan(pressure_angle) * add)
    ref_point_d = LEFT * ded + DOWN * (pitch / 4 + np.tan(pressure_angle) * ded)
    ref_point_c = DOWN * pitch / 2 + LEFT * ded

    ref_point_mid = (ref_point_a * 0.75 + ref_point_d * 0.25) + DOWN * 0.2 + RIGHT * 0.2

    # ref_point_c_mirror = ref_point_c * np.array([1, -1, 1])

    ref_point_e = ref_point_c + RIGHT * (1 + add + ded)
    ref_point_f = ref_point_a_0 + RIGHT * (1)

    ref_line_crv = PolyLineCurve(
        [
            ref_point_a_0,
            ref_point_a,
            ref_point_mid,
            ref_point_d,
            ref_point_c,
            # ref_point_e,
            # ref_point_f,
            # ref_point_a_0,
        ]
    )
    ref_line_crv = ref_line_crv.fillet_at_index(0.4, 2)

    return ref_line_crv


def make_wavy_cut_profile(
    addendum_coefficient,
    dedendum_coefficient,
    pitch,
    pressure_angle,
    module,
    sagitta=0.05,
    amplitude=0.03,
    periods=3.5,
):
    """Stress test cutter: rack flank is a slight arc with a sine wave superimposed.

    sagitta: height of the arc over the straight flank, in module units,
        positive bulges into the gap (away from the cutter material)
    amplitude: sine wave amplitude, in module units
    periods: number of sine periods along the flank, half-integers keep the
        flank end points on the straight rack corners
    """
    add = addendum_coefficient * module
    ded = dedendum_coefficient * module
    ref_point_a_0 = RIGHT * add
    ref_point_a = RIGHT * add + DOWN * (pitch / 4 - np.tan(pressure_angle) * add)
    ref_point_d = LEFT * ded + DOWN * (pitch / 4 + np.tan(pressure_angle) * ded)
    ref_point_c = DOWN * pitch / 2 + LEFT * ded

    chord = ref_point_d - ref_point_a
    chord_len = np.linalg.norm(chord)
    # away from the cutter material
    normal_gap = -normal_of(chord / chord_len)
    sag = sagitta * module
    amp = amplitude * module
    if sag != 0:
        radius = (chord_len**2 / 4 + sag**2) / (2 * sag)

    def flank_func(u):
        offset = amp * np.sin(2 * np.pi * periods * u)
        if sag != 0:
            x = (u - 0.5) * chord_len
            offset += np.sign(radius) * np.sqrt(radius**2 - x**2) - (radius - sag)
        return ref_point_a + chord * u + normal_gap * offset

    return CurveChain(
        LineCurve(ref_point_a_0, ref_point_a),
        Curve(flank_func),
        LineCurve(ref_point_d, ref_point_c),
    )


path_curve_offset = 0  # profile shift
# outer radius of the blank, nothing beyond it gets cut
r_blank = r_pitch + addendum_coefficient * module + path_curve_offset


def normal_of(tangent):
    # rotated +90deg, points into the cutter material (material left of travel)
    normal = np.array([-tangent[1], tangent[0], 0])
    return normal / np.linalg.norm(normal)


def segment_point_normal(segment: Curve, t):
    # one-sided derivative at the ends, stay within the segment
    direction = 1 if t == 0 else (-1 if t == 1 else 0)
    return segment(t), normal_of(segment.derivative(t, direction=direction))


def calc_roll_angle(CutPoint, CutNormal):
    """Roll angle where the normal line of CutPoint hits the pole on the rolling line.

    Returns nan if the normal is parallel to the rolling line (point never cuts).
    """
    if abs(CutNormal[0]) < 1e-12:
        return np.nan
    # distance along normal to rolling line (x = -profile shift)
    lambda_n = (-path_curve_offset - CutPoint[0]) / CutNormal[0]
    y_center = CutPoint[1] + lambda_n * CutNormal[1]
    # roll angle is negative if y displacement is positive
    return -y_center / r_pitch


def transform_to_blank(CutPoint, Phi):
    """Map a cutter-frame point into the blank frame at roll angle Phi.

    Result is shifted by -r_pitch in x so it overlays the cutter in plots.
    """
    blank_center = -RIGHT * (r_pitch + path_curve_offset) + UP * (-r_pitch * Phi)
    R = Rotation.from_euler("z", -Phi)
    return R.apply(-blank_center + CutPoint) - RIGHT * r_pitch


def blank_radius_of(SweptPoint):
    # undo the plotting shift of transform_to_blank()
    return np.linalg.norm((SweptPoint + RIGHT * r_pitch)[:2])


def point_phi_range(CutPoint):
    """Roll angle range where a fixed cutter point is within r_blank.

    |P|^2 = (r_pitch + shift + x)^2 + (y + r_pitch * phi)^2 <= r_blank^2
    Returns None if the point never reaches the blank.
    """
    dx = r_pitch + path_curve_offset + CutPoint[0]
    if abs(dx) > r_blank:
        return None
    h = np.sqrt(r_blank**2 - dx**2)
    return (-CutPoint[1] - h) / r_pitch, (-CutPoint[1] + h) / r_pitch


def calc_angle_limits(cut_profile: CurveChain, n_samples=100, margin=0.05):
    """Roll angle range that covers the complete cut.

    Outside of it no point of the cutter is within the blank.
    Points are sampled per segment, margin (fraction of the range) covers sampling.
    """
    ranges = []
    for segment in cut_profile.get_curves():
        for t in np.linspace(0, 1, n_samples):
            result = point_phi_range(segment(t))
            if result is not None:
                ranges.append(result)
    ranges = np.array(ranges)
    lo, hi = np.min(ranges[:, 0]), np.max(ranges[:, 1])
    return np.array([lo, hi]) + np.array([-1, 1]) * margin * (hi - lo)


def find_roots(func, t_scan, values):
    """Roots of func on the scanned interval, bracketed by sign changes."""
    roots = []
    for i in range(len(t_scan) - 1):
        if values[i] == 0:
            roots.append(t_scan[i])
        elif values[i] * values[i + 1] < 0:
            roots.append(brentq(func, t_scan[i], t_scan[i + 1], xtol=1e-14))
    if values[-1] == 0:
        roots.append(t_scan[-1])
    return roots


def envelope_point(segment: Curve, t):
    CutPoint, CutNormal = segment_point_normal(segment, t)
    Phi = calc_roll_angle(CutPoint, CutNormal)
    if np.isnan(Phi):
        return np.nan, np.full(3, np.nan)
    return Phi, transform_to_blank(CutPoint, Phi)


def clip_at_blank_radius(segment: Curve, t_root, t_other, n_scan=100):
    """Parameter between t_other and an asymptote at t_root where the envelope leaves r_blank.

    Near the asymptote phi -> +-inf, so the envelope radius grows without bound.
    Returns None if the envelope is outside the blank on the whole interval.
    """

    def radius_error(t):
        return blank_radius_of(envelope_point(segment, t)[1]) - r_blank

    # stop just short of the asymptote, the root itself has no envelope point
    t_near = t_root + (t_other - t_root) * 1e-12
    t_scan = np.linspace(t_other, t_near, n_scan)
    values = np.array([radius_error(t) for t in t_scan])
    inside = np.nonzero(values < 0)[0]
    if len(inside) == 0:
        return None
    j = inside[-1]
    if j == n_scan - 1:
        return t_near
    return brentq(radius_error, t_scan[j], t_scan[j + 1], xtol=1e-14)


def gen_point_trajectory(CutPoint, n_points):
    """A single cutter point cutting over its whole roll angle range within the blank."""
    result = point_phi_range(CutPoint)
    if result is None:
        return None
    phis = np.linspace(*result, n_points)
    return phis, np.array([transform_to_blank(CutPoint, Phi) for Phi in phis])


def gen_segment_points(
    segment: Curve,
    n_points,
    n_scan=100,
    tol_parallel=1e-9,
    tol_on_line=1e-9,
    tol_end=1e-5,
):
    """Envelope of a smooth cutter segment.

    Normals parallel to the rolling line (n_x = 0) split the segment into pieces:
        - point off the rolling line: no contact, phi -> +-inf at the split,
          the piece is clipped where the envelope leaves r_blank
        - point on the rolling line: the normal contains the pole for every phi,
          this point cuts on its own like a vertex
    A straight segment with n_x = 0 everywhere only cuts with its point on the
    rolling line.

    Returns a list of (kind, param, phis, points), kind is "segment" or "point",
    param is segment t (n_points for segment pieces, single t for points).
    """
    branches = []

    def normal_x(t):
        return segment_point_normal(segment, t)[1][0]

    def offset_x(t):
        return segment(t)[0] + path_curve_offset

    def add_point_branch(t):
        result = gen_point_trajectory(segment(t), n_points)
        if result is not None:
            branches.append(("point", t, *result))

    t_scan = np.linspace(0, 1, n_scan)
    nx_values = np.array([normal_x(t) for t in t_scan])
    if np.all(np.abs(nx_values) < tol_parallel):
        offset_values = np.array([offset_x(t) for t in t_scan])
        for t_root in find_roots(offset_x, t_scan, offset_values):
            add_point_branch(t_root)
        return branches

    roots = find_roots(normal_x, t_scan, nx_values)
    # roots at the segment ends show no sign change, and are blurred by the
    # numeric derivative noise (~1e-6), e.g. a fillet starting on a radial wall
    if abs(nx_values[0]) < tol_end and (not roots or roots[0] != 0):
        roots.insert(0, 0.0)
    if abs(nx_values[-1]) < tol_end and (not roots or roots[-1] != 1):
        roots.append(1.0)

    # bounds of pieces: (t, kind), kind is "end", "asymptote" or "on_line"
    bounds = [(0.0, "end")]
    for t_root in roots:
        kind = "on_line" if abs(offset_x(t_root)) < tol_on_line else "asymptote"
        if kind == "on_line":
            add_point_branch(t_root)
        if t_root == 0:
            bounds[0] = (0.0, kind)
        else:
            bounds.append((t_root, kind))
    if bounds[-1][0] < 1:
        bounds.append((1.0, "end"))

    for (ta, kind_a), (tb, kind_b) in zip(bounds[:-1], bounds[1:]):
        if tb - ta < 1e-12:
            continue
        t_start, t_end = ta, tb
        if kind_a == "asymptote":
            t_start = clip_at_blank_radius(segment, ta, tb)
        if kind_b == "asymptote":
            t_end = clip_at_blank_radius(segment, tb, ta)
        if t_start is None or t_end is None:
            continue
        t_vals = np.linspace(t_start, t_end, n_points)
        phis = np.full(n_points, np.nan)
        points = np.full((n_points, 3), np.nan)
        for i, t in enumerate(t_vals):
            # on_line ends are 0/0, the point trajectory branch covers them
            if (i == 0 and kind_a == "on_line") or (
                i == n_points - 1 and kind_b == "on_line"
            ):
                continue
            Phi, SweptPoint = envelope_point(segment, t)
            # sanity check, outside the limits the cutter is out of the blank
            if angle_limits[0] <= Phi <= angle_limits[1]:
                phis[i] = Phi
                points[i] = SweptPoint
        branches.append(("segment", t_vals, phis, points))
    return branches


def gen_vertex_points(
    segment_in: Curve, segment_out: Curve, n_points, tol_smooth=1e-5, tol_on_line=1e-9
):
    """Trajectory of a sharp convex vertex between 2 segments.

    The normals of the vertex form a fan from the incoming to the outgoing normal,
    each normal in the fan gives a roll angle, the vertex cuts over that range.
    If the fan contains the rolling line direction, phi -> +-inf in the fan:
        - vertex off the rolling line: 2 branches, each clipped at r_blank
        - vertex on the rolling line: the vertex cuts over its whole range in the blank
    Returns a list of (phis, points), empty if the vertex is smooth, concave
    or never reaches the blank.
    tol_smooth: sin of the tangent angle below which the vertex counts as smooth,
        must stay above the numeric derivative noise (~1e-6).
    """
    Vertex = segment_in(1)
    tangent_in = segment_in.derivative(1, direction=-1)
    tangent_out = segment_out.derivative(0, direction=1)
    tangent_in = tangent_in / np.linalg.norm(tangent_in)
    tangent_out = tangent_out / np.linalg.norm(tangent_out)
    # material on the left: left turn = convex corner
    turn = tangent_in[0] * tangent_out[1] - tangent_in[1] * tangent_out[0]
    if turn < tol_smooth:
        return []
    point_range = point_phi_range(Vertex)
    if point_range is None:
        return []
    lo, hi = point_range

    normal_in, normal_out = normal_of(tangent_in), normal_of(tangent_out)
    phi_in = calc_roll_angle(Vertex, normal_in)
    phi_out = calc_roll_angle(Vertex, normal_out)
    offset = Vertex[0] + path_curve_offset

    # fan rotates CCW by less than 180deg, n_x changes sign iff it contains +-y
    if normal_in[0] * normal_out[0] > 0:
        ranges = [(phi_in, phi_out)]
    elif abs(offset) < tol_on_line:
        ranges = [(lo, hi)]
    else:
        # n_y of the fan normal parallel to the rolling line
        if abs(normal_in[0]) < 1e-12:
            ny_cross = normal_in[1]
        elif abs(normal_out[0]) < 1e-12:
            ny_cross = normal_out[1]
        else:
            # CCW from +x side reaches +y, from -x side reaches -y
            ny_cross = np.sign(normal_in[0])
        # near the crossing phi ~ offset * n_y / n_x / r_pitch
        ranges = []
        if abs(normal_in[0]) >= 1e-12:
            far = hi if offset * ny_cross * normal_in[0] > 0 else lo
            ranges.append((phi_in, far))
        if abs(normal_out[0]) >= 1e-12:
            far = hi if offset * ny_cross * normal_out[0] > 0 else lo
            ranges.append((far, phi_out))

    branches = []
    for phi_a, phi_b in ranges:
        phi_a, phi_b = np.clip(phi_a, lo, hi), np.clip(phi_b, lo, hi)
        if abs(phi_b - phi_a) < 1e-12:
            continue
        phis = np.linspace(phi_a, phi_b, n_points)
        points = np.array([transform_to_blank(Vertex, Phi) for Phi in phis])
        branches.append((phis, points))
    return branches


def fixed_point_eval(CutPoint, Phi):
    return Phi, transform_to_blank(CutPoint, Phi)


def gen_swept_points(cut_profile: CurveChain, n_points):
    """Generate candidate boundary points of the cut, grouped by source.

    Returns a list of dicts in profile order:
        source: ("segment", k), ("point", k) or ("vertex", k),
            vertex k is the end of segment k,
            point k is a point of segment k on the rolling line
        param: segment t (segment) or roll angle (point, vertex)
        phi: roll angle of each point
        points: (n_points, 3) array in the blank frame, nan where no contact
        eval: exact evaluator, eval(param) -> (phi, point)
    A source can have multiple branches.
    """
    segments = cut_profile.get_curves()
    n_vertex = len(segments) if cut_profile.is_closed else len(segments) - 1
    branches = []
    for k, segment in enumerate(segments):
        for kind, param, phis, points in gen_segment_points(segment, n_points):
            if kind == "segment":
                branch = {"param": param, "eval": partial(envelope_point, segment)}
            else:
                # param is the t of the point on the segment
                branch = {
                    "param": phis,
                    "eval": partial(fixed_point_eval, segment(param)),
                }
            branch.update(source=(kind, k), phi=phis, points=points)
            branches.append(branch)
        if k < n_vertex:
            vertex_branches = gen_vertex_points(
                segment, segments[(k + 1) % len(segments)], n_points
            )
            for phis, points in vertex_branches:
                branches.append(
                    {
                        "source": ("vertex", k),
                        "param": phis,
                        "phi": phis,
                        "points": points,
                        "eval": partial(fixed_point_eval, segment(1)),
                    }
                )
    return branches


# --------------------------------------------------------------------------------------
# Trimming
# Branches are split at every event where their validity could change: leaving the
# domain, cusps, intersections. Between events a piece is either entirely valid or
# entirely invalid, so each piece is classified once, at its midpoint.
# The domain is half a tooth space: within r_blank, between the gear tooth
# centerline (polar angle 0) and the tooth space centerline (polar angle psi_space).
# Assumes the cutter profile runs from the rack gap center (y=0) to the rack tooth
# center (y=-pitch/2).
# --------------------------------------------------------------------------------------

psi_space = -pitch / 2 / r_pitch
DOMAIN_NAMES = ("blank", "tooth_center", "space_center")


def gear_polar(points):
    """Radius and polar angle around the gear center, points in plot coords."""
    x = points[..., 0] + r_pitch
    y = points[..., 1]
    return np.hypot(x, y), np.arctan2(y, x)


def domain_functions(points):
    """Distance-like functions of the domain boundaries, all negative inside.

    Order follows DOMAIN_NAMES, angles are scaled by the radius to get lengths.
    """
    r, psi = gear_polar(points)
    return np.stack([r - r_blank, psi * r, (psi_space - psi) * r], axis=-1)


def eval_point(branch, param):
    return branch["eval"](param)[1]


def param_span(branch):
    params = branch["param"]
    return abs(params[-1] - params[0])


def find_domain_events(branch):
    """Parameters where the branch crosses a domain boundary."""
    params = branch["param"]
    values = domain_functions(branch["points"])
    events = []
    for k, name in enumerate(DOMAIN_NAMES):

        def func(param, k=k):
            return domain_functions(eval_point(branch, param))[k]

        v = values[:, k]
        for i in range(len(params) - 1):
            # nan samples are outside the angle limits, out of the blank as well
            if np.isfinite(v[i]) and np.isfinite(v[i + 1]) and v[i] * v[i + 1] < 0:
                events.append(
                    (brentq(func, params[i], params[i + 1], xtol=1e-14), name)
                )
    return events


def find_cusps(branch):
    """Parameters where the branch reverses direction (envelope folds back)."""
    params = branch["param"]
    chords = np.diff(branch["points"][:, :2], axis=0)
    h = 1e-7 * param_span(branch)
    events = []
    for i in range(len(chords) - 1):
        if not np.all(np.isfinite(chords[i : i + 2])):
            continue
        if np.dot(chords[i], chords[i + 1]) >= 0:
            continue
        ref = chords[i]

        def speed_along_ref(param, ref=ref):
            delta = eval_point(branch, param + h) - eval_point(branch, param - h)
            return np.dot(delta[:2], ref)

        a, b = params[i], params[i + 2]
        if speed_along_ref(a) * speed_along_ref(b) < 0:
            events.append((brentq(speed_along_ref, a, b, xtol=1e-14), "cusp"))
        else:
            events.append((params[i + 1], "cusp"))
    return events


def cross2(a, b):
    return a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]


def polyline_crossings(points_a, points_b, same=False):
    """Index pairs and local parameters of crossing polyline segments.

    nan points never cross (comparisons with nan are False).
    """
    a0, da = points_a[:-1, :2], np.diff(points_a[:, :2], axis=0)
    b0, db = points_b[:-1, :2], np.diff(points_b[:, :2], axis=0)
    w = b0[np.newaxis, :] - a0[:, np.newaxis]
    with np.errstate(divide="ignore", invalid="ignore"):
        denom = cross2(da[:, np.newaxis], db[np.newaxis, :])
        s = cross2(w, db[np.newaxis, :]) / denom
        u = cross2(w, da[:, np.newaxis]) / denom
    mask = (s >= 0) & (s < 1) & (u >= 0) & (u < 1)
    if same:
        i, j = np.indices(mask.shape)
        mask &= j > i + 1
    i, j = np.nonzero(mask)
    return zip(i, j, s[i, j], u[i, j])


def refine_intersection(branch_a, branch_b, param_a, param_b, tol=1e-9):
    def residual(x):
        return eval_point(branch_a, x[0])[:2] - eval_point(branch_b, x[1])[:2]

    sol = root(residual, [param_a, param_b], method="hybr", options={"xtol": 1e-14})
    # judged by the residual only, the solver reports no progress when it reaches
    # the noise floor of the envelope (~1e-10, numeric derivative of the cutter)
    if not np.all(np.isfinite(sol.x)) or np.linalg.norm(residual(sol.x)) > tol:
        return None
    return sol.x


def is_param_within(branch, param, tol=1e-9):
    params = branch["param"]
    lo, hi = min(params[0], params[-1]), max(params[0], params[-1])
    return lo - tol * param_span(branch) <= param <= hi + tol * param_span(branch)


def is_param_at_end(branch, param, tol=1e-9):
    params = branch["param"]
    span = param_span(branch)
    return min(abs(param - params[0]), abs(param - params[-1])) < tol * span


def find_intersections(branches):
    """Intersections between all branch pairs (and self-intersections).

    Returns per-branch event lists. Endpoints shared by both branches are
    joins created by construction (eg. segment end - vertex trajectory), not events.
    """
    events = [[] for _ in branches]
    for ia, branch_a in enumerate(branches):
        for ib in range(ia, len(branches)):
            branch_b = branches[ib]
            same = ia == ib
            crossings = polyline_crossings(
                branch_a["points"], branch_b["points"], same=same
            )
            for i, j, s, u in crossings:
                pa, pb = branch_a["param"], branch_b["param"]
                guess_a = pa[i] + s * (pa[i + 1] - pa[i])
                guess_b = pb[j] + u * (pb[j + 1] - pb[j])
                sol = refine_intersection(branch_a, branch_b, guess_a, guess_b)
                if sol is None:
                    continue
                param_a, param_b = sol
                if not (
                    is_param_within(branch_a, param_a)
                    and is_param_within(branch_b, param_b)
                ):
                    continue
                if is_param_at_end(branch_a, param_a) and is_param_at_end(
                    branch_b, param_b
                ):
                    continue
                if same and abs(param_a - param_b) < 1e-9 * param_span(branch_a):
                    continue
                events[ia].append((param_a, "intersection"))
                events[ib].append((param_b, "intersection"))
    return events


def make_rack_polygon(cut_profile: CurveChain, n_samples=2000):
    """Closed polygon of the full rack: profile mirrored to a full tooth, repeated.

    The rack is long enough that its closed-off ends never reach the blank
    within the angle limits.
    Line segments only need their end points, other curves are sampled densely.
    """
    travel = r_pitch * np.max(np.abs(angle_limits))
    n_teeth_side = int(np.ceil((travel + r_blank) / pitch)) + 1
    half = []
    for segment in cut_profile.get_curves():
        n = 2 if isinstance(segment, LineCurve) else n_samples
        # ArcCurve does not take array input
        half.append(np.array([segment(t) for t in np.linspace(0, 1, n)])[:, :2])
    half = np.concatenate(half)
    # mirror on the rack tooth center: y -> -pitch - y
    mirrored = half[::-1] * np.array([1, -1]) + np.array([0, -pitch])
    tooth = np.vstack([half, mirrored])
    outline = np.vstack(
        [tooth + np.array([0, -k * pitch]) for k in range(-n_teeth_side, n_teeth_side)]
    )
    x_body = np.max(outline[:, 0]) + 2 * module
    outline = np.vstack(
        [outline, [x_body, outline[-1, 1]], [x_body, outline[0, 1]]]
    )
    # duplicate points at the mirror and tooth joins turn into self-intersections
    # if the polygon is transformed
    keep = np.append(True, np.linalg.norm(np.diff(outline, axis=0), axis=1) > 1e-12)
    polygon = Polygon(outline[keep])
    shapely.prepare(polygon)
    return polygon


def penetration_depth(point, rack_polygon, n_phi=4000):
    """Max depth a blank point gets inside the rack over the angle limits.

    0 for points on the boundary of the cut, positive for points removed.
    """
    blank_point = point[:2] + np.array([r_pitch, 0])
    rack_boundary = rack_polygon.exterior

    def signed_depth(phis):
        c, s = np.cos(phis), np.sin(phis)
        x = c * blank_point[0] - s * blank_point[1] - r_pitch - path_curve_offset
        y = s * blank_point[0] + c * blank_point[1] - r_pitch * phis
        inside = shapely.contains_xy(rack_polygon, x, y)
        depth = np.zeros_like(phis)
        if np.any(inside):
            depth[inside] = shapely.distance(
                rack_boundary, shapely.points(x[inside], y[inside])
            )
        return depth

    phis = np.linspace(*angle_limits, n_phi)
    depths = signed_depth(phis)
    k = np.argmax(depths)
    if depths[k] == 0:
        return 0.0
    result = minimize_scalar(
        lambda phi: -signed_depth(np.array([phi]))[0],
        bounds=(phis[max(k - 1, 0)], phis[min(k + 1, n_phi - 1)]),
        method="bounded",
        options={"xatol": 1e-12},
    )
    return max(depths[k], -result.fun)


def split_pieces(branch, events, end_types=("end", "end"), tol=1e-9):
    """Split a branch at its events, in the branch's own direction.

    Events closer than tol (fraction of the branch) to the ends or each other
    are merged, eg. a branch clipped at r_blank also has a blank event at its end.

    Returns a list of pieces: dict(branch, p0, p1, end_types).
    """
    params = branch["param"]
    finite = np.nonzero(np.isfinite(branch["points"][:, 0]))[0]
    if len(finite) < 2:
        return []
    p_start, p_end = params[finite[0]], params[finite[-1]]
    span = p_end - p_start
    cuts = [(p_start, end_types[0])]
    for param, kind in sorted(events, key=lambda e: (e[0] - p_start) / span):
        progress = (param - p_start) / span
        if tol < progress < 1 - tol and (param - cuts[-1][0]) / span > tol:
            cuts.append((param, kind))
    cuts.append((p_end, end_types[1]))
    return [
        {"branch": branch, "p0": pa, "p1": pb, "end_types": (ka, kb)}
        for (pa, ka), (pb, kb) in zip(cuts[:-1], cuts[1:])
    ]


def piece_end_point(piece, which):
    """End point of a piece, stepping slightly inside if the end has no point (0/0)."""
    p0, p1 = piece["p0"], piece["p1"]
    param, inward = (p0, p1) if which == 0 else (p1, p0)
    point = eval_point(piece["branch"], param)
    if not np.all(np.isfinite(point)):
        point = eval_point(piece["branch"], param + (inward - param) * 1e-9)
    return point


def classify_pieces(pieces, rack_polygon, tol_domain=1e-9, tol_depth=1e-5):
    """Mark pieces as 'outside' (domain), 'removed' or 'valid'."""
    for piece in pieces:
        mid = eval_point(piece["branch"], (piece["p0"] + piece["p1"]) / 2)
        if not np.all(np.isfinite(mid)) or np.any(domain_functions(mid) > tol_domain):
            piece["status"] = "outside"
            continue
        piece["depth"] = penetration_depth(mid, rack_polygon)
        piece["status"] = "valid" if piece["depth"] < tol_depth * module else "removed"


def chain_pieces(valid_pieces, tol=1e-7):
    """Order the valid pieces from the tip to the root by connecting end points.

    Starts at the end point on the blank circle (or tooth centerline) closest to
    the tooth center. Each piece gets a 'reversed' flag for its direction in the chain.
    """
    tol = tol * r_blank
    for piece in valid_pieces:
        piece["ends"] = (piece_end_point(piece, 0), piece_end_point(piece, 1))

    start, best_psi = None, -np.inf
    for piece in valid_pieces:
        for which, point in enumerate(piece["ends"]):
            f_blank, f_tooth, _ = domain_functions(point)
            psi = gear_polar(point)[1]
            on_start_boundary = abs(f_blank) < tol or abs(f_tooth) < tol
            if on_start_boundary and psi > best_psi:
                start, best_psi = (piece, which), psi
    if start is None:
        raise RuntimeError("no valid piece reaches the blank circle")

    chain = []
    unused = list(valid_pieces)
    piece, which = start
    while True:
        unused.remove(piece)
        piece["reversed"] = which == 1
        chain.append(piece)
        current = piece["ends"][1 - which]
        candidates = [
            (other, w)
            for other in unused
            for w in (0, 1)
            if np.linalg.norm(other["ends"][w] - current) < tol
        ]
        if not candidates:
            break
        if len(candidates) > 1:
            print("warning: ambiguous chain at", current[:2], len(candidates))
        piece, which = candidates[0]

    end_distance = abs(domain_functions(current)[2])
    if end_distance > tol:
        print(f"warning: chain ends {end_distance:.2e} off the tooth space centerline")
    if unused:
        print(f"warning: {len(unused)} valid pieces not in the chain")
    return chain


def piece_to_curve(piece):
    """Exact Curve of a piece in the gear frame (gear center at origin, tooth on +x)."""
    p_from, p_to = piece["p0"], piece["p1"]
    if piece["reversed"]:
        p_from, p_to = p_to, p_from
    branch = piece["branch"]

    def func(u):
        return eval_point(branch, p_from + (p_to - p_from) * u) + RIGHT * r_pitch

    return Curve(func)


def is_root_arc(piece, tol=1e-9, n=11):
    """True if the piece lies on a circle around the gear center."""
    params = np.linspace(piece["p0"], piece["p1"], n)
    radii = gear_polar(np.array([eval_point(piece["branch"], p) for p in params]))[0]
    return np.ptp(radii) < tol * r_blank


def trim_cut(branches, cut_profile, n_piece_samples=300):
    """Trim candidate branches to the tooth flank and the addendum-dedendum arcs.

    Returns dict:
        flank: CurveChain from the tip to the root, gear frame
        addendum: ArcCurve on r_blank from the tooth centerline to the flank, or None
        dedendum: ArcCurve from the flank to the tooth space centerline, or None
        chain: the valid pieces in order (debug info: source, params, end types)
        pieces: all pieces with their status
    """
    # branches on the blank circle (eg. from the rack root line) would intersect
    # the blank circle everywhere, the blank circle itself takes their place
    candidates = []
    for branch in branches:
        radius_error = domain_functions(branch["points"])[:, 0]
        radius_error = radius_error[np.isfinite(radius_error)]
        if len(radius_error) > 1 and np.all(np.abs(radius_error) < 1e-9 * r_blank):
            continue
        candidates.append(branch)

    # split at domain events and cusps first, the pieces are smooth and resampled
    # over their exact ranges, so their polylines reach right into cusp tips,
    # where intersections hide within the chord error of the branch samples
    smooth_pieces = []
    for branch in candidates:
        events = find_domain_events(branch) + find_cusps(branch)
        for piece in split_pieces(branch, events):
            params = np.linspace(piece["p0"], piece["p1"], n_piece_samples)
            points = np.array([eval_point(branch, p) for p in params])
            smooth_pieces.append(
                dict(branch, param=params, points=points, end_types=piece["end_types"])
            )

    intersection_events = find_intersections(smooth_pieces)
    pieces = []
    for smooth_piece, events in zip(smooth_pieces, intersection_events):
        pieces += split_pieces(smooth_piece, events, smooth_piece["end_types"])

    rack_polygon = make_rack_polygon(cut_profile)
    classify_pieces(pieces, rack_polygon)
    chain = chain_pieces([p for p in pieces if p["status"] == "valid"])

    flank_pieces = list(chain)
    dedendum = None
    root_pieces = []
    while flank_pieces and is_root_arc(flank_pieces[-1]):
        root_pieces.insert(0, flank_pieces.pop())
    if root_pieces:
        first = root_pieces[0]
        p_root_start = first["ends"][1 if first["reversed"] else 0]
        last = root_pieces[-1]
        p_root_end = last["ends"][0 if last["reversed"] else 1]
        dedendum = ArcCurve.from_2_point_center(
            p0=p_root_start + RIGHT * r_pitch,
            p1=p_root_end + RIGHT * r_pitch,
            center=ORIGIN,
        )

    first = chain[0]
    p_tip = first["ends"][1 if first["reversed"] else 0]
    addendum = None
    if abs(domain_functions(p_tip)[1]) > 1e-7 * r_blank:
        addendum = ArcCurve.from_2_point_center(
            p0=RIGHT * r_blank, p1=p_tip + RIGHT * r_pitch, center=ORIGIN
        )

    return {
        "flank": CurveChain(*[piece_to_curve(p) for p in flank_pieces]),
        "addendum": addendum,
        "dedendum": dedendum,
        "chain": chain,
        "pieces": pieces,
    }


def print_chain(chain):
    for piece in chain:
        end_types = piece["end_types"][::-1] if piece["reversed"] else piece["end_types"]
        print(
            f"  {str(piece['branch']['source']):16s}"
            f" param [{piece['p0']:.6f}, {piece['p1']:.6f}]"
            f" ends {end_types[0]} -> {end_types[1]}"
            f" depth {piece['depth']:.1e}"
        )


def plot_result(result, branches, cut_profile):
    for branch in branches:
        points = branch["points"]
        plt.plot(points[:, 0], points[:, 1], "-", color="0.8", linewidth=0.8)
    for piece in result["pieces"]:
        if piece["status"] == "removed":
            params = np.linspace(piece["p0"], piece["p1"], 100)
            points = np.array([eval_point(piece["branch"], p) for p in params])
            plt.plot(points[:, 0], points[:, 1], ":", color="tab:red", linewidth=0.8)
    for piece in result["chain"]:
        params = np.linspace(piece["p0"], piece["p1"], 200)
        points = np.array([eval_point(piece["branch"], p) for p in params])
        plt.plot(points[:, 0], points[:, 1], "-", linewidth=2,
                 label=str(piece["branch"]["source"]))
        for point, kind in zip(piece["ends"], piece["end_types"]):
            marker = {"intersection": "rx", "cusp": "m^"}.get(kind, "k.")
            plt.plot(point[0], point[1], marker)
    for name in ("addendum", "dedendum"):
        arc = result[name]
        if arc is not None:
            points = np.array([arc(u) for u in np.linspace(0, 1, 100)])
            points = points - RIGHT * r_pitch
            plt.plot(points[:, 0], points[:, 1], "k-", linewidth=2.5)
    points = cut_profile(np.linspace(0, 1, 500))
    plt.plot(points[:, 0], points[:, 1], "k:", label="cutter")
    # zoom to the cutter, the blank is much larger
    plt.gca().set_aspect("equal")
    plt.xlim(np.min(points[:, 0]) - module, np.max(points[:, 0]) + module)
    plt.ylim(np.min(points[:, 1]) - module, np.max(points[:, 1]) + module)
    plt.legend()


use_wavy_cutter = True
if use_wavy_cutter:
    cut_profile = make_wavy_cut_profile(
        addendum_coefficient, dedendum_coefficient, pitch, pressure_angle, module
    )
else:
    cut_profile = make_cut_profile(
        addendum_coefficient, dedendum_coefficient, pitch, pressure_angle, module
    )

angle_limits = calc_angle_limits(cut_profile)
print("angle limits:", angle_limits)
branches = gen_swept_points(cut_profile, n_points=400)
result = trim_cut(branches, cut_profile)
print("flank chain, tip to root:")
print_chain(result["chain"])

plot_result(result, branches, cut_profile)
plt.show()
