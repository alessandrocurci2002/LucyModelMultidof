"""Reading and validation of the mechanical configuration (``config/gimbal.yaml``), conversion to SI.

This is the only place where the geometry of the gimbal enters the package and where the units of the input
file (mm, kg m^2, N m/rad) are converted to SI. Everything downstream receives an immutable, validated
:class:`GimbalConfig`.

Supported topology (the one of the thesis model): on the polar axis the telescope element has
``mode: telescope`` and the hub model is enabled (``polar.hub.enabled: true``); on the azimuth axis the
telescope hangs from the load node through ``azimuth.telescope_attachment``. Both axes therefore carry the
telescope as a node P off the shaft chain, which is what :func:`lucy_sim.structure.lagrange.axis_chains`
requires.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal, Mapping

import numpy as np
import yaml

from .inertia import build_cad_frames, direction

MM = 1e-3
"""mm to m."""

AxisName = Literal["azimuth", "polar"]


class ConfigError(ValueError):
    """The configuration file does not follow the expected schema."""


# ---------------------------------------------------------------------------------------------- data (SI)


@dataclass(frozen=True)
class Material:
    name: str
    rho: float
    """Density, kg/m^3."""
    G: float
    """Shear modulus, Pa."""


@dataclass(frozen=True)
class Layer:
    """Coaxial annular layer of a segment (``d_in = 0`` for a solid cylinder)."""

    shape: Literal["solid", "hollow"]
    d_out: float
    d_in: float
    material: Material
    carries_torque: bool
    """``True`` if the layer carries torque (adds stiffness); otherwise it only adds inertia."""


@dataclass(frozen=True)
class Segment:
    """Constant-section stretch of an element; the segments of an element are in series."""

    length: float
    layers: tuple[Layer, ...]


@dataclass(frozen=True)
class ExtraInertia:
    """Concentrated inertia attached to an element (split half and half on its two nodes)."""

    name: str
    J: float
    m: float | None
    I_t: float | None


@dataclass(frozen=True)
class Element:
    """Shaft element ``id`` between node ``id`` and node ``id + 1``.

    ``geometry``: inertia and stiffness from the layers (``k_manual``, if given, replaces the stiffness only);
    ``manual``: ``J_manual`` and ``k_manual`` from CAD; ``telescope``: no shaft, the rigid telescope joins the
    two nodes (its inertia is in the ``payload`` section)."""

    id: int
    mode: Literal["geometry", "manual", "telescope"]
    length: float
    segments: tuple[Segment, ...]
    k_joint: float | None
    extra: tuple[ExtraInertia, ...]
    J_manual: float | None
    k_manual: float | None
    m_manual: float | None
    I_t_manual: float | None
    y_c_manual: float | None


@dataclass(frozen=True)
class NodeInertia:
    """Concentrated inertia at a node. With ``frame: cad_cg`` (Onshape export: tensor about the centre of mass,
    axes parallel to N) ``J = I_zz,G + m d^2`` about the azimuth axis."""

    node: int
    name: str
    J: float
    m: float | None
    I_t: float | None
    includes_polar_shafts: bool = False
    """``True`` if the CAD export of the azimuth load already contains the polar shafts (not the telescope)."""


@dataclass(frozen=True)
class HubSection:
    """Axisymmetric stretch of a side hub, or the telescope tube clamped in the ring (SI)."""

    d_out: float
    d_in: float
    length: float
    material: Material


@dataclass(frozen=True)
class HubRing:
    """Telescope support ring: its axis is the pointing axis, the side hubs are fixed at two opposite points."""

    d_out: float
    d_in: float
    width: float
    material: Material
    poisson: float
    stiffness_factor: float = 1.0


@dataclass(frozen=True)
class Hub:
    """Polar hub (two side hubs, tube + flange, and the support ring) and compliant telescope attachment.

    The two halves of the hub replace the elements ``side_elements`` (massless springs), the hub centre H takes
    the nodes of ``center_element`` (the telescope element) and the telescope is a node P off the chain, attached
    to H by ``k_b``."""

    side_elements: tuple[int, int]
    center_element: int
    m: float
    r_G_N: np.ndarray
    I_G_N: np.ndarray
    in_azimuth_load: bool
    tube: HubSection
    flange: HubSection
    ring: HubRing
    k_b: float
    telescope_section: HubSection | None = None
    J_polar: float = 0.0
    """Inertia of the hub about the polar axis (telescope section clamped in the ring included), kg m^2."""
    J_azimuth: float = 0.0
    """Inertia of the hub about z_hat at theta = 0 (same section included), kg m^2."""


@dataclass(frozen=True)
class AzimuthAttachment:
    """Telescope hanging from the azimuth load L through a spring ``k`` about z_hat."""

    k: float


@dataclass(frozen=True)
class AxisConfig:
    """Data of one axis; nodes and elements numbered from 1."""

    name: AxisName
    n_nodes: int
    encoder_node: int
    torque_node: int
    end_node: int
    """``load_node`` on the azimuth axis, ``end_node`` on the polar axis."""
    masters: tuple[int | str, ...]
    elements: tuple[Element, ...]
    node_inertias: tuple[NodeInertia, ...]
    bearing_b: tuple[tuple[int, float], ...]
    """``(element, b)`` viscous bearing friction, N m s/rad."""
    payload_element: int | None = None
    payload_xi: float | None = None
    y_node1_from_o2: float | None = None
    nodes_direction: float = 1.0
    hub: Hub | None = None
    attachment: AzimuthAttachment | None = None

    def element(self, id_: int) -> Element:
        return self.elements[id_ - 1]


@dataclass(frozen=True)
class PayloadConfig:
    """Rigid telescope: Onshape tensor about its centre of mass, axes parallel to N (``frame: cad_cg``)."""

    m: float
    I_G_N: np.ndarray
    """Tensor as written in the file, kg m^2."""
    r_G_N: np.ndarray
    products_sign: Literal["tensor", "positive"]


R_P2 = np.array([[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]])
"""Columns: I_hat, J_hat, K_hat in the internal pointing frame P (x_P pointing, y_P = j_hat, z_P = x_P x y_P up
at theta = 0)."""


@dataclass(frozen=True)
class Frames:
    """Reference frames (De Donatis convention), SI.

    N inertial (origin O_N, z_hat up along the azimuth axis, the Onshape assembly frame); body 1 (fork, origin o_1
    on the azimuth axis); body 2 (telescope, origin o_2 = centre of rotation of the polar shaft, K_hat pointing,
    J_hat polar axis, theta > 0 raises the pointing). The dynamics uses internally the pointing frame P of the
    2-DOF model (x_P = K_hat, y_P = j_hat = -J_hat, z_P = I_hat)."""

    azimuth_axis_xy: np.ndarray
    z_o1: float
    z_o2: float
    theta_cad: float
    R_N1_0: np.ndarray
    R_NP_0: np.ndarray
    """Attitude of the pointing frame P in the CAD pose."""

    @property
    def r_o2(self) -> np.ndarray:
        return np.array([self.azimuth_axis_xy[0], self.azimuth_axis_xy[1], self.z_o2])


@dataclass(frozen=True)
class GimbalConfig:
    """Mechanical model of the gimbal, SI units."""

    azimuth: AxisConfig
    polar: AxisConfig
    payload: PayloadConfig
    frames: Frames
    zeta: float
    """Damping ratio of the first flexible mode of each axis (stiffness-proportional damping)."""
    g: float
    source: Path | None = None

    def axis(self, name: AxisName) -> AxisConfig:
        return self.azimuth if name == "azimuth" else self.polar


# ---------------------------------------------------------------------------------------------- readers


def _require(d: Mapping[str, Any], key: str, ctx: str) -> Any:
    if not isinstance(d, Mapping) or key not in d or d[key] is None:
        raise ConfigError(f"{ctx}: missing required field '{key}'")
    return d[key]


def _positive(value: Any, name: str, ctx: str) -> float:
    try:
        v = float(value)
    except (TypeError, ValueError) as e:
        raise ConfigError(f"{ctx}: '{name}' = {value!r} is not a number") from e
    if not np.isfinite(v) or v <= 0.0:
        raise ConfigError(f"{ctx}: '{name}' = {v} must be positive")
    return v


def _optional(value: Any, name: str, ctx: str, *, scale: float = 1.0, positive: bool = True) -> float | None:
    if value is None:
        return None
    try:
        v = float(value)
    except (TypeError, ValueError) as e:
        raise ConfigError(f"{ctx}: '{name}' = {value!r} is not a number") from e
    if positive and v < 0.0:
        raise ConfigError(f"{ctx}: '{name}' = {v} cannot be negative")
    return v * scale


def _vector3(v: Any, name: str, ctx: str) -> np.ndarray:
    try:
        a = np.array([float(x) for x in v], dtype=float)
    except (TypeError, ValueError) as e:
        raise ConfigError(f"{ctx}: '{name}' must be a list of 3 numbers") from e
    if a.shape != (3,):
        raise ConfigError(f"{ctx}: '{name}' must have 3 components")
    return a


def _tensor(v: Any, ctx: str) -> np.ndarray:
    """Symmetric 3x3 tensor as written in the file (kg m^2)."""
    try:
        I = np.array([[float(x) for x in row] for row in v], dtype=float)
    except (TypeError, ValueError) as e:
        raise ConfigError(f"{ctx}: 'I_kgm2' must be a 3x3 matrix of numbers") from e
    if I.shape != (3, 3):
        raise ConfigError(f"{ctx}: 'I_kgm2' must be 3x3")
    scale = max(float(np.max(np.abs(I))), 1e-300)
    if np.max(np.abs(I - I.T)) > 1e-9 * scale:
        raise ConfigError(f"{ctx}: 'I_kgm2' is not symmetric")
    return 0.5 * (I + I.T)


def _materials(data: Mapping[str, Any]) -> dict[str, Material]:
    raw = _require(data, "materials", "config")
    out = {}
    for name, m in raw.items():
        ctx = f"materials.{name}"
        out[str(name)] = Material(str(name), _positive(_require(m, "rho", ctx), "rho", ctx),
                                  _positive(_require(m, "G", ctx), "G", ctx))
    return out


def _material(name: Any, materials: Mapping[str, Material], ctx: str) -> Material:
    if str(name) not in materials:
        raise ConfigError(f"{ctx}: material '{name}' not defined in materials ({sorted(materials)})")
    return materials[str(name)]


def _layer(d: Mapping[str, Any], materials: Mapping[str, Material], ctx: str) -> Layer:
    shape = str(_require(d, "shape", ctx)).strip().lower()
    d_out = _positive(_require(d, "d_out_mm", ctx), "d_out_mm", ctx) * MM
    d_in_raw = d.get("d_in_mm")
    if shape == "solid":
        if d_in_raw is not None and float(d_in_raw) != 0.0:
            raise ConfigError(f"{ctx}: solid cylinder with d_in_mm = {d_in_raw}; use shape: hollow")
        d_in = 0.0
    elif shape == "hollow":
        if d_in_raw is None:
            raise ConfigError(f"{ctx}: hollow cylinder without d_in_mm")
        d_in = float(d_in_raw) * MM
        if not 0.0 < d_in < d_out:
            raise ConfigError(f"{ctx}: hollow cylinder needs 0 < d_in < d_out")
    else:
        raise ConfigError(f"{ctx}: shape = '{shape}', allowed 'solid' and 'hollow'")
    return Layer(shape, d_out, d_in, _material(_require(d, "material", ctx), materials, ctx),  # type: ignore[arg-type]
                 bool(d.get("carries_torque", True)))


def _extra(items: Any, ctx: str) -> tuple[ExtraInertia, ...]:
    out = []
    for j, e in enumerate(items or []):
        c = f"{ctx}.extra_inertias[{j}]"
        out.append(ExtraInertia(str(e.get("name", f"extra {j + 1}")), float(_optional(_require(e, "J_kgm2", c), "J_kgm2", c) or 0.0),
                                _optional(e.get("m_kg"), "m_kg", c), _optional(e.get("I_t_kgm2"), "I_t_kgm2", c)))
    return tuple(out)


def _element(d: Mapping[str, Any], materials: Mapping[str, Material], ctx: str) -> Element:
    id_ = int(_require(d, "id", ctx))
    ctx = f"{ctx} (element {id_})"
    mode = str(d.get("mode", "geometry")).strip().lower()
    k_joint = _optional(d.get("k_joint_Nm_rad"), "k_joint_Nm_rad", ctx)
    if k_joint is not None and k_joint <= 0.0:
        raise ConfigError(f"{ctx}: k_joint_Nm_rad must be positive or null")
    extra = _extra(d.get("extra_inertias"), ctx)
    k_man = _optional(d.get("k_Nm_rad"), "k_Nm_rad", ctx)
    if k_man is not None and k_man <= 0.0:
        raise ConfigError(f"{ctx}: k_Nm_rad must be positive")
    if mode == "geometry":
        segments = []
        for s, seg in enumerate(_require(d, "segments", ctx)):
            cs = f"{ctx}.segments[{s}]"
            L = _positive(_require(seg, "length_mm", cs), "length_mm", cs) * MM
            layers = tuple(_layer(x, materials, f"{cs}.layers[{j}]") for j, x in enumerate(_require(seg, "layers", cs)))
            if k_man is None and not any(x.carries_torque for x in layers):
                raise ConfigError(f"{cs}: no layer with carries_torque: true and no k_Nm_rad: zero stiffness")
            segments.append(Segment(L, layers))
        length = float(sum(s.length for s in segments))
        if d.get("length_mm") is not None and abs(float(d["length_mm"]) * MM - length) > 1e-9:
            raise ConfigError(f"{ctx}: length_mm differs from the sum of the segments")
        return Element(id_, "geometry", length, tuple(segments), k_joint, extra, None, k_man, None, None, None)
    if mode == "manual":
        J = _positive(_require(d, "J_kgm2", ctx), "J_kgm2", ctx)
        if k_man is None:
            raise ConfigError(f"{ctx}: manual mode without k_Nm_rad")
        L = _positive(_require(d, "length_mm", ctx), "length_mm", ctx) * MM
        return Element(id_, "manual", L, (), k_joint, extra, J, k_man, _optional(d.get("m_kg"), "m_kg", ctx),
                       _optional(d.get("I_t_kgm2"), "I_t_kgm2", ctx), _optional(d.get("y_c_mm"), "y_c_mm", ctx, scale=MM))
    if mode == "telescope":
        if d.get("segments") or d.get("J_kgm2") is not None or k_man is not None:
            raise ConfigError(f"{ctx}: mode: telescope has no shaft (no segments, J_kgm2 or k_Nm_rad)")
        L = _positive(_require(d, "length_mm", ctx), "length_mm", ctx) * MM
        return Element(id_, "telescope", L, (), None, (), None, None, None, None, None)
    raise ConfigError(f"{ctx}: mode = '{mode}', allowed 'geometry', 'manual' and 'telescope'")


def _node_inertias(d: Mapping[str, Any], name: AxisName, n_nodes: int, end_node: int, axis_xy: np.ndarray,
                   ctx: str) -> tuple[NodeInertia, ...]:
    out = []
    for j, n in enumerate(d.get("node_inertias") or []):
        c = f"{ctx}.node_inertias[{j}]"
        node = int(_require(n, "node", c))
        if not 1 <= node <= n_nodes:
            raise ConfigError(f"{c}: node {node} outside 1..{n_nodes}")
        shafts = bool(n.get("includes_polar_shafts", False))
        if shafts and not (name == "azimuth" and node == end_node):
            raise ConfigError(f"{c}: includes_polar_shafts applies only to the azimuth load (node {end_node})")
        if str(n.get("frame", "axis")).strip() == "cad_cg":
            if name != "azimuth":
                raise ConfigError(f"{c}: frame: cad_cg is supported on the azimuth axis only")
            m = _positive(_require(n, "m_kg", c), "m_kg", c)
            I_G = _tensor(_require(n, "I_kgm2", c), c)
            r_G = _vector3(_require(n, "r_G_mm", c), "r_G_mm", c) * MM
            dist = float(np.hypot(r_G[0] - axis_xy[0], r_G[1] - axis_xy[1]))
            J_G = float(I_G[2, 2])  # diagonal moment: independent of the sign of the products
            if J_G <= 0.0:
                raise ConfigError(f"{c}: I_zz about the centre of mass is not positive")
            out.append(NodeInertia(node, str(n.get("name", f"node {node}")), J_G + m * dist**2, m, None, shafts))
            continue
        J_raw = n.get("J_kgm2")
        if J_raw is None:
            raise ConfigError(f"{c}: missing J_kgm2")
        out.append(NodeInertia(node, str(n.get("name", f"node {node}")), float(_optional(J_raw, "J_kgm2", c) or 0.0),
                               _optional(n.get("m_kg"), "m_kg", c), _optional(n.get("I_t_kgm2"), "I_t_kgm2", c), shafts))
    if name == "azimuth" and not any(x.node == end_node for x in out):
        raise ConfigError(f"{ctx}: the inertia of the load (fork) at node {end_node} is required")
    return tuple(out)


def _hub_section(d: Mapping[str, Any], materials: Mapping[str, Material], ctx: str, solid: bool = False) -> HubSection:
    d_out = _positive(_require(d, "d_out_mm", ctx), "d_out_mm", ctx) * MM
    d_in = 0.0 if solid else float(_require(d, "d_in_mm", ctx)) * MM
    if not 0.0 <= d_in < d_out:
        raise ConfigError(f"{ctx}: needs 0 <= d_in_mm < d_out_mm")
    return HubSection(d_out, d_in, _positive(_require(d, "length_mm", ctx), "length_mm", ctx) * MM,
                      _material(_require(d, "material", ctx), materials, ctx))


def _hub(d: Any, materials: Mapping[str, Material], elements: tuple[Element, ...], payload_element: int) -> Hub:
    ctx = "polar.hub"
    if not d or not bool(d.get("enabled", False)):
        raise ConfigError(f"{ctx}: the hub model is required (enabled: true)")
    side = tuple(int(e) for e in _require(d, "side_elements", ctx))
    center = int(_require(d, "center_element", ctx))
    if center != payload_element or elements[center - 1].mode != "telescope":
        raise ConfigError(f"{ctx}: center_element must be the telescope element (payload_element) in mode: telescope")
    if sorted(side) != [center - 1, center + 1]:
        raise ConfigError(f"{ctx}: side_elements must be the elements next to {center}")
    ci = f"{ctx}.inertia"
    inertia = _require(d, "inertia", ctx)
    if str(inertia.get("frame", "cad_cg")) != "cad_cg":
        raise ConfigError(f"{ci}: only frame: cad_cg (Onshape export) is supported")
    ring = _require(d, "ring", ctx)
    cr = f"{ctx}.ring"
    ring_ = HubRing(_positive(_require(ring, "d_out_mm", cr), "d_out_mm", cr) * MM, _positive(_require(ring, "d_in_mm", cr), "d_in_mm", cr) * MM,
                    _positive(_require(ring, "width_mm", cr), "width_mm", cr) * MM, _material(_require(ring, "material", cr), materials, cr),
                    float(ring.get("poisson", 0.33)), _positive(ring.get("stiffness_factor", 1.0), "stiffness_factor", cr))
    if ring_.d_in >= ring_.d_out:
        raise ConfigError(f"{cr}: needs d_in_mm < d_out_mm")
    sides = _require(d, "side", ctx)
    tr = d.get("telescope_on_ring")
    section = _hub_section(tr, materials, f"{ctx}.telescope_on_ring") if tr and bool(tr.get("enabled", False)) else None
    att = _require(d, "telescope_attachment", ctx)
    I_G = np.array(_require(inertia, "I_kgm2", ci), dtype=float)
    if I_G.shape != (3, 3) or not np.allclose(I_G, I_G.T, atol=1e-9):
        raise ConfigError(f"{ci}: I_kgm2 must be a symmetric 3x3 matrix")
    return Hub(side_elements=(side[0], side[1]), center_element=center, m=_positive(_require(inertia, "m_kg", ci), "m_kg", ci),
               r_G_N=np.array(_require(inertia, "r_G_mm", ci), dtype=float) * MM, I_G_N=I_G,
               in_azimuth_load=bool(inertia.get("in_azimuth_load", False)),
               tube=_hub_section(_require(sides, "tube", f"{ctx}.side"), materials, f"{ctx}.side.tube"),
               flange=_hub_section(_require(sides, "flange", f"{ctx}.side"), materials, f"{ctx}.side.flange", solid=True),
               ring=ring_, k_b=_positive(_require(att, "k_b_Nm_rad", f"{ctx}.telescope_attachment"), "k_b_Nm_rad", ctx),
               telescope_section=section)


def _axis(d: Mapping[str, Any], name: AxisName, materials: Mapping[str, Material], axis_xy: np.ndarray) -> AxisConfig:
    ctx = name
    n_nodes = int(_require(d, "n_nodes", ctx))
    enc, torque = int(_require(d, "encoder_node", ctx)), int(_require(d, "torque_node", ctx))
    end = int(_require(d, "load_node" if name == "azimuth" else "end_node", ctx))
    for label, n in (("encoder_node", enc), ("torque_node", torque), ("load_node/end_node", end)):
        if not 1 <= n <= n_nodes:
            raise ConfigError(f"{ctx}: {label} = {n} outside 1..{n_nodes}")
    elements = tuple(_element(e, materials, f"{ctx}.elements[{j}]") for j, e in enumerate(_require(d, "elements", ctx)))
    if [e.id for e in elements] != list(range(1, n_nodes)):
        raise ConfigError(f"{ctx}: element ids must be 1..{n_nodes - 1}, contiguous and sorted")
    nodes = _node_inertias(d, name, n_nodes, end, axis_xy, ctx)
    bearing = []
    for k, b in (d.get("bearing_b") or {}).items():
        if not 1 <= int(k) <= n_nodes - 1 or float(b) < 0.0:
            raise ConfigError(f"{ctx}.bearing_b: invalid entry {k}: {b}")
        bearing.append((int(k), float(b)))
    masters = []
    for m in _require(d, "masters", ctx):
        if isinstance(m, str) and m.strip().upper() == "P":
            if name != "polar":
                raise ConfigError(f"{ctx}: master node P exists only on the polar axis")
            masters.append("P")
        else:
            masters.append(int(m))
    required = [enc, torque, end] + (["P"] if name == "polar" else [])
    if not all(r in masters for r in required):
        raise ConfigError(f"{ctx}: encoder, motor, load/end node (and P on the polar axis) must be master nodes")
    if name == "azimuth":
        att = d.get("telescope_attachment")
        if not att or not bool(att.get("enabled", False)):
            raise ConfigError("azimuth.telescope_attachment: required (enabled: true)")
        attachment = AzimuthAttachment(_positive(_require(att, "k_Nm_rad", f"{ctx}.telescope_attachment"), "k_Nm_rad", ctx))
        return AxisConfig(name, n_nodes, enc, torque, end, tuple(masters), elements, nodes, tuple(bearing), attachment=attachment)
    pe = int(_require(d, "payload_element", ctx))
    if not 1 <= pe <= n_nodes - 1 or elements[pe - 1].mode != "telescope":
        raise ConfigError(f"{ctx}: payload_element must be the element in mode: telescope")
    if any(e.mode == "telescope" and e.id != pe for e in elements):
        raise ConfigError(f"{ctx}: only the telescope element can have mode: telescope")
    centre = d.get("payload_center_from_node3_mm")
    xi = 0.5 if centre is None else float(centre) * MM / elements[pe - 1].length
    direction_ = str(d.get("nodes_increase_along", "+J_hat")).strip()
    if direction_ not in ("+J_hat", "-J_hat"):
        raise ConfigError(f"{ctx}: nodes_increase_along must be '+J_hat' or '-J_hat'")
    return AxisConfig(name, n_nodes, enc, torque, end, tuple(masters), elements, nodes, tuple(bearing), payload_element=pe,
                      payload_xi=xi, y_node1_from_o2=_optional(d.get("y_node1_from_o2_mm"), "y_node1_from_o2_mm", ctx, scale=MM, positive=False),
                      nodes_direction=1.0 if direction_ == "+J_hat" else -1.0, hub=_hub(d.get("hub"), materials, elements, pe))


def _payload(pl: Mapping[str, Any]) -> PayloadConfig:
    ctx = "payload"
    if str(pl.get("frame", "cad_cg")).strip() != "cad_cg":
        raise ConfigError(f"{ctx}.frame: only cad_cg (Onshape export about the centre of mass) is supported")
    sign = str(pl.get("products_sign", "tensor")).strip().lower()
    if sign not in ("tensor", "positive"):
        raise ConfigError(f"{ctx}.products_sign = '{sign}': allowed 'tensor' and 'positive'")
    return PayloadConfig(_positive(_require(pl, "m_kg", ctx), "m_kg", ctx), _tensor(_require(pl, "I_kgm2", ctx), ctx),
                         _vector3(_require(pl, "r_G_mm", ctx), "r_G_mm", ctx) * MM, sign)  # type: ignore[arg-type]


def _frames(f: Mapping[str, Any]) -> Frames:
    ctx = "frames"
    xy = np.array([float(v) for v in f.get("azimuth_axis_xy_mm", [0.0, 0.0])]) * MM
    if xy.shape != (2,):
        raise ConfigError(f"{ctx}.azimuth_axis_xy_mm must have 2 components")
    theta_cad = float(np.deg2rad(float(f.get("theta_cad_deg", 0.0))))
    try:
        R1, RP = build_cad_frames(direction(f.get("pointing_cad", "-Y")), direction(f.get("polar_axis_cad", "+X")), theta_cad)
    except (AssertionError, ValueError) as e:
        raise ConfigError(f"{ctx}: pointing_cad and polar_axis_cad are not consistent: {e}") from e
    return Frames(xy, float(_require(f, "o1_z_mm", ctx)) * MM, float(_require(f, "o2_z_mm", ctx)) * MM, theta_cad, R1, RP)


def load_config(path: str | Path) -> GimbalConfig:
    """Reads and validates the mechanical configuration; derived hub inertias are computed here."""
    from .hub import hub_inertia, telescope_section_inertia

    path = Path(path)
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ConfigError(f"{path}: not a YAML mapping")
    materials = _materials(data)
    frames = _frames(_require(data, "frames", "config"))
    azimuth = _axis(_require(data, "azimuth", "config"), "azimuth", materials, frames.azimuth_axis_xy)
    polar = _axis(_require(data, "polar", "config"), "polar", materials, frames.azimuth_axis_xy)
    hub = polar.hub
    assert hub is not None
    J_pol, J_az = hub_inertia(hub, frames)
    section = telescope_section_inertia(hub)
    dJ = section[2] if section is not None else 0.0
    from dataclasses import replace
    polar = replace(polar, hub=replace(hub, J_polar=J_pol + dJ, J_azimuth=J_az + dJ))
    damping = data.get("damping", {}) or {}
    if str(damping.get("model", "stiffness_proportional")) != "stiffness_proportional":
        raise ConfigError("damping.model: only 'stiffness_proportional' is supported")
    zeta = float(damping.get("zeta_first_flexible", 0.01))
    if not 0.0 <= zeta < 1.0:
        raise ConfigError("damping.zeta_first_flexible outside [0, 1)")
    return GimbalConfig(azimuth, polar, _payload(_require(data, "payload", "config")), frames, zeta,
                        float(data.get("g", 9.80665)), path)
