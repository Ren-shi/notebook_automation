"""The interactive scene of the planner app: the solids of :mod:`physim.nuclear.scene` drawn with NiceGUI's 3D
scene (three.js), with selection and constrained dragging.

A detector is selected by a click and can then be dragged. The limits of :func:`physim.nuclear.scene.limits` hold
the drag in the browser; every position is checked again here, so no drag can end in the beam or inside another
solid. Needs NiceGUI (``pip install physim-engine[app]``).
"""

from __future__ import annotations

import math
from typing import Callable, Optional

import numpy as np

from . import scene as _scene

#: The scene's colours in each theme of the page.
PALETTES = {
    "light": {"ground": "#FFFFFF", "ink": "#16191C", "label": "rgba(255,255,255,0.85)", "active": "#4C78A8",
              "active2": "#7FA3CC", "board": "#3F7D4E", "crystal": "#B58B4C", "housing": "#8C959C",
              "foil": "#C9A227", "beam": "#D1495B", "selected": "#F28E2B", "element": "#E15759", "unsafe": "#B07AA1",
              "blocked": "#D62728", "wall": "#9AA3AA", "line": "#22272B"},
    "dark": {"ground": "#171B1E", "ink": "#E9EBE9", "label": "rgba(23,27,30,0.85)", "active": "#5B8FCC",
             "active2": "#86B0E0", "board": "#4E9A61", "crystal": "#C9A163", "housing": "#A9B2B9",
             "foil": "#E0BC45", "beam": "#F0707F", "selected": "#FFA94D", "element": "#FF7B7D", "unsafe": "#C99BD0",
             "blocked": "#FF5A52", "wall": "#6F7A82", "line": "#DADDDB"},
}

HEIGHT = 540


def _name(key: str, element: Optional[int] = None) -> str:
    return key if element is None else f"{key}|{element}"


def parse_name(name: Optional[str]) -> tuple:
    """(key, element) from the name of a scene object; (None, None) for anything that cannot be selected."""
    if not name or not name.split("|")[0].split(":")[0] in ("detector", "gamma", "target"):
        return None, None
    key, _, element = name.partition("|")
    return key, int(element) if element else None


def dimensions(solid) -> str:
    """A solid's size in a line, for its label in the scene."""
    hull = solid.hull()
    if solid.kind == "detector":
        p = hull[0]
        size = (f"{p.w:g} × {p.h:g} mm" if p.kind == "box" else
                f"⌀ {2 * p.r_in:g}–{2 * p.r_out:g} mm" if p.r_in else f"⌀ {2 * p.r_out:g} mm")
        return f"{size}, {p.depth * 1e3:g} µm thick"
    if solid.kind == "gamma":
        crystals = [p for p in solid.parts if p.role == "crystal"]
        c = crystals[0]
        text = f"{len(crystals)} × " if len(crystals) > 1 else ""
        return f"{text}⌀ {2 * c.r_out:g} × {c.depth:g} mm"
    return f"⌀ {2 * hull[0].r_out:g} mm (as drawn)"


class SceneView:
    """The scene of one page. ``planner`` returns the page's planner; ``on_select(key, element)`` is called when
    the selection changes, ``on_live(values)`` while a detector is dragged (see
    :meth:`physim.nuclear.planner.Planner.live`) and ``on_moved(key, position, mode)`` when it is dropped at an
    allowed place."""

    def __init__(self, planner: Callable, theme: str = "light", on_select: Optional[Callable] = None,
                 on_live: Optional[Callable] = None, on_moved: Optional[Callable] = None):
        from nicegui import ui

        self._ui = ui
        self.planner = planner
        self.theme = theme
        self.colours = PALETTES[theme]
        self.on_select, self.on_live, self.on_moved = on_select, on_live, on_moved
        self.selected: Optional[str] = None
        self.element: Optional[int] = None
        self.mode = "angle"
        self.limits: Optional[_scene.Limits] = None
        self._groups: dict = {}
        self._meshes: dict = {}
        self._solids: dict = {}
        self._notes: list = []
        self._good: Optional[np.ndarray] = None
        self._tinted = False
        self.drawn = None
        #: True while a drop is being handled (the page then leaves the redrawing to the view).
        self.busy = False
        self.extent = 100.0
        self.scene = ui.scene(
            width=800, height=HEIGHT, grid=False, background_color=self.colours["ground"],
            camera=ui.scene.perspective_camera(fov=40, near=1.0, far=1e5), on_click=self._clicked,
            click_events=["click"], on_drag_end=self._dropped).classes("w-full").style(f"height: {HEIGHT}px")
        self.scene.on("drag", self._dragging, throttle=0.12)
        self.scene.on("init", self._clicks_only)
        self.draw()
        self.look("default", duration=0)

    def _clicks_only(self) -> None:
        """Turning the view ends in a click too; only a press and release at one place should select."""
        self._ui.run_javascript(f"""
            const el = getElement({self.scene.id}).$el;
            let down = null;
            el.addEventListener('pointerdown', (e) => {{ down = [e.clientX, e.clientY]; }}, true);
            el.addEventListener('click', (e) => {{
                if (down && Math.hypot(e.clientX - down[0], e.clientY - down[1]) > 4) e.stopImmediatePropagation();
            }}, true);
        """)

    # -- drawing ------------------------------------------------------------------------------------------------

    @property
    def alive(self) -> bool:
        return not self.scene.is_deleted

    def draw(self) -> None:
        """Draw the experiment as it is now (the camera stays where it is)."""
        exp = self.planner().experiment
        #: The experiment the scene shows.
        self.drawn = exp
        c = self.colours
        solids = _scene.solids(exp)
        self._solids = {s.key: s for s in solids}
        self.extent = _scene.extent_mm(exp, solids)
        try:
            unsafe = self.planner().safety()
        except Exception:  # noqa: BLE001 -- the scene must draw even where the rates cannot be computed
            unsafe = {}
        if self.selected not in self._solids:
            self.selected, self.element, self.limits = None, None, None
        self.scene.clear()
        self._groups, self._meshes, self._notes = {}, {}, []
        e = self.extent
        with self.scene as sc:
            r = _scene.beam_radius_mm(exp)
            sc.cylinder(r, r, e, 16).rotate(math.pi / 2, 0, 0).move(0, 0, -e / 2).material(c["beam"], 0.9)
            sc.cylinder(r, r, e, 16).rotate(math.pi / 2, 0, 0).move(0, 0, e / 2).material(c["beam"], 0.25)
            style = f"color: {c['ink']}; font-size: 11px; background: {c['label']}; padding: 0 3px; border-radius: 3px"
            sc.text("beam", style).move(0, 0.03 * e, -0.95 * e)
            step = 10 if e < 80 else 50 if e < 400 else 100
            for k in range(-int(e // step), int(e // step) + 1):
                if k:
                    z = k * step
                    sc.line([0, -0.012 * e, z], [0, 0.012 * e, z]).material(c["line"])
                    sc.text(f"{z:g}", style + "; opacity: 0.75").move(0, -0.04 * e, z)
            if exp.chamber is not None:
                sc.sphere(_scene._mm(exp.chamber.radius), 32, 16, wireframe=True).material(c["wall"], 0.25)
            for s in solids:
                self._draw_solid(sc, s, unsafe.get(s.name, ()))
        self._paint()
        self._annotate()
        self._limit()

    def _draw_solid(self, sc, s, unsafe) -> None:
        c = self.colours
        with sc.group() as group:
            group.with_name(s.key)
            elements = any(p.element is not None and p.role == "active" for p in s.parts)
            for p in s.parts:
                if p.kind == "line":
                    sc.line(p.points[0], p.points[1]).material(c["line"]).with_name(s.key)
                    continue
                if p.role == "active" and p.hull and elements:
                    continue  # its rings or strips are drawn instead
                colour = c[p.role] if p.role in c else c["active"]
                if p.role == "active" and p.element is not None:
                    colour = c["unsafe"] if p.element in unsafe else c["active2" if p.element % 2 else "active"]
                opacity = 0.22 if p.role == "housing" else 1.0
                name = _name(s.key, p.element)
                mid = p.z0 - p.depth / 2
                meshes = []
                if p.kind == "box":
                    meshes.append(sc.box(p.w, p.h, p.depth).move(p.cx, p.cy, mid))
                elif p.r_in <= 0:
                    meshes.append(sc.cylinder(p.r_out, p.r_out, p.depth, 40).rotate(math.pi / 2, 0, 0).move(
                        p.cx, p.cy, mid))
                else:
                    for z in ((p.z0, p.z0 - p.depth) if p.hull else (p.z0,)):
                        meshes.append(sc.ring(p.r_in, p.r_out, 64).move(p.cx, p.cy, z))
                for m in meshes:
                    m.material(colour, opacity, side="both").with_name(name)
                    self._meshes.setdefault(s.key, []).append((m, p, colour, opacity))
        group.move(*s.centre).rotate_R(s.axes.tolist())
        self._groups[s.key] = group

    def _paint(self) -> None:
        """Colour the selected detector, and the selected ring, strip or crystal within it."""
        c = self.colours
        for key, meshes in self._meshes.items():
            for m, p, colour, opacity in meshes:
                if key == self.selected and p.role in ("active", "crystal"):
                    colour = c["element"] if (self.element is not None and p.element == self.element) else (
                        c["blocked"] if self._tinted else c["selected"])
                m.material(colour, opacity, side="both")
        for key, group in self._groups.items():
            group.draggable(key == self.selected and key != "target")

    def _annotate(self) -> None:
        """The selected solid's name, place and size, written beside it, with a line from the target."""
        for obj in self._notes:
            obj.delete()
        self._notes = []
        s = self._solids.get(self.selected)
        if s is None:
            return
        c = self.colours
        style = (f"color: {c['ink']}; font-size: 12px; background: {c['label']}; padding: 2px 6px; "
                 "border-radius: 4px; white-space: nowrap; margin-top: -30px")
        with self.scene as sc:
            if s.kind == "target":
                text = f"{s.name} · {dimensions(s)}"
            else:
                d = float(np.linalg.norm(s.centre))
                theta = math.degrees(math.acos(max(-1.0, min(1.0, s.centre[2] / d))))
                text = f"{s.name} · {d:.1f} mm from the target · θ {theta:.1f}° · {dimensions(s)}"
                self._notes.append(sc.line([0, 0, 0], list(s.centre)).material(c["selected"]))
            self._notes.append(sc.text(text, style).move(*s.centre))

    def look(self, view: str = "default", duration: float = 0.4) -> None:
        """Point the camera: "default" (from upstream, above and to the side), "side", "top" or "beam" (along the
        beam, from upstream). The vertical axis (y) is up, except from the top."""
        e = 2.6 * self.extent
        x, y, z = {"default": (0.55 * e, 0.4 * e, -0.75 * e), "side": (e, 0, 0), "top": (0, e, 0),
                   "beam": (0, 0, -e)}[view]
        up = (1, 0, 0) if view == "top" else (0, 1, 0)
        self.scene.move_camera(x=x, y=y, z=z, look_at_x=0, look_at_y=0, look_at_z=0, up_x=up[0], up_y=up[1],
                               up_z=up[2], duration=duration)

    # -- selection and dragging ---------------------------------------------------------------------------------

    def select(self, key: Optional[str], element: Optional[int] = None, notify: bool = True) -> None:
        """Select a solid by key (``None`` for the whole experiment), and one of its rings, strips or crystals."""
        if key not in self._solids:
            key, element = None, None
        changed = (key, element) != (self.selected, self.element)
        self.selected, self.element, self._tinted = key, element, False
        self._limit()
        self._paint()
        self._annotate()
        if changed and notify and self.on_select:
            self.on_select(key, element)

    def set_mode(self, mode: str) -> None:
        """What a drag changes: the "angle" (at a fixed distance) or the "distance" (along a fixed direction)."""
        self.mode = mode
        self._limit()

    def _limit(self) -> None:
        key = self.selected
        if key is None or key == "target":
            self.limits = None
            return
        exp = self.planner().experiment
        self.limits = _scene.limits(exp, key, self.mode)
        self._good = _scene.position_of(exp, key)
        self.scene._props["drag-constraints"] = self.limits.constraints()
        self.scene.update()

    def _clicked(self, e) -> None:
        for hit in e.hits:
            key, element = parse_name(hit.object_name)
            if key is not None:
                self.select(key, element)
                return
        self.select(None)

    def _trial(self, x: float, y: float, z: float) -> tuple:
        pos = self.limits.apply((x, y, z))
        fields = _scene.move_fields(self.planner().experiment, self.selected, pos, self.limits.mode)
        return pos, fields

    def _dragging(self, e) -> None:
        key, _ = parse_name(e.args.get("object_name"))
        if key is None or key != self.selected or self.limits is None:
            return
        pos, fields = self._trial(e.args["x"], e.args["y"], e.args["z"])
        live = self.planner().live(key, fields)
        blocked = bool(live.get("problem"))
        if not blocked:
            self._good = pos
            if self.limits.mode == "angle":
                # Keep it facing the target as it moves over the sphere.
                trial = _scene.moved(self.planner().experiment, key, fields)
                turned = next(s for s in _scene.solids(trial) if s.key == key)
                self._groups[key].rotate_R(turned.axes.tolist())
        if blocked != self._tinted:
            self._tinted = blocked
            self._paint()
        if self.on_live:
            self.on_live(live)

    def _dropped(self, e) -> None:
        key, _ = parse_name(e.object_name)
        if key is None or key != self.selected or self.limits is None:
            self.draw()
            return
        pos, fields = self._trial(e.x, e.y, e.z)
        exp = self.planner().experiment
        why = _scene.blocked(exp, key, fields)
        if why:
            self._ui.notify(f"{why} It stays at the last place that is allowed.", type="warning")
            pos = self._good if self._good is not None else _scene.position_of(exp, key)
            fields = _scene.move_fields(exp, key, pos, self.limits.mode)
            if _scene.blocked(exp, key, fields):
                pos = _scene.position_of(exp, key)
        else:
            held = self.limits.stopped(pos)
            if held:
                self._ui.notify(f"Stopped here: {held}", type="info")
        self._tinted = False
        mode = self.limits.mode
        if self.on_moved and float(np.linalg.norm(pos - _scene.position_of(exp, key))) > 1e-6:
            self.busy = True
            try:
                self.on_moved(key, pos, mode)
            finally:
                self.busy = False
        elif self.on_select:
            self.on_select(key, self.element)  # nothing moved: the side panel goes back to the detector's numbers
        if self.alive:
            self.draw()


__all__ = ["HEIGHT", "PALETTES", "SceneView", "dimensions", "parse_name"]
