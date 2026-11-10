r"""Tuning of the LQG-PI with the procedure of Gawronski (2008, par. 8.2), as in the thesis.

Internal model without bearings, :math:`Q = V`. Transfer functions and 1 degree step on the reduced (design) model; the
doublets of LUCY on the plant with Dahl bearings and the closed-loop poles with the bearings sliding and stuck are
constraints. Step 1 is the first try of the book (small weights). Step 2: :math:`q_p` is raised until the bandwidth
reaches the fundamental. Step 3: :math:`q_i` is raised until oscillations appear, within the critical gain (7.9). Step 4:
the weights of the modes with too high a peak are raised. Then fine tuning on the index (8.27) reduced to settling time
and bandwidth, :math:`q_p` raised again, and :math:`V = \rho Q`: :math:`T_{ry}` on the design model does not depend on
V, so V is judged on the Dahl plant too. The thresholds not given by the book are in :func:`Tuner.violations` and in
``config/tuning.yaml``; the low-authority constraint :math:`q_p/q_f \ge 10` (par. 8.1.4) is also there.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import yaml
from scipy.linalg import expm
from scipy.signal import butter, sosfiltfilt

from ..control import lqg
from ..model import GimbalModel
from ..paths import ROOT
from ..simulation.loops import run_lqg
from .doublets import Window, doublet_metrics, load_typical_commands, simulate_doublet


def db(G: np.ndarray) -> np.ndarray:
    return 20 * np.log10(np.maximum(np.abs(G), 1e-300))


def zoh(f: np.ndarray, T: float) -> np.ndarray:
    return np.sinc(f * T) * np.exp(-1j * np.pi * f * T)


def bandwidth(f: np.ndarray, T: np.ndarray) -> float:
    """First frequency at which |T| falls 3 dB below its low-frequency value."""
    k = np.flatnonzero(np.abs(T) < 10 ** (-3 / 20) * abs(T[0]))
    return float(f[k[0]]) if k.size else float(f[-1])


def frf(A, B, C, f) -> np.ndarray:
    return np.array([C @ np.linalg.solve(2j * np.pi * x * np.eye(len(A)) - A, B) for x in f])


def local_maximum(f, d):
    """Highest local maximum in the band, with its height above the higher of the two side minima; ``None`` if monotonic."""
    k = [i for i in range(1, len(d) - 1) if d[i] > d[i - 1] and d[i] >= d[i + 1]]
    if not k:
        return None
    i = max(k, key=lambda j: d[j])
    return {"f_Hz": float(f[i]), "dB": float(d[i]), "height_dB": float(d[i] - max(d[:i].min(), d[i:].min()))}


class _RecordEstimate:
    """Controller wrapper that keeps the position estimate."""

    def __init__(self, c):
        self.c, self.y_hat = c, 0.0

    def initial_state(self):
        return self.c.initial_state()

    def output(self, st, r):
        self.y_hat = float(st["x"][0])
        return self.c.output(st, r)

    def __getattr__(self, name):
        return getattr(self.c, name)

    def update(self, st, y, u):
        self.c.update(st, y, u)


def plain(o):
    """JSON-serialisable copy (numpy scalars to Python, keys to strings)."""
    if isinstance(o, dict):
        return {str(k): plain(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [plain(v) for v in o]
    if isinstance(o, (np.floating, np.integer, np.bool_)):
        return o.item()
    return o


class Tuner:
    """One axis: design model, typical doublet commands, full plants and the evaluation of a set of weights."""

    def __init__(self, model: GimbalModel, axis: str, settings: dict, commands: list[Window]):
        self.model, self.axis, self.s = model, axis, settings
        self.a = a = model.axes[axis]
        self.dt = model.dt_master
        self.n = len(a.Ap)
        fg = settings["frequency_grid_Hz"]
        self.F = np.geomspace(float(fg["start"]), float(fg["stop"]), int(fg["n"]))
        self.P = np.array([(a.Cp @ np.linalg.solve(2j * np.pi * x * np.eye(self.n) - a.Ap, a.Bp))[0, 0] for x in self.F])
        st = settings["step"]
        self.step_amplitude, self.step_start, self.step_duration = np.deg2rad(float(st["amplitude_deg"])), float(st["start_s"]), float(st["duration_s"])
        bands = [tuple(map(float, b)) for b in settings["resonance_bands_Hz"][axis]]
        self.modes = []
        for m in a.lightly_damped_modes:
            lo, hi = next((x, y) for x, y in bands if x <= m.f_hz <= y)
            self.modes.append({"f_Hz": m.f_hz, "zeta": m.zeta, "states": list(m.indices), "band": (lo, hi)})
        self.f_fund = self.modes[0]["f_Hz"]
        self.W = commands
        an = a.velocity_loop
        self.plants = {k: (A, an.B[:, 0], an.out_position) for k, A in a.bearing_variants().items()}
        lim = settings["limits"]
        self.I_lim, self.mod_min, self.peak_max = float(lim["current_A"]), float(lim["min_modulus"]), float(lim["peak_dB"])
        self.low_authority = float(lim["low_authority_ratio"])

    # ----------------------------------------------------------------------------- evaluation of one set of weights
    def controller(self, q_p, q_i, q_f=None, rho=1.0):
        return self.model.lqg_controller(self.axis, q_p, q_i, q_f, rho, base=float(self.s["base_flexible_weight"]))

    def loop(self, c) -> tuple[np.ndarray, np.ndarray]:
        """``(L, N)`` on the design model: :math:`L = F_y\\,ZOH\\,P`, :math:`N = F_r\\,ZOH\\,P`, :math:`T_{ry} = N/(1 + L)`."""
        z = np.exp(2j * np.pi * self.F * self.dt)
        Fr, Fy = c.transfer_functions(z, self.model.delay)
        h = zoh(self.F, self.dt)
        return Fy * h * self.P, Fr * h * self.P

    def evaluate(self, c) -> dict:
        L, N = self.loop(c)
        T = N / (1 + L)
        d = db(T)
        lam = np.log(c.regulator_poles.astype(complex)) / self.dt
        lam = lam[lam.imag > 0]
        modes = []
        for m in self.modes:
            b = (self.F >= m["band"][0]) & (self.F <= m["band"][1])
            x = lam[np.argmin(np.abs(np.abs(lam) / (2 * np.pi) - m["f_Hz"]))]
            modes.append({"f_Hz": m["f_Hz"], "band_max_dB": float(d[b].max()), "peak": local_maximum(self.F[b], d[b]), "zeta_regulator": float(-x.real / abs(x)),
                          "K_f_norm": float(np.linalg.norm(c.K_f[[j - 1 for j in m["states"]]])), "K_e_norm": float(np.linalg.norm(c.K_e[m["states"], 0]))})
        km = int(np.argmin(np.abs(1 + L)))
        a = self.a
        return {"k_p": c.k_p, "k_i": c.k_i, "bandwidth_Hz": bandwidth(self.F, T), "modulus": float(np.abs(1 + L[km])), "modulus_Hz": float(self.F[km]),
                "T_ry_max_dB": float(d.max()), "modes": modes, "K_e_position_per_s": float(c.K_e[0, 0] / self.dt),
                "largest_pole": lqg.worst_pole(lqg.closed_loop_poles(c, a.Ap, a.Bp[:, 0], a.Cp[0], self.model.delay), self.dt)}

    def step_metrics(self, t, y, fs) -> dict:
        G = self.step_amplitude
        out = np.flatnonzero(np.abs(y - G)[t >= self.step_start] > 0.02 * G)
        vib = sosfiltfilt(butter(4, (30.0, min(150.0, 0.45 * fs)), "bandpass", fs=fs, output="sos"), y) if np.abs(y).max() < 1e3 else np.full_like(y, np.nan)
        return {"overshoot_pct": float(100 * (y.max() / G - 1)), "t_set_2pct_s": float(out[-1] * (t[1] - t[0])) if out.size else 0.0,
                "vib_peak_urad": float(1e6 * np.nanmax(np.abs(vib))), "vib_after_100ms_urad": float(1e6 * np.nanmax(np.abs(vib[t >= self.step_start + 0.1])))}

    def step(self, c, duration: float | None = None) -> dict:
        """1 degree step on the design model at the master period, with the same delay as the simulators."""
        a, N = self.a, self.n
        E = expm(np.block([[a.Ap, a.Bp], [np.zeros((1, N + 1))]]) * self.dt)
        Phi, Gam, C = E[:N, :N], E[:N, N], a.Cp[0]
        t = np.arange(0.0, self.step_duration if duration is None else duration, self.dt)
        x, st = np.zeros(N), c.initial_state()
        y, u = np.empty(len(t)), np.empty(len(t))
        queue = [0.0] * self.model.delay
        for k, tk in enumerate(t):
            y[k] = C @ x
            if hasattr(c, "measure"):
                c.measure(st, y[k])
            queue.append(c.output(st, self.step_amplitude if tk >= self.step_start else 0.0))
            u_app = queue.pop(0)
            u[k] = u_app
            x = Phi @ x + Gam * u_app
            c.update(st, y[k], u_app)
        return {"t": t, "y": y, "u": u, "metrics": self.step_metrics(t, y, 1 / self.dt)}

    def oscillations(self, r: dict) -> int:
        """Extrema of the error outside the 2 % band after the first crossing of the command, below 20 Hz: the
        oscillations of the integral gain are slower than the fundamental (par. 8.1.7), the vibration of the modes is
        step 4."""
        if not np.all(np.isfinite(r["y"])) or np.abs(r["y"]).max() > 1e3:
            return 999
        y = sosfiltfilt(butter(4, 20.0, "lowpass", fs=1 / self.dt, output="sos"), r["y"])
        t, e = r["t"], (y - self.step_amplitude) / self.step_amplitude
        e = e[t >= self.step_start]
        i0 = np.flatnonzero(e >= 0)
        if not i0.size:
            return 0
        e = e[i0[0]:]
        d = np.diff(e)
        ext = np.flatnonzero(d[:-1] * d[1:] < 0) + 1
        return int(np.sum(np.abs(e[ext]) > 0.02))

    def doublets(self, c) -> dict:
        ms, peak = [], 0.0
        for w in self.W:
            r = simulate_doublet(self.model, self.axis, c, w)
            ms.append(doublet_metrics(w, r["p"], r["i"]))
            peak = max(peak, r["max_i_cmd"])
        return {"i_max_A": float(peak), "t_set_median_s": float(np.median([m["t_set_s"] for m in ms])), "t_set_max_s": float(max(m["t_set_s"] for m in ms)),
                "err_rms_median_mdeg": float(np.median([m["err_rms_mdeg"] for m in ms])), "over_max_mdeg": float(max(m["over_mdeg"] for m in ms))}

    def point(self, q_p, q_i, q_f=None, rho=1.0, with_doublets=True) -> dict:
        c = self.controller(q_p, q_i, q_f, rho)
        v = self.evaluate(c)
        r = self.step(c)
        out = {"q_p": float(q_p), "q_i": float(q_i), "q_f": {str(k): x for k, x in (q_f or {}).items()}, "rho": rho} | v | {
            "step": r["metrics"], "extrema_outside_2pct": self.oscillations(r), "k_i_over_critical": float(v["k_i"] / (0.25 * v["k_p"] ** 2)),
            "full_plants": {k: lqg.worst_pole(lqg.closed_loop_poles(c, *x, self.model.delay), self.dt) for k, x in self.plants.items()}}
        if with_doublets:
            out["doublets"] = self.doublets(c)
        return out

    def violations(self, x: dict, integral: bool = False) -> list[str]:
        """Violated conditions: bandwidth at the fundamental and current (book), min |1 + L| and stability on the full
        plants, low authority, plus the critical gain (7.9) and the oscillations for step 3."""
        v = []
        if x["bandwidth_Hz"] >= self.f_fund:
            v.append("bandwidth at the fundamental")
        if "doublets" in x and x["doublets"]["i_max_A"] > self.I_lim:
            v.append("current")
        if x["modulus"] < self.mod_min:
            v.append("min |1 + L|")
        if any(q["abs_z"] > 1 + 1e-9 for q in x["full_plants"].values()):
            v.append("unstable on a full plant")
        if integral and x["k_i_over_critical"] > 1:
            v.append("critical gain")
        if integral and x["extrema_outside_2pct"] >= 2:
            v.append("oscillations")
        if x.get("q_f") and x["q_p"] / max(x["q_f"].values()) < self.low_authority:
            v.append("low authority")
        return v

    # ----------------------------------------------------------------------------- steps of the procedure
    def raise_along(self, scale, fun, integral: bool) -> tuple[list, float, list[str]]:
        """Climbs ``scale``, then bisects between the last good point and the first that violates a limit."""
        points, good, bad, reason = [], None, None, []
        for q in scale:
            x = fun(q)
            x["limits"] = self.violations(x, integral)
            points.append(x)
            print(self.axis, "%.4g" % q, round(x["k_p"], 2), round(x["k_i"], 2), round(x["bandwidth_Hz"], 3), round(x["step"]["t_set_2pct_s"], 3),
                  round(x["doublets"]["i_max_A"], 1), round(x["modulus"], 3), x["limits"], flush=True)
            if x["limits"]:
                bad, reason = q, x["limits"]
                break
            good = q
        if good is None:
            return points, float(scale[0]), reason
        if bad is None:
            return points, float(good), []
        lo, hi = good, bad
        for _ in range(int(self.s["bisection_steps"])):
            m = float(np.sqrt(lo * hi))
            x = fun(m)
            x["limits"] = self.violations(x, integral)
            points.append(x)
            print(self.axis, " bis %.4g" % m, round(x["k_p"], 2), round(x["doublets"]["i_max_A"], 1), round(x["modulus"], 3), x["limits"], flush=True)
            if x["limits"]:
                hi, reason = m, x["limits"]
            else:
                lo = m
        return points, float(lo), reason

    def step_2(self):
        q_i0 = float(self.s["first_try"]["q_i"])
        return self.raise_along([float(x) for x in self.s["scales"]["q_p"]], lambda q: self.point(q, q_i0), False)

    def step_3(self, q_p, q_f=None):
        q_i0 = float(self.s["first_try"]["q_i"])
        kp = self.point(q_p, q_i0, q_f, with_doublets=False)["k_p"]
        q_c = (0.25 * kp ** 2) ** 2
        n = int(np.ceil(np.log10(max(q_c, 1.0) / q_i0) * 2)) + 2
        return self.raise_along([q_i0 * 10 ** (0.5 * k) for k in range(n)], lambda q: self.point(q_p, q, q_f), True)

    def step_4(self, q_p, q_i, q_f=None) -> tuple[list, dict]:
        """Raises the weight of the modes with the peak above the limit only, one at a time, until the peak falls or a
        limit stops it."""
        q_f = dict(q_f or {})
        history = []
        x = self.point(q_p, q_i, q_f)
        x["limits"] = self.violations(x)
        history.append(x)
        for im in range(len(self.modes)):
            if x["modes"][im]["band_max_dB"] < self.peak_max:
                continue
            for q in self.s["scales"]["q_f"]:
                trial = q_f | {im: float(q)}
                y = self.point(q_p, q_i, trial)
                y["limits"] = self.violations(y)
                y["raised_mode"] = im
                history.append(y)
                print(self.axis, "step 4 mode", im, "%.3g" % float(q), round(y["modes"][im]["band_max_dB"], 1), round(y["modulus"], 3), y["limits"], flush=True)
                if y["limits"]:
                    break
                q_f, x = trial, y
                if y["modes"][im]["band_max_dB"] < self.peak_max:
                    break
        return history, q_f

    @staticmethod
    def index(x: dict, ref: dict) -> float:
        """(8.27) reduced to settling time and bandwidth, each relative to the reference point."""
        return float(x["step"]["t_set_2pct_s"] / ref["step"]["t_set_2pct_s"] - x["bandwidth_Hz"] / ref["bandwidth_Hz"])

    def fine(self, ref, q_f) -> tuple[list, dict]:
        grid = []
        q_p0, q_i0, kp = ref["q_p"], ref["q_i"], ref["k_p"]
        q_c = (0.25 * kp ** 2) ** 2
        axis_i = sorted(set([float(x) for x in np.geomspace(q_i0, max(q_c, 10 * q_i0), 7)]))
        for fq in self.s["scales"]["fine_q_p"]:
            for qi in axis_i:
                x = self.point(q_p0 * float(fq), qi, q_f, with_doublets=False)
                lim = self.violations(x, True)
                if not lim:
                    x["doublets"] = self.doublets(self.controller(x["q_p"], qi, q_f))
                    lim = self.violations(x, True)
                x["limits"] = lim
                x["index"] = self.index(x, ref)
                grid.append(x)
        ok = [x for x in grid if not x["limits"]]
        return grid, min(ok, key=lambda x: x["index"]) if ok else ref

    def raise_again(self, q_p, q_i, q_f, rho=1.0) -> tuple[list, float | None, list[str]]:
        """From q_p up until a limit stops it; the peak of the modes counts too, because step 4 is done."""
        points, qs, reason = [], None, []
        for f in self.s["scales"]["q_p_raise"]:
            x = self.point(q_p * float(f), q_i, q_f, rho)
            x["limits"] = self.violations(x) + (["resonance peak"] if any(m["band_max_dB"] >= self.peak_max for m in x["modes"]) else [])
            points.append(x)
            if x["limits"]:
                reason = x["limits"]
                break
            qs = q_p * float(f)
        return points, qs, reason

    def try_V(self, q_p, q_i, q_f) -> list[dict]:
        """V = rho Q, at fixed q_p and raising it again for every rho."""
        out = []
        for rho in self.s["scales"]["rho"]:
            rho = float(rho)
            points, qs, reason = self.raise_again(q_p, q_i, q_f, rho)
            out.append({"rho": rho, "at_same_q_p": points[0], "q_p_max": qs, "stop": reason, "at_q_p_max": next((x for x in points if x["q_p"] == qs), None)})
            p0 = points[0]
            print(self.axis, "rho", rho, "base", round(p0["modulus"], 3), round(p0["doublets"]["i_max_A"], 1), round(p0["doublets"]["t_set_median_s"], 2),
                  round(p0["doublets"]["err_rms_median_mdeg"], 1), p0["limits"], "q_p max", qs, reason, flush=True)
        return out

    def choose_V(self, V: list[dict], ref: dict) -> tuple[float, float, list[dict]]:
        """V changes only if the index on the model improves and on the Dahl plant neither the settling time nor the
        rms error get worse."""
        imp = float(self.s["improvement"])
        res = []
        for q in V:
            x = q["at_q_p_max"]
            if q["rho"] == 1.0 or x is None:
                continue
            d, d0 = x["doublets"], ref["doublets"]
            res.append({"rho": q["rho"], "index": self.index(x, ref), "t_set_ratio": d["t_set_median_s"] / d0["t_set_median_s"],
                        "err_ratio": d["err_rms_median_mdeg"] / d0["err_rms_median_mdeg"],
                        "accepted": bool(self.index(x, ref) < -imp and d["t_set_median_s"] <= d0["t_set_median_s"] and d["err_rms_median_mdeg"] <= d0["err_rms_median_mdeg"])})
        ok = [e for e in res if e["accepted"]]
        if not ok:
            return 1.0, ref["q_p"], res
        b = min(ok, key=lambda e: e["index"])
        return b["rho"], next(q["q_p_max"] for q in V if q["rho"] == b["rho"]), res

    # ----------------------------------------------------------------------------- checks outside the procedure
    def cascade_reference(self) -> dict:
        """The current cascade on the full model: bandwidth, step and doublets for comparison."""
        cas = self.a.cascade
        A, B, C = cas.A, cas.B[:, 0], cas.out_position
        T = frf(A, B, C, self.F)
        N = len(A)
        E = expm(np.block([[A, B[:, None]], [np.zeros((1, N + 1))]]) * self.dt)
        t = np.arange(0.0, self.step_duration, self.dt)
        x, y = np.zeros(N), np.empty(len(t))
        for k, tk in enumerate(t):
            y[k] = C @ x
            x = E[:N, :N] @ x + E[:N, N] * (self.step_amplitude if tk >= self.step_start else 0.0)
        return {"bandwidth_Hz": bandwidth(self.F, T), "step": self.step_metrics(t, y, 1 / self.dt), "doublets": self.doublets(None), "K_pp": self.a.drive.p.Kpp}

    def dahl_check(self, c, duration: float = 120.0) -> dict:
        """Poles with the bearings sliding and stuck, and a 1 degree step on the Dahl plant (no current limit)."""
        t = self.model.time(duration)
        rec = _RecordEstimate(c)
        r = run_lqg(self.a, rec, np.where(t >= self.step_start, self.step_amplitude, 0.0), 0.0, self.dt, self.model.delay)
        y = r["p"]
        k = t <= self.step_duration
        return {"largest_pole": {n: lqg.worst_pole(lqg.closed_loop_poles(c, *v, self.model.delay), self.dt) for n, v in self.plants.items()},
                "y_100_120s_mdeg": [float(np.rad2deg(y[t >= 100].min()) * 1e3), float(np.rad2deg(y[t >= 100].max()) * 1e3)],
                "y_minus_y_hat_end_mdeg": float(np.rad2deg(y[-1] - rec.y_hat) * 1e3), "step_20s": self.step_metrics(t[k], y[k], 1 / self.model.dt),
                "i_max_A": float(np.abs(r["i"]).max())}

    # ----------------------------------------------------------------------------- whole procedure
    def run(self) -> dict:
        imp = float(self.s["improvement"])
        p2, qp, m2 = self.step_2()
        p3, qi, m3 = self.step_3(qp)
        st4, qf = self.step_4(qp, qi)
        after = self.point(qp, qi, qf)
        grid, best = self.fine(after, qf)
        second = bool(best is not after and best["index"] < -imp)
        qp2, qi2 = (best["q_p"], best["q_i"]) if second else (qp, qi)
        st4b, qf2 = self.step_4(qp2, qi2, qf) if second else ([], qf)
        pts_r, qp_r, stop_r = self.raise_again(qp2, qi2, qf2)
        qp_r = qp_r if qp_r is not None else qp2
        at_Q_V = self.point(qp_r, qi2, qf2)
        V = self.try_V(qp_r, qi2, qf2)
        rho, qp3, res_v = self.choose_V(V, at_Q_V)
        final = self.point(qp3, qi2, qf2, rho)
        final["limits"] = self.violations(final, True)
        check = self.dahl_check(self.controller(final["q_p"], final["q_i"], final["q_f"], final["rho"]))
        return {"fundamental_Hz": self.f_fund, "typical_commands": [w.name for w in self.W], "q_p": p2, "q_p_selected": qp, "q_p_stop": m2,
                "q_i": p3, "q_i_selected": qi, "q_i_stop": m3, "step_4": st4, "q_f_after_step_4": {str(k): v for k, v in qf.items()},
                "after_steps": after, "fine_grid": grid, "fine_best": best if best is not after else None, "second_iteration": second,
                "step_4_second": st4b, "q_f_final": {str(k): v for k, v in qf2.items()}, "at_Q_equal_V": at_Q_V, "q_p_raised": pts_r,
                "q_p_raised_selected": qp_r, "q_p_raised_stop": stop_r, "V": V, "V_decision": res_v, "rho_selected": rho, "final": final,
                "cascade": self.cascade_reference(), "dahl_check": check}


def load_settings(root: str | Path = ROOT, path: str = "config/tuning.yaml") -> dict:
    return yaml.safe_load((Path(root) / path).read_text(encoding="utf-8"))


def tune(model: GimbalModel, axes=("azimuth", "polar"), root: str | Path = ROOT, settings: dict | None = None) -> dict:
    """Runs the procedure on the axes; returns the summary per axis (``summary[axis]["final"]`` are the tuned weights)."""
    s = settings or load_settings(root)
    return {ax: plain(Tuner(model, ax, s, load_typical_commands(Path(root) / "data" / "doublets" / f"{ax}.npz")).run()) for ax in axes}


def weights(summary: dict) -> dict:
    """Weights of ``config/controllers.yaml`` from the summary of :func:`tune`."""
    return {ax: {"q_p": r["final"]["q_p"], "q_i": r["final"]["q_i"], "q_f": {int(k): v for k, v in r["final"]["q_f"].items()}, "rho": r["final"]["rho"]}
            for ax, r in summary.items()}
