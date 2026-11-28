#!/usr/bin/env python3
# The cheap body of the agent fish (fish_agents.py) as a model to fit (plan 2, work package E): the speed through the
# water along the axis U, the turn rate and the slip across the axis, driven by the gait (f, a, b).
#
#   U relaxes to s_eff f a L r(t) over n_up tail beats (faster gait) or n_down tail beats (slower), with the tail
#   still (a < 0.5) it coasts with tau_coast (towards the gait as it fades), with the body bent (|b| > 0.5) while the
#   tail stills (a < a_brake: 0.5 as fish_agents.py, 0.95 for the fit: from the start of the brake of the brain) it
#   brakes, dU/dt = -C_brake U^2; r(t) = min(1, t/t_ramp), the ramp of the gait
#   with C_d > 0 (the drag law) the glide with the tail still follows dU/dt = C_d/L (u|u| - U|U|), u = s_eff f a L
#   r(t), instead of tau_coast (validation 30: the 2D swimmer glides that way, 1/U grows linearly in time)
#   s_eff = s (1 + s_f (f/f0 - 1)) (1 - s_b b^2)        the stride, its slope with the frequency, its loss in a bend
#   omega = k_eff b U / L, k_eff = k (1 - k_b |b|)        swimming; with the tail still k_brake b U / L
#   v_slip = beta_b b U                                  the velocity across the axis through the water
# With s_f = s_b = k_b = beta_b = 0 it is the body of fish_agents.py with n_up, n_down (body_step), the same numbers.
# The power of the deformation into the water per tail beat: P = c3 (f a)^3 + c1 (f a) U^2 (per the depth of a 2D
# body), fitted apart.
#
# Architect: Hans Bihs
import math
import numpy as np

PARAMS = ("s", "s_f", "s_b", "n_up", "n_down", "tau_coast", "C_d", "C_brake", "k", "k_b", "k_brake", "beta_b")
DEFAULT = {"L": 1.0, "f0": 1.0, "t_ramp": 1.0, "s": 0.5, "s_f": 0.0, "s_b": 0.0, "n_up": 0.8, "n_down": 0.65,
           "tau_coast": 11.0, "C_d": 0.0, "C_brake": 3.5, "k": 1.3, "k_b": 0.0, "k_brake": 1.77, "beta_b": 0.0, "a_brake": 0.5}


def stride(p, f, b):
    return p["s"] * (1.0 + p["s_f"] * (f / p["f0"] - 1.0)) * (1.0 - p["s_b"] * b * b)


def step(U, f, a, b, t, dt, p):
    """U after dt; the turn rate and the slip velocity with it (as fish_agents.py turns with the new speed)"""
    L = p["L"]
    if a < p["a_brake"] and abs(b) > 0.5:
        U1 = U - dt * p["C_brake"] * U * abs(U)
        return U1, p["k_brake"] * b * U1 / L, p["beta_b"] * b * U1
    tgt = stride(p, f, b) * f * a * L * min(1.0, t / p["t_ramp"])
    cd = p.get("C_d", 0.0)
    if cd > 0.0 and a < 0.5:
        # the drag law of the glide: with the tail still dU/dt = C_d/L (tgt|tgt| - U|U|), not past tgt
        U1 = U + dt * cd / L * (tgt * abs(tgt) - U * abs(U))
        U1 = max(U1, tgt) if tgt < U else min(U1, tgt)
    elif a < 0.5:
        U1 = U + dt / p["tau_coast"] * (tgt - U)
    else:
        n = p["n_up"] if tgt >= U else p["n_down"]
        U1 = U + dt * max(f, 0.1 * p["f0"]) / n * (tgt - U)
    if a < 0.5:
        return U1, p["k_brake"] * b * U1 / L, p["beta_b"] * b * U1
    kk = p["k"] * (1.0 - p["k_b"] * abs(b))
    return U1, kk * b * U1 / L, p["beta_b"] * b * U1


def simulate(t, f, a, b, p, U0=0.0):
    """the model along the gait of a run (arrays at the times t): U, omega, v_slip at t"""
    U = np.zeros_like(t); W = np.zeros_like(t); V = np.zeros_like(t)
    u = U0
    for i in range(len(t)):
        if i > 0:
            dt = t[i] - t[i - 1]
            u, _, _ = step(u, f[i - 1], a[i - 1], b[i - 1], t[i - 1], dt, p)
        _, w, v = step(u, f[i], a[i], b[i], t[i], 0.0, p)
        U[i], W[i], V[i] = u, w, v
    return U, W, V


def power(p, f, a, U):
    return p.get("c3", 0.0) * (f * a) ** 3 + p.get("c1", 0.0) * (f * a) * U * U
