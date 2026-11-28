#!/usr/bin/env python3
# The swimming physiology of an Atlantic salmon for the agent fish: what a speed costs, the oxygen debt and the fatigue.
#
# The model (a fish's version of the critical-power / W' balance model of exercise physiology):
#  * The oxygen demand of swimming at U through the water, in units of the maximum metabolic rate MMR, rises
#    exponentially with the speed (the exponential regressions of MO2 against U, Hvas & Oppedal 2019):
#        m(U) = r_s exp(beta U/U_crit),   r_s = SMR/MMR (0.25: SMR 123-137, MMR 511-531 mg O2/kg/h, post-smolts),
#    with beta = ln(1/r_s)/u_a: the demand reaches MMR at the aerobic limit U_a = u_a U_crit.
#  * Above it (m > 1) the fish swims on anaerobic power and builds an oxygen debt D [s at MMR]:
#        dD/dt = m - 1.
#  * Below it the fish pays the debt back from its spare aerobic scope 1 - m: all of it while the debt is large, a
#    first-order recovery below tau_r:
#        dD/dt = -(1 - m) min(1, D/tau_r),
#    and its oxygen uptake is MO2 = MMR (m + (1 - m) min(1, D/tau_r)) (the excess post-exercise oxygen consumption).
#  * The fish is exhausted when its debt reaches its anaerobic capacity C. It can then swim no faster than U_a until it
#    has paid the debt back to f_resume C (it falls back against the rear grid of a tunnel, the rear net of a cage).
#  * The fatigue F = D/C (0 rested, 1 exhausted).
# U_crit is the critical swimming speed as measured by the protocol of Brett (1964). By default from the fit of 390
# measurements in one tunnel (Hvas et al. 2021, Reviews in Aquaculture 13, Fig. 3; data of Remen et al. 2016a, Hvas &
# Oppedal 2017, Hvas et al. 2017a, c):
#        U_crit [cm/s] = 52.2020 + 0.0203 W + 4.0963 T - 6.1904e-6 W^2 - 0.1331 T^2   (W [g], T [C]),
# held at its maximum in W beyond 1640 g. The individual fish differ: their U_crit with a coefficient of variation
# cv_ucrit, their capacity C with cv_C (log-normal).
# The constants u_a, C, tau_r, cv_C and cv_ucrit are fitted (calibrate.py, validation 32) to the fluctuating currents
# of Athammer et al. (2024, Aquaculture Environment Interactions 16, 103-114), with the U_crit test of the same fish
# held to its measured value, and to a steady 80 % U_crit sustained for 4 h (Hvas & Oppedal 2017): u_a 0.785,
# C 1300 s (the debt of 22 min at MMR, about 190 mg O2/kg, on average), tau_r 450 s. The fish differ in their capacity
# (cv_C 0.70), not in their U_crit (cv_ucrit 0): the spread of the U_crit test (7 %; 9 % in the model) is that of
# the capacity.
# "gamma" > 0: the demand as SMR plus a power of the speed, m = r_s + (1 - r_s) (U/U_a)^gamma, instead of the
# exponential (fitted: gamma 2.7, a poorer fit).
#
#   import physiology as ph
#   p = ph.Physiology(n=100, ucrit=ph.ucrit(W=800, T=9) , seed=1)   # or ucrit=[...] per fish
#   for each step:  st = p.step(U, dt)        # U: the speed through the water of each fish [m/s]
#                   st["F"], st["MO2"], st["exhausted"], p.u_max()   (the speed each fish can reach now)
#
# Architect: Hans Bihs
import math
import numpy as np

# the constants of calibrate.py (validation 32): fitted to Athammer et al. (2024)
DEFAULT = {"r_s": 0.25, "u_a": 0.785, "C": 1300.0, "tau_r": 450.0, "cv_ucrit": 0.0, "cv_C": 0.70, "f_resume": 0.5,
           "MMR": 520.0, "gamma": 0.0}


def ucrit_cm(W, T):
    """the critical swimming speed [cm/s] of an Atlantic salmon of W grams at T C (Hvas et al. 2021, Fig. 3)"""
    W = np.minimum(np.asarray(W, float), 0.0203 / (2 * 6.1904e-6))
    return 52.2020 + 0.0203 * W + 4.0963 * T - 6.1904e-6 * W ** 2 - 0.1331 * np.asarray(T, float) ** 2


def ucrit(W, T):
    """the same in m/s"""
    return ucrit_cm(W, T) / 100.0


class Physiology:
    """the oxygen debt and the fatigue of n fish (see the header)"""

    def __init__(self, n=1, ucrit=1.0, seed=None, vary=True, **kw):
        self.p = dict(DEFAULT, **kw)
        p = self.p
        rnd = np.random.default_rng(seed)
        uc = np.broadcast_to(np.asarray(ucrit, float), (n,)).copy()
        if vary and p["cv_ucrit"] > 0.0 and np.ndim(ucrit) == 0:
            uc *= np.clip(1.0 + p["cv_ucrit"] * rnd.standard_normal(n), 0.5, 1.5)
        self.ucrit = uc
        s = math.sqrt(math.log(1.0 + p["cv_C"] ** 2)) if (vary and p["cv_C"] > 0.0) else 0.0
        self.C = p["C"] * np.exp(s * rnd.standard_normal(n) - 0.5 * s * s)
        self.beta = math.log(1.0 / p["r_s"]) / p["u_a"]
        self.D = np.zeros(n)
        self.exhausted = np.zeros(n, bool)
        self.t_exhausted = np.full(n, np.nan)       # the first time each fish was exhausted
        self.t = 0.0

    def demand(self, U):
        """the oxygen demand of swimming at U [MMR]: exponential (gamma 0) or SMR plus a power of the speed"""
        u = np.abs(U) / self.ucrit
        if self.p.get("gamma", 0.0) > 0.0:
            return self.p["r_s"] + (1.0 - self.p["r_s"]) * (u / self.p["u_a"]) ** self.p["gamma"]
        return self.p["r_s"] * np.exp(self.beta * u)

    def u_aerobic(self):
        return self.p["u_a"] * self.ucrit

    def u_max(self):
        """the fastest each fish can swim now: U_a when exhausted, else no limit (inf)"""
        return np.where(self.exhausted, self.u_aerobic(), np.inf)

    def step(self, U, dt):
        p = self.p
        m = self.demand(U)
        up = m > 1.0
        # above the aerobic limit the debt grows; below it is paid back at (1 - m) min(1, D/tau_r): all the spare scope
        # while D > tau_r, then first order (exact over dt)
        sp = np.clip(1.0 - m, 0.0, None)
        lin = self.D - sp * dt
        Dn = np.where(up, self.D + (m - 1.0) * dt,
                      np.where(lin >= p["tau_r"], lin, np.minimum(self.D, p["tau_r"]) * np.exp(-sp * dt / p["tau_r"])))
        mo2 = np.where(up, 1.0, m + sp * np.minimum(self.D / p["tau_r"], 1.0))
        self.D = Dn
        self.t += dt
        new = (~self.exhausted) & (self.D >= self.C)
        self.t_exhausted = np.where(new & np.isnan(self.t_exhausted), self.t, self.t_exhausted)
        self.exhausted = (self.exhausted | new) & ~(self.exhausted & (self.D <= p["f_resume"] * self.C))
        return {"m": m, "MO2": p["MMR"] * mo2, "F": np.clip(self.D / self.C, 0.0, None), "exhausted": self.exhausted.copy()}


# ---------------------------------------------------------------------------------------------------------------
# the tunnel protocols
def brett(n=1, ucrit=0.945, u0=0.15, du=0.15, dt_step=1800.0, dt=1.0, t_max=6 * 3600.0, seed=1, vary=True, **kw):
    """the U_crit test of Brett (1964) on the model: from u0 [m/s], du every dt_step [s] until each fish is
    exhausted; U_crit = U_f + du t_f/dt_step (U_f the last speed completed, t_f the time at the speed of fatigue)"""
    ph = Physiology(n, ucrit, seed=seed, vary=vary, **kw)
    t = 0.0
    while t < t_max and np.isnan(ph.t_exhausted).any():
        k = int(t // dt_step)
        ph.step(np.full(n, u0 + k * du), dt)
        t += dt
    te = np.where(np.isnan(ph.t_exhausted), t_max, ph.t_exhausted)
    k = np.floor(te / dt_step)
    return (u0 + (k - 1) * du) + du * (te - k * dt_step) / dt_step, ph


def steady(u, T=4 * 3600.0, ucrit=0.945, n=1000, dt=1.0, seed=1, vary=True, **kw):
    """a steady speed u [m/s] for T [s]: the time each fish was exhausted (nan: it endured)"""
    ph = Physiology(n, ucrit, seed=seed, vary=vary, **kw)
    for _ in range(int(round(T / dt))):
        ph.step(np.full(n, u), dt)
    return ph.t_exhausted.copy(), ph


def wave_exponent(mean_rel):
    """p such that the mean of s^p, s = (1 - cos 2 pi phi)/2, is mean_rel (0 < mean_rel <= 0.5)"""
    lo, hi = 1e-3, 50.0
    g = lambda p: math.exp(math.lgamma(p + 0.5) - math.lgamma(p + 1.0)) / math.sqrt(math.pi)
    for _ in range(80):
        mid = math.sqrt(lo * hi)
        if g(mid) > mean_rel:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi)


def fluctuating(u_peak, u_min, u_mean, period, ucrit=0.945, n=1000, T=4 * 3600.0, dt=0.5, seed=1, vary=True, **kw):
    """fluctuating currents (Athammer et al. 2024): the speed u_min + (u_peak - u_min) s^p over each period, p set by
    the mean speed the fish met; the fish holds its station (swims at the current). Returns the time each fish was
    exhausted (nan: it endured T) and the model"""
    p_ = wave_exponent((u_mean - u_min) / (u_peak - u_min))
    ph = Physiology(n, ucrit, seed=seed, vary=vary, **kw)
    t = 0.0
    while t < T - 1e-9:
        s = 0.5 * (1.0 - math.cos(2.0 * math.pi * (t + 0.5 * dt) / period))
        ph.step(np.full(n, u_min + (u_peak - u_min) * s ** p_), dt)
        t += dt
    return ph.t_exhausted.copy(), ph


# the trials of Athammer et al. (2024), Table 1 (the speeds the fish met, m/s), the period [s], and the result:
# the share that endured 4 h, the mean time of fatigue of the others [min] (n = 18 per group)
ATHAMMER = [
    {"name": "80 %, 1 min", "peak": 0.79, "min": 0.23, "mean": 0.49, "period": 60.0, "endured": 1.0, "t_fat": None},
    {"name": "100 %, 1 min", "peak": 1.09, "min": 0.25, "mean": 0.63, "period": 60.0, "endured": 1.0, "t_fat": None},
    {"name": "120 %, 1 min", "peak": 1.20, "min": 0.24, "mean": 0.67, "period": 60.0, "endured": 4 / 18, "t_fat": 119.0},
    {"name": "140 %, 1 min", "peak": 1.41, "min": 0.23, "mean": 0.75, "period": 60.0, "endured": 0.0, "t_fat": 51.0},
    {"name": "120 %, 0.5 min", "peak": 1.25, "min": 0.23, "mean": 0.68, "period": 30.0, "endured": 3 / 18, "t_fat": 113.0},
    {"name": "120 %, 2 min", "peak": 1.24, "min": 0.19, "mean": 0.68, "period": 120.0, "endured": 2 / 18, "t_fat": 123.0},
]
# a steady speed: post-smolts sustained 80 % of their U_crit for 4 h (Hvas & Oppedal 2017; Hvas et al. 2021)
SUSTAINED = {"name": "80 % steady, 4 h", "u_rel": 0.80, "T": 4 * 3600.0, "endured": 1.0}
ATHAMMER_UCRIT = 0.945          # the U_crit test of the same fish (+0.15 m/s every 30 min from 0.15 m/s), SD 0.066
