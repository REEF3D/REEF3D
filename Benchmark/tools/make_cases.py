#!/usr/bin/env python3
"""
Writes the benchmark case directories (cases/<name>/case.json, control.txt, ctrl.txt).

The case files in cases/ are the source of truth and can be edited by hand; this script documents
how they were set up (gauge lists, geometry) and can regenerate them:  ./tools/make_cases.py
"""

import json
import math
import os
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CASES = os.path.join(ROOT, "cases")
TUT = os.environ.get("REEF3D_TUTORIALS", os.path.join(ROOT, "..", "Tutorials"))


def write_case(name, control, ctrl, case, files=None):
    d = os.path.join(CASES, name)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "control.txt"), "w") as f:
        f.write(control.strip() + "\n")
    with open(os.path.join(d, "ctrl.txt"), "w") as f:
        f.write(ctrl.strip() + "\n")
    with open(os.path.join(d, "case.json"), "w") as f:
        json.dump(case, f, indent=1)
        f.write("\n")
    for src, dst in (files or {}).items():
        shutil.copy(src, os.path.join(d, dst))


def gauges(xs, y, fmt="P 51 %.4f %.4f"):
    return "\n".join(fmt % (x, y) for x in xs)


def frange(a, b, h):
    n = int(round((b - a) / h))
    return [a + i * h for i in range(n + 1)]


# references (short label, full citation)
REF = {
    "beji": ("Beji & Battjes (1993)",
             "Beji, S. & Battjes, J.A. (1993) Experimental investigation of wave propagation over a bar. Coastal "
             "Engineering 19, 151-162; Luth, Klopman & Kitou (1994) / Dingemans (1994). Measured series: secondary, "
             "digitised (Basilisk test 'bar', gauge-4..11), see refdata/waves/PROVENANCE.md"),
    "berkhoff": ("Berkhoff, Booy & Radder (1982)",
                 "Berkhoff, J.C.W., Booy, N. & Radder, A.C. (1982) Verification computations with linear wave "
                 "propagation models. Coastal Engineering 6, 255-279. Sections 2, 3, 5, 7: secondary, digitised "
                 "(Basilisk example 'shoal'), see refdata/waves/PROVENANCE.md"),
    "synolakis": ("Synolakis (1987)",
                  "Synolakis, C.E. (1987) The runup of solitary waves. J. Fluid Mech. 185, 523-545; run-up law as "
                  "quoted in Synolakis et al. (2008) PAGEOPH 165 (NOAA NCTR benchmark 'canonical bathymetry')"),
    "ting": ("Ting & Kirby (1995)",
             "Ting, F.C.K. & Kirby, J.T. (1995) Dynamics of surf-zone turbulence in a strong plunging breaker. "
             "Coastal Engineering 24, 177-204 (plunging case H = 0.128 m, T = 5 s, 1:35 slope); breaking point "
             "x_b = 7.795 m from Derakhti et al. (2016) Table 1, see refdata/waves/ting_kirby_breaking/setup.txt"),
    "fenton": ("Fenton (1985) 5th-order Stokes",
               "Fenton, J.D. (1985) A fifth-order Stokes theory for steady waves. J. Waterway, Port, Coastal and "
               "Ocean Eng. 111(2), 216-234 (the generating theory, B 92 5, recorded with P 50)"),
    "stokes2": ("2nd-order Stokes theory",
                "Stokes second-order wave theory (e.g. Dean & Dalrymple 1991), recorded with P 50"),
    "ritter": ("Ritter (1892), SWASHES",
               "Ritter, A. (1892) Die Fortpflanzung der Wasserwellen. Z. Verein Deutscher Ing. 36; exact solution "
               "as in Delestre et al. (2013) SWASHES, Int. J. Numer. Meth. Fluids 72, 269-300"),
    "stoker": ("Stoker (1957), SWASHES",
               "Stoker, J.J. (1957) Water Waves. Interscience; exact solution as in Delestre et al. (2013) SWASHES"),
    "solitary": ("Solitary wave celerity",
                 "Boussinesq solitary wave, c = sqrt(g (d + H)) (e.g. Synolakis 1987; Laitone 1960)"),
    "martin": ("Martin & Moyce (1952)",
               "Martin, J.C. & Moyce, W.J. (1952) An experimental study of the collapse of liquid columns on a "
               "rigid horizontal plane. Phil. Trans. R. Soc. A 244, 312-324; n^2 = 2, a = 2.25 in, surge front "
               "Z(T) digitised from Fig. 3 (secondary: PySPH db_exp_data.py), see refdata/dambreak/PROVENANCE.md"),
    "kleefsman": ("Kleefsman et al. (2005)",
                  "Kleefsman, K.M.T. et al. (2005) A volume-of-fluid based simulation method for wave impact "
                  "problems. J. Comput. Phys. 206, 363-393 (MARIN, SPHERIC Test 2). P1/P3 pressures: secondary, "
                  "digitised (PySPH), see refdata/dambreak/PROVENANCE.md"),
    "cylinder": ("Williamson St(Re); Qu et al. (2013)",
                 "Williamson (1988) St-Re fit as given by Beaudan & Moin (1994) Eq. 1 (50 < Re < 180); mean C_D at "
                 "Re = 100: Park et al. (1998) 1.33, Posdziech & Grundmann (2007) 1.325, Qu et al. (2013) 1.319 "
                 "(Qu et al. 2013 J. Fluids Struct. 39, Table 1), see refdata/misc/cylinder_re100_unconfined.txt"),
    "chen": ("Chen et al. (2014)",
             "Chen, L.F., Zang, J., Hillis, A.J., Morgan, G.C.J. & Plummer, A.R. (2014) Numerical investigation of "
             "wave-structure interaction using OpenFOAM. Ocean Engineering 88, 91-109; measured data from the "
             "REEF3D tutorial 11_14 (Bihs et al. 2016, Computers & Fluids 140)"),
    "irschik": ("Irschik et al. (2002)",
                "Irschik, K., Sparboom, U. & Oumeraci, H. (2002) Breaking wave characteristics for the loading of "
                "a slender pile. Proc. 28th ICCE; GWK Hannover; measured data from the REEF3D tutorial 11_15 "
                "(Kamath et al. 2016, Ocean Engineering)"),
    "hulme": ("Hulme (1982) linear theory",
              "Hulme, A. (1982) The wave forces acting on a floating hemisphere undergoing forced periodic "
              "oscillations. J. Fluid Mech. 121, 443-463 (added mass and damping, table via LHEEA "
              "Hulme_Heaving_Sphere), see refdata/misc/hulme_hemisphere_heave.txt; equilibrium: Archimedes"),
    "scour": ("HEC-18 / Sumer et al. (1992)",
              "Arneson et al. (2012) FHWA HEC-18 (CSU pier scour equation); Sumer, Fredsoe & Christiansen (1992) "
              "J. Waterway, Port, Coastal, Ocean Eng. 118, S/D = 1.3 (sigma 0.7) for live-bed current scour; see "
              "refdata/misc/pier_scour.txt"),
}


def ref(k):
    return {"short": REF[k][0], "citation": REF[k][1]}


# ==============================================================================================
# SFLOW
# ==============================================================================================

def sflow_dambreak(name, hr, key):
    """1D dam break, flat frictionless bed, hydrostatic SFLOW. hl = 0.5 m, dam at x0 = 20 m."""
    hl, x0 = 0.5, 20.0
    gx = [12.0, 16.0, 19.0, 21.0, 24.0, 28.0, 32.0, 36.0]
    control = """
C 11 21
C 12 3
C 13 3
C 14 21
C 15 21
C 16 21

B 1 0.01
B 10 0.0 40.0 0.0 0.01 0.0 1.0

M 10 1
M 20 2
"""
    ctrl = """
A 10 2
A 210 3
A 211 4
A 220 0
A 217 1
F 60 %(f60)g
F 72 0.0 %(x0)g -1.0 1.0 %(hl)g
N 41 4.0
N 47 0.25
M 10 1
P 12 100
W 22 -9.81
%(g)s
""" % {"f60": max(hr, 1.0e-4), "hl": hl, "x0": x0, "g": gauges(gx, 0.005)}
    case = {
        "description": "1D dam break on a flat frictionless bed, h_l = %g m, h_r = %g m, dam at x = %g m, hydrostatic "
                       "SFLOW (A 220 0), WENO reconstruction; depth time series at 8 gauges for 0 < t < 4 s against "
                       "the exact %s solution." % (hl, hr, x0, "Ritter (dry bed)" if hr == 0 else "Stoker (wet bed)"),
        "reference": ref(key),
        "tags": ["sflow", "analytical", "1d", "quick"],
        "time": 4.0,
        "check": {"type": "swe_dambreak", "x0": x0, "hl": hl, "hr": hr, "z_bed": -max(hr, 1.0e-4), "datum": 0.0,
                  "t_start": 0.0, "t_end": 4.0, "tol": 0.025 if hr == 0 else 0.002},
        "levels": {
            "nightly": {"control_set": {"B 1": "0.05", "B 10": "0.0 40.0 0.0 0.05 0.0 1.0"},
                        "ctrl_set": {"P 51": [("%.4f 0.025" % x) for x in gx]},
                        "check_set": {"tol": 0.03 if hr == 0 else 0.005}},
            "release": {"np": 2}},
        "baseline": {"nightly": 0.0212 if hr == 0 else 0.0026, "release": 0.0198 if hr == 0 else 0.0005,
                     "commit": "ff1bb68"},
    }
    write_case(name, control, ctrl, case)


def sflow_solitary():
    """Tutorial 8_3: solitary wave H/d = 0.1 on d = 0.5 m, non-hydrostatic SFLOW."""
    control = """
C 11 6
C 12 3
C 13 3
C 14 7
C 15 21
C 16 3

B 1 0.01
B 10 0.0 40.0 0.0 0.01 0.0 1.0

M 10 1
M 20 2
"""
    ctrl = """
A 10 2
A 210 3
A 211 4
A 216 2
A 220 2
A 223 0.5
A 240 1
A 241 1
A 243 0
A 246 0
B 90 1
B 92 9
B 91 0.05 4.0
B 96 4.0 8.0
B 98 2
B 99 1
F 60 0.5
N 41 18.0
N 47 0.2
M 10 1
P 12 100
P 51 8.0 0.005
P 51 16.0 0.005
P 51 24.0 0.005
W 22 -9.81
"""
    case = {
        "description": "Solitary wave H/d = 0.1 (d = 0.5 m) propagating over constant depth with the non-hydrostatic "
                       "SFLOW (A 220 2), set-up of tutorial 8_3 extended to 40 m: celerity between x = 8 and 24 m "
                       "against c = sqrt(g(d+H)) and crest height at x = 24 m.",
        "reference": ref("solitary"),
        "tags": ["sflow", "analytical", "1d", "quick"],
        "time": 18.0,
        "check": {"type": "solitary", "x1": 8.0, "x2": 24.0, "d": 0.5, "H": 0.05, "tol_c": 0.01, "tol_H": 0.03},
        "levels": {
            "nightly": {"control_set": {"B 1": "0.02", "B 10": "0.0 40.0 0.0 0.02 0.0 1.0"},
                        "ctrl_set": {"P 51": ["8.0 0.01", "16.0 0.01", "24.0 0.01"]},
                        "check_set": {"tol_c": 0.01, "tol_H": 0.03}},
            "release": {"np": 2}},
    }
    write_case("sflow_solitary_celerity", control, ctrl, case)


# Beji & Battjes: measured gauges WG4..WG11 in the Basilisk frame (bar toe at x = 6 m)
BB_GAUGES = [("WG4", 10.5), ("WG5", 12.5), ("WG6", 13.5), ("WG7", 14.5), ("WG8", 15.7), ("WG9", 17.3),
             ("WG10", 19.0), ("WG11", 21.0)]


def bb_signals(shift):
    return [{"name": n, "source": "gauge", "index": i + 1, "data": "waves/beji_battjes_bar/%s_eta_timeseries.txt" % n,
             "data_y_factor": 0.01} for i, (n, x) in enumerate(BB_GAUGES)]


def bb_check(tol_rms, tol_h):
    return {"type": "timeseries", "signals": bb_signals(0), "align": [n for n, _ in BB_GAUGES], "lag_range": [-1.2, 1.2],
            "lag_step": 0.005, "local_lag": 0.1, "tol_rms": tol_rms, "tol_height": tol_h,
            "note": "WG4..WG11 at x = 10.5..21 m from the wave board (bar toe at 6 m); the measured time origin "
                    "is arbitrary, one common lag is fitted on all gauges, plus at most +-0.1 s (T/20) per gauge for the gauge "
                    "position and digitising uncertainty"}


def sflow_beji():
    shift = 5.0   # tutorial 8_4 geometry: bar toe at x = 11 m
    control = """
C 11 6
C 12 3
C 13 3
C 14 7
C 15 21
C 16 3

B 1 0.01
B 10 0.0 40.0 0.0 0.01 0.0 0.8

S 61 11.0 17.0 0.0 0.01 0.0 0.3
S 10 17.0 19.0 0.0 0.01 0.0 0.3
S 61 19.0 22.0 0.0 0.01 0.3 0.0

M 10 1
M 20 2
"""
    ctrl = """
A 10 2
A 220 3
B 90 1
B 92 4
B 93 0.02 2.02
B 96 5.0 9.0
B 98 2
B 99 1
F 60 0.40
N 41 42.0
N 47 0.2
M 10 1
P 12 100
W 10 0.0
W 22 -9.81
%s
""" % gauges([x + shift for _, x in BB_GAUGES], 0.005)
    case = {
        "description": "Regular waves (H = 2 cm, T = 2.02 s) over the Beji & Battjes submerged bar, non-hydrostatic "
                       "SFLOW with the improved-dispersion quadratic pressure (A 220 3), geometry of tutorial 8_4 (bar toe at x = 11 m, i.e. +5 m against the experiment). "
                       "Surface elevation at WG4..WG11 against the measured series (shape incl. the released higher "
                       "harmonics on the lee side).",
        "reference": ref("beji"),
        "tags": ["sflow", "experiment", "1d", "dispersion"],
        "time": 42.0,
        "check": bb_check(0.35, 0.2),
        "levels": {
            "nightly": {"control_set": {"B 1": "0.025", "B 10": "0.0 40.0 0.0 0.025 0.0 0.8",
                                        "S 61": ["11.0 17.0 0.0 0.025 0.0 0.3", "19.0 22.0 0.0 0.025 0.3 0.0"],
                                        "S 10": "17.0 19.0 0.0 0.025 0.0 0.3"},
                        "ctrl_set": {"P 51": ["%.4f 0.0125" % (x + shift) for _, x in BB_GAUGES]},
                        "check_set": {"tol_rms": 0.35, "tol_height": 0.2}},
            "release": {"np": 2}},
    }
    write_case("sflow_beji_battjes_bar", control, ctrl, case)


# ==============================================================================================
# NHFLOW
# ==============================================================================================

def nhflow_beji():
    control = """
C 11 6
C 12 3
C 13 3
C 14 7
C 15 21
C 16 3

B 2 1500 1 10
B 10 0.0 35.0 0.0 0.01 0.0 1.0

S 61 6.0 12.0 0.0 0.01 0.0 0.3
S 10 12.0 14.0 0.0 0.01 0.0 0.3
S 61 14.0 17.0 0.0 0.01 0.3 0.0

B 103 5
B 113 2.5
B 116 1.0

M 10 4
"""
    ctrl = """
A 10 5
B 90 1
B 92 4
B 93 0.02 2.02
B 96 3.73 8.73
B 98 2
B 99 2
F 60 0.4
N 41 42.0
N 47 0.5
M 10 4
P 12 100
W 22 -9.81
%s
""" % gauges([x for _, x in BB_GAUGES], 0.005)
    case = {
        "description": "Regular waves (H = 2 cm, T = 2.02 s) over the Beji & Battjes submerged bar with NHFLOW, "
                       "geometry of tutorial NHFLOW/2 (bar toe at x = 6 m as in the experiment). Surface elevation "
                       "at WG4..WG11 against the measured series.",
        "reference": ref("beji"),
        "tags": ["nhflow", "experiment", "2dv", "dispersion"],
        "np": 4,
        "time": 42.0,
        "check": bb_check(0.35, 0.3),
        "levels": {
            "nightly": {"np": 2, "control_set": {"B 2": "1000 1 8", "M 10": "2"},
                        "check_set": {"tol_rms": 0.45, "tol_height": 0.3}},
            "release": {}},
    }
    write_case("nhflow_beji_battjes_bar", control, ctrl, case)


def nhflow_stokes5():
    control = """
C 11 6
C 12 3
C 13 3
C 14 7
C 15 21
C 16 3

B 2 800 1 10
B 10 0.0 200.0 0.0 0.05 0.0 1.0

B 103 5
B 113 2.5
B 116 1.0

M 10 4
M 20 2
"""
    ctrl = """
A 10 5
B 90 1
B 92 5
B 93 1.0 4.5
B 96 25.0 50.0
B 98 2
B 99 1
F 60 4.01
I 12 1
N 41 80.0
N 47 0.5
M 10 4
P 12 100
P 50 50.0 0.025
P 50 75.0 0.025
P 50 100.0 0.025
P 50 125.0 0.025
P 51 50.0 0.025
P 51 75.0 0.025
P 51 100.0 0.025
P 51 125.0 0.025
W 22 -9.81
"""
    case = {
        "description": "5th-order Stokes waves (Fenton), H = 1 m, T = 4.5 s, d = 4.01 m (kd ~ 1, H/L ~ 0.04) "
                       "propagating 100 m in NHFLOW (tutorial NHFLOW/1 extended). The simulated surface elevation "
                       "must follow the generating theory (P 50) at x = 50..125 m: no loss of height, no phase "
                       "drift, the crest/trough asymmetry kept.",
        "reference": ref("fenton"),
        "tags": ["nhflow", "analytical", "2dv", "quick"],
        "np": 4,
        "time": 80.0,
        "check": {"type": "theory_gauges", "t_start": 50.0, "t_end": 80.0, "tol_rms": 0.10, "tol_height": 0.03},
        "levels": {
            "nightly": {"np": 2, "time": 60.0, "control_set": {"B 2": "400 1 8", "M 10": "2"},
                        "check_set": {"t_start": 45.0, "t_end": 60.0, "tol_rms": 0.12, "tol_height": 0.04}},
            "release": {}},
    }
    write_case("nhflow_stokes5_propagation", control, ctrl, case)


def synolakis_geometry(d, Hd, dx, x_toe, x_end, z_top):
    """plane beach 1:19.85, toe at x_toe, still water depth d; returns (x_shore, S 61 line)"""
    cot = 19.85
    x_shore = x_toe + cot * d
    z_end = (x_end - x_toe) / cot
    return x_shore, "S 61 %.3f %.3f 0.0 %g 0.0 %.4f" % (x_toe, x_end, dx, z_end)


def nhflow_synolakis():
    d, Hd = 1.0, 0.0185
    x_toe, x_end = 50.0, 79.0
    x_shore, s61 = synolakis_geometry(d, Hd, 0.05, x_toe, x_end, 1.6)
    control = """
C 11 6
C 12 3
C 13 3
C 14 21
C 15 21
C 16 3

B 1 0.05
B 10 0.0 %(xe)g 0.0 0.05 0.0 1.6
%(s61)s

M 10 4
""" % {"xe": x_end, "s61": s61}
    ctrl = """
A 10 5
B 90 1
B 92 9
B 93 %(H)g 10.0
B 96 10.0 0.0
B 98 2
B 99 0
F 60 %(d)g
N 41 45.0
N 47 0.5
M 10 4
P 12 100
P 51 20.0 0.025
P 51 %(xt)g 0.025
P 133 0.025
P 134 0.025
W 22 -9.81
""" % {"H": Hd * d, "d": d, "xt": x_toe}
    case = {
        "description": "Non-breaking solitary wave H/d = 0.0185 (d = 1 m) running up a 1:19.85 plane beach "
                       "(Synolakis 1987 canonical case), NHFLOW with wetting and drying. Maximum run-up (P 134) "
                       "against the run-up law R/d = 2.831 sqrt(cot b)(H/d)^1.25 = 0.0862 (valid for non-breaking "
                       "waves; breaking starts at H/d = 0.045 on this slope).",
        "reference": ref("synolakis"),
        "tags": ["nhflow", "analytical", "2dv", "runup"],
        "np": 4,
        "time": 45.0,
        "check": {"type": "runup_law", "d": d, "H_over_d": Hd, "cot_beta": 19.85, "swl": d, "tol": 0.08},
        "levels": {
            "nightly": {"np": 2, "control_set": {"B 1": "0.1", "B 10": "0.0 %g 0.0 0.1 0.0 1.6" % x_end,
                                                 "S 61": "%.3f %.3f 0.0 0.1 0.0 %.4f" % (x_toe, x_end, (x_end - x_toe) / 19.85),
                                                 "M 10": "2"},
                        "ctrl_set": {"P 51": ["20.0 0.05", "%g 0.05" % x_toe], "P 133": "0.05", "P 134": "0.05"},
                        "check_set": {"tol": 0.12}},
            "release": {"control_set": {"B 1": "0.025", "B 10": "0.0 %g 0.0 0.025 0.0 1.6" % x_end,
                                        "S 61": "%.3f %.3f 0.0 0.025 0.0 %.4f" % (x_toe, x_end, (x_end - x_toe) / 19.85)},
                        "ctrl_set": {"P 51": ["20.0 0.0125", "%g 0.0125" % x_toe], "P 133": "0.0125",
                                     "P 134": "0.0125"}}},
        "x_shore": x_shore,
    }
    write_case("nhflow_synolakis_runup_nonbreaking", control, ctrl, case)

    # breaking case H/d = 0.28: profiles
    Hd = 0.28
    ctrl2 = """
A 10 5
A 550 1
B 90 1
B 92 9
B 93 %(H)g 10.0
B 96 10.0 0.0
B 98 2
B 99 0
F 60 %(d)g
N 41 40.0
N 47 0.5
M 10 4
P 12 100
P 52 0.025
P 55 0.02
P 133 0.025
P 134 0.025
W 22 -9.81
""" % {"H": Hd * d, "d": d}
    profiles = [[t, "waves/synolakis_runup/profile_Hd0.28_t%d.txt" % t] for t in range(10, 70, 5)]
    case2 = {
        "description": "Breaking solitary wave H/d = 0.28 (d = 1 m) on the 1:19.85 beach (Synolakis 1987), NHFLOW "
                       "with breaking (A 550 1) and wetting/drying: free-surface profiles at t sqrt(g/d) = 10..65 "
                       "(shoaling, breaking, run-up and run-down) against the measured profiles. The time origin is "
                       "fitted on the first profile.",
        "reference": ref("synolakis"),
        "tags": ["nhflow", "experiment", "2dv", "runup", "breaking"],
        "np": 4,
        "time": 40.0,
        "check": {"type": "profiles", "d": d, "H_over_d": Hd, "swl": d, "x_shore": x_shore, "profiles": profiles,
                  "tol": 0.22},
        "levels": {
            "nightly": {"np": 2, "control_set": {"B 1": "0.1", "B 10": "0.0 %g 0.0 0.1 0.0 1.6" % x_end,
                                                 "S 61": "%.3f %.3f 0.0 0.1 0.0 %.4f" % (x_toe, x_end, (x_end - x_toe) / 19.85),
                                                 "M 10": "2"},
                        "ctrl_set": {"P 52": "0.05", "P 133": "0.05", "P 134": "0.05"},
                        "check_set": {"tol": 0.2}},
            "release": {"control_set": {"B 1": "0.025", "B 10": "0.0 %g 0.0 0.025 0.0 1.6" % x_end,
                                        "S 61": "%.3f %.3f 0.0 0.025 0.0 %.4f" % (x_toe, x_end, (x_end - x_toe) / 19.85)},
                        "ctrl_set": {"P 52": "0.0125", "P 133": "0.0125", "P 134": "0.0125"}}},
    }
    write_case("nhflow_synolakis_breaking_profiles", control, ctrl2, case2)


BERKHOFF_SECTIONS = [
    {"name": "section 2", "x": 3.0, "sign": -1.0, "data": "waves/berkhoff_shoal/section2_H_over_H0.txt"},
    {"name": "section 3", "x": 5.0, "sign": -1.0, "data": "waves/berkhoff_shoal/section3_H_over_H0.txt"},
    {"name": "section 5", "x": 9.0, "sign": -1.0, "data": "waves/berkhoff_shoal/section5_H_over_H0.txt"},
    {"name": "section 7", "y": 0.0, "sign": 1.0, "data": "waves/berkhoff_shoal/section7_H_over_H0.txt"},
]


def berkhoff_gauges():
    """P 51 lines at the measured points; domain frame = shoal frame + (10, 10)"""
    out = []
    for sec in BERKHOFF_SECTIONS:
        rows = [l.split() for l in open(os.path.join(ROOT, "refdata", sec["data"])) if l.strip() and l[0] != "#"]
        for r in rows:
            s = float(r[0])
            x = sec["x"] if "x" in sec else sec["sign"] * s
            y = sec["y"] if "y" in sec else sec["sign"] * s
            out.append("P 51 %.4f %.4f" % (x + 10.0, y + 10.0))
    return "\n".join(out)


def berkhoff(model):
    """elliptic shoal, 3D: domain x' = -10..15, y' = -10..10 (shoal frame) -> 0..25, 0..20"""
    control = """
C 11 6
C 12 21
C 13 21
C 14 7
C 15 21
C 16 3

B 1 0.04
B 10 0.0 25.0 0.0 20.0 0.0 0.6
B 103 5
B 113 2.0
B 116 0.6

G 10 1
G 15 2
G 31 0

M 10 16
M 20 2
"""
    if model == "nhflow":
        mk = "A 10 5\nB 92 4\nN 47 0.5"
        sig = "B 2 625 500 10"
    else:
        mk = "A 10 3\nA 310 3\nA 311 4\nA 343 0\nB 92 4\nN 47 1.0"
        sig = "B 2 625 500 10"
    control = control.replace("B 1 0.04", sig)
    ctrl = """
%(mk)s
B 90 1
B 93 0.0464 1.0
B 96 1.5 3.0
B 98 2
B 99 1
F 60 0.45
N 41 30.0
M 10 16
P 12 100
W 22 -9.81
%(g)s
""" % {"mk": mk, "g": berkhoff_gauges()}
    case = {
        "description": "Berkhoff et al. (1982) elliptic shoal on a 1:50 slope rotated by 20 deg, regular waves "
                       "T = 1 s, H = 4.64 cm (3D refraction, diffraction, nonlinear focusing behind the shoal) with "
                       "%s. Wave heights H/H0 along sections 2, 3, 5 (across the tank at x = 3, 5, 9 m behind the "
                       "shoal centre) and 7 (centre line) against the measurements. Bathymetry from geo.dat written "
                       "by tools/berkhoff_geo.py (Basilisk shoal.c formula, minimum depth 0.07 m)." % model.upper(),
        "reference": ref("berkhoff"),
        "tags": [model, "experiment", "3d", "refraction"],
        "np": 16,
        "time": 30.0,
        "generate": {"script": "berkhoff_geo.py", "args": [0.04]},
        "check": {"type": "section_heights", "sections": BERKHOFF_SECTIONS, "dx": 10.0, "dy": 10.0, "H0": 0.0464,
                  "t_start": 22.0, "t_end": 30.0, "tol": 0.15,
                  "note": "section coordinates as in Basilisk shoal.c (sections 2-5 plotted against -y)"},
        "levels": {
            "nightly": {"np": 4, "control_set": {"B 2": "250 200 8", "M 10": "4"},
                        "ctrl_set": {"M 10": "4"}, "generate": {"script": "berkhoff_geo.py", "args": [0.1]},
                        "check_set": {"tol": 0.3}},
            "release": {}},
    }
    if model == "nhflow":
        case["levels"]["nightly"]["control_set"]["B 2"] = "250 200 5"
        case["levels"]["nightly"]["xfail"] = (
            "hans_dev ff1bb68, dx = 0.1 m (L/15): the wave heights are about 60 % of the measured ones on all "
            "sections, already at the first point of section 7 (H/H0 0.67 vs 1.07); 5 or 8 sigma layers give "
            "the same. FNPF on the same grid is within 0.17-0.28. Cause not investigated (horizontal numerical "
            "damping or 3D wave generation)")
    write_case("%s_berkhoff_shoal" % model, control, ctrl, case)


# ==============================================================================================
# FNPF
# ==============================================================================================

def fnpf_beji():
    shift = 5.0   # tutorial 9_2 geometry: bar toe at x = 11 m
    control = """
C 11 6
C 12 3
C 13 3
C 14 7
C 15 21
C 16 3

B 2 2000 1 10
B 10 0.0 40.0 0.0 0.02 0.0 1.0

S 61 11.0 17.0 0.0 0.02 0.0 0.3
S 10 17.0 19.0 0.0 0.02 0.0 0.3
S 61 19.0 22.0 0.0 0.02 0.3 0.0

B 103 5
B 113 3.0
B 116 1.0

M 10 2
M 20 2
"""
    ctrl = """
A 10 3
A 310 3
A 311 4
A 343 0
B 90 1
B 92 4
B 93 0.02 2.02
B 96 5.0 10.0
B 98 2
B 99 1
F 60 0.4
N 41 42.0
N 47 1.0
M 10 2
P 12 100
W 22 -9.81
%s
""" % gauges([x + shift for _, x in BB_GAUGES], 0.01)
    case = {
        "description": "Regular waves (H = 2 cm, T = 2.02 s) over the Beji & Battjes submerged bar with FNPF, "
                       "geometry of tutorial 9_2 (bar toe at x = 11 m, i.e. +5 m against the experiment). Surface "
                       "elevation at WG4..WG11 against the measured series.",
        "reference": ref("beji"),
        "tags": ["fnpf", "experiment", "2dv", "dispersion", "quick"],
        "np": 2,
        "time": 42.0,
        "check": bb_check(0.3, 0.25),
        "levels": {
            "nightly": {"np": 1, "control_set": {"B 2": "800 1 8", "M 10": "1", "B 10": "0.0 40.0 0.0 0.05 0.0 1.0",
                                                 "S 61": ["11.0 17.0 0.0 0.05 0.0 0.3", "19.0 22.0 0.0 0.05 0.3 0.0"],
                                                 "S 10": "17.0 19.0 0.0 0.05 0.0 0.3"},
                        "ctrl_set": {"M 10": "1", "P 51": ["%.4f 0.025" % (x + shift) for _, x in BB_GAUGES]},
                        "check_set": {"tol_rms": 0.3, "tol_height": 0.2}},
            "release": {}},
    }
    write_case("fnpf_beji_battjes_bar", control, ctrl, case)


def fnpf_stokes2():
    control = """
C 11 6
C 12 3
C 13 3
C 14 7
C 15 21
C 16 3

B 2 800 1 10
B 10 0.0 40.0 0.0 0.05 0.0 2.0
B 103 5
B 113 2.5
B 116 2.0

M 10 1
M 20 2
"""
    ctrl = """
A 10 3
A 310 3
A 311 4
A 343 0
B 90 1
B 92 4
B 91 0.04 4.0
B 96 4.0 8.0
B 98 2
B 99 1
F 60 2.0
N 41 50.0
N 47 1.0
M 10 1
P 12 100
P 50 15.0 0.025
P 50 20.0 0.025
P 50 25.0 0.025
P 51 15.0 0.025
P 51 20.0 0.025
P 51 25.0 0.025
W 22 -9.81
"""
    case = {
        "description": "2nd-order Stokes waves H = 0.04 m, L = 4 m in deep water (d = 2 m), FNPF, tutorial 9_1: "
                       "the simulated surface elevation must follow the generating theory (P 50) at x = 15..25 m.",
        "reference": ref("stokes2"),
        "tags": ["fnpf", "analytical", "2dv", "quick"],
        "check": {"type": "theory_gauges", "t_start": 35.0, "t_end": 50.0, "tol_rms": 0.15, "tol_height": 0.02,
                  "note": "hans_dev ff1bb68: the simulated wave is about 0.4 % faster than the theory (phase lead "
                          "0.017 s at x = 15 m, 0.032 s at 25 m, at both levels), height within 1 %"},
        "time": 50.0,
        "levels": {
            "nightly": {"time": 40.0, "control_set": {"B 2": "400 1 8"},
                        "check_set": {"t_start": 30.0, "t_end": 40.0, "tol_rms": 0.17, "tol_height": 0.03}},
            "release": {}},
    }
    write_case("fnpf_stokes2_propagation", control, ctrl, case)


def fnpf_ting_kirby():
    # tutorial 9_3: slope from x = 5.8 m (depth 0.4 m), depth 0.38 m at x = 6.5 m (Ting & Kirby x = 0)
    x0 = 5.8 + 0.02 * 35.0
    xs = frange(8.0, 18.0, 0.1)
    control = """
C 11 6
C 12 3
C 13 3
C 14 7
C 15 21
C 16 3

B 2 600 1 10
B 10 0.0 30.0 0.0 0.025 0.0 0.748
B 103 5
B 113 1.0
B 116 0.4

S 61 5.8 32.0 0.0 0.025 0.0 0.748

M 10 2
M 20 2
"""
    ctrl = """
A 10 3
A 341 2.0
A 343 1
A 350 1
A 351 3
A 352 3
A 365 0.0025
B 90 1
B 92 8
B 93 0.128 5.0
B 96 9.5 9.0
B 98 3
B 99 1
F 60 0.4
N 41 100.0
N 47 1.0
M 10 2
P 12 100
W 22 -9.81
%s
""" % gauges(xs, 0.0125)
    case = {
        "description": "Plunging breaker of Ting & Kirby (1995): cnoidal waves H = 0.128 m, T = 5 s on a 1:35 slope, "
                       "FNPF with the viscosity breaking model (tutorial 9_3 at its resolution dx = 0.05 m, Dirichlet generation). Breaking point "
                       "= location of the maximum wave height from 101 gauges between x = 8 and 18 m, against the "
                       "measured breaking point x_b = 7.795 m from the x = 0 station (d = 0.38 m, here x = %.2f m)."
                       % x0,
        "reference": ref("ting"),
        "tags": ["fnpf", "experiment", "2dv", "breaking"],
        "np": 2,
        "time": 100.0,
        "check": {"type": "breaking_point", "x_b": 7.795, "x_offset": x0, "t_start": 70.0, "t_end": 100.0,
                  "tol": 0.4},
        "levels": {
            "nightly": {"np": 1, "time": 70.0, "control_set": {"M 10": "1"},
                        "ctrl_set": {"M 10": "1"}, "check_set": {"t_start": 45.0, "t_end": 70.0, "tol": 0.8}},
            "release": {}},
        "xfail": "on hans_dev ff1bb68 (dx = 0.05 m, tutorial resolution) the maximum wave height is reached at "
                 "x = 12.4 m (d = 0.21 m, H = 0.20 m), 1.9 m before the measured breaking point (d = 0.157 m); "
                 "with dx = 0.025 m the run stops at t = 91 s (N 61 velocity limit) after the wave train becomes "
                 "irregular at t ~ 75 s",
    }
    write_case("fnpf_ting_kirby_plunging", control, ctrl, case)


# ==============================================================================================
# CFD
# ==============================================================================================

def cfd_martin_moyce():
    a = 0.05715  # 2.25 in
    Zs = [1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    L = 14 * a
    control = """
C 11 21
C 12 3
C 13 3
C 14 21
C 15 21
C 16 21

B 1 %(dx)g
B 10 0.0 %(L).5f 0.0 %(dx)g 0.0 %(H).5f

M 10 2
M 20 2
""" % {"dx": a / 40, "L": L, "H": 3 * a}
    ctrl = """
D 10 4
D 20 2
D 30 1
F 30 3
F 40 3
F 54 %(a).5f
F 56 %(h).5f
N 40 3
N 41 0.5
N 47 0.2
M 10 2
P 12 100
T 10 0
W 22 -9.81
%(g)s
""" % {"a": a, "h": 2 * a, "g": gauges([Z * a for Z in Zs], a / 80)}
    case = {
        "description": "2D collapse of a water column a x 2a (n^2 = 2, a = 2.25 in = 0.05715 m) in a closed tank, "
                       "two-phase CFD with level set (set-up of tutorial 11_1 scaled to the experiment). Surge-front "
                       "arrival at gauges Z = x/a = 1.5..10 as T = t sqrt(2g/a) against Martin & Moyce (1952). "
                       "The front has arrived when the surface is 0.05 a above the bed.",
        "reference": ref("martin"),
        "tags": ["cfd", "experiment", "2dv", "two-phase", "quick"],
        "np": 2,
        "time": 0.5,
        "check": {"type": "front_arrival", "a": a, "x_wall": 0.0, "z_bed": 0.0, "h_thr": 0.05 * a, "datum": 0.0,
                  "data": "dambreak/martin_moyce_1952_surge_n2_2_a2.25in_pysph.dat", "tol_max": 0.25,
                  "tol_mean": 0.08,
                  "note": "the largest deviation is at Z = 1.5 (front 20 % early): the experiment includes the "
                          "opening of the column, the simulation starts from rest without it"},
        "levels": {
            "nightly": {"np": 1, "control_set": {"B 1": "%g" % (a / 20), "B 10": "0.0 %.5f 0.0 %g 0.0 %.5f" % (L, a / 20, 3 * a),
                                                 "M 10": "1"},
                        "ctrl_set": {"M 10": "1", "P 51": ["%.5f %g" % (Z * a, a / 40) for Z in Zs]},
                        "check_set": {"tol_max": 0.25, "tol_mean": 0.1}},
            "release": {}},
    }
    write_case("cfd_dambreak_2d_martin_moyce", control, ctrl, case)


def cfd_kleefsman():
    H = 0.55
    tf = 1.0 / math.sqrt(9.81 / H)
    pf = 1000.0 * 9.81 * H
    control = """
C 11 21
C 12 21
C 13 21
C 14 21
C 15 21
C 16 3

B 1 0.01
B 10 0.0 3.22 0.0 1.0 0.0 1.0
S 10 0.67 0.83 0.3 0.7 0.0 0.16

M 10 16
M 20 2
"""
    ctrl = """
D 10 4
D 20 2
D 30 1
F 30 3
F 40 3
F 51 2.0
F 56 0.55
N 40 3
N 41 4.0
N 47 0.2
M 10 16
P 12 100
T 10 0
W 1 1000.0
W 22 -9.81
P 64 0.835 0.5 0.021
P 64 0.835 0.5 0.101
P 51 2.725 0.5
P 51 1.0 0.5
"""
    sig = lambda n, f: {"name": "P%d" % n, "source": "pressure", "n": 1 if n == 1 else 2,
                        "data": "dambreak/kleefsman_2005_P%d_pysph.dat" % n,
                        "data_t_factor": tf, "data_y_factor": pf}
    case = {
        "description": "MARIN 3D dam break with a box obstacle (Kleefsman et al. 2005, SPHERIC Test 2): tank "
                       "3.22 x 1 x 1 m, water column 1.22 x 1 x 0.55 m behind the door at x = 2 m, box 0.16 x 0.40 "
                       "x 0.16 m at x = 0.67..0.83 m. Pressure at P1 (z = 0.021 m) and P3 (z = 0.101 m) on the "
                       "water-facing box face against the measurements: impact peak value and time, rms error. "
                       "Sensor heights follow PySPH/LS-DYNA (0.021/0.101 m; the ComFLOW page gives 0.025/0.099 m).",
        "reference": ref("kleefsman"),
        "tags": ["cfd", "experiment", "3d", "two-phase", "impact"],
        "np": 16,
        "time": 4.0,
        "check": {"type": "timeseries", "signals": [sig(1, 1), sig(3, 2)], "align": [], "lag": 0.0,
                  "tol_rms": 0.35, "tol_height": None, "rms_about_mean": False,
                  "peak": {"signals": ["P1"], "window": [0.3, 0.8], "tol_value": 0.30, "tol_time": 0.03}},
        "levels": {
            "nightly": {"np": 4, "time": 2.0, "control_set": {"B 1": "0.025", "M 10": "4"},
                        "ctrl_set": {"M 10": "4", "P 64": ["0.85 0.5 0.021", "0.85 0.5 0.101"]},
                        "check_set": {"tol_rms": 0.5, "window": [0.0, 2.0],
                                      "peak": {"signals": ["P1"], "window": [0.3, 0.8], "tol_value": 0.3,
                                               "tol_time": 0.1}}},
            "release": {}},
    }
    write_case("cfd_dambreak_3d_kleefsman", control, ctrl, case)


def cfd_cylinder():
    """laminar 2D cylinder Re = 100, D = 0.1 m, U = 1 m/s, nu = 1e-3; domain 30D x 16D"""
    D, U, nu = 0.1, 1.0, 1.0e-3
    Hz = 1.6
    control = """
C 11 1
C 12 3
C 13 3
C 14 2
C 15 3
C 16 3

B 1 0.002
B 10 0.0 3.0 0.0 0.002 0.0 1.6
B 101 11
B 127 0.002 0.03 0.9 0.4 1.05
B 103 11
B 129 0.002 0.03 0.79 0.3 1.05

S 32 0.8 0.79 0.05

M 10 4
M 20 2
"""
    ctrl = """
B 10 0
B 20 2
B 60 1
B 61 1
D 10 4
D 20 2
D 22 2
D 30 1
F 30 0
F 40 0
I 11 1
N 40 3
N 41 50.0
N 47 0.3
M 10 4
P 12 100
P 81 0.7 0.9 -1.0 1.0 0.69 0.89
T 10 0
W 1 1000.0
W 2 %(nu)g
W 10 %(Q)g
"""
    case = {
        "description": "Laminar flow past a circular cylinder at Re = U D / nu = 100 (D = 0.1 m, U = 1 m/s, nu = "
                       "1e-3 m2/s), single-phase CFD in 2D, domain 30 D x 16 D (8 D upstream, blockage 6 %, cylinder 0.1 D below "
                       "the centre line so that the shedding starts), "
                       "cell-size stretching around the cylinder (D/50 release, D/25 nightly), implicit diffusion with the wall ghost values at the solid (D 22 2: with the default D 22 1 the cylinder wall is free of shear in the implicit diffusion, C_D drops by a third and the wake does not shed). Strouhal number from "
                       "the lift and mean drag coefficient after the shedding has settled.",
        "reference": ref("cylinder"),
        "tags": ["cfd", "analytical", "2d", "single-phase", "quick"],
        "np": 4,
        "time": 50.0,
        "check": {"type": "vortex_shedding", "U": U, "D": D, "W": 0.002, "nu": nu, "rho": 1000.0, "lift": "Fz",
                  "Cd_ref": 1.33, "t_start": 30.0, "t_end": 50.0, "tol_St": 0.05, "tol_Cd": 0.10,
                  "note": "blockage D/H = 1/16 raises St and C_D by a few per cent against the unconfined values"},
        "levels": {
            "nightly": {"np": 2, "time": 40.0,
                        "control_set": {"B 1": "0.004", "B 10": "0.0 3.0 0.0 0.004 0.0 1.6",
                                        "B 127": "0.004 0.04 0.9 0.4 1.08", "B 129": "0.004 0.04 0.79 0.3 1.08",
                                        "M 10": "2"},
                        "ctrl_set": {"M 10": "2", "W 10": "%g" % (U * Hz * 0.004)},
                        "check": {"type": "vortex_shedding", "U": U, "D": D, "W": 0.004, "nu": nu, "rho": 1000.0,
                                  "lift": "Fz", "Cd_ref": 1.33, "t_start": 25.0, "t_end": 40.0, "tol_St": 0.08,
                                  "tol_Cd": 0.15},
                        "xfail": "hans_dev ff1bb68, D/25: St = 0.1745 (+6 %) and the lift amplitude C_L = 0.31 "
                                 "(rms 0.22, literature 0.225-0.235) agree, but the P 81 drag gives C_D = 1.01 "
                                 "instead of 1.33 (-24 %)"},
            "release": {}},
    }
    ctrl = ctrl % {"nu": nu, "Q": U * Hz * 0.002}
    write_case("cfd_cylinder_re100", control, ctrl, case)


def cfd_beji():
    control = """
C 11 6
C 12 3
C 13 3
C 14 7
C 15 21
C 16 3

B 1 0.01
B 10 0.0 30.0 0.0 0.01 0.0 0.6

S 61 6.0 12.0 0.0 0.01 0.0 0.3
S 10 12.0 14.0 0.0 0.01 0.0 0.3
S 61 14.0 17.0 0.0 0.01 0.3 0.0

M 10 8
M 20 2
"""
    ctrl = """
B 10 1
B 50 0.0001
B 90 1
B 92 2
B 93 0.02 2.02
B 96 3.73 5.0
B 98 2
B 99 2
D 10 4
D 20 2
D 30 1
F 30 3
F 40 3
F 42 0.6
F 60 0.4
I 12 1
N 40 3
N 41 42.0
N 47 0.2
M 10 8
P 12 100
W 22 -9.81
%s
""" % gauges([x for _, x in BB_GAUGES], 0.005)
    case = {
        "description": "Regular waves (H = 2 cm, T = 2.02 s) over the Beji & Battjes submerged bar, two-phase CFD in "
                       "2D (tutorial 11_10 with the flume extended to 30 m so that all gauges lie outside the "
                       "numerical beach). Surface elevation at WG4..WG11 against the measured series.",
        "reference": ref("beji"),
        "tags": ["cfd", "experiment", "2dv", "two-phase", "dispersion"],
        "np": 8,
        "time": 42.0,
        "check": bb_check(0.35, 0.2),
        "levels": {
            "nightly": {"np": 4, "control_set": {"B 1": "0.02", "B 10": "0.0 30.0 0.0 0.02 0.0 0.6",
                                                 "S 61": ["6.0 12.0 0.0 0.02 0.0 0.3", "14.0 17.0 0.0 0.02 0.3 0.0"],
                                                 "S 10": "12.0 14.0 0.0 0.02 0.0 0.3", "M 10": "4"},
                        "ctrl_set": {"M 10": "4", "P 51": ["%.4f 0.01" % x for _, x in BB_GAUGES]},
                        "check_set": {"tol_rms": 0.4, "tol_height": 0.2}},
            "release": {}},
    }
    write_case("cfd_beji_battjes_bar", control, ctrl, case)


def cfd_ting_kirby():
    # tutorial 11_11: slope from x = 13.8 m, depth 0.38 m at x = 14.5 m
    x0 = 13.8 + 0.02 * 35.0
    xs = frange(18.0, 25.0, 0.1)
    control = """
C 11 6
C 12 3
C 13 3
C 14 7
C 15 21
C 16 3

B 1 0.005
B 10 0.0 40.0 0.0 0.005 0.0 1.00
S 61 13.8 40.0 0.0 0.005 0.0 0.748

M 10 32
M 20 2
"""
    ctrl = """
B 10 1
B 11 1
B 50 0.0001
B 90 1
B 92 8
B 93 0.128 5.0
B 96 9.8 0.0
B 98 2
B 99 2
D 10 4
D 20 2
D 30 1
F 30 3
F 40 3
F 42 1.0
F 60 0.4
I 12 1
N 40 3
N 41 60.0
N 47 0.1
M 10 32
P 12 100
T 10 2
T 36 2
W 10 0.0
W 22 -9.81
%s
""" % gauges(xs, 0.0025)
    case = {
        "description": "Plunging breaker of Ting & Kirby (1995), cnoidal waves H = 0.128 m, T = 5 s on a 1:35 "
                       "slope, two-phase CFD with k-omega (tutorial 11_11). Breaking point = location of the "
                       "maximum wave height from gauges every 0.1 m between x = 18 and 25 m, against the measured "
                       "x_b = 7.795 m from the d = 0.38 m station (here x = %.2f m)." % x0,
        "reference": ref("ting"),
        "tags": ["cfd", "experiment", "2dv", "two-phase", "breaking"],
        "np": 32,
        "time": 60.0,
        "check": {"type": "breaking_point", "x_b": 7.795, "x_offset": x0, "t_start": 40.0, "t_end": 60.0,
                  "tol": 0.4},
        "levels": {
            "nightly": {"np": 4, "time": 45.0,
                        "control_set": {"B 1": "0.02", "B 10": "0.0 40.0 0.0 0.02 0.0 1.00",
                                        "S 61": "13.8 40.0 0.0 0.02 0.0 0.748", "M 10": "4"},
                        "ctrl_set": {"M 10": "4", "P 51": ["%.4f 0.01" % x for x in xs]},
                        "check_set": {"t_start": 30.0, "t_end": 45.0, "tol": 1.0}},
            "release": {}},
    }
    write_case("cfd_ting_kirby_plunging", control, ctrl, case)


def cfd_chen():
    t = os.path.join(TUT, "REEF3D_CFD", "11_14 Non-Breaking Wave Forces")
    control = """
C 11 6
C 12 21
C 13 3
C 14 7
C 15 21
C 16 3

B 1 0.025
B 10 0.0 18.0 1.5 3.0 0.0 1.0

S 33 7.50 1.50 0.125

M 10 16
M 20 2
"""
    ctrl = """
B 10 1
B 11 1
B 50 0.0001
B 90 1
B 92 4
B 93 0.07 1.22
B 96 2.11 4.22
B 98 2
B 99 2
D 10 4
D 20 2
D 30 1
F 30 3
F 40 3
F 42 1.0
F 60 0.505
I 12 1
N 40 3
N 41 33.0
N 47 0.1
M 10 16
P 12 100
T 10 2
T 36 2
W 22 -9.81
P 51 5.00 1.5125
P 81 7.325 7.675 1.4 1.675 0.0 1.0
"""
    case = {
        "description": "Regular waves (2nd-order Stokes, H = 0.07 m, T = 1.22 s, d = 0.505 m) on a vertical "
                       "cylinder D = 0.25 m (Chen et al. 2014), two-phase CFD with k-omega, tutorial 11_14 on the "
                       "half domain y = 1.5..3 m (symmetry plane through the cylinder axis, force x 2). Surface "
                       "elevation at x = 5 m and inline force against the measurements; the time lag is fitted on "
                       "the gauge.",
        "reference": ref("chen"),
        "tags": ["cfd", "experiment", "3d", "two-phase", "wave-force"],
        "np": 16,
        "time": 33.0,
        "check": {"type": "timeseries",
                  "signals": [{"name": "eta x=5 m", "source": "gauge", "index": 1,
                               "data": "wave_forces/chen2014_eta_x5.txt"},
                              {"name": "Fx", "source": "force", "n": 1, "component": "Fx", "factor": 2.0,
                               "data": "wave_forces/chen2014_force.txt"}],
                  "align": ["eta x=5 m"], "lag_range": [-1.3, 1.3], "lag_step": 0.005,
                  "tol_rms": 0.25, "tol_height": 0.10},
        "levels": {
            "nightly": {"np": 4, "control_set": {"B 1": "0.05", "M 10": "4"}, "ctrl_set": {"M 10": "4",
                        "P 51": "5.0 1.525"},
                        "check_set": {"tol_rms": 0.4, "tol_height": 0.2}},
            "release": {}},
    }
    write_case("cfd_wave_force_chen2014", control, ctrl, case)


def cfd_irschik():
    control = """
C 11 6
C 12 21
C 13 3
C 14 8
C 15 21
C 16 3

B 1 0.05
B 10 0.0 54.0 2.5 5.0 0.0 7.0

S 33 44.0 2.50 0.35
S 61 21.0 44.0 2.5 5.0 0.0 2.30
S 10 44.0 54.0 2.5 5.0 0.0 2.30

M 10 64
M 20 2
"""
    ctrl = """
B 10 1
B 50 0.0001
B 90 1
B 92 5
B 93 1.3 4.0
B 96 21.0 0.0
B 98 2
B 99 3
D 10 4
D 20 2
D 30 1
F 30 3
F 40 3
F 42 7.0
F 60 3.80
I 12 1
N 40 3
N 41 27.0
N 47 0.1
M 10 64
P 12 100
T 10 2
T 36 2
W 22 -9.81
P 51 43.65 2.55
P 81 43.30 44.70 2.5 3.20 2.30 8.0
"""
    case = {
        "description": "Breaking wave impact on a vertical pile D = 0.7 m in the GWK (Irschik et al. 2002): 5th-order "
                       "Stokes waves H = 1.3 m, T = 4 s on d = 3.8 m, 1:10 slope, the wave breaks at the pile "
                       "(tutorial 11_15, half domain with a symmetry plane, force x 2). Release level only: the "
                       "slamming force (peak and rms over the measured impact) against the measured force; the "
                       "lag is fitted on the force itself because the position of the measured gauge in the tutorial "
                       "data is not documented.",
        "reference": ref("irschik"),
        "tags": ["cfd", "experiment", "3d", "two-phase", "wave-force", "breaking"],
        "np": 64,
        "time": 27.0,
        "check": {"type": "timeseries",
                  "signals": [{"name": "Fx", "source": "force", "n": 1, "component": "Fx", "factor": 2.0,
                               "data": "wave_forces/irschik2002_force.txt"}],
                  "align": ["Fx"], "lag_range": [-4.0, 4.0], "lag_step": 0.005,
                  "tol_rms": 0.4, "tol_height": 0.25},
        "levels": {"release": {}},
    }
    write_case("cfd_breaking_wave_force_irschik2002", control, ctrl, case)


def hulme_values():
    import sys
    sys.path.insert(0, ROOT)
    import benchmark
    T, zeta, kR = benchmark.hulme_prediction(benchmark.read_refdata("misc/hulme_hemisphere_heave.txt"), 1.0, 0.5)
    return T, zeta


def cfd_sphere():
    t = os.path.join(TUT, "REEF3D_CFD", "11_16 Heave Decay of Sphere")
    control = """
C 11 21
C 12 21
C 13 21
C 14 21
C 15 21
C 16 3

B 1 0.025
B 10 0.0 6.0 0.0 6.0 0.0 4.0

B 101 11
B 127 0.025 0.2 3.0 1.5 1.08
B 102 11
B 128 0.025 0.2 3.0 1.5 1.08
B 103 11
B 129 0.025 0.2 2.0 1.5 1.08

M 10 32
M 20 2
"""
    ctrl = """
B 10 1
B 50 0.00001
B 90 1
B 99 2
B 107 0.0 0.0 0.0 6.0 0.5
B 107 6.0 6.0 0.0 6.0 0.5
B 107 0.0 6.0 0.0 0.0 0.5
B 107 0.0 6.0 6.0 6.0 0.5
D 10 4
D 20 2
D 30 1
F 30 3
F 35 5
F 40 3
F 60 2.0
I 10 1
N 40 4
N 41 8.0
N 47 0.3
M 10 32
P 12 100
W 1 1000.0
W 22 -9.81
X 10 1
X 11 0 0 1 0 0 0
X 21 500.0
X 180 1
X 182 3.0 3.0 2.1
"""
    case = {
        "description": "Free heave decay of a sphere R = 0.5 m of half the water density (tutorial 11_16, heave "
                       "only, initial offset 0.1 m = 0.2 R, damping zones at the walls). Equilibrium from Archimedes "
                       "(centre at the still water level), damped period and damping ratio of the first periods "
                       "against linear potential theory with Hulme's (1982) hemisphere coefficients "
                       "(T = %.3f s, zeta = %.3f). The offset is moderate, so a few per cent nonlinear deviation is "
                       "expected and allowed in the tolerances." % hulme_values(),
        "reference": ref("hulme"),
        "tags": ["cfd", "semi-analytical", "3d", "two-phase", "6dof"],
        "np": 32,
        "time": 8.0,
        "files": ["floating.stl"],
        "check": {"type": "heave_decay", "R": 0.5, "mass_ratio": 1.0, "z_eq": 2.0, "n_periods": 2,
                  "hulme": "misc/hulme_hemisphere_heave.txt", "tol_T": 0.05, "tol_zeta": 0.25, "tol_eq": 0.02,
                  "eq_window": 1.0},
        "levels": {
            "nightly": {"np": 4, "control_set": {"B 1": "0.05", "B 127": "0.05 0.3 3.0 1.5 1.1",
                                                 "B 128": "0.05 0.3 3.0 1.5 1.1", "B 129": "0.05 0.3 2.0 1.5 1.1",
                                                 "M 10": "4"},
                        "ctrl_set": {"M 10": "4"},
                        "check_set": {"tol_T": 0.06, "tol_zeta": 0.25, "tol_eq": 0.02}},
            "release": {}},
    }
    write_case("cfd_sphere_heave_decay", control, ctrl, case, {os.path.join(t, "floating.stl"): "floating.stl"})


def cfd_pier_scour():
    D, h, U = 0.2, 0.3, 0.3
    Fr = U / math.sqrt(9.81 * h)
    hec = 2.0 * 1.0 * 1.0 * 1.1 * (h / D) ** 0.35 * Fr ** 0.43
    control = """
C 11 1
C 12 21
C 13 21
C 14 2
C 15 21
C 16 3

B 1 0.02
B 10 1.0 5.5 0.0 2.0 -0.3 0.5

S 33 2.5 1.0 0.1

M 10 32
M 13 0
"""
    ctrl = """
B 10 1
B 11 1
B 50 0.006
B 60 1
B 61 2
D 10 4
D 20 2
D 30 1
F 30 3
F 40 3
F 60 0.3
I 10 1
N 40 2
N 47 0.3
M 10 32
P 12 100
P 122 1
T 10 2
T 36 1
W 10 0.18
W 22 -9.81
S 10 1
S 11 3
S 12 0
S 13 5.0
S 14 0.3
S 16 4
S 19 75600
S 20 0.00097
S 21 5.0
S 22 2570.0
S 24 0.75
S 30 0.035
S 41 1
S 42 1
S 43 100
S 44 1
S 50 1
S 57 0.0
S 60 0
S 73 0.00 0.5 0.0 1.0 0.0
S 80 4
S 81 35.0
S 82 10.0
S 90 1
S 91 10
"""
    case = {
        "description": "Local scour at a vertical circular pier D = 0.2 m in a steady current (h = 0.3 m, U = 0.3 "
                       "m/s, d50 = 0.97 mm), CFD with k-omega and the Exner model (tutorial 11_18). Release: "
                       "21 h of sediment time (S 19 75600); the equilibrium scour depth S/D must lie in the range "
                       "of established estimates: HEC-18/CSU %.2f (K1 = 1, K2 = 1, K3 = 1.1), Sumer et al. 1.3 "
                       "(sigma 0.7). Nightly: 1 h of sediment time on a coarse grid as a smoke test (scour must "
                       "start, S/D between 0.05 and 2)." % hec,
        "reference": ref("scour"),
        "tags": ["cfd", "empirical", "3d", "sediment"],
        "np": 32,
        "check": {"type": "scour_depth", "D": D, "band": [0.9, 2.0], "hec18": round(hec, 3), "bed0": 0.0},
        "levels": {
            "nightly": {"np": 4, "control_set": {"B 1": "0.04", "M 10": "4"},
                        "ctrl_set": {"M 10": "4", "S 19": "3600", "S 17": "5.0"},
                        "check_set": {"band": [0.05, 2.0]}},
            "release": {}},
    }
    write_case("cfd_pier_scour", control, ctrl, case)


def main():
    if os.path.isdir(CASES):
        shutil.rmtree(CASES)
    sflow_dambreak("sflow_dambreak_ritter", 0.0, "ritter")
    sflow_dambreak("sflow_dambreak_stoker", 0.05, "stoker")
    sflow_solitary()
    sflow_beji()
    nhflow_beji()
    nhflow_stokes5()
    nhflow_synolakis()
    berkhoff("nhflow")
    berkhoff("fnpf")
    fnpf_beji()
    fnpf_stokes2()
    fnpf_ting_kirby()
    cfd_martin_moyce()
    cfd_kleefsman()
    cfd_cylinder()
    cfd_beji()
    cfd_ting_kirby()
    cfd_chen()
    cfd_irschik()
    cfd_sphere()
    cfd_pier_scour()
    print("\n".join(sorted(os.listdir(CASES))))


if __name__ == "__main__":
    main()
