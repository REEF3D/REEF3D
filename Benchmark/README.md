# REEF3D literature benchmark suite

Established test cases from the literature for **REEF3D::CFD** (incl. 6DOF and sediment),
**REEF3D::NHFLOW**, **REEF3D::FNPF** and **REEF3D::SFLOW**: laboratory experiments, exact and
semi-analytical solutions. Each case runs REEF3D, compares the output with the published reference
data and gives PASS/FAIL against a tolerance.

The suite complements the regression suite in `Regression/`:

| set | question | cases |
|---|---|---|
| `Regression/` | does the code give the **same** numbers as before? (bitwise / tolerance) | short runs, stored references |
| **`Benchmark/` (this)** | does REEF3D reproduce the **established literature benchmarks**? | 33 cases, experiments + exact solutions |

It uses the same case format (`case.json` + `control.txt` + `ctrl.txt`, `"base"` + overrides) and
is self-contained (Python 3 standard library; matplotlib only for plots; it does not need the
regression dump). It lives in the repository as `Benchmark/` next to `Regression/`; reports and
plots of verification runs are kept out of the repository (`~/Dropbox/Claude/REEF3D_benchmark_suite/results/`).

## Files

| path | content |
|---|---|
| `benchmark.py` | runner, checkers, report, plots |
| `cases/<name>/` | `case.json` (description, reference, check, levels), `control.txt`, `ctrl.txt`, extra files |
| `refdata/` | reference data with provenance (`*/PROVENANCE.md`, `PROVENANCE_git_sources.md`) |
| `tools/make_cases.py` | writes `cases/` (documents the gauge lists and geometry; the case files are the source of truth) |
| `tools/berkhoff_geo.py`, `conical_geo.py`, `thacker_geo.py` | `geo.dat` bathymetries (run automatically by `write`/`run`, `"generate"` in case.json) |
| `tools/mase_kirby_waverecon.py` | `waverecon.dat` wave components from the measured Mase & Kirby record |
| `tools/run_suite.sh` | run one level + plots |

## Two levels

Every case has two levels, selected with `--level`:

| level | purpose | grids | cost | CI tier (REEF3D CI plan) |
|---|---|---|---|---|
| `nightly` | catch degradations of accuracy | coarse (about half the release resolution, 2D cases 1-4 ranks) | seconds to ~1 h per case | Tier 2 (nightly, NTNU workstation) |
| `release` | reproduce the published agreement | the resolution of the papers / tutorials | minutes to many hours, up to 64 ranks | Tier 3 (release, Idun) |

The levels are overrides in `case.json` (`"levels": {"nightly": {...}, "release": {...}}`) of the
grid (`control_set`), the inputs (`ctrl_set`), the simulated time (`time`), the ranks (`np`) and the
tolerances (`check_set`). A case without a `nightly` entry is release-only.

## Quick start

```bash
cd Benchmark                  # or ~/Dropbox/Claude/REEF3D_benchmark_suite
./benchmark.py list --level nightly

export REEF3D_MPIRUN="mpirun"  # "mpirun --oversubscribe" with fewer cores than ranks
./benchmark.py run --level nightly --reef3d ~/Codelite/REEF3D/bin/REEF3D \
    --divemesh ~/Codelite/DIVEMesh/bin/DiveMESH --out ~/reef3d_bench/nightly --check

./benchmark.py check ~/reef3d_bench/nightly        # evaluate again (e.g. after changing a tolerance)
./benchmark.py plot  ~/reef3d_bench/nightly        # one PNG per case: simulation vs. reference
```

Options: `--cases 'nhflow_*'`, `--tags quick` (fast subset), `--max-np 8` (cap ranks on a small
machine), `--keep` (keep grids and VTU/VTP), `--timeout` (s per case, default 6 h).
`./benchmark.py write --out DIR --level release` writes the resolved input files only, e.g. to
submit the release cases as cluster jobs, then `./benchmark.py check DIR` evaluates them (put a
`suite.json` with `{"level": "release"}` into DIR and a `run.json` with `{"status": "ok"}` into each
case folder, or run through `benchmark.py run`).

The report goes to `<out>/benchmark.md` (+ `benchmark.json`): one line per case with the error
measure and the tolerance, and per case the details (per gauge, per profile, ...).
Return code 0 if all cases pass (or fail as expected, see XFAIL).

## Cases

| case | model | reference | compared | nightly | release |
|---|---|---|---|---|---|
| `sflow_dambreak_ritter` | SFLOW | Ritter (1892), exact (SWASHES) | depth h(t) at 8 gauges, rel. L1 | dx 0.05, 1 rank, s | dx 0.01 |
| `sflow_dambreak_stoker` | SFLOW | Stoker (1957), exact (SWASHES) | depth h(t) at 8 gauges, rel. L1 | dx 0.05 | dx 0.01 |
| `sflow_solitary_celerity` | SFLOW | solitary wave c = sqrt(g(d+H)) | celerity, crest height after 16 m | dx 0.02 | dx 0.01 |
| `sflow_beji_battjes_bar` | SFLOW | Beji & Battjes (1993) | eta(t) at WG4-WG11 | dx 0.025 | dx 0.01 |
| `nhflow_stokes5_propagation` | NHFLOW | Fenton (1985) 5th-order Stokes | P 51 vs. theory P 50 over 125 m | 400 x 8 | 800 x 10 |
| `nhflow_beji_battjes_bar` | NHFLOW | Beji & Battjes (1993) | eta(t) at WG4-WG11 | 700 x 5 | 1500 x 10 |
| `nhflow_synolakis_runup_nonbreaking` | NHFLOW | Synolakis (1987) run-up law and profiles | max. run-up R/d and 5 measured profiles, H/d = 0.0185 | dx 0.1 | dx 0.025 |
| `nhflow_synolakis_breaking_profiles` | NHFLOW | Synolakis (1987) | 12 surface profiles, H/d = 0.28 | dx 0.1 | dx 0.025 |
| `nhflow_berkhoff_shoal` | NHFLOW | Berkhoff et al. (1982) | H/H0 on sections 1-8 (3D) | dx 0.1 | dx 0.04 |
| `fnpf_stokes2_propagation` | FNPF | 2nd-order Stokes | P 51 vs. theory P 50 | 400 x 8 | 800 x 10 |
| `fnpf_beji_battjes_bar` | FNPF | Beji & Battjes (1993) | eta(t) at WG4-WG11 | 800 x 8 | 2000 x 10 |
| `fnpf_ting_kirby_plunging` | FNPF | Ting & Kirby (1995) | breaking point (max. H) | dx 0.05 | dx 0.05, 100 s |
| `fnpf_berkhoff_shoal` | FNPF | Berkhoff et al. (1982) | H/H0 on sections 1-8 (3D) | dx 0.1 | dx 0.04 |
| `cfd_dambreak_2d_martin_moyce` | CFD | Martin & Moyce (1952) | surge front Z(T) | a/20 | a/40 |
| `cfd_dambreak_3d_kleefsman` | CFD | Kleefsman et al. (2005), SPHERIC 2 | pressures P1-P8, water heights at 4 gauges; P1 peak | dx 0.025, 2 s | dx 0.01, 6 s |
| `cfd_cylinder_re100` | CFD | Williamson St(Re); C_D compilation | St, mean C_D at Re = 100 | D/25 | D/50 |
| `cfd_beji_battjes_bar` | CFD | Beji & Battjes (1993) | eta(t) at WG4-WG11 | dx 0.02 | dx 0.01 |
| `cfd_ting_kirby_plunging` | CFD | Ting & Kirby (1995) | breaking point (max. H) | dx 0.02 | dx 0.005 |
| `cfd_wave_force_chen2014` | CFD | Chen et al. (2014) | eta(t), inline force on cylinder | dx 0.05 | dx 0.025 |
| `cfd_breaking_wave_force_irschik2002` | CFD | Irschik et al. (2002), GWK | slamming force | – | dx 0.05 |
| `cfd_sphere_heave_decay` | CFD 6DOF | Hulme (1982) linear theory, Archimedes | damped period, damping ratio, equilibrium | dx 0.05 | dx 0.025 |
| `cfd_pier_scour` | CFD sediment | HEC-18 (CSU), Sumer et al. (1992) | equilibrium scour depth S/D | smoke test, 1 h sediment time | dx 0.02, 21 h sediment time |
| `sflow_thacker_parabola` | SFLOW | Thacker (1981), exact (SWASHES) | depth h(t) at 9 gauges incl. wetting/drying, rel. L1 | dx 0.02 | dx 0.005 |
| `sflow_conical_island_A`, `_C` | SFLOW | Briggs et al. (1995) | gauges 6, 9, 16, 22; run-up on 12 radial lines | dx 0.1 | dx 0.05 |
| `nhflow_conical_island_C` | NHFLOW | Briggs et al. (1995) | gauges 6, 9, 16, 22; run-up on 12 radial lines | dx 0.1 x 5 | dx 0.05 x 5 |
| `sflow_mase_kirby_irregular` | SFLOW | Mase & Kirby (1992) | Hm0 and skewness at 12 gauges on a 1:20 slope | dx 0.04, 150 s | dx 0.02, 420 s |
| `nhflow_mase_kirby_irregular` | NHFLOW | Mase & Kirby (1992) | Hm0 and skewness at 12 gauges | 362 x 5, 200 s | 725 x 8, 420 s |
| `fnpf_mase_kirby_irregular` | FNPF | Mase & Kirby (1992) | Hm0 and skewness at 12 gauges | 362 x 8, 200 s | 725 x 10, 420 s |
| `fnpf_jonswap_spectrum` | FNPF | JONSWAP (DNV-RP-C205) | Hm0, Tp, spectral shape at 3 gauges | 400 x 8, 600 s | 800 x 10, 1100 s |
| `cfd_sloshing_linear` | CFD | linear sloshing theory | first-mode period | dx 0.02 | dx 0.01 |
| `nhflow_sloshing_linear` | NHFLOW | linear sloshing theory | first-mode period | 50 x 5 | 100 x 10 |
| `cfd_porous_dambreak_lin1998` | CFD VRANS | Lin (1998), Liu et al. (1999) | 5 free-surface profiles | dx 0.01 | dx 0.005 |

Each `case.json` holds the full description, the citation (`"reference"`), the check parameters and
the tolerances of both levels. The inputs are based on the user-guide tutorials where one exists;
the changes against the tutorial are given in the description (e.g. flume extended so that all
gauges are outside the numerical beach, half domain with a symmetry plane for the cylinder cases,
smaller initial offset of the sphere for linear theory).

## Checkers (`benchmark.py`)

| type | used for | error measure |
|---|---|---|
| `timeseries` | gauges, forces, pressure probes vs. measured series | rel. rms error rms(sim - data)/rms(data) and height ratio per signal; common time lag fitted when the measured time origin is arbitrary, plus an optional small per-signal lag (`local_lag`, gauge position / digitising uncertainty); optional peak value and time |
| `theory_gauges` | wave propagation | P 51 vs. P 50 (generating theory): rel. rms error and height ratio |
| `front_arrival` | dam-break surge front | arrival time T = t sqrt(2g/a) at Z = x/a vs. measured T(Z) |
| `vortex_shedding` | cylinder | Strouhal number (lift zero crossings), mean C_D |
| `heave_decay` | floating sphere | damped period, damping ratio (log. decrement), final equilibrium |
| `runup_law` | solitary run-up | R/d vs. 2.831 sqrt(cot b)(H/d)^1.25 |
| `profiles` | surface profiles at given times | rms(eta - eta_data)/H per profile, time origin fitted on the first profile |
| `section_heights` | 3D wave fields | H/H0 along measurement sections, rms error |
| `breaking_point` | breaking on slopes | location of max. wave height vs. measured x_b |
| `swe_dambreak` | shallow-water dam break | exact Ritter/Stoker depth, rel. L1 |
| `solitary` | solitary propagation | celerity and crest height |
| `scour_depth` | pier scour | S/D within the band of established estimates |
| `runup_angles` | run-up around an island | highest wetted gauge on radial lines vs. measured run-up, rms/d |
| `wave_stats` | irregular waves in the surf zone | Hm0 = 4 sigma and skewness vs. the measured records (same window) |
| `spectrum` | irregular waves | Hm0 = 4 sqrt(m0), Tp and S(f) shape vs. the target JONSWAP spectrum (FFT in pure Python) |
| `oscillation_period` | sloshing | period from up-crossings vs. linear theory |
| `thacker1d` | parabolic bowl | exact Thacker depth, rel. L1 |
| `profiles_abs` | free-surface profiles at physical times | rms(z - z_data)/h_ref per profile |

A case can have `"extra_checks"` (a list of further check dicts, e.g. time series plus run-up);
all of them must pass.

New case: copy a case folder, edit, `./benchmark.py run --cases <new> --level nightly`, look at the
plot, then set the tolerances. A case that deviates from the reference for a known reason gets
`"xfail": "<reason>"` (report shows XFAIL and the suite still passes; XPASS once it is fixed).

## Reference data

All reference data are in `refdata/`, each file with a header giving the source (URL, or git
repository + commit + path), the original reference, units and scaling.
`refdata/*/PROVENANCE.md` and `refdata/PROVENANCE_git_sources.md` summarise sources and confidence
(A = original instrument records, B = digitised or tabulated by a third party, C = unclear).
Values are copied unchanged; conversions (cm -> m, amplitude -> H/H0, non-dimensional scaling) are
done in the checkers and documented in `case.json`.

| data | source | confidence |
|---|---|---|
| SPHERIC Test 2 (Kleefsman et al. 2005): P1-P8, 4 water heights, dt 1 ms | AQUAgpusph repository | A |
| Briggs et al. (1995) conical island: gauges 6, 9, 16, 22 (cases A-C), dt 0.04 s | FUNWAVE-TVD benchmark | A |
| Briggs et al. (1995) run-up around the island (24 angles) | NHWAVE repository | A/B |
| Mase & Kirby (1992): 12 gauge records on a 1:20 slope, dt 0.05 s | FUNWAVE-TVD benchmark | A |
| Berkhoff et al. (1982) sections 1-8 (amplitudes) | FUNWAVE-TVD benchmark | B |
| Synolakis (1987) H/d = 0.0185 profiles | FUNWAVE-TVD (NTHMP BP4) | B |
| Lin (1998) porous dam break, 6 profiles | olaFlow tutorial CR35_dambreak | B (digitised) |
| Beji & Battjes WG4-WG11 eta(t) | Basilisk `test/bar.c` data files | B (digitised) |
| Berkhoff sections 2, 3, 5, 7 (first round, no longer used by a case) | Basilisk `examples/shoal.c` | B (digitised) |
| Synolakis H/d = 0.28 profiles | Basilisk `test/beach.c` data files | B (digitised) |
| Martin & Moyce surge front (a = 2.25 in) | PySPH `db_exp_data.py` | B (digitised) |
| Kleefsman P1, P3 (first round, no longer used by a case) | PySPH | B (digitised) |
| Chen et al. eta + force, Irschik et al. force | REEF3D tutorials 11_14, 11_15 | B (digitised) |
| Hulme hemisphere added mass / damping | LHEEA `Hulme_Heaving_Sphere` | tabulated |
| Ritter, Stoker, Thacker, solitary, Stokes, sloshing, JONSWAP, run-up law, St(Re), C_D | formulas / published values | exact / tabulated |

The public test suites of FUNWAVE-TVD, NHWAVE, AQUAgpusph and olaFlow were cloned with git (the web
pages of GitHub are blocked in the sandbox, `git clone` works); `refdata/PROVENANCE_git_sources.md`
lists the commits. Further measured data found there and not yet used by a case: OSU shelf/island
run-up (9 gauges, 3 ADVs), composite beach (NTHMP BP5), Monai valley (BP7), Synolakis H/d = 0.3
profiles, Ghia et al. (1982) cavity tables (REEF3D has no moving-lid wall to run it).

Because many series are digitised, the time-series comparisons allow a fitted time lag and use rms
errors of 0.1-0.6 as tolerances rather than tighter values.

## Status of the reference data (what is missing)

- **Beji & Battjes original records** (Dingemans 1994, Delft Hydraulics H1684, 0.05 s sampling): not
  in any repository found; the digitised Basilisk series are used.
- **Ting & Kirby (1995)** measured wave height and set-up along the flume: not found as data; the
  cases check only the breaking point.
- **Whalin (1971)** harmonic amplitudes: not found as data (OceanWave3D ships only the set-up).
- **Conical island gauge positions**: taken from the FUNWAVE-TVD grid indices (gauge radii 3.60, 2.45,
  2.55, 2.55 m from the island centre); the lab coordinates are not in the repositories, and a
  per-gauge lag of up to 0.2 s is allowed. The run-up angle convention (270 deg = facing the wave)
  is inferred from the data.
- **SPHERIC Test 2 height labels**: the four height columns follow AQUAgpusph's reading (column 10 =
  gauge nearest the reservoir); the cases name the gauges by position.

## Findings on hans_dev while setting up the suite

Observed while building and running the cases (hans_dev `ff1bb6819`, later `4bb681b`, Linux,
MPICH 4.2.3, REEFMG, `make all`). The fixes are delivered as patches in
`~/Dropbox/Claude/REEF3D_benchmark_suite/patches/` (0001-0005 on `4bb681b`, applied in `d80a861`;
0006-0007 on `d80a861`).

1. (Fixed on hans_dev in the meantime: `src/momentum_FCC3.cpp`, which still included the deleted
   `momentum_FCC3.h` and broke `make`, has been removed after `ff1bb68`.)
2. **CFD implicit diffusion at solids (default `D 22 1`)**: the laminar cylinder at Re = 100 with
   `D 20 2` and the default `D 22 1` gave C_D = 0.62 and no vortex shedding after 40 s (400 D/U);
   explicit diffusion (`D 20 1`) or `D 22 2` give the expected behaviour. **Cause:** in
   `idiff2_FS[uvw]*.cpp` a fluid face next to a direct-forcing solid face (`DF < 0`) got a Neumann
   coupling (`M.p += M.s`), i.e. a free-slip wall in the implicit diffusion. **Fixed in patch 0003**
   (Dirichlet u = 0 at the forced solid face; domain walls unchanged): the cylinder sheds with
   `D 22 1` (St 0.172, C_L 0.30, C_D 1.19 at 8-15 s; `D 22 2`: St 0.175, C_D 1.33). The remaining
   10 % in C_D is the staircase position of the wall (one face instead of the level-set distance);
   `D 22 2` (already forced for `X 10` bodies) is the more accurate choice and a candidate default.
3. **`B 20 3` with a DIVEMesh solid**: with `B 20 3` the flow did not see the solid (cylinder force
   ~0, u_max 1.15 instead of 1.28 m/s). **Cause:** `ghostcell::solid_forcing`
   (`gc_solid_forcing.cpp`) built the direct forcing of solids only for `B 20 1` and `B 20 2` (with
   `B 21 1`); for `B 20 3/4` no branch ran and the forcing stayed zero. `B 20` only selects the
   ghost-cell condition at domain walls and topo, and the `B 20 1` block was an identical copy of
   the `B 20 2` block. **Fixed in patch 0007** (forcing depends on `B 21` only): `B 20 3` gives
   C_D = 1.165 at t = 1-2 s, the same as `B 20 2`.
4. **SFLOW `A 243 1`**: cells that are dry at the start never become wet, so a dam break onto a
   dry bed does not move at all (u_max ~ 1e-11). Default `A 243 2` works. By design option 1 is a
   pure depth criterion ("depth criterion"; option 2 adds "wetting from higher wet neighbours"):
   continuity and the HLL fluxes skip dry cells, so a dry cell can never gain water. Not changed;
   option 1 is only suitable for drying (worth a note in the user guide, or an input warning).
5. **3D dam break with box (Kleefsman set-up, dx = 0.025 m, N 40 3, level set)**: the water volume
   drops from 0.67 to 0.55 m3 (-18 %) until t = 1.7 s and the run stops at t = 1.77 s with the N 61
   velocity limit (u_max 365 m/s). Up to t = 1.7 s the front pressures agree (P1 peak 10.6 vs
   11.2 kPa, 0.05 s late; rms P1-P3 0.41-0.45, P4 0.64), the top probes P5-P8 read close to zero
   (rms 0.83-0.95), heights at x = 1.456 / 0.960 m rms 0.17 / 0.12, behind the box (x = 0.464 m) the
   water arrives late (rms 0.46). **Cause:** the level set is not mass conserving (coarse grid
   dx = 0.05: -33 % with `N 40 3`, -26 % with `N 40 13`, so not the per-stage reinitialisation of the
   RK scheme), and the volume correction `F 46` was broken with `N 40 2/3/4`: `momentum_rk` never
   called `volcalc`, and `picard_lsm` took its reference at `count == 1` while the scheme is already
   called at `count 0`, so `F 46 3` removed all the water in the first step. **Fixed in patch 0004.**
   With the patch and `F 46 3` (`F 47 10`): volume 0.673 -> 0.668 m3 (-0.7 %), the run passes
   t = 1.77 s (stopped at 1.85 s for time, u_max 30-50 m/s in the air after the impact), and up to
   1.8 s all signals are within the nightly tolerances: P1 peak 10.8 vs 11.2 kPa (0.05 s late),
   rms P1-P3 0.38-0.44, P4 0.59, P5-P8 0.63-0.74, heights at x = 2.606 / 1.456 / 0.960 / 0.464 m
   0.02 / 0.23 / 0.20 / 0.39. The case
   itself keeps the defaults (no `F 46`) until the patch is applied; then `F 46 3` is recommended.
6. **FNPF plunging breaker (tutorial 9_3 set-up)**: at the tutorial resolution (dx = 0.05 m) the
   maximum wave height is reached 1.9 m before the measured breaking point; at dx = 0.025 m the
   wave train becomes irregular at t ~ 75 s and the run stops at t = 91 s (N 61).
7. **FNPF regular waves**: 2nd-order Stokes waves in deep water run about 0.4 % faster than the
   theory (phase lead 0.017 s at x = 15 m, 0.032 s at x = 25 m), independent of the grid level.
8. **P 81 drag of the cylinder (Re = 100, nightly grid D/25, `D 22 2`)**: shedding frequency
   (St = 0.1745, +6 % with 6 % blockage) and lift amplitude (C_L = 0.31, i.e. rms 0.22 vs. 0.225-0.235
   in the literature) are right, but the mean drag was C_D = 1.01 instead of 1.33. **Cause:** the
   viscous force in `force_force.cpp` was `mu*A*(du*ny+du*nz)`, the wall shear multiplied by signed
   normal components, so the friction drag cancelled between the upper and lower half of the body
   (the same in `6DOF_obj_forces_lsm.cpp`). **Fixed in patch 0002** (tau = mu u_t/dn with the
   tangential velocity, relative to the body for 6DOF): C_D = 1.33 at 4bb681b + patch, lift unchanged;
   6DOF sphere heave decay unchanged (T and zeta within 0.1 %).
9. **`B 105` origin inside the domain mirrors the waves**: with irregular-wave reconstruction
   (`B 92 51`) and `B 105` at the toe gauge (x0 = 2.75 m) the generated waves were 13x too small
   (Hm0 15 % of measured at 4bb681b). **Cause:** `xgen`/`ygen` in `iowave_dist.cpp` (the coordinates
   the wave theories are evaluated at) were unsigned distances to the `B 105` line, so the wave
   field was mirrored there and the generation zone upstream of it got reversed phases (the zone
   itself is not moved, as first assumed). `ygen` was also not perpendicular to `xgen` for oblique
   `B 105_1`. **Fixed in patch 0001** (signed wave-frame coordinates): the Mase & Kirby case with
   `B 105 0 2.75 0` then gives the same result as with the default origin; regression case
   `nhflow_2d_origin` (origin at x = 10 m) run to 60 s: wave height 0.992 of theory (0.964 before).
10. **`B 92 20` (wave maker from a measured time series) in SFLOW**: the SFLOW wave heights came out
   1.7x too large for the Mase & Kirby record. **Cause** (code reading): `wave_lib_piston_eta`
   derives the velocity with the shallow-water transfer u = eta sqrt(g/h), uniform in depth; for the
   intermediate-depth spectrum (kh = 1.4-5) that overestimates the depth-averaged velocity by
   sqrt(kh/tanh kh) = 1.26 at the peak and 1.7-2.3 at 1.5-2 fp. The inflow then acts as a reflecting
   velocity boundary (the ghost eta is overwritten by the Neumann update), so outgoing long waves are
   re-reflected. Fix proposal (not implemented, it belongs to the iowave redesign): a linear transfer
   per frequency component (FFT of the record, u = eta omega/(k h)) and a Riemann/Flather-type
   generating-absorbing inflow; also `wave_fi` of that class returns an uninitialised value. The cases
   use `B 92 51` with relaxation generation instead.
11. **CFD initial free surface `F 57` / `F 70`**: a sloping plane (`F 57`) or a cosine built from
   `F 70` strips is resolved only to the cell (staircase): `ini_phi.cpp` sets the level set to an
   inside/outside value per cell, so the interface lies on cell faces. That gives a 60 % error in
   the linear sloshing period check; the exact linear surface with `F 60` / `F 62` / `F 63` gives
   +0.4 %. Not changed (a signed-distance initialisation would change all dam-break cases; for
   grid-aligned boxes it makes no difference).
12. **SFLOW wave gauges read the neighbouring cell**: in the Thacker parabola the gauge level was off
   by about one cell of bed rise near the shoreline. **Cause:** `sflow_print_wsf::ini_location`
   rounded (x - x0)/dx to the nearest integer, so a gauge in the upper half of cell i read cell
   i+1 (and a uniform dx was assumed). **Fixed in patch 0006** (the containing cell, `posc_i/j`, as
   the other SFLOW probes): Thacker error 0.065 -> 0.051 (gauges x = 0.75 / 3.25 m: 0.139 -> 0.093,
   0.499 -> 0.100). A dry gauge reports bed + `A 244` (relative to `F 60`). The checkers still treat
   a gauge as wet only when its level rises above its initial value (conical-island run-up) and mask
   dry gauges (Thacker, `h_dry`).
13. **FNPF in the inner surf zone (Mase & Kirby)**: Hm0 is too high in the shallowest gauges and the
   deviation grows with grid refinement: +9 % / +26 % at h = 5 / 2.5 cm on the nightly grid
   (362 x 8), +18 % / +39 % on the release grid (725 x 10). **Not the breaking model**: switching
   the breaking filter off (`A 352 0`), an 8x / 40x larger breaking viscosity (`A 365`), a lower
   slope criterion or a wider coastline sponge change it by a few per cent only. Split into bands
   (now in the report), FNPF is within 6 % in the wind-wave band (f > 0.3 Hz) at every gauge; the
   excess is infragravity energy (f < 0.3 Hz): 1.6x measured at the toe, 1.9x at h = 2.5 cm. The
   linear (first-order) wave generation radiates free long waves that the static coastline
   reflects. Second-order (subharmonic) generation or an absorbing shoreline would address it.
14. **NHFLOW numerical damping at coarse resolution (Berkhoff shoal)**: the nightly Berkhoff case
   (dx = 0.1 m) gave about 60 % of the measured wave heights. **Not a 3D effect**: a 2D flume with
   the same wave (T = 1 s, h = 0.45 m, flat bed) loses 16 % of the height over 14 m at dx = 0.1 m
   (L/15), 2 % at dx = 0.05 m and nothing at dx = 0.025 m, independent of the number of sigma
   layers, HLL/HLLC (`A 511`) and RK2/RK3 (`A 510`). The default WENO-JS reconstruction damps;
   **WENO-Z (`A 527 1`)** keeps 0.95 along the flume at dx = 0.1 m and reduces the Berkhoff error
   from 0.44 to 0.31 (centreline section 7: 0.72 -> 0.49), the rest being the coarse grid over the
   shoal (L/10). The same damping shows in Mase & Kirby: NHFLOW is 10 % low in the wind-wave band
   already at the toe. WENO-Z as the NHFLOW default is optional patch 0005: the NHFLOW nightly
   benchmarks pass with it with practically unchanged errors (Beji & Battjes 0.355 -> 0.37, Mase &
   Kirby 0.217 -> 0.203, Stokes 5 / Synolakis / sloshing within 0.002).
15. Housekeeping: my first `git status` on the repo left an empty `.git/index.lock` (the sandbox
   cannot delete files there); I renamed it to `.git/stale_index.lock_from_claude` so git works.
   That file can be deleted.

## Verification run

Run in the Claude cloud workspace on hans_dev `ff1bb6819` (+ the orphan `momentum_FCC3.cpp`
removed), `make all` (-O3, REEFMG, no hypre), MPICH 4.2.3, **2 cores** (`--max-np 2`; the release
runs of the newer cases with `--max-np 1`). Nightly: 29 of 33 cases run, 25 PASS, 3 XFAIL, 1 FAIL
(Kleefsman). Release: 19 cases run, 18 PASS, 1 XFAIL. Reports and
one plot per case are in `~/Dropbox/Claude/REEF3D_benchmark_suite/results/` (`nightly_ff1bb68/`,
`release_ff1bb68/`). The tolerances in `case.json` were set from these runs with a margin; where no
run was possible they are provisional (marked below).

| case | nightly | release |
|---|---|---|
| `sflow_dambreak_ritter` | PASS, L1 0.021 | PASS, 0.020 |
| `sflow_dambreak_stoker` | PASS, L1 0.0026 | PASS, 0.0005 |
| `sflow_solitary_celerity` | PASS, c -0.2 %, H -1.1 % | PASS, c -0.2 %, H -1.1 % |
| `sflow_beji_battjes_bar` | PASS, rms <= 0.28 | PASS, rms <= 0.29 |
| `sflow_thacker_parabola` | PASS, L1 0.065 (wet gauges) | PASS, 0.041 |
| `sflow_conical_island_A` | PASS, gauges rms 0.36-0.43, run-up rms 0.055 d | not run (about 4 h on 1 core; tol provisional) |
| `sflow_conical_island_C` | PASS, gauges rms 0.32-0.48, run-up -0.11..+0.08 d | not run (tol provisional) |
| `sflow_mase_kirby_irregular` | PASS, Hm0 ratio 0.89-1.14, skewness within 0.26 | not run (tol provisional) |
| `nhflow_stokes5_propagation` | PASS, rms <= 0.067 | PASS, rms <= 0.057 |
| `nhflow_beji_battjes_bar` | PASS, rms <= 0.36 | PASS, rms <= 0.29 |
| `nhflow_synolakis_runup_nonbreaking` | PASS, R/d -8 %; profiles rms/H 0.09-0.15 (t* 30-60), 0.32 (t* 70) | PASS, R/d -5 %; profiles 0.09-0.16, 0.33 (t* 70) |
| `nhflow_synolakis_breaking_profiles` | PASS, rms/H <= 0.12 | PASS, rms/H <= 0.16 |
| `nhflow_berkhoff_shoal` | XFAIL, H/H0 ~60 % of measured (rms 0.34-0.72 on sections 1-7) | not run (tol provisional) |
| `nhflow_conical_island_C` | PASS, gauges rms 0.33-0.48, run-up within 0.06 d | not run (tol provisional) |
| `nhflow_mase_kirby_irregular` | PASS, Hm0 ratio 0.78-0.93 (low), skewness within 0.31 | PASS, Hm0 ratio 0.83-0.95 (low), skewness within 0.14 |
| `nhflow_sloshing_linear` | PASS, T -0.6 % | PASS, T -0.3 % |
| `fnpf_stokes2_propagation` | PASS, rms <= 0.13 (phase lead) | PASS, rms <= 0.12 |
| `fnpf_beji_battjes_bar` | PASS, rms <= 0.24 | PASS, rms <= 0.25 |
| `fnpf_ting_kirby_plunging` | XFAIL, max. H 1.7 m before x_b | xfail entry; dx = 0.025 m stopped at t = 91 s, release now at dx = 0.05 m (not run) |
| `fnpf_berkhoff_shoal` | PASS, rms 0.07-0.27 on sections 1-6, 8; 0.31 on section 7 (focus) | PASS, rms 0.09-0.18 on sections 1-6, 8; 0.21 on section 7 (focus 7 % high) |
| `fnpf_mase_kirby_irregular` | PASS, Hm0 ratio 0.98-1.09, 1.26 at h = 2.5 cm | XFAIL, Hm0 ratio 0.99-1.09 for h >= 7.5 cm, 1.18 / 1.39 at h = 5 / 2.5 cm (finding 13) |
| `fnpf_jonswap_spectrum` | PASS, Hm0 -3..-4 %, Tp 1.44 vs 1.50 s, shape rms 0.09 | PASS, Hm0 -0.2..-0.8 %, Tp 1.44-1.48 s, shape rms 0.06 |
| `cfd_dambreak_2d_martin_moyce` | PASS, mean 8 %, max 18 % (Z = 1.5) | PASS, mean 7 %, max 20 % |
| `cfd_dambreak_3d_kleefsman` | ff1bb68 without `F 46`: **FAIL**, N 61 stop at t = 1.77 s (volume -18 %). With `F 46 3` (patch 0004, now in the case): volume -0.7 %, passes 1.77 s, all signals within tolerance up to 1.8 s (P1 peak -4 %, 0.05 s late); run to 2 s not completed (small time step after the impact, > 1 h on 2 cores) | not run |
| `cfd_cylinder_re100` | ff1bb68: XFAIL, C_D 1.01 (P 81 sign error); with patch 0002: PASS, C_D 1.33, St +6 % (blockage), C_L 0.31 | not run (tol provisional) |
| `cfd_beji_battjes_bar` | PASS, rms <= 0.33, heights -15 % | not run (tol provisional) |
| `cfd_sphere_heave_decay` | PASS, T +4.1 %, zeta +1 %, equilibrium 0.001 R | not run (tol provisional) |
| `cfd_sloshing_linear` | PASS, T +0.44 % | PASS, T +0.13 % |
| `cfd_porous_dambreak_lin1998` | PASS, profiles rms 0.023-0.060 h | PASS, 0.022-0.047 h |
| `cfd_ting_kirby_plunging` | smoke test only (20 steps); a nightly run reached t = 0.7 s in 12 min on 2 shared cores (about 12 h for 45 s), stopped - run on the cluster | not run |
| `cfd_wave_force_chen2014` | smoke test only (20 steps); 3D, not run (cluster) | not run |
| `cfd_pier_scour` | smoke test only (20 steps) | not run |
| `cfd_breaking_wave_force_irschik2002` | (release only) | not run |

Smoke test = inputs accepted, grid and output files written. The 3D CFD cases (and the 2D CFD
plunging breaker) need several hours each at the nightly level, more than the 2-core sandbox
allows; their nightly and release levels are meant for the cluster.
