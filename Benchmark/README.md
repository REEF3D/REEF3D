# REEF3D literature benchmark suite

Established test cases from the literature for **REEF3D::CFD** (incl. 6DOF and sediment),
**REEF3D::NHFLOW**, **REEF3D::FNPF** and **REEF3D::SFLOW**: laboratory experiments, exact and
semi-analytical solutions. Each case runs REEF3D, compares the output with the published reference
data and gives PASS/FAIL against a tolerance.

The suite complements the regression suite in `Regression/`:

| set | question | cases |
|---|---|---|
| `Regression/` | does the code give the **same** numbers as before? (bitwise / tolerance) | short runs, stored references |
| **`Benchmark/` (this)** | does REEF3D reproduce the **established literature benchmarks**? | 22 cases, experiments + exact solutions |

It uses the same case format (`case.json` + `control.txt` + `ctrl.txt`, `"base"` + overrides) and
is self-contained (Python 3 standard library; matplotlib only for plots; it does not need the
regression dump). The intended place in the repository is `Benchmark/` next to `Regression/`
(patch `261004_REEF3D_Benchmark-on-8714770.patch`).

## Files

| path | content |
|---|---|
| `benchmark.py` | runner, checkers, report, plots |
| `cases/<name>/` | `case.json` (description, reference, check, levels), `control.txt`, `ctrl.txt`, extra files |
| `refdata/` | reference data with provenance (`*/PROVENANCE.md`) |
| `tools/make_cases.py` | writes `cases/` (documents the gauge lists and geometry; the case files are the source of truth) |
| `tools/berkhoff_geo.py` | `geo.dat` bathymetry of the Berkhoff shoal (run automatically by `write`/`run`) |
| `tools/run_suite.sh` | run one level + plots |
| `results/` | reports and plots of the verification runs below |

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
| `nhflow_synolakis_runup_nonbreaking` | NHFLOW | Synolakis (1987) run-up law | max. run-up R/d, H/d = 0.0185 | dx 0.1 | dx 0.025 |
| `nhflow_synolakis_breaking_profiles` | NHFLOW | Synolakis (1987) | 12 surface profiles, H/d = 0.28 | dx 0.1 | dx 0.025 |
| `nhflow_berkhoff_shoal` | NHFLOW | Berkhoff et al. (1982) | H/H0 on sections 2, 3, 5, 7 (3D) | dx 0.1 | dx 0.04 |
| `fnpf_stokes2_propagation` | FNPF | 2nd-order Stokes | P 51 vs. theory P 50 | 400 x 8 | 800 x 10 |
| `fnpf_beji_battjes_bar` | FNPF | Beji & Battjes (1993) | eta(t) at WG4-WG11 | 800 x 8 | 2000 x 10 |
| `fnpf_ting_kirby_plunging` | FNPF | Ting & Kirby (1995) | breaking point (max. H) | dx 0.05 | dx 0.05, 100 s |
| `fnpf_berkhoff_shoal` | FNPF | Berkhoff et al. (1982) | H/H0 on sections 2, 3, 5, 7 (3D) | dx 0.1 | dx 0.04 |
| `cfd_dambreak_2d_martin_moyce` | CFD | Martin & Moyce (1952) | surge front Z(T) | a/20 | a/40 |
| `cfd_dambreak_3d_kleefsman` | CFD | Kleefsman et al. (2005), SPHERIC 2 | pressure P1, P3: peak, rms | dx 0.025 | dx 0.01 |
| `cfd_cylinder_re100` | CFD | Williamson St(Re); C_D compilation | St, mean C_D at Re = 100 | D/25 | D/50 |
| `cfd_beji_battjes_bar` | CFD | Beji & Battjes (1993) | eta(t) at WG4-WG11 | dx 0.02 | dx 0.01 |
| `cfd_ting_kirby_plunging` | CFD | Ting & Kirby (1995) | breaking point (max. H) | dx 0.02 | dx 0.005 |
| `cfd_wave_force_chen2014` | CFD | Chen et al. (2014) | eta(t), inline force on cylinder | dx 0.05 | dx 0.025 |
| `cfd_breaking_wave_force_irschik2002` | CFD | Irschik et al. (2002), GWK | slamming force | – | dx 0.05 |
| `cfd_sphere_heave_decay` | CFD 6DOF | Hulme (1982) linear theory, Archimedes | damped period, damping ratio, equilibrium | dx 0.05 | dx 0.025 |
| `cfd_pier_scour` | CFD sediment | HEC-18 (CSU), Sumer et al. (1992) | equilibrium scour depth S/D | smoke test, 1 h sediment time | dx 0.02, 21 h sediment time |

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

New case: copy a case folder, edit, `./benchmark.py run --cases <new> --level nightly`, look at the
plot, then set the tolerances. A case that deviates from the reference for a known reason gets
`"xfail": "<reason>"` (report shows XFAIL and the suite still passes; XPASS once it is fixed).

## Reference data

All reference data are in `refdata/`, each file with a header giving the source URL, the original
reference, units and scaling; `refdata/*/PROVENANCE.md` summarise sources and confidence.
**Important: most measured series are secondary, digitised copies** (from the test suites of
Basilisk, PySPH, Lethe and from the REEF3D tutorials), not the original instrument records:

| data | source | confidence |
|---|---|---|
| Beji & Battjes WG4-WG11 eta(t) | Basilisk `test/bar.c` data files | digitised (Yamazaki et al. 2009 figure) |
| Berkhoff sections 2, 3, 5, 7 | Basilisk `examples/shoal.c` data files | digitised |
| Synolakis H/d = 0.28 profiles | Basilisk `test/beach.c` data files | digitised |
| Synolakis run-up law | Synolakis et al. (2008) PAGEOPH | formula |
| Martin & Moyce surge front (a = 2.25 in) | PySPH `db_exp_data.py` (Fig. 3 of the paper) | digitised |
| Kleefsman P1, P3 | PySPH | digitised, coarse |
| Chen et al. eta + force, Irschik et al. force | REEF3D tutorials 11_14, 11_15 | digitised |
| Hulme hemisphere added mass / damping | LHEEA `Hulme_Heaving_Sphere` (Hulme's table) | tabulated |
| Ritter, Stoker, solitary, Stokes, run-up law, St(Re), C_D | formulas / published values | exact / tabulated |

Because of the digitising, the time series comparisons allow a fitted time lag and use rms errors
of 0.1-0.4 as tolerances rather than tighter values.

## Status of the reference data (what is missing)

- **Kleefsman/SPHERIC Test 2**: the water heights H1-H4 and pressures P2, P4-P8 are only in
  `SPHERIC_Test2.zip` (`test_case_2_exp_data.xls`) on spheric-sph.org (Test 2), which could not be
  downloaded from the sandbox. With the file, add `refdata/dambreak/kleefsman_2005_<sensor>.txt`
  (t [s], value) and further signals in `cases/cfd_dambreak_3d_kleefsman/case.json`.
- **Beji & Battjes original records** (Dingemans 1994, Delft Hydraulics H1684, 0.05 s sampling):
  likely in FUNWAVE `BENCHMARK_FUNWAVE/car_luth_AC` and ERDC `air-water-vv`; would replace the
  digitised series.
- **Ting & Kirby (1995)** measured wave height and set-up along the flume (only figures found): the
  cases check only the breaking point.
- **Synolakis H/d = 0.0185 and 0.3 profiles** (NOAA NCTR `.xls`), **Berkhoff sections 1, 4, 6, 8**,
  **Whalin (1971)** harmonics: not found as text. Whalin's geometry is in `refdata/waves/whalin_shoal`
  for a future case.

## Findings on hans_dev while setting up the suite

Observed while building and running the cases (hans_dev `ff1bb6819`, Linux, MPICH 4.2.3, REEFMG,
`make all`). Worth a look; none of them was changed in the code.

1. (Fixed on hans_dev in the meantime: `src/momentum_FCC3.cpp`, which still included the deleted
   `momentum_FCC3.h` and broke `make`, has been removed after `ff1bb68`.)
2. **CFD implicit diffusion at solids (default `D 22 1`)**: the laminar cylinder at Re = 100 with
   `D 20 2` and the default `D 22 1` gives C_D = 0.62 and no vortex shedding after 40 s (400 D/U);
   explicit diffusion (`D 20 1`) or `D 22 2` give C_D = 0.92 at t = 0.8 s and the expected
   behaviour. Zero normal gradient towards solid cells removes the wall shear in the implicit
   diffusion (the earlier `tests/validation` README noted the same for walls). Consider `D 22 2` as the default.
3. **`B 20 3` with a DIVEMesh solid**: with `B 20 3` + `D 22 2` the cylinder force is ~0
   (-5e-5 N instead of 0.18 N) and u_max drops from 1.28 to 1.15 m/s, i.e. the flow does not see
   the solid. `B 20 3` seems to act on domain walls only.
4. **SFLOW `A 243 1`**: cells that are dry at the start never become wet, so a dam break onto a
   dry bed does not move at all (u_max ~ 1e-11). Default `A 243 2` works.
   `sflow_eta_wetdry.cpp`: option 1 only sets `wet = 1` where the water level is already above the
   threshold.
5. **3D dam break with box (Kleefsman set-up, dx = 0.025 m, N 40 3, level set)**: the water volume
   drops from 0.656 to 0.532 m3 (-19 %) until t = 1.7 s and the run stops at t = 1.71 s with the
   N 61 velocity limit (u_max 365 m/s). The impact at P1 comes 0.07 s late (9.4 vs 10.9 kPa).
6. **FNPF plunging breaker (tutorial 9_3 set-up)**: at the tutorial resolution (dx = 0.05 m) the
   maximum wave height is reached 1.9 m before the measured breaking point; at dx = 0.025 m the
   wave train becomes irregular at t ~ 75 s and the run stops at t = 91 s (N 61).
7. **FNPF regular waves**: 2nd-order Stokes waves in deep water run about 0.4 % faster than the
   theory (phase lead 0.017 s at x = 15 m, 0.032 s at x = 25 m), independent of the grid level.
8. **P 81 drag of the cylinder (Re = 100, nightly grid D/25, `D 22 2`)**: shedding frequency
   (St = 0.1745, +6 % with 6 % blockage) and lift amplitude (C_L = 0.31, i.e. rms 0.22 vs. 0.225-0.235
   in the literature) are right, but the mean drag is C_D = 1.01 instead of 1.33. The viscous part
   in `force_force.cpp` uses u/dx one cell off the wall and the pressure is sampled P 91 cells off the
   surface; one of the two probably underestimates the force.
9. Housekeeping: my first `git status` on the repo left an empty `.git/index.lock` (the sandbox
   cannot delete files there); I renamed it to `.git/stale_index.lock_from_claude` so git works.
   That file can be deleted.

## Verification run

Run in the Claude cloud workspace on hans_dev `ff1bb6819` (+ the orphan `momentum_FCC3.cpp`
removed), `make all` (-O3, REEFMG, no hypre), MPICH 4.2.3, **2 cores** (`--max-np 2`). Reports and
one plot per case are in `results/nightly_ff1bb68/` and `results/release_ff1bb68/`. The tolerances
in `case.json` were set from these runs with a margin; where no run was possible they are
provisional (marked below).

| case | nightly | release |
|---|---|---|
| `sflow_dambreak_ritter` | PASS, L1 0.021 | PASS, 0.020 |
| `sflow_dambreak_stoker` | PASS, L1 0.0026 | PASS, 0.0005 |
| `sflow_solitary_celerity` | PASS, c -0.2 %, H -1.1 % | PASS, c -0.2 %, H -1.1 % |
| `sflow_beji_battjes_bar` | PASS, rms <= 0.28 | PASS, rms <= 0.29 |
| `nhflow_stokes5_propagation` | PASS, rms <= 0.067 | PASS, rms <= 0.057 |
| `nhflow_beji_battjes_bar` | PASS, rms <= 0.36 | PASS, rms <= 0.29 |
| `nhflow_synolakis_runup_nonbreaking` | PASS, R/d -8 % | not run (tol provisional) |
| `nhflow_synolakis_breaking_profiles` | PASS, rms/H <= 0.12 | PASS, rms/H <= 0.16 |
| `nhflow_berkhoff_shoal` | XFAIL, H/H0 ~60 % of measured | not run (tol provisional) |
| `fnpf_stokes2_propagation` | PASS, rms <= 0.13 (phase lead) | PASS, rms <= 0.12 |
| `fnpf_beji_battjes_bar` | PASS, rms <= 0.24 | PASS, rms <= 0.25 |
| `fnpf_ting_kirby_plunging` | XFAIL, max. H 1.7 m before x_b | xfail entry; dx = 0.025 m stopped at t = 91 s, release now at dx = 0.05 m (not run) |
| `fnpf_berkhoff_shoal` | PASS, rms 0.22 | not run (tol provisional) |
| `cfd_dambreak_2d_martin_moyce` | PASS, mean 8 %, max 18 % (Z = 1.5) | PASS, mean 7 %, max 20 % |
| `cfd_dambreak_3d_kleefsman` | **FAIL**: N 61 stop at t = 1.71 s; impact 0.07 s late | not run |
| `cfd_cylinder_re100` | XFAIL: St +6 %, C_L ok, C_D 1.01 (-24 %) | not run (tol provisional) |
| `cfd_beji_battjes_bar` | PASS, rms <= 0.33, heights -15 % | not run (tol provisional) |
| `cfd_sphere_heave_decay` | PASS, T +4.1 %, zeta +1 %, equilibrium 0.001 R | not run (tol provisional) |
| `cfd_ting_kirby_plunging` | smoke test only (20 steps) | not run |
| `cfd_wave_force_chen2014` | smoke test only (20 steps) | not run |
| `cfd_pier_scour` | smoke test only (20 steps) | not run |
| `cfd_breaking_wave_force_irschik2002` | (release only) | not run |

Smoke test = inputs accepted, grid and output files written (the 3D CFD cases need about 2-4 h each
on 2 cores at the nightly level). The release level of the 3D cases needs the cluster.
