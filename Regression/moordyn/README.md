# REEF3D::CORAL mooring against MoorDyn v2

The dynamic mooring lines of CORAL (`src/coral/mooring`, `mooring_system`) against **MoorDyn v2.7.2**
(MoorDyn-C, linked as a library): single lines, line networks and a moored floating body, in still water,
currents and waves, the same case in both codes. Two cases also compare with published results: the flexible
jumper of Lowe and Langley (2006) and the chain tests of Azcona et al. (2017) in the basin of Ecole Centrale
de Nantes.

| file | content |
|---|---|
| `coral_moordyn_cases.h` | the cases and the CORAL side (also used by the unit test) |
| `coral_moordyn_compare.cpp` | runs every case in CORAL and in MoorDyn, writes the time series and the differences; `--ref` writes the references |
| `build_moordyn.sh` | clones and builds MoorDyn (tag `v2.7.2`), builds `coral_moordyn_compare` |
| `plot_compare.py` | figures of the time series (`out/fig/<case>.png`, `out/fig/overview.png`) |
| `ref/<case>.csv` | MoorDyn's time series of the unit cases, the references of `Regression/unit/coral_moordyn_test.cpp` |
| `ext/lowe_langley_2006_case<n>.csv` | the published top tension of the Lowe and Langley jumper (digitized for MoorDyn's tests, BSD-3-Clause) |

```bash
cd Regression/moordyn
./build_moordyn.sh                    # MoorDyn into ./moordyn_build, then coral_moordyn_compare
./coral_moordyn_compare               # all cases (or: ./coral_moordyn_compare bridle spar_decay)
./plot_compare.py                     # figures in out/fig
./coral_moordyn_compare --ref         # new references after a deliberate change of CORAL or of the cases
```

The unit test needs no MoorDyn: `cd Regression/unit && ./coral_tests.sh` (or build `coral_moordyn_test.cpp`
alone) runs the CORAL side of the unit cases (`unit` in `coral_moordyn_cases.h`, 26 of 33) against `ref/`,
the published and measured values included; `./coral_moordyn_test <case>` runs any case.

## Set-up

* **Channels:** the mooring force on the body (Fx, Fy, Fz) and the moment My about its origin, the
  tensions of the end segments of every line (with the axial damping), the positions of the free points,
  and for a free body its surge, heave and pitch.
* **Same model in both codes:** the line types (d, mass, EA or a stress-strain table, BA or a damping
  ratio, Cdn, Can, Cdt, Cat), points (mass, volume, CdA, Ca), seabed (flat, kb, cb, no friction), still
  water and an optional current profile. CORAL uses MoorDyn's axial damping law (`BA_slack`).
* **Time stepping:** the same sub-step in both codes (CORAL's automatic one) and RK4. MoorDyn starts from
  its stationary initial-condition solver (dynamic relaxation for the heavy chain). Its free points start
  where CORAL's static start put them, because MoorDyn's solver stalls from a rough guess.
* **Body ends:** they sit on a coupled MoorDyn body, so MoorDyn's forces include the inertia of the end
  nodes, like CORAL's attachments.
* **Prescribed motion:** a smooth ramp into a harmonic surge, heave and/or pitch. CORAL steps 0.01 s and
  moves the ends on Hermite curves. MoorDyn takes 1 ms coupling steps: it moves a coupled body linearly
  from its start velocity, which costs a·h²/2 of fairlead position, and stiff chain segments
  (EA/l ≈ 1e8 N/m) turn that into about 15 kN at h = 0.01 s.
* **Free body:** a planar rigid body (surge, heave, pitch) with linear stiffness and damping and a
  harmonic surge force. CORAL integrates it with generalized-α (planar joint, `linear6`). The MoorDyn
  driver integrates the same equations with semi-implicit Euler and MoorDyn's forces (1 ms).
* **Waves:** regular components of CORAL's linear wave field (finite depth). MoorDyn takes the same
  kinematics at its nodes through its API (`WaveKin 1`, at the middle of each coupling step), so both codes see
  the same water; MoorDyn's surface then stays at z = 0, CORAL's nodes see the elevation.
* **Error metric:** the largest difference of a channel over the run, relative to its range or to a tenth
  of its peak, whichever is larger. Snap loads (a light chain going slack and snapping taut) are compared by
  the mean of the per-period peaks of the force magnitude and the tensions over the second part of the run,
  because the instants of individual snaps are sensitive.

## Results (MoorDyn 2.7.2, October 2026)

| case | | largest difference | unit test tolerance |
|---|---|---|---|
| `oc3_static` | OC3 chain at rest (static start) | 0.00 % (0.03 N of 45 kN) | 0.2 % |
| `oc3_surge` | OC3 chain, fairlead surge 8 m, T 20 s | 0.03 % | 0.5 % |
| `oc3_heave` | OC3 chain, fairlead heave 3 m, T 10 s | 0.06 % | 0.5 % |
| `oc3_zeta` | OC3 chain, BA = -0.8 (damping ratio), heave 3 m | 0.05 % | 0.5 % |
| `taut_poly` | taut polyester, surge 4 m and pitch 3°, T 12 s | 0.00 % | 0.5 % |
| `poly_curve` | taut polyester with a stress-strain table, surge 6 m | 0.01 % | 0.5 % |
| `chain_snap` | heavy chain in 50 m, heave 4 m, T 6 s (slack, snap loads) | 1.68 % | 2.5 % |
| `chain_surface` | heavy chain through the free surface, fairlead at +5 m, heave 3 m | 0.84 % | 1.2 % |
| `oc3_current` | OC3 chain in a uniform cross current of 1 m/s (steady state) | 0.00 % | 0.2 % |
| `oc3_shear` | OC3 chain in a sheared current profile, surge 5 m | 0.03 % | 0.5 % |
| `bridle` | two wires to a connection point, chain to the anchor; surge and pitch | 0.14 % | 0.5 % |
| `clump` | clump weight on the seabed between two chains; heave and surge | 1.66 % | 2.5 % |
| `buoy` | submerged buoy between two chains; surge 5 m | 2.52 % | 3.5 % |
| `chain_wire_chain` | chain-wire-chain at massless points; surge 6 m | 0.04 % | 0.5 % |
| `spar_decay` | OC3-like spar on three lines, free surge decay from 10 m | 0.03 % (7 mm of 16 m) | 0.5 % |
| `spar_forced` | the spar under a harmonic surge force 0.5 MN, T 30 s | 0.02 % | 0.5 % |
| `oc3_waves` | OC3 chain in regular waves H 6 m, T 10 s, surge 3 m | 0.04 % | 2 % |
| `chain_irregular` | heavy chain in 50 m in five oblique wave components, heave 1 m | 1.0 % | 3 % |
| `spar_waves` | the spar in waves H 8 m, T 12 s (on the lines) with a surge force | 0.52 % | 1 % |
| `lowe_langley_1..6` | Lowe and Langley jumper (EI, BA -0.5): fairlead motion x, y, z; waves from 0, 45, 90° | 0.6 - 3.5 % (Fx, Fy of the jumper vary little: ≤ 0.15 % of its tension) | 1 - 5 % |
| `azcona_c1/2_static` | Azcona chain at rest | 0.6 % (MoorDyn's stationary solver stops short) | 1 % |
| `azcona_c1/2_T158 .. T474` | Azcona chain, snap loads, mean peaks | ≤ 0.3 % | 5 % |

## Published and measured results

**Lowe and Langley (OMAE 2006), flexible jumper.** A published time-domain solution, not a measurement:
170 m of line with bending stiffness (EI 120.8 kN m², BA -0.5) between an anchor 55 m down and a fairlead
5 m below the surface; the fairlead moves 10 m at T 27 s along x, y, z (cases 1-3), or waves of 5 m amplitude,
T 10 s pass the fixed jumper from 0°, 45°, 90° (cases 4-6). Their 'top tension' is the magnitude of the force
on the fairlead (channel `F`; the end segment's tension is 0.5 - 2 kN lower). Largest difference / rms over the
published record (kN, tension about 40 - 60 kN):

| case | CORAL | MoorDyn |
|---|---|---|
| 1 (x) | 0.44 / 0.10 | 0.45 / 0.11 |
| 2 (y) | 0.28 / 0.08 | 0.29 / 0.09 |
| 3 (z) | 0.63 / 0.16 | 0.66 / 0.19 |
| 4 (waves 0°) | 0.81 / 0.33 | 0.82 / 0.33 |
| 5 (waves 45°) | 0.97 / 0.32 | 0.97 / 0.32 |
| 6 (waves 90°) | 0.46 / 0.20 | 0.47 / 0.20 |

The unit test asks 1 kN (as MoorDyn's own test of these cases).

**Azcona et al. (Ocean Engineering 2017), chain in the ECN basin.** Measurements: 21 m of studless chain
(69 g/m, EA 3.4e5 N, equivalent diameter 3.4 mm) in 5 m of water, the fairlead at the water plane moved
horizontally in the plane of the line, anchor 19.364 m (configuration 1) or 19.870 m (2) away; seabed
20 N/m² and 0.1 N s/m² per length, Cdn 1.4, Cdt 0.67, Can 1, fresh water (the paper's own model data).

* Static fairlead tension: measured 8.13 N and 14.48 N; CORAL 8.114 N (-0.2 %) and 14.87 N (+2.7 %), MoorDyn
  the same (the paper's code: -0.4 % and +1.5 %; configuration 2 is near taut and sensitive to the 2 mm
  uncertainty of the distance and to the load cell in the 21 m).
* Dynamic tests: the text gives ratios only (the time series are in figures): at T 1.58 s the largest tension
  is "between 2 (Configuration 1) and 3 (Configuration 2) times" the quasi-static one, and three times larger
  in configuration 2 than in 1. The table's 0.25 m is taken as the **stroke** (amplitude 0.125 m): then CORAL
  gives 1.76, 2.74 and 3.07 (MoorDyn the same within 0.3 %); with an amplitude of 0.25 m both codes give 3.8,
  6.9 and 3.9. The unit test checks the mean peak against 2 and 3 times the quasi-static peak within 20 %
  (-12 % and -9 % now). To be confirmed against the paper's figures 5 and 6.

## Known differences (modelling, not errors)

* **Seabed contact.** MoorDyn's seabed damping also pulls a node down while it rises from the seabed;
  CORAL's contact only pushes. This explains `chain_snap` (1.68 %, and 0.17 % when CORAL is given
  MoorDyn's law) and `clump` (1.66 %, and 0.15 % with MoorDyn's law). It is also part of `buoy`, together
  with short axial ringing of the anchor chain where it lifts off: 2.5 %, and 1.7 % with MoorDyn's law,
  with an rms of 6 N on 6 kN.
* **Free surface.** CORAL takes a node's submerged fraction over its own piece of line; MoorDyn takes it
  per segment. Switching CORAL to MoorDyn's per-segment fraction changed `chain_surface` by less than
  0.1 %; its remaining 0.8 % is the fairlead motion within a coupling step.
* **Points at the surface.** MoorDyn's points always carry their whole volume as buoyancy, even above
  the water. CORAL takes only the submerged part, so a buoy floats at the surface instead of flying out.
* **Initial current.** MoorDyn's initial condition leaves the current out and the line then drifts into
  it; CORAL's static start already includes it. `oc3_current` compares the steady state, and `oc3_shear`
  starts its motion after the drift.
* **Stress-strain table.** Beyond the last point, MoorDyn holds the last tension and CORAL continues
  with the last slope. The cases stay inside the table.
* **Water depth.** MoorDyn-C lowers the water depth to the lowest point of its input, so the harness
  keeps seeded points at the seabed or above.
* **Waves at the surface.** With waves given through the API MoorDyn's free surface stays at z = 0; CORAL's
  nodes take the elevation (a node above the crest is dry). The wave cases keep their lines submerged.
