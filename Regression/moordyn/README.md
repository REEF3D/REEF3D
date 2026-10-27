# REEF3D::CORAL mooring against MoorDyn v2

The dynamic mooring lines of CORAL (`src/coral/mooring`, `mooring_system`) against **MoorDyn v2.7.2**
(MoorDyn-C, linked as a library): single lines, line networks and a moored floating body, the same case
in both codes.

| file | content |
|---|---|
| `coral_moordyn_cases.h` | the cases and the CORAL side (also used by the unit test) |
| `coral_moordyn_compare.cpp` | runs every case in CORAL and in MoorDyn, writes the time series and the differences; `--ref` writes the references |
| `build_moordyn.sh` | clones and builds MoorDyn (tag `v2.7.2`), builds `coral_moordyn_compare` |
| `plot_compare.py` | figures of the time series (`out/fig/<case>.png`, `out/fig/overview.png`) |
| `ref/<case>.csv` | MoorDyn's time series, the references of `Regression/unit/coral_moordyn_test.cpp` |

```bash
cd Regression/moordyn
./build_moordyn.sh                    # MoorDyn into ./moordyn_build, then coral_moordyn_compare
./coral_moordyn_compare               # all cases (or: ./coral_moordyn_compare bridle spar_decay)
./plot_compare.py                     # figures in out/fig
./coral_moordyn_compare --ref         # new references after a deliberate change of CORAL or of the cases
```

The unit test needs no MoorDyn: `cd Regression/unit && ./coral_tests.sh` (or build `coral_moordyn_test.cpp`
alone) runs the CORAL side of every case against `ref/`.

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
* **Error metric:** the largest difference of a channel over the run, relative to its range or to a tenth
  of its peak, whichever is larger.

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
