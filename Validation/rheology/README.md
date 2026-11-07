# Rheology validation (SFLOW, rheology core W 120)

These cases check the rheology core (`src/rheology_core.{h,cpp}`) and its SFLOW adapter (`src/sflow_rheology_core.{h,cpp}`, W 120 - W 128). They compare against analytical solutions, an ODE reference with the same closure, and the Carbopol experiments of Cochard & Ancey (2009).

- **Repo:** REEF3D `hans_dev` 1e67465d + the rheology core patch, DIVEMesh `hans_dev` 2ee2a7e.
- **Build:** Linux sandbox, g++ 13, `make release LTO=0` (`-O3 -march=native`), no HYPRE, Open MPI.
- **Unit test of the core:** `Regression/unit/rheology_core_test.cpp`. All checks pass.

## Model

The flow is a single depth-averaged layer of the rheological material: SFLOW, hydrostatic (`A 220 0`), RK3, HLL.

- **Basal friction.** After the explicit update of each stage, the speed s of every wet cell is solved from
  `s + dt tau_b(s)/(rho h) = |U*|`, which is backward Euler. The direction of U* is kept.
  - The layer stops when `|U*| <= dt tau_static/(rho h)`. This is the yield criterion: a layer below the friction angle stays exactly at rest, and the friction never changes sign.
  - Closed form for Coulomb and Voellmy. Bisection on the stress for Herschel-Bulkley and Herschel-Bulkley-Coulomb (whose plug-flow law gives U(tau) explicitly). Bisection on the speed for the others.
- **Inclined reference plane.** `W 20`/`W 21` act as along-slope gravity (body force h g_x). `W 22` is the gravity normal to the plane. With W 20 = W 21 = 0, the bed slope drives the flow as usual.
- **Pore pressure.** The basal ratio lambda sets the effective normal stress `sigma_n' = (1 - lambda) rho g_n h`.
  - It is constant, or advected (first-order upwind) and relaxed once per step: `lambda_h + (lambda - lambda_h) exp(-pi^2 D dt/(4 h^2))` (Iverson & Denlinger 2001).
  - The friction step uses the ratio at the middle of the step. Using the ratio at the start instead gave a 2.5 % error in case 04.

| key | type | meaning (default) |
|---|---|---|
| W 120 | int | basal law: 0 off, 1 Coulomb, 2 Voellmy, 3 Pouliquen-Forterre μ(h,Fr), 4 μ(I) (Bagnold profile), 5 Herschel-Bulkley plug flow (Bingham n = 1), 6 Herschel-Bulkley-Coulomb (0) |
| W 121 | double | basal friction angle δ [deg] (Coulomb, Voellmy, HBC) (30) |
| W 122 | double | Voellmy ξ [m/s²] (1000) |
| W 123 | 5 × double | Pouliquen-Forterre ζ1 ζ2 ζ3 [deg], β, L [m] (21, 30.7, 22.2, 0.136, 0.00065) |
| W 124 | 5 × double | μ(I): δ_s δ_2 [deg], I0, d [m], φ_s (20.9, 32.76, 0.279, 0.0005, 0.6) |
| W 125 | 3 × double | Herschel-Bulkley τ_y [Pa] (cohesion c for W 120 6), K [Pa sⁿ], n (0, 1, 1) |
| W 126 | int | pore pressure: 0 off, 1 constant ratio, 2 advected and relaxing (0) |
| W 127 | 3 × double | λ (initial or constant), λ_h (value the excess relaxes to), D [m²/s] (0, 0, 0) |
| W 128 | int | bed-normal gravity g cos²θ from the bed slope in the normal stress (0) |

The density of the material is `W 1`. W 90 1 in SFLOW, the legacy explicit Bingham source, still runs when W 120 is 0. W 120 does not yet act in the semi-implicit integrator (`A 210 12/13`); the code prints a note when both are set.

**Settings for fast granular and debris flows.** Every case here sets these two keys:

- `A 531 1.0e6`: switches off the Froude limiter. Its default, Fr ≤ 3, holds back the fronts of fast flows on slopes.
- `A 244 1.0e-5`: the wetting threshold. The default of 1 mm made the front of a dry dam break on a 35° slope lag by about 20 % (frictionless check, t = 2 s). With 1e-5 m the lag is 4 %.

## Running

```bash
cd Validation/rheology
tools/run_case.sh <REEF3D> <DiveMESH> 01_coulomb_dambreak_mangeney/cases/mangeney <run>/mangeney
python3 tools/rheo_eval.py mangeney <run>/mangeney 01_coulomb_dambreak_mangeney/results
python3 tools/rheo_eval.py arrest <run>/hb_core <run>/legacy 02_viscoplastic_arrest/results
python3 tools/rheo_eval.py cochard <run>/cochard 03_cochard_ancey_12deg/results         # 2 ranks: add "mpirun -np 2"
python3 tools/rheo_eval.py pore <run>/layer 04_pore_pressure_layer/results
g++ -O2 -std=c++17 -I../../src tools/uniform_ode.cpp ../../src/rheology_core.cpp -o uniform_ode
python3 tools/rheo_eval.py uniform <run>/u_ ./uniform_ode 05_uniform_layer_laws/results  # runs <run>/u_<case>
```

`rheo_eval.py` reads the free-surface VTP files with `../turbulence/tools/r3read.py` (Python 3, numpy, matplotlib).

## Results

### 01 Coulomb dam break on a 35° incline (Mangeney et al. 2000)

- **Set-up:** h0 = 1 m, reservoir 25 m long, δ = 25°, dx = 5 cm, 1D, slope-aligned frame (W 20 = g sin35°).
- **Reference:** the exact solution `h = (2c0 - x/t + mt/2)²/(9 g_n)`, `u = 2/3 (c0 + x/t + mt)`, with c0 = √(g_n h0) and m = g(sinθ − cosθ tanδ).
- **Error measure:** L1 over the region the solution covers (ahead of the rarefaction from the upstream wall).

| t (s) | L1(h) | L1(u) | front h > 1 mm, SFLOW / exact (m) |
|---|---|---|---|
| 1 | 0.11 % | 0.18 % | 5.65 / 6.35 |
| 2 | 0.14 % | 0.15 % | 13.40 / 14.55 |
| 3 | 0.19 % | 0.15 % | 22.95 / 24.65 |

The thin tail of the front (h < ~2 cm) lags by about 7 %, as it does in the frictionless SFLOW dam break. This comes from the wet/dry front, not the rheology. Serial and 2 ranks give identical results.

![](01_coulomb_dambreak_mangeney/results/01_mangeney_profiles.png)

### 02 Viscoplastic arrest on a horizontal plane

- **Set-up:** Carbopol 0.30 % (Cochard & Ancey 2009: τ_y = 89 Pa, K = 47.68 Pa sⁿ, n = 0.415), ρ = 1000. Lock 0.5 m × 0.3 m against a wall, 1D, 60 s.
- **Reference:** a deposit that has yielded everywhere satisfies `h² = 2 τ_y (x_f − x)/(ρ g)`. Volume conservation then gives x_f = 1.408 m and h(0) = 0.160 m.

| run | x_f (m) | h(0) (m) | volume (m²) | max \|u\| at 60 s |
|---|---|---|---|---|
| exact (fully yielded) | 1.408 | 0.160 | 0.15 | 0 |
| W 120 5, dx = 1 cm | 1.460 | 0.145 | 0.14997 | 3.0e-4 m/s |
| W 120 5, dx = 5 mm | 1.445 | 0.145 | 0.14997 | 2.9e-4 m/s |
| legacy W 90 1 (explicit 1.5τ_y + 3Ku/h) | — | — | — | emergency stop at step 3 (N 61) |

- **Front:** within 2.6–3.7 % of the fully yielded profile.
- **Wall region:** h(0) is 9 % lower than the analytical value. The part near the wall stops below yield after the dynamic phase, so the deposit is not fully yielded there, and more material reaches the front.
- **Creep:** the front region still creeps at about 1e-4 m/s, with h dh/dx about 15 % above τ_y/(ρg). It slows down slowly and is independent of the grid.

![](02_viscoplastic_arrest/results/02_arrest.png)

### 03 Carbopol on a 12° incline (Cochard & Ancey 2009)

- **Set-up:** Carbopol 0.30 %, reservoir 0.51 m × 0.30 m with walls, 33 cm deep at the gate, free surface horizontal (43 kg). Gate removed instantly, 2D, dx = 2 cm, 2 ranks, 30 s. dx = 1 cm (10 s) gives the same front to 1 cm. Volume conserved to 0.1 %.
- **Reference:** the paper gives the front only in figures. The text gives a front speed of about 0.7 m/s in the inertial phase (t < 1 s), an end of the inertial phase at 1.20–1.45 m (12°, 0.25–0.40 %), and slow creep afterwards with x_f ∝ t^0.275 to 1.9 m over hours (no arrest).

| | SFLOW | Cochard & Ancey |
|---|---|---|
| mean front speed 0–1 s | 0.78 m/s | ~0.7 m/s (gate opened over 0.8 s) |
| x_f at 2 s / 10 s / 30 s | 0.86 / 0.94 / 1.00 m | 1.20–1.45 m at the end of the inertial phase |
| late regime | x_f ∝ t^0.06, arrests | x_f ∝ t^0.275, creeps for hours |

- **Inertial phase:** the speed matches.
- **Late regime:** the model stops about 25–30 % short and arrests. That is what a Herschel-Bulkley material must do: the deposit thins to the critical depth h_c = τ_y/(ρ g sinθ) = 4.4 cm. 0.043 m³ at h_c covers about 1 m², which is the area the model deposit covers.
- **Why the experiment goes further:** a front at 1.9 m means a mean depth well below h_c. A pure Herschel-Bulkley yield stress can't hold that. Wall slip and thixotropy of Carbopol (not in the model) are the likely reasons, as are the 3D non-shallow reservoir and the gate opening.

![](03_cochard_ancey_12deg/results/03_cochard_front.png)
![](03_cochard_ancey_12deg/results/03_cochard_map.png)

### 04 Uniform layer with dissipating pore pressure

- **Set-up:** h = 1 m on 20°, δ = 25°, λ from 0.5 to 0 with T = 4h²/(π²D) = 5 s, 1D, 300 m (speed taken at x = 150 m).
- **Reference:** `u = (g sinθ − g cosθ tanδ (1−λ_h)) t + g cosθ tanδ (λ0 − λ_h) T (1 − e^{−t/T})` until u = 0. After that the layer stays at rest.

| | SFLOW | analytic |
|---|---|---|
| peak speed | 2.1448 m/s | 2.1455 m/s |
| stop | first output at rest 9.80 s | 9.781 s |
| max \|u − u_exact\| | 1e-4 m/s | |

Serial and 2 ranks give identical results.

![](04_pore_pressure_layer/results/04_pore_pressure.png)

### 05 Uniform layers, all basal laws

An infinite uniform layer on an incline, starting from rest, compared with `du/dt = g_x − τ_b(u)/(ρh)` integrated with the same closure in steps of 1e-6 s (`tools/uniform_ode.cpp`). The domains are long enough that the walls don't reach the probe.

| law | set-up | steady u, SFLOW / ODE (m/s) | max error / max u |
|---|---|---|---|
| Voellmy | h = 1 m, 30°, δ = 21.8°, ξ = 200 | 5.4004 / 5.3965 (t = 8 s) | 0.13 % |
| Pouliquen-Forterre | h = 5 mm, 25°, glass beads | 0.14247 / 0.14246 | 0.7 % |
| μ(I) | h = 5 mm, 25°, glass beads | 0.08688 / 0.08688 | 2.0 % |
| Herschel-Bulkley | h = 0.1 m, 12°, Carbopol 0.30 % | 0.1202 / 0.1202 | 0.8 % |
| Herschel-Bulkley-Coulomb | h = 0.1 m, 12°, c = 20 Pa, δ = 5.7° | 0.04998 / 0.04998 | 1.2 % |

The largest errors are in the first steps of the acceleration (RK stage weights and output times). The steady states agree to 4–5 digits.

![](05_uniform_layer_laws/results/05_uniform_layers.png)

## Not covered yet

- Mesh refinement (G 1): the friction step runs in every level's stage. The pore pressure is advected per level and not exchanged between levels.
- Semi-implicit SFLOW (A 210 12/13): no friction step yet.
- Restart from checkpoints: the pore pressure ratio is not in the checkpoints yet and starts again from W 127.
- The local (3D) laws of the core are only covered by the unit test. They are for the CFD and NHFLOW adapters (later stages).
