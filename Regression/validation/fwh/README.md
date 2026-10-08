# FW-H underwater acoustics: validation runs

Long validation runs of the permeable-surface FW-H module (`U` keywords, `src/acoustics_fwh*`,
kernel `src/acoustics_fwh_kernel.*`). The short regression versions are
`cases/cfd_3d_fwh_oscillating_sphere` and `cases/cfd_3d_fwh_sphere_endcaps`, the kernel unit test is
`unit/fwh_test.cpp`; `cases/cfd_3d_fwh_propeller_al` is the short version of 12.

Each case directory holds `ctrl.txt`, `control.txt` and `results.txt` (the output of the analysis
scripts for the runs below). The runs used `ahmet_dev`: 01–05 a build before the Box.dat diagnostic
(same FW-H), 06–09 commit `14c892c41`, 10 the `U 22` end caps on top of it.

## Running a case

```
mkdir run && cp cases/<case>/ctrl.txt cases/<case>/control.txt run/ && cd run
python3 ../scripts/make_inputs.py oscillating        # 01-05: floating.stl and 6DOF_motion.dat
python3 ../scripts/make_inputs.py shedding           # 06-10: 6DOF_motion.dat
cp ../../../cases/cfd_3d_heave_sphere_6dof/floating.stl .   # 06-10: sphere D = 1
DiveMESH && mpirun -np 8 REEF3D                      # DIVEMesh with grid format v2
python3 ../scripts/<script> <t0> <t1>                # as listed in results.txt
```

8 ranks; 01–05 ~2.4 s per step (120³ cells, 300–400 steps), 06–10 ~1 s per step (1.05 M cells,
~1100 steps for 50 s, ~1900 for 80 s).

## 01–05 Oscillating sphere (potential-flow dipole)

Sphere a = 0.1 m in surge x(t) = X0 (1 − cos ωt), X0 = 0.01 m, f = 1 Hz, closed tank of 8 m,
120³ sinh-stretched grid (Δx ≈ 0.02 m at the centre), fixed Δt = 5 ms. Reference
p = ρ a³ U̇ cosθ / (2 r²). Probes are taken relative to a probe at θ = 90° (the pressure level of the
closed tank floats).

| case | variant | script |
|---|---|---|
| 01 | 5 observers (directions, r = 5 m outside the domain) | `compare_oscillating.py 0.5 2.0` |
| 02 | radial line r = 0.3–1.5 m, probes off the rank borders | `compare_oscillating_radial.py 0.5 1.5` |
| 03 | velocity probes P 61, N40 = 4 (loop_cfd) | `compare_oscillating_velocity.py` |
| 04 | as 03, N40 = 14 (loop_cfd_df) | `compare_oscillating_velocity.py` |
| 05 | as 03, level set on (F 30 3, F 60 10) | `compare_oscillating_velocity.py` |

FW-H / analytic = 1.00–1.04 for r = 0.3–5 m. The CFD pressure probes give ≈ 0.75 × analytic and
the CFD velocity 1.27–1.64 × potential flow; the probes do not satisfy ρ ∂u/∂t = −∂p/∂x. 03/04 agree
to 7 digits and 03/05 bitwise, so neither the time scheme nor the single-phase φ is the cause. With
the moving direct-forcing body this case is not a clean reference for FW-H against probes.

## 06–07 Vortex-shedding sphere, box-size independence

Fixed sphere D = 1 (6DOF, `X 11 0 2 0 0 0 0` with one lateral kick in the first 4 s to break the
symmetry), U = 1 (`B 60 1`, `W 10 144`), ν = 0.002 (Re = 500), single phase, 240×66×66 grid
(Δx = 0.1 at the sphere). Box A [−1, 1.5]×[−0.95, 0.95]², box B [−1.6, 2.6]×[−1.58, 1.58]²
(B stopped at t ≈ 26 s).

- `boxforce.py`: box momentum balance from `REEF3D-CFD-FWH-Box.dat`, F_box = −(Lp + Lm − V + dM/dt),
  and the dipole radiated by FW-H, −(Lp + Lm + dM/dt), against the REEF3D body force.
- `curle.py`: FW-H at r = 50 D against Curle from the REEF3D body force, and at r = 3 D against probes.

| t = 10–25 s | −Lp | −Lm | F_box | F_REEF3D |
|---|---|---|---|---|
| A, Fx | 472 | −165 | 312 | 217 |
| B, Fx | 224 | +70 | 307 | 217 |
| A / B, Fy rms | | | 16.8 / 17.2 | 14.3 |

The FW-H far field equals the box momentum balance; the box force is independent of the box size
(1–2 %); the viscous traction on the box is negligible. The direct-forcing body force reported by
REEF3D (STL surface integral) is lower than the momentum the fluid loses (drag 217 vs ~310 N,
lateral ~0.85×) at this resolution (D/10); FW-H against Curle with the REEF3D force is ≈ 1.18 for
that reason. At r = 3 D FW-H and the probes disagree (pseudo-sound of the wake outside the box).

## 08–10 End caps

Same CFD as 06 (bitwise), 50 s, analysis window 25–50 s. 08: box A with the downstream face open
(`U 21 2`); 09: long box [−1, 3]×[−0.95, 0.95]², closed; 10: long box with 4 averaged end caps at
x = 3.0, 2.5, 2.0, 1.5 (`U 21 2`, `U 22 4 1.5`). `selfdipole.py` compares FW-H with the compact dipole
of the run's own box momentum balance, which removes the bias of the REEF3D body force and leaves
the end-cap / non-compactness error; `endcap.py` compares two runs of the same CFD.

| r = 50 D | lateral ratio | lateral corr. | lateral rms diff. | streamwise corr. |
|---|---|---|---|---|
| 06 box A, closed (cap at x = 1.5) | 0.97–1.05 | 0.995 | 11 % | 0.79 |
| 09 long box, closed (cap at x = 3) | 0.97–1.12 | 0.96 | 27–30 % | 0.37 |
| 10 long box, 4 averaged caps | 0.96–1.07 | 0.98 | 18–20 % | 0.57 |
| 08 box A, open | ~100 | ~0 | — | — |

- An open face does not work with incompressible input: on an open surface ∫ ρ u̇_n dS ≠ 0, a
  spurious 1/r monopole with the same amplitude at every far observer.
- The error of a closed cap grows with the wake that crosses it (better close to the body).
- Averaging 4 caps reduces the error of the long box from 30 to 19 %; the caps span 1.5 m, much
  less than the convected wavelength U·T, so the cancellation is partial.
- The weighted box momentum balance equals the closed box (295 / 295 N, Fy rms 7.59 / 7.66),
  which checks the cap weights.

## 11–12 Open-water propeller as actuator line

Ship module (`X 350`, `ship.dat`) on a fixed hub sphere (radius 0.1, `make_inputs.py hub`) held in a
current of 0.8 m/s (`B 60 1`, `X 101 0 0 180`): D = 1, n = 1 rps, Z = 4 (`propeller_blades 4`),
KT = 0.2, KQ = 0.03 constant (T = 200 N, Q = 30 Nm, J = 0.8), line width ε = 0.1 m, Δx = 0.05 at the
rotor (240×88×88), fixed Δt = 0.01 s (25 steps per BPF period), 5 s. Box [−0.75, 1]×[−0.75, 0.75]².
11 without, 12 with the mean flow `U 50 0.8 0 0`.

References (`al_reference.py`, `al_reference_smeared.py`): the linear incompressible pressure of the
momentum source, ∇²p = ∇·f, i.e. the sum of the rotating dipoles p = 1/(4π) Σ f·(x − y)/|x − y|³,
for thin lines and for the smeared distribution exactly as the code applies it.

| BPF amplitude [Pa] | thin lines | smeared (as applied) | FW-H 11 (no U 50) | FW-H 12 (U 50) | probe |
|---|---|---|---|---|---|
| (0, 1, 0), r = 2R | 0.656 | 0.444 | 0.496 | 0.451 | 0.470 |
| (0, 0, 1.5) | 0.083 | 0.057 | 0.057 | 0.058 | 0.068 |
| (0.5, 1.2, 0.3) | 0.125 | 0.089 | 0.083 | 0.089 | 0.105 |

- With the mean flow the BPF of FW-H agrees with the applied force distribution within 0–2 %; the
  phase lags by ~6° (the blade angle is taken at the start of the time step).
- Without `U 50` the steady part is wrong by O(100 Pa) (ρ U U_n on the box and the linear terms
  ρ U u' outside it); with `U 50` it is within ~1 Pa (the slipstream leaving the box remains).
- The Gaussian smearing reduces the BPF (azimuthal mode Z) by ~30 % against thin lines at ε = 0.1,
  R = 0.5: the line width is part of the model.

## Notes

- `X 120` (analytic sphere) gives no surface triangles with the current grid code; STL spheres used.
- `F 30 3` with inflow/outflow and the free surface above the domain (`F 60 10`) gave NaN in the
  first step; 06–10 are single phase.
