# REEF3D turbulence model review: SFLOW, NHFLOW, CFD

- **Source:** local working copy `~/Codelite/REEF3D/src`, HEAD `d83fe5287`, with uncommitted changes. Snapshot taken 2026-10-03.
- **Method:** read-only static review; no files were changed. Each module was audited in full. The items marked ✔ were also re-checked line by line during consolidation.
- **Severity:** **C** = wrong results or crash in common setups; **M** = wrong in some configurations; **m** = minor, dead code or inconsistency.
- **Confidence:** Confirmed / Likely / Needs-check.

---

## 0. Fix first

| # | Module | Location | Issue | Sev |
|---|---|---|---|---|
| 1 ✔ | SFLOW | `ini.cpp:62` | `A263=2.7` should be `A264=2.7`. As written, A264 (ce_γ) is never initialised and the A263 limiter default becomes 2.7 instead of 10. | C |
| 2 ✔ | SFLOW | `sflow_ifou.cpp:58-82` | `ivel1`/`jvel1` are never computed (`u_flux` fills `iadvec`, not `ivel1`), and `ul,ur,vl,vr` are never reset. The k-ω convection is broken. | C |
| 3 ✔ | SFLOW | `sflow_turb_kw_IM1.cpp:137` | The ω bed source uses the ε form u*⁴/h², which is dimensionally wrong for ω. | C |
| 4 ✔ | NHFLOW | `nhflow_idiff_w.cpp:106` (also `nhflow_ediff.cpp:147`) | Copy-paste: the second VH bracket should be `(VH[IJm1Kp1]-VH[IJm1Km1])`. | C |
| 5 ✔ | NHFLOW | `nhflow_idiff_u/v/w(_2D).cpp` cross terms | σz factor is missing in u/w and divides instead of multiplies in v. The ∂²/∂x∂z terms are wrong by D or D². | C |
| 6 ✔ | NHFLOW | `nhflow_idiff_w.cpp:99`, `nhflow_idiff_u_2D.cpp:100` | A stray `/p->W1` in the σxx term of `M.t`. | M |
| 7 ✔ | CFD | `kepsilon_IM1.cpp:69-71` | The ε BC/ghost elimination (`bckeps_start(...,30)`) runs **after** `psolv->start`. During the solve, ghost couplings act as ε = 0. | C |
| 8 ✔ | CFD | `komega_func_PLIC.cpp:135` | The T41 limiter has no `Sij2>1e-20` guard (the non-PLIC version has one), so it returns 0/0 = NaN in still flow. | C |
| 9 ✔ | CFD | `idiff2_PLICu.cpp:103` | `(1.0-H_ddy_p)` should be `(1.0-H_ddy_m)`. | M |
| 10 ✔ | NHFLOW | `nhflow_idiff_scalar(_2D).cpp` | `EV0` has no ghost exchange or MPI update and is not initialised. Halo diffusivity is 0, and the first step has PK0 = 0. | M |
| 11 ✔ | CFD | `rans_ini.cpp:35-42` | `uref` is read uninitialised when B60 = 0 and B90 = 0. | M |

---

## 1. SFLOW

The model is selected with `A 260`: 0 void, 1 k-ε, 2 k-ω, 3 Prandtl, 4 parabolic, 5 `kw_IM1_v1`.

### Critical

**S-C1 ✔ `ini.cpp:61-62`, A264 never initialised**
- The code reads `A263=10.0; A263=2.7;`. The second line was meant to set `A264`.
- `ceg(p->A264)` in `sflow_turb_ke_IM1.cpp:35` is therefore indeterminate, usually 0. With ceg = 0 the ε bed source vanishes. ε then stays at about 0, k has no sink, and ν_t sits on the floor.
- The limiter default is also wrong (see S-M2).

**S-C2 ✔ `sflow_ifou.cpp`, k-ω convection broken (A260 = 2 and 5)**
- `pflux->u_flux(ipol,uvel,iadvec,ivel2)` should pass `ivel1,ivel2`. The same applies to v.
- `ul/ur/vl/vr` are never reset. Once one is set to 1 it stays 1, so reversed flow is upwinded wrongly.
- Net effect: no inflow flux from the upwind cell, only a sink u·k/dx.
- Fix: pass `ivel1`/`jvel1`, and set `ul=ur=vl=vr=0` at the top of `aij`.

**S-C3 ✔ `sflow_turb_kw_IM1.cpp:137`, ω bed source has wrong dimensions**
- The term `3.456/cf^0.75 * cmu * u*⁴/h²` has units m²/s⁴, but the ω equation needs s⁻².
- Derive it from the Rastogi–Rodi equilibrium with ω = ε/(Cμk):
  - k = u*²/(3.6 √Cμ cf^¼)
  - ω = 3.6 u*/(√Cμ cf^¼ h)
  - P_ωv = βω² = C_ω u*²/h², with C_ω = 12.96 β/(Cμ √cf) ≈ 10.8/√cf
- For typical values the current term is about 10⁻³–10⁻⁴ too small. ω ends up far too small and ν_t too large.

**S-C4 k-ε production uses no-slip wall gradients wherever `flagslice4<0` or a neighbour is dry** (`sflow_turb_ke_IM1.cpp:168-181`)
- In 1D, all y-ghost cells have flag −10. Every cell then gets dudy = 2u/dx and a spurious P_k ≈ 4ν_t u²/dx². This can be orders of magnitude larger than bed production.
- In 2D, SFLOW walls are free-slip in momentum (`sflow_momentum_func.cpp:672-684`). The same treatment is applied at inflow/outflow ghosts and at every shoreline cell.
- Fix: add a `j_dir` guard and use zero-gradient (slip) there. Confirmed.

**S-C5 The v1 model (A260 = 5) is inconsistent**
- Line 117: the k sink `cmu/sqrt(W)` should be `cmu*sqrt(W)`.
- Line 139: `3.5*eddyv*Qw` has the wrong units. Qw is also always 0; see m-list.
- Line 141: the bed term needs Cμ^−1.5, but the code uses Cμ^+1.5 (Likely).
- Line 302: the wall value should use Cμ^0.5.
- Output: W is written to vtp under the name "omega".
- Recommendation: remove the model or label it experimental.

### Major

- **S-M1 Bed friction is inconsistent between turbulence and momentum.**
  - Turbulence uses `n = ks^(1/6)/26` and g = 9.81. Momentum (`sflow_rough_manning.cpp:44`) uses `/20` and `W22`. So cf_turb = 0.59·cf_mom and u* is about 0.77× too small.
  - With defaults (`A218=0`) momentum friction is off, yet the turbulence models still generate bed production from ks = B50.
  - Sediment `bedshear_sflow.cpp` uses a third, log-law definition.
  - Recommendation: compute cf in one shared place.
- **S-M2 The k-ε limiter is scaled by Cμ** (`sflow_turb_ke_IM1.cpp:109-111`). The code computes `cmu*MAX(MIN(k²/ε, A263·k/S), 1e-4ν)`.
  - The effective limit is 0.243 k/S. That is below the equilibrium value of 0.30 k/S, so standard k-ε is clipped by about 20%.
  - The 3D code applies Cμ only to k²/ε and uses T31·k/S.
  - SFLOW k-ω uses T31 while k-ε uses A263.
- **S-M3 The k-ω side-wall function is applied wherever any neighbour has flag < 0** (`kw:251-255`). This includes inflow/outflow cells and every cell in 1D. `wallf` also removes bed production, shear production and dissipation from those cells (`kw:116-122`).
- **S-M4 The k-ω wall ω value is overwritten by the solve** (`kw:298-300`). It is assigned before `psolv->start`. The 3D code uses a 1e20 penalty instead.
- **S-M5 NaN when ks = 0** (`kw:243-253`): log(30y/0) = ∞, and with k = 0 this gives 0·∞ = NaN. Also, `dist` is modified in the loop and never reset, and there is no `MAX(0.01, log())` guard as in 3D.
- **S-M6 k-ω has no wet/dry handling.**
  - `HP` uses 1e-20 instead of A244.
  - Dry cells get no identity rows.
  - ν_t is not zeroed in dry cells.
  - k-ε does all three (lines 115-120, 339-355).
- **S-M7 With defaults, ν_t never reaches the momentum.** `A212=0` selects `sflow_diffusion_void`, and no warning is printed. With explicit diffusion (`A212=1`), `sflow_etimestep` has no viscous dt limit, and `p->viscmax` is never computed (Likely).
- **S-M8 `b->kin` is only filled by k-ε.** With any other model, sediment `S16==4` (τ = 0.3ρk) gives τ = 0.
- **S-M9 The parabolic model uses face velocities** `P(i,j)`, `Q(i,j)` (`sflow_turb_parabolic.cpp:52`) instead of cell-centred averages. This underestimates u* at walls and shorelines.
- **S-M10 Robustness.**
  - k-ε: the k sink is explicit (`-MAX(eps,0)` on the rhs). The 3D code treats it implicitly as (ε/k)·k. In shallow cells, P_kv ∝ h^−4/3 and k can go negative.
  - k-ω: α(ω/k)P_k has no bound (3D `komega_func.cpp:248` has one). (Likely)

### Minor

- **k-ε:** `wallf` is computed but unused. The condition `wet[...]<0` can never be true.
- **k-ω:** `Qw` (line 163) is always 0. `Vw` and `gcval_omega` are unused. The headers redeclare `gcval_kin/eps`, which shadows the base class members.
- **k-ω `clearrhs` does not clear M.** It only works because `sflow_ifou` assigns with `=`. Switching to `sflow_iweno_hj` would break it.
- **Hard-coded constants:** the HP threshold differs (A244 vs 1e-20), and `9.81` is used instead of `W22`.
- **Spacing:** turbulence uses `DXM` everywhere, while momentum uses DXP/DXN.
- **Boundary conditions:** there is no turbulence inflow BC or relaxation, and k/ε are Neumann at inflow.
- **Memory:** convection and diffusion objects are never deleted, and `sflow_turbulence` has no virtual destructor.
- **NaN reliance:** `A263*k/S` with S = 0 relies on the MIN/MAX operand order to discard the NaN. This is unsafe with `-ffast-math`.
- **Restart:** k-ε/k-ω are not in the state files (SFLOW restart is currently disabled).
- **Outside scope, but found:** `sflow_vtp_fsf.cpp:113,120` has `p->count<<p->P46_ie`. `<<` should be `<`.

### Checked OK

- **Constants:** k-ε 1.44/1.92/1.0/1.3. k-ω α 5/9, β 3/40, σ = 2 (applied as ν + ν_t/σ).
- **Rastogi–Rodi k-ε sources:** correct form. c_εγ defaults to 2.7, versus 3.6 in the original paper.
- **2D production / strain:** the staggered stencils have correct offsets.
- **Discretisation:** WENO-HJ coefficients and `sflow_idiff::diff_scalar` are correct.
- **Ghost cells / MPI:** `gcsl_start4` calls and the flagslice4 MPI exchange are in place.
- **Output:** vtp offsets are consistent.

---

## 2. NHFLOW

The model is selected with `A 560`: 0 void, 1/21 k-ε, 2/22 k-ω, 31 LES.
- 21/22 differ from 1/2 only by the URANS filter ν_t = min(ν_t, A568·Cμ·Δ·√k).
- `control_parse.cpp:82` forces A512 = 2 (implicit diffusion) whenever A560 > 0.

### Critical

**N-C1 ✔ `nhflow_idiff_w.cpp:106` (also `nhflow_ediff.cpp:147`), VH cross derivative**
- The code is `(VH[IJp1Kp1]-VH[IJp1Km1]) - (VH[IJp1Kp1]-VH[IJm1Km1])`. The second bracket should be `(VH[IJm1Kp1]-VH[IJm1Km1])`.
- As written, the term becomes a first y-derivative divided by Δσ. It acts as a spurious source in every 3D run with ν_t > 0, and it grows with vertical refinement.

**N-C2 ✔ σz factor in the momentum cross derivatives**
- The correct form is ∂²f/∂x∂z = σz ∂²f/∂x∂σ (+σx terms).
- The code is inconsistent:
  - `nhflow_idiff_u.cpp:118` (WH), `nhflow_idiff_u_2D.cpp:106`, `nhflow_idiff_w.cpp:105` (UH) and `nhflow_idiff_w_2D.cpp:92` have **no** σz factor: wrong by a factor D.
  - `nhflow_idiff_v.cpp:104` **divides** by `sigz`: wrong by a factor D², e.g. ×400 at D = 20 m.
- For constant ν and ∇·u = 0, these terms should cancel against the doubled normal stresses. With D ≠ 1 m they don't, which leaves spurious stresses ∝ ν_t.

### Major

- **N-M1 ✔ σxx term in `M.t` divided by W1** (`nhflow_idiff_w.cpp:99`, `nhflow_idiff_u_2D.cpp:100`). `M.b` and the other operators have no `/W1`, so the upper half of σxx·∂f/∂σ is lost wherever the surface or bed is curved.
- **N-M2 Swapped DZP in `nhflow_idiff_w_2D.cpp:85/88`.** `M.t` uses `DZP[KM1]` and `M.b` uses `DZP[KP]`; `M.p` and all other operators use the opposite. The row sum is ≠ 0 on stretched σ grids.
- **N-M3 ✔ `EV0` is never ghost-exchanged or initialised.**
  - It is used by `nhflow_idiff_scalar(_2D)` for the k/ε/ω face diffusivity and by `PK0`.
  - The halo is 0 at MPI and periodic interfaces.
  - `nhflow_rans_io::ini` sets only `EV`, so on the first step PK0 = 0 and diffusion is molecular only.
  - CFD does `pgc->start4(p,eddyv0,24)`.
- **N-M4 `nhflow_scalar_advec_CDS2.cpp:52-56` blocks the y flux at MPI boundaries.** `j==0` and `j==knoy-1` are local indices. This affects k/ε/ω and suspended sediment in y-decomposed runs.
- **N-M5 Inflow k/ε/ω (B60 1) never enters the implicit system** (`nhflow_kepsilon_bc.cpp:220,308`; `nhflow_komega_bc.cpp:218,306`). The ghost coupling is eliminated using the cell's own value, i.e. zero gradient. The inflow profile written into the ghosts is ignored.
- **N-M6 Inflow / initial profile uses the wrong height** (`nhflow_rans_ini.cpp:96-98`).
  - `H=beddist` (distance to bed) is used in the Manning slope, where the flow depth is needed (CFD uses `depth+DXM`).
  - The profile is therefore wrong over the column: k is about 4.6× too high near the bed for H/D = 0.01.
  - Lines 71/74 divide by `EV` without a guard.
- **N-M7 k-ε wall function** (`nhflow_kepsilon_bc.cpp:113`).
  - There is no `MAX(0.01, log())` guard (k-ω and CFD have it). After the `dist=ks/30` clamp, u⁺ = 0 in shallow or rough wall cells. WALLF has already removed P_k and dissipation, so k then has no source or sink at all.
  - ks ≤ 0 is not clamped (CFD `ks_val` clamps it to 1e-4), so ks = 0 gives NaN.
- **N-M8 ω/ε production is unbounded when the EV0 floor is active** (`nhflow_komega_func.cpp:195`, `nhflow_kepsilon_func.cpp:196`). CFD `komega_func.cpp:237ff` has the cap; NHFLOW does not. (Likely)
- **N-M9 No default initialisation of k/ε/ω unless B60 1.** The B90 block is commented out. For wave cases, k-ω then starts with ω → 0⁺. (Needs-check with a test case)
- **N-M10 The VRANS ω source is multiplied by porosity** (`vrans_nhflow_source_omega.cpp:48`). The comment says it should balance βω∞². The two models then give different equilibria:
  - k-ω: ω = √n·ω∞
  - k-ε: ω = ω∞
  - At n ≈ 0.45, k differs by about 50% between models.
- **N-M11 LES Δ in 2D includes the dummy DYN** (`nhflow_LES_Smagorinsky.cpp:45`).

### Minor

- **Output:** A560 = 21 writes ε without a name declaration (`nhflow_rans_io.cpp:286,298,321,333` test only `A560==1`).
- **Vertex interpolation:** in `print_3D` (rans_io and les_io) it averages only over j and k, not i. The output is shifted by half a cell in x.
- **`print_2D`:** broken (wrong array size, missing write) and dead.
- **B11 handling:** k-ε tests `>0`, k-ω tests `==1`. `wallf_update` marks WALLF regardless, so with B11 = 0 wall cells lose P_k and ε with no replacement.
- **Wall roughness:** k-ε uses B50 and k-ω uses B57. At the bed, S20·S21 is used and S28 is ignored.
- **A567 meaning differs by model:**
  - k-ε treats T37 as a fraction of WL.
  - k-ω A567 = 1/2 treats T37 as an absolute length; A567 = 3 treats it as a fraction.
  - k-ω A567 = 2 is identical to A567 = 1.
- **`turb_relax_nhflow`:** relaxes KIN and EV but not EPS or EV0. `kinupdate` runs before relaxation.
- **ν_t in production:** k production uses EV (limited) while ε/ω production uses EV0. CFD k-ε uses the limited value for both. Possibly intentional.
- **Wet/dry fronts:** they act as k = ε = 0 Dirichlet (dry rows are identity with rhs 0). Needs-check.
- **Source timing:** `isource/jsource` use `d->WL` at time n instead of the RK-stage WL.
- **Negative k:** never clipped.
- **Momentum diffusion details:**
  - Factor 2 is applied to σx² + σy² in diff_w; only σz² should be doubled.
  - σx averaging differs between 3D and 2D.
  - Central-difference denominators use DZN instead of DZP.
- **`nhflow_ediff`** (laminar only): uses `DXP[IP]`/`DZP[KP]` for both faces, uses UH instead of VH in the diff_v σ term, contains the N-C1 copy-paste bug, and has no viscous dt limit.
- **LES:**
  - The s11…s23 loop is dead code.
  - EV is not zeroed in DF < 0 cells, and the VRANS eddyv treatment is not applied.
  - EV0 is not set.
  - Cs = 0.2 is hard-coded and there is no wall damping.
- **Driver robustness:** an unknown A560 value leaves `pnhfturb` uninitialised. With knoz = 1, the bed penalty and `epsfsf` hit the same cell.
- **Dead code:**
  - `flowdepth_inflow` iterates over the wrong list.
  - `kepsini_default` is declared but has no definition.
  - `plain_wallfunc` is empty.
  - The `sst_*` constants are unused.
  - A569 is read but unused.
  - The `kinval/epsval` stubs return 0.
- **Memory:** KN, EN, KIN, EPS, WALLF, PK, PK0 and PK_b are never freed.

### Checked OK

- **Constants:** k-ε and Wilcox-88 k-ω constants. SST a1 = 0.31, and F2 matches Menter.
- **Limiter:** T31 = 1/√3 (Durbin).
- **Strain and production:** the strain/production with the full σ chain rule, including the 2D reduction.
- **Implicit linearisations:** correct.
- **`diff_scalar`:** the σ-Laplacian is correct.
- **Wall-law forms:** correct, apart from N-M7.
- **Ghost cells:** `start20V/30V/24V`.
- **WL:** clamped to A544, so 1/WL is safe.
- **VRANS:** the Nakayama–Kuwahara k∞ and ε∞ forms.

---

## 3. CFD

The model is selected with `T 10`: 0 void, 1/21 k-ε, 2/22 k-ω, PLIC k-ω (F80 4), 31 Smagorinsky, 33 WALE.

### Critical

**F-C1 ✔ `kepsilon_IM1.cpp:69-71`, ε BC applied after the solve**
- The k block calls `bckeps_start` before `psolv` (line 54). The ε block calls `psolv` first, then `epsfsf`, then `bckeps_start`.
- `bicgstab_ijk::fillxvec` fills `x` only in FLEXLOOP, so non-MPI ghosts are 0 during the solve. As a result:
  - Inflow, outflow and symmetry boundaries act as ε = 0 Dirichlet.
  - The convective ε inflow is lost.
  - The wall ε is imposed only weakly, by overwriting after the solve.
  - The gcval == 30 elimination loop (`kepsilon_bc.cpp:94-134`) is effectively dead.
- Consequence: ε is too low and ν_t too high near every boundary.
- Fix:
  1. Call `bckeps_start(...,gcval_eps)` before `psolv->start`.
  2. Impose the wall ε by penalty, as in `komega_bc::wall_law_omega`.

**F-C2 ✔ `komega_func_PLIC.cpp:135`, T41 limiter NaN**
- `Qij2/Sij2` has no guard, so u = 0 gives 0/0 = NaN.
- `MIN(eddyv0,NaN)` returns NaN, which then propagates into the momentum.
- Copy the `Sij2_val>1e-20` guard from `komega_func.cpp:139-146`.

### Major

- **F-M1 Inflow and outflow k/ε/ω are effectively zero.**
  - k-ω x faces (`komega_bc.cpp:157-193`): IO==1 sets `M.s=0` with no rhs, i.e. a ghost value of 0. IO==2 is Neumann.
  - k-ω y/z faces: the condition requires `IO==0`, so both inflow and outflow are Dirichlet 0.
  - k-ε: the inflow ghost is never set (`gc_epol4.cpp:113-119` is commented out).
  - `kinval/epsval/eddyval` in `ioflow_f.cpp:52-69` and `iowave.cpp:51-62` are never used.
  - Those values are also wrong: the ω coefficient is too small by a factor of about 1/Cμ, and `kinval=0` gives `eddyval=0/0`.
- **F-M2 The PLIC ω production has no bound** (`komega_func_PLIC.cpp:234`). The non-PLIC version (`komega_func.cpp:246-255`) has one.
- **F-M3 k-ε lacks the k-ω handling for DF/SF bodies, VRANS and NWT.**
  - `wallf=1` at gcdf4 cells removes P_k and dissipation, but no `QGCDF4LOOP` wall law is added, so k has no source or sink there.
  - Also missing in k-ε: the `flagsf4` turn-off, `solid_forcing_lsm`, `vrans_wall_law_*`, `turb_relax`, T36 = 3, T44 and T45.
- **F-M4 ✔ `uref` is read uninitialised** in `rans_ini.cpp:35-42` (`turbulence.h:68`). The garbage value enters `tau_calc` and the initial k/ε/ω.
- **F-M5 `ifou.cpp:128-138` divides by the wrong cell width.** The positive-flow branch uses `DX[IM1]`. This is non-conservative on stretched grids (also affects suspended sediment).
- **F-M6 T36 free-surface damping overwrites ε/ω with a δ-weighted value.** This drives ε/ω towards 0 at the band edges, so ν_t spikes there (only T31/T34 catch it).
  - δ has units 1/m, so the result is grid-dependent.
  - T36 = 1/2 treat T37 as a length, but its default of 0.07 is dimensionless.
  - Suggestion: use `MAX(eps, eps_s)` or a Heaviside blend. (Likely; this is a design issue)
- **F-M7 Vegetation sources have the wrong dimensions** (`vrans_veg_k.cpp:47,75`, `_eps.cpp:47`, `_omega.cpp:47`). For example, `sqrt(uu*kin)` should be `sqrt(uu)*kin`. They also give NaN for k < 0.
- **F-M8 LES `T 21 2` silently disables the SGS model.** `LES_filter_f2::start` is empty, so ν_sgs = 0.
- **F-M9 The momentum wall law at DF bodies is applied to the wall-normal component, with a stale `deltaZ`** (`bcmom.cpp:60,70,80`, `bcmomPLIC.cpp:41,51,61`). The `cs = 1/4` cases are not excluded, as they are in the gcb loop.
- **F-M10 T45 buoyancy is an explicit, unbounded sink.** It can drive k strongly negative in air cells. Suggestion: Patankar-linearise it onto the diagonal when G_b < 0. (Likely)
- **F-M11 ✔ `idiff2_PLICu.cpp:103`:** `(1.0-H_ddy_p)` should be `(1.0-H_ddy_m)`.
- **F-M12 The 2D filter width includes DYN.** Affected: `LES_smagorinsky.cpp:67`, `LES_WALE.cpp:~78`, and the k-ω T39 floor (`komega_func.cpp:192` + PLIC). In 2D use √(DXN·DZN).
- **F-M13 `etimestep.cpp:160,174,181,188`:** the viscosity term inside the square root is `|u|/dx + visc`, which adds m²/s to 1/s. It needs the 1/dx² factor. This affects only explicit diffusion.
- **F-M14 The VRANS wall law double-sources k** (`komega_bc.cpp:221-276`). Those cells keep `wallf=0`, so the bulk P_k and dissipation are added on top. (Needs-check)

### Minor

- **k-ω inflow ν_t override** (`komega_func.cpp:150-167` + PLIC):
  - hard-codes `0.212` instead of T31;
  - writes up to 4 cells beyond the local subdomain, so results depend on the decomposition;
  - goes out of range if `knox==1`.
- **`roughness.cpp:62,65`:** `topo(i+1,j,k-1)` should be `topo(i+1,j,k)`.
- **T44 sampling height:** velocity is taken at `bed+T44*DZN`, but u⁺ is evaluated at 0.5·DZN. (Needs-check)
- **Restart:** `eddyv0` is not in the state files. The first k-ω step after a restart has no turbulent diffusion and no ω production.
- **`rans_ini.cpp`:**
  - line 94 `eps=kin/eddyv` has no guard;
  - lines 109 and 115 omit T10 = 21/22;
  - global `DXM` is used.
- **Strain limiter:** `T31*k/S` with S = 0 relies on the MIN/MAX operand order (`kepsilon_func.cpp:83`, `komega_func.cpp:104,127,164`).
- **Limiter consistency:** k-ε always applies T31, k-ω only with T34. k-ε URANS uses global DXM; k-ω uses local Δ and T23.
- **Two gradient routines:** `strain::pk` uses `gradient_momscalar` (not solid-aware), while `strainterm/Sij2/Qij2` are flag-aware.
- **`wallf` also covers bc 6/7/8** (wave inflow/outflow, AWA), so those cells get no P_k or dissipation. (Needs-check)
- **ε/ω wall values:** the roughness clamp `30y<ks` is not applied.
- **T39 SGS floor:** an extra √2 makes the effective Cs 0.283 instead of 0.2.
- **No validation:** unknown T10, T12 or T21 values leave pointers uninitialised.
- **Dead code:**
  - the gcval == 30 loop (see F-C1);
  - `ioflow_turbulence.cpp`;
  - the SST constants;
  - T43 (read but unused);
  - the T42 comment says "λ1" but it is λ2;
  - `kepsini_default` and `les_io::tau_calc` are declared but have no definition;
  - `LES_filter_box` is an identity copy.

### Checked OK

- **k-ε:** constants, implicit sources and ν_t form.
- **k-ω:** Wilcox-88 sources and the Larsen–Fuhrman limiter (Sij2 = 2SijSij, Qij2 = 2ΩijΩij).
- **Diffusion:** `idiff2` uses ν + ν_t/σ.
- **Strain:** √(2SijSij) from the staggered velocities.
- **Wall values:** ε* = Cμ^¾k^{3/2}/(κy) and ω* = √k/(Cμ^¼κy).
- **WALE:** term-by-term match with Nicoud & Ducros, Cw = 0.6.
- **LES filter f1:** correct.
- **URANS filter form:** correct.
- **Buoyancy:** sign of G_b is correct.
- **Paraview output:** offsets are consistent.

---

## 4. Cross-module consistency

| Topic | SFLOW | NHFLOW | CFD |
|---|---|---|---|
| ν_t limiter | k-ε: Cμ·min(k²/ε, A263 k/S) ✗; k-ω: T31 | T31·k/S, A564 | k-ε always T31; k-ω only with T34 |
| ν_t used in production | same ν_t | k: EV; ε/ω: EV0 | k-ω: k limited, ω EV0; k-ε: limited for both |
| ω-production bound (EV0 floor) | none | none ✗ | yes (non-PLIC); none in PLIC ✗ |
| u⁺ log guard `MAX(0.01,…)` | none (k-ω) ✗ | k-ω yes, k-ε no ✗ | yes |
| ks ≤ 0 clamp | no ✗ | no ✗ | yes (`ks_val`) |
| Inflow k/ε/ω | none (Neumann) | ghosts written but eliminated ✗ | effectively 0 ✗ |
| k sink | explicit (k-ε) ✗ | implicit | implicit |
| Wall roughness parameter | B50 + Manning/26 | k-ε B50, k-ω B57, bed S20·S21 | `roughness::ks_val` (S28 aware) |
| 2D LES/SGS Δ | n/a | uses DYN ✗ | uses DYN ✗ |
| Free-surface ε/ω parameter | n/a | A567 (meaning differs by model) | T36/T37 (δ-weighted) |
| EV0 ghost update | n/a | missing ✗ | yes |
| Turbulence in restart | no | no | k/ε/ν_t yes, eddyv0 no |

Two refactors would remove most of these inconsistencies:

1. **One shared wall-law helper.** It should compute u⁺ with the guard and the ks clamp, and return τ, the k source/sink, ε* and ω*.
2. **One shared bed cf / u\* routine per module.** It should use the same constants as the momentum friction.

---

## 5. Design notes (deliberate choices, flagged for awareness)

- **Wall functions** use the rough-wall log law only, with no smooth-wall E/y⁺ branch. The wall distance is vertical (0.5·Δz), not bed-normal.
- **NHFLOW k/ε/ω** are solved in non-conservative advective form, not as D·k.
- **SFLOW turbulence** transport is not h-weighted. Momentum diffusion has no ∇ν_t term.
- **SFLOW c_εγ** defaults to 2.7, versus 3.6 in Rastogi & Rodi.
- **Default initial state.** CFD I13 = 0 and NHFLOW (without B60) both start with k = ε = ω = 0, and k-ε then relies on the limiter.
- **Two-phase CFD diffusion** uses kinematic ν + ν_t per cell rather than (μ + ρν_t)/ρ.
- **AMR** is disabled when turbulence is on, in both SFLOW and NHFLOW. Correct.
