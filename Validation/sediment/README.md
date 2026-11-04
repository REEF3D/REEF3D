# Sediment transport validation (SFLOW, NHFLOW)

These are cases measured in the laboratory, used to validate:
- the waves in the sediment transport (S 110, S 11 6–9);
- the depth-averaged suspended load of SFLOW (S 120);
- the NHFLOW suspended load (S 12);
- the sediment supply at the inflow (S 129, S 130, S 124).

| folder | experiment | result |
|---|---|---|
| `01_vanrijn_trench` | van Rijn (1986) trench migration in a current, 3 trenches | NHFLOW RMS 1.9 – 3.3 cm (validated, after calibrating the grain roughness to the measured supply); SFLOW migrates too fast (RMS 5–6 cm) |
| `02_lip11d_1b` | LIP 11D test 1B, erosive breaker bar (Delta Flume) | no skill yet (BSS about 0), see the README for the next steps |

`tools/`:
- `mk_trench.py`, `mk_lip.py`: write control.txt, ctrl.txt and geo.dat; extra keys as `"S 21=30"`, and `"A 534=del"` removes a key.
- `run_case.sh`: runs DiveMESH and REEF3D in a run folder.
- `sed_eval.py`: scores (RMS / BSS against `refdata`).
- `flux.py`: transport along the flume.
- `vtp.py`, `vtu.py`: readers for the raw appended VTK files.

Needs Python 3 and numpy.
