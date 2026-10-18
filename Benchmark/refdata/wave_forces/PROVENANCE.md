# Wave force reference data: provenance

Copied unchanged (numeric rows only, the one-line header of the original dropped) from the REEF3D
repository, branch `hans_dev`, commit `ff1bb6819`, folder `Tutorials/REEF3D_CFD/`. These are the
measured series the user guide tutorials are compared with.

| File | Content | Original file | Experiment | Notes |
|---|---|---|---|---|
| `chen2014_eta_x5.txt` (286 rows: t [s], eta [m]) | surface elevation at x = 5 m | `11_14 Non-Breaking Wave Forces/expt_data_wsf.txt` | Chen et al. (2014), Ocean Eng. 88 | time base of the tutorial comparison; the checker fits the lag |
| `chen2014_force.txt` (274 rows: t [s], Fx [N]) | inline force on the D = 0.25 m cylinder | `11_14 .../expt_data_force.txt` | Chen et al. (2014) | full cylinder |
| `irschik2002_force.txt` (106 rows: t [s], Fx [N]) | breaking-wave slamming force on the D = 0.7 m pile, one impact | `11_15 Breaking Wave Forces/expt_data_force.txt` (header `ti2 force`) | Irschik, Sparboom & Oumeraci (2002), GWK | full pile |
| `irschik2002_eta_wg2.txt` (177 rows: t [s], eta [m]) | surface elevation at gauge "wg2" | `11_15 .../expt_data_wsf.txt` (header `wg2 wg_y`) | Irschik et al. (2002) | gauge position not documented, not used by the checker |

Confidence: the files are the data the REEF3D group has used in its publications (Bihs et al. 2016,
Computers & Fluids; Kamath et al. 2016, Ocean Engineering). The irregular time steps
(0.003-0.12 s) indicate that all four series were digitised from figures, not exported from the
instrument records; treat them like the other confidence-B data sets (digitised, secondary).
