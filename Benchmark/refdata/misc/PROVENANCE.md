# Misc validation reference data: provenance

Collected 2026-10-03 using WebSearch and WebFetch only. The shell had no route to these hosts: curl, pip and gh were all refused by the proxy.

**Caveat on the fetch path:** WebFetch passes each page through a summarising model, and that model makes things up. In this session it invented two tables that were not in the content it received:
- A K1 table "from" HEC-18 FHWA-HIF-12-003. The PDF had been truncated after Chapter 2.
- A Kramer et al. (2021) sphere table. A second fetch gave different numbers and a third said the table was absent.

Neither appears in any file. To guard against this, every number kept in these files was handled in one of two ways:
- (a) Fetched twice with different prompts, with identical results (flagged "x2").
- (b) Taken from raw extracted text or code (verbatim strings).

Everything else is marked UNVERIFIED in the file itself.

| File | Content | Source (fetched) | Cites | Confidence |
|---|---|---|---|---|
| `hulme_hemisphere_heave.txt` | 26-row table: kA, A = a33/(ρV), B, for kA 0 to 10 | raw.githubusercontent.com/LHEEA/Hulme_Heaving_Sphere/master/comparison/HulmeResult.dat (x2, identical) | Hulme 1982 JFM 121:443 | **B+**. This is a third party's transcription of Hulme's table. The added-mass normalisation ρ(2/3)πa³ and the kA→0 limit 0.830951 are confirmed by arXiv math/0110302. The damping normalisation b33/(ρVω) is **not** confirmed verbatim, because the LHEEA code labels it "rho g V omega" (a typo). |
| `ito1977_free_decay.txt` | Ito tests a **2D horizontal cylinder**, not a sphere: r = 0.0762 m, ρb = 500, d = 1.22 m, initial offset 0.02454 m, NWT 10 × 1.8 m | Bihs et al. MARINE 2017 (REEF3D paper, upcommons.upc.edu); MIT DSpace metadata | Ito 1977 MIT MS thesis | **C**. The text values come from REEF3D's own paper. The thesis PDF could not be read (HTTP 405). There is no decay data, which exists only as a plot. |
| `cylinder_re100_unconfined.txt` | Williamson fit St = A/Re + B + C·Re, with A = -3.3265, B = 0.1816, C = 1.6e-4, for 50 ≤ Re ≤ 180. Re = 100 table of Cd, St, CL'(rms) and more for 8 studies, with their domains. | Beaudan & Moin 1994 Stanford TF-62 Eq. (1); Basilisk cylinder-strouhal.c; Qu et al. 2013 JFS 39 Table 1 (Chalmers full text, x2) | Williamson 1988 / 1991; Park et al. 1998; Posdziech & Grundmann 2007 and others | **A-** for the table. The fit's range is as Beaudan & Moin state it. Lift is given as **rms**, not amplitude. |
| `schaefer_turek_2d2.txt` | Problem definition; FEATFLOW Q2/P1 tables of C_D, C_L and St per level and dt; 1996 intervals | FEATFLOW 2D-2 page (wwwold.mathematik.tu-dortmund.de, x2); 1996 intervals from github AdebanjiAdelowo/fem-cylinder-flow README | Schäfer & Turek 1996 Table 4 | FEATFLOW tables **A-**. The 1996 intervals (C_D,max 3.22–3.24, C_L,max 0.99–1.01, St 0.295–0.305, Δp 2.46–2.50) are **UNVERIFIED secondary**: that repo says itself they were not re-checked. The original PDF returned 404. |
| `pier_scour.txt` | HEC-18 equation (FHWA-HRT-16-045 Fig. 3 alt text); K1, K2, K3 and limits (HEC-RAS 6.5 tech ref); Sumer wave scour S/D = 1.3{1 − exp[−0.03(KC−6)]}, KC ≥ 6; steady S/D = 1.3, σ = 0.7; Melville 1997 equations 9–17 and Tables 5–6 | fhwa.dot.gov, hec.usace.army.mil, mdpi JMSE 12:963 Eq. (6), Mostafa & Agamy 2011, Frontiers 2024, ponce.sdsu.edu/bridge_scour_1.html | HEC-18 5th ed.; Sumer et al. 1992; Breusers et al. 1977; Melville 1997 | **C** (verbatim equations). HEC-18 5th-ed equation **number** was not verified. The K tables come from HEC-RAS, which still carries K4 from the 4th edition. Melville Table 6 numbers come from one character-exact fetch. |
| `swe_analytic_solutions.txt` | Stoker, Ritter, Thacker 1D, and two Thacker 2D cases, with SWASHES parameters; solitary wave sech² profile, γ = √(3H/4), Xs, and c = √(g(d+H)) | arXiv 1110.0288v7 (SWASHES) raw pdf text; NOAA PAGEOPH 2008 + NCTR page; Basilisk beach.c | Ritter 1892, Stoker 1957, Thacker 1981, Synolakis 1987 | **A-**. Rebuilt from raw line-broken pdf text. The h-region conditions lost in extraction are marked [inferred]. SWASHES has no solitary wave. |

## Not obtained
- Hulme 1982 original Table: surge values, and explicit verbatim damping normalisation.
- Any **sphere** free-decay data attributed to Ito (1977). Ito's thesis is about a cylinder, and no tabulated decay curve, period or damping was found. Kramer et al. 2021 (Energies 14:269) holds a public sphere decay dataset, but its numbers were not transcribed because the fetches were inconsistent.
- Williamson 1988 paper text, which was blocked by ADS robots. Its coefficients were confirmed only via secondary sources.
- Lift **amplitude** at Re = 100 for Park et al. 1998. Only the rms value, 0.235, was found.
- Schäfer & Turek 1996 original PDF (404 and archive blocked), the Δp time-instant definition, and the FEATFLOW Δp table for 2D-2.
- HEC-18 5th edition (HIF-12-003) body text, equation numbers and Table 7.x. The PDF was truncated in the fetch.
- Melville & Sutherland 1988 individual K-factor definitions.
- Stoker and Ritter h-interval labels as printed. They were inferred from the u-intervals in the same extraction.
