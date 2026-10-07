# Wave validation reference data: provenance

Collected 2026-10-03 with WebSearch/WebFetch only (shell curl had no route to these hosts).
Every number in these files was copied from a page fetched in that session. Nothing was recalled from memory or
read off a plot by eye. Each file's `#` header gives its source URL.

**Caveat on the fetch path:** WebFetch passes pages through a summarising model. Each numeric file was fetched as
raw text and asked for verbatim output. As a spot check, `gauge-9`, `section-5` and `t-50` were fetched a second time
with a different prompt, and both fetches matched exactly.

## Confidence scale
- **A**: original tabulated experimental values.
- **B**: experimental values digitised from a published figure by a third party and shipped with a model test case.
- **C**: geometry or parameters quoted from text.

| File(s) | Content | Source | Cites | Confidence |
|---|---|---|---|---|
| `beji_battjes_bar/WG4..WG11_eta_timeseries.txt` (8 files, 34-35 rows x 2 cols: t[s], eta[cm]) | Measured eta(t), T = 2.02 s, about 33-39 s | `http://basilisk.fr/src/test/gauge-N?raw`, used by `http://basilisk.fr/src/test/bar.c` | Beji & Battjes 1993, Luth et al. 1994, Dingemans 1994, Yamazaki et al. 2009 | **B**. Irregular dt (about 0.18 s) and quantised values point to digitising, probably from Yamazaki et al. 2009 Fig. 4. The raw file holds the series twice, and the duplicate copy was dropped. Model forcing amplitude was tuned to WG4, and the time origin is arbitrary. |
| `beji_battjes_bar/geometry.txt`, `gauges.txt` | Bar geometry and gauge x | Beji & Battjes 1993 PDF (serdarbeji.com), basilisk bar.c, NHWAVE paper (udel), arXiv 1004.3436 | | **C**. h0 0.4, 1:20 up, crest 0.1 m deep and 2 m long, 1:10 down, toe 6 m from the paddle. Confirmed by 2 or more sources. |
| `berkhoff_shoal/section{2,3,5,7}_H_over_H0.txt` (26, 26, 26, 19 rows x 2 cols) | Normalised wave height along sections 2, 3, 5, 7 | `https://basilisk.fr/src/examples/section-N?raw`, used by `shoal.c` | Berkhoff et al. 1982 via Lannes & Marche 2014 | **B** (digitised). Sections 1, 4, 6 and 8 were not found. |
| `berkhoff_shoal/geometry.txt` | Bathymetry formula, T = 1 s, H = 4.64 cm | basilisk shoal.c code, NHWAVE paper, FUNWAVE manual, XBeach-NH report | | **C**. The 5.82 vs 5.84 offset is unresolved. |
| `synolakis_runup/profile_Hd0.28_tNN.txt` (12 files, t = 10..65, 37-57 rows x 2 cols: x/d, eta/d) | Surface profiles, breaking case H/d = 0.28 | `https://basilisk.fr/src/test/t-NN?raw`, used by `beach.c` | Synolakis 1987 | **B** (digitised). Exact duplicate pairs were removed and rows sorted by x. The x origin is the shoreline, positive onshore. The t origin follows beach.c's definition. |
| `synolakis_runup/runup_law_and_setup.txt` | R = 2.831 sqrt(cot b) H^(5/4), valid for sqrt(H) <= 0.288 tan b; slope 1:19.85; toe 14.95 m | Synolakis et al. 2008 PAGEOPH (nctr.pmel.noaa.gov), NOAA NCTR benchmark page | | **C** |
| `whalin_shoal/geometry.txt` | Shoal formula, tank 25.6 x 6.096 m, depths 0.4572 to 0.1524 m | Celeris paper arXiv 1611.05984 | Whalin 1971 | **C**. The slope term's sign is garbled in the PDF extraction and must be checked. There is no measured data. |
| `ting_kirby_breaking/setup.txt` | Slope 1:35, h = 0.4 m, H and T for both cases, break points | Derakhti et al. (udel), ERDC air-water-vv docs, ICCE paper | Ting & Kirby 1994 | **C**. There is no measured data. |

## Not obtained
- Original Dingemans (1994) gauge records sampled at 0.05 s. Repos that may hold them could not be listed: FUNWAVE `BENCHMARK_FUNWAVE/car_luth_AC` and ERDC `air-water-vv`. GitHub tree pages are blocked by robots.txt, and the API, codeload and jsDelivr returned 403.
- Berkhoff sections 1, 4, 6 and 8.
- Whalin harmonic amplitudes for all cases.
- Ting & Kirby measured H and set-up.
- Synolakis H/d = 0.0185 and 0.3 profiles. They exist only as NOAA `.xls` files, which are binary and unreadable here.
- Measured run-up vs H/d tables.
