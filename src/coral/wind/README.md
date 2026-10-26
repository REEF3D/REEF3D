# CORAL: wind/

Wind inflow.

## Content

* `coral_wind.h/.cpp`: `wind_inflow` (velocity at a point and time, hub speed); `wind_field`: steady wind with power-law or logarithmic shear, veer, upflow, the hub's horizontal position and a start ramp; uniform wind in time (hub speed, direction, vertical speed, shear exponent, linear horizontal and vertical shear, gust; InflowWind's uniform wind file); IEC 61400-1 ed. 3 turbine classes (`iec_turbine`: V_ref, I_ref, Lambda_1, NTM, ETM and EWM sigma_1, V_e50, V_e1) and the deterministic events EOG, EDC+-, ECD+-, EWS vertical and horizontal +-; an optional fluctuation grid on the mean wind
* `coral_wind_grid.h/.cpp`: `wind_grid`: u v w on a lateral-vertical grid in time (turbulence box) with tower points, frozen and convected with the hub speed, bilinear in space and linear in time, periodic or held; TurbSim full fields (.bts) read and written; IEC Kaimal turbulence with the exponential coherence (u, optionally v and w) by Veers' method (Cholesky of the cross-spectral matrix per frequency, diagonal where the coherence of neighbours has died out, a Fourier series of 2^n levels, portable seeds), optional scaling to sigma at the hub
* `coral_wind_probe.h/.cpp`: output channels of the inflow at a point fixed in space or on a body

Keywords (`io/coral_input.h`): `rho_air`, `wind steady`, `wind file`, `wind iec`, `wind event`, `wind turbulence kaimal`, `wind turbsim`, `wind_probe`.

## Planned

* Mann turbulence
* height-dependent turbulence intensity and other spectra (von Karman, API/NPD for offshore)
* wind from REEF3D's resolved air phase at the blade and tower points

## May include from

`core`
