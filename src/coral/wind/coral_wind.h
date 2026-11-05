/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Architect: Hans Bihs
--------------------------------------------------------------------*/
// REEF3D::CORAL: wind inflow for the aerodynamic loads (core layer only).
//
// Coordinates: CORAL's frame, z up, z = 0 the still water level (offshore) or the ground; heights are z.
// The wind frame: xi along the mean direction of propagation (direction 0 = +x, counterclockwise from above
// as the waves of CORAL), eta to the left looking downwind, z up.
//
// wind_field: the mean wind and its deterministic changes, plus an optional fluctuation field:
//   U(z, t) = U_ref(t) s(z) + gust(t) + event,   s(z) = (z/z_ref)^alpha  or  ln(z/z0)/ln(z_ref/z0)
//   veer: the direction turns with height, dtheta/dz [rad/m] (positive counterclockwise from above)
//   upflow: the inflow inclined by an angle (positive upwards)
//   series: hub speed, direction, vertical speed, shear exponent, linear shears and gust in time (uniform wind)
//   IEC 61400-1 (ed. 3) deterministic events of a turbine class: EOG, EDC, ECD, EWS (vertical and horizontal)
//   fluctuations: a wind_grid of u v w in the wind frame (turbulence), frozen and convected with U_hub
// wind_grid (coral_wind_grid.h) is also an inflow of its own: the full field of a TurbSim file.

#ifndef CORAL_WIND_H_
#define CORAL_WIND_H_

#include"coral/core/coral_types.h"
#include<memory>
#include<string>
#include<vector>

namespace coral
{
    // source of the wind velocity at a point (global frame)
    class wind_inflow
    {
    public:
        virtual ~wind_inflow() = default;
        virtual Vec3 velocity(const Vec3 &p, double t) const = 0;
        // mean speed at the reference (hub) height at time t (convection of frozen turbulence, controllers)
        virtual double hub_speed(double t) const = 0;
    };

    // IEC 61400-1 ed. 3 wind turbine class and the turbulence models
    struct iec_turbine
    {
        double V_ref = 50.0;      // reference wind speed (class I 50, II 42.5, III 37.5 m/s)
        double I_ref = 0.16;      // turbulence intensity at 15 m/s (A+ 0.18, A 0.16, B 0.14, C 0.12)
        double z_hub = 90.0;
        double D = 126.0;         // rotor diameter

        // cls "I", "II", "III"; cat "A+", "A", "B", "C"; false if unknown
        bool set_class(const std::string &cls, const std::string &cat);

        double Lambda1() const {return z_hub <= 60.0 ? 0.7*z_hub : 42.0;}     // longitudinal scale parameter
        double sigma_ntm(double V_hub) const {return I_ref*(0.75*V_hub + 5.6);}
        double sigma_etm(double V_hub) const;
        double sigma_ewm(double V_hub) const {return 0.11*V_hub;}
        double V_e50(double z) const;             // steady extreme wind model, recurrence period 50 years
        double V_e1(double z) const {return 0.8*V_e50(z);}
    };

    enum class iec_event {none, eog, edc_pos, edc_neg, ecd_pos, ecd_neg, ews_v_pos, ews_v_neg, ews_h_pos, ews_h_neg};
    // "eog", "edc+", "edc-", "ecd+", "ecd-", "ews_v+", "ews_v-", "ews_h+", "ews_h-"; none if unknown
    iec_event iec_event_from(const std::string&);

    class wind_grid;

    class wind_field : public wind_inflow
    {
    public:
        // U_ref at the height z_ref, direction of propagation [rad]
        wind_field(double U_ref, double z_ref, double direction=0.0);

        void set_power_law(double alpha) {alpha_ = alpha; z0_ = 0.0;}
        void set_log_law(double z0) {z0_ = z0;}
        void set_veer(double dtheta_dz) {veer_ = dtheta_dz;}
        // horizontal position of the hub (lateral shear, events, frozen turbulence)
        void set_origin(double x, double y) {ox_ = x; oy_ = y;}
        void set_upflow(double angle) {upflow_ = angle;}
        // smooth start 0.5 (1 - cos(pi t / t_ramp)) of the whole field
        void set_ramp(double t) {t_ramp_ = t;}

        // uniform wind in time (linear between the rows): hub speed, direction (rad, counterclockwise), vertical
        // speed, shear exponent, horizontal and vertical linear shear (per L_ref), gust speed; replaces U_ref,
        // the direction and alpha
        struct series_row {double t, U, dir, w, alpha, hshear, vshear, gust;};
        void set_series(const std::vector<series_row> &rows, double L_ref);
        // InflowWind's uniform wind file (.hh/.wnd): t U dir(deg, positive clockwise from above) [w] [h-shear]
        // [alpha] [lin. v-shear] [gust] [upflow]; comment lines start with '!' or '#'
        static bool read_uniform(const std::string &file, std::vector<series_row> &rows, std::string &msg);

        void set_event(iec_event e, const iec_turbine &tc, double t0) {event_ = e; iec_ = tc; t_event_ = t0;}

        // fluctuations u v w in the wind frame (has_mean false), convected with the hub speed
        void set_turbulence(std::shared_ptr<const wind_grid> g) {turb_ = std::move(g);}

        Vec3 velocity(const Vec3 &p, double t) const override;
        double hub_speed(double t) const override;

        // the mean speed (along the local direction) at height z and lateral position eta, the local direction
        double speed(double eta, double z, double t) const;
        double direction(double z, double t) const;

        double z_ref() const {return z_ref_;}
        const iec_turbine &iec() const {return iec_;}

    private:
        series_row row(double t) const;
        double shear(double z, double alpha) const;
        double ramp(double t) const;

        double U_ref_, z_ref_, dir_;
        double ox_ = 0.0, oy_ = 0.0, alpha_ = 0.0, z0_ = 0.0, veer_ = 0.0, upflow_ = 0.0, t_ramp_ = 0.0, L_ref_ = 1.0;
        std::vector<series_row> series_;
        iec_event event_ = iec_event::none;
        iec_turbine iec_;
        double t_event_ = 0.0;
        std::shared_ptr<const wind_grid> turb_;
    };
}

#endif
