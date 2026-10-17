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
// REEF3D::CORAL: incident waves as a sum of linear (Airy) components, for running CORAL without a fluid
// grid (excitation from the hydrodynamic database, Morison kinematics). Regular waves, JONSWAP and
// Pierson-Moskowitz spectra with seeded random phases, or user components; finite or infinite depth.
//
// Coordinates: z up, still water level z = 0, seabed z = -depth. Component i:
//   eta = a cos(theta),  theta = omega t - k (x cos(beta) + y sin(beta)) + phase
// (the excitation of the hydrodynamic database refers to the wave elevation at the origin).
// Above the still water level the kinematics of z = 0 are used.

#ifndef CORAL_WAVES_H_
#define CORAL_WAVES_H_

#include"coral/core/coral_types.h"
#include<vector>

namespace coral
{
    struct wave_component
    {
        double a = 0.0;        // amplitude [m]
        double w = 0.0;        // circular frequency [rad/s]
        double k = 0.0;        // wave number [1/m]
        double beta = 0.0;     // direction of propagation [rad], 0 = +x
        double phase = 0.0;    // [rad]
    };

    enum class wave_spectrum {jonswap, pierson_moskowitz};

    // Source of the wave kinematics of the Morison members and the mooring lines: CORAL's linear waves
    // (wave_field) or the resolved fluid of REEF3D sampled at the loads' points (sampled_kinematics).
    // A load that asks at one of its own points passes itself and the index of the point (owner, index),
    // so that a sampled source can return the values sampled there; the load lists these points, in
    // index order, in load::kinematics_points. Elevation: z of the free surface (CORAL frame).
    class wave_kinematics
    {
    public:
        virtual ~wave_kinematics() = default;
        virtual double elevation(double x, double y, double t) const = 0;
        // velocity u, acceleration du (local), dynamic pressure / rho at the point p
        virtual void kinematics(const Vec3 &p, double t, Vec3 &u, Vec3 &du, double &p_rho) const = 0;
        // the same at point index of owner (default: at p)
        virtual double elevation_at(const void*, int, const Vec3 &p, double t) const {return elevation(p(0), p(1), t);}
        virtual void kinematics_at(const void*, int, const Vec3 &p, double t, Vec3 &u, Vec3 &du, double &p_rho) const
        {kinematics(p, t, u, du, p_rho);}
    };

    class wave_field : public wave_kinematics
    {
    public:
        // depth <= 0: infinite depth
        explicit wave_field(double depth=-1.0, double g=9.81);

        void add_component(double a, double w, double beta=0.0, double phase=0.0);
        void add_regular(double H, double T, double beta=0.0, double phase=0.0);
        // n components between wmin and wmax (equal spacing, frequencies shifted randomly inside their
        // bin against the repetition of the series), random phases from seed; gamma for JONSWAP
        void add_spectrum(wave_spectrum, double Hs, double Tp, double gamma, double beta,
                          double wmin, double wmax, int n, unsigned seed);

        // smooth start: 0.5 (1 - cos(pi t / t_ramp)) for t < t_ramp
        void set_ramp(double t_ramp) {t_ramp_ = t_ramp;}
        double ramp(double t) const;

        double elevation(double x, double y, double t) const override;
        // velocity u, acceleration du (local, linear), dynamic pressure / rho at the point p
        void kinematics(const Vec3 &p, double t, Vec3 &u, Vec3 &du, double &p_rho) const override;

        const std::vector<wave_component> &components() const {return comp_;}
        double depth() const {return depth_;}
        double g() const {return g_;}
        // significant wave height of the components, 4 sqrt(m0)
        double hs() const;

        // spectral density S(w) [m^2 s] of the spectrum (normalized to Hs)
        static double spectrum(wave_spectrum, double w, double Hs, double Tp, double gamma);
        // dispersion relation w^2 = g k tanh(k d)
        static double wavenumber(double w, double depth, double g);

    private:
        double depth_, g_;
        double t_ramp_ = 0.0;
        std::vector<wave_component> comp_;
    };
}

#endif
