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
// REEF3D::CORAL: linear potential-flow loads on the bodies of a hydrodynamic database (Cummins
// equation), the engineering hydrodynamics of WEC-Sim / HydroDyn type:
//   F = rho g V e_z - C xi + F_exc(t) - A_inf xi'' - int K(t - tau) xi'(tau) dtau - B_lin xi' - Q |xi'| xi'
// xi = [p - p0; theta] of the reference point of each database body (theta = log(R R0^T)), forces at
// the reference point in the global frame. Buoyancy acts at the reference point (the hydrostatic
// matrix of WAMIT-type databases holds the buoyancy terms about that point); gravity is a separate
// load at the centre of gravity. Unknown volume: buoyancy = weight of the body (equilibrium).
// Radiation memory by convolution with the committed velocity history or by state-space.
// Excitation from the incident waves (wave_field) and X(w, beta) of the database.

#ifndef CORAL_HYDRO_LOAD_H_
#define CORAL_HYDRO_LOAD_H_

#include"coral/core/coral_load.h"
#include"coral/core/coral_state.h"
#include"coral/hydro/coral_hydro_db.h"
#include"coral/hydro/coral_radiation.h"
#include"coral/hydro/coral_waves.h"
#include<memory>
#include<string>
#include<vector>

namespace coral
{
    struct hydro_settings
    {
        bool hydrostatics = true;
        bool buoyancy = true;
        bool excitation = true;
        bool radiation = true;
        bool state_space = false;     // radiation memory as states (default: convolution)
        double t_memory = 20.0;       // length of the radiation memory [s]
        double dt_memory = 0.01;      // sampling of K(t) [s]
        double ss_tol = 1.0e-4;       // state-space: discarded Hankel singular values
        int ss_order = 12;            // state-space: max order per DOF pair
        MatX linear_damping;          // 6N x 6N (empty: none)
        MatX quadratic_drag;          // 6N x 6N, F_i = -Q_ij |xi'_j| xi'_j (empty: none)
    };

    class hydro_load : public load
    {
    public:
        // bodies[i]: model body of database body i
        hydro_load(const model&, std::shared_ptr<const hydro_db>, const std::vector<int> &bodies,
                   std::shared_ptr<const wave_field>, const hydro_settings& = hydro_settings(),
                   const std::string &name="hydro");

        std::string name() const override {return name_;}
        coupling_kind coupling() const override {return coupling_kind::local;}
        void coupled_bodies(std::vector<int> &b) const override {b.insert(b.end(), map_.begin(), map_.end());}
        void initialize(const model&, double t0) override;
        void eval(const model&, double t, VecX &f) const override;
        void commit(const model&, double t, double dt) override;
        // restart from a checkpoint (coral_state.h)
        void save_state(std::vector<double>&) const override;
        bool load_state(const std::vector<double>&, size_t &k) override;
        void channel_names(std::vector<std::string>&) const override;
        void channel_values(const model&, double t, std::vector<double>&) const override;

        // 6N vectors of the database bodies in the current state of the model
        VecX displacement(const model&) const;
        VecX velocity(const model&) const;
        VecX acceleration(const model&) const;

        // force parts (6N) for the current state at time t
        VecX hydrostatic_force(const model&) const;
        VecX excitation_force(double t) const;
        VecX radiation_force(const model&, double t) const;
        VecX drag_force(const model&) const;

        const MatX &ainf() const {return Ainf_;}
        const retardation &memory() const {return K_;}
        const radiation_ss &state_space() const {return ss_;}

    private:
        VecX memory_force(const VecX &v1, double t) const;

        std::string name_;
        std::shared_ptr<const hydro_db> db_;
        std::shared_ptr<const wave_field> waves_;
        hydro_settings set_;
        std::vector<int> map_;
        std::vector<Vec3> s_, p0_;        // reference points (body frame), at equilibrium (global)
        std::vector<Mat3> R0_;
        std::vector<double> buoy_;
        MatX Ainf_;
        retardation K_;
        radiation_ss ss_;
        std::vector<VecXc> Xc_;           // excitation of the wave components, a X exp(i phase)

        // committed history (convolution) or states (state-space)
        std::vector<double> th_;
        std::vector<VecX> vh_;
        std::vector<Eigen::VectorXcd> xs_;
        VecX v_last_;
        double t_last_ = 0.0;
        bool started_ = false;

        // caches for the Newton iterations at one time
        mutable double cache_t_ = -1.0e300;
        mutable size_t cache_n_ = 0;
        mutable VecX cache_hist_;
        mutable double exc_t_ = -1.0e300;
        mutable VecX exc_f_;
    };
}

#endif
