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
// REEF3D::CORAL: loads of a resolved fluid (REEF3D FNPF; NHFLOW and CFD to follow) on rigid bodies of
// the model. The fluid solver integrates the pressure on the body surfaces and hands the loads over once
// per fluid time step (set_sample): per body the force and the moment about the CoG (global frame,
// without gravity, which is CORAL's gravity load) and, where the fluid solver has one, the instantaneous
// added mass A (6N x 6N, global, moments about the CoGs, block (i,j): load on body i from the
// accelerations of body j). Between the samples the load is a polynomial in time through the last
// order + 1 samples (Lagrange; order 1: linear, the default; 0: the last sample):
//   committed samples: one per fluid time level t_n, at the state the fluid and CORAL share there
//   provisional sample: the fluid load at a later stage time t_p of the same fluid step (computed with
//     a predicted body pose), replaced by the next committed sample
// With a provisional sample the CORAL step interpolates instead of extrapolating (predictor-corrector
// coupling without a re-step of CORAL). The added mass acts on the accelerations of the current
// (trial) state:
//   F_total(t) = F(t) - A [a_1; alpha_1; ...; a_N; alpha_N]
// a: acceleration of the CoG, alpha: angular acceleration (global). F is the part of the fluid load
// without the acceleration term (FNPF: psi_0 = phi_t or chi without the body accelerations), so the
// added mass is implicit in CORAL's Newton iterations and the coupling is stable for light bodies.
// The loads act at the CoG of each body, or at a point of application given by set_attachment.

#ifndef CORAL_FLUID_LOAD_H_
#define CORAL_FLUID_LOAD_H_

#include"coral/core/coral_load.h"
#include<string>
#include<vector>

namespace coral
{
    class fluid_load : public load
    {
    public:
        fluid_load(const model&, const std::vector<int> &bodies, const std::string &name="fluid");

        std::string name() const override {return name_;}
        coupling_kind coupling() const override {return coupling_kind::local;}
        void coupled_bodies(std::vector<int> &b) const override {b.insert(b.end(), map_.begin(), map_.end());}
        void eval(const model&, double t, VecX &f) const override;
        void channel_names(std::vector<std::string>&) const override;
        void channel_values(const model&, double t, std::vector<double>&) const override;

        int nbody() const {return (int)map_.size();}
        const std::vector<int> &bodies() const {return map_;}

        // fluid loads at the fluid time t (6N) and the added mass (6N x 6N; empty: keep the last);
        // provisional: a stage of the current fluid step (see above)
        void set_sample(double t, const VecX &F, const MatX &A=MatX(), bool provisional=false);
        bool has_sample() const {return !ts_.empty();}
        double sample_time() const {return ts_.empty() ? 0.0 : ts_.back();}
        bool has_provisional() const {return prov_;}

        // the extrapolated fluid load without the added-mass part, and the added-mass part of the state
        VecX force(double t) const;
        VecX added_mass_force(const model&) const;
        VecX acceleration(const model&) const;
        const MatX &added_mass() const {return A_;}
        // point of application of the load of fluid body k (global, now; default its CoG): the loads are
        // about that point (one-way loads of a fluid body on a CORAL body that is not the fluid body itself)
        void set_attachment(const model&, int k, const Vec3 &p);

        int order = 1;      // polynomial order in time through the samples (0, 1, 2)

    private:
        std::string name_;
        std::vector<int> map_;
        std::vector<std::string> bname_;
        std::vector<Vec3> s_;           // points of application, body frame
        std::vector<double> ts_;        // sample times (the committed ones, then the provisional one)
        std::vector<VecX> Fs_;
        bool prov_ = false;
        MatX A_;
    };
}

#endif
