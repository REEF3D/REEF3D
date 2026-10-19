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
// REEF3D::CORAL inside REEF3D (Z 40 1): runs the CORAL model of coral.txt in the case folder alongside the
// fluid solver (CFD, NHFLOW, FNPF), in step with the fluid time: each fluid step from simtime to
// simtime + dt is covered by CORAL steps of at most the dt of coral.txt. Every rank holds the same model
// and computes the same motion; rank 0 reads coral.txt and broadcasts it, and writes the time series to
// REEF3D_CORAL/. The CORAL side of the coupling is coral::coupling (coral/io/coral_coupling.h); this
// class adds REEF3D: MPI, the sampling of the fluid and the interface of the 6DOF bodies:
// * fluid_body: CORAL is the external dynamics of the 6DOF bodies (6DOF_external.h, FNPF two-way; FNPF
//   and NHFLOW one-way); two-way, the CORAL step of a fluid step happens inside the fluid step
// * fluid_member: the hull loads of a one-way 6DOF body in strips along a CORAL member
// * fluid_kinematics: the Morison members and mooring lines of coral.txt take the velocity and free
//   surface of the resolved fluid (FNPF, NHFLOW), sampled at their points at every fluid time level
//   (after the fluid step, or at its first stage with two-way bodies)
// * fluid_coupling record: a journal of all of it for the standalone CORAL (replay)

#ifndef CORAL_REEF3D_H_
#define CORAL_REEF3D_H_

#include"6DOF_external.h"
#include<memory>
#include<string>
#include<vector>

class lexer;
class ghostcell;
class coral_fluid_sampler;

namespace coral
{
    class coral_case;
    class coupling;
}

class coral_reef3d : public sixdof_external
{
public:
    coral_reef3d(lexer*, ghostcell*);
    ~coral_reef3d();

    // set up at p->simtime (consistent accelerations, first output)
    void ini(lexer*, ghostcell*);
    // advance from p->simtime to p->simtime + p->dt (nothing to do when the fluid bodies already did)
    void step(lexer*, ghostcell*);
    // sampling of the fluid for fluid_kinematics (set by the driver before ini, FNPF and NHFLOW)
    void set_sampler(coral_fluid_sampler*);

    // checkpoints (P 48, I 45): the state of CORAL, the same on all ranks (that of rank 0 is written);
    // restart before ini: ini continues the outputs and takes the state instead of setting up
    void ckpt_state(std::vector<double>&) const;
    void restart(const std::vector<double>&);
    // the time series and the journal on disk
    void flush();

    // external dynamics of the 6DOF bodies of the fluid solver (fluid_body, fluid_member)
    bool drives(int n) const override;
    bool one_way(int n) const override;
    bool strips(int n, Eigen::Vector3d &a, Eigen::Vector3d &e, double &L, int &ns) const override;
    void fluid_strips(int n, double t, const std::vector<double> &F) override;
    void connect(int n, double mass, const Eigen::Vector3d &c, const Eigen::Matrix3d &R, const Eigen::Matrix3d &I) override;
    bool predictor() const override;
    void fluid_loads(double t, const std::vector<vec6> &F, const std::vector<mat6> &A, bool provisional) override;
    void advance(double t0, double t1) override;
    void predict(int n, double t, Eigen::Vector3d &c, Eigen::Matrix3d &R, Eigen::Vector3d &v, Eigen::Vector3d &w) override;
    void state(int n, double t, Eigen::Vector3d &c, Eigen::Matrix3d &R, Eigen::Vector3d &v, Eigen::Vector3d &w) override;

private:
    void fail(const std::string&);
    bool load_state(const std::vector<double>&);
    void sample_kinematics(double t);   // the fluid at the points of the sampled kinematics at time t
    // several ranks: every rank advances the same model; its input from the fluid is broadcast from
    // rank 0 (bit for bit the same on all ranks), and after each step the state of the model is
    // compared over the ranks (a warning at the first difference)
    void bcast(double *v, int n);
    void check_ranks(const char *where);
    bool rank_diff_told_ = false;

    lexer *p_ = nullptr;
    ghostcell *pgc_ = nullptr;
    std::unique_ptr<coral::coral_case> case_;
    std::unique_ptr<coral::coupling> cp_;
    std::unique_ptr<coral_fluid_sampler> sampler_;
    bool sampler_told_ = false;
    std::vector<double> restart_state_;
    bool restart_ = false;
};

#endif
