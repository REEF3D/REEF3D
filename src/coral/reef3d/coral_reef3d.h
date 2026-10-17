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
// and computes the same motion (no communication); rank 0 reads coral.txt and broadcasts it, and writes
// the time series to REEF3D_CORAL/.
// CORAL's own engineering hydrodynamics (Cummins, Morison, mooring, PTO) work as in the standalone
// CORAL. Bodies marked fluid_body in coral.txt are the 6DOF bodies of the fluid solver (FNPF so far):
// CORAL is their external dynamics (6DOF_external.h): the fluid hands over its loads and added mass,
// CORAL steps the whole model and returns the bodies' poses and velocities at the stage times of the
// fluid (cubic Hermite interpolation over the CORAL steps, coral/core/coral_track.h). The CORAL step of
// a fluid step then happens inside the fluid step; step() after the fluid step finds it done.

#ifndef CORAL_REEF3D_H_
#define CORAL_REEF3D_H_

#include"6DOF_external.h"
#include"coral/core/coral_track.h"
#include<memory>
#include<string>

class lexer;
class ghostcell;

namespace coral
{
    class coral_case;
    class timeint;
    class time_series;
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

    // external dynamics of the 6DOF bodies of the fluid solver (fluid_body in coral.txt)
    bool drives(int n) const override;
    void connect(int n, double mass, const Eigen::Vector3d &c, const Eigen::Matrix3d &R, const Eigen::Matrix3d &I) override;
    bool predictor() const override;
    void fluid_loads(double t, const std::vector<vec6> &F, const std::vector<mat6> &A, bool provisional) override;
    void advance(double t0, double t1) override;
    void predict(int n, double t, Eigen::Vector3d &c, Eigen::Matrix3d &R, Eigen::Vector3d &v, Eigen::Vector3d &w) override;
    void state(int n, double t, Eigen::Vector3d &c, Eigen::Matrix3d &R, Eigen::Vector3d &v, Eigen::Vector3d &w) override;

private:
    void fail(ghostcell*, const std::string&);
    void advance_to(double t1);
    int fluid_index(int n) const;       // index in the fluid load of 6DOF body n, -1: not driven

    lexer *p_ = nullptr;
    ghostcell *pgc_ = nullptr;
    coral::body_track track_;
    bool fluid_started_ = false;

    std::unique_ptr<coral::coral_case> case_;
    std::unique_ptr<coral::timeint> ti_;
    std::unique_ptr<coral::time_series> out_;
    double t_ = 0.0;
    double next_out_ = 0.0;
    long nstep_ = 0;
};

#endif
