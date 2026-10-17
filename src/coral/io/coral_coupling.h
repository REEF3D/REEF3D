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
// REEF3D::CORAL: the CORAL side of the coupling with a resolved fluid, without the fluid solver: the time
// stepping of the model in step with the fluid, the bodies driven by the fluid (fluid_body, two-way and
// one-way), the distributed loads (fluid_member), the sampled kinematics (fluid_kinematics), and a
// journal of all inputs from the fluid (fluid_coupling record <file>).
//
// REEF3D (coral/reef3d/coral_reef3d) feeds it with the loads and samples of its fluid solver; the
// standalone CORAL replays a journal (CORAL coral.txt replay <journal>) and repeats the coupled run
// without the fluid: a different PTO, controller or structure on the same fluid input (one-way coupling;
// for two-way bodies the replay only repeats the recorded run).
//
// Positions passed in and out are in the frame of the fluid solver (CORAL's frame plus fluid_coupling
// origin); loads, velocities and rotations are the same in both frames.

#ifndef CORAL_COUPLING_H_
#define CORAL_COUPLING_H_

#include"coral/io/coral_input.h"
#include"coral/io/coral_output.h"
#include"coral/core/coral_timeint.h"
#include"coral/core/coral_track.h"
#include<fstream>
#include<memory>
#include<ostream>
#include<string>
#include<vector>

namespace coral
{
    class coupling
    {
    public:
        // log: messages (only the writing rank should pass a real stream); out: whether this instance writes
        // the time series and the journal
        coupling(coral_case&, std::ostream &log, bool out);
        ~coupling();

        // start of the journal (the file of fluid_coupling record, if any, in dir)
        void record(const std::string &dir);

        // set-up at t0 in dir (static equilibrium and modes if asked, consistent accelerations, output)
        void ini(double t0, const std::string &dir);
        // after the fluid step to t1; with only one-way loads CORAL follows the newest loads (returns the
        // time reached); nothing when two-way bodies have stepped already
        double step_to(double t1);
        // two-way bodies: the CORAL step of the fluid step from t0 to t1
        void advance(double t0, double t1);

        // REEF3D 6DOF body n (see 6DOF_external.h)
        bool drives(int n) const;
        bool one_way(int n) const;
        bool predictor() const {return c_.fluid_predictor;}
        bool strips(int n, Vec3 &a, Vec3 &e, double &L, int &ns) const;
        void connect(int n, double mass, const Vec3 &c, const Mat3 &R, const Mat3 &I);
        void fluid_loads(double t, const std::vector<Vec6> &F, const std::vector<Mat6> &A, bool provisional);
        // strip forces of body n summed over the ranks
        void fluid_strips(int n, double t, const VecX &F);
        void predict(int n, double t, Vec3 &c, Mat3 &R, Vec3 &v, Vec3 &w) const;
        void state(int n, double t, Vec3 &c, Mat3 &R, Vec3 &v, Vec3 &w) const;

        // sampled kinematics: the points now (fluid frame), then the sample at time t (eta: z of the surface
        // in the fluid frame; inside 0/1)
        bool has_kinematics() const {return c_.fluid_kin != nullptr;}
        void kinematics_points(std::vector<Vec3> &pts) const;
        void kinematics(double t, const std::vector<Vec3> &u, const std::vector<double> &eta, const std::vector<int> &in);

        double time() const {return t_;}
        long steps() const {return nstep_;}
        bool started() const {return fluid_started_;}
        coral_case &model_case() {return c_;}

        // replays the journal file on the case, output to dir; returns the number of events
        static long replay(coral_case&, const std::string &journal, const std::string &dir, std::ostream &log);

    private:
        void advance_to(double t1);
        int fluid_index(int n) const;
        int member_index(int n) const;
        void jrec(int tag, const std::vector<double> &v);

        coral_case &c_;
        std::ostream &log_;
        bool out_on_;
        std::unique_ptr<timeint> ti_;
        std::unique_ptr<time_series> out_;
        body_track track_;
        double t_ = 0.0, next_out_ = 0.0;
        long nstep_ = 0;
        bool fluid_started_ = false;
        bool twoway_ = false;
        double t_load_ = -1.0e300;
        std::vector<VecX> F_ref_;
        std::ofstream jout_;
    };
}

#endif
