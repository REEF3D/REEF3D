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

#ifndef STEP_REPEAT_H_
#define STEP_REPEAT_H_

//  N 66: the repeat of a time step.
//
//  The adaptive time step of a step comes from the velocities at its start.  When they jump within the
//  step (an impact, a closing gap), the step was too large for what it produced, and the criterion says
//  so only for the step after it.  With N 66 r the main loop keeps the state of the start of every step
//  in memory (reefrun_ckpt::mem_save: what a checkpoint holds, so a repeated step starts from exactly
//  the state the first try started from), and once the time step of the next step is known:
//
//      dt of the step > r * dt the criterion gives after it   ->   the state is taken back, the step
//                                                                   repeated with half its time step
//
//  and also when the velocities after the step are above N 61 (the run would stop),
//  at most MAXTRY times per step.  N 66 -k is the test of it: every k. step is taken back and done again
//  with the same time step, and the run has to be the run without N 66, bit for bit.
//  The outputs of a step are written after this, so a repeated step
//  leaves no trace in them.  Collective: the time steps are those of the run, the same on all ranks.

#include"lexer.h"
#include"reefrun_ckpt.h"
#include"reef3d_mpmd.h"
#include<string>
#include<iostream>
#include<iomanip>
#include<algorithm>
#include<mpi.h>

class step_repeat
{
public:
    // the start of a step, after the count: the state kept.  refuse: why the repeat is not there for
    // this run (empty: it is); key: the ratio r (N 66), 0 off; adaptive: the run has an adaptive time step
    void begin(lexer *p, reefrun_ckpt &ck, double key, bool adaptive, const std::string &refuse)
    {
        on = false;
        if(key==0.0 || off)
        return;

        std::string why = (adaptive || key<0.0) ? refuse : std::string("the time step is fixed");
        if(why.empty())
        ck.mem_save(why);
        if(!why.empty())
        {
            if(p->mpirank==0)
            std::cout<<"N 66: the repeat of a time step is not available here ("<<why<<"): not used"<<std::endl;
            off = true;
            return;
        }

        if(!told)
        {
            long long b = ck.mem_bytes(), bmax = 0;
            MPI_Reduce(&b,&bmax,1,MPI_LONG_LONG,MPI_MAX,0,reef3d_comm());
            if(p->mpirank==0 && key>0.0)
            std::cout<<"N 66 "<<key<<": a step is repeated with half the time step when the time step after it is below 1/"<<key
                     <<" of it; the state of the start of a step in memory, "<<std::setprecision(3)<<double(bmax)/1048576.0<<" MB on the largest rank"<<std::endl;
            if(p->mpirank==0 && key<0.0)
            std::cout<<"N 66 "<<key<<": the test of the repeat, every "<<std::max(1,int(-key))<<". step is taken back and done again with the same time step"
                     <<" (the run is that without N 66); the state of the start of a step in memory, "<<std::setprecision(3)<<double(bmax)/1048576.0<<" MB on the largest rank"<<std::endl;
        }
        told = true;

        r = key;
        on = true;
        tries = 0;
        dtry = p->dt;
    }

    // after the time step of the next step is set (p->dt): true when the state of the start of the
    // step is back and p->dt is half the time step of the try before: the step is to be done again
    // stop: the run would stop after this step (the velocities above N 61): repeated whatever the ratio
    bool again(lexer *p, reefrun_ckpt &ck, bool stop=false)
    {
        if(!on || tries>=MAXTRY)
        return false;

        const double dtnext = p->dt;
        const bool test = r<0.0;
        if(test && (tries>0 || p->count%std::max(1,int(-r))!=0))
        return false;
        if(!test && !stop && !(dtry > r*dtnext))
        return false;

        const int cnt = p->count;
        std::string why;
        if(!ck.mem_restore(why))
        {
            if(p->mpirank==0)
            std::cout<<std::endl<<"!!! N 66: the state of the start of step "<<cnt<<" could not be taken back ("<<why<<") !!!"<<std::endl<<std::endl;
            MPI_Abort(reef3d_comm(),-66);
        }

        const double dold = dtry;
        if(!test)
        dtry = 0.5*dtry;
        p->dt = dtry;
        ++tries;
        ++total;

        if(p->mpirank==0 && !test)
        std::cout<<"N 66: step "<<p->count<<" repeated with the time step "<<std::setprecision(5)<<dtry<<" (it was "<<dold
                 <<", the criterion after it gave "<<dtnext<<(stop ? ", velocities above N 61" : "")<<"); "<<total<<" repeats so far"<<std::endl;
        return true;
    }

    int repeats() const { return total; }
    // the state of this step is kept (begin)
    bool active() const { return on; }
    // the time step of the last try of this step (the one it was taken with)
    double dt_taken() const { return dtry; }

private:
    enum { MAXTRY=3 };
    bool on = false, off = false, told = false;
    double r = 0.0, dtry = 0.0;
    int tries = 0, total = 0;
};

#endif
