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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"fnpf_prof.h"
#include"lexer.h"
#include"ghostcell.h"
#include<iostream>
#include<iomanip>

namespace fnpf_prof
{
    double t[N] = {0.0};
    double tstep = 0.0;
    int steps = 0;
    
    namespace
    {
        int reported = -1;
        
        const char *label[N] =
        {
            "wave generation (precalc)", "RK stages (total)", "printer + LAGOON poll", "time step",
            "kfsfbc + dfsfbc", "fsfdisc (WENO, metrics)", "sigma_update", "fsfbc_sig + bedbc_sig",
            "Fi halo exchange (start7V)", "Laplace solve", "fsfwvel", "velcalc_sig (+U,V,W halos)",
            "eta/Fifsf relaxation", "damping (vb: breaking, coastline)", "wetdry", "coastline", "breaking", "inflow_fnpf",
            "slice halo exchange (gcsl)", "body (6DOF/ALE)"
        };
    }
    
    void step_end(double steptime)
    {
        tstep += steptime;
        ++steps;
    }
    
    void report(lexer *p, ghostcell *pgc, bool final)
    {
        if(steps==0 || steps==reported)
        return;
        
        reported = steps;
        
        // RK parts timed at the top level, the rest of the stage and of the step
        double rkparts = 0.0;
        for(int k=fsfbc; k<N; ++k)
        rkparts += t[k];
        
        const double rkrest = t[rk_total] - rkparts;
        const double steprest = tstep - t[wavegen] - t[rk_total] - t[print] - t[timestep];
        
        const double ns = double(steps);
        const double nr = double(p->mpi_size);
        
        auto mean = [&](double v) {return pgc->globalsum(v)/nr;};
        
        double avg[N], mx[N];
        for(int k=0; k<N; ++k)
        {
        avg[k] = mean(t[k]);
        mx[k] = pgc->globalmax(t[k]);
        }
        
        const double avg_step = mean(tstep), mx_step = pgc->globalmax(tstep);
        const double avg_rkrest = mean(rkrest), mx_rkrest = pgc->globalmax(rkrest);
        const double avg_rest = mean(steprest), mx_rest = pgc->globalmax(steprest);
        
        if(p->mpirank!=0)
        return;
        
        const std::streamsize pr = std::cout.precision();
        const std::ios::fmtflags fl = std::cout.flags();
        
        auto row = [&](const char *name, double a, double m, bool indent)
        {
            std::cout<<(indent ? "    " : "  ")<<std::left<<std::setw(indent ? 30 : 32)<<name<<std::right
                     <<std::fixed<<std::setprecision(1)
                     <<std::setw(10)<<1000.0*a/ns
                     <<std::setw(10)<<1000.0*m/ns
                     <<std::setw(8)<<100.0*a/avg_step<<" %"<<std::endl;
        };
        
        std::cout<<std::endl<<"FNPF step profile, "<<steps<<" steps"<<(final ? " (end of run)" : "")
                 <<"  [ms/step: mean over ranks, slowest rank, % of the mean step]"<<std::endl;
        
        row("wave generation (precalc)",avg[wavegen],mx[wavegen],false);
        row("RK stages (total)",avg[rk_total],mx[rk_total],false);
        
        for(int k=fsfbc; k<N; ++k)
        row(label[k],avg[k],mx[k],true);
        
        row("rest of the stages",avg_rkrest,mx_rkrest,true);
        row("printer + LAGOON poll",avg[print],mx[print],false);
        row("time step",avg[timestep],mx[timestep],false);
        row("rest of the step",avg_rest,mx_rest,false);
        row("step",avg_step,mx_step,false);
        std::cout<<std::endl;
        
        std::cout.precision(pr);
        std::cout.flags(fl);
    }
}
