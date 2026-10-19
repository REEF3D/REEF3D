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
// REEF3D::CORAL standalone executable: runs a CORAL model without a fluid grid.
//   CORAL [input file]                       (default coral.txt)
//   CORAL <input file> replay <journal>      the fluid input recorded by a coupled REEF3D run
//                                            (fluid_coupling record), output in this folder

#include"coral/core/coral_version.h"
#include"coral/core/coral_timeint.h"
#include"coral/io/coral_input.h"
#include"coral/io/coral_output.h"
#include"coral/io/coral_vtp.h"
#include"coral/core/coral_static.h"
#include"coral/linear/coral_modes.h"
#include"coral/io/coral_modes_out.h"
#include"coral/io/coral_coupling.h"
#include<chrono>
#include<sys/stat.h>
#include<cmath>
#include<cstdio>
#include<exception>
#include<iostream>
#include<string>

int main(int argc, char **argv)
{
    using namespace coral;
    const std::string file = argc > 1 ? argv[1] : "coral.txt";
    std::printf("%s\n", banner());
    try
    {
        coral_case c;
        c.read(file);
        if(argc > 3 && std::string(argv[2]) == "replay")
        {
            const auto w0 = std::chrono::steady_clock::now();
            const long nev = coupling::replay(c, argv[3], ".", std::cout);
            const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - w0).count();
            std::printf("replay of %s: %ld events, %.2f s\n", argv[3], nev, wall);
            return 0;
        }
        std::printf("input %s: %d bodies, %d constraint equations, %d loads, %s waves\n", file.c_str(),
                    c.m.nbody(), c.m.ncon(), (int)c.m.load_list.size(), c.waves ? std::to_string(c.waves->components().size()).c_str() : "no");
        if(c.fluid)
            std::printf("fluid_body: no resolved fluid in the standalone run, the fluid loads are zero\n");
        const run_settings &r = c.run;
        const auto w0 = std::chrono::steady_clock::now();
        if(r.static_ini)
        {
            const static_result sr = static_solve(c.m, r.t_start);
            std::printf("static equilibrium: %s after %d iterations, max |residual| %.2e\n",
                        sr.converged ? "converged" : "NOT converged", sr.iterations, sr.residual);
        }
        if(r.modes > 0)
        {
            const mode_result mr = modes(c.m, r.t_start, r.modes, r.modes_damped);
            print_modes(c.m, mr, r.modes_damped);
            write_modes("", c.m, mr);
        }
        if(r.vtp_dt >= 0.0)
            mkdir("coral_vtp", 0777);
        int nvtp = 0;
        auto vtp = [&]() {char b[64]; std::snprintf(b, sizeof(b), "coral_vtp/structure_%05d.vtp", nvtp++); write_vtp(c.m, b);};
        if(r.t_end <= r.t_start)
        {
            if(r.vtp_dt >= 0.0) vtp();
            return 0;
        }
        timeint ti(r.ti);
        ti.initialize(c.m, r.t_start);
        if(r.vtp_dt >= 0.0) vtp();
        const long nvtp_every = r.vtp_dt > 0.0 ? std::max(1L, std::lround(r.vtp_dt/r.dt)) : 1L;
        time_series out;
        out.open(r.out_file, c.m);
        out.write(c.m, r.t_start);
        const long n = std::lround((r.t_end - r.t_start)/r.dt);
        const long nout = r.out_dt > 0.0 ? std::max(1L, std::lround(r.out_dt/r.dt)) : 1L;
        long next_report = n/10;
        for(long i=0; i<n; ++i)
        {
            const double t = r.t_start + i*r.dt;
            if(!ti.step(c.m, t, r.dt))
            {
                // sub-steps of the failed step
                bool ok = false;
                for(int level=1; level<=6 && !ok; ++level)
                {
                    const int k = 1<<level;
                    ok = true;
                    for(int s=0; s<k && ok; ++s)
                        ok = ti.step(c.m, t + s*r.dt/k, r.dt/k);
                }
                if(!ok)
                {
                    std::fprintf(stderr, "CORAL: the time step at t = %g s failed, stopping\n", t);
                    out.close();
                    return 2;
                }
            }
            if((i+1)%nout==0)
                out.write(c.m, r.t_start + (i+1)*r.dt);
            if(r.vtp_dt >= 0.0 && (i+1)%nvtp_every==0)
                vtp();
            if(n >= 10 && i+1 >= next_report)
            {
                std::printf("  t = %10.4f s  (%3.0f%%)\n", r.t_start + (i+1)*r.dt, 100.0*(i+1)/n);
                std::fflush(stdout);
                next_report += n/10;
            }
        }
        out.close();
        const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - w0).count();
        std::printf("done: %ld steps, %ld Newton iterations, %ld Jacobians, %ld retries, %ld failures, max |g| %.1e, %.2f s\n",
                    ti.stats.steps, ti.stats.iterations, ti.stats.jacobians, ti.stats.retries, ti.stats.failures,
                    ti.stats.max_constraint, wall);
        std::printf("results: %s\n", r.out_file.c_str());
    }
    catch(const std::exception &e)
    {
        std::fprintf(stderr, "CORAL: %s\n", e.what());
        return 1;
    }
    return 0;
}
