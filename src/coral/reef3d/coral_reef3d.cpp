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

#include"coral/reef3d/coral_reef3d.h"
#include"coral/core/coral_timeint.h"
#include"coral/core/coral_version.h"
#include"coral/io/coral_input.h"
#include"coral/io/coral_output.h"
#include"lexer.h"
#include"ghostcell.h"
#include<mpi.h>
#include<cmath>
#include<fstream>
#include<iostream>
#include<sstream>
#include<sys/stat.h>

coral_reef3d::coral_reef3d(lexer *p, ghostcell *pgc)
{
    // rank 0 reads coral.txt, all ranks build the same model
    std::string text;
    int ok = 1;
    if(p->mpirank==0)
    {
        std::ifstream in("coral.txt");
        if(!in)
            ok = 0;
        else
        {
            std::stringstream ss;
            ss << in.rdbuf();
            text = ss.str();
        }
    }
    MPI_Bcast(&ok, 1, MPI_INT, 0, pgc->mpi_comm);
    if(!ok)
        fail(pgc, "Z 40: coral.txt not found in the case directory");
    int n = (int)text.size();
    MPI_Bcast(&n, 1, MPI_INT, 0, pgc->mpi_comm);
    text.resize(n);
    if(n > 0)
        MPI_Bcast(&text[0], n, MPI_CHAR, 0, pgc->mpi_comm);

    case_ = std::make_unique<coral::coral_case>();
    try
    {
        case_->read_string(text, ".");
    }
    catch(const std::exception &e)
    {
        fail(pgc, e.what());
    }
    ti_ = std::make_unique<coral::timeint>(case_->run.ti);
    if(p->mpirank==0)
    {
        std::cout<<coral::banner()<<": "<<case_->m.nbody()<<" bodies, "<<case_->m.ncon()<<" constraint equations, "
                 <<case_->m.load_list.size()<<" loads, max. step "<<case_->run.dt<<" s"<<std::endl;
    }
}

coral_reef3d::~coral_reef3d()
{
    if(out_)
        out_->close();
}

void coral_reef3d::fail(ghostcell *pgc, const std::string &msg)
{
    std::cerr<<"\n!!! REEF3D::CORAL: "<<msg<<" !!!\n"<<std::endl;
    MPI_Abort(pgc->mpi_comm, 1);
}

void coral_reef3d::ini(lexer *p, ghostcell*)
{
    t_ = p->simtime;
    ti_->initialize(case_->m, t_);
    if(p->mpirank==0)
    {
        mkdir("./REEF3D_CORAL", 0777);
        std::string f = case_->run.out_file;
        const size_t s = f.find_last_of('/');
        if(s != std::string::npos)
            f = f.substr(s+1);
        out_ = std::make_unique<coral::time_series>();
        out_->open("./REEF3D_CORAL/" + f, case_->m);
        out_->write(case_->m, t_);
    }
    next_out_ = t_ + case_->run.out_dt;
}

void coral_reef3d::step(lexer *p, ghostcell *pgc)
{
    const double t1 = p->simtime + p->dt;
    const double hmax = case_->run.dt > 0.0 ? case_->run.dt : p->dt;
    const int n = std::max(1, (int)std::ceil((t1 - t_)/hmax - 1.0e-9));
    const double h = (t1 - t_)/n;
    for(int k=0; k<n; ++k)
    {
        const double t = t_ + k*h;
        if(!ti_->step(case_->m, t, h))
        {
            // sub-steps of the failed step
            bool ok = false;
            for(int level=1; level<=6 && !ok; ++level)
            {
                const int m = 1<<level;
                ok = true;
                for(int s=0; s<m && ok; ++s)
                    ok = ti_->step(case_->m, t + s*h/m, h/m);
            }
            if(!ok)
                fail(pgc, "time step failed at t = " + std::to_string(t) + " s");
        }
        ++nstep_;
        if(out_ && (case_->run.out_dt <= 0.0 || t + h >= next_out_ - 1.0e-9*h))
        {
            out_->write(case_->m, t + h);
            if(case_->run.out_dt > 0.0)
                while(next_out_ <= t + h + 1.0e-9*h)
                    next_out_ += case_->run.out_dt;
        }
    }
    t_ = t1;
}
