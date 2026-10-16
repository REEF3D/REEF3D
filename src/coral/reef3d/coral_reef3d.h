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
// Phase 1: CORAL's own engineering hydrodynamics (Cummins, Morison, mooring, PTO); the loads from the
// resolved fluid and the body motion back to the fluid follow with the fluid interfaces.

#ifndef CORAL_REEF3D_H_
#define CORAL_REEF3D_H_

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

class coral_reef3d
{
public:
    coral_reef3d(lexer*, ghostcell*);
    ~coral_reef3d();

    // set up at p->simtime (consistent accelerations, first output)
    void ini(lexer*, ghostcell*);
    // advance from p->simtime to p->simtime + p->dt
    void step(lexer*, ghostcell*);

private:
    void fail(ghostcell*, const std::string&);

    std::unique_ptr<coral::coral_case> case_;
    std::unique_ptr<coral::timeint> ti_;
    std::unique_ptr<coral::time_series> out_;
    double t_ = 0.0;
    double next_out_ = 0.0;
    long nstep_ = 0;
};

#endif
