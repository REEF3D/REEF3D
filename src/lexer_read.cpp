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

#include"lexer.h"
#include"ghostcell.h"

void lexer::lexer_read(ghostcell *pgc)
{
    ii_recv=ii_send=0;
    dd_recv=dd_send=0;

    if(mpirank==0)
    control::read_control(this);

    // SFLOW: A 212 not given in ctrl.txt -> explicit horizontal diffusion when a turbulence model is on (A 260 > 0),
    // so that the eddy viscosity reaches the momentum equations; no diffusion otherwise (the previous default)
    if(mpirank==0 && A212<0)
    A212 = (A260>0) ? 1 : 0;

#ifndef REEF3D_USE_HYPRE
    // built without hypre (the default; hypre is opt-in via make HYPRE=1): the hypre solvers N 10 10-39 are replaced by REEFMG (N 10 1)
    if(mpirank==0 && N10>=10)
    {
        std::cout<<"N 10 "<<N10<<": this build has no hypre (hypre is opt-in: make HYPRE=1), using REEFMG (N 10 1)"<<std::endl;
        N10=1;
    }
#endif

    Iarray(ictrl,ctrlsize);
    Darray(dctrl,ctrlsize);

    if(mpirank==0)
    control::ctrlsend();

    pgc->globalctrl(this);

    if(mpirank>0)
    control::ctrlrecv();

    del_Iarray(ictrl,ctrlsize);
    del_Darray(dctrl,ctrlsize);

    read_grid();

    control::parse(this);

    lexer_ini();
}
