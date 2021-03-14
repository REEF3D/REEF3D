/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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
--------------------------------------------------------------------*/

#include"VOF_PLIC.h"
#include"gradient.h"
#include"initialize.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"solver.h"
#include"ghostcell.h"
#include"freesurface_header.h"
#include"ioflow.h"
#include"fluid_update_vof.h"
#include"heat.h"
#include"hires.h"
#include"weno_hj.h"
#include"hric.h"

VOF_PLIC::VOF_PLIC
(
    lexer* p,
    fdm *a,
    ghostcell* pgc,
    heat *pheat
):gradient(p),norm_vec(p),alpha(p),nx(p),ny(p),nz(p),vof1(p),vof2(p),vof3(p)
{
    if(p->F50==1)
    gcval_frac=71;

    if(p->F50==2)
    gcval_frac=72;

    if(p->F50==3)
    gcval_frac=73;

    if(p->F50==4)
    gcval_frac=74;

    pupdate = new fluid_update_vof(p,a,pgc);

    reini_ = new reini_RK3(p,1);

    sSweep = -1;

    ininorVecLS(p);
}

VOF_PLIC::~VOF_PLIC()
{
}


void VOF_PLIC::start
(
    fdm* a,
    lexer* p,
    convection* pconvec,
    solver* psolv,
    ghostcell* pgc,
    ioflow* pflow,
    reini* preini,
    particlecorr* ppart,
    field &F
)
{
    pflow->fsfinflow(p,a,pgc);

    starttime=pgc->timer();

    int sweep = 0;
    if (sSweep < 2)
    {
        sSweep++;
        sweep = sSweep;
    }
    else
    {
        sSweep = 0;
    }

    double Q1, Q2;

    // x-sweep (0), y-sweep (1), z-sweep (2)
    for (int nSweep = 0; nSweep < 3; nSweep++)
    {
        LOOP
        {
            //- Calculate left and right fluxes Q1 and Q2
            calcFlux(a, p, Q1, Q2, sweep);


            //- PLIC loop
            vof1(i, j, k) = 0.0;
            vof2(i, j, k) = 0.0;
            vof3(i, j, k) = 0.0;

            if (a->vof(i, j, k) >= 0.999)
            {
                // Fluxes leave and enter cell in a straight manner
                vof1(i, j, k) = max(-Q1, 0.0);
                vof2(i, j, k) = 1.0 - max(Q1, 0.0) + min(Q2, 0.0);
                vof3(i, j, k) = max(Q2, 0.0);
            }
            else
            {
                // Reconstruct plane in cell
                reconstructPlane(a, p);

                // Advect interface using Lagrangian approach
                advectPlane(a, p, Q1, Q2, sweep);

                // Update volume fraction
                updateVolumeFraction(a, p, Q1, Q2, sweep);
            }
        }


        //- Distribute volume fractions
        pgc->start4(p,vof1,gcval_frac);
        pgc->start4(p,vof2,gcval_frac);
        pgc->start4(p,vof3,gcval_frac);


        //- Calculate updated vof from volume fractions and distribute
        updateVOF(a, p, sweep);
        pgc->start4(p,a->vof,gcval_frac);


        //- Change sweep
        if (sweep < 2)
        {
            sweep++;
        }
        else
        {
            sweep = 0;
        }
    }

    //- Redistance distance function from updated plane equations
    //redistance(a, p, pdisc, pgc, pflow, 20);
    //- Distribute ls function
    //pgc->start4(p,a->phi,gcval_frac);

    pflow->vof_relax(p,pgc,a->vof);
	pgc->start4(p,a->vof,gcval_frac);
    pupdate->start(p,a,pgc);

    p->lsmtime=pgc->timer()-starttime;

    if(p->mpirank==0)
    cout<<"vofplictime: "<<setprecision(3)<<p->lsmtime<<endl;

    LOOP
    {
        //a->test(i,j,k) = a->vof(i,j,k);

        if (a->vof(i,j,k) > 0.5)
        {
            a->phi(i,j,k) = 1.0;
        }
        else
        {
            a->phi(i,j,k) = -1.0;
        }
    }

    for (int tt = 0; tt < 10; tt++)
    {
        reini_->start(a,p,a->phi,pgc,pflow);
    }

}


void VOF_PLIC::ltimesave(lexer* p, fdm *a, field &F)
{
}


void VOF_PLIC::update
(
    lexer *p,
    fdm *a,
    ghostcell *pgc,
    field &F
)
{
    pupdate->start(p,a,pgc);
}
