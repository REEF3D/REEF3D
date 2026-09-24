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

#include "fnpf_fsfbc_wd.h"
#include "lexer.h"
#include "fdm_fnpf.h"
#include "ghostcell.h"
#include "fnpf_coastline.h"

void fnpf_fsfbc_wd::wetdry(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &Fifsf)
{
    if(p->count==0)
    {
        SLICELOOP4
        {
            wetcoast(i,j)=1;
            p->wet[IJ]=1;
        }

        SLICELOOP4
        if(p->wd - c->bed(i,j) < c->wd_criterion)
        wetcoast(i,j)=0;

        if(p->A343>=1)
        SLICELOOP4
        if(p->wd - c->bed(i,j) < c->wd_criterion)
        p->wet[IJ]=0;

        SLICELOOP4
        c->WL(i,j) = MAX(c->wd_criterion, eta(i,j) + p->wd - c->bed(i,j));

        // initially wet cells may dry right away
        SLICELOOP4
        wetage(i,j) = p->A330;

        if(p->A343==2)
        {
            // runup: every cell below the maximum runup elevation A335 (above
            // still water) may become wet, not only the initially wet ones
            SLICELOOP4
            wetcoast(i,j) = (p->flagslice4[IJ]>0 && c->bed(i,j) - p->wd <= p->A335) ? 1 : 0;

            // dry cells carry a thin film: eta = bed + A344. This is the level
            // the rewetting test compares with, so it must be set before the
            // first time step (initial land eta is 0)
            SLICELOOP4
            if(p->wet[IJ]==0 && wetcoast(i,j)==1)
            eta(i,j) = c->wd_criterion - c->depth(i,j);

            pgc->gcsl_start4(p,eta,gcval_eta);
        }

        if(p->A343>=2)
        {
            pgc->gcsl_start4Vint(p,p->wet,50);
            wd_front_mask(p,pgc);
        }
    }
    else if(p->count>=1)
    {
        SLICELOOP4
        c->WL(i,j) = eta(i,j) + p->wd - c->bed(i,j);

        pgc->gcsl_start4(p,c->WL,50);

        // dynamic wetting-drying: 2 runup and rundown, 3 rundown only
        if(p->A343>=2)
        wetdry_dynamic(p,c,pgc,eta,Fifsf);

        if(p->A343==1)
        SLICELOOP4
        if(c->WL(i,j)<=c->wd_criterion && wetcoast(i,j)==1)
        {
            eta(i,j) = 1.1*c->wd_criterion - c->depth(i,j);
            c->WL(i,j) = eta(i,j) + c->depth(i,j);

            if(p->j_dir)
                Fifsf(i,j) = 0.25*(Fifsf(i-1,j) + Fifsf(i+1,j) + Fifsf(i,j-1) + Fifsf(i,j+1));
            else
                Fifsf(i,j) = 0.5*(Fifsf(i-1,j) + Fifsf(i+1,j));
        }

        pgc->gcsl_start4Vint(p,p->wet,50);
        pgc->gcsl_start4(p,eta,gcval_eta);
        pgc->gcsl_start4(p,c->WL,gcval_eta);
    }

    pgc->gcsl_start4Vint(p,p->wet,50);

    if(coastline_count==0)
    {
        pcoast->start(p,c,pgc,c->coastline,p->wet,c->wet_n);
        ++coastline_count;
    }

}


// ---------------------------------------------------------------------------
// Dynamic wetting-drying
//   A343 2: runup and rundown. All cells below the elevation A335 may wet and
//           dry, no coastline damping (only the front viscosity band A332)
//   A343 3: rundown only. Only initially wet cells may dry and rewet, the
//           coastline damping A341/A342 around the initial coastline is kept
//
// Compared with the former A343 2 implementation:
//  - the wet/dry flags change only once per time step (first call with a new
//    p->count, i.e. the first RK stage); later stages only enforce the
//    constraints, so the RK stages see one consistent mask
//  - drying keeps Fifsf; a rewetted cell gets Fifsf from its wet face
//    neighbours instead of the stale value/zero, which otherwise produces a
//    jump of the size of the Bernoulli drift of Fi (O(100) m^2/s) and thus
//    O(10) m/s surface velocities in the first step after rewetting
//  - hysteresis: rewetting needs a wet neighbour whose surface is at least
//    A331*A344 above the dry cell's level, and a rewetted cell stays wet for
//    at least A330 time steps (below the criterion it is clamped as for A343 1)
//  - dry cells are held at WL = A344 in every stage
//  - optional (A334 1): the volume added/removed by the clamps is taken
//    from/given to the wet face neighbours, so wetting-drying conserves volume
// ---------------------------------------------------------------------------
double fnpf_fsfbc_wd::wet_nb_average(lexer *p, slice &f)
{
    double sum=0.0;
    int n=0;

    if(p->wet[Im1J]==1 && p->flagslice4[Im1J]>0)
    {sum+=f(i-1,j); ++n;}

    if(p->wet[Ip1J]==1 && p->flagslice4[Ip1J]>0)
    {sum+=f(i+1,j); ++n;}

    if(p->j_dir==1)
    {
        if(p->wet[IJm1]==1 && p->flagslice4[IJm1]>0)
        {sum+=f(i,j-1); ++n;}

        if(p->wet[IJp1]==1 && p->flagslice4[IJp1]>0)
        {sum+=f(i,j+1); ++n;}
    }

    return n>0 ? sum/double(n) : f(i,j);
}

void fnpf_fsfbc_wd::wetdry_dynamic(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &Fifsf)
{
    const double crit   = c->wd_criterion;
    const double margin = p->A331*crit;

    // neighbour values are read below: make the halos current for this stage
    pgc->gcsl_start4(p,eta,gcval_eta);
    pgc->gcsl_start4(p,Fifsf,gcval_fifsf);

    // ---------------------------------------------------------------
    // 1. flag update, once per time step
    // ---------------------------------------------------------------
    if(p->count!=wd_flagcount)
    {
        wd_flagcount = p->count;

        int nrewet=0, ndry=0;

        SLICELOOP4
        {
            p->wet_n[IJ] = p->wet[IJ];
            temp[IJ] = p->wet[IJ];

            if(p->wet[IJ]==1 && wetage(i,j)<p->A330)
            ++wetage(i,j);
        }

        SLICELOOP4
        if(wetcoast(i,j)==1)
        {
            if(p->wet[IJ]==0)
            {
                int rewet=0;

                if(p->wet[Ip1J]==1 && p->flagslice4[Ip1J]>0 && eta(i+1,j)>eta(i,j)+margin && c->WL(i+1,j)>crit+margin)
                rewet=1;

                if(p->wet[Im1J]==1 && p->flagslice4[Im1J]>0 && eta(i-1,j)>eta(i,j)+margin && c->WL(i-1,j)>crit+margin)
                rewet=1;

                if(p->j_dir==1)
                {
                    if(p->wet[IJp1]==1 && p->flagslice4[IJp1]>0 && eta(i,j+1)>eta(i,j)+margin && c->WL(i,j+1)>crit+margin)
                    rewet=1;

                    if(p->wet[IJm1]==1 && p->flagslice4[IJm1]>0 && eta(i,j-1)>eta(i,j)+margin && c->WL(i,j-1)>crit+margin)
                    rewet=1;
                }

                if(rewet==1)
                {
                    temp[IJ] = 1;
                    wetage(i,j) = 0;
                    Fifsf(i,j) = wet_nb_average(p,Fifsf);   // uses the old mask: wet neighbours only
                    ++nrewet;
                }
            }
            else if(c->WL(i,j)<=crit && wetage(i,j)>=p->A330)
            {
                temp[IJ] = 0;
                ++ndry;
            }
        }

        SLICELOOP4
        if(wetcoast(i,j)==1)
        p->wet[IJ] = temp[IJ];

        pgc->gcsl_start4Vint(p,p->wet,50);

        nrewet = pgc->globalisum(nrewet);
        ndry   = pgc->globalisum(ndry);

        if(p->mpirank==0 && (p->count%p->P12==0))
        cout<<"wetdry: rewetted "<<nrewet<<"  dried "<<ndry<<endl;

        wd_front_mask(p,pgc);
    }

    // ---------------------------------------------------------------
    // 2. constraints, every stage
    // ---------------------------------------------------------------
    SLICELOOP4
    wd_dvol(i,j) = 0.0;

    SLICELOOP4
    if(wetcoast(i,j)==1)
    {
        if(p->wet[IJ]==0)
        {
            // dry: water level held at the criterion, Fifsf untouched.
            // Only cells that were wet at the start of this step exchange
            // volume; cells that stay dry have K=0 and any change there is
            // not water (relaxation zone, initial state)
            const double etad = crit - c->depth(i,j);

            if(p->wet_n[IJ]==1)
            wd_dvol(i,j) = etad - eta(i,j);

            eta(i,j) = etad;
        }
        else if(eta(i,j) + c->depth(i,j) < crit)
        {
            // wet but below the criterion (young cell or later RK stage): clamp as for A343 1
            const double etaw = 1.1*crit - c->depth(i,j);
            wd_dvol(i,j) = etaw - eta(i,j);
            eta(i,j) = etaw;
            Fifsf(i,j) = wet_nb_average(p,Fifsf);
        }
    }

    // ---------------------------------------------------------------
    // 3. volume redistribution to the wet face neighbours
    //    (gather formulation, MPI-safe through the halo exchange)
    // ---------------------------------------------------------------
    if(p->A334==1)
    {
        SLICELOOP4
        {
            int nw=0;

            if(wd_dvol(i,j)!=0.0)
            {
                if(p->wet[Im1J]==1 && p->flagslice4[Im1J]>0) ++nw;
                if(p->wet[Ip1J]==1 && p->flagslice4[Ip1J]>0) ++nw;

                if(p->j_dir==1)
                {
                    if(p->wet[IJm1]==1 && p->flagslice4[IJm1]>0) ++nw;
                    if(p->wet[IJp1]==1 && p->flagslice4[IJp1]>0) ++nw;
                }
            }

            wd_nwet(i,j) = double(nw);
        }

        pgc->gcsl_start4(p,wd_dvol,50);
        pgc->gcsl_start4(p,wd_nwet,50);

        SLICELOOP4
        if(p->wet[IJ]==1)
        {
            double take=0.0;

            if(p->flagslice4[Im1J]>0 && wd_nwet(i-1,j)>0.5)
            take += wd_dvol(i-1,j)/wd_nwet(i-1,j);

            if(p->flagslice4[Ip1J]>0 && wd_nwet(i+1,j)>0.5)
            take += wd_dvol(i+1,j)/wd_nwet(i+1,j);

            if(p->j_dir==1)
            {
                if(p->flagslice4[IJm1]>0 && wd_nwet(i,j-1)>0.5)
                take += wd_dvol(i,j-1)/wd_nwet(i,j-1);

                if(p->flagslice4[IJp1]>0 && wd_nwet(i,j+1)>0.5)
                take += wd_dvol(i,j+1)/wd_nwet(i,j+1);
            }

            // a donor is never pushed below the criterion
            if(take>0.0)
            take = MIN(take, MAX(0.0, eta(i,j) + c->depth(i,j) - crit));

            eta(i,j) -= take;
        }
    }

    SLICELOOP4
    c->WL(i,j) = eta(i,j) + c->depth(i,j);
}

// wet cells with a dry cell inside the +-3 cell WENO stencil (x and y lines);
// used for the wet-only eta gradients at the front (A336)
void fnpf_fsfbc_wd::wd_front_mask(lexer *p, ghostcell *pgc)
{
    SLICELOOP4
    {
        int front=0;

        if(p->wet[IJ]==1)
        for(int q=-3; q<=3; ++q)
        {
            if(p->wet[(i-p->imin+q)*p->jmax + (j-p->jmin)]==0)
            front=1;

            if(p->j_dir==1 && p->wet[(i-p->imin)*p->jmax + (j-p->jmin+q)]==0)
            front=1;
        }

        wdfront(i,j) = front;
    }
}
