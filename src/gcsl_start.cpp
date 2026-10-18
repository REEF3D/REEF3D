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
#include"slice.h"

// The distro routines write margin cells deep behind a boundary cell. Per entry the depth is
// limited to the solid depth of column [2] of the gcbsl lists (gcb_seed; 0: not set, full
// margin): a solid one or two cells thick inside the domain is filled from both sides, but
// the fluid cells behind it are not overwritten.
static inline int gcsl_depth(int d, int m0)
{
    return (d>0 && d<m0) ? d : m0;
}

void ghostcell::gcsl_start1(lexer *p, slice &f, int gcv)
{
    const int m0 = margin;
    starttime=timer();

    QQGCSL1LOOP
    {
    margin = gcsl_depth(p->gcbsl1[qq][2],m0);
    gcsldistro1(p,f,p->gcbsl1[qq][0], p->gcbsl1[qq][1], p->gcbsl1[qq][5], gcv, p->gcbsl1[qq][4], p->gcbsl1[qq][3]);
    margin = m0;
    }

    endtime=timer();
    p->gctime+=endtime-starttime;

    //  MPI Boundary Swap
    if(do_comms)
    {
        starttime=timer();
        gcslparax(p,f,1);
        gcslparacox(p,f,gcv);
        endtime=timer();
        p->xtime+=endtime-starttime;
    }

    if(p->Y40==1 || p->Y40==3)
    dgcslpol1(p,f);

}

void ghostcell::gcsl_start2(lexer *p, slice &f, int gcv)
{
    const int m0 = margin;
    starttime=timer();

    QQGCSL2LOOP
    {
    margin = gcsl_depth(p->gcbsl2[qq][2],m0);
    gcsldistro2(p,f,p->gcbsl2[qq][0], p->gcbsl2[qq][1], p->gcbsl2[qq][5], gcv, p->gcbsl2[qq][4], p->gcbsl2[qq][3]);
    margin = m0;
    }

    endtime=timer();
    p->gctime+=endtime-starttime;

    //  MPI Boundary Swap
    if(do_comms)
    {
        starttime=timer();
        gcslparax(p,f,2);
        gcslparacox(p,f,gcv);
        endtime=timer();
        p->xtime+=endtime-starttime;
    }

    if(p->Y40==1  || p->Y40==3)
    dgcslpol2(p,f);

}

void ghostcell::gcsl_start4(lexer *p, slice &f, int gcv)
{
    const int m0 = margin;
    starttime=timer();
    QQGCSL4LOOP
    {
    margin = gcsl_depth(p->gcbsl4[qq][2],m0);
    gcsldistro4(p,f,p->gcbsl4[qq][0],p->gcbsl4[qq][1], p->gcbsl4[qq][5], gcv, p->gcbsl4[qq][4], p->gcbsl4[qq][3]);
    margin = m0;
    }
    endtime=timer();
    p->gctime+=endtime-starttime;
    
    solid_forcing_bed(p,f);

    //  MPI Boundary Swap
    if(do_comms)
    {
        starttime=timer();
        gcslparax(p,f,4);
        gcslparacox(p,f,gcv);
        endtime=timer();
        p->xtime+=endtime-starttime;
    }

    if(p->Y40==1  || p->Y40==3)
    dgcslpol4(p,f);
}

void ghostcell::gcsl_start1int(lexer *p, sliceint &f, int gcv)
{
    const int m0 = margin;
    starttime=timer();
    QQGCSL1LOOP
    {
    margin = gcsl_depth(p->gcbsl1[qq][2],m0);
    gcsldistro1int(p,f,p->gcbsl1[qq][0],p->gcbsl1[qq][1], p->gcbsl1[qq][5], gcv, p->gcbsl1[qq][4], p->gcbsl1[qq][3]);
    margin = m0;
    }
    endtime=timer();
    p->gctime+=endtime-starttime;

    //  MPI Boundary Swap
    if(do_comms)
    {
        starttime=timer();
        gcslparax_int(p,f,1);
        gcslparacox_int(p,f,gcv);
        endtime=timer();
        p->xtime+=endtime-starttime;
    }
}

void ghostcell::gcsl_start2int(lexer *p, sliceint &f, int gcv)
{
    const int m0 = margin;
    starttime=timer();
    QQGCSL2LOOP
    {
    margin = gcsl_depth(p->gcbsl2[qq][2],m0);
    gcsldistro2int(p,f,p->gcbsl2[qq][0],p->gcbsl2[qq][1], p->gcbsl2[qq][5], gcv, p->gcbsl2[qq][4], p->gcbsl2[qq][3]);
    margin = m0;
    }
    endtime=timer();
    p->gctime+=endtime-starttime;

    //  MPI Boundary Swap
    if(do_comms)
    {
        starttime=timer();
        gcslparax_int(p,f,1);
        gcslparacox_int(p,f,gcv);
        endtime=timer();
        p->xtime+=endtime-starttime;
    }
}

void ghostcell::gcsl_start4int(lexer *p, sliceint &f, int gcv)
{
    const int m0 = margin;
    starttime=timer();
    QQGCSL4LOOP
    {
    margin = gcsl_depth(p->gcbsl4[qq][2],m0);
    gcsldistro4int(p,f,p->gcbsl4[qq][0],p->gcbsl4[qq][1], p->gcbsl4[qq][5], gcv, p->gcbsl4[qq][4], p->gcbsl4[qq][3]);
    margin = m0;
    }
    endtime=timer();
    p->gctime+=endtime-starttime;

    //  MPI Boundary Swap
    if(do_comms)
    {
        starttime=timer();
        gcslparax_int(p,f,4);
        gcslparacox_int(p,f,gcv);
        endtime=timer();
        p->xtime+=endtime-starttime;
    }
}

void ghostcell::gcsl_start4Vint(lexer *p, int *f, int gcv)
{
    const int m0 = margin;
    starttime=timer();
    QQGCSL4LOOP
    {
    margin = gcsl_depth(p->gcbsl4[qq][2],m0);
    gcsldistro4Vint(p,f,p->gcbsl4[qq][0],p->gcbsl4[qq][1], p->gcbsl4[qq][5], gcv, p->gcbsl4[qq][4], p->gcbsl4[qq][3]);
    margin = m0;
    }
    endtime=timer();
    p->gctime+=endtime-starttime;

    SLICELOOP4
    {
        if(p->flagslice4[Im1Jm1]<0 && p->flagslice4[Im1J]<0 && p->flagslice4[IJm1]<0)
        f[Im1Jm1] = f[IJ];

        if(p->flagslice4[Im1Jp1]<0 && p->flagslice4[Im1J]<0 && p->flagslice4[IJp1]<0)
        f[Im1Jp1] = f[IJ];

        if(p->flagslice4[Ip1Jm1]<0 && p->flagslice4[Ip1J]<0 && p->flagslice4[IJm1]<0)
        f[Ip1Jm1] = f[IJ];

        if(p->flagslice4[Ip1Jp1]<0 && p->flagslice4[Ip1J]<0 && p->flagslice4[IJp1]<0)
        f[Ip1Jp1] = f[IJ];
    }

    //  MPI Boundary Swap
    if(do_comms)
    {
        starttime=timer();
        gcslparaxV_int(p,f,4);
        gcslparacoxV_int(p,f,gcv);
        endtime=timer();
        p->xtime+=endtime-starttime;
    }
}
