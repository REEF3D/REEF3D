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
#include"sliceint.h"

// The distro routines write margin cells deep behind a boundary cell. Per entry the depth is
// limited to the solid depth of column [2] of the gcbsl lists (gcb_seed; 0: not set, full
// margin): a solid one or two cells thick inside the domain is filled from both sides, but
// the fluid cells behind it are not overwritten.
static inline int gcsl_depth(int d, int m0)
{
    return (d>0 && d<m0) ? d : m0;
}

namespace
{
// the label of an entry (gcsleval*) depends on gcv, its bc and side cs and the settings of the run,
// not on the cell: evaluated once per (bc, cs) of a call
struct gcsl_labels
{
    int n = 0;
    int bc[8], cs[8], lab[8];
    template<class EV>
    int get(int b, int c, EV ev)
    {
        for(int k=0; k<n; ++k)
        if(bc[k]==b && cs[k]==c)
        return lab[k];
        const int l = ev();
        if(n<8)
        {
            bc[n]=b; cs[n]=c; lab[n]=l; ++n;
        }
        return l;
    }
};

// gcsl_neumann / gcsl_noslip / gcsl_neumann_int at the cell c (sx: the stride in i), mg cells deep,
// inline (the same writes in the same order)
template<class T>
inline void gcsl_neumann_at(T *c, int sx, int cs, int mg)
{
    if(cs==1)
    for(int q=0;q<mg;++q)
    c[-(q+1)*sx] = c[0];

    if(cs==2)
    for(int q=0;q<mg;++q)
    c[q+1] = c[0];

    if(cs==3)
    for(int q=0;q<mg;++q)
    c[-(q+1)] = c[0];

    if(cs==4)
    for(int q=0;q<mg;++q)
    c[(q+1)*sx] = c[0];
}

inline void gcsl_zero_at(double *c, int sx, int cs, int mg)
{
    if(cs==1)
    for(int q=0;q<mg;++q)
    c[-(q+1)*sx] = 0.0;

    if(cs==2)
    for(int q=0;q<mg;++q)
    c[q+1] = 0.0;

    if(cs==3)
    for(int q=0;q<mg;++q)
    c[-(q+1)] = 0.0;

    if(cs==4)
    for(int q=0;q<mg;++q)
    c[(q+1)*sx] = 0.0;
}
}

void ghostcell::gcsl_start1(lexer *p, slice &f, int gcv)
{
    const int m0 = margin;
    starttime=timer();

    // (labels 4 and 5 inline, the others through gcsldistro1)
    gcsl_labels L;
    QQGCSL1LOOP
    {
    const int ii = p->gcbsl1[qq][0], jj = p->gcbsl1[qq][1], bc = p->gcbsl1[qq][4], cs = p->gcbsl1[qq][3];
    const int mg = gcsl_depth(p->gcbsl1[qq][2],m0);
    const int lab = L.get(bc,cs,[&]{ return gcsleval1(p,gcv,bc,cs); });
    if(lab==4)
    gcsl_neumann_at(&f(ii,jj),p->jmax,cs,mg);
    else if(lab==5)
    gcsl_zero_at(&f(ii,jj),p->jmax,cs,mg);
    else if(lab!=-1)
    {
    margin = mg;
    gcsldistro1(p,f,ii,jj,p->gcbsl1[qq][5],gcv,bc,cs);
    margin = m0;
    }
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

    // (labels 4 and 5 inline, the others through gcsldistro2)
    gcsl_labels L;
    QQGCSL2LOOP
    {
    const int ii = p->gcbsl2[qq][0], jj = p->gcbsl2[qq][1], bc = p->gcbsl2[qq][4], cs = p->gcbsl2[qq][3];
    const int mg = gcsl_depth(p->gcbsl2[qq][2],m0);
    const int lab = L.get(bc,cs,[&]{ return gcsleval2(p,gcv,bc,cs); });
    if(lab==4)
    gcsl_neumann_at(&f(ii,jj),p->jmax,cs,mg);
    else if(lab==5)
    gcsl_zero_at(&f(ii,jj),p->jmax,cs,mg);
    else if(lab!=-1)
    {
    margin = mg;
    gcsldistro2(p,f,ii,jj,p->gcbsl2[qq][5],gcv,bc,cs);
    margin = m0;
    }
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
    // (labels 4, 14 and 5 inline, the others through gcsldistro4)
    gcsl_labels L;
    QQGCSL4LOOP
    {
    const int ii = p->gcbsl4[qq][0], jj = p->gcbsl4[qq][1], bc = p->gcbsl4[qq][4], cs = p->gcbsl4[qq][3];
    const int mg = gcsl_depth(p->gcbsl4[qq][2],m0);
    const int lab = L.get(bc,cs,[&]{ return gcsleval4(p,gcv,bc,cs); });
    if(lab==4 || (lab==14 && (cs==1 || cs==4)))
    gcsl_neumann_at(&f(ii,jj),p->jmax,cs,mg);
    else if(lab==5)
    gcsl_zero_at(&f(ii,jj),p->jmax,cs,mg);
    else if(lab==7 || lab==8)
    {
    margin = mg;
    gcsldistro4(p,f,ii,jj,p->gcbsl4[qq][5],gcv,bc,cs);
    margin = m0;
    }
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
    // (gcsl_neumann_int inline)
    QQGCSL1LOOP
    gcsl_neumann_at(&f(p->gcbsl1[qq][0],p->gcbsl1[qq][1]),p->jmax,p->gcbsl1[qq][3],gcsl_depth(p->gcbsl1[qq][2],m0));
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
    // (gcsl_neumann_int inline)
    QQGCSL2LOOP
    gcsl_neumann_at(&f(p->gcbsl2[qq][0],p->gcbsl2[qq][1]),p->jmax,p->gcbsl2[qq][3],gcsl_depth(p->gcbsl2[qq][2],m0));
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
    // (gcsl_neumann_int inline)
    QQGCSL4LOOP
    gcsl_neumann_at(&f(p->gcbsl4[qq][0],p->gcbsl4[qq][1]),p->jmax,p->gcbsl4[qq][3],gcsl_depth(p->gcbsl4[qq][2],m0));
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
    // (gcsl_neumann_V_int inline)
    QQGCSL4LOOP
    gcsl_neumann_at(f + (p->gcbsl4[qq][0]-p->imin)*p->jmax + (p->gcbsl4[qq][1]-p->jmin),p->jmax,p->gcbsl4[qq][3],gcsl_depth(p->gcbsl4[qq][2],m0));
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
