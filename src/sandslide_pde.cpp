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

#include"sandslide_pde.h"
#include"sediment_fdm.h"
#include"sediment_mixture.h"
#include"lexer.h"
#include"ghostcell.h"

sandslide_pde::sandslide_pde(lexer *p) : norm_vec(p), bedslope(p), fh(p)
{
    if(p->S50==1)
	gcval_topo=151;

	if(p->S50==2)
	gcval_topo=152;

	if(p->S50==3)
	gcval_topo=153;
	
	if(p->S50==4)
	gcval_topo=154;

	dxs=sqrt(2.0*p->DXM*p->DXM);
	fac1 = (1.0/6.0);
	fac2 = (1.0/12.0);
}

sandslide_pde::~sandslide_pde()
{
}

void sandslide_pde::start(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    SLICEBASELOOP
    s->slide_fh(i,j)=0.0;

    // converged when no neighbour pair is steeper than phi by more than tol (height)
    tol = 1.0e-4*MAX(p->S20,1.0e-6);

    // the pair uses the mean angle of repose of its two cells: neighbour values across ranks
    pgc->gcsl_start4(p,s->phi,1);

    // mainloop
    for(int qn=0; qn<p->S91; ++qn)
    {
        count=0;

        // fill
        SLICEBASELOOP
        fh(i,j)=0.0;

        pgc->gcsl_start4(p,fh,1);

        if(s->pmix!=nullptr)
        s->pmix->slide_zero(p,pgc);

        // slide loop: sediment cells in the S 77 window only
        SEDSLICELOOP
        if(p->pos_x()>p->S77_xs && p->pos_x()<p->S77_xe)
        {
            slide(p,pgc,s);
        }

        pgc->gcslparax_fh(p,fh,4);

        // fill back
        SEDSLICELOOP
        {
        s->slide_fh(i,j)+=fh(i,j);
        s->bedzh(i,j)+=fh(i,j);
        }

        pgc->gcsl_start4(p,s->bedzh,1);

        // multi-fraction bed: sorting of the slid material
        if(s->pmix!=nullptr)
        s->pmix->slide_finish(p,pgc,s);

        count=pgc->globalimax(count);

        p->slidecells=count;

        if(p->slidecells==0)
        break;

        if(p->mpirank==0)
        cout<<"sandslide_pde corrections: "<<p->slidecells<<endl;
    }
}

void sandslide_pde::slide(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    // finite-volume relaxation of the bed towards the angle of repose: between the cell and each
    // of its 8 neighbours (4 faces, 4 diagonals) the volume
    //      V = K*min(A_i,A_n)*sign(z_n - z_i)*max(|z_n - z_i| - d*tan(phi), 0)
    // is exchanged, d the distance of the cell centres, phi the mean of the two cells.
    // Only the excess over the angle of repose moves, so the iteration converges to phi in
    // every direction; the exchange is symmetric (conservative on non-uniform grids).
    // (was: diffusion of the full height difference, switched on by the central-difference
    // cell slope; it stopped at about 37 deg for phi = 35 deg, steeper on the faces and diagonals)
    // Pairs to cells without an erodible bed (structures, DFBED<0, solids, domain boundary) and
    // to cells outside the S 77 window are closed; a diagonal pair needs at least one of the two
    // cells beside it open.
    const double K = 0.1;
    const int ni[8] = {1,-1,0,0,1,1,-1,-1};
    const int nj[8] = {0,0,1,-1,1,-1,1,-1};
    double flux[8];
    const double Ai = p->DXN[IP]*p->DYN[JP];
    int over=0;

    auto open = [&](int di, int dj)
    {
        int ii=i+di;
        int jj=j+dj;

        return SLIDE_NB(di,dj) && p->flagslice4[(ii-p->imin)*p->jmax + jj-p->jmin]>0 && p->DFBED[(ii-p->imin)*p->jmax + jj-p->jmin]>0
            && p->XP[IP+di]>p->S77_xs && p->XP[IP+di]<p->S77_xe;
    };

    fh(i,j) = 0.0;

    for(int f=0;f<8;++f)
    {
    flux[f] = 0.0;

    int di=ni[f];
    int dj=nj[f];

        if(!open(di,dj))
        continue;

        if(di!=0 && dj!=0 && !open(di,0) && !open(0,dj))
        continue;

    double ddx = di<0?p->DXP[IM1]:(di>0?p->DXP[IP]:0.0);
    double ddy = dj<0?p->DYP[JM1]:(dj>0?p->DYP[JP]:0.0);
    double d = sqrt(ddx*ddx + ddy*ddy);
    double An = p->DXN[IP+di]*p->DYN[JP+dj];

    double dz = s->bedzh(i+di,j+dj)-s->bedzh(i,j);
    double e = fabs(dz) - d*tan(0.5*(s->phi(i,j)+s->phi(i+di,j+dj)));

        if(e>0.0)
        {
        flux[f] = K*(MIN(Ai,An)/Ai)*(dz>0.0?e:-e);
        fh(i,j) += flux[f];

        if(e>tol)
        over=1;
        }
    }

    count += over;

    // multi-fraction bed: upwind composition for each pair
    if(s->pmix!=nullptr)
    s->pmix->slide_pde(p,s,i,j,ni,nj,flux,8);
}
