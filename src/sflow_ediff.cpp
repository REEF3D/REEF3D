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

#include"sflow_ediff.h"
#include"lexer.h"
#include"fdm2D.h"

sflow_ediff::sflow_ediff(lexer* p)
{
}

sflow_ediff::~sflow_ediff()
{
}

void sflow_ediff::diff_u(lexer* p, fdm2D *b, ghostcell *pgc, solver2D *psolv, slice &UHdiff, slice &UH, slice &U, slice &V, slice &WL, double alpha)
{
	SLICELOOP4
    {
    UHdiff(i,j) = UH(i,j);
    
    WETDRY
	b->F(i,j) += WL(i,j)*viscosity(p,b)*(2.0*dxx(p,U) + dyy(p,U) + dxy(p,V));
    }
}

void sflow_ediff::diff_v(lexer* p, fdm2D *b, ghostcell *pgc, solver2D *psolv, slice &VHdiff, slice &VH, slice &U, slice &V, slice &WL, double alpha)
{
	SLICELOOP4
    {
    VHdiff(i,j) = VH(i,j);
    
    WETDRY
	b->G(i,j) += WL(i,j)*viscosity(p,b)*(dxx(p,V) + 2.0*dyy(p,V) + dxy(p,U))*p->y_dir;
    }
}

void sflow_ediff::diff_w(lexer* p, fdm2D *b, ghostcell *pgc, solver2D *psolv, slice &WHdiff, slice &WH, slice &U, slice &V, slice &W, slice &WL, double alpha)
{
	SLICELOOP4
    {
    WHdiff(i,j) = WH(i,j);
    
    WETDRY
	b->H(i,j) += WL(i,j)*viscosity(p,b)*(dxx(p,W) + dyy(p,W));
    }
}

void sflow_ediff::diff_scalar(lexer* p, fdm2D *b, ghostcell *pgc, solver2D *psolv, slice &f, double sig, double alpha)
{
}

double sflow_ediff::viscosity(lexer *p, fdm2D *b)
{
    double visc = p->W2 + b->eddyv(i,j);
    
    if(p->A246==2 && b->breaking(i,j)==1)
    visc += p->A250;
    
    return visc;
}

// dry neighbours (inside the domain, wet 0) are free-slip, zero normal gradient: the cell's own value is
// used instead of the dry cell's (U = 0 there would act as a no-slip wall at the shoreline); walls and
// open boundaries keep their ghost values
#define SFLOW_DRY(X) (p->flagslice4[X]>0 && p->wet[X]==0)

double sflow_ediff::dxx(lexer *p, slice &f)
{
    const double fp = SFLOW_DRY(Ip1J) ? f(i,j) : f(i+1,j);
    const double fm = SFLOW_DRY(Im1J) ? f(i,j) : f(i-1,j);
    
    return ((fp - f(i,j))/p->DXP[IP] - (f(i,j) - fm)/p->DXP[IM1])/p->DXN[IP];
}

double sflow_ediff::dyy(lexer *p, slice &f)
{
    if(p->j_dir==0)
    return 0.0;
    
    const double fp = SFLOW_DRY(IJp1) ? f(i,j) : f(i,j+1);
    const double fm = SFLOW_DRY(IJm1) ? f(i,j) : f(i,j-1);
    
    return ((fp - f(i,j))/p->DYP[JP] - (f(i,j) - fm)/p->DYP[JM1])/p->DYN[JP];
}

double sflow_ediff::dxy(lexer *p, slice &f)
{
    if(p->j_dir==0)
    return 0.0;
    
    if(p->flagslice4[Ip1Jp1]<0 || p->flagslice4[Ip1Jm1]<0 || p->flagslice4[Im1Jp1]<0 || p->flagslice4[Im1Jm1]<0)
    return 0.0;
    
    if(p->wet[Ip1Jp1]==0 || p->wet[Ip1Jm1]==0 || p->wet[Im1Jp1]==0 || p->wet[Im1Jm1]==0)   // no cross term next to dry cells
    return 0.0;
    
    return (f(i+1,j+1) - f(i+1,j-1) - f(i-1,j+1) + f(i-1,j-1))/((p->DXP[IP]+p->DXP[IM1])*(p->DYP[JP]+p->DYP[JM1]));
}
