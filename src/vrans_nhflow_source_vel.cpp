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

#include"vrans_nhflow_f.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"slice.h"

#define VRANS_NHFLOW_DRAG_EXPLICIT (p->B204==0)

// Darcy-Forchheimer resistance of the static porous structures (B 201, B 202 layers, B 210-230).
// The coefficients use the static-structure porosity n_s = 1 - H(1 - n), not d->POR, so that
// porous floating bodies (X 16), which carry their own resistance, are not counted twice.
// (1 - n_s) -> 0 outside the structure, so no additional Heaviside weighting is needed.
void vrans_nhflow_f::porous_coeff(lexer *p, fdm_nhf *d)
{
    H = Hporface(p,d,0,0,0);
    
    // porous layers (B 202): cell values, B 201 otherwise
    if(p->B202>0)
    {
    porval   = LN[IJK];
    partval  = LD50[IJK];
    alphaval = LALPHA[IJK];
    betaval  = LBETA[IJK];
    }
    
    else
    {
    porval   = p->B201_n;
    partval  = p->B201_d50;
    alphaval = p->B201_alpha;
    betaval  = p->B201_beta;
    }
    
    porval  = 1.0 - H*(1.0 - porval);
    Aporval = Apor(porval,partval,alphaval,d->VISC[IJK]);
    Bporval = Bpor(porval,partval,betaval);
}

void vrans_nhflow_f::u_source(lexer *p, fdm_nhf *d, slice &WL)
{
    // B 204 1: the drag is applied point-implicitly after the stage update (implicit_drag)
    if(p->B200>0 && VRANS_NHFLOW_DRAG_EXPLICIT)
    LOOP
	{
        porous_coeff(p,d);
        
        if(H<1.0e-12)
        continue;

        porousterm = Aporval*d->U[IJK] + Bporval*d->U[IJK]*fabs(d->U[IJK]); 
    	
        d->Fext[IJK] -= WL(i,j)*porousterm;
    
        d->maxF = MAX(fabs(d->maxF),fabs(WL(i,j)*porousterm));
	}
}

void vrans_nhflow_f::v_source(lexer *p, fdm_nhf *d, slice &WL)
{
    // B 204 1: the drag is applied point-implicitly after the stage update (implicit_drag)
    if(p->B200>0 && VRANS_NHFLOW_DRAG_EXPLICIT)
    LOOP
	{
        porous_coeff(p,d);
        
        if(H<1.0e-12)
        continue;

        porousterm = Aporval*d->V[IJK] + Bporval*d->V[IJK]*fabs(d->V[IJK]); 
    	
        d->Gext[IJK] -= WL(i,j)*porousterm;
    
        d->maxG = MAX(fabs(d->maxG),fabs(WL(i,j)*porousterm));
	}
}

void vrans_nhflow_f::w_source(lexer *p, fdm_nhf *d, slice &WL)
{
    // B 204 1: the drag is applied point-implicitly after the stage update (implicit_drag)
    if(p->B200>0 && VRANS_NHFLOW_DRAG_EXPLICIT)
    LOOP
	{
        porous_coeff(p,d);
        
        if(H<1.0e-12)
        continue;

        porousterm = Aporval*d->W[IJK] + Bporval*d->W[IJK]*fabs(d->W[IJK]); 
    	
        d->Hext[IJK] -= WL(i,j)*porousterm;
    
        d->maxH = MAX(fabs(d->maxH),fabs(WL(i,j)*porousterm));
	}
}

// Point-implicit Darcy-Forchheimer drag (B 204 1, default).
// The explicit drag -(A + B|u|) u is unstable for dt CPORNH (A + 2B|u|) > 2 (RK2). This is reached in a
// fine-grained core (B 201) below a coarse, permeable armour (B 202): at the run-up / infiltration front
// the projection pushes the armour flow into the core cells, where the core resistance acts on it, the
// drag overshoots and the non-hydrostatic pressure grows to MPa within a few steps.
// Here the drag is linearised in the stage input velocity U (Picard) and solved with the stage update,
//   UH_new (1 + fac dt CPORNH (A + B|U|)) = UH_stage,   fac = RK stage weight of the explicit terms,
// which is unconditionally stable and keeps the Darcy-Forchheimer balance in steady flow.
void vrans_nhflow_f::implicit_drag(lexer *p, fdm_nhf *d, double *UH, double fac, int comp)
{
    if(p->B200==0 || VRANS_NHFLOW_DRAG_EXPLICIT)
    return;
    
    double *VEL = comp==0 ? d->U : (comp==1 ? d->V : d->W);
    
    LOOP
	{
        porous_coeff(p,d);
        
        if(H<1.0e-12)
        continue;
        
        UH[IJK] /= 1.0 + fac*p->dt*CPORNH*(Aporval + Bporval*fabs(VEL[IJK]));
	}
}

double vrans_nhflow_f::Apor(double por, double part, double alpha, double visc)
{
	val = alpha*(pow(1.0-por,2.0)/pow(por,3.0))*(visc/pow(part,2.0));
    
    if(val!=val)
    val=0.0;
	
	return val;
}

double vrans_nhflow_f::Bpor(double por, double part, double beta)
{
	val = beta*(1.0 + 7.5/Cval)*((1.0-por)/pow(por,4.0))*(1.0/part);
    
    if(val!=val)
    val=0.0;
	
	return val;
}
