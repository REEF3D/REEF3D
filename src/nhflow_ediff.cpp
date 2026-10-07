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

#include"nhflow_ediff.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"solver.h"

// neighbour value in the diffusion stencil, the same boundary rule as nhflow_idiff: faces to solid or ghost
// cells (flag4 < 0 or DF < 0) have zero gradient with A 513 1 and take the ghost value with A 513 2. Without it
// the bed ghost of A 518 2 (zero) acted as an additional no-slip wall on top of the bed shear (A 519)
static inline double ediff_nb(lexer *p, const double *f, int c, int n)
{
    if(p->A513==1 && (p->flag4[n]<0 || p->DF[n]<0))
    return f[c];
    
    return f[n];
}

nhflow_ediff::nhflow_ediff(lexer* p)
{
	gcval_u=10;
	gcval_v=11;
	gcval_w=12;
    
    gcval_uh=14;
	gcval_vh=15;
	gcval_wh=16;
}

nhflow_ediff::~nhflow_ediff()
{
}

// explicit momentum diffusion (A 512 1), the same operator as nhflow_idiff: backward differences with the
// backward spacing (was DXP[IP] on both sides), face metrics, the sigxx term, factor 2 only on the normal
// stress (in w only on sigz^2), and the v equation uses VH in its sigma cross terms (it used UH)
void nhflow_ediff::diff_u(lexer *p, fdm_nhf *d, ghostcell *pgc, ioflow *pflow, solver *psolv, double *UHdiff, double *UHin, double *UH, double *VH, double *WH, slice &WL, double alpha)
{
    LOOP
    UHdiff[IJK] = UHin[IJK];
    
    pgc->start4V(p,UHdiff,gcval_uh);
    
    pflow->rkinflow_nhflow(p,d,pgc,UHdiff,UHin);

    LOOP
    {
    if(p->wet[IJ]==0 || p->DF[IJK]<0)
    continue;

    visc = d->VISC[IJK] + d->EV[IJK];

    // vertical coefficient at the faces k-1/2 and k+1/2
    const double s2b = pow(p->sigx[FIJK],2.0) + pow(p->sigy[FIJK],2.0) + pow(p->sigz[IJ],2.0);
    const double s2t = pow(p->sigx[FIJKp1],2.0) + pow(p->sigy[FIJKp1],2.0) + pow(p->sigz[IJ],2.0);

    const double fc = UH[IJK];
    const double fim = ediff_nb(p,UH,IJK,Im1JK), fip = ediff_nb(p,UH,IJK,Ip1JK);
    const double fjm = ediff_nb(p,UH,IJK,IJm1K), fjp = ediff_nb(p,UH,IJK,IJp1K);
    const double fkm = ediff_nb(p,UH,IJK,IJKm1), fkp = ediff_nb(p,UH,IJK,IJKp1);

    d->F[IJK] += 2.0*visc*((fip-fc)/p->DXP[IP] - (fc-fim)/p->DXP[IM1])/p->DXN[IP]

                   + visc*((fjp-fc)/p->DYP[JP] - (fc-fjm)/p->DYP[JM1])/p->DYN[JP]*p->y_dir

                   + visc*(s2t*(fkp-fc)/p->DZP[KP] - s2b*(fc-fkm)/p->DZP[KM1])/p->DZN[KP]

                   + visc*p->sigxx[FIJK]*(fkp-fkm)/(p->DZP[KP]+p->DZP[KM1])
        // transpose stress (leading order, as nhflow_idiff)
         + visc*((VH[Ip1Jp1K]-VH[Im1Jp1K]) - (VH[Ip1Jm1K]-VH[Im1Jm1K]))/((p->DXP[IP]+p->DXP[IM1])*(p->DYP[JP]+p->DYP[JM1]))*p->y_dir
         + visc*((WH[Ip1JKp1]-WH[Im1JKp1]) - (WH[Ip1JKm1]-WH[Im1JKm1]))/((p->DXP[IP]+p->DXP[IM1])*(p->DZP[KP]+p->DZP[KM1]))*p->sigz[IJ]

        + visc*2.0*0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(UH[Ip1JKp1] - UH[Im1JKp1] - UH[Ip1JKm1] + UH[Im1JKm1])
                            /((p->DXP[IP]+p->DXP[IM1])*(p->DZP[KP]+p->DZP[KM1]))
                        
        + visc*2.0*0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(UH[IJp1Kp1] - UH[IJm1Kp1] - UH[IJp1Km1] + UH[IJm1Km1])
                            /((p->DYP[JP]+p->DYP[JM1])*(p->DZP[KP]+p->DZP[KM1]))*p->y_dir;
	}
}

void nhflow_ediff::diff_v(lexer *p, fdm_nhf *d, ghostcell *pgc, ioflow *pflow, solver *psolv, double *VHdiff, double *VHin, double *UH, double *VH, double *WH, slice &WL, double alpha)
{
    LOOP
    VHdiff[IJK] = VHin[IJK];
    
    pgc->start4V(p,VHdiff,gcval_vh);
    
    pflow->rkinflow_nhflow(p,d,pgc,VHdiff,VHin);

    LOOP
    {
    if(p->wet[IJ]==0 || p->DF[IJK]<0)
    continue;

    visc = d->VISC[IJK] + d->EV[IJK];

    // vertical coefficient at the faces k-1/2 and k+1/2
    const double s2b = pow(p->sigx[FIJK],2.0) + pow(p->sigy[FIJK],2.0) + pow(p->sigz[IJ],2.0);
    const double s2t = pow(p->sigx[FIJKp1],2.0) + pow(p->sigy[FIJKp1],2.0) + pow(p->sigz[IJ],2.0);

    const double fc = VH[IJK];
    const double fim = ediff_nb(p,VH,IJK,Im1JK), fip = ediff_nb(p,VH,IJK,Ip1JK);
    const double fjm = ediff_nb(p,VH,IJK,IJm1K), fjp = ediff_nb(p,VH,IJK,IJp1K);
    const double fkm = ediff_nb(p,VH,IJK,IJKm1), fkp = ediff_nb(p,VH,IJK,IJKp1);

    d->G[IJK] += visc*((fip-fc)/p->DXP[IP] - (fc-fim)/p->DXP[IM1])/p->DXN[IP]

                   + 2.0*visc*((fjp-fc)/p->DYP[JP] - (fc-fjm)/p->DYP[JM1])/p->DYN[JP]*p->y_dir

                   + visc*(s2t*(fkp-fc)/p->DZP[KP] - s2b*(fc-fkm)/p->DZP[KM1])/p->DZN[KP]

                   + visc*p->sigxx[FIJK]*(fkp-fkm)/(p->DZP[KP]+p->DZP[KM1])
        // transpose stress (leading order, as nhflow_idiff)
         + visc*((UH[Ip1Jp1K]-UH[Ip1Jm1K]) - (UH[Im1Jp1K]-UH[Im1Jm1K]))/((p->DYP[JP]+p->DYP[JM1])*(p->DXP[IP]+p->DXP[IM1]))
         + visc*((WH[IJp1Kp1]-WH[IJm1Kp1]) - (WH[IJp1Km1]-WH[IJm1Km1]))/((p->DYP[JP]+p->DYP[JM1])*(p->DZP[KP]+p->DZP[KM1]))*p->sigz[IJ]

        + visc*2.0*0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(VH[Ip1JKp1] - VH[Im1JKp1] - VH[Ip1JKm1] + VH[Im1JKm1])
                            /((p->DXP[IP]+p->DXP[IM1])*(p->DZP[KP]+p->DZP[KM1]))
                        
        + visc*2.0*0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(VH[IJp1Kp1] - VH[IJm1Kp1] - VH[IJp1Km1] + VH[IJm1Km1])
                            /((p->DYP[JP]+p->DYP[JM1])*(p->DZP[KP]+p->DZP[KM1]))*p->y_dir;
	}
}

void nhflow_ediff::diff_w(lexer *p, fdm_nhf *d, ghostcell *pgc, ioflow *pflow, solver *psolv, double *WHdiff, double *WHin, double *UH, double *VH, double *WH, slice &WL, double alpha)
{
    LOOP
    WHdiff[IJK] = WHin[IJK];
    
    pgc->start4V(p,WHdiff,gcval_wh);
    
    pflow->rkinflow_nhflow(p,d,pgc,WHdiff,WHin);

    LOOP
    {
    if(p->wet[IJ]==0 || p->DF[IJK]<0)
    continue;

    visc = d->VISC[IJK] + d->EV[IJK];

    // vertical coefficient at the faces k-1/2 and k+1/2; w: + sigz^2 once more from 2 dw/dz
    const double s2b = pow(p->sigx[FIJK],2.0) + pow(p->sigy[FIJK],2.0) + 2.0*pow(p->sigz[IJ],2.0);
    const double s2t = pow(p->sigx[FIJKp1],2.0) + pow(p->sigy[FIJKp1],2.0) + 2.0*pow(p->sigz[IJ],2.0);

    const double fc = WH[IJK];
    const double fim = ediff_nb(p,WH,IJK,Im1JK), fip = ediff_nb(p,WH,IJK,Ip1JK);
    const double fjm = ediff_nb(p,WH,IJK,IJm1K), fjp = ediff_nb(p,WH,IJK,IJp1K);
    const double fkm = ediff_nb(p,WH,IJK,IJKm1), fkp = ediff_nb(p,WH,IJK,IJKp1);

    d->H[IJK] += visc*((fip-fc)/p->DXP[IP] - (fc-fim)/p->DXP[IM1])/p->DXN[IP]

                   + visc*((fjp-fc)/p->DYP[JP] - (fc-fjm)/p->DYP[JM1])/p->DYN[JP]*p->y_dir

                   + visc*(s2t*(fkp-fc)/p->DZP[KP] - s2b*(fc-fkm)/p->DZP[KM1])/p->DZN[KP]

                   + visc*p->sigxx[FIJK]*(fkp-fkm)/(p->DZP[KP]+p->DZP[KM1])
        // transpose stress (leading order, as nhflow_idiff)
         + visc*((UH[Ip1JKp1]-UH[Ip1JKm1]) - (UH[Im1JKp1]-UH[Im1JKm1]))/((p->DZP[KP]+p->DZP[KM1])*(p->DXP[IP]+p->DXP[IM1]))*p->sigz[IJ]
         + visc*((VH[IJp1Kp1]-VH[IJp1Km1]) - (VH[IJm1Kp1]-VH[IJm1Km1]))/((p->DYP[JP]+p->DYP[JM1])*(p->DZP[KP]+p->DZP[KM1]))*p->sigz[IJ]*p->y_dir

        + visc*2.0*0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(WH[Ip1JKp1] - WH[Im1JKp1] - WH[Ip1JKm1] + WH[Im1JKm1])
                            /((p->DXP[IP]+p->DXP[IM1])*(p->DZP[KP]+p->DZP[KM1]))
                        
        + visc*2.0*0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(WH[IJp1Kp1] - WH[IJm1Kp1] - WH[IJp1Km1] + WH[IJm1Km1])
                            /((p->DYP[JP]+p->DYP[JM1])*(p->DZP[KP]+p->DZP[KM1]))*p->y_dir;
	}
}

void nhflow_ediff::diff_scalar(lexer *p, fdm_nhf *d, ghostcell *pgc, solver *psolv, double *F, double sig, double alpha)
{
}