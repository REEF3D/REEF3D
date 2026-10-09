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

#include"CPM.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"
#include<cstdlib>

/*--------------------------------------------------------------------
hybrid suspension (Q 58 3 and 4, S 10 1): the bed and the bedload layer are parcels, the suspended
load is the Eulerian concentration of REEF3D (suspended_IM1: convection with w - w_s, eddy
diffusion, implicit exchange with the bed in the first fluid cell above the bed)

The resolved suspension of Q 58 2 needs many small parcels: at moderate Shields numbers a column
holds a fraction of a suspended parcel, the suspended sand falls back within a few columns and
the profile is noisy. The concentration carries the suspension without this sampling noise, and
the parcels keep the bed: pickup into the bedload layer, hops, slopes and avalanches, the packed
bed, the seepage.

Exchange per bed column, in the first fluid cell above the fluid's bed (as the Eulerian model):
    erosion     E = w_s c_be,   c_be: van Rijn (1984) c_a at a = max(k_s, 2 d50, 0.01 h) <= 0.5 h,
                                transferred to the cell centre with the Rouse profile (bedconc_VR);
                                bed shear stress of the bedload layer (log law, Q 68), onset of
                                suspension after van Rijn as Q 58 2 (Q 58 4: none, pickup for T > 0
                                as the Eulerian model); zero where the column has no erodible parcels
    deposition  D = w_s c_1,    c_1: concentration of that cell (implicit in the solver)
The parcels take the net flux: per column the account blCs gains Q 65 (E - D) A dt (solid volume;
MORFAC as for the bedload layer). Above one parcel volume a parcel of the layer (moving, else
exposed) leaves the bed into the concentration; below minus one parcel volume a parcel is
deposited on the bed (placed at the bed level of the layer, at rest). The fluid's bed follows the
parcels: the concentration of a cell taken into the bed goes to the account of its column, a cell
released from the bed takes its concentration from it. The faces of the flow cells to the bed, solid
bodies, air and walls are closed for the concentration (suspended_IM1::bcsusp_start), so sediment
volume of parcels + concentration + accounts is conserved (log: suspended volume; CPM_SUSPDBG=1 prints
the balance of the exchange and the near-bed concentrations).
--------------------------------------------------------------------*/

// first fluid cell above the fluid's bed in the column (i,j), -1 if none: the cell of the
// bed exchange in suspended_IM1::suspsource
int CPM::susp_cell(lexer *p, fdm *a, int ii, int jj)
{
    if(p->DFBED[(ii-p->imin)*p->jmax + jj-p->jmin]<=0 || p->XP[ii+marge]<p->S71 || p->XP[ii+marge]>p->S72)
    return -1;
    
    for(int kk=1; kk<p->knoz; ++kk)
    {
        int ijk = (ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin;
        
        if(p->flag4[ijk]>0 && p->DF[ijk]>0 && a->topo(ii,jj,kk)>0.0 && a->topo(ii,jj,kk-1)<0.0)
        return a->phi(ii,jj,kk)>=0.0 ? kk : -1;
    }
    
    return -1;
}

// settling velocity and near-bed concentration of the erosion flux, before the concentration step
void CPM::susp_cbe(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s)
{
    const double gmag = sqrt(p->W20*p->W20 + p->W21*p->W21 + p->W22*p->W22);
    const double R = (p->S22 - p->W1)/p->W1;
    const double Dst = p->S20*pow(R*gmag/(p->W2*p->W2), 1.0/3.0);
    const double ws = settling_velocity(p,p->S20);
    const double ucs = (Dst<=10.0 ? 4.0/MAX(Dst,1.0) : 0.4)*ws;
    
    s->ws = ws;
    afdm = a;
    
    if(getenv("CPM_SUSPDBG")!=nullptr)
    dbg_v0 = susp_volume(p,a,pgc);
    
    // cells that changed between the bed and the flow since the last concentration step (the
    // fluid's bed follows the parcels): the concentration of a cell taken into the bed is deposited,
    // a cell released from the bed takes its concentration from the bed (column mass difference)
    double dsw=0.0;
    
    if(blMc_ok==1)
    for(i=0;i<p->knox;++i)
    for(j=0;j<p->knoy;++j)
    {
    double dm = susp_column(p,a,i,j) - blMc(i,j);
    blCs(i,j) += p->Q65*dm;
    dsw += dm;
    }
    
    if(getenv("CPM_SUSPDBG")!=nullptr)
    {
        dsw = pgc->globalsum(dsw);
        dbg_sw += dsw;
    }
    
    SLICELOOP4
    s->cbe(i,j) = 0.0;
    
    for(i=0;i<p->knox;++i)
    for(j=0;j<p->knoy;++j)
    {
        int k1 = susp_cell(p,a,i,j);
        
        // columns without erodible parcels (no bed, fixed floor, structure) do not erode
        if(k1<0 || blNc(i,j)<=0.0)
        continue;
        
        double tm = sqrt(blTx(i,j)*blTx(i,j) + blTy(i,j)*blTy(i,j));
        double ustar = sqrt(tm/p->W1);
        double T = tm/((p->S22-p->W1)*gmag*p->S20*p->Q60) - 1.0;
        
        // onset of suspension after van Rijn (1984) as Q 58 2; Q 58 4: none (pickup for T > 0, as the Eulerian model)
        if((p->Q58==3 && ustar<=ucs) || T<=0.0)
        continue;
        
        // water depth above the bed: free surface of the column (top cell in water + phi)
        double zb = p->ZP[k1+marge] - a->topo(i,j,k1);
        double zs = -1.0e20;
        
        for(k=k1; k<p->knoz; ++k)
        if(a->phi(i,j,k)>=0.0)
        zs = p->ZP[KP] + a->phi(i,j,k);
        
        zs = MIN(zs, p->global_zmax);   // lid (free surface above the domain)
        double h = zs>zb ? zs-zb : 1.0e20;
        
        double adist = MAX(p->S21*p->S20, 2.0*p->S20);
        adist = MAX(adist, h<1.0e19 ? 0.01*h : 0.0);
        adist = MIN(adist, 0.5*h);
        
        double ca = MIN(0.05, 0.015*p->S20*pow(T,1.5)/(adist*pow(Dst,0.3)));
        
        // Rouse transfer from z = a to the first cell centre
        double z1 = MAX(a->topo(i,j,k1), adist);
        z1 = MIN(z1, 0.99*h);
        double P = MIN(ws/(0.4*MAX(ustar,1.0e-6)), 5.0);
        double ratio = h<1.0e19 ? (adist/(h-adist))*((h-z1)/z1) : adist/z1;
        
        s->cbe(i,j) = ca*pow(ratio, P);
    }
    
    pgc->gcsl_start4(p,s->cbe,1);
}

// net exchange with the bed after the concentration step: rate of the parcel account per column
void CPM::susp_flux(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s)
{
    SLICELOOP4
    blSr(i,j) = 0.0;
    
    for(i=0;i<p->knox;++i)
    for(j=0;j<p->knoy;++j)
    {
        int k1 = susp_cell(p,a,i,j);
        
        if(k1<0)
        continue;
        
        double c1 = MAX(a->conc(i,j,k1), 0.0);
        
        blSr(i,j) = p->Q65*s->ws*(s->cbe(i,j) - c1)*p->DXN[IP]*p->DYN[JP];
    }
    
    // column mass of the concentration after the step
    for(i=0;i<p->knox;++i)
    for(j=0;j<p->knoy;++j)
    blMc(i,j) = susp_column(p,a,i,j);
    
    blMc_ok=1;
    
    if(getenv("CPM_SUSPDBG")!=nullptr)
    {
        double ex=0.0;
        for(i=0;i<p->knox;++i)
        for(j=0;j<p->knoy;++j)
        ex += blSr(i,j)*p->dt/MAX(p->Q65,1.0e-12);
        ex = pgc->globalsum(ex);
        double v1 = susp_volume(p,a,pgc);
        dbg_ex += ex;
        dbg_dv += v1 - dbg_v0;
        double cmin=0.0, sc=0.0, s1=0.0, sz=0.0, sn=0.0;
        LOOP
        cmin = MIN(cmin, a->conc(i,j,k));
        cmin = pgc->globalmin(cmin);
        for(i=0;i<p->knox;++i)
        for(j=0;j<p->knoy;++j)
        {
            int k1 = susp_cell(p,a,i,j);
            if(k1<0) continue;
            sc += s->cbe(i,j); s1 += a->conc(i,j,k1); sz += a->topo(i,j,k1); sn += 1.0;
        }
        sc = pgc->globalsum(sc); s1 = pgc->globalsum(s1); sz = pgc->globalsum(sz); sn = pgc->globalsum(sn);
        if(p->mpirank==0 && p->count%50==0 && sn>0.0)
        cout<<"CPM susp near bed: mean c_be "<<sc/sn<<" c_1 "<<s1/sn<<" z_1 "<<sz/sn<<endl;
        if(p->mpirank==0 && (p->count%50==0 || fabs(v1-dbg_v0-ex)>1.0e-3*P.ParcelFactor*Vp))
        cout<<"CPM susp balance: step "<<p->count<<" exchange "<<dbg_ex<<" conc change in steps "<<dbg_dv<<" switch "<<dbg_sw<<" step: "<<ex<<" "<<v1-dbg_v0<<" cmin "<<cmin<<endl;
    }
}

// sediment volume in suspension in the column (i,j)
double CPM::susp_column(lexer *p, fdm *a, int ii, int jj)
{
    double vol=0.0;
    
    for(int kk=0; kk<p->knoz; ++kk)
    {
        int ijk = (ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin;
        
        // flow cells only: a cell taken into the fluid's bed (topo < 0) keeps its value until the
        // concentration step sets it to zero
        if(p->flag4[ijk]>0 && a->phi(ii,jj,kk)>=0.0 && a->topo(ii,jj,kk)>0.0)
        vol += MAX(a->conc(ii,jj,kk),0.0)*p->DXN[ii+marge]*p->DYN[jj+marge]*p->DZN[kk+marge];
    }
    
    return vol;
}

// sediment volume in suspension (for the log)
double CPM::susp_volume(lexer *p, fdm *a, ghostcell *pgc)
{
    double vol=0.0;
    
    if(p->Q58>=3 && p->S10==1)
    LOOP
    if(a->phi(i,j,k)>=0.0 && a->topo(i,j,k)>0.0)
    vol += MAX(a->conc(i,j,k),0.0)*p->DXN[IP]*p->DYN[JP]*p->DZN[KP];
    
    return pgc->globalsum(vol);
}
