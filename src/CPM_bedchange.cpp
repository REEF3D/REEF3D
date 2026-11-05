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
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"CPM.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

// bed interface from the particles: iso-surface of the solid volume fraction theta_bed = Q 26 (1-S 24)
// first order level set estimate, reinitialised afterwards by reinitopo
void CPM::topo_update(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s)
{
    // sub-grid bedload layer (Q 58): the parcels keep the iso-surface (exposure, near-bed closure),
    // the fluid and the layer get the single valued bed surface from the solid volume of each column
    if(p->Q58>0 && p->S10==1 && zsplit==0)
    {
        topo_iso(p,pgc,Tiso);
        topo_column(p,a,pgc);
    }
    
    else
    topo_iso(p,pgc,a->topo);
}

void CPM::topo_iso(lexer *p, ghostcell *pgc, field &f)
{
    double gx,gy,gz,grad,h,val;
    double Tm,Tp;
    
    BASELOOP
    {
        h = p->DZN[KP];
        
        Tm = wallcell(p,i-1,j,k) ? Ts(i,j,k) : Ts(i-1,j,k);
        Tp = wallcell(p,i+1,j,k) ? Ts(i,j,k) : Ts(i+1,j,k);
        gx = (Tp-Tm)/(p->DXP[IM1]+p->DXP[IP]);
        
        gy=0.0;
        if(p->j_dir==1)
        {
        Tm = wallcell(p,i,j-1,k) ? Ts(i,j,k) : Ts(i,j-1,k);
        Tp = wallcell(p,i,j+1,k) ? Ts(i,j,k) : Ts(i,j+1,k);
        gy = (Tp-Tm)/(p->DYP[JM1]+p->DYP[JP]);
        }
        
        Tm = wallcell(p,i,j,k-1) ? (1.0-p->S24) : Ts(i,j,k-1);
        Tp = wallcell(p,i,j,k+1) ? Ts(i,j,k) : Ts(i,j,k+1);
        gz = (Tp-Tm)/(p->DZP[KM1]+p->DZP[KP]);
        
        grad = sqrt(gx*gx + gy*gy + gz*gz);
        grad = MAX(grad, theta_bed/(2.0*h));
        
        val = (theta_bed - Ts(i,j,k))/grad;
        
        val = MAX(val,-3.0*h);
        val = MIN(val, 3.0*h);
        
        f(i,j,k) = val;
    }

    pgc->start4a(p,f,150);
}

// bed level set seen by the parcels: the iso-surface of the solid fraction
double CPM::ptopo(lexer *p, fdm *a, double xp, double yp, double zp)
{
    if(p->Q58>0 && p->S10==1 && zsplit==0)
    return cip4a(p,Tiso,xp,yp,zp);

    return p->ccipol4_b(a->topo,xp,yp,zp);
}

double CPM::ptopo_p(lexer *p, fdm *a, int q, double xp, double yp, double zp)
{
    if(q<0)
    return ptopo(p,a,xp,yp,zp);
    
    if(size_t(q)>=tpV.size())
    {
        const size_t N = MAX(size_t(P.index), size_t(q)+1);
        tpX.resize(N,0.0); tpY.resize(N,0.0); tpZ.resize(N,0.0); tpV.resize(N,0.0);
        tpS.resize(N,-1);
    }
    
    if(tpS[q]==p->count && tpX[q]==xp && tpY[q]==yp && tpZ[q]==zp)
    return tpV[q];
    
    const double v = ptopo(p,a,xp,yp,zp);
    
    tpX[q]=xp; tpY[q]=yp; tpZ[q]=zp; tpV[q]=v; tpS[q]=p->count;
    
    return v;
}

// bed elevation from the topo level set
void CPM::bedzh_update(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s)
{
    double h;
    
    ILOOP
    JLOOP
    {
        h = s->bedzh(i,j);
        
        if(a->topo(i,j,0)>=0.0 && p->nb5<0)
        h = p->ZN[0+marge];
        
        KLOOP
        PBASECHECK
        if(k>0 && a->topo(i,j,k-1)<0.0 && a->topo(i,j,k)>=0.0)
            h = -(a->topo(i,j,k-1)*p->DZP[KM1])/(a->topo(i,j,k)-a->topo(i,j,k-1)) + p->pos_z()-p->DZP[KM1];
            
        s->bedzh(i,j)=h;
        a->bed(i,j)=h;
    }
    
    pgc->gcsl_start4(p,s->bedzh,1);
    pgc->gcsl_start4(p,a->bed,50);
}

/*--------------------------------------------------------------------
bed level from the solid volume of the columns (Q 58, S 10 1, no vertical decomposition)

The iso-surface of the solid fraction depends on where the parcels sit inside the top cell
of the bed, so the bed level jumps by up to a cell between neighbouring columns and the
fluid sees a staircase (form drag, sheltered columns without transport). With the bedload
layer the bed surface is single valued: from the bottom face of the highest cell of the bed
(theta >= theta_bed), that cell and the one above count with their solid volume (parcels in
their nearest cell, without the bedload layer and the suspension),

    z_b = z_face + (theta_k/theta_0) dz_k + (theta_k+1/theta_0) dz_k+1,     topo = z - z_b.

Loose cells further down do not lower the surface, parcels further up (suspension) do not count.
--------------------------------------------------------------------*/
void CPM::topo_column(lexer *p, fdm *a, ghostcell *pgc)
{
    const int nj = p->knoy+2;
    const int nk = p->knoz;
    const double vpar = P.ParcelFactor*Vp;
    std::vector<double> zb((p->knox+2)*nj, 0.0);
    
    // solid volume of the bed per cell: parcels in their nearest cell (no kernel spread of the bed
    // into the cell above), without the parcels of the bedload layer and the suspension
    // (moving faster than 0.1 w_s above the bed surface of the parcels)
    std::vector<double> occ(size_t(p->knox)*p->knoy*nk, 0.0);
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]>=ACTIVE && P.Hop[n]<=0.0)
    {
        int ic = pci(p,n,P.X[n]);
        int jc = p->j_dir==1 ? pcj(p,n,P.Y[n]) : 0;
        int kc = pck(p,n,P.Z[n]);
        
        if(ic<0 || ic>=p->knox || jc<0 || jc>=p->knoy || kc<0 || kc>=nk)
        continue;
        
        double sp = sqrt(P.U[n]*P.U[n] + P.V[n]*P.V[n] + P.W[n]*P.W[n]);
        
        if(sp>0.1*settling_velocity(p,P.D[n]) && ptopo(p,a,P.X[n],P.Y[n],P.Z[n])>=0.0)
        continue;
        
        occ[(size_t(ic)*p->knoy + jc)*nk + kc] += vpar*P.Wt[n];
    }

    for(i=0;i<p->knox;++i)
    for(j=0;j<p->knoy;++j)
    {
        const double A = p->DXN[IP]*p->DYN[JP];
        const double *oc = &occ[(size_t(i)*p->knoy + j)*nk];
        
        // highest cell of the bed
        int ktop=-1;

        for(k=0;k<nk;++k)
        {
            // the bed ends at a solid body above it (parcels at or in a pipeline are no bed level);
            // solid cells below the bed (a fixed floor) are skipped
            if(p->solidread>0 && a->solid(i,j,k)<0.0)
            {
                if(ktop>=0)
                break;
                
                continue;
            }
            
            if(oc[k] >= theta_bed*A*p->DZN[KP])
            ktop=k;
            
            // the bed ends at the first almost empty cell above it: single parcels resting higher up in
            // the flow (on a structure, in a recirculation, at rest for a moment in the wake) are no bed
            // level; on fine grids two slow parcels fill a cell to theta_bed (a pipeline wake on 2.5 mm
            // cells gave a bed 7 cm above the real one).
            // An almost empty cell with a packed cell within the next two cells above it is a cavity in the
            // bed (the grid-limited step leaves single cells of a scoured bed with one or two parcels), not
            // its end: the column went on (pipeline, 2.5 mm: single columns of the bed level lay up to 8
            // cells below their neighbours, under 4-10 packed cells above such a cavity)
            else if(ktop>=0 && oc[k] < 0.25*theta_bed*A*p->DZN[KP])
            {
                bool cavity=false;
                
                for(int kq=k+1; kq<=MIN(k+2,nk-1); ++kq)
                {
                    if(p->solidread>0 && a->solid(i,j,kq)<0.0)
                    break;
                    
                    if(oc[kq] >= theta_bed*A*p->DZN[kq+marge])
                    cavity=true;
                }
                
                if(!cavity)
                break;
            }
        }

        // the cells below the highest cell of the bed count as full (dilated or loose cells inside
        // the bed do not lower the surface), the highest cell and the one above with their solid volume
        // no bed cell (only a solid floor): the top of the floor
        int kb = MAX(ktop,0);
        
        if(ktop<0 && p->solidread>0)
        for(k=0;k<nk && a->solid(i,j,k)<0.0;++k)
        kb = k+1;
        
        kb = MIN(kb,nk-1);
        double z = p->ZN[kb+marge];

        for(k=kb;k<=MIN(kb+1,nk-1);++k)
        if(!(p->solidread>0 && a->solid(i,j,k)<0.0))
        z += MIN(oc[k]/(theta_0*A*p->DZN[KP]), 1.0)*p->DZN[KP];

        // the fluid sees the bed level relaxed in time (Q 63, default 10 s): single parcels moving in and
        // out of the top cells (pickup, deposition from the suspension) make the level jump by
        // V_p/(theta_0 A) each time; a bed boundary that crosses the cell centres back and forth
        // switches fluid cells on and off and stalls the near-bed flow (bed shear stress down to 1/5)
        blZb(i,j) = z;
        zbl_ok = 1;
        
        if(zbf_ini==0)
        zbf(i,j) = z;
        
        // Q 63 is a time of the bed (morphological time, Q 65 t): with the flow time the level seen by the fluid
        // lagged the parcels by Q 65 Q 63 of bed time (pipeline, Q 65 50, Q 63 2: 100 s; the fluid bed 1 cm above
        // the parcels under the pipe and up to 1 cm below them on the mound, where hundreds of deposited grains
        // stood in the resolved flow and were carried into the wake)
        const double mfz = p->Q65>0.0 ? p->Q65 : 1.0;
        
        zbf(i,j) += MIN(1.0, mfz*p->dt/MAX(p->Q63,1.0e-12))*(z - zbf(i,j));
        
        zb[(i+1)*nj + j+1] = zbf(i,j);
    }
    
    // Q 72: the level seen by the fluid is smoothed over the neighbouring columns (1-2-1 filter, Q 72 passes).
    // A fluid cell is bed or fluid by the sign of topo at its centre: a single column whose level passes the
    // centre by a fraction of a cell (parcel noise) blocks the whole cell, and the one-cell bump stalls the
    // near-bed flow behind it (periodic channel, bed level at a cell centre: bed shear stress down by 40 %
    // for 2 s, w = 0.13 m/s at the bump; with Q 72 2 no cell switched and no stall)
    if(p->Q72>0)
    {
        for(i=0;i<p->knox;++i)
        for(j=0;j<p->knoy;++j)
        zbt(i,j) = zbf(i,j);
        
        std::vector<double> tmp(size_t(p->knox)*p->knoy);
        
        for(int qn=0;qn<p->Q72;++qn)
        {
            pgc->gcsl_start4(p,zbt,1);
            
            for(i=0;i<p->knox;++i)
            for(j=0;j<p->knoy;++j)
            {
                // the domain sides (no neighbour subdomain, not periodic) mirror the column: the ghost values of
                // the slice are no bed level there (an inflow boundary pulled the first column down by 1/4 of
                // the bed height per pass)
                const double c  = zbt(i,j);
                const double xm = (i==0 && p->nb1<0 && perx!=1) ? c : zbt(i-1,j);
                const double xp = (i==p->knox-1 && p->nb4<0 && perx!=1) ? c : zbt(i+1,j);
                
                double v = 0.25*xm + 0.5*c + 0.25*xp;
                
                if(p->j_dir==1)
                {
                    const double ym = (j==0 && p->nb3<0 && pery!=1) ? c : zbt(i,j-1);
                    const double yp = (j==p->knoy-1 && p->nb2<0 && pery!=1) ? c : zbt(i,j+1);
                    
                    v = 0.5*(v + 0.25*ym + 0.5*c + 0.25*yp);
                }
                
                tmp[size_t(i)*p->knoy + j] = v;
            }
            
            for(i=0;i<p->knox;++i)
            for(j=0;j<p->knoy;++j)
            zbt(i,j) = tmp[size_t(i)*p->knoy + j];
        }
        
        // the filter removes the noise, not the bed forms: the level of a column moves by at most half a cell
        // (the edge of the initial hole under a pipeline, a step of several cells, was smeared over the
        // neighbouring columns, opened cells next to the pipe and made the flow unstable)
        for(i=0;i<p->knox;++i)
        for(j=0;j<p->knoy;++j)
        {
            const int kc = MAX(0, MIN(nk-1, p->posc_k(zbf(i,j))));
            const double dmax = 0.5*p->DZN[kc+marge];
            
            zb[(i+1)*nj + j+1] = zbf(i,j) + MAX(-dmax, MIN(dmax, zbt(i,j) - zbf(i,j)));
        }
    }
    
    zbf_ini=1;
    
    // layer bed level of the neighbour subdomains: the slope in CPM_bedload reads blZb(-1,j) etc.
    // at subdomain boundaries (was never exchanged, 0 there: spurious slope at every boundary)
    pgc->gcsl_start4(p,blZb,1);
    
    BASELOOP
    {
        double h = p->DZN[KP];
        double val = p->ZP[KP] - zb[(i+1)*nj + j+1];

        a->topo(i,j,k) = MAX(-3.0*h, MIN(3.0*h, val));
    }

    pgc->start4a(p,a->topo,150);
}
