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


#include"nhflow_thinbody.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"slice.h"
#include<cmath>

nhflow_thinbody::nhflow_thinbody(lexer *p, fdm_nhf *d, ghostcell *ppgc) : etaL(p), fp(p), cbx(p), cby(p), phi(p), rL(p), Hx(p), Hy(p), first(true), nlowg(0), nlow(0)
{
    pgc = ppgc;
    ncell = p->imax*p->jmax*(p->kmax+2);

    
    p->Darray(sideC,ncell);
    p->Darray(bx,ncell);
    p->Darray(by,ncell);
    p->Darray(bz,ncell);
    p->Darray(ux,ncell);
    p->Darray(uy,ncell);
    p->Darray(uz,ncell);
    p->Darray(low,ncell);
    p->Darray(fz,ncell);
    
    SLICELOOP4
    {
    etaL(i,j) = d->eta(i,j);
    fp(i,j) = 0.0;
    cbx(i,j) = 0.0;
    cby(i,j) = 0.0;
    phi(i,j) = 0.0;
    }
}

nhflow_thinbody::~nhflow_thinbody()
{
    delete [] sideC;
    delete [] bx;
    delete [] by;
    delete [] bz;
    delete [] ux;
    delete [] uy;
    delete [] uz;
    delete [] low;
    delete [] fz;
}

void nhflow_thinbody::begin(lexer *p)
{
    for(int q=0; q<ncell; ++q)
    bx[q]=by[q]=bz[q]=ux[q]=uy[q]=uz[q]=0.0;
}

void nhflow_thinbody::finish(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    pgc->start4V(p,bx,1);
    pgc->start4V(p,by,1);
    pgc->start4V(p,bz,1);
    pgc->start4V(p,ux,1);
    pgc->start4V(p,uy,1);
    pgc->start4V(p,uz,1);
    
    // vertical momentum faces between cell centres on opposite sides; lower cells: odd number of them above the
    // cell in its column (sigma columns are not split between the subdomains)
    for(int q=0; q<ncell; ++q)
    low[q]=fz[q]=0.0;
    
    int n=0;
    
    ILOOP
    JLOOP
    {
        fp(i,j) = 0.0;
        
        if(p->wet[IJ]==0)
        continue;
        
        for(k=0; k<p->knoz-1; ++k)
        if(sideC[IJK]*sideC[IJKp1]<0.0)
        fz[IJK]=1.0;
        
        double par=0.0;
        
        for(k=p->knoz-1; k>=0; --k)
        {
            if(k<p->knoz-1 && fz[IJK]>0.5)
            par = 1.0-par;
            
            low[IJK] = par;
            
            if(par>0.5)
            {
            fp(i,j) = 1.0;
            ++n;
            }
        }
    }
    
    nlow = n;
    nlowg = pgc->globalisum(n);
    
    // columns with a blocked link to their +x / +y neighbour (at any level)
    SLICELOOP4
    {
    cbx(i,j)=0.0;
    cby(i,j)=0.0;
    }
    
    LOOP
    {
        if(bx[IJK]>0.5)
        cbx(i,j)=1.0;
        
        if(by[IJK]>0.5)
        cby(i,j)=1.0;
    }
    
    pgc->gcsl_start4(p,cbx,1);
    pgc->gcsl_start4(p,cby,1);
    
    pgc->start4V(p,low,1);
    pgc->start4V(p,fz,1);
    pgc->gcsl_start4(p,fp,1);
}

int nhflow_thinbody::nlower(lexer *p, ghostcell *pgc)
{
    return pgc->globalisum(nlow);
}

double nhflow_thinbody::head(lexer *p, fdm_nhf *d, int i, int j, int k)
{
    return low[IJK]>0.5 ? etaL(i,j) : d->eta(i,j);
}

void nhflow_thinbody::update_etaL(lexer *p, fdm_nhf *d, int iter)
{
    // harmonic extension of eta into the footprint columns (Gauss-Seidel, warm start, Dirichlet eta around it)
    SLICELOOP4
    if(fp(i,j)<0.5)
    etaL(i,j) = d->eta(i,j);
    
    pgc->gcsl_start4(p,etaL,50);
    
    for(int it=0; it<iter; ++it)
    {
        SLICELOOP4
        if(fp(i,j)>0.5 && p->wet[IJ]==1)
        {
            double s=0.0, w=0.0;
            
            auto add = [&](int ii, int jj, double dl)
            {
                const int qq = (ii-p->imin)*p->jmax + (jj-p->jmin);
                
                if(p->wet[qq]==0)
                return;
                
                const double v = fp(ii,jj)>0.5 ? etaL(ii,jj) : d->eta(ii,jj);
                s += v/(dl*dl);
                w += 1.0/(dl*dl);
            };
            
            add(i+1,j,p->DXP[IP]);
            add(i-1,j,p->DXP[IM1]);
            
            if(p->j_dir==1)
            {
            add(i,j+1,p->DYP[JP]);
            add(i,j-1,p->DYP[JM1]);
            }
            
            if(w>0.0)
            etaL(i,j) = s/w;
        }
        
        pgc->gcsl_start4(p,etaL,50);
    }
}

void nhflow_thinbody::flux_hook(lexer *p, fdm_nhf *d, int id, int ipol, double *Fx, double *Fy)
{
    // Faces handled here: blocked faces (wall flux), faces next to lower cells (head eta_L) and faces whose
    // reconstruction stencil (WENO5: two cells to either side, eta per column) reaches across the body ("near"
    // faces). The latter get a local flux from states reconstructed only from cells on the same side (van Leer,
    // slope 0 next to a blocked link), central hydrostatic flux of the cell heads and HLL dissipation: gravity wave
    // speed for the continuity and the normal momentum, flow speed for the tangential momentum (as HLLC). Below a
    // floor (no free surface there): flow speed only, continuity the Rhie-Chow flux of the projection with the cell
    // water depths and the lid correction.
    const double g = fabs(p->W22);
    double *R = ipol==1 ? d->F : (ipol==2 ? d->G : d->H);
    
    if(ipol==1)
    {
        update_etaL(p,d,first?2000:20);
        first=false;
    }
    
    auto Dc = [&](int ii, int jj) {return MAX(d->eta(ii,jj) + d->depth(ii,jj), p->A544);};
    
    auto vl = [](double a, double b)
    {
        const double den = fabs(a) + fabs(b);
        return den>1.0e-20 ? (a*fabs(b) + fabs(a)*b)/den : 0.0;
    };
    
    // advected quantity: U, V or W
    const double *Q = ipol==1 ? d->U : (ipol==2 ? d->V : d->W);
    
    // ---------------------------------------------------------------- x faces (i+1/2), i = -1 .. knox-1
    for(i=-1; i<p->knox; ++i)
    for(j=0; j<p->knoy; ++j)
    {
        bool near=false;
        
        for(int o=-2; o<=2; ++o)
        if(cbx(i+o,j)>0.5)
        near=true;
        
        for(k=0; k<p->knoz; ++k)
        {
            if(p->flag4[IJK]<=0 || p->flag4[Ip1JK]<=0 || p->wet[IJ]==0 || p->wet[Ip1J]==0)
            continue;
            
            const bool blk = bx[IJK]>0.5;
            const bool L = low[IJK]>0.5 || low[Ip1JK]>0.5;
            
            if(!blk && !L && !near)
            continue;
            
            const double Da = Dc(i,j), Db = Dc(i+1,j);
            const double hf = 0.5*(d->depth(i,j) + d->depth(i+1,j));
            const double ea = head(p,d,i,j,k), eb = head(p,d,i+1,j,k);
            const double Ua = d->U[IJK], Ub = d->U[Ip1JK];
            
            if(ipol==4)
            {
                if(blk)
                Fx[IJK] = 0.5*(Da+Db)*ux[IJK];
                else
                if(L)
                {
                    // below a floor: Rhie-Chow face velocity of the projection (link form, nhflow_membrane_rc_flux)
                    const double dxs = p->DXP[IM1], dxc = p->DXP[IP], dxn = p->DXP[IP1];
                    const double dUi  = (dxc*d->MRCX[IJK]   + dxs*d->MRCX[Im1JK])/(dxc+dxs);
                    const double dUi1 = (dxn*d->MRCX[Ip1JK] + dxc*d->MRCX[IJK])/(dxn+dxc);
                    
                    Fx[IJK] = 0.5*(Da+Db)*(0.5*(Ua+Ub) + d->MRCX[IJK] - 0.5*(dUi+dUi1));
                }
                else
                {
                    // near the body: HLL with states from the own side, gravity wave speed
                    const bool om = bx[Im1JK]<0.5 && p->flag4[Im1JK]>0 && p->wet[Im1J]==1;
                    const bool op = bx[Ip1JK]<0.5 && p->flag4[Ip2JK]>0 && p->wet[Ip2J]==1;
                    const double dc = (Ub-Ua)/p->DXP[IP];
                    const double UL = Ua + 0.5*p->DXN[IP]*(om ? vl(dc,(Ua-d->U[Im1JK])/p->DXP[IM1]) : 0.0);
                    const double UR = Ub - 0.5*p->DXN[IP1]*(op ? vl((d->U[Ip2JK]-Ub)/p->DXP[IP1],dc) : 0.0);
                    const double c = sqrt(g*MAX(Da,Db)) + MAX(fabs(UL),fabs(UR));
                    
                    Fx[IJK] = 0.5*(Da*UL + Db*UR) - 0.5*c*(eb - ea);
                }
                continue;
            }
            
            const double Qa = Q[IJK], Qb = Q[Ip1JK];
            
            if(blk)
            {
                if(ipol==1)
                {
                    const double Fa = g*(0.5*ea*ea + ea*hf) + Da*ux[IJK]*ux[IJK];
                    const double Fb = g*(0.5*eb*eb + eb*hf) + Db*ux[IJK]*ux[IJK];
                    
                    Fx[IJK] = Fa;
                    
                    if(i+1<p->knox)
                    R[Ip1JK] += (Fb - Fa)/p->DXN[IP1];
                }
                else
                Fx[IJK] = 0.5*(Da+Db)*ux[IJK]*0.5*(Qa+Qb);
                
                continue;
            }
            
            // states from the cells on either side, slopes not across blocked links
            const bool om = bx[Im1JK]<0.5 && p->flag4[Im1JK]>0 && p->wet[Im1J]==1;
            const bool op = bx[Ip1JK]<0.5 && p->flag4[Ip2JK]>0 && p->wet[Ip2J]==1;
            
            auto states = [&](const double *F, double &fl, double &fr)
            {
                const double dc = (F[Ip1JK]-F[IJK])/p->DXP[IP];
                const double sa = om ? vl(dc,(F[IJK]-F[Im1JK])/p->DXP[IM1]) : 0.0;
                const double sb = op ? vl((F[Ip2JK]-F[Ip1JK])/p->DXP[IP1],dc) : 0.0;
                fl = F[IJK]   + 0.5*p->DXN[IP]*sa;
                fr = F[Ip1JK] - 0.5*p->DXN[IP1]*sb;
            };
            
            double UL,UR,QL,QR;
            states(d->U,UL,UR);
            states(Q,QL,QR);
            
            // dissipation: flow speed; normal momentum near the body (not below a floor) also the gravity wave speed
            const double sx = (ipol==1 && !L) ? sqrt(g*MAX(Da,Db)) + MAX(fabs(UL),fabs(UR)) : fabs(0.5*(UL+UR));
            
            Fx[IJK] = 0.5*(Da*UL*QL + Db*UR*QR) - 0.5*sx*(Db*QR - Da*QL);
            
            if(ipol==1)
            {
                const double ef = 0.5*(ea+eb);
                Fx[IJK] += g*(0.5*ef*ef + ef*hf);
            }
        }
    }
    
    // ---------------------------------------------------------------- y faces (j+1/2), j = -1 .. knoy-1
    if(p->j_dir==1)
    for(i=0; i<p->knox; ++i)
    for(j=-1; j<p->knoy; ++j)
    {
        bool near=false;
        
        for(int o=-2; o<=2; ++o)
        if(cby(i,j+o)>0.5)
        near=true;
        
        for(k=0; k<p->knoz; ++k)
        {
            if(p->flag4[IJK]<=0 || p->flag4[IJp1K]<=0 || p->wet[IJ]==0 || p->wet[IJp1]==0)
            continue;
            
            const bool blk = by[IJK]>0.5;
            const bool L = low[IJK]>0.5 || low[IJp1K]>0.5;
            
            if(!blk && !L && !near)
            continue;
            
            const double Da = Dc(i,j), Db = Dc(i,j+1);
            const double hf = 0.5*(d->depth(i,j) + d->depth(i,j+1));
            const double ea = head(p,d,i,j,k), eb = head(p,d,i,j+1,k);
            const double Va = d->V[IJK], Vb = d->V[IJp1K];
            
            if(ipol==4)
            {
                if(blk)
                Fy[IJK] = 0.5*(Da+Db)*uy[IJK];
                else
                if(L)
                {
                    const double dys = p->DYP[JM1], dyc = p->DYP[JP], dyn = p->DYP[JP1];
                    const double dVj  = (dyc*d->MRCY[IJK]   + dys*d->MRCY[IJm1K])/(dyc+dys);
                    const double dVj1 = (dyn*d->MRCY[IJp1K] + dyc*d->MRCY[IJK])/(dyn+dyc);
                    
                    Fy[IJK] = 0.5*(Da+Db)*(0.5*(Va+Vb) + d->MRCY[IJK] - 0.5*(dVj+dVj1));
                }
                else
                {
                    const bool om = by[IJm1K]<0.5 && p->flag4[IJm1K]>0 && p->wet[IJm1]==1;
                    const bool op = by[IJp1K]<0.5 && p->flag4[IJp2K]>0 && p->wet[IJp2]==1;
                    const double dc = (Vb-Va)/p->DYP[JP];
                    const double VL = Va + 0.5*p->DYN[JP]*(om ? vl(dc,(Va-d->V[IJm1K])/p->DYP[JM1]) : 0.0);
                    const double VR = Vb - 0.5*p->DYN[JP1]*(op ? vl((d->V[IJp2K]-Vb)/p->DYP[JP1],dc) : 0.0);
                    const double c = sqrt(g*MAX(Da,Db)) + MAX(fabs(VL),fabs(VR));
                    
                    Fy[IJK] = 0.5*(Da*VL + Db*VR) - 0.5*c*(eb - ea);
                }
                continue;
            }
            
            const double Qa = Q[IJK], Qb = Q[IJp1K];
            
            if(blk)
            {
                if(ipol==2)
                {
                    const double Fa = g*(0.5*ea*ea + ea*hf) + Da*uy[IJK]*uy[IJK];
                    const double Fb = g*(0.5*eb*eb + eb*hf) + Db*uy[IJK]*uy[IJK];
                    
                    Fy[IJK] = Fa;
                    
                    if(j+1<p->knoy)
                    R[IJp1K] += (Fb - Fa)/p->DYN[JP1];
                }
                else
                Fy[IJK] = 0.5*(Da+Db)*uy[IJK]*0.5*(Qa+Qb);
                
                continue;
            }
            
            const bool om = by[IJm1K]<0.5 && p->flag4[IJm1K]>0 && p->wet[IJm1]==1;
            const bool op = by[IJp1K]<0.5 && p->flag4[IJp2K]>0 && p->wet[IJp2]==1;
            
            auto states = [&](const double *F, double &fl, double &fr)
            {
                const double dc = (F[IJp1K]-F[IJK])/p->DYP[JP];
                const double sa = om ? vl(dc,(F[IJK]-F[IJm1K])/p->DYP[JM1]) : 0.0;
                const double sb = op ? vl((F[IJp2K]-F[IJp1K])/p->DYP[JP1],dc) : 0.0;
                fl = F[IJK]   + 0.5*p->DYN[JP]*sa;
                fr = F[IJp1K] - 0.5*p->DYN[JP1]*sb;
            };
            
            double VL,VR,QL,QR;
            states(d->V,VL,VR);
            states(Q,QL,QR);
            
            const double sy = (ipol==2 && !L) ? sqrt(g*MAX(Da,Db)) + MAX(fabs(VL),fabs(VR)) : fabs(0.5*(VL+VR));
            
            Fy[IJK] = 0.5*(Da*VL*QL + Db*VR*QR) - 0.5*sy*(Db*QR - Da*QL);
            
            if(ipol==2)
            {
                const double ef = 0.5*(ea+eb);
                Fy[IJK] += g*(0.5*ef*ef + ef*hf);
            }
        }
    }
    
    if(ipol==4)
    {
        if(nlowg>0)
        lid_correction(p,d,Fx,Fy);
        
        return;
    }
    
    // ---------------------------------------------------------------- vertical faces, bed slope term below floors
    // vertical faces within two faces of a blocked one: first-order upwind (the reconstruction in z reaches
    // across the floor)
    LOOP
    {
        if(k<p->knoz-1)
        {
            if(fz[IJK]>0.5)
            d->Fz[IJK] = 0.0;
            else
            {
                bool nearz=false;
                
                for(int o=-2; o<=2; ++o)
                if(k+o>=0 && k+o<p->knoz-1 && fz[IJK+o]>0.5)
                nearz=true;
                
                if(nearz)
                {
                    const double om = d->omegaF[FIJKp1];
                    d->Fz[IJK] = om*(om>=0.0 ? Q[IJK] : Q[IJKp1]);
                }
            }
        }
        
        // pjm::upgrad / vpgrad add eta g (dfx_i - dfx_i-1)/dx with eta; below a floor the head is eta_L
        if(low[IJK]>0.5 && p->wet[IJ]==1)
        {
            const double de = etaL(i,j) - d->eta(i,j);
            
            if(ipol==1)
            R[IJK] += de*g*(0.5*(d->depth(i+1,j)-d->depth(i-1,j)))/p->DXN[IP];
            
            if(ipol==2 && p->j_dir==1)
            R[IJK] += de*g*(0.5*(d->depth(i,j+1)-d->depth(i,j-1)))/p->DYN[JP];
        }
    }
}

void nhflow_thinbody::projection_rhs(lexer *p, fdm_nhf *d, double *U, double *V, double *W, double alpha)
{
    // the face velocity of a blocked link in the divergence of nhflow_pjm::rhs is the wall velocity: row n (node k at
    // the bottom of cell k) averages the faces of the cells k (weight 1 - fac) and k-1 (fac)
    const double adt = alpha*p->dt;
    int n=0;
    
    LOOP
    {
        WETDRYDEEP
        {
            for(int c=0; c<2; ++c)
            {
                const int kk = k-c;
                
                if(kk<0)
                continue;
                
                double facx = p->DZN[KM1]/(p->DZN[KP]+p->DZN[KM1]);
                double facy = facx;
                
                if(k==0)
                {
                facx = MAX((1.0 - p->A522*fabs(d->Bx(i,j))),0.0)*p->DZN[KM1]/(p->DZN[KP]+p->DZN[KM1]);
                facy = MAX((1.0 - p->A522*fabs(d->By(i,j))),0.0)*p->DZN[KM1]/(p->DZN[KP]+p->DZN[KM1]);
                }
                
                const double wx = c==0 ? 1.0-facx : facx;
                const double wy = c==0 ? 1.0-facy : facy;
                const double cx = 2.0/((p->DXP[IP]+p->DXP[IM1])*adt);
                const double cy = 2.0/((p->DYP[JP]+p->DYP[JM1])*adt);
                
                const int q   = (i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + kk-p->kmin;
                const int qe  = (i+1-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + kk-p->kmin;
                const int qw  = (i-1-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + kk-p->kmin;
                
                if(bx[q]>0.5)
                d->rhsvec.V[n] += wx*(0.5*(U[q]+U[qe]) - ux[q])*cx;
                
                if(bx[qw]>0.5)
                d->rhsvec.V[n] += wx*(ux[qw] - 0.5*(U[qw]+U[q]))*cx;
                
                if(p->j_dir==1)
                {
                const int qn  = (i-p->imin)*p->jmax*p->kmax + (j+1-p->jmin)*p->kmax + kk-p->kmin;
                const int qs  = (i-p->imin)*p->jmax*p->kmax + (j-1-p->jmin)*p->kmax + kk-p->kmin;
                
                if(by[q]>0.5)
                d->rhsvec.V[n] += wy*(0.5*(V[q]+V[qn]) - uy[q])*cy;
                
                if(by[qs]>0.5)
                d->rhsvec.V[n] += wy*(uy[qs] - 0.5*(V[qs]+V[q]))*cy;
                }
            }
        }
        
        ++n;
    }
}

void nhflow_thinbody::cut_forcing(lexer *p, fdm_nhf *d, double *WH, slice &WL)
{
    LOOP
    if(bz[IJK]>0.5 && p->wet[IJ]==1)
    {
        d->W[IJK] = uz[IJK];
        WH[IJK] = uz[IJK]*WL(i,j);
    }
    
    pgc->start4V(p,d->W,12);
}

double nhflow_thinbody::link_force(lexer *p, fdm_nhf *d, int dir, int i, int j, int k)
{
    const double rg = p->W1*fabs(p->W22);
    const double *P = d->P;
    
    if(dir==0)
    {
        const double pa = rg*head(p,d,i,j,k)   + 0.5*(P[FIJK]+P[FIJKp1]);
        const double pb = rg*head(p,d,i+1,j,k) + 0.5*(P[FIp1JK]+P[FIp1JKp1]);
        const double D = 0.5*(d->eta(i,j)+d->depth(i,j) + d->eta(i+1,j)+d->depth(i+1,j));
        
        return (pa-pb)*p->DYN[JP]*p->DZN[KP]*D;
    }
    
    if(dir==1)
    {
        const double pa = rg*head(p,d,i,j,k)   + 0.5*(P[FIJK]+P[FIJKp1]);
        const double pb = rg*head(p,d,i,j+1,k) + 0.5*(P[FIJp1K]+P[FIJp1Kp1]);
        const double D = 0.5*(d->eta(i,j)+d->depth(i,j) + d->eta(i,j+1)+d->depth(i,j+1));
        
        return (pa-pb)*p->DXN[IP]*p->DZN[KP]*D;
    }
    
    // vertical link through the cut cell k: the cell centres on either side of the floor
    int kb=k-1, ka=k+1;
    
    if(k>0 && fz[IJKm1]>0.5)
    {
    kb=k-1;
    ka=k;
    }
    else
    if(fz[IJK]>0.5)
    {
    kb=k;
    ka=k+1;
    }
    
    kb = MAX(kb,0);
    ka = MIN(ka,p->knoz-1);
    
    const double pb = rg*head(p,d,i,j,kb) + P[FIJK];
    const double pa = rg*head(p,d,i,j,ka) + P[FIJKp1];
    
    return (pb-pa)*p->DXN[IP]*p->DYN[JP];
}

void nhflow_thinbody::lid_correction(lexer *p, fdm_nhf *d, double *Fx, double *Fy)
{
    // Rigid lid below a closed floor. The column of a footprint carries the inner free surface, the cells below its
    // floor are outside water under a lid: their net horizontal flux has to vanish (fixed floor), otherwise water
    // moves between the lower region and the inner free surface. The continuity fluxes of the lower cells are not
    // exactly divergence free in the column sense (the projection works on node control volumes), so they are
    // corrected with the gradient of a lid potential phi, the same on all open lower layers of a face:
    //     div_h( H_f grad phi ) = r      in the footprint columns,  phi = 0 in the columns around it,
    // r: net outflow of the lower cells of the column (sum over layers, DZN weighted), H_f: sum of DZN of the open
    // lower layers of a face. The imbalance goes to the outside columns, the inner level only sees the inner fluxes.
    auto isL = [&](int q) {return low[q]>0.5;};
    
    SLICELOOP4
    {
    rL(i,j)=0.0;
    Hx(i,j)=0.0;
    Hy(i,j)=0.0;
    }
    
    // open lower layers of the faces i+1/2 (i = -1 .. knox-1) and j+1/2
    for(i=-1; i<p->knox; ++i)
    for(j=0; j<p->knoy; ++j)
    for(k=0; k<p->knoz; ++k)
    if((isL(IJK) || isL(Ip1JK)) && bx[IJK]<0.5 && p->flag4[IJK]>0 && p->flag4[Ip1JK]>0 && p->wet[IJ]==1 && p->wet[Ip1J]==1)
    Hx(i,j) += p->DZN[KP];
    
    if(p->j_dir==1)
    for(i=0; i<p->knox; ++i)
    for(j=-1; j<p->knoy; ++j)
    for(k=0; k<p->knoz; ++k)
    if((isL(IJK) || isL(IJp1K)) && by[IJK]<0.5 && p->flag4[IJK]>0 && p->flag4[IJp1K]>0 && p->wet[IJ]==1 && p->wet[IJp1]==1)
    Hy(i,j) += p->DZN[KP];
    
    // net outflow of the lower cells
    LOOP
    if(isL(IJK) && p->wet[IJ]==1)
    rL(i,j) += p->DZN[KP]*((Fx[IJK] - Fx[Im1JK])/p->DXN[IP] + (p->j_dir==1 ? (Fy[IJK] - Fy[IJm1K])/p->DYN[JP] : 0.0));
    
    // SOR, warm start
    const double w = 1.6;
    auto ph = [&](int ii, int jj) {return fp(ii,jj)>0.5 ? phi(ii,jj) : 0.0;};
    
    for(int it=0; it<(first_lid?400:60); ++it)
    {
        SLICELOOP4
        if(fp(i,j)>0.5 && p->wet[IJ]==1)
        {
            const double ae = Hx(i,j)/(p->DXP[IP]*p->DXN[IP]);
            const double aw = Hx(i-1,j)/(p->DXP[IM1]*p->DXN[IP]);
            const double an = p->j_dir==1 ? Hy(i,j)/(p->DYP[JP]*p->DYN[JP]) : 0.0;
            const double as = p->j_dir==1 ? Hy(i,j-1)/(p->DYP[JM1]*p->DYN[JP]) : 0.0;
            const double sa = ae+aw+an+as;
            
            if(sa<=0.0)
            {
                phi(i,j)=0.0;
                continue;
            }
            
            const double pn = (ae*ph(i+1,j) + aw*ph(i-1,j) + (p->j_dir==1 ? an*ph(i,j+1) + as*ph(i,j-1) : 0.0) - rL(i,j))/sa;
            phi(i,j) += w*(pn - phi(i,j));
        }
        
        pgc->gcsl_start4(p,phi,1);
    }
    
    first_lid=false;
    
    // correction on the open lower layers
    for(i=-1; i<p->knox; ++i)
    for(j=0; j<p->knoy; ++j)
    if(Hx(i,j)>0.0)
    {
        const double c = -(ph(i+1,j) - ph(i,j))/p->DXP[IP];
        
        for(k=0; k<p->knoz; ++k)
        if((isL(IJK) || isL(Ip1JK)) && bx[IJK]<0.5 && p->flag4[IJK]>0 && p->flag4[Ip1JK]>0)
        Fx[IJK] += c;
    }
    
    if(p->j_dir==1)
    for(i=0; i<p->knox; ++i)
    for(j=-1; j<p->knoy; ++j)
    if(Hy(i,j)>0.0)
    {
        const double c = -(ph(i,j+1) - ph(i,j))/p->DYP[JP];
        
        for(k=0; k<p->knoz; ++k)
        if((isL(IJK) || isL(IJp1K)) && by[IJK]<0.5 && p->flag4[IJK]>0 && p->flag4[IJp1K]>0)
        Fy[IJK] += c;
    }
}
