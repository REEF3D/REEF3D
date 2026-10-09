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

nhflow_thinbody::nhflow_thinbody(lexer *p, fdm_nhf *d, ghostcell *ppgc) : etaL(p), fp(p), cbx(p), cby(p), phi(p), rL(p), Hx(p), Hy(p), cr(p), cd(p), cq(p), wf(p), rt(p), first(true), nlowg(0), nlow(0)
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
    p->Darray(pw,ncell);
    p->Darray(side0,ncell);
    p->Darray(tu,ncell);
    p->Darray(tv,ncell);
    p->Darray(sideN,p->imax*p->jmax*p->kmaxF);
    
    for(int q=0; q<ncell; ++q)
    pw[q]=0.5;
    p->Darray(fz,ncell);
    
    SLICELOOP4
    {
    etaL(i,j) = d->eta(i,j);
    fp(i,j) = 0.0;
    cbx(i,j) = 0.0;
    cby(i,j) = 0.0;
    phi(i,j) = 0.0;
    wf(i,j) = 0.0;
    rt(i,j) = 0.0;
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
    delete [] pw;
    delete [] side0;
    delete [] tu;
    delete [] tv;
    delete [] sideN;
    delete [] fz;
}

void nhflow_thinbody::begin(lexer *p)
{
    for(int q=0; q<ncell; ++q)
    bx[q]=by[q]=bz[q]=ux[q]=uy[q]=uz[q]=0.0;
    
    for(int q=0; q<p->imax*p->jmax*p->kmaxF; ++q)
    sideN[q]=0.0;
}

void nhflow_thinbody::finish(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    pgc->start4V(p,bx,1);
    pgc->start4V(p,by,1);
    pgc->start4V(p,bz,1);
    pgc->start4V(p,ux,1);
    pgc->start4V(p,uy,1);
    pgc->start4V(p,uz,1);
    
    // node sides of the neighbouring subdomains (node links across the subdomain border: node_cut_x/y)
    pgc->gcparax7(p,sideN,7);
    pgc->gcparax7co(p,sideN,7);
    
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
    
    // vertical wall velocity of the floor above the lower cells of a column (the vertical node link at the topmost floor
    // crossing, cell just below or just above it): with the normal wall velocities, the inflow through the cut face of the
    // column (the blocked side faces of the lower cells carry their own wall fluxes)
    ILOOP
    JLOOP
    {
        wf(i,j) = 0.0;
        
        if(fp(i,j)<0.5)
        continue;
        
        for(k=p->knoz-2; k>=0; --k)
        if(fz[IJK]>0.5)
        {
            if(bz[IJK]>0.5)
            wf(i,j) = uz[IJK];
            else
            if(bz[IJKp1]>0.5)
            wf(i,j) = uz[IJKp1];
            break;
        }
    }
    
    pgc->gcsl_start4(p,wf,1);
    
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
    
    // cell pressure next to the body: only the node on the cell's own side
    for(int q=0; q<ncell; ++q)
    pw[q]=0.5;
    
    LOOP
    if(sideC[IJK]!=0.0)
    {
        const double s = sideC[IJK];
        const bool ok0 = !(sideN[FIJK]*s<0.0), ok1 = !(sideN[FIJKp1]*s<0.0);
        
        if(!ok0 && ok1)
        pw[IJK]=0.0;
        
        if(ok0 && !ok1)
        pw[IJK]=1.0;
    }
    
    pgc->start4V(p,pw,1);
    
    pgc->start4V(p,low,1);
    pgc->start4V(p,fz,1);
    pgc->gcsl_start4(p,fp,1);
    
    top_region(p,d);
}

void nhflow_thinbody::wall_faces(lexer *p, int i, int j, int k, int *w) const
{
    w[0] = bx[Im1JK]>0.5;
    w[1] = bx[IJK]>0.5;
    w[2] = p->j_dir==1 && by[IJm1K]>0.5;
    w[3] = p->j_dir==1 && by[IJK]>0.5;
    w[4] = k>0 && fz[IJKm1]>0.5;
    w[5] = k<p->knoz-1 && fz[IJK]>0.5;
}

void nhflow_thinbody::matrix_walls(lexer *p, fdm_nhf *d, const double *F)
{
    int w[6];
    int n=0;
    
    LOOP
    {
        if(p->wet[IJ]==1)
        {
            wall_faces(p,i,j,k,w);
            
            if(w[0])
            {
            d->rhsvec.V[n] -= d->M.s[n]*F[IJK];
            d->M.s[n] = 0.0;
            }
            
            if(w[1])
            {
            d->rhsvec.V[n] -= d->M.n[n]*F[IJK];
            d->M.n[n] = 0.0;
            }
            
            if(w[2])
            {
            d->rhsvec.V[n] -= d->M.e[n]*F[IJK];
            d->M.e[n] = 0.0;
            }
            
            if(w[3])
            {
            d->rhsvec.V[n] -= d->M.w[n]*F[IJK];
            d->M.w[n] = 0.0;
            }
            
            if(w[4])
            {
            d->rhsvec.V[n] -= d->M.b[n]*F[IJK];
            d->M.b[n] = 0.0;
            }
            
            if(w[5])
            {
            d->rhsvec.V[n] -= d->M.t[n]*F[IJK];
            d->M.t[n] = 0.0;
            }
        }
        
        ++n;
    }
}

double nhflow_thinbody::inner_volume(lexer *p, fdm_nhf *d)
{
    double V=0.0;
    
    LOOP
    if(p->wet[IJ]==1 && ((rt(i,j)>0.5) != (low[IJK]>0.5)))
    V += p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*d->WL(i,j);
    
    return pgc->globalsum(V);
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
    // head of the cells below a body crossing ("lower" cells, the region other than the one of the column's free surface):
    // columns with the free surface inside the body (rt = 1, bag interior): the outer level, harmonic extension of eta of
    // the outside columns (Dirichlet) over the inside columns; columns with the free surface outside and lower cells (a wall
    // that leans outwards, a dent or fold of a flexible membrane: the lower cells hold inner water): the inner level,
    // harmonic extension of eta of the inside columns over these columns (no flux to other outside columns; isolated:
    // own eta). Gauss-Seidel, warm start.
    SLICELOOP4
    if(!(rt(i,j)>0.5 || fp(i,j)>0.5))
    etaL(i,j) = d->eta(i,j);
    
    pgc->gcsl_start4(p,etaL,50);
    
    for(int it=0; it<iter; ++it)
    {
        SLICELOOP4
        if((rt(i,j)>0.5 || fp(i,j)>0.5) && p->wet[IJ]==1)
        {
            const bool in = rt(i,j)>0.5;
            double s=0.0, w=0.0;
            
            auto add = [&](int ii, int jj, double dl)
            {
                const int qq = (ii-p->imin)*p->jmax + (jj-p->jmin);
                
                if(p->wet[qq]==0)
                return;
                
                const bool nin = rt(ii,jj)>0.5;
                double v;
                
                if(in)
                v = nin ? etaL(ii,jj) : d->eta(ii,jj);
                else
                {
                    if(nin)
                    v = d->eta(ii,jj);
                    else
                    if(fp(ii,jj)>0.5)
                    v = etaL(ii,jj);
                    else
                    return;
                }
                
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
            
            etaL(i,j) = w>0.0 ? s/w : d->eta(i,j);
        }
        
        pgc->gcsl_start4(p,etaL,50);
    }
}

void nhflow_thinbody::top_region(lexer *p, fdm_nhf *d)
{
    // region of the free surface of each column: 1 inside the body (the top cell on the inner side, or connected to such a
    // column by open links of the top layer), 0 outside
    SLICELOOP4
    {
        k = p->knoz-1;
        rt(i,j) = sideC[IJK]<0.0 ? 1.0 : (sideC[IJK]>0.0 ? 0.0 : -1.0);
    }
    
    pgc->gcsl_start4(p,rt,1);
    
    for(int it=0; it<200; ++it)
    {
        int ch=0;
        k = p->knoz-1;
        
        SLICELOOP4
        if(rt(i,j)<-0.5)
        {
            if((rt(i+1,j)>0.5 && bx[IJK]<0.5) || (rt(i-1,j)>0.5 && bx[Im1JK]<0.5)
            || (p->j_dir==1 && ((rt(i,j+1)>0.5 && by[IJK]<0.5) || (rt(i,j-1)>0.5 && by[IJm1K]<0.5))))
            {
                rt(i,j) = 1.0;
                ++ch;
            }
        }
        
        pgc->gcsl_start4(p,rt,1);
        
        if(pgc->globalisum(ch)==0)
        break;
    }
    
    SLICELOOP4
    if(rt(i,j)<-0.5)
    rt(i,j) = 0.0;
    
    pgc->gcsl_start4(p,rt,1);
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
                
                const double sN = nside(p,i,j,k);
                
                if(bx[q]>0.5)
                d->rhsvec.V[n] += wx*(0.5*(U[q]+U[qe]) - ux[q])*cx;
                else
                if(other_side(q,qe,sN))
                d->rhsvec.V[n] += wx*0.5*(U[q]+U[qe])*cx;
                
                if(bx[qw]>0.5)
                d->rhsvec.V[n] += wx*(ux[qw] - 0.5*(U[qw]+U[q]))*cx;
                else
                if(other_side(qw,q,sN))
                d->rhsvec.V[n] -= wx*0.5*(U[qw]+U[q])*cx;
                
                if(p->j_dir==1)
                {
                const int qn  = (i-p->imin)*p->jmax*p->kmax + (j+1-p->jmin)*p->kmax + kk-p->kmin;
                const int qs  = (i-p->imin)*p->jmax*p->kmax + (j-1-p->jmin)*p->kmax + kk-p->kmin;
                
                if(by[q]>0.5)
                d->rhsvec.V[n] += wy*(0.5*(V[q]+V[qn]) - uy[q])*cy;
                else
                if(other_side(q,qn,sN))
                d->rhsvec.V[n] += wy*0.5*(V[q]+V[qn])*cy;
                
                if(by[qs]>0.5)
                d->rhsvec.V[n] += wy*(uy[qs] - 0.5*(V[qs]+V[q]))*cy;
                else
                if(other_side(qs,q,sN))
                d->rhsvec.V[n] -= wy*0.5*(V[qs]+V[q])*cy;
                }
            }
        }
        
        ++n;
    }
}

void nhflow_thinbody::side_change(lexer *p, fdm_nhf *d, double *UH, double *VH, double *WH, slice &WL)
{
    // a cell whose centre crossed the body since the last stage keeps the velocity of the other side; it is replaced
    // by the mean of its face neighbours on its new side that did not cross (no neighbour: at rest)
    auto flipped = [&](int q) {return sideC[q]*side0[q]<0.0;};
    int n=0;
    
    LOOP
    if(p->wet[IJ]==1 && flipped(IJK))
    {
        const double s = sideC[IJK];
        const int nb[6] = {Im1JK,Ip1JK,IJm1K,IJp1K,IJKm1,IJKp1};
        double su=0.0, sv=0.0, sw=0.0;
        int c=0;
        
        for(int r=0; r<6; ++r)
        {
            const int q = nb[r];
            
            if(r>=2 && r<4 && p->j_dir==0)
            continue;
            
            if(r==4 && k==0)
            continue;
            
            if(r==5 && k==p->knoz-1)
            continue;
            
            if(p->flag4[q]<=0 || flipped(q) || sideC[q]*s<0.0)
            continue;
            
            // neighbours far from the body (side 0) are on the same side only if not across a blocked link
            if(sideC[q]==0.0 && ((r==0 && bx[Im1JK]>0.5) || (r==1 && bx[IJK]>0.5) || (r==2 && by[IJm1K]>0.5) || (r==3 && by[IJK]>0.5)))
            continue;
            
            su += d->U[q];
            sv += d->V[q];
            sw += d->W[q];
            ++c;
        }
        
        const double u = c>0 ? su/c : 0.0, v = c>0 ? sv/c : 0.0, w = c>0 ? sw/c : 0.0;
        
        d->U[IJK] = u;
        d->V[IJK] = v;
        d->W[IJK] = w;
        UH[IJK] = u*WL(i,j);
        VH[IJK] = v*WL(i,j);
        WH[IJK] = w*WL(i,j);
        ++n;
    }
    
    for(int q=0; q<ncell; ++q)
    side0[q] = sideC[q];
    
    pgc->start4V(p,d->U,10);
    pgc->start4V(p,d->V,11);
    pgc->start4V(p,d->W,12);
}

void nhflow_thinbody::cut_forcing(lexer *p, fdm_nhf *d, double *UH, double *VH, double *WH, slice &WL)
{
    LOOP
    if(bz[IJK]>0.5 && p->wet[IJ]==1)
    {
        d->W[IJK] = uz[IJK];
        WH[IJK] = uz[IJK]*WL(i,j);
    }
    
    // pockets of the staircase: three or more blocked faces, the cell moves with the body (mean wall velocity of its
    // blocked links)
    LOOP
    if(p->wet[IJ]==1)
    {
        int n=0;
        double u=0.0, v=0.0, w=0.0, nu=0.0, nv=0.0, nw=0.0;
        
        if(bx[Im1JK]>0.5) {++n; u+=ux[Im1JK]; nu+=1.0;}
        if(bx[IJK]>0.5)   {++n; u+=ux[IJK];   nu+=1.0;}
        
        if(p->j_dir==1)
        {
        if(by[IJm1K]>0.5) {++n; v+=uy[IJm1K]; nv+=1.0;}
        if(by[IJK]>0.5)   {++n; v+=uy[IJK];   nv+=1.0;}
        }
        
        if(k>0 && fz[IJKm1]>0.5) {++n; w+=wf(i,j); nw+=1.0;}
        if(k<p->knoz-1 && fz[IJK]>0.5) {++n; w+=wf(i,j); nw+=1.0;}
        
        if(n<3)
        continue;
        
        if(bz[IJK]>0.5) {w+=uz[IJK]; nw+=1.0;}
        
        u = nu>0.0 ? u/nu : 0.0;
        v = nv>0.0 ? v/nv : 0.0;
        w = nw>0.0 ? w/nw : 0.0;
        
        d->U[IJK] = u;
        d->V[IJK] = v;
        d->W[IJK] = w;
        UH[IJK] = u*WL(i,j);
        VH[IJK] = v*WL(i,j);
        WH[IJK] = w*WL(i,j);
    }
    
    // cells next to a blocked face: the velocity component normal to the face lies between the wall velocity and the
    // velocity on the other side of the cell (monotone). The projection reaches this component only through the open
    // face (half of the wide gradient), a mode with the cell flowing into the body is not controlled and can grow
    // (flexible bag, a lower cell below a dented floor: 0.25 -> 11 m/s within 0.7 s)
    pgc->start4V(p,d->U,10);
    pgc->start4V(p,d->V,11);
    
    LOOP
    if(p->wet[IJ]==1)
    {
        if(bx[Im1JK]>0.5 || bx[IJK]>0.5)
        {
            const double uw = bx[Im1JK]>0.5 ? ux[Im1JK] : d->U[Im1JK];
            const double ue = bx[IJK]>0.5   ? ux[IJK]   : d->U[Ip1JK];
            const double u = MAX(MIN(uw,ue), MIN(MAX(uw,ue), d->U[IJK]));
            
            UH[IJK] += (u - d->U[IJK])*WL(i,j);
            tu[IJK] = u;
        }
        else
        tu[IJK] = d->U[IJK];
        
        if(p->j_dir==1)
        {
            if(by[IJm1K]>0.5 || by[IJK]>0.5)
            {
                const double vs = by[IJm1K]>0.5 ? uy[IJm1K] : d->V[IJm1K];
                const double vn = by[IJK]>0.5   ? uy[IJK]   : d->V[IJp1K];
                const double v = MAX(MIN(vs,vn), MIN(MAX(vs,vn), d->V[IJK]));
                
                VH[IJK] += (v - d->V[IJK])*WL(i,j);
                tv[IJK] = v;
            }
            else
            tv[IJK] = d->V[IJK];
        }
    }
    
    LOOP
    if(p->wet[IJ]==1)
    {
        d->U[IJK] = tu[IJK];
        
        if(p->j_dir==1)
        d->V[IJK] = tv[IJK];
    }
    
    pgc->start4V(p,d->U,10);
    pgc->start4V(p,d->V,11);
    pgc->start4V(p,d->W,12);
}

double nhflow_thinbody::link_force(lexer *p, fdm_nhf *d, int dir, int i, int j, int k)
{
    const double rg = p->W1*fabs(p->W22);
    const double *P = d->P;
    
    if(dir==0)
    {
        const double pa = rg*head(p,d,i,j,k)   + pcell(p,P,i,j,k);
        const double pb = rg*head(p,d,i+1,j,k) + pcell(p,P,i+1,j,k);
        const double D = 0.5*(d->eta(i,j)+d->depth(i,j) + d->eta(i+1,j)+d->depth(i+1,j));
        
        return (pa-pb)*p->DYN[JP]*p->DZN[KP]*D;
    }
    
    if(dir==1)
    {
        const double pa = rg*head(p,d,i,j,k)   + pcell(p,P,i,j,k);
        const double pb = rg*head(p,d,i,j+1,k) + pcell(p,P,i,j+1,k);
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
    
    // (also the ghost columns i = -1, j = -1: the faces i-1/2, j-1/2 of the first local columns are accumulated below)
    for(i=-1; i<=p->knox; ++i)
    for(j=-1; j<=p->knoy; ++j)
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
    
    // conjugate gradients on the footprint columns (symmetric positive definite -L phi = -r; parallel: the
    // operator only needs the ghost columns of the search direction), warm start from the last stage, relative
    // residual 1e-8
    auto ph = [&](int ii, int jj) {return fp(ii,jj)>0.5 ? phi(ii,jj) : 0.0;};
    auto act = [&](int ii, int jj) {return fp(ii,jj)>0.5 && p->wet[(ii-p->imin)*p->jmax + (jj-p->jmin)]==1;};
    
    // A s = sa s - sum a_nb s_nb on the active columns (s = 0 elsewhere)
    auto apply = [&](slice4 &s, slice4 &As)
    {
        pgc->gcsl_start4(p,s,1);
        
        SLICELOOP4
        {
            As(i,j) = 0.0;
            
            if(!act(i,j))
            continue;
            
            const double ae = Hx(i,j)/(p->DXP[IP]*p->DXN[IP]);
            const double aw = Hx(i-1,j)/(p->DXP[IM1]*p->DXN[IP]);
            const double an = p->j_dir==1 ? Hy(i,j)/(p->DYP[JP]*p->DYN[JP]) : 0.0;
            const double as = p->j_dir==1 ? Hy(i,j-1)/(p->DYP[JM1]*p->DYN[JP]) : 0.0;
            const double sa = ae+aw+an+as;
            
            auto sv = [&](int ii, int jj) {return act(ii,jj) ? s(ii,jj) : 0.0;};
            
            As(i,j) = sa>0.0 ? sa*s(i,j) - (ae*sv(i+1,j) + aw*sv(i-1,j) + (p->j_dir==1 ? an*sv(i,j+1) + as*sv(i,j-1) : 0.0))
                             : s(i,j);
        }
    };
    
    SLICELOOP4
    if(!act(i,j))
    phi(i,j)=0.0;
    
    // r = b - A phi, b = -rL (sa = 0: phi = 0)
    apply(phi,cq);
    
    double rr=0.0, bb=0.0;
    
    SLICELOOP4
            {
        cr(i,j) = 0.0;
        
        if(act(i,j))
        {
            const double sa = Hx(i,j)+Hx(i-1,j)+(p->j_dir==1 ? Hy(i,j)+Hy(i,j-1) : 0.0);
            const double b = sa>0.0 ? -(rL(i,j) + wf(i,j)) : 0.0;
            cr(i,j) = b - cq(i,j);
            rr += cr(i,j)*cr(i,j);
            bb += b*b;
            }
            
        cd(i,j) = cr(i,j);
    }
    
    rr = pgc->globalsum(rr);
    bb = pgc->globalsum(bb);
    
    const double tol2 = 1.0e-16*MAX(bb,1.0e-60);
    
    for(int it=0; it<500 && rr>tol2; ++it)
    {
        apply(cd,cq);
        
        double dq=0.0;
        SLICELOOP4
        if(act(i,j))
        dq += cd(i,j)*cq(i,j);
        
        dq = pgc->globalsum(dq);
        
        if(dq<=0.0)
        break;
        
        const double al = rr/dq;
        double rn=0.0;
        
        SLICELOOP4
        if(act(i,j))
        {
            phi(i,j) += al*cd(i,j);
            cr(i,j)  -= al*cq(i,j);
            rn += cr(i,j)*cr(i,j);
        }
        
        rn = pgc->globalsum(rn);
        
        const double be = rn/rr;
        rr = rn;
        
        SLICELOOP4
        cd(i,j) = act(i,j) ? cr(i,j) + be*cd(i,j) : 0.0;
        }
        
        pgc->gcsl_start4(p,phi,1);
    
    
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
