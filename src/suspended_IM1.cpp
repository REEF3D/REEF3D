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

#include"suspended_IM1.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"diffusion.h"
#include"ioflow.h"
#include"solver.h"
#include"sediment_fdm.h"

suspended_IM1::suspended_IM1(lexer* p) : concn(p),wvel(p)
{
	gcval_susp=60;
}

suspended_IM1::~suspended_IM1()
{
}

void suspended_IM1::start(fdm* a, lexer* p, convection* pconvec, diffusion* pdiff, solver* psolv, ghostcell* pgc, ioflow* pflow, sediment_fdm *s)
{
    starttime=pgc->timer();
    clearrhs(p,a);
    fill_wvel(p,a,pgc,s);
    pconvec->start(p,a,a->conc,4,a->u,a->v,wvel);
	pdiff->idiff_scalar(p,a,pgc,psolv,a->conc,a->eddyv,1.0,1.0);
	suspsource(p,a,a->conc,s);
	timesource(p,a,a->conc);
    bcsusp_start(p,a,pgc,s,a->conc);
	psolv->start(p,a,pgc,a->conc,a->rhsvec,4);
	sedfsf(p,a,a->conc);
	pgc->start4(p,a->conc,gcval_susp);
    fillconc(p,a,pgc,s);
	p->susptime=pgc->timer()-starttime;
	p->suspiter=p->solveriter;
	if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"suspiter: "<<p->suspiter<<"  susptime: "<<setprecision(3)<<p->susptime<<endl;
}

void suspended_IM1::timesource(lexer* p, fdm* a, field& fn)
{
    int count=0;
    int q;

    LOOP
    {
        a->M.p[count]+= 1.0/p->dt;

        a->rhsvec.V[count] += a->L(i,j,k) + a->conc(i,j,k)/p->dt;

	++count;
    }
}

void suspended_IM1::ctimesave(lexer *p, fdm* a)
{
    LOOP
    concn(i,j,k)=a->conc(i,j,k);
}

void suspended_IM1::fill_wvel(lexer *p, fdm* a, ghostcell *pgc, sediment_fdm *s)
{
    double ws_eff,nval,Re_p;
    
    Re_p = s->ws*p->S20/p->W2;
    nval = (4.7 + 0.41*pow(Re_p,0.75))/(1 + 0.175*pow(Re_p,0.75));
    
    WLOOP
    {
    
    ws_eff = s->ws * pow(MAX(1.0 - a->conc(i,j,k)/0.635, 0.0), nval);
    wvel(i,j,k) = a->w(i,j,k) - ws_eff;
    }
    
    pgc->start3(p,wvel,12);
}

void suspended_IM1::suspsource(lexer* p,fdm* a,field& conc, sediment_fdm *s)
{    
    double zdist;
    
    count=0;
    LOOP
    {
        // exchange with the bed only where the Exner equation applies it (DFBED>0, erodible window
        // S 71 - S 72) and in water: erosion into an air cell is deleted by sedfsf()
        if(p->DF[IJK]>0 && s->DFBED[IJ]>0 && p->XP[IP]>=p->S71 && p->XP[IP]<=p->S72 && a->phi(i,j,k)>=0.0)
        if(a->topo(i,j,k)>0.0 && a->topo(i,j,k-1)<0.0)
        {
        zdist = p->DZN[KP];
        
        a->rhsvec.V[count]  += (-s->ws)*(-s->cbe(i,j))/zdist;
        a->M.p[count] += (s->ws)/zdist;
        
        
        //a->rhsvec.V[count]  += s->ws*s->cbe(i,j)/(zdist);
        }
	++count;
    }
}

void suspended_IM1::bcsusp_start(lexer* p, fdm* a,ghostcell *pgc, sediment_fdm *s, field& conc)
{
    double cval;
    
    // The bed exchanges with the concentration only through suspsource (erosion w_s c_be, deposition
    // w_s c_1), which the Exner equation (susp_ED) and the CPM hybrid (Q 58 3, 4) apply to the bed.
    // The faces to the cells of the bed (topo < 0), of solid bodies, to the air and to the walls of the
    // domain are closed: no settling, convective or diffusive flux (the outflow and diffusion terms of
    // the face are taken out of the diagonal); in- and outflow boundaries stay open. With the bed face
    // open, the deposition left the water a second time by settling through the face, and the bed never
    // received it: most of the sand eroded into suspension was lost (pipeline scour: 73-91 %).
    // Above a fixed, non-erodible bed (no exchange) the sand stays in suspension.
    
    auto closed = [&](int ii, int jj, int kk)
    {
        // walls at the bottom, top and the y sides of the domain (in- and outflow in x stay open)
        if((kk<0 && p->nb5<0 && p->periodic3==0) || (kk>=p->knoz && p->nb6<0 && p->periodic3==0))
        return true;
        
        if(p->j_dir==1 && ((jj<0 && p->nb3<0 && p->periodic2==0) || (jj>=p->knoy && p->nb2<0 && p->periodic2==0)))
        return true;
        
        // bed, solid bodies, air above the free surface
        return a->topo(ii,jj,kk)<0.0 || (p->solidread>0 && a->solid(ii,jj,kk)<0.0) || a->phi(ii,jj,kk)<0.0;
    };
    
        n=0;
        LOOP
        {
            {
            const double dx=p->DXN[IP], dy=p->DYN[JP], dz=p->DZN[KP];
            double vel;
            
            if(closed(i-1,j,k))
            {
            vel = a->u(i-1,j,k);
            a->M.p[n] += (a->M.s[n] + MAX(vel,0.0)/dx) + MIN(vel,0.0)/dx;
            a->M.s[n] = 0.0;
            }
            
            if(closed(i+1,j,k))
            {
            vel = a->u(i,j,k);
            a->M.p[n] += (a->M.n[n] - MIN(vel,0.0)/dx) - MAX(vel,0.0)/dx;
            a->M.n[n] = 0.0;
            }
            
            if(p->j_dir==1 && closed(i,j-1,k))
            {
            vel = a->v(i,j-1,k);
            a->M.p[n] += (a->M.e[n] + MAX(vel,0.0)/dy) + MIN(vel,0.0)/dy;
            a->M.e[n] = 0.0;
            }
            
            if(p->j_dir==1 && closed(i,j+1,k))
            {
            vel = a->v(i,j,k);
            a->M.p[n] += (a->M.w[n] - MIN(vel,0.0)/dy) - MAX(vel,0.0)/dy;
            a->M.w[n] = 0.0;
            }
            
            if(closed(i,j,k-1))
            {
            vel = wvel(i,j,k-1);
            a->M.p[n] += (a->M.b[n] + MAX(vel,0.0)/dz) + MIN(vel,0.0)/dz;
            a->M.b[n] = 0.0;
            }
            
            if(closed(i,j,k+1))
            {
            vel = wvel(i,j,k);
            a->M.p[n] += (a->M.t[n] - MIN(vel,0.0)/dz) - MAX(vel,0.0)/dz;
            a->M.t[n] = 0.0;
            }
            }
            
            if(p->flag4[Im1JK]<0 || (p->DF[IJK]>0 && p->DF[Im1JK]<0))
            {
            a->rhsvec.V[n] -= a->M.s[n]*conc(i-1,j,k);
            a->M.s[n] = 0.0;
            }
            
            if(p->flag4[Ip1JK]<0 || (p->DF[IJK]>0 && p->DF[Ip1JK]<0))
            {
            a->rhsvec.V[n] -= a->M.n[n]*conc(i+1,j,k);
            a->M.n[n] = 0.0;
            }
            
            if((p->flag4[IJm1K]<0 || (p->DF[IJK]>0 && p->DF[IJm1K]<0)) && p->j_dir==1)
            {
            a->rhsvec.V[n] -= a->M.e[n]*conc(i,j-1,k);
            a->M.e[n] = 0.0;
            }
            
            if((p->flag4[IJp1K]<0 || (p->DF[IJK]>0 && p->DF[IJp1K]<0)) && p->j_dir==1)
            {
            a->rhsvec.V[n] -= a->M.w[n]*conc(i,j+1,k);
            a->M.w[n] = 0.0;
            }
            
            if(p->flag4[IJKm1]<0 || (p->DF[IJK]>0 && p->DF[IJKm1]<0))
            {
            a->rhsvec.V[n] -= a->M.b[n]*conc(i,j,k-1);
            a->M.b[n] = 0.0;
            }
            
            if(p->flag4[IJKp1]<0 || (p->DF[IJK]>0 && p->DF[IJKp1]<0))
            {
            a->rhsvec.V[n] -= a->M.t[n]*conc(i,j,k+1);
            a->M.t[n] = 0.0;
            }

        ++n;
        }
        
        
        n=0;
        BASELOOP
        {
            if(p->DF[IJK]<0)
            {
            a->M.p[n] = 1.0;

            a->M.n[n] = 0.0;
            a->M.s[n] = 0.0;

            a->M.w[n] = 0.0;
            a->M.e[n] = 0.0;

            a->M.t[n] = 0.0;
            a->M.b[n] = 0.0;
            
            a->rhsvec.V[n] = 0.0;
            }
        ++n;
        }
}

void suspended_IM1::fillconc(lexer* p, fdm* a, ghostcell *pgc, sediment_fdm *s)
{
    // near-bed concentration of the first fluid cell above the bed (s->bedk), the cell of the
    // bed exchange in suspsource(). (The loop over the solid-forcing cells gcdf4 kept the last
    // entry of each column: next to structures and walls a cell high up the structure.)
    // The reference concentration cbe is given at the same cell centre (bedconc_VR), so cb is the
    // cell value; S 61 2 (Rouse transfer to z = 2 d50) would compare two different levels.
    SLICELOOP4
    {
    s->cb(i,j) = 0.0;
    
    k = s->bedk(i,j);
    
    if(k>=0 && k<p->knoz)
    if(a->phi(i,j,k)>=0.0)
    s->cb(i,j) = MAX(a->conc(i,j,k),0.0);
    }
    
    pgc->gcsl_start4(p,s->cb,1);
}

double suspended_IM1::Rouse_formula(lexer* p, fdm *a, sediment_fdm *s, double Cc)
{
    double Ca;    
    double za,zc,P;
    
    za = 2.0*p->S20;
    
    zc = 0.5*p->DZN[KP]*p->WL[IJ];
    
    P = s->ws/(0.4* (s->shearvel_eff(i,j)>0.0?s->shearvel_eff(i,j):1.0e-6) );
    
    P = MAX(P,0.8);
    P = MIN(P,2.5);
    
    
    Ca = Cc * pow( ((p->WL[IJ]-za)/za) / ((p->WL[IJ]-zc)/zc), P);
    
    Ca = MIN(Ca,0.1);
    
    //cout<<"Cc: "<<Cc<<" Ca: "<<Ca<<" | P: "<<P<<" "<<s->shearvel_eff(i,j)<<endl;

    return Ca;
}

void suspended_IM1::sedfsf(lexer* p,fdm* a,field& conc)
{
    LOOP
    if(a->phi(i,j,k)<0.0)
    conc(i,j,k)=0.0;
}

void suspended_IM1::clearrhs(lexer* p, fdm* a)
{
    n=0;
    LOOP
    {
    a->M.p[n] = 0.0;

    a->M.n[n] = 0.0;
    a->M.s[n] = 0.0;

    a->M.w[n] = 0.0;
    a->M.e[n] = 0.0;

    a->M.t[n] = 0.0;
    a->M.b[n] = 0.0;
            
    a->rhsvec.V[n] = 0.0;
    a->L(i,j,k)=0.0;
    
	++n;
    }
}
