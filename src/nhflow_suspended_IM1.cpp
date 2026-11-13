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

#include"nhflow_suspended_IM1.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"sediment_inflow.h"
#include"nhflow_scalar_convection.h"
#include"nhflow_diffusion.h"
#include"ioflow.h"
#include"solver.h"
#include"sediment_fdm.h"
#include"ckpt_registry.h"

nhflow_suspended_IM1::nhflow_suspended_IM1(lexer* p) : ED(p),RT(p)
{
	gcval_susp=60;

    p->Darray(WVEL,p->imax*p->jmax*(p->kmax+2));
    p->Darray(WLN,p->imax*p->jmax);
    wl_ini=0;
    
    SLICEBASELOOP
    ED(i,j)=RT(i,j)=0.0;
    
    acc[0]=acc[1]=0.0;
    ckid = ckpt_registry::add_raw(p,acc,2,false);
}

nhflow_suspended_IM1::~nhflow_suspended_IM1()
{
    ckpt_registry::remove(ckid);
}

void nhflow_suspended_IM1::start(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_scalar_convection *pconvec, nhflow_diffusion *pdiff, solver *psolv, ioflow *pflow, sediment_fdm *s)
{
    starttime=pgc->timer();
    drysave(p,d,s);
    clearrhs(p,d);
    fill_wvel(p,d,pgc,s);
    inflow_conc(p,d,s);
    pconvec->start(p,d,d->CONC,4,d->U,d->V,WVEL);
    pdiff->diff_scalar(p,d,pgc,psolv,d->CONC,1.0,1.0);
    settling_fit(p,d,s);
	suspsource(p,d,d->CONC,s);
	timesource(p,d,d->CONC);
    bcsusp_start(p,d,pgc,s,d->CONC);
    psolv->startV(p,pgc,d->CONC,d->rhsvec,d->M,4);
	pgc->start60V(p,d->CONC,gcval_susp);
    fillconc(p,d,pgc,s);
    
    // the exchange with the bed of this solve, for the bed update (S 62 1)
    exchange_sum(p,d,s);
    
    // depth belonging to the new concentration, for the D^n/D^(n+1) ratio of the next time step
    SLICELOOP4
    WLN[IJ] = d->WL(i,j);
    wl_ini=1;
    
	p->susptime=pgc->timer()-starttime;
	p->suspiter=p->solveriter;
	if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"suspiter: "<<p->suspiter<<"  susptime: "<<setprecision(3)<<p->susptime<<endl;
    
    // sediment balance (every P 12 steps): suspended volume, bed change (solids) since the start, the sediment
    // of dried columns not yet deposited and the exchange since the last bed update, which the bed takes in the
    // next sediment step (S 62 1); closed domains: the sum stays constant when the sediment time step is the flow
    // time since the last bed update (S 15 4, or S 15 3 with S 42 1, S 44 1) and S 17 1, S 35 1
    if(p->count%p->P12==0)
    {
    double vol=0.0, bed=0.0, dry=0.0, pend=0.0;
    
    SLICELOOP4
    pend += ED(i,j)*p->DXN[IP]*p->DYN[JP];
    
    LOOP
    if(p->DF[IJK]>0 && p->wet[IJ]==1)
    vol += MAX(d->CONC[IJK],0.0)*p->DZN[KP]*MAX(d->WL(i,j),0.0)*p->DXN[IP]*p->DYN[JP];
    
    SEDSLICELOOP
    {
    bed += (1.0-p->S24)*(s->bedzh(i,j)-s->bedzh0(i,j))*p->DXN[IP]*p->DYN[JP];
    dry += s->dryd(i,j)*p->DXN[IP]*p->DYN[JP];
    }
    
    vol = pgc->globalsum(vol);
    bed = pgc->globalsum(bed);
    dry = pgc->globalsum(dry);
    pend = pgc->globalsum(pend);
    
    if(p->mpirank==0)
    cout<<"suspended load (S 12): volume "<<setprecision(8)<<vol<<" m^3, bed change "<<bed<<" m^3, dried columns "<<dry
        <<" m^3, pending exchange "<<-pend<<" m^3, balance "<<vol+bed+dry-pend<<" m^3"<<endl;
    }
}

void nhflow_suspended_IM1::timesource(lexer* p, fdm_nhf *d, double *FN)
{
    // conservative form for D*C with the transport coefficients divided by D^(n+1):
    //   (D^(n+1) C^(n+1) - D^n C^n)/(D^(n+1) dt) = C^(n+1)/dt - (D^n/D^(n+1)) C^n/dt
    // D^n: depth at the end of the last solve (includes bed changes and the flow step since then)
    int count=0;
    double hn,ho;

    LOOP
    {
        hn = MAX(d->WL(i,j),1.0e-20);
        ho = wl_ini==1?MAX(WLN[IJ],0.0):hn;
        
        d->M.p[count]+= 1.0/p->dt;

        d->rhsvec.V[count] += d->L[IJK] + (ho/hn)*d->CONC[IJK]/p->dt;

	++count;
    }
}

void nhflow_suspended_IM1::ctimesave(lexer *p, fdm_nhf *d)
{
}

void nhflow_suspended_IM1::drysave(lexer *p, fdm_nhf *d, sediment_fdm *s)
{
    // columns that fell dry since the last solve: the solve sets their concentration to zero,
    // so the sediment they still hold (D^n C^n) is handed to the bed (deposited by the Exner step)
    if(wl_ini==0)
    return;
    
    SLICELOOP4
    if(p->wet[IJ]==0)
    {
        KLOOP
        PCHECK
        s->dryd(i,j) += MAX(d->CONC[IJK],0.0)*p->DZN[KP]*MAX(WLN[IJ],0.0);
    }
}

void nhflow_suspended_IM1::fill_wvel(lexer *p, fdm_nhf *d, ghostcell *pgc, sediment_fdm *s)
{
    // WVEL: vertical volume flux across the sigma faces for the conservative (form=2) ifou scheme:
    // omegaF is the face volume flux D*dsigma/dt [m/s],
    // settling across a sigma face is -ws (D*dsigma/dt of -ws = -ws).
    // Face k lies between cells k-1 and k; bed (k=0) and surface (k=knoz) faces stay closed,
    // bed exchange is handled by suspsource().
    // S 127 1: hindered settling (Richardson-Zaki), in the bed exchange as well (suspsource); S 127 0: ws
    double ws_eff,nval,Re_p,cface;

    Re_p = s->ws*p->S20/p->W2;
    nval = (4.7 + 0.41*pow(Re_p,0.75))/(1 + 0.175*pow(Re_p,0.75));

    FLOOP
    {
    WVEL[FIJK] = 0.0;

        if(k>0 && k<p->knoz && p->DF[IJK]>0 && p->wet[IJ]==1)
        {
        cface = 0.5*(d->CONC[IJK] + d->CONC[IJKm1]);
        ws_eff = s->ws;
        
        if(p->S127==1)
        ws_eff *= pow(MAX(1.0 - cface/0.635, 0.0), nval);
        WVEL[FIJK] = d->omegaF[FIJK] - ws_eff;
        }
    }
}

void nhflow_suspended_IM1::settling_fit(lexer *p, fdm_nhf *d, sediment_fdm *s)
{
    // Vertical settling and diffusion as one exponentially fitted flux (Allen & Southwell 1955, Il'in 1969; the
    // exponential scheme of Patankar). The first-order upwind flux of WVEL carries the numerical diffusivity
    // |w| dz/2; near the bed that is as large as the eddy diffusivity kappa u* z itself (ws dz/2 against
    // kappa u* dz/2 in the bottom cell), and the profile above the bed came out several times fuller than the
    // Rouse profile on every practical grid. The diffusivity of each face between two cells is reduced by
    //     eps (1 + x - x coth x),   x = |w| dz/(2 eps),
    // so that upwind plus diffusion is the exact flux of the steady settling-diffusion balance with eps and w of
    // the face: eps - |w| dz/2 for x << 1 (second order), no diffusion left for x >> 1 (upwind). The reduced
    // diffusivity is never negative: the matrix stays an M-matrix. Bed and free surface faces are closed.
    //
    // The face between the bottom cell and the cell above: the bottom cell carries the mean of the near-bed
    // profile (bedconc_VR), which is well above the profile at its centre (a factor 2 for a Rouse number of 1.8
    // and a cell of five reference heights). Taken as the value at the centre it fed the cells above with that
    // factor. The flux of this face sees the bottom cell as a node at rbz = max(centre, a) with the value
    // rbf c_0, the Rouse profile at rbz: a (c_1 - rbf c_0), with a over the distance from rbz to the centre
    // above. The same flux leaves the bottom cell and enters the cell above (conservative).
    auto fit = [](double eps, double wd)
    {
        if(eps<=1.0e-20 || wd<=0.0)
        return 0.0;
        
        const double x = 0.5*wd/eps;
        
        if(x>30.0)
        return -eps;
        
        if(x<1.0e-4)
        return -eps*x;
        
        const double e = exp(-2.0*x);
        
        return eps*(x*(1.0+e)/(1.0-e) - 1.0 - x);
    };
    
    double eps,c;
    
    n=0;
    LOOP
    {
        if(p->DF[IJK]>0 && p->wet[IJ]==1)
        {
        const double epsP = d->VISC[IJK] + d->EV0[IJK];
        const double sz2 = p->sigz[IJ]*p->sigz[IJ];
        const double sxyz2 = pow(0.5*(p->sigx[FIJK]+p->sigx[FIJKp1]),2.0) + pow(0.5*(p->sigy[FIJK]+p->sigy[FIJKp1]),2.0) + sz2;
        const bool bedcell = p->DFBED[IJ]>0 && p->knoz>1;
        
            // the face between the bottom cell and the cell above, seen from either cell
            if(bedcell && (k==0 || k==1) && p->DF[k==0 ? IJKp1 : IJKm1]>0)
            {
            const int kb = k;
            k=0;
            const double zc = p->ZP[KP]*d->WL(i,j);
            const double dzs = p->DZP[KP];
            const double wface = fabs(WVEL[FIJKp1]);
            eps = 0.5*(d->VISC[IJK] + d->EV0[IJK] + d->VISC[IJKp1] + d->EV0[IJKp1]);
            k=kb;
            
            // sigma distance from the node of the bottom cell to the centre above
            const double dzn = MAX(dzs - MAX(s->rbz(i,j)-zc,0.0)/MAX(d->WL(i,j),1.0e-20), 0.25*dzs);
            const double phi = s->rbf(i,j);
            
            const double a_old = eps*sxyz2/(dzs*p->DZN[KP]);
            const double a_new = (eps*sxyz2 + fit(eps, wface*dzn*d->WL(i,j))*sz2)/(dzn*p->DZN[KP]);
            
                if(k==0)
                {
                d->M.t[n] -= a_new - a_old;
                d->M.p[n] += a_new*phi - a_old;
                }
                
                if(k==1)
                {
                d->M.b[n] -= a_new*phi - a_old;
                d->M.p[n] += a_new - a_old;
                }
            }
        
            // face above
            if(k<p->knoz-1 && p->DF[IJKp1]>0 && !(bedcell && k==0))
            {
            eps = 0.5*(epsP + d->VISC[IJKp1] + d->EV0[IJKp1]);
            c = fit(eps, fabs(WVEL[FIJKp1])*p->DZP[KP]*d->WL(i,j))*sz2/(p->DZP[KP]*p->DZN[KP]);
            
            d->M.t[n] -= c;
            d->M.p[n] += c;
            }
            
            // face below
            if(k>0 && p->DF[IJKm1]>0 && !(bedcell && k==1))
            {
            eps = 0.5*(epsP + d->VISC[IJKm1] + d->EV0[IJKm1]);
            c = fit(eps, fabs(WVEL[FIJK])*p->DZP[KM1]*d->WL(i,j))*sz2/(p->DZP[KM1]*p->DZN[KP]);
            
            d->M.b[n] -= c;
            d->M.p[n] += c;
            }
        }
    ++n;
    }
}

void nhflow_suspended_IM1::suspsource(lexer* p, fdm_nhf *d, double *CONC, sediment_fdm *s)
{    
    double zdist,wsb;
    
    const double Re_p = s->ws*p->S20/p->W2;
    const double nval = (4.7 + 0.41*pow(Re_p,0.75))/(1 + 0.175*pow(Re_p,0.75));
    
    count=0;
    LOOP
    {   
        // exchange with the bed only where the Exner equation applies it (erodible region S71-S72, DFBED>0),
        // otherwise erosion from or deposition onto a fixed bed creates or deletes sediment
        // E - D = ws (c_a - c_0/rbe) = wsx (cbe - c_0), wsx = ws/rbe: pick-up ws c_a at the reference level z = a
        // and deposition with the near-bed concentration c_0/rbe of the Rouse profile in the bottom cell
        // (bedconc_VR: cbe = rbe c_a is the cell mean). The equilibrium c_0 = cbe is that of ws (cbe - c_0); away
        // from it the exchange does not slow down with the thickness of the bottom cell (rbe << 1 in a thick
        // cell). S 127 1: the fall velocity hindered with the concentration of the cell, as at the faces above.
        if(k==0 && p->DF[IJK]>0 && p->wet[IJ]==1 && p->DFBED[IJ]>0 && p->XP[IP]>=p->S71 && p->XP[IP]<=p->S72)
        {
        zdist = p->DZN[KP]*d->WL(i,j);
        
        wsb = s->ws;
        
        if(p->S127==1)
        wsb *= pow(MAX(1.0 - MAX(d->CONC[IJK],0.0)/0.635, 0.0), nval);
        
        s->wsx(i,j) = wsb/MAX(s->rbe(i,j),1.0e-8);
        
        d->rhsvec.V[count]  += s->wsx(i,j)*s->cbe(i,j)/zdist;
        d->M.p[count] += s->wsx(i,j)/zdist;
        }
        
	++count;
    }
}

void nhflow_suspended_IM1::inflow_conc(lexer *p, fdm_nhf *d, sediment_fdm *s)
{
    // prescribed inflow concentration (S 124 >= 0, or the time series of S 129 2): the inflow ghost columns (IO 1)
    // at the x- edge carry the depth-averaged concentration C_in with the profile of the first column (uniform
    // while that column is clear); bcsusp_start takes them as the upwind (Dirichlet) value of the inflow face.
    // C_in < 0: zero gradient (the face coefficient is folded into the diagonal). All four edges of the domain.
    cin_on = 0;
    
    const double cin = (s->pinflow!=nullptr) ? s->pinflow->conc(p) : p->S124;
    
    if(cin<0.0)
    return;
    
    cin_on = 1;
    
    double cm;
    
    // x-
    if(p->origin_i==0)
    for(j=0; j<p->knoy; ++j)
    {
    i=0;
    cm=0.0;
    
        for(k=0; k<p->knoz; ++k)
        cm += MAX(d->CONC[IJK],0.0)*p->DZN[KP];
    
        for(k=0; k<p->knoz; ++k)
        if(p->IO[Im1JK]==1)
        d->CONC[Im1JK] = d->CONC[Im2JK] = d->CONC[Im3JK] = (cm>1.0e-12) ? cin*MAX(d->CONC[IJK],0.0)/cm : cin;
    }
    
    // x+
    if(p->origin_i+p->knox==p->gknox)
    for(j=0; j<p->knoy; ++j)
    {
    i=p->knox-1;
    cm=0.0;
    
        for(k=0; k<p->knoz; ++k)
        cm += MAX(d->CONC[IJK],0.0)*p->DZN[KP];
    
        for(k=0; k<p->knoz; ++k)
        if(p->IO[Ip1JK]==1)
        d->CONC[Ip1JK] = d->CONC[Ip2JK] = d->CONC[Ip3JK] = (cm>1.0e-12) ? cin*MAX(d->CONC[IJK],0.0)/cm : cin;
    }
    
    // y-
    if(p->j_dir==1 && p->origin_j==0)
    for(i=0; i<p->knox; ++i)
    {
    j=0;
    cm=0.0;
    
        for(k=0; k<p->knoz; ++k)
        cm += MAX(d->CONC[IJK],0.0)*p->DZN[KP];
    
        for(k=0; k<p->knoz; ++k)
        if(p->IO[IJm1K]==1)
        d->CONC[IJm1K] = d->CONC[IJm2K] = d->CONC[IJm3K] = (cm>1.0e-12) ? cin*MAX(d->CONC[IJK],0.0)/cm : cin;
    }
    
    // y+
    if(p->j_dir==1 && p->origin_j+p->knoy==p->gknoy)
    for(i=0; i<p->knox; ++i)
    {
    j=p->knoy-1;
    cm=0.0;
    
        for(k=0; k<p->knoz; ++k)
        cm += MAX(d->CONC[IJK],0.0)*p->DZN[KP];
    
        for(k=0; k<p->knoz; ++k)
        if(p->IO[IJp1K]==1)
        d->CONC[IJp1K] = d->CONC[IJp2K] = d->CONC[IJp3K] = (cm>1.0e-12) ? cin*MAX(d->CONC[IJK],0.0)/cm : cin;
    }
}

void nhflow_suspended_IM1::bcsusp_start(lexer *p, fdm_nhf *d, ghostcell *pgc, sediment_fdm *s, double *CONC)
{
        n=0;
        LOOP
        {
            if(p->DF[IJK]>0 && p->wet[IJ]==1)
            {
            // closed faces (walls, domain edges, solid and dry neighbours, bed, free surface):
            // implicit zero gradient, the off-diagonal is folded into the diagonal, so the
            // diffusive flux through the face is exactly zero (explicit C^n left a flux ~ C^(n+1)-C^n);
            // at inflow edges the inflow concentration is the cell value
            if(cin_on==1 && i+p->origin_i==0 && p->IO[Im1JK]==1 && p->DF[Im1JK]>0)
            {
            d->rhsvec.V[n] -= d->M.s[n]*CONC[Im1JK];
            d->M.s[n] = 0.0;
            }
            
            if(cin_on==1 && i+p->origin_i==p->gknox-1 && p->IO[Ip1JK]==1 && p->DF[Ip1JK]>0)
            {
            d->rhsvec.V[n] -= d->M.n[n]*CONC[Ip1JK];
            d->M.n[n] = 0.0;
            }
            
            if(cin_on==1 && p->j_dir==1 && j+p->origin_j==0 && p->IO[IJm1K]==1 && p->DF[IJm1K]>0)
            {
            d->rhsvec.V[n] -= d->M.e[n]*CONC[IJm1K];
            d->M.e[n] = 0.0;
            }
            
            if(cin_on==1 && p->j_dir==1 && j+p->origin_j==p->gknoy-1 && p->IO[IJp1K]==1 && p->DF[IJp1K]>0)
            {
            d->rhsvec.V[n] -= d->M.w[n]*CONC[IJp1K];
            d->M.w[n] = 0.0;
            }
            
            if(p->flag4[Im1JK]<0 || p->DF[Im1JK]<0 || p->wet[Im1J]==0)
            {
            d->M.p[n] += d->M.s[n];
            d->M.s[n] = 0.0;
            }
            
            if(p->flag4[Ip1JK]<0 || p->DF[Ip1JK]<0 || p->wet[Ip1J]==0)
            {
            d->M.p[n] += d->M.n[n];
            d->M.n[n] = 0.0;
            }
            
            if(p->j_dir==1)
            if(p->flag4[IJm1K]<0 || p->DF[IJm1K]<0 || p->wet[IJm1]==0)
            {
            d->M.p[n] += d->M.e[n];
            d->M.e[n] = 0.0;
            }
            
            if(p->j_dir==1)
            if(p->flag4[IJp1K]<0 || p->DF[IJp1K]<0 || p->wet[IJp1]==0)
            {
            d->M.p[n] += d->M.w[n];
            d->M.w[n] = 0.0;
            }
            
            // bed: no diffusive flux, the exchange with the bed is in suspsource (implicit zero gradient)
            if(k==0)
            {
            d->M.p[n] += d->M.b[n];
            d->M.b[n] = 0.0;
            }
            
            if((p->flag4[IJKm1]<0 || p->DF[IJKm1]<0) && k>0)
            {
            d->M.p[n] += d->M.b[n];
            d->M.b[n] = 0.0;
            }
            
            if((p->flag4[IJKp1]<0 || p->DF[IJKp1]<0) && k<p->knoz-1)
            {
            d->M.p[n] += d->M.t[n];
            d->M.t[n] = 0.0;
            }
            
            // free surface: impermeable, implicit zero gradient (no diffusive loss)
            if(k==p->knoz-1)
            {
            d->M.p[n] += d->M.t[n];
            d->M.t[n] = 0.0;
            }
            }

        ++n;
        }
        
        
    // turn off inside direct forcing body
        n=0;
        LOOP
        {
            if(p->DF[IJK]<0 || p->wet[IJ]==0)
            {
            d->M.p[n] = 1.0;

            d->M.n[n] = 0.0;
            d->M.s[n] = 0.0;

            d->M.w[n] = 0.0;
            d->M.e[n] = 0.0;

            d->M.t[n] = 0.0;
            d->M.b[n] = 0.0;
            
            d->rhsvec.V[n] = 0.0;
            }
        ++n;
        }
}

void nhflow_suspended_IM1::fillconc(lexer* p, fdm_nhf *d, ghostcell *pgc, sediment_fdm *s)
{
    k=0;
    SLICELOOP4
    {
    if(p->DF[IJK]<0 || p->wet[IJ]==0)
    s->cb(i,j) = 0.0;
    
        if(p->DF[IJK]>0 && p->wet[IJ]==1)
        {
            // cell value at the level of cbe (bedconc_VR gives cbe at the first cell centre); no upper
            // cap, the bed exchange in suspsource() uses the unclipped concentration
            s->cb(i,j) = MAX(d->CONC[IJK],0.0);
        }
    }    
    pgc->gcsl_start4(p,s->cb,1);
    
    
    double Uh;
    
    SEDSLICELOOP
    s->qbs(i,j) = 0.0;
    
    SEDSLICELOOP
    KLOOP
    {
        Uh = sqrt(d->U[IJK]*d->U[IJK] + d->V[IJK]*d->V[IJK]);
        
        s->qbs(i,j) += Uh*d->CONC[IJK]*p->DZN[KP]*d->WL(i,j);
    }
    
    pgc->gcsl_start4(p,s->qbs,1);
}

void nhflow_suspended_IM1::clearrhs(lexer* p, fdm_nhf *d)
{
    n=0;
    LOOP
    {    
    d->rhsvec.V[n]=0.0;
    d->L[IJK]=0.0;
    
    
            d->M.p[n] = 0.0;

            d->M.n[n] = 0.0;
            d->M.s[n] = 0.0;

            d->M.w[n] = 0.0;
            d->M.e[n] = 0.0;

            d->M.t[n] = 0.0;
            d->M.b[n] = 0.0;
	++n;
    }
}

void nhflow_suspended_IM1::exchange_sum(lexer *p, fdm_nhf *d, sediment_fdm *s)
{
    // the implicit exchange of the solve, wsx (cbe - c_0) on the cells of suspsource; the bed takes it from the
    // sediment step on (susp_ED: count > S 43)
    if(p->S62!=1 || p->count<=p->S43)
    return;
    
    k=0;
    SLICELOOP4
    if(p->DF[IJK]>0 && p->wet[IJ]==1 && p->DFBED[IJ]>0 && p->XP[IP]>=p->S71 && p->XP[IP]<=p->S72)
    ED(i,j) += s->wsx(i,j)*(s->cbe(i,j) - d->CONC[IJK])*p->dt;
    
    acc[0] += p->dt;
    acc[1] += 1.0;
}

void nhflow_suspended_IM1::bed_rate(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    // one flow step since the last bed update: the bed takes ws (cbn - cb) of that step (susp_ED);
    // several flow steps (S 44 > 1, S 42 2, 3): the mean rate of all of them, so that with S 15 4 the bed takes
    // exactly the exchange of the water column
    s->edrate = nullptr;
    
    if(acc[1]<1.5 || acc[0]<=0.0)
    return;
    
    SLICELOOP4
    RT(i,j) = ED(i,j)/acc[0];
    
    pgc->gcsl_start4(p,RT,1);
    
    s->edrate = &RT;
}

void nhflow_suspended_IM1::bed_done(lexer *p, sediment_fdm *s)
{
    SLICEBASELOOP
    ED(i,j) = 0.0;
    
    acc[0] = acc[1] = 0.0;
    s->edrate = nullptr;
}
