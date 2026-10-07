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

#include"sediment_mixture.h"
#include"sediment_fdm.h"
#include"lexer.h"
#include"ghostcell.h"
#include"bedload.h"
#include<cmath>
#include<cstdlib>
#include<iostream>
#include<iomanip>
#include<sys/stat.h>

sediment_mixture::sediment_mixture(lexer *p) : Hs(p),d50(p),d90(p),dm(p),qbe_raw(p),qbe_tot(p),fac(p),
                                               shields_eff0(p),shields_crit0(p),tau_crit0(p),shearvel_crit0(p)
{
    nf = p->S51;
    rho_s = p->S22;
    hide_type = p->S54;
    hide_m = p->S55;
    S20 = p->S20;

    if(nf<1 || nf>20)
    {
        if(p->mpirank==0)
        cout<<"sediment_mixture: number of fractions (S 51 lines) must be between 1 and 20, found "<<nf<<endl;
        exit(1);
    }

    if(p->S10==12)
    {
        if(p->mpirank==0)
        cout<<"sediment_mixture: multi-fraction bed (S 51) is not available with the RK2 sediment solver S 10 12"<<endl;
        exit(1);
    }

    d = new double[nf];
    order = new int[nf];
    Fe = new double[nf];
    V = new double[nf];
    FI = new double[nf];
    dhc = new double[nf];
    inS = new int[nf];
    clip_k = new double[nf];
    clip_sum = new double[nf];
    noneq_ini_k = new int[nf];

    F = new slice4*[nf];
    Fs = new slice4*[nf];
    qbe_k = new slice4*[nf];
    qb_k = new slice4*[nf];
    qbn_k = new slice4*[nf];
    vz_k = new slice4*[nf];
    dh_k = new slice4*[nf];
    fh_k = new slice4*[nf];

    for(int q=0;q<nf;++q)
    {
        F[q] = new slice4(p);
        Fs[q] = new slice4(p);
        qbe_k[q] = new slice4(p);
        qb_k[q] = new slice4(p);
        qbn_k[q] = new slice4(p);
        vz_k[q] = new slice4(p);
        dh_k[q] = new slice4(p);
        fh_k[q] = new slice4(p);

        d[q] = p->S51_d[q];
        clip_k[q] = 0.0;
        clip_sum[q] = 0.0;
        noneq_ini_k[q] = 0;

        if(d[q]<=0.0)
        {
            if(p->mpirank==0)
            cout<<"sediment_mixture: grain diameter of fraction "<<q+1<<" must be positive"<<endl;
            exit(1);
        }
    }

    // sandslide transfer to the 8 neighbours, index (di+1)*3+(dj+1)
    out = new slice4*[9];
    for(int q=0;q<9;++q)
    out[q] = new slice4(p);

    // fractions sorted by diameter for the grain size statistics
    for(int q=0;q<nf;++q)
    order[q]=q;

    for(int q=1;q<nf;++q)
    for(int r=q;r>0 && d[order[r]]<d[order[r-1]];--r)
    {
        int tmp=order[r];
        order[r]=order[r-1];
        order[r-1]=tmp;
    }
}

sediment_mixture::~sediment_mixture()
{
    for(int q=0;q<nf;++q)
    {
        delete F[q];
        delete Fs[q];
        delete qbe_k[q];
        delete qb_k[q];
        delete qbn_k[q];
        delete vz_k[q];
        delete dh_k[q];
        delete fh_k[q];
    }

    for(int q=0;q<9;++q)
    delete out[q];

    delete [] F;
    delete [] Fs;
    delete [] qbe_k;
    delete [] qb_k;
    delete [] qbn_k;
    delete [] vz_k;
    delete [] dh_k;
    delete [] fh_k;
    delete [] out;

    delete [] d;
    delete [] order;
    delete [] Fe;
    delete [] V;
    delete [] FI;
    delete [] dhc;
    delete [] clip_k;
    delete [] clip_sum;
    delete [] noneq_ini_k;
}

// --------------------------------------------------------------------
// initialisation
// --------------------------------------------------------------------

void sediment_mixture::ini(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    const int size = p->imax*p->jmax;
    double suma=0.0, sums=0.0;

    for(int q=0;q<nf;++q)
    {
    suma += MAX(p->S51_fa[q],0.0);
    sums += MAX(p->S51_fs[q],0.0);
    }

    if(suma<=0.0)
    {
        if(p->mpirank==0)
        cout<<"sediment_mixture: active layer fractions (S 51, 2nd value) sum to zero"<<endl;
        exit(1);
    }

    // an empty substrate definition means: substrate has the composition of the active layer
    for(int q=0;q<nf;++q)
    {
    Fe[q] = MAX(p->S51_fa[q],0.0)/suma;
    V[q]  = sums>0.0?MAX(p->S51_fs[q],0.0)/sums:Fe[q];
    }

    if(p->mpirank==0 && (fabs(suma-1.0)>1.0e-6 || (sums>0.0 && fabs(sums-1.0)>1.0e-6)))
    cout<<"sediment_mixture: fractions normalised to sum 1 (active: "<<suma<<", substrate: "<<sums<<")"<<endl;

    for(int q=0;q<nf;++q)
    for(int nn=0;nn<size;++nn)
    {
    F[q]->V[nn]  = Fe[q];
    Fs[q]->V[nn] = V[q];
    qbe_k[q]->V[nn] = 0.0;
    qb_k[q]->V[nn]  = 0.0;
    qbn_k[q]->V[nn] = 0.0;
    vz_k[q]->V[nn]  = 0.0;
    dh_k[q]->V[nn]  = 0.0;
    fh_k[q]->V[nn]  = 0.0;
    }

    for(int nn=0;nn<size;++nn)
    {
    Hs.V[nn] = MAX(p->S53,0.0);
    qbe_raw.V[nn] = 0.0;
    qbe_tot.V[nn] = 0.0;
    fac.V[nn] = 1.0;
    }

    grain_stats(p,pgc);

    // active layer thickness
    if(p->S52>0.0)
    La = p->S52;

    if(p->S52<=0.0)
    {
        // d90 of the initial active layer
        double c0,c1,cum=0.0,d90ini=d[order[nf-1]];
        double *c = V;

        for(int r=0;r<nf;++r)
        {
        c[r] = cum + 0.5*Fe[order[r]];
        cum += Fe[order[r]];
        }

        if(0.9<=c[0])
        d90ini = d[order[0]];

        for(int r=0;r<nf-1;++r)
        if(0.9>c[r] && 0.9<=c[r+1])
        {
        c0=c[r];
        c1=c[r+1];
        d90ini = exp(log(d[order[r]]) + (0.9-c0)/MAX(c1-c0,1.0e-20)*(log(d[order[r+1]])-log(d[order[r]])));
        }

        La = 2.0*d90ini;
    }

    // roughness from the local grain size
    if(p->S56>0)
    SLICELOOP4
    {
    s->ks(i,j) = p->S21*ks_diameter(p,i,j);
    s->ks_eff(i,j) = p->S21*ks_diameter(p,i,j);
    }

    if(p->mpirank==0)
    {
    cout<<"sediment mixture: "<<nf<<" fractions, active layer La = "<<La<<" m, substrate Hs = "<<p->S53<<" m"<<endl;
    cout<<"  hiding/exposure: "<<(p->S54==1?"Wu, Wang & Jia (2000), m = ":"off")<<(p->S54==1?p->S55:0.0)<<endl;

    for(int q=0;q<nf;++q)
    cout<<"  fraction "<<q+1<<"  d = "<<d[q]<<" m   F_active = "<<Fe[q]<<"   F_substrate = "<<Fs[q]->V[0]<<endl;

    mkdir("./REEF3D_Log",0777);
    
    // ini() runs twice with a hotstart: a second open() on an open stream fails all later writes
    if(!mixlog.is_open())
    mixlog.open("./REEF3D_Log/REEF3D_sediment_mixture.dat");
    mixlog<<"# multi-fraction sediment bed, bed volume per fraction (incl. pores) and cumulative limiter volume"<<endl;
    mixlog<<"# sediter \t sedtime";
    for(int q=0;q<nf;++q)
    mixlog<<" \t V_"<<q+1;
    for(int q=0;q<nf;++q)
    mixlog<<" \t Vlim_"<<q+1;
    mixlog<<endl;
    }
}

// --------------------------------------------------------------------
// hiding/exposure factor for the critical Shields parameter
// Wu, Wang & Jia (2000), J. Hydraul. Res. 38(6)
// --------------------------------------------------------------------

double sediment_mixture::hiding(int ii, int jj, int q)
{
    if(hide_type==0)
    return 1.0;

    double ph=0.0;

    for(int r=0;r<nf;++r)
    ph += (*F[r])(ii,jj)*d[r]/(d[q]+d[r]);

    ph = MIN(MAX(ph,1.0e-6),1.0-1.0e-6);

    return pow((1.0-ph)/ph, -hide_m);
}

// --------------------------------------------------------------------
// bedload per fraction
// --------------------------------------------------------------------

void sediment_mixture::save_base(lexer *p, sediment_fdm *s)
{
    const int size = p->imax*p->jmax;

    for(int nn=0;nn<size;++nn)
    {
    shields_eff0.V[nn]   = s->shields_eff.V[nn];
    shields_crit0.V[nn]  = s->shields_crit.V[nn];
    tau_crit0.V[nn]      = s->tau_crit.V[nn];
    shearvel_crit0.V[nn] = s->shearvel_crit.V[nn];
    }
}

void sediment_mixture::restore_base(lexer *p, sediment_fdm *s)
{
    const int size = p->imax*p->jmax;

    for(int nn=0;nn<size;++nn)
    {
    s->shields_eff.V[nn]   = shields_eff0.V[nn];
    s->shields_crit.V[nn]  = shields_crit0.V[nn];
    s->tau_crit.V[nn]      = tau_crit0.V[nn];
    s->shearvel_crit.V[nn] = shearvel_crit0.V[nn];
    }

    s->dk = p->S20;
}

void sediment_mixture::load_fraction(lexer *p, ghostcell *pgc, sediment_fdm *s, int q)
{
    double xi;
    const double r = S20/d[q];


    s->dk = d[q];

    // bedshear computes shields_eff, tau_crit etc. with d = S20 and theta_cr = S30*reduce
    SLICELOOP4
    {
    xi = hiding(i,j,q);

    s->shields_eff(i,j)   = shields_eff0(i,j)*r;
    s->shields_crit(i,j)  = shields_crit0(i,j)*xi;
    s->tau_crit(i,j)      = tau_crit0(i,j)*xi/r;
    s->shearvel_crit(i,j) = shearvel_crit0(i,j)*sqrt(xi/r);
    }
}

void sediment_mixture::bedload_fractions(lexer *p, ghostcell *pgc, sediment_fdm *s, bedload *pbed)
{
    const int size = p->imax*p->jmax;

    save_base(p,s);

    for(int nn=0;nn<size;++nn)
    qbe_raw.V[nn] = 0.0;

    for(int q=0;q<nf;++q)
    {
        load_fraction(p,pgc,s,q);

        pbed->start(p,pgc,s);

        for(int nn=0;nn<size;++nn)
        {
        qbe_k[q]->V[nn] = F[q]->V[nn]*s->qbe.V[nn];
        qbe_raw.V[nn] += qbe_k[q]->V[nn];
        }
    }

    restore_base(p,s);

    for(int nn=0;nn<size;++nn)
    s->qbe.V[nn] = qbe_raw.V[nn];
}

// --------------------------------------------------------------------
// Exner
// bedload direction and relaxation zones act on the total qbe as
// cell-local factors; the same factor is applied to every fraction
// --------------------------------------------------------------------

void sediment_mixture::exner_begin(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    const int size = p->imax*p->jmax;

    for(int nn=0;nn<size;++nn)
    {
    qbe_tot.V[nn] = s->qbe.V[nn];
    fac.V[nn] = fabs(qbe_raw.V[nn])>1.0e-20?qbe_tot.V[nn]/qbe_raw.V[nn]:0.0;
    }

    save_base(p,s);
}

void sediment_mixture::exner_end(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    const int size = p->imax*p->jmax;

    restore_base(p,s);

    for(int nn=0;nn<size;++nn)
    s->qbe.V[nn] = qbe_tot.V[nn];
}

// --------------------------------------------------------------------
// active layer / substrate update for the bed change dhk[q] (m of bed,
// incl. pores) of each fraction
// --------------------------------------------------------------------

void sediment_mixture::bedchange(lexer *p, ghostcell *pgc, sediment_fdm *s, slice4 **dhk)
{
    double dz,dz0,Vs,area,corr;
    int neg;

    SEDSLICELOOP
    {
        area = p->DXN[IP]*p->DYN[JP];

        dz0=0.0;
        for(int q=0;q<nf;++q)
        {
        dhc[q] = (*dhk[q])(i,j);
        dz0 += dhc[q];
        }

        // limit the erosion of fractions that are not available:
        // for the set S of fractions that would become negative, V_q = 0 is solved exactly,
        //   dhc_q = dhk_q + x_q,  x_q = X FI_q - V_q,  X = -sum_S V_q/(1 - sum_S FI_q)
        // (V_q at dhk, FI fixed). S only grows; FI changes with the sign of dz, so the
        // solve is repeated until S and the sign are consistent (a few passes).
        // (The previous fixed-point update reduced the deficit only by FI_q per pass and
        // could stop after 20 passes with a remaining deficit, e.g. armoured beds.)
        for(int q=0;q<nf;++q)
        inS[q]=0;
        
        for(int it=0;it<2*nf+4;++it)
        {
            dz=0.0;
            for(int q=0;q<nf;++q)
            dz += dhc[q];

            for(int q=0;q<nf;++q)
            FI[q] = dz>=0.0?(*F[q])(i,j):(*Fs[q])(i,j);

            neg=0;
            for(int q=0;q<nf;++q)
            {
            V[q] = La*(*F[q])(i,j) + dhc[q] - dz*FI[q];

                if(V[q]<-1.0e-12*La && inS[q]==0)
                {
                inS[q]=1;
                neg=1;
                }
            }

            if(neg==0)
            break;
            
            // exact solve for S, starting from the unlimited bed change
            double sumV=0.0, sumFI=0.0, dzk=0.0;
            
            for(int q=0;q<nf;++q)
            dzk += (*dhk[q])(i,j);
            
            for(int q=0;q<nf;++q)
            FI[q] = dzk>=0.0?(*F[q])(i,j):(*Fs[q])(i,j);
            
            for(int q=0;q<nf;++q)
            if(inS[q]==1)
            {
            sumV  += La*(*F[q])(i,j) + (*dhk[q])(i,j) - dzk*FI[q];
            sumFI += FI[q];
            }
            
            if(sumFI>1.0-1.0e-12)
            break;
            
            double X = -sumV/(1.0-sumFI);
            
            for(int q=0;q<nf;++q)
            {
            dhc[q] = (*dhk[q])(i,j);
            
            if(inS[q]==1)
            dhc[q] += X*FI[q] - (La*(*F[q])(i,j) + (*dhk[q])(i,j) - dzk*FI[q]);
            }
        }
        
        // final volumes with the limited bed change
        dz=0.0;
        for(int q=0;q<nf;++q)
        dz += dhc[q];
        
        for(int q=0;q<nf;++q)
        {
        FI[q] = dz>=0.0?(*F[q])(i,j):(*Fs[q])(i,j);
        V[q] = La*(*F[q])(i,j) + dhc[q] - dz*FI[q];
        }

        // bookkeeping of the limited volume
        corr=0.0;
        for(int q=0;q<nf;++q)
        {
        corr += dhc[q] - (*dhk[q])(i,j);
        clip_k[q] += (dhc[q] - (*dhk[q])(i,j))*area;

            if(V[q]<0.0)
            {
            clip_k[q] += -V[q]*area;
            V[q] = 0.0;
            }
        }

        s->bedzh(i,j) += corr;

        dz = dz0 + corr;

        // substrate
        if(dz>0.0)
        {
            for(int q=0;q<nf;++q)
            (*Fs[q])(i,j) = (Hs(i,j)*(*Fs[q])(i,j) + dz*(*F[q])(i,j))/(Hs(i,j)+dz);

        Hs(i,j) += dz;
        }

        // degradation: the active layer takes dz*Fs from the substrate; beyond its thickness
        // there is nothing left, the missing volume is booked as limiter volume
        if(dz<=0.0)
        {
            if(Hs(i,j)+dz<0.0)
            for(int q=0;q<nf;++q)
            clip_k[q] += -(Hs(i,j)+dz)*(*Fs[q])(i,j)*area;
            
        Hs(i,j) = MAX(Hs(i,j)+dz,0.0);
        }

        // active layer
        Vs=0.0;
        for(int q=0;q<nf;++q)
        Vs += V[q];

        if(Vs>1.0e-30)
        for(int q=0;q<nf;++q)
        (*F[q])(i,j) = V[q]/Vs;
    }

    for(int q=0;q<nf;++q)
    {
    pgc->gcsl_start4(p,*F[q],1);
    pgc->gcsl_start4(p,*Fs[q],1);
    }

    pgc->gcsl_start4(p,Hs,1);
    pgc->gcsl_start4(p,s->bedzh,1);

    grain_stats(p,pgc);
}

// --------------------------------------------------------------------
// sandslide
// --------------------------------------------------------------------

void sediment_mixture::slide_zero(lexer *p, ghostcell *pgc)
{
    const int size = p->imax*p->jmax;

    for(int q=0;q<nf;++q)
    for(int nn=0;nn<size;++nn)
    fh_k[q]->V[nn] = 0.0;

    for(int q=0;q<9;++q)
    for(int nn=0;nn<size;++nn)
    out[q]->V[nn] = 0.0;
}

// called inside the slide loops: must not touch the static loop indices
void sediment_mixture::slide_transfer(int i0, int j0, int i1, int j1, double x)
{
    const int di=i1-i0;
    const int dj=j1-j0;

    if(abs(di)<=1 && abs(dj)<=1)
    {
    (*out[(di+1)*3+(dj+1)])(i0,j0) += x;
    return;
    }

    // not a direct neighbour: move with the active layer composition
    for(int q=0;q<nf;++q)
    {
    (*fh_k[q])(i0,j0) -= x*(*F[q])(i0,j0);
    (*fh_k[q])(i1,j1) += x*(*F[q])(i0,j0);
    }
}

// flux form sandslide (S90 5): upwind composition at each face
void sediment_mixture::slide_pde(lexer *p, sediment_fdm *s, slice &ci, int ii, int jj, double *fc)
{
    const int ni[4] = {1,-1,0,0};
    const int nj[4] = {0,0,1,-1};
    double flux;

    for(int f=0;f<4;++f)
    {
    int i1 = ii+ni[f];
    int j1 = jj+nj[f];

    flux = fc[f]*(s->bedzh(i1,j1)-s->bedzh(ii,jj))*0.5*(ci(i1,j1)+ci(ii,jj));

        for(int q=0;q<nf;++q)
        (*fh_k[q])(ii,jj) += flux*(flux>0.0?(*F[q])(i1,j1):(*F[q])(ii,jj));
    }
}

void sediment_mixture::composition_eroded(int ii, int jj, double L, double *C)
{
    if(L<=La)
    {
    for(int q=0;q<nf;++q)
    C[q] = (*F[q])(ii,jj);

    return;
    }

    for(int q=0;q<nf;++q)
    C[q] = (La*(*F[q])(ii,jj) + (L-La)*(*Fs[q])(ii,jj))/L;
}

// after the sandslide fill back: per-fraction bed change and sorting update
void sediment_mixture::slide_finish(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    double L,x;

    SEDSLICELOOP
    {
        L=0.0;
        for(int dir=0;dir<9;++dir)
        L += (*out[dir])(i,j);

        if(L>0.0)
        {
            composition_eroded(i,j,L,Fe);

            for(int dir=0;dir<9;++dir)
            {
            x = (*out[dir])(i,j);

                if(x>0.0)
                for(int q=0;q<nf;++q)
                {
                // receiver height with the area ratio (volume conserving, as the bed sand slide)
                (*fh_k[q])(i,j) -= x*Fe[q];
                (*fh_k[q])(i+dir/3-1,j+dir%3-1) += x*Fe[q]*(p->DXN[IP]*p->DYN[JP])/(p->DXN[IP+dir/3-1]*p->DYN[JP+dir%3-1]);
                }
            }
        }
    }

    for(int q=0;q<nf;++q)
    pgc->gcslparax_fh(p,*fh_k[q],4);

    bedchange(p,pgc,s,fh_k);
}

// --------------------------------------------------------------------
// grain size statistics of the active layer
// --------------------------------------------------------------------

void sediment_mixture::grain_stats(lexer *p, ghostcell *pgc)
{
    double cum,c0,c1,val;
    const double P[2] = {0.5,0.9};

    SLICELOOP4
    {
        dm(i,j)=0.0;
        for(int q=0;q<nf;++q)
        dm(i,j) += (*F[q])(i,j)*d[q];

        // cumulative distribution at the class centres, log-linear interpolation
        cum=0.0;
        for(int r=0;r<nf;++r)
        {
        FI[r] = cum + 0.5*(*F[order[r]])(i,j);
        cum += (*F[order[r]])(i,j);
        }

        for(int m=0;m<2;++m)
        {
            val = d[order[nf-1]];

            if(P[m]<=FI[0])
            val = d[order[0]];

            for(int r=0;r<nf-1;++r)
            if(P[m]>FI[r] && P[m]<=FI[r+1])
            {
            c0=FI[r];
            c1=FI[r+1];
            val = exp(log(d[order[r]]) + (P[m]-c0)/MAX(c1-c0,1.0e-20)*(log(d[order[r+1]])-log(d[order[r]])));
            }

            if(m==0)
            d50(i,j)=val;

            if(m==1)
            d90(i,j)=val;
        }
    }

    pgc->gcsl_start4(p,dm,1);
    pgc->gcsl_start4(p,d50,1);
    pgc->gcsl_start4(p,d90,1);
}

double sediment_mixture::ks_diameter(lexer *p, int ii, int jj)
{
    if(p->S56==1)
    return d50(ii,jj);

    if(p->S56==2)
    return dm(ii,jj);

    if(p->S56==3)
    return d90(ii,jj);

    return p->S20;
}

// --------------------------------------------------------------------
// log: bed volume per fraction and limiter volume
// --------------------------------------------------------------------

void sediment_mixture::print_log(lexer *p, ghostcell *pgc)
{
    double area;

    for(int q=0;q<nf;++q)
    {
    V[q]=0.0;
    clip_sum[q] = pgc->globalsum(clip_k[q]);
    }

    SEDSLICELOOP
    {
    area = p->DXN[IP]*p->DYN[JP];

        for(int q=0;q<nf;++q)
        V[q] += (La*(*F[q])(i,j) + Hs(i,j)*(*Fs[q])(i,j))*area;
    }

    for(int q=0;q<nf;++q)
    V[q] = pgc->globalsum(V[q]);

    if(p->mpirank==0)
    {
    mixlog<<p->sediter<<" \t "<<setprecision(8)<<p->sedtime;

    for(int q=0;q<nf;++q)
    mixlog<<" \t "<<setprecision(12)<<V[q];

    for(int q=0;q<nf;++q)
    mixlog<<" \t "<<setprecision(6)<<clip_sum[q];

    mixlog<<endl;
    }
}

int sediment_mixture::nfields()
{
    return nf+4;
}
