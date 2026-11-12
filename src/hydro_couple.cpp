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

#include"hydro_couple.h"
#include"hydro_core.h"
#include"hydro_sflow.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"reefckpt.h"
#include<algorithm>
#include<climits>
#include<cmath>
#include<iomanip>
#include<iostream>
#include<sys/stat.h>
#include<sys/types.h>

hydro_couple::hydro_couple(lexer *p, fdm2D *b, ghostcell *pgc) : hc(nullptr), t0(0.0), tnext(0.0), warm(0.0), qw(0.0), ret(-1),
    Vqo0(0.0), Vleave(0.0), Vinfout(0.0), S0(0.0), Vadj(0.0), Vwarm(0.0), nexch(0), started(false), opened(false)
{
    for(double &v : b0)
    v=0.0;
}

hydro_couple::~hydro_couple()
{
    delete hc;
}

void hydro_couple::stop(lexer *p, ghostcell *pgc, const string &msg) const
{
    if(p->mpirank==0)
    cout<<endl<<"!!! HYDRO (R 100 2): "<<msg<<" !!!"<<endl<<endl;

    pgc->final(true);
}

// world coordinates of the centre of an SFLOW cell
double hydro_couple::xw(lexer *p, int ii, int jj) const
{
    const double x = p->geo_curv==1 ? p->XC2D[(ii-p->imin)*p->jmax + (jj-p->jmin)] : p->XP[ii+p->marge];
    const double y = p->geo_curv==1 ? p->YC2D[(ii-p->imin)*p->jmax + (jj-p->jmin)] : p->YP[jj+p->marge];
    return p->Xout(x,y);
}

double hydro_couple::yw(lexer *p, int ii, int jj) const
{
    const double x = p->geo_curv==1 ? p->XC2D[(ii-p->imin)*p->jmax + (jj-p->jmin)] : p->XP[ii+p->marge];
    const double y = p->geo_curv==1 ? p->YC2D[(ii-p->imin)*p->jmax + (jj-p->jmin)] : p->YP[jj+p->marge];
    return p->Yout(x,y);
}

void hydro_couple::ini(lexer *p, fdm2D *b, ghostcell *pgc)
{
    if(p->mpirank==0)
    cout<<endl<<"REEF3D::HYDRO catchment model coupled with SFLOW (R 100 2)"<<endl;

    if(b->hydro==nullptr)
    stop(p,pgc,"no HYDRO object of SFLOW");

    hc = new hydro_core(p,pgc);
    hc->ini(p,pgc);

    const int N=hc->N;
    if(hc->poff.back()==0)
    stop(p,pgc,"hydro_grid.dat without the flow paths (version 2): run DIVEMesh W 10 again");

    t0=hc->t;
    warm = p->R109;
    if(warm>0.0 && p->R108!=1)
    {
        if(p->mpirank==0)
        cout<<"!!! HYDRO: R 109 (river start state) needs the hot start R 108 1, no start state"<<endl;
        warm=0.0;
    }
    tnext=p->simtime+warm;
    b->hydro->cpl=true;

    // overlap: the hydrology cell of every SFLOW cell
    vector<int> map(size_t(hc->ni)*hc->nj,-1);
    for(int n=0; n<N; ++n)
    map[hc->I[n]+size_t(hc->ni)*hc->J[n]]=n;

    aov.assign(N,0.0);
    hcell.clear();
    SLICELOOP4
    {
        const int II=int(floor((xw(p,i,j)-hc->x0)/hc->D)), JJ=int(floor((yw(p,i,j)-hc->y0)/hc->D));
        int h=-1;
        if(II>=0 && JJ>=0 && II<hc->ni && JJ<hc->nj)
        h=map[II+size_t(hc->ni)*JJ];

        hcell.push_back(h);
        if(h>=0)
        aov[h]+=p->cellarea(i,j);
    }
    MPI_Allreduce(MPI_IN_PLACE,aov.data(),N,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);

    int nov=0;
    double atot=0.0;
    for(int n=0; n<N; ++n)
    {
        hc->share[n] = aov[n]>0.0 ? min(1.0,aov[n]/hc->area[n]) : 0.0;
        if(aov[n]>0.0)
        ++nov;
        atot+=aov[n];
    }

    // entry points: the first point of the flow path of a cell in a fluid SFLOW cell (any rank)
    vector<int> kf(N,INT_MAX), own(N,INT_MAX);
    vector<entry> cand;
    for(int n=0; n<N; ++n)
    for(int k=hc->poff[n]; k<hc->poff[n+1]; ++k)
    {
        const double xm=p->Xin(hc->px[k],hc->py[k]), ym=p->Yin(hc->px[k],hc->py[k]);
        int ci, cj;
        if(b->hydro->locate(p,xm,ym,ci,cj))
        {
            kf[n]=k;
            cand.push_back({n,ci,cj,p->cellarea(ci,cj)});
            break;
        }
    }
    vector<int> kmin(kf);
    MPI_Allreduce(MPI_IN_PLACE,kmin.data(),N,MPI_INT,MPI_MIN,pgc->mpi_comm);
    for(int n=0; n<N; ++n)
    if(kf[n]==kmin[n] && kf[n]<INT_MAX)
    own[n]=p->mpirank;
    MPI_Allreduce(MPI_IN_PLACE,own.data(),N,MPI_INT,MPI_MIN,pgc->mpi_comm);

    ent.clear();
    for(const entry &e : cand)
    if(own[e.n]==p->mpirank)
    ent.push_back(e);

    int nent=0, emax=-1;
    for(int n=0; n<N; ++n)
    {
        hc->enter[n] = (kmin[n]<INT_MAX) ? 1 : 0;
        if(hc->enter[n])
        {
            ++nent;
            emax=n;
        }
    }

    // return: the first cell at or below the hydrology cell of an SFLOW cell whose flow path does not enter SFLOW
    rcell.clear();
    for(int h : hcell)
    {
        int r=h;
        while(r>=0 && hc->enter[r])
        r=hc->down[r];
        rcell.push_back(r);
    }
    ret = (emax>=0) ? hc->down[emax] : -1;

    // the infiltration of SFLOW so far
    SLICELOOP4
    b->hydro->Fp(i,j) = (1.0-b->hydro->sc(i,j))*b->hydro->F(i,j);

    // river start state (R 109): the steady discharges of the coupled network, SFLOW with them for warm seconds
    if(warm>0.0 && p->I45!=1)
    start_state(p,b,pgc);

    if(p->I45!=1)
    hc->begin(p,pgc);

    if(p->mpirank==0)
    {
        cout<<"HYDRO coupling: "<<nov<<" of "<<N<<" hydrology cells under SFLOW ("<<atot*1.0e-6<<" km2), "<<nent
            <<" flow paths enter SFLOW; the outflow of SFLOW returns into the cells below the entering paths"
            <<(p->A226==1 ? " (per boundary cell)" : " (explicit SFLOW: all into cell "+to_string(ret)+")")
            <<"; exchange every "<<hc->dt<<" s, catchment time at the start "<<t0<<" s"<<endl;
        if(warm>0.0)
        cout<<"HYDRO coupling: river start state, SFLOW runs "<<warm<<" s with the steady sources of the network ("
            <<b->hydro->qcpl<<" m3/s); the channels under SFLOW gave up "<<Vadj<<" m3"<<endl;
    }
}

// the steady discharges of the coupled network for the runoff of the catchment state, and the sources of SFLOW
void hydro_couple::start_state(lexer *p, fdm2D *b, ghostcell *pgc)
{
    const int N=hc->N;
    vector<double> Qin(N,0.0), Qs(N,0.0), qs(N,0.0), qe(N,0.0);

    for(int n=0; n<N; ++n)
    {
        qs[n]=hc->share[n]*hc->qgen[n];
        Qs[n]=Qin[n]+hc->qgen[n]-qs[n];

        if(hc->enter[n])
        qe[n]=Qs[n];
        else if(hc->down[n]>=0)
        Qin[hc->down[n]]+=Qs[n];
    }

    Vadj=0.0;
    if(p->R106==1)
    for(int n=0; n<N; ++n)
    if(hc->rescell[n]<0)
    Vadj+=hc->len[n]*hc->alpha[n]*(pow(hc->Q[n],0.6)-pow(Qs[n],0.6));

    hc->Q=Qs;
    sources(p,b,pgc,qs,qe);

    // the steady sources, raised over the first half of the start state (step)
    Cw.clear();
    SLICELOOP4
    Cw.push_back(b->hydro->C(i,j));
    qw=b->hydro->qcpl;
}

// the start of the coupled budget: the first exchange (after the river start state)
void hydro_couple::start_budget(lexer *p, fdm2D *b, ghostcell *pgc)
{
    hydro_sflow &H=*b->hydro;

    Vqo0=H.Vqo;
    SLICELOOP4
    {
    H.Bq(i,j)=0.0;
    H.Fp(i,j)=(1.0-H.sc(i,j))*H.F(i,j);
    }
    fill(hc->qret.begin(),hc->qret.end(),0.0);
    fill(hc->addv.begin(),hc->addv.end(),0.0);

    Vleave=0.0;
    Vinfout=0.0;
    Vwarm=H.Vcpl;
    const double v[6]={H.Vrain,H.Vsrc,H.Vcpl,H.Vqi,H.Vevap,H.Vabs};
    for(int k=0; k<6; ++k)
    b0[k]=v[k];

    S0=hc->total_storage(pgc)+H.Vstore;
    started=true;

    if(p->mpirank==0 && warm>0.0)
    cout<<"HYDRO coupling: the river start state took "<<Vwarm<<" m3 from the network, the coupled run starts"<<endl;
}

// SFLOW -> HYDRO: the outflow and the infiltration of the last interval
void hydro_couple::exchange_in(lexer *p, fdm2D *b, ghostcell *pgc)
{
    hydro_sflow &H=*b->hydro;
    const int N=hc->N;
    const double dVo=H.Vqo-Vqo0;
    Vqo0=H.Vqo;

    vector<double> &qret=hc->qret;
    fill(qret.begin(),qret.end(),0.0);
    double out=0.0;

    if(p->A226==1)
    {
        int q=0;
        SLICELOOP4
        {
            const double v=H.Bq(i,j);
            if(v!=0.0)
            {
                if(rcell[q]>=0)
                qret[rcell[q]]+=v;
                else
                out+=v;
            }
            H.Bq(i,j)=0.0;
            ++q;
        }
        MPI_Allreduce(MPI_IN_PLACE,qret.data(),N,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);
        Vleave+=pgc->globalsum(out);

        for(int n=0; n<N; ++n)
        qret[n]/=hc->dt;
    }
    else
    {
        if(ret>=0)
        qret[ret]=dVo/hc->dt;
        else
        Vleave+=dVo;
    }

    fill(hc->addv.begin(),hc->addv.end(),0.0);
    out=0.0;
    int q=0;
    SLICELOOP4
    {
        const double Fn=(1.0-H.sc(i,j))*H.F(i,j);
        const double dv=(Fn-H.Fp(i,j))*p->cellarea(i,j);
        H.Fp(i,j)=Fn;

        if(hcell[q]>=0)
        hc->addv[hcell[q]]+=dv;
        else
        out+=dv;

        ++q;
    }
    MPI_Allreduce(MPI_IN_PLACE,hc->addv.data(),N,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);
    Vinfout+=pgc->globalsum(out);
}

// HYDRO -> SFLOW: the source rates of the interval (qs runoff share, qe discharge of the entering cells) [m3/s]
void hydro_couple::sources(lexer *p, fdm2D *b, ghostcell *pgc, const vector<double> &qs, const vector<double> &qe)
{
    slice4 &C=b->hydro->C;
    int q=0;
    double v=0.0;

    SLICELOOP4
    {
        const int h=hcell[q];
        C(i,j) = (h>=0 && aov[h]>0.0) ? qs[h]/aov[h] : 0.0;
        ++q;
    }

    for(const entry &e : ent)
    C(e.ci,e.cj) += qe[e.n]/e.area;

    SLICELOOP4
    v+=C(i,j)*p->cellarea(i,j);

    b->hydro->qcpl=pgc->globalsum(v);
}

void hydro_couple::budget(lexer *p, fdm2D *b, ghostcell *pgc)
{
    double s[3]={hc->Vp,hc->Vea,hc->Vgm};
    MPI_Allreduce(MPI_IN_PLACE,s,3,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);

    double pend=0.0;
    for(int n=0; n<hc->N; ++n)
    pend+=hc->qret[n]*hc->dt+hc->addv[n];

    const double Sh=hc->total_storage(pgc);
    const hydro_sflow &H=*b->hydro;
    const double Stot=Sh+H.Vstore+pend;
    const double rain=H.Vrain-b0[0], rec=H.Vcpl-b0[2], usr=H.Vsrc-b0[1]-rec, qi=H.Vqi-b0[3], ev=H.Vevap-b0[4], ab=H.Vabs-b0[5];
    const double res=Stot-S0-(s[0]+s[2]+rain+usr+qi-s[1]-ev-ab-Vinfout-hc->Vout-Vleave);

    if(p->mpirank==0)
    {
        if(!opened)
        {
            mkdir("./REEF3D_HYDRO",0777);
            if(reefckpt_open(fout,"./REEF3D_HYDRO/REEF3D-HYDRO-Coupled.dat",p->I45))
            {
            fout<<"// REEF3D::HYDRO coupled budget [m3], cumulative since the start of the coupling; storage = catchment model + SFLOW + returns of the last interval"<<endl;
            fout<<"// time_SFLOW[s] time_HYDRO[s] precipitation glacier_melt evaporation_HYDRO outflow_catchment rain_SFLOW sources_SFLOW inflow_SFLOW"
                <<" evaporation_SFLOW abstraction_SFLOW infiltration_outside outflow_SFLOW_leaving storage_HYDRO storage_SFLOW returns dS residual"
                <<" to_SFLOW received_SFLOW returned"<<endl;
            }
            fout<<setprecision(14);
            opened=true;
        }

        fout<<p->simtime<<" "<<hc->t<<" "<<s[0]<<" "<<s[2]<<" "<<s[1]<<" "<<hc->Vout<<" "<<rain<<" "<<usr<<" "<<qi<<" "
            <<ev<<" "<<ab<<" "<<Vinfout<<" "<<Vleave<<" "<<Sh<<" "<<H.Vstore<<" "<<pend<<" "<<Stot-S0<<" "<<res<<" "
            <<hc->Vto<<" "<<rec<<" "<<hc->Vfrom<<endl;
    }

    if(p->mpirank==0 && p->count%p->P12==0)
    cout<<"HYDRO coupled: t = "<<hc->t<<" s, to SFLOW "<<hc->Vto<<" m3, residual "<<res<<" m3"<<endl;
}

void hydro_couple::step(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // river start state: the steady sources with a cosine ramp over the first half
    if(!started && !Cw.empty())
    {
        const double PI=3.14159265358979323846;
        const double f = (p->simtime<0.5*warm) ? 0.5*(1.0-cos(PI*p->simtime/(0.5*warm))) : 1.0;
        int q=0;
        SLICELOOP4
        {
        b->hydro->C(i,j)=f*Cw[q];
        ++q;
        }
        b->hydro->qcpl=f*qw;
    }

    if(p->simtime>=tnext-1.0e-9*hc->dt)
    {
        if(!started)
        start_budget(p,b,pgc);
        else
        exchange_in(p,b,pgc);

        budget(p,b,pgc);

        hc->step(p,pgc);
        sources(p,b,pgc,hc->qsf,hc->Qent);

        tnext+=hc->dt;
        ++nexch;
    }

    // the SFLOW steps end on the next exchange: the rest of the interval in equal steps no longer than the step
    // of the time step control (one short step at the end jolts the semi-implicit scheme)
    const double rem=tnext-p->simtime;
    if(rem>0.0 && p->dt>0.0)
    {
        const double n=ceil(rem/p->dt-1.0e-9);
        p->dt=rem/max(n,1.0);
    }
}

void hydro_couple::final(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // the run ends on an exchange: the last budget with the returns of the last interval
    if(started && fabs(p->simtime-tnext)<=1.0e-9*hc->dt)
    {
        exchange_in(p,b,pgc);
        budget(p,b,pgc);
    }

    hc->end(p,pgc);
}

/*
checkpoint (rank 0's blob): t0 tnext warm Vqo0 Vleave Vinfout S0 Vadj Vwarm b0[6] nexch started, then the
catchment model (hydro_core::ckpt_state)
*/
void hydro_couple::ckpt_state(lexer *p, ghostcell *pgc, vector<double> &v)
{
    vector<double> hv;
    hc->ckpt_state(p,pgc,hv);

    v={t0,tnext,warm,Vqo0,Vleave,Vinfout,S0,Vadj,Vwarm};
    v.insert(v.end(),b0,b0+6);
    v.push_back(double(nexch));
    v.push_back(started ? 1.0 : 0.0);
    v.insert(v.end(),hv.begin(),hv.end());
}

bool hydro_couple::ckpt_restore(lexer *p, fdm2D *b, const vector<double> &v)
{
    const size_t m=17;
    if(v.size()<m)
    return false;

    t0=v[0]; tnext=v[1]; warm=v[2]; Vqo0=v[3]; Vleave=v[4]; Vinfout=v[5]; S0=v[6]; Vadj=v[7]; Vwarm=v[8];
    for(int k=0; k<6; ++k)
    b0[k]=v[9+k];
    nexch=int(v[15]);
    started=(v[16]>0.5);

    return hc->ckpt_restore(p,vector<double>(v.begin()+m,v.end()));
}

void hydro_couple::flush()
{
    if(fout.is_open())
    fout.flush();
    hc->flush();
}
