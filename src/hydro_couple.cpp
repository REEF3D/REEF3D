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
#include<algorithm>
#include<climits>
#include<cmath>
#include<iomanip>
#include<iostream>
#include<sys/stat.h>
#include<sys/types.h>

hydro_couple::hydro_couple(lexer *p, fdm2D *b, ghostcell *pgc) : hc(nullptr), t0(0.0), tnext(0.0),
    Vqo0(0.0), Vret(0.0), Vleave(0.0), Vinfout(0.0), Vinf0(0.0), S0(0.0), nexch(0)
{
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

    if(p->I45>0)
    stop(p,pgc,"the restart of a coupled run (I 45) is not available yet");

    hc = new hydro_core(p,pgc);
    hc->ini(p,pgc);
    hc->begin(p,pgc);

    const int N=hc->N;
    if(hc->poff.back()==0)
    stop(p,pgc,"hydro_grid.dat without the flow paths (version 2): run DIVEMesh W 10 again");

    t0=hc->t;
    tnext=p->simtime;
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
    hc->ret = (emax>=0) ? hc->down[emax] : -1;

    // the infiltration of SFLOW so far
    Fp.clear();
    SLICELOOP4
    Fp.push_back((1.0-b->hydro->sc(i,j))*b->hydro->F(i,j));

    Vqo0=b->hydro->Vqo;
    S0=hc->total_storage(pgc)+b->hydro->Vstore;

    if(p->mpirank==0)
    {
        cout<<"HYDRO coupling: "<<nov<<" of "<<N<<" hydrology cells under SFLOW ("<<atot*1.0e-6<<" km2), "<<nent
            <<" flow paths enter SFLOW, the outflow of SFLOW returns into cell "<<hc->ret
            <<(hc->ret<0 ? " (none: it leaves the catchment)" : "")<<"; exchange every "<<hc->dt<<" s, catchment time at the start "<<t0<<" s"<<endl;

        mkdir("./REEF3D_HYDRO",0777);
        fout.open("./REEF3D_HYDRO/REEF3D-HYDRO-Coupled.dat");
        fout<<"// REEF3D::HYDRO coupled budget [m3], cumulative since the start; storage = catchment model + SFLOW + returns of the last interval"<<endl;
        fout<<"// time_SFLOW[s] time_HYDRO[s] precipitation glacier_melt evaporation_HYDRO outflow_catchment rain_SFLOW sources_SFLOW inflow_SFLOW"
            <<" evaporation_SFLOW abstraction_SFLOW infiltration_outside outflow_SFLOW_leaving storage_HYDRO storage_SFLOW returns dS residual"
            <<" to_SFLOW received_SFLOW returned"<<endl;
        fout<<setprecision(14);
    }
}

// SFLOW -> HYDRO: the outflow and the infiltration of the last interval
void hydro_couple::exchange_in(lexer *p, fdm2D *b, ghostcell *pgc)
{
    const double dVo=b->hydro->Vqo-Vqo0;
    Vqo0=b->hydro->Vqo;

    if(hc->ret>=0)
    hc->Qret=dVo/hc->dt;
    else
    {
        hc->Qret=0.0;
        Vleave+=dVo;
    }

    fill(hc->addv.begin(),hc->addv.end(),0.0);
    double out=0.0;
    int q=0;
    SLICELOOP4
    {
        const double Fn=(1.0-b->hydro->sc(i,j))*b->hydro->F(i,j);
        const double dv=(Fn-Fp[q])*p->cellarea(i,j);
        Fp[q]=Fn;

        if(hcell[q]>=0)
        hc->addv[hcell[q]]+=dv;
        else
        out+=dv;

        ++q;
    }
    MPI_Allreduce(MPI_IN_PLACE,hc->addv.data(),hc->N,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);
    Vinfout+=pgc->globalsum(out);
}

// HYDRO -> SFLOW: the source rates of the interval
void hydro_couple::sources(lexer *p, fdm2D *b, ghostcell *pgc)
{
    slice4 &C=b->hydro->C;
    int q=0;
    double v=0.0;

    SLICELOOP4
    {
        const int h=hcell[q];
        C(i,j) = (h>=0 && aov[h]>0.0) ? hc->qsf[h]/aov[h] : 0.0;
        ++q;
    }

    for(const entry &e : ent)
    C(e.ci,e.cj) += hc->Qent[e.n]/e.area;

    SLICELOOP4
    v+=C(i,j)*p->cellarea(i,j);

    b->hydro->qcpl=pgc->globalsum(v);
}

void hydro_couple::budget(lexer *p, fdm2D *b, ghostcell *pgc)
{
    double s[3]={hc->Vp,hc->Vea,hc->Vgm};
    MPI_Allreduce(MPI_IN_PLACE,s,3,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);

    double pend=hc->Qret*hc->dt;
    for(double a : hc->addv)
    pend+=a;

    const double Sh=hc->total_storage(pgc);
    const hydro_sflow &H=*b->hydro;
    const double Stot=Sh+H.Vstore+pend;
    const double usr=H.Vsrc-H.Vcpl;
    const double res=Stot-S0-(s[0]+s[2]+H.Vrain+usr+H.Vqi-s[1]-H.Vevap-H.Vabs-Vinfout-hc->Vout-Vleave);

    if(p->mpirank==0)
    fout<<p->simtime<<" "<<hc->t<<" "<<s[0]<<" "<<s[2]<<" "<<s[1]<<" "<<hc->Vout<<" "<<H.Vrain<<" "<<usr<<" "<<H.Vqi<<" "
        <<H.Vevap<<" "<<H.Vabs<<" "<<Vinfout<<" "<<Vleave<<" "<<Sh<<" "<<H.Vstore<<" "<<pend<<" "<<Stot-S0<<" "<<res<<" "
        <<hc->Vto<<" "<<H.Vcpl<<" "<<hc->Vfrom<<endl;

    if(p->mpirank==0 && p->count%p->P12==0)
    cout<<"HYDRO coupled: t = "<<hc->t<<" s, to SFLOW "<<hc->Vto<<" m3, residual "<<res<<" m3"<<endl;
}

void hydro_couple::step(lexer *p, fdm2D *b, ghostcell *pgc)
{
    if(p->simtime>=tnext-1.0e-9*hc->dt)
    {
        if(nexch>0)
        exchange_in(p,b,pgc);

        budget(p,b,pgc);

        hc->step(p,pgc);
        sources(p,b,pgc);

        tnext+=hc->dt;
        ++nexch;
    }

    // the SFLOW steps end on the next exchange
    const double rem=tnext-p->simtime;
    if(rem>0.0 && p->dt>rem)
    p->dt=rem;
}

void hydro_couple::final(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // the run ends on an exchange: the last budget with the returns of the last interval
    if(fabs(p->simtime-(tnext))<=1.0e-9*hc->dt)
    {
        exchange_in(p,b,pgc);
        budget(p,b,pgc);
    }

    hc->end(p,pgc);
}
