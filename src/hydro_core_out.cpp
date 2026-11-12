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

#include"hydro_core.h"
#include"lexer.h"
#include"ghostcell.h"
#include"regression_dump.h"
#include"reefckpt.h"
#include<algorithm>
#include<cmath>
#include<iomanip>
#include<iostream>
#include<sstream>
#include<sys/stat.h>
#include<sys/types.h>

void hydro_core::step(lexer *p, ghostcell *pgc)
{
    forcing_mean(p,t,t+dt);
    vertical(p);
    route(p,pgc);
    t+=dt;

    budget(p,pgc,false);

    if(p->R103>0.0 && t>=tvtk-1.0e-6*dt)
    {
        print_vtk(p,pgc);
        tvtk+=p->R103;
    }
}

void hydro_core::begin(lexer *p, ghostcell *pgc)
{
    budget(p,pgc,true);

    if(p->R103>0.0)
    {
        print_vtk(p,pgc);
        tvtk=t+p->R103;
    }
}

void hydro_core::end(lexer *p, ghostcell *pgc)
{
    budget(p,pgc,true);

    if(p->R107==1)
    write_state(p,pgc);
}

void hydro_core::run(lexer *p, ghostcell *pgc)
{
    if(p->mpirank==0)
    cout<<endl<<"REEF3D::HYDRO catchment model (R 100 1)"<<endl;

    ini(p,pgc);
    begin(p,pgc);

    regression_dump rdump(p);
    vector<string> rn;
    vector<vector<double>> rd;
    p->count=0;
    p->simtime=t;
    p->dt=dt;
    rdump.model_ini(p,"sum(SP^2) sum(SM^2) sum(SUZ^2) sum(SLZ^2) sum(Q^2)");
    regression_fields(rn,rd);
    rdump.model_state(p,rn,rd);

    const double t0=pgc->timer();
    int count=0;
    while(t<p->N41-1.0e-6*dt && count<p->N45)
    {
        step(p,pgc);
        ++count;

        p->count=count;
        p->simtime=t;
        rdump.model_step(p,regression_sums());

        if(p->mpirank==0 && p->P12>0 && count%p->P12==0)
        cout<<"HYDRO step "<<count<<"  t = "<<t/86400.0<<" day  Q outlet = "<<(gcell.empty() ? 0.0 : Q[gcell[0]])<<" m3/s"<<endl;
    }

    end(p,pgc);

    regression_fields(rn,rd);
    rdump.model_state(p,rn,rd);

    if(p->mpirank==0)
    cout<<"HYDRO: "<<count<<" steps, "<<pgc->timer()-t0<<" s"<<endl;
}

/*
REEF3D-HYDRO-Budget.dat, cumulative volumes [m3] since t = 0:
    precipitation (after PCORR and SFCF), glacier melt, actual evaporation, outflow,
    storage (snow, soil, upper and lower zone, MAXBAS buffer, channels), dS = storage - storage(0),
    residual = dS - (P + glacier melt - Ea - outflow)
REEF3D-HYDRO-Gauges.dat: discharge [m3/s] at the outlet (gauge 0) and the gauges
REEF3D-HYDRO-Means.dat: catchment means of the step: P, T, PET, Ea, runoff [mm/h] and states [mm]
*/
void hydro_core::budget(lexer *p, ghostcell *pgc, bool force)
{
    if((!force && p->R102>0.0 && t<tout-1.0e-6*dt) || t==tlast)
    return;

    tlast=t;

    if(p->R102>0.0)
    while(tout<=t+1.0e-6*dt)
    tout+=p->R102;

    // channel volume (identical on every rank)
    double vch=0.0;
    if(p->R106==1)
    for(int n=0; n<N; ++n)
    if(rescell[n]<0)
    vch+=len[n]*alpha[n]*pow(Q[n],0.6);

    // owned sums
    enum {SS,VP,VEA,VGM,AA,MP,MT,ME,MEA,MRUN,MSWE,MSM,MSUZ,MSLZ,NS};
    double s[NS]={0.0};
    s[SS]=owned_storage();
    s[VP]=Vp;
    s[VEA]=Vea;
    s[VGM]=Vgm;
    for(int n=n0; n<n1; ++n)
    {
        const int o=n-n0;
        const double A=area[n];
        double swe=0.0, sm=0.0;
        for(int b=0; b<nb; ++b)
        {
            swe += bf[n*nb+b]*(SP[o*nb+b]+WC[o*nb+b]);
            sm  += bf[n*nb+b]*(1.0-bg[n*nb+b])*SM[o*nb+b];
        }
        s[AA]+=A;
        s[MP]+=A*cP[o];
        s[MT]+=A*cT[o];
        s[ME]+=A*cE[o];
        s[MEA]+=A*cEa[o];
        s[MRUN]+=A*cRun[o];
        s[MSWE]+=A*swe;
        s[MSM]+=A*sm;
        s[MSUZ]+=A*SUZ[o];
        s[MSLZ]+=A*SLZ[o];
    }
    MPI_Allreduce(MPI_IN_PLACE,s,NS,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);

    const double store=s[SS]+vch+reservoir_storage();
    if(first)
    S0=store;
    first=false;

    if(p->mpirank!=0)
    return;

    // the files (continued on a restart from a checkpoint of a coupled run, I 45)
    if(!opened)
    {
        mkdir("./REEF3D_HYDRO",0777);
        const int I45 = (p->R100==2) ? p->I45 : 0;

        if(reefckpt_open(fb,"./REEF3D_HYDRO/REEF3D-HYDRO-Budget.dat",I45))
        {
        fb<<"// REEF3D::HYDRO volume budget [m3], cumulative since t = 0"<<endl;
        fb<<"// time[s] precipitation glacier_melt evaporation outflow storage channels dS residual relative_residual"
          <<(p->R100==2 ? " to_SFLOW from_SFLOW" : "")<<endl;
        }

        if(reefckpt_open(fg,"./REEF3D_HYDRO/REEF3D-HYDRO-Gauges.dat",I45))
        {
        fg<<"// REEF3D::HYDRO discharge [m3/s]"<<endl<<"// time[s]";
        for(size_t g=0; g<gid.size(); ++g)
        fg<<" Q_"<<gid[g];
        fg<<endl;
        }

        if(reefckpt_open(fm,"./REEF3D_HYDRO/REEF3D-HYDRO-Means.dat",I45))
        {
        fm<<"// REEF3D::HYDRO catchment means of the last step"<<endl;
        fm<<"// time[s] P[mm/h] T[C] PET[mm/h] Ea[mm/h] runoff[mm/h] SWE[mm] SM[mm] SUZ[mm] SLZ[mm]"<<endl;
        }

        fb<<setprecision(12);
        fg<<setprecision(12);
        fm<<setprecision(10);
        opened=true;
    }

    const double dS=store-S0;
    const double res=dS-(s[VP]+s[VGM]-s[VEA]-Vout-Vto+Vfrom);
    const double ref=max(1.0e-30,s[VP]+s[VGM]);
    fb<<t<<" "<<s[VP]<<" "<<s[VGM]<<" "<<s[VEA]<<" "<<Vout<<" "<<store<<" "<<vch<<" "<<dS<<" "<<res<<" "<<res/ref;
    if(p->R100==2)
    fb<<" "<<Vto<<" "<<Vfrom;
    fb<<endl;

    fg<<t;
    for(size_t g=0; g<gid.size(); ++g)
    fg<<" "<<Q[gcell[g]];
    fg<<endl;

    reservoir_output(p);

    const double A=s[AA];
    fm<<t<<" "<<s[MP]/A<<" "<<s[MT]/A<<" "<<s[ME]/A<<" "<<s[MEA]/A<<" "<<s[MRUN]/A<<" "<<s[MSWE]/A<<" "<<s[MSM]/A<<" "<<s[MSUZ]/A<<" "<<s[MSLZ]/A<<endl;
}

// cell fields on the hydrology raster (legacy VTK, cells outside the catchment: mask 0)
void hydro_core::print_vtk(lexer *p, ghostcell *pgc)
{
    enum {SWE,SMf,SUZf,SLZf,RUN,PR,TE,EA,NF};
    const char *name[NF]={"SWE_mm","SM_mm","SUZ_mm","SLZ_mm","runoff_mm_h","P_mm_h","T_C","Ea_mm_h"};

    vector<double> loc(size_t(n1-n0)*NF,0.0);
    for(int n=n0; n<n1; ++n)
    {
        const int o=n-n0;
        double *v=&loc[size_t(o)*NF];
        for(int b=0; b<nb; ++b)
        {
            v[SWE] += bf[n*nb+b]*(SP[o*nb+b]+WC[o*nb+b]);
            v[SMf] += bf[n*nb+b]*(1.0-bg[n*nb+b])*SM[o*nb+b];
        }
        v[SUZf]=SUZ[o];
        v[SLZf]=SLZ[o];
        v[RUN]=cRun[o];
        v[PR]=cP[o];
        v[TE]=cT[o];
        v[EA]=cEa[o];
    }

    vector<int> rc(cnt.size()), rd(cnt.size());
    for(size_t r=0; r<cnt.size(); ++r)
    {
        rc[r]=cnt[r]*NF;
        rd[r]=displ[r]*NF;
    }

    vector<double> all(p->mpirank==0 ? size_t(N)*NF : 1);
    MPI_Gatherv(loc.data(),int(loc.size()),MPI_DOUBLE,all.data(),rc.data(),rd.data(),MPI_DOUBLE,0,pgc->mpi_comm);

    if(p->mpirank==0)
    {
        vector<int> map(size_t(ni)*nj,-1);
        for(int n=0; n<N; ++n)
        map[I[n]+size_t(ni)*J[n]]=n;

        ostringstream fn;
        fn<<"./REEF3D_HYDRO/REEF3D-HYDRO-"<<setw(6)<<setfill('0')<<nvtk<<".vtk";
        ofstream f(fn.str());

        f<<"# vtk DataFile Version 3.0"<<endl<<"REEF3D::HYDRO t = "<<t<<" s"<<endl<<"ASCII"<<endl<<"DATASET STRUCTURED_POINTS"<<endl;
        f<<"DIMENSIONS "<<ni+1<<" "<<nj+1<<" 1"<<endl;
        f<<setprecision(12)<<"ORIGIN "<<x0<<" "<<y0<<" 0"<<endl<<"SPACING "<<D<<" "<<D<<" 1"<<endl;
        f<<"FIELD FieldData 1"<<endl<<"TIME 1 1 double"<<endl<<t<<endl;
        f<<"CELL_DATA "<<size_t(ni)*nj<<endl<<setprecision(8);

        f<<"SCALARS mask int 1"<<endl<<"LOOKUP_TABLE default"<<endl;
        for(size_t m=0; m<map.size(); ++m)
        f<<(map[m]>=0 ? 1 : 0)<<"\n";

        for(int k=0; k<NF; ++k)
        {
            f<<"SCALARS "<<name[k]<<" double 1"<<endl<<"LOOKUP_TABLE default"<<endl;
            for(size_t m=0; m<map.size(); ++m)
            f<<(map[m]>=0 ? all[size_t(map[m])*NF+k] : 0.0)<<"\n";
        }

        f<<"SCALARS Q_m3_s double 1"<<endl<<"LOOKUP_TABLE default"<<endl;
        for(size_t m=0; m<map.size(); ++m)
        f<<(map[m]>=0 ? Q[map[m]] : 0.0)<<"\n";
    }

    ++nvtk;
}

void hydro_core::regression_fields(vector<string> &names, vector<vector<double>> &data) const
{
    names={"SP","WC","SM","SUZ","SLZ","buf","qgen","Q"};
    data={SP,WC,SM,SUZ,SLZ,buf,vector<double>(qgen.begin()+n0,qgen.begin()+n1),vector<double>(Q.begin()+n0,Q.begin()+n1)};
}

vector<double> hydro_core::regression_sums() const
{
    vector<double> s(5,0.0);
    for(double v : SP)  s[0]+=v*v;
    for(double v : SM)  s[1]+=v*v;
    for(double v : SUZ) s[2]+=v*v;
    for(double v : SLZ) s[3]+=v*v;
    for(int n=n0; n<n1; ++n)
    s[4]+=Q[n]*Q[n];
    return s;
}

// water in the catchment model [m3]: snow, soil, zones, MAXBAS buffers, channels (collective)
double hydro_core::total_storage(ghostcell *pgc) const
{
    double s=owned_storage();
    MPI_Allreduce(MPI_IN_PLACE,&s,1,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);

    if(kw)
    for(int n=0; n<N; ++n)
    if(rescell[n]<0)
    s+=len[n]*alpha[n]*pow(Q[n],0.6);

    return s+reservoir_storage();
}

/*
the state of all cells, gathered on rank 0: per cell SP[nb] WC[nb] SM[nb] SUZ SLZ buf[nmb]
*/
static void cell_pack(const hydro_core *h, int nb, int nmb, int n0, int n1, const vector<double> &SP, const vector<double> &WC,
                      const vector<double> &SM, const vector<double> &SUZ, const vector<double> &SLZ, const vector<double> &buf,
                      vector<double> &loc)
{
    const int m=3*nb+2+nmb;
    loc.assign(size_t(n1-n0)*m,0.0);
    for(int n=n0; n<n1; ++n)
    {
        const int o=n-n0;
        double *v=&loc[size_t(o)*m];
        for(int b=0; b<nb; ++b)
        {
            v[b]=SP[o*nb+b];
            v[nb+b]=WC[o*nb+b];
            v[2*nb+b]=SM[o*nb+b];
        }
        v[3*nb]=SUZ[o];
        v[3*nb+1]=SLZ[o];
        for(int k=0; k<nmb; ++k)
        v[3*nb+2+k]=buf[size_t(o)*nmb+k];
    }
}

void hydro_core::gather_cells(lexer *p, ghostcell *pgc, vector<double> &all)
{
    const int m=3*nb+2+nmb;
    vector<double> loc;
    cell_pack(this,nb,nmb,n0,n1,SP,WC,SM,SUZ,SLZ,buf,loc);

    vector<int> rc(cnt.size()), rd(cnt.size());
    for(size_t r=0; r<cnt.size(); ++r)
    {
        rc[r]=cnt[r]*m;
        rd[r]=displ[r]*m;
    }

    all.assign(p->mpirank==0 ? size_t(N)*m : 1,0.0);
    MPI_Gatherv(loc.data(),int(loc.size()),MPI_DOUBLE,all.data(),rc.data(),rd.data(),MPI_DOUBLE,0,pgc->mpi_comm);
}

void hydro_core::scatter_cells(const double *all)
{
    const int m=3*nb+2+nmb;
    for(int n=n0; n<n1; ++n)
    {
        const int o=n-n0;
        const double *v=&all[size_t(n)*m];
        for(int b=0; b<nb; ++b)
        {
            SP[o*nb+b]=v[b];
            WC[o*nb+b]=v[nb+b];
            SM[o*nb+b]=v[2*nb+b];
        }
        SUZ[o]=v[3*nb];
        SLZ[o]=v[3*nb+1];
        for(int k=0; k<nmb; ++k)
        buf[size_t(o)*nmb+k]=v[3*nb+2+k];
    }
}

/*
hot start (R 107 / R 108), binary: "R3DHYST3", int N, nb, nmb, double t, the cells (gather_cells), Q[N], qgen[N],
int reservoirs, their levels ("R3DHYST2": without the reservoirs, "R3DHYST1": without qgen either)
*/
void hydro_core::write_state(lexer *p, ghostcell *pgc)
{
    vector<double> all;
    gather_cells(p,pgc,all);

    if(p->mpirank==0)
    {
        mkdir("./REEF3D_HYDRO",0777);
        ofstream f("./REEF3D_HYDRO/hydro_state.dat",ios::binary);
        const char magic[8]={'R','3','D','H','Y','S','T','3'};
        f.write(magic,8);
        int32_t h[3]={N,nb,nmb};
        f.write((char*)h,sizeof(h));
        f.write((char*)&t,sizeof(double));
        f.write((char*)all.data(),all.size()*sizeof(double));
        f.write((char*)Q.data(),Q.size()*sizeof(double));
        f.write((char*)qgen.data(),qgen.size()*sizeof(double));
        const int32_t nr=int32_t(resv.size());
        f.write((char*)&nr,sizeof(nr));
        for(const reservoir &r : resv)
        f.write((char*)&r.h,sizeof(double));
        cout<<"HYDRO: state at t = "<<t<<" s written to REEF3D_HYDRO/hydro_state.dat"<<endl;
    }
}

void hydro_core::read_state(lexer *p, ghostcell *pgc)
{
    ifstream f("hydro_state.dat",ios::binary);
    char magic[8]={0};
    int32_t h[3]={0,0,0};
    double ts=0.0;
    const int m=3*nb+2+nmb;
    vector<double> all(size_t(N)*m), q(N), g(N,0.0);

    f.read(magic,8);
    f.read((char*)h,sizeof(h));
    f.read((char*)&ts,sizeof(double));
    f.read((char*)all.data(),all.size()*sizeof(double));
    f.read((char*)q.data(),q.size()*sizeof(double));
    const string mg(magic,8);
    if(mg=="R3DHYST2" || mg=="R3DHYST3")
    f.read((char*)g.data(),g.size()*sizeof(double));

    // the levels of the reservoirs (a state without them: the levels of hydro_reservoir.dat)
    int32_t nr=0;
    vector<double> lev;
    if(mg=="R3DHYST3")
    {
        f.read((char*)&nr,sizeof(nr));
        lev.resize(max(nr,0));
        f.read((char*)lev.data(),lev.size()*sizeof(double));
    }

    const int ok = (f && (mg=="R3DHYST1" || mg=="R3DHYST2" || mg=="R3DHYST3") && h[0]==N && h[1]==nb && h[2]==nmb
                    && (mg!="R3DHYST3" || nr==int32_t(resv.size()))) ? 1 : 0;
    if(pgc->globalimin(ok)==0)
    stop(p,pgc,"R 108 1: hydro_state.dat is missing or does not fit hydro_grid.dat and MAXBAS (copy REEF3D_HYDRO/hydro_state.dat of an R 107 1 run)");

    t=ts;
    scatter_cells(all.data());
    Q=q;
    qgen=g;
    for(size_t k=0; k<lev.size(); ++k)
    resv[k].h=lev[k];

    if(p->mpirank==0)
    cout<<"HYDRO: hot start from hydro_state.dat at t = "<<t<<" s"<<endl;
}

/*
checkpoint of a coupled run (rank 0's blob): N nb nmb, t tout tvtk tlast nvtk S0 Vp Vea Vgm Vout Vto Vfrom,
the cells (gather_cells), Q, qgen, qret, addv
*/
void hydro_core::ckpt_state(lexer *p, ghostcell *pgc, vector<double> &v)
{
    double s[3]={Vp,Vea,Vgm};
    MPI_Allreduce(MPI_IN_PLACE,s,3,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);

    vector<double> all;
    gather_cells(p,pgc,all);

    v={double(N),double(nb),double(nmb),t,tout,tvtk,tlast,double(nvtk),S0,s[0],s[1],s[2],Vout,Vto,Vfrom};
    if(p->mpirank==0)
    v.insert(v.end(),all.begin(),all.end());
    else
    v.insert(v.end(),size_t(N)*(3*nb+2+nmb),0.0);
    for(const vector<double> *a : {&Q,&qgen,&qret,&addv})
    v.insert(v.end(),a->begin(),a->end());
    for(const reservoir &r : resv)
    {
        const double rv[4]={r.h,r.qin,r.qrel,r.qspill};
        v.insert(v.end(),rv,rv+4);
    }
}

bool hydro_core::ckpt_restore(lexer *p, const vector<double> &v)
{
    const int m=3*nb+2+nmb;
    if(v.size()!=15+size_t(N)*m+4*size_t(N)+4*resv.size() || int(v[0])!=N || int(v[1])!=nb || int(v[2])!=nmb)
    return false;

    t=v[3]; tout=v[4]; tvtk=v[5]; tlast=v[6]; nvtk=int(v[7]); S0=v[8];
    Vp = p->mpirank==0 ? v[9] : 0.0;
    Vea = p->mpirank==0 ? v[10] : 0.0;
    Vgm = p->mpirank==0 ? v[11] : 0.0;
    Vout=v[12]; Vto=v[13]; Vfrom=v[14];
    first=false;

    size_t k=15;
    scatter_cells(&v[k]);
    k+=size_t(N)*m;
    for(vector<double> *a : {&Q,&qgen,&qret,&addv})
    {
        a->assign(v.begin()+k,v.begin()+k+N);
        k+=N;
    }
    for(reservoir &r : resv)
    {
        r.h=v[k]; r.qin=v[k+1]; r.qrel=v[k+2]; r.qspill=v[k+3];
        k+=4;
    }
    return true;
}

void hydro_core::flush()
{
    fb.flush();
    fg.flush();
    fm.flush();
    if(fr.is_open())
    fr.flush();
}
