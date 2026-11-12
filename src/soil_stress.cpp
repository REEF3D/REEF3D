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

#include"soil_stress.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<cmath>
#include<iostream>
#include<iomanip>
#include<sys/stat.h>

using namespace std;

double soil_stress::DXN(int ii) const {return p->DXN[ii+p->marge];}
double soil_stress::DYN(int jj) const {return p->DYN[jj+p->marge];}
double soil_stress::DZN(int kk) const {return p->DZN[kk+p->marge];}
double soil_stress::DXP(int ii) const {return p->DXP[ii+p->marge];}
double soil_stress::DYP(int jj) const {return p->DYP[jj+p->marge];}
double soil_stress::DZP(int kk) const {return p->DZP[kk+p->marge];}

soil_stress::soil_stress(lexer *pp, fdm *a, ghostcell *pgc) : p(pp)
{
    model = p->S200;
    knox = p->knox;
    knoy = p->knoy;
    knoz = p->knoz;
    dim3 = p->j_dir;
    NY = knoy+2;
    NZ = knoz+2;
    N = (knox+2)*NY*NZ;

    // material: cross-anisotropic skeleton, Jeng (1997)
    Ev = p->S201;
    const double mhh = p->S202_hh;
    const double mvh = p->S202_vh;
    const double om  = p->S203_O;
    const double den = (1.0+mhh)*(1.0 - mhh - 2.0*om*mvh*mvh);

    if(den<=0.0 || Ev<=0.0 || om<=0.0)
    {
        if(p->mpirank==0)
        cout<<"S 200: the elastic constants S 201 - S 203 are not admissible (Ev > 0, Eh/Ev > 0, (1+mhh)(1-mhh-2 Eh/Ev mvh^2) > 0)"<<endl;
        pgc->final();
        exit(1);
    }

    C11 = om*(1.0 - om*mvh*mvh)/den;
    C12 = om*(mhh + om*mvh*mvh)/den;
    C13 = om*mvh*(1.0+mhh)/den;
    C33 = (1.0 - mhh*mhh)/den;
    Gh  = om/(2.0*(1.0+mhh));
    Gv  = p->S203_L>0.0 ? p->S203_L : 0.5/(1.0+mvh);

    rhof = p->W1;
    grav = sqrt(p->W20*p->W20 + p->W21*p->W21 + p->W22*p->W22);
    const double gamw = rhof*grav;
    Kx = p->S204_x/gamw;
    Ky = p->S204_y/gamw;
    Kz = p->S204_z/gamw;
    rho = (1.0-p->S205)*p->S209 + p->S205*rhof;
    gsub = (1.0-p->S205)*(p->S209-rhof)*grav;
    K0 = p->S210;
    patm = 101325.0;

    // arrays
    ctype.assign(N,double(NOCELL));
    for(int c=0;c<3;++c)
    kind[c].assign(N,double(NODOF));
    X.assign(4*N,0.0);
    Xn.assign(4*N,0.0);
    Xnn.assign(4*N,0.0);
    pb.assign(N,0.0);
    pref.assign(N,0.0);
    stor.assign(N,0.0);
    D.assign(4*N,1.0);
    post.assign(size_t(npost)*N,0.0);

    t_last = dt = dt_old = 0.0;
    a0 = a1 = a2 = c0 = c1 = c2 = 0.0;
    nsteps = nsolve = 0;
    iter = iter_max = soilcells = 0;
    resid = soiltime = 0.0;
    printtime = 0.0;
    printcount = 0;

    if(p->mpirank==0)
    {
    cout<<"REEF3D::CFD soil stress model S 200 "<<model<<(model==2?" (dynamic u-p)":" (quasi-static u-p)")<<endl;
    cout<<"  Ev "<<Ev<<" Pa, C11 "<<C11<<" C12 "<<C12<<" C13 "<<C13<<" C33 "<<C33<<" G "<<Gh<<" L "<<Gv<<endl;
    cout<<"  k "<<p->S204_x<<" "<<p->S204_y<<" "<<p->S204_z<<" m/s, n "<<p->S205<<", Sr "<<p->S206<<", Kw "<<p->S207<<" Pa, start "<<p->S213<<" s"<<endl;
    }

    if(p->toporead==0 && p->S10==0 && p->mpirank==0)
    cout<<"S 200: no sediment bed (S 10 with S 57 or a DIVEMesh topography): no soil cells"<<endl;

    if(p->mpirank==0 && (p->S215>0 || p->S216>0))
    mkdir("./REEF3D_CFD_Soil",0777);

    MPI_Comm_split(pgc->mpi_comm,p->origin_i*(p->gknoy+1)+p->origin_j,p->origin_k,&colcomm);

    // start now (the initial CFD pressure the reference) or at the first step after S 213
    if(p->simtime>=p->S213-1.0e-12)
    ini(p,a,pgc);
}

void soil_stress::step(lexer *pp, fdm *a, ghostcell *pgc)
{
    p = pp;

    if(!active)
    {
        if(p->simtime>=p->S213-1.0e-12)
        ini(p,a,pgc);
        return;
    }

    ++nsteps;
    if(nsteps%p->S212!=0)
    return;

    dt = p->simtime - t_last;
    if(dt<=0.0)
    return;

    const double starttime = pgc->timer();

    classify(p,a,pgc);
    bed_pressure(p,a,pgc);

    // time discretisation: BDF2 with variable steps, BDF1 for the first step
    if(nsolve==0)
    {
        a0 = 1.0/dt;
        a1 = -1.0/dt;
        a2 = 0.0;
        Xnn = Xn;
        dt_old = dt;
    }
    else
    {
        const double r = dt/dt_old;
        a0 = (1.0+2.0*r)/(dt*(1.0+r));
        a1 = -(1.0+r)/dt;
        a2 = r*r/(dt*(1.0+r));
    }
    c0 = 2.0/(dt*(dt+dt_old));
    c1 = -2.0/(dt*dt_old);
    c2 = 2.0/(dt_old*(dt+dt_old));

    solve(p,pgc);

    Xnn = Xn;
    Xn = X;
    dt_old = dt;
    t_last = p->simtime;
    ++nsolve;

    soiltime = pgc->timer() - starttime;

    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<"soil: cells "<<soilcells<<"  iter "<<iter<<"  res "<<setprecision(3)<<resid<<"  time "<<setprecision(3)<<soiltime<<endl;

    // output
    bool printnow = false;
    if(p->P30>0.0 && p->simtime>printtime)
    {
        printnow = true;
        printtime += p->P30;
    }
    if(p->P30<=0.0 && p->P20>0 && p->count%p->P20==0)
    printnow = true;

    if(p->S216>0 || (p->S215>0 && printnow))
    postprocess(p);

    if(p->S216>0)
    print_probes(p,pgc);

    if(p->S215>0 && printnow)
    print_vtu(p);
}

void soil_stress::ini(lexer *p, fdm *a, ghostcell *pgc)
{
    // the reference of the dynamic bed pressure: the CFD pressure at the start time
    for(int i=0;i<knox;++i)
    for(int j=0;j<knoy;++j)
    for(int k=0;k<knoz;++k)
    pref[id(i,j,k)] = a->press(i,j,k);
    halo(pgc,pref.data());

    // storage n beta: Pw0 the absolute pore pressure, S 208 or local (atmosphere + CFD pressure)
    for(int i=-1;i<=knox;++i)
    for(int j=-1;j<=knoy;++j)
    for(int k=-1;k<=knoz;++k)
    {
        const int q = id(i,j,k);
        const double pw0 = p->S208>0.0 ? p->S208 : patm + max(pref[q],0.0);
        stor[q] = p->S205*(1.0/p->S207 + (1.0-p->S206)/pw0);
    }

    classify(p,a,pgc);
    bed_pressure(p,a,pgc);

    t_last = p->simtime;
    active = true;
    printtime = p->simtime + (p->S215>0 && p->P30>0.0 ? p->P30 : 0.0);

    if(p->mpirank==0)
    cout<<"soil: start at "<<p->simtime<<" s, soil cells "<<soilcells<<endl;

    if(p->S216>0 && p->mpirank==0)
    {
        probefile.open("./REEF3D_CFD_Soil/REEF3D-CFD-Soil-Probe.dat");
        probefile<<"REEF3D::CFD soil stress probes: excess pore pressure p (Pa), displacement ux uy uz (m), liquefaction potential LP"<<endl;
        for(int n=0;n<p->S216;++n)
        probefile<<"probe "<<n+1<<"  x "<<p->S216_x[n]<<"  y "<<p->S216_y[n]<<"  z "<<p->S216_z[n]<<endl;
        probefile<<endl<<"time";
        for(int n=0;n<p->S216;++n)
        probefile<<"\tp_"<<n+1<<"\tux_"<<n+1<<"\tuy_"<<n+1<<"\tuz_"<<n+1<<"\tLP_"<<n+1;
        probefile<<endl;
    }

    if(p->S215>0)
    {
        postprocess(p);
        print_vtu(p);
    }
}

void soil_stress::classify(lexer *p, fdm *a, ghostcell *pgc)
{
    const bool topo = p->toporead==1 || p->S10>0;
    const int lateral = p->S211==2 ? WALL : ROLLER;

    // the cells of the rank
    soilcells = 0;
    for(int i=-1;i<=knox;++i)
    for(int j=-1;j<=knoy;++j)
    for(int k=-1;k<=knoz;++k)
    {
        const bool inside = i>=0 && i<knox && j>=0 && j<knoy && k>=0 && k<knoz;
        int t = NOCELL;

        if(inside)
        {
            if(p->solidread==1 && a->solid(i,j,k)<0.0)
            t = WALL;
            else if(topo && a->topo(i,j,k)<0.0)
            {
            t = SOIL;
            ++soilcells;
            }
            else
            t = FREE;
        }
        else
        {
            // outside the domain: bottom fixed, lateral S 211, top free
            // (serial periodic directions, periodic 1, are filled by the halo)
            const bool xout = p->periodic1!=1 && ((i<0 && p->nb1==-2) || (i>=knox && p->nb4==-2));
            const bool yout = dim3 && p->periodic2!=1 && ((j<0 && p->nb3==-2) || (j>=knoy && p->nb2==-2));
            const bool zlow = p->periodic3!=1 && k<0 && p->nb5==-2;
            const bool zup  = p->periodic3!=1 && k>=knoz && p->nb6==-2;

            if(zlow)
            t = WALL;
            else if(xout || yout)
            t = lateral;
            else if(zup)
            t = FREE;
        }
        ctype[id(i,j,k)] = double(t);
    }
    halo(pgc,ctype.data());

    // 2D: the cells of the side rows are those of the row
    if(!dim3)
    for(int i=-1;i<=knox;++i)
    for(int k=-1;k<=knoz;++k)
    ctype[id(i,-1,k)] = ctype[id(i,1,k)] = ctype[id(i,0,k)];

    soilcells = int(pgc->globalsum(double(soilcells)));

    // dof kinds of the faces
    auto facekind = [&](int ta, int tb) -> int
    {
        if(ta==SOIL && tb==SOIL)
        return UNK;
        if((ta==SOIL && tb==FREE) || (ta==FREE && tb==SOIL))
        return UNK;
        if(ta==SOIL || tb==SOIL)
        return ZERO;
        return NODOF;
    };

    for(int c=0;c<3;++c)
    {
        std::fill(kind[c].begin(),kind[c].end(),double(NODOF));
        if(c==1 && !dim3)
        continue;

        for(int i=-1;i<=knox;++i)
        for(int j=-1;j<=knoy;++j)
        for(int k=-1;k<=knoz;++k)
        {
            int i2=i, j2=j, k2=k;
            if(c==0) ++i2;
            if(c==1) ++j2;
            if(c==2) ++k2;
            if(i2>knox || j2>knoy || k2>knoz)
            continue;
            kind[c][id(i,j,k)] = double(facekind(type(i,j,k),type(i2,j2,k2)));
        }
        halo(pgc,kind[c].data());
    }

    // owned unknowns
    for(int c=0;c<4;++c)
    own[c].clear();

    for(int i=0;i<knox;++i)
    for(int j=0;j<knoy;++j)
    for(int k=0;k<knoz;++k)
    {
        const int q = id(i,j,k);
        for(int c=0;c<3;++c)
        if(int(kind[c][q])==UNK)
        own[c].push_back(q);

        if(int(ctype[q])==SOIL)
        own[3].push_back(q);
    }

    // values without unknowns are zero (a moving bed: new soil starts from rest)
    for(int q=0;q<N;++q)
    {
        for(int c=0;c<3;++c)
        if(int(kind[c][q])!=UNK)
        X[c*N+q] = Xn[c*N+q] = Xnn[c*N+q] = 0.0;

        if(int(ctype[q])!=SOIL)
        X[3*N+q] = Xn[3*N+q] = Xnn[3*N+q] = 0.0;
    }
}

void soil_stress::bed_pressure(lexer *p, fdm *a, ghostcell *pgc)
{
    // the pore pressure at the bed: the dynamic pressure of the fluid cells, or S 217 p0 cos(kx - wt)
    const double t = p->simtime;
    const bool wave = p->S217>0 && p->S217_L>0.0 && p->S217_T>0.0;
    const double pi = 3.14159265358979323846;
    const double wk = wave ? 2.0*pi/p->S217_L : 0.0;
    const double ww = wave ? 2.0*pi/p->S217_T : 0.0;

    std::fill(pb.begin(),pb.end(),0.0);

    for(int i=0;i<knox;++i)
    for(int j=0;j<knoy;++j)
    for(int k=0;k<knoz;++k)
    if(type(i,j,k)==FREE)
    {
        const int q = id(i,j,k);
        if(wave)
        pb[q] = p->S217_p0*cos(wk*p->XP[i+p->marge] - ww*t);
        else
        pb[q] = a->press(i,j,k) - pref[q];
    }
    halo(pgc,pb.data());

    if(!dim3)
    for(int i=-1;i<=knox;++i)
    for(int k=-1;k<=knoz;++k)
    pb[id(i,-1,k)] = pb[id(i,1,k)] = pb[id(i,0,k)];
}

// one-layer halo with edges and corners: x, then y, then z, each with the halo of the directions before
void soil_stress::halo(ghostcell *pgc, double *f)
{
    std::vector<double> sbuf, rbuf;
    MPI_Comm comm = pgc->mpi_comm;

    auto exchange = [&](int dir)
    {
        int nbm, nbp, nm, per;
        if(dir==0) {nbm=p->nb1; nbp=p->nb4; nm=knox; per=p->periodic1;}
        else if(dir==1) {nbm=p->nb3; nbp=p->nb2; nm=knoy; per=p->periodic2;}
        else {nbm=p->nb5; nbp=p->nb6; nm=knoz; per=p->periodic3;}

        const int rm = nbm==-2 ? MPI_PROC_NULL : nbm;
        const int rp = nbp==-2 ? MPI_PROC_NULL : nbp;

        // the plane: the two other directions with their halo
        const int n1 = dir==0 ? knoy+2 : knox+2;
        const int n2 = dir==2 ? knoy+2 : knoz+2;
        sbuf.resize(n1*n2);
        rbuf.resize(n1*n2);

        auto index = [&](int l, int m1, int m2) -> int
        {
            if(dir==0) return id(l,m1-1,m2-1);
            if(dir==1) return id(m1-1,l,m2-1);
            return id(m1-1,m2-1,l);
        };

        // serial periodic direction (periodic 1, one rank across): the copy of the own layers
        if(per==1 && nbm==-2 && nbp==-2)
        {
            for(int m1=0;m1<n1;++m1)
            for(int m2=0;m2<n2;++m2)
            {
                f[index(nm,m1,m2)] = f[index(0,m1,m2)];
                f[index(-1,m1,m2)] = f[index(nm-1,m1,m2)];
            }
            return;
        }

        // to the minus side, from the plus side
        for(int m1=0;m1<n1;++m1)
        for(int m2=0;m2<n2;++m2)
        sbuf[m1*n2+m2] = f[index(0,m1,m2)];
        MPI_Sendrecv(sbuf.data(),n1*n2,MPI_DOUBLE,rm,710+dir,rbuf.data(),n1*n2,MPI_DOUBLE,rp,710+dir,comm,MPI_STATUS_IGNORE);
        if(rp!=MPI_PROC_NULL)
        for(int m1=0;m1<n1;++m1)
        for(int m2=0;m2<n2;++m2)
        f[index(nm,m1,m2)] = rbuf[m1*n2+m2];

        // to the plus side, from the minus side
        for(int m1=0;m1<n1;++m1)
        for(int m2=0;m2<n2;++m2)
        sbuf[m1*n2+m2] = f[index(nm-1,m1,m2)];
        MPI_Sendrecv(sbuf.data(),n1*n2,MPI_DOUBLE,rp,720+dir,rbuf.data(),n1*n2,MPI_DOUBLE,rm,720+dir,comm,MPI_STATUS_IGNORE);
        if(rm!=MPI_PROC_NULL)
        for(int m1=0;m1<n1;++m1)
        for(int m2=0;m2<n2;++m2)
        f[index(-1,m1,m2)] = rbuf[m1*n2+m2];
    };

    exchange(0);
    if(dim3)
    exchange(1);
    exchange(2);
}

void soil_stress::halo_all(ghostcell *pgc, std::vector<double> &x)
{
    for(int c=0;c<4;++c)
    if(c!=1 || dim3)
    halo(pgc,x.data()+size_t(c)*N);
}
