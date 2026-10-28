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


#include"seastate_remote.h"
#include"reef3d_mpmd.h"
#include"seastate_f.h"
#include"seastate_vtp.h"
#include"fdm_seastate.h"
#include"lexer.h"
#include"ghostcell.h"
#include<algorithm>
#include<cmath>
#include<iostream>
#include<vector>

MPI_Comm seastate_remote::ic = MPI_COMM_NULL;
int seastate_remote::role = 0;
seastate_coupling *seastate_remote::link = nullptr;
lexer *seastate_remote::plex = nullptr;

void seastate_remote::roles(lexer *p, ghostcell *pgc)
{
    plex = p;

    const int napps = reef3d_mpmd::napps;
    const int me = reef3d_mpmd::app;
    const int mine_role = (p->A750==2) ? (p->A10==7 ? 2 : ((p->A10==2 || p->A10==5) ? 1 : 3)) : 0;

    // one model: a wave app alone has no host (a host alone stops in remote_ini with the launch line)
    if(napps==1)
    {
        if(mine_role==2 || mine_role==3)
        {
            if(p->mpirank==0)
            cout<<endl<<"!!! A 750 2: "<<(mine_role==2 ? "the remote waves (A 10 7) need a host model, SFLOW (A 10 2) or NHFLOW (A 10 5), in the same MPI job, e.g."
                                                   : "the remote waves are coupled to SFLOW (A 10 2) or NHFLOW (A 10 5), e.g.")<<endl
                <<"    mpirun -np 4 -wdir surge REEF3D : -np 2 -wdir waves REEF3D !!!"<<endl<<endl;
            pgc->final(true);
        }
        return;
    }

    // the role of every app (rank 0 of each app speaks for it)
    std::vector<int> all(reef3d_mpmd::world_size);
    MPI_Allgather(&mine_role,1,MPI_INT,all.data(),1,MPI_INT,MPI_COMM_WORLD);

    int host=-1, wave=-1, nhost=0, nwave=0, nbad=0;
    for(int a=0; a<napps; ++a)
    {
    const int r = all[reef3d_mpmd::leader[a]];
        if(r==1)
        {
        host = a;
        ++nhost;
        }
        if(r==2)
        {
        wave = a;
        ++nwave;
        }
        if(r==3)
        ++nbad;
    }

    if(nhost+nwave+nbad==0)
    return;

    if(nhost!=1 || nwave!=1 || nbad>0)
    {
        if(reef3d_mpmd::world_rank==0)
        cout<<endl<<"!!! A 750 2: one host model (SFLOW A 10 2 or NHFLOW A 10 5, A 750 2) and one wave model (SEASTATE A 10 7,"
            <<" A 750 2) per job; found "<<nhost<<" host(s), "<<nwave<<" wave model(s)"<<(nbad>0 ? " and A 750 2 in other models" : "")<<" !!!"<<endl<<endl;
        MPI_Abort(MPI_COMM_WORLD,1);
    }

    if(reef3d_mpmd::world_rank==0)
    cout<<endl<<"REEF3D coupled run: SEASTATE remote waves (A 750 2), host app "<<host<<", wave app "<<wave<<endl;

    if(me!=host && me!=wave)
    return;

    role = (me==host) ? 1 : 2;

    // the rebalance stop (G 46) would end one model of the pair
    if(p->G46>0.0)
    {
        if(p->mpirank==0)
        cout<<"coupled run: G 46 (rebalance stop) is switched off"<<endl;
        p->G46 = 0.0;
    }

    MPI_Intercomm_create(reef3d_comm(),0,MPI_COMM_WORLD,reef3d_mpmd::leader[me==host ? wave : host],tag_ic,&ic);
}

void seastate_remote::finish()
{
    if(link!=nullptr && plex!=nullptr)
    link->remote_stop(plex);

    link = nullptr;
}

// ------------------------------------------------------------------ wave side

seastate_remote_wave::seastate_remote_wave(lexer *p) : seastate_coupling(p), etaW(p), WLW(p), UW(p), VW(p), wetd(p)
{
}

void seastate_remote_wave::run(lexer *p, ghostcell *pgc)
{
    const double t0 = pgc->timer();
    remote = true;

    if(ic_null())
    {
        if(p->mpirank==0)
        cout<<endl<<"!!! A 750 2: no host model in the job !!!"<<endl<<endl;
        pgc->final(true);
    }

    // set-up: the keys both sides must share and the grids (mirror of seastate_coupling::remote_ini)
    std::vector<double> xs, ys;
    centres(p,pgc,xs,ys);

    const double mine[9] = {double(xs.size()), double(ys.size()), p->A706, double(p->A751), double(p->A753),
                            0.0, p->W203, p->A780, double(p->A10)};
    double theirs[9] = {0,0,0,0,0,0,0,0,0};

    if(p->mpirank==0)
    {
    MPI_Sendrecv(mine,9,MPI_DOUBLE,0,4952,theirs,9,MPI_DOUBLE,0,4951,seastate_remote::ic,MPI_STATUS_IGNORE);

    std::vector<double> g(xs);
    g.insert(g.end(),ys.begin(),ys.end());
    std::vector<double> h(size_t(theirs[0])+size_t(theirs[1]));
    MPI_Sendrecv(g.data(),int(g.size()),MPI_DOUBLE,0,4954,h.data(),int(h.size()),MPI_DOUBLE,0,4953,seastate_remote::ic,MPI_STATUS_IGNORE);
    rx.assign(h.begin(),h.begin()+size_t(theirs[0]));
    ry.assign(h.begin()+size_t(theirs[0]),h.end());
    }

    MPI_Bcast(theirs,9,MPI_DOUBLE,0,pgc->mpi_comm);
    rx.resize(size_t(theirs[0]));
    ry.resize(size_t(theirs[1]));
    MPI_Bcast(rx.data(),int(rx.size()),MPI_DOUBLE,0,pgc->mpi_comm);
    MPI_Bcast(ry.data(),int(ry.size()),MPI_DOUBLE,0,pgc->mpi_comm);

    // the host aborts on a mismatch of A 706, A 751, A 753 (it reports both sides)
    if(theirs[2]!=p->A706 || int(theirs[3])!=p->A751 || int(theirs[4])!=p->A753)
    MPI_Abort(MPI_COMM_WORLD,1);

    host_wind = theirs[5]>0.5;

    ini_wave(p,pgc,int(theirs[8])==5 ? "NHFLOW" : "SFLOW");

    fdm_seastate *e = pwave->e;

    // the first wave fields: own environment, no current
    IMALOOP
    JMALOOP
    {
    etaW(i,j) = 0.0;
    WLW(i,j)  = e->depth(i,j);
    UW(i,j) = VW(i,j) = 0.0;
    }

    forces(p,pgc,UW,VW);

    if(p->P10>0)
    pwave->printer()->print2D(p,e,pgc);

    const size_t Nh = rx.size()*ry.size();
    rbuf.assign(2+4*Nh,0.0);

    if(p->mpirank==0)
    cout<<"SEASTATE remote waves (A 750 2): host grid "<<rx.size()<<" x "<<ry.size()<<" cells, this grid "<<xs.size()<<" x "<<ys.size()
        <<"; wave steps of "<<p->A706<<" s between the exchanges, until the host stops"<<endl<<"starting mainloop.SEASTATE (remote)"<<endl;

    const double t_ini = pgc->timer() - t0;
    double t_step = 0.0;

    for(;;)
    {
        IMALOOP
        JMALOOP
        wetd(i,j) = e->wet(i,j);

        // send the wave fields of this wave time, receive the host fields of the host's exchange time (or the stop)
        gather(p,pgc,{&Fw,&Gw,&divM,&Hsw,&Tpw,&wetd},sbuf,1.0);
        exchange(p,pgc,sbuf,rbuf,4962,4961);

        if(rbuf[0]<0.5)
        break;

        const double ts = pgc->timer();

        map_host(p);
        environment(p,etaW,WLW,UW,VW);

        // the next exchange is A 706 s of host time later: the wave state of that time
        p->simtime = rbuf[1] + p->A706;
        ++p->count;

        wave_step(p,pgc,p->simtime - tlast);
        forces(p,pgc,UW,VW);
        print(p,pgc);

        t_step += pgc->timer() - ts;
    }

    if(p->mpirank==0)
    {
    cout<<endl<<"******************************"<<endl<<endl;
    cout<<"SEASTATE remote waves: "<<nstep<<" wave steps to simtime "<<p->simtime<<" s; set-up "<<t_ini<<" s, wave steps "<<t_step
        <<" s, waiting for the host "<<t_wait<<" s"<<endl<<endl;
    }

    pgc->final();
}

bool seastate_remote_wave::ic_null() const
{
    return seastate_remote::ic==MPI_COMM_NULL;
}

void seastate_remote_wave::map_host(lexer *p)
{
    fdm_seastate *e = pwave->e;
    const int nxh = int(rx.size()), nyh = int(ry.size());
    const size_t Nh = size_t(nxh)*size_t(nyh);
    const double *Eh = &rbuf[2], *Wh = &rbuf[2+Nh], *Uh = &rbuf[2+2*Nh], *Vh = &rbuf[2+3*Nh];

    IMALOOP
    JMALOOP
    {
    const double xc = p->XP[IP], yc = p->YP[JP];

    // outside the host grid (or no wave cell): the own still water
    etaW(i,j) = 0.0;
    WLW(i,j)  = std::max(p->wd - e->bed(i,j),0.0);
    UW(i,j) = VW(i,j) = 0.0;

        if(e->wet0(i,j)==0 || !covers(rx,ry,xc,yc))
        continue;

    // the host cells inside this wave cell
    const int a0 = int(std::lower_bound(rx.begin(),rx.end(),p->XN[IP])  - rx.begin());
    const int a1 = int(std::lower_bound(rx.begin(),rx.end(),p->XN[IP1]) - rx.begin());
    int b0 = 0, b1 = 1;

        if(nyh>1)
        {
        b0 = int(std::lower_bound(ry.begin(),ry.end(),p->YN[JP])  - ry.begin());
        b1 = int(std::lower_bound(ry.begin(),ry.end(),p->YN[JP1]) - ry.begin());
        }

    double se=0.0, su=0.0, sv=0.0, sw=0.0, nw=0.0, n=0.0;

        if(a1>a0 && b1>b0)
        {
            // mean over the host cells of the box
            for(int a=a0; a<a1; ++a)
            for(int b=b0; b<b1; ++b)
            {
            const size_t q = size_t(a)*nyh + b;

            sw += Wh[q];
            n += 1.0;

                if(Wh[q]>p->A705)
                {
                se += Eh[q];
                su += Uh[q];
                sv += Vh[q];
                nw += 1.0;
                }
            }
        }
        else
        {
            // host grid coarser than the wave grid: bilinear at the cell centre
            int a, b;
            double wx, wy;
            bilinear(rx,ry,xc,yc,a,b,wx,wy);

            for(int da=0; da<2; ++da)
            for(int db=0; db<2; ++db)
            {
            const double w = (da==0 ? 1.0-wx : wx)*(db==0 ? 1.0-wy : wy);
            const size_t q = size_t(std::min(a+da,nxh-1))*nyh + std::min(b+db,nyh-1);

                if(w<=0.0)
                continue;

            sw += w*Wh[q];
            n += w;

                if(Wh[q]>p->A705)
                {
                se += w*Eh[q];
                su += w*Uh[q];
                sv += w*Vh[q];
                nw += w;
                }
            }
        }

        if(!(n>0.0))
        continue;

    WLW(i,j) = sw/n;

        if(nw>0.0)
        {
        etaW(i,j) = se/nw;
        UW(i,j) = su/nw;
        VW(i,j) = sv/nw;
        }
    }
}
