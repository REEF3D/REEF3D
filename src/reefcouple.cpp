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

#include"reefcouple.h"
#include"reef3d_mpmd.h"
#include"fnpf_couple_send.h"
#include"wave_lib_fnpf_couple.h"
#include"lexer.h"
#include"ghostcell.h"
#include<iostream>

MPI_Comm reefcouple::src_ic = MPI_COMM_NULL;
std::vector<MPI_Comm> reefcouple::recv_ic;
std::vector<int> reefcouple::recv_app;
std::vector<int> reefcouple::recv_model;
fnpf_couple_send *reefcouple::sender = nullptr;
wave_lib_fnpf_couple *reefcouple::recv = nullptr;

void reefcouple::roles(lexer *p, ghostcell *pgc)
{
    const int napps = reef3d_mpmd::napps;
    const int me = reef3d_mpmd::app;
    const bool want = (p->B92==62 && p->B90>0);
    
    // one model: nothing to couple
    if(napps==1)
    {
        if(want)
        {
            if(p->mpirank==0)
            cout<<endl<<"!!! B 92 62: the FNPF coupling needs an FNPF model in the same MPI job, e.g."<<endl
                <<"    mpirun -np 2 -wdir fnpf REEF3D : -np 8 -wdir cfd REEF3D !!!"<<endl<<endl;
            pgc->final(true);
        }
        return;
    }
    
    // model and role of every app (rank 0 of each app speaks for it)
    int mine[2] = {p->A10, want ? 1 : 0};
    std::vector<int> all(2*reef3d_mpmd::world_size);
    MPI_Allgather(mine,2,MPI_INT,all.data(),2,MPI_INT,MPI_COMM_WORLD);
    
    std::vector<int> model(napps), recv(napps);
    for(int a=0; a<napps; ++a)
    {
        model[a] = all[2*reef3d_mpmd::leader[a]];
        recv[a]  = all[2*reef3d_mpmd::leader[a]+1];
    }
    
    // the source: the lowest FNPF app that is no receiver
    int src=-1;
    for(int a=0; a<napps; ++a)
    if(model[a]==3 && recv[a]==0)
    {
    src=a;
    break;
    }
    
    int nrecv=0;
    for(int a=0; a<napps; ++a)
    nrecv += recv[a];
    
    if(reef3d_mpmd::world_rank==0)
    {
        cout<<endl<<"REEF3D coupled run: "<<napps<<" models"<<endl;
        for(int a=0; a<napps; ++a)
        cout<<"  app "<<a<<": A 10 "<<model[a]<<(a==src && nrecv>0 ? "  FNPF source" : "")<<(recv[a] ? "  receiver (B 92 62)" : "")<<endl;
        cout<<endl;
    }
    
    for(int a=0; a<napps; ++a)
    if(recv[a]==1 && model[a]!=6 && model[a]!=5)
    {
        if(reef3d_mpmd::world_rank==0)
        cout<<"!!! B 92 62: app "<<a<<" (A 10 "<<model[a]<<") cannot receive, only CFD (A 10 6) and NHFLOW (A 10 5) !!!"<<endl;
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    
    if(nrecv>0 && src<0)
    {
        if(reef3d_mpmd::world_rank==0)
        cout<<"!!! B 92 62: no FNPF model (A 10 3) in the job to receive from !!!"<<endl;
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    
    if(nrecv==0)
    return;
    
    // the rebalance stop (G 46) would end one model of the pair: off in coupled models
    if((me==src || recv[me]==1) && p->G46>0.0)
    {
        if(p->mpirank==0)
        cout<<"coupled run: G 46 (rebalance stop) is switched off"<<endl;
        p->G46 = 0.0;
    }
    
    // intercommunicators, source <-> each receiver in ascending order
    if(me==src)
    {
        for(int a=0; a<napps; ++a)
        if(recv[a]==1)
        {
            MPI_Comm ic;
            MPI_Intercomm_create(reef3d_comm(),0,MPI_COMM_WORLD,reef3d_mpmd::leader[a],4700+a,&ic);
            recv_ic.push_back(ic);
            recv_app.push_back(a);
            recv_model.push_back(model[a]);
        }
    }
    
    if(recv[me]==1)
    MPI_Intercomm_create(reef3d_comm(),0,MPI_COMM_WORLD,reef3d_mpmd::leader[src],4700+me,&src_ic);
}

void reefcouple::finish()
{
    if(sender!=nullptr)
    sender->finish();
    
    if(recv!=nullptr)
    recv->finish();
    
    sender = nullptr;
    recv = nullptr;
}
