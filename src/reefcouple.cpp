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
#include<string>

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
    
    // model, role and chosen source (I 246) of every app (rank 0 of each app speaks for it)
    int mine[3] = {p->A10, want ? 1 : 0, p->I246};
    std::vector<int> all(3*reef3d_mpmd::world_size);
    MPI_Allgather(mine,3,MPI_INT,all.data(),3,MPI_INT,MPI_COMM_WORLD);
    
    std::vector<int> model(napps), recv(napps), wsrc(napps), src(napps,-1);
    for(int a=0; a<napps; ++a)
    {
        model[a] = all[3*reef3d_mpmd::leader[a]];
        recv[a]  = all[3*reef3d_mpmd::leader[a]+1];
        wsrc[a]  = all[3*reef3d_mpmd::leader[a]+2];
    }
    
    auto fail = [&](const std::string &msg)
    {
        if(reef3d_mpmd::world_rank==0)
        cout<<"!!! "<<msg<<" !!!"<<endl;
        MPI_Abort(MPI_COMM_WORLD,1);
    };
    
    // default source: the lowest FNPF app that is no receiver
    int dsrc=-1;
    for(int a=0; a<napps; ++a)
    if(model[a]==3 && recv[a]==0)
    {
    dsrc=a;
    break;
    }
    
    int nrecv=0;
    for(int a=0; a<napps; ++a)
    if(recv[a]==1)
    {
        ++nrecv;
        
        if(model[a]!=6 && model[a]!=5 && model[a]!=3)
        fail("B 92 62: app "+std::to_string(a)+" (A 10 "+std::to_string(model[a])+") cannot receive, only CFD (A 10 6), NHFLOW (A 10 5) and FNPF (A 10 3)");
        
        src[a] = (wsrc[a]>=0) ? wsrc[a] : dsrc;
        
        if(src[a]<0)
        fail("B 92 62: no FNPF model (A 10 3) in the job to receive from (app "+std::to_string(a)+")");
        if(src[a]>=napps)
        fail("I 246 "+std::to_string(src[a])+" of app "+std::to_string(a)+": the job has the apps 0.."+std::to_string(napps-1));
        if(src[a]==a)
        fail("I 246 of app "+std::to_string(a)+": a model cannot receive from itself");
        if(model[src[a]]!=3)
        fail("I 246 of app "+std::to_string(a)+": app "+std::to_string(src[a])+" is no FNPF model (A 10 3)");
    }
    
    // no cycles in the chains (an FNPF receiver can send on)
    for(int a=0; a<napps; ++a)
    if(recv[a]==1)
    {
        int b=src[a], n=0;
        while(b>=0 && recv[b]==1 && n<=napps)
        {
            if(b==a)
            fail("B 92 62 / I 246: the coupling chain of app "+std::to_string(a)+" is a cycle");
            b=src[b];
            ++n;
        }
    }
    
    bool issrc=false;
    for(int a=0; a<napps; ++a)
    if(recv[a]==1 && src[a]==me)
    issrc=true;
    
    if(reef3d_mpmd::world_rank==0)
    {
        cout<<endl<<"REEF3D coupled run: "<<napps<<" models"<<endl;
        for(int a=0; a<napps; ++a)
        {
            bool s=false;
            for(int r=0; r<napps; ++r)
            if(recv[r]==1 && src[r]==a)
            s=true;
            
            cout<<"  app "<<a<<": A 10 "<<model[a]<<(s ? "  FNPF source" : "");
            if(recv[a]==1)
            cout<<"  receiver (B 92 62) of app "<<src[a];
            cout<<endl;
        }
        cout<<endl;
    }
    
    // feedback keys: I 244 (overlap box) only for CFD, I 243/I 244 not for an FNPF receiver
    if(want && p->A10==3 && (p->I243>0 || p->I244>0))
    {
        if(p->mpirank==0)
        cout<<"!!! B 92 62: I 243 / I 244 (feedback) are not available for an FNPF receiver !!!"<<endl;
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    
    if(want && p->A10==5 && p->I244>0)
    {
        if(p->mpirank==0)
        cout<<"!!! B 92 62: I 244 (overlap box) is available for CFD receivers only !!!"<<endl;
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    
    if(nrecv==0)
    return;
    
    // the rebalance stop (G 46) would end one model of the pair: off in coupled models
    if((issrc || recv[me]==1) && p->G46>0.0)
    {
        if(p->mpirank==0)
        cout<<"coupled run: G 46 (rebalance stop) is switched off"<<endl;
        p->G46 = 0.0;
    }
    
    // intercommunicators source <-> receiver, all pairs in the order of the receivers: every app
    // takes its pairs in the same global order (a model in a chain is receiver and source)
    for(int a=0; a<napps; ++a)
    if(recv[a]==1)
    {
        if(me==src[a])
        {
            MPI_Comm ic;
            MPI_Intercomm_create(reef3d_comm(),0,MPI_COMM_WORLD,reef3d_mpmd::leader[a],4700+a,&ic);
            recv_ic.push_back(ic);
            recv_app.push_back(a);
            recv_model.push_back(model[a]);
        }
        
        if(me==a)
        MPI_Intercomm_create(reef3d_comm(),0,MPI_COMM_WORLD,reef3d_mpmd::leader[src[a]],4700+a,&src_ic);
    }
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
