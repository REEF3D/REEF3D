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

#include"reef3d_mpmd.h"
#include<cstring>
#include<cstdlib>
#include<iostream>
#include<unistd.h>

MPI_Comm reef3d_mpmd::model = MPI_COMM_WORLD;
int reef3d_mpmd::app = 0;
int reef3d_mpmd::napps = 1;
int reef3d_mpmd::world_rank = 0;
int reef3d_mpmd::world_size = 1;
std::vector<int> reef3d_mpmd::leader(1,0);

void reef3d_mpmd::split(int &argc, char **argv)
{
    MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&world_size);

    // app number from the launcher (MPI_APPNUM: the app context of mpirun a : b)
    int *val=nullptr, flag=0;
    MPI_Comm_get_attr(MPI_COMM_WORLD,MPI_APPNUM,&val,&flag);
    app = (flag && val!=nullptr) ? *val : 0;

    // command line: -dir <path>, -app <n>
    for(int n=1; n<argc; ++n)
    {
        if(strcmp(argv[n],"-app")==0 && n+1<argc)
        app = atoi(argv[++n]);

        else
        if(strcmp(argv[n],"-dir")==0 && n+1<argc)
        {
            if(chdir(argv[++n])!=0)
            {
            std::cout<<"!!! REEF3D: cannot change to the directory "<<argv[n]<<" (world rank "<<world_rank<<") !!!"<<std::endl;
            MPI_Abort(MPI_COMM_WORLD,1);
            }
        }
    }

    if(app<0)
    app=0;

    int amax=0;
    MPI_Allreduce(&app,&amax,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    napps = amax+1;

    leader.assign(napps,-1);

    if(napps==1)
    {
    model = MPI_COMM_WORLD;
    leader[0] = 0;
    return;
    }

    MPI_Comm_split(MPI_COMM_WORLD,app,world_rank,&model);

    std::vector<int> all(world_size);
    MPI_Allgather(&app,1,MPI_INT,all.data(),1,MPI_INT,MPI_COMM_WORLD);

    for(int r=world_size-1; r>=0; --r)
    leader[all[r]] = r;

    for(int a=0; a<napps; ++a)
    if(leader[a]<0)
    {
        if(world_rank==0)
        std::cout<<"!!! REEF3D MPMD: app "<<a<<" has no processes (apps 0.."<<napps-1<<" expected) !!!"<<std::endl;
        MPI_Abort(MPI_COMM_WORLD,1);
    }
}
