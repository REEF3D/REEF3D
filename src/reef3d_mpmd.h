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

#ifndef REEF3D_MPMD_H_
#define REEF3D_MPMD_H_

#include<mpi.h>
#include<vector>

/*--------------------------------------------------------------------
Several REEF3D models in one MPI job (MPMD), for the runtime coupling
(reefcouple, B 92 62):

  mpirun -np 4  -wdir fnpf REEF3D : -np 16 -wdir cfd REEF3D

Each app context (MPI_APPNUM) is one model with its own working
directory (ctrl.txt, control.txt, grid, output) and its own
communicator, reef3d_comm(). A run with one app is unchanged
(reef3d_comm() is MPI_COMM_WORLD).

Command line (for launchers without -wdir or without MPI_APPNUM,
e.g. srun --multi-prog):
  -dir <path>   change to <path> before anything is read
  -app <n>      the app number of this process
--------------------------------------------------------------------*/

struct reef3d_mpmd
{
    static MPI_Comm model;          // communicator of this model
    static int app;                 // app number of this process
    static int napps;               // number of apps in the job
    static int world_rank, world_size;
    static std::vector<int> leader; // world rank of rank 0 of each app

    // after MPI_Init: -dir/-app, split of MPI_COMM_WORLD by app
    static void split(int &argc, char **argv);
};

// the communicator of this model: use it instead of MPI_COMM_WORLD
inline MPI_Comm reef3d_comm() {return reef3d_mpmd::model;}

#endif
