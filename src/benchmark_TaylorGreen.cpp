/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

#include"benchmark_TaylorGreen.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<sys/stat.h>

benchmark_TaylorGreen::benchmark_TaylorGreen(lexer *p, fdm *a) : gradient(p)
{
    double L = 1.0;
    double U = 1600.0;
    
    double x,y,z;
    
    ULOOP
    {
        x = p->pos1_x();
        y = p->pos1_y();
        z = p->pos1_z();
        a->u(i,j,k) = U*(sin(x/L)*cos(y/L)*cos(z/L));
    }   
    
    VLOOP
    {
        x = p->pos2_x();
        y = p->pos2_y();
        z = p->pos2_z();
        a->v(i,j,k) = U*(cos(x/L)*sin(y/L)*cos(z/L));
    }
    
    WLOOP
    a->w(i,j,k) = 0.0;
    
    if(p->mpirank==0)
    {
        mkdir("./REEF3D_CFD_TaylorGreen_Diss", 0777);
	
        ofstream print;
        print.open("./REEF3D_CFD_TaylorGreen_Diss/REEF3D_CFD_TG_Diss.dat");
        print<<"time \t epsilon "<<endl;
        print.close();
    }
}

benchmark_TaylorGreen::~benchmark_TaylorGreen()
{
}

void benchmark_TaylorGreen::start(lexer* p, fdm *a, ghostcell *pgc, convection *pconvec )
{
    double vx, vy, vz, vmag2, vol;

    double vVolAvg = 0.0;
    double volTot = 0.0;

    // Local calculation of volume averaged vorticity
    LOOP
    {
        if(p->j_dir==1)
        {
            vx = pvdz(p,a) - pwdy(p,a); 
            vy = pudz(p,a) - pwdx(p,a);
            vz = pvdx(p,a) - pudy(p,a); 
        }
        else
        {
            vx = 0.0;
            vy = pudz(p,a) - pwdx(p,a);
            vz = 0.0;
        }
       
        vmag2 = vx*vx + vy*vy + vz*vz;
        
        vol = p->DXN[IP]*p->DYN[JP]*p->DZN[KP];

        vVolAvg += vmag2*vol;
        volTot += vol; 
    }

    // Sum up vorticity and volume and send to rank 0 
    double dissipation = 0.0;
    double volume = 0.0;
    MPI_Reduce(&vVolAvg, &dissipation, 1, MPI_DOUBLE, MPI_SUM, 0, pgc->mpi_comm);
    MPI_Reduce(&volTot, &volume, 1, MPI_DOUBLE, MPI_SUM, 0, pgc->mpi_comm);

    if(p->mpirank == 0)
    {
        dissipation = p->W2*dissipation/volume;
        
        ofstream print;
        print.open("./REEF3D_CFD_TaylorGreen_Diss/REEF3D_CFD_TG_Diss.dat", ofstream::app);
        print<<p->simtime<<" \t "<<dissipation<<endl;
        print.close();
    }
}
