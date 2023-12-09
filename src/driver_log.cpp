/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"driver.h"
#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"
#include"freesurface_header.h"
#include"turbulence_header.h"
#include"momentum_header.h"
#include"pressure_header.h"
#include"fdm_header.h"
#include"sediment_header.h"
#include"convection_header.h"
#include"solver_header.h"
#include"field_header.h"
#include"heat_header.h"
#include"concentration_header.h"
#include"benchmark_header.h"
#include"6DOF_header.h"
#include"lexer.h"
#include"cart1.h"
#include"cart2.h"
#include"cart3.h"
#include"cart4.h"
#include"cart4a.h"
#include<sys/stat.h>
#include<sys/types.h>

void driver::log_ini()
{

	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_Log",0777);

    if(p->mpirank==0)
    {
    if(p->P14==0)
    mainlogout.open("REEF3D_mainlog.dat");
    if(p->P14==1)
    mainlogout.open("./REEF3D_Log/REEF3D_mainlog.dat");
    mainlogout<<"REEF3D version:  "<<version<<endl<<endl;
    mainlogout<<"number of cells:  "<<p->cellnumtot<<endl<<endl;
    mainlogout<<"#iteration \t #timestep \t #simtime \t #itertime \t #piter \t #ptime \t #Volume 1 \t #Volume2 \t #Inflow \t #Outflow \t #Ui \t #Phimean \t #Phiout "<<endl;
    }

    if(p->mpirank==0)
    {
    if(p->P14==0)
    maxlogout.open("REEF3D_maxlog.dat");
    if(p->P14==1)
    maxlogout.open("./REEF3D_Log/REEF3D_maxlog.dat");

    maxlogout<<"number of cells:  "<<p->cellnumtot<<endl<<endl;
    maxlogout<<"#iteration \t #umax \t\t #vmax \t\t #wmax \t\t #viscmax \t\t #kinmax \t\t #epsmax \t\t #pressmax "<<endl;
    }

    if(p->mpirank==0)
    {
    if(p->P14==0)
    solvlogout.open("REEF3D_solverlog.dat");
    if(p->P14==1)
    solvlogout.open("./REEF3D_Log/REEF3D_solverlog.dat");

    solvlogout<<"number of cells:  "<<p->cellnumtot<<endl<<endl;
    solvlogout<<"#iteration \t #itertime \t #totaltime \t |#piter \t #ptime \t| #uiter \t #utime \t| #viter \t #vtime \t| #witer \t #wtime \t|";
    solvlogout<<"#kiter \t #ktime \t| #eiter \t #etime \t|";
    solvlogout<<"#liter \t #ltime \t| #reiniiter \t #reinitime"<<endl;
    }

}

void driver::mainlog(lexer *p)
{
	 if(p->count%p->P12==0)
	 {
     mainlogout<<fixed<<p->count<<" \t "<<setprecision(5)<<p->dt<<" \t "<<setprecision(5)<<p->simtime<<" \t ";
	 mainlogout<<fixed<<setprecision(4)<<p->itertime<<" \t ";
	 mainlogout<<p->poissoniter<<" \t "<<setprecision(4)<<p->poissontime<<" \t ";
     mainlogout<<fixed<<setprecision(4)<<p->volume1<<" \t "<<setprecision(4)<<p->volume2<<" \t ";
     mainlogout<<fixed<<setprecision(6)<<p->Qi<<" \t "<<setprecision(6)<<p->Qo<<" \t ";
	 mainlogout<<fixed<<setprecision(4)<<p->Ui<<" \t "<<setprecision(6)<<p->phimean<<" \t "<<setprecision(6)<<p->phiout;
	 mainlogout<<endl;
	 }
}

void driver::maxlog(lexer *p)
{
	 if(p->count%p->P12==0)
	 {
     maxlogout<<p->count<<"\t \t"<<p->umax<<" \t "<<setprecision(4)<<p->vmax<<" \t "<<setprecision(4)<<p->wmax<<" \t ";
     maxlogout<<setprecision(4)<<p->viscmax<<" \t "<<setprecision(4)<<p->kinmax<<" \t "<<setprecision(4)<<p->epsmax<<" \t ";
     maxlogout<<setprecision(4)<<p->pressmax<<endl;
	 }
}

void driver::solverlog(lexer* p)
{
	 if(p->count%p->P12==0)
	 {
     solvlogout<<p->count<<"\t \t";
	 solvlogout<<setprecision(4)<<p->itertime<<" \t ";
	 solvlogout<<setprecision(4)<<p->totaltime<<" \t ";
	 solvlogout<<p->poissoniter<<" \t "<<setprecision(4)<<p->poissontime<<" \t ";
     solvlogout<<p->uiter<<" \t "<<setprecision(4)<<p->utime<<" \t ";
     solvlogout<<p->viter<<" \t "<<setprecision(4)<<p->vtime<<" \t ";
     solvlogout<<p->witer<<" \t "<<setprecision(4)<<p->wtime<<" \t ";

     solvlogout<<p->kiniter<<" \t "<<setprecision(4)<<p->kintime<<" \t ";
     solvlogout<<p->epsiter<<" \t "<<setprecision(4)<<p->epstime<<" \t ";
     
     solvlogout<<p->lsmiter<<" \t "<<setprecision(4)<<p->lsmtime<<" \t ";
     solvlogout<<p->F44<<" \t "<<setprecision(4)<<p->reinitime<<" \t "<<endl;
	 }
}





