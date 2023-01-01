/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2023 Tobias Martin

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

#include"mooring_Spring.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

mooring_Spring::mooring_Spring(int number):line(number){}

mooring_Spring::~mooring_Spring(){}


void mooring_Spring::initialize(lexer *p, fdm *a, ghostcell *pgc)
{   
    curr_time = p->simtime;

    dx = p->X311_xe[line] - p->X311_xs[line];			
    dy = p->X311_ye[line] - p->X311_ys[line];				
    dz = p->X311_ze[line] - p->X311_zs[line];	
	L0 = sqrt(dx*dx + dy*dy + dz*dz);	
   
	k = p->X312_k[line];
    T0 = p->X312_T0[line];

	if(p->mpirank==0 && p->P14==1)
	{
		char str[1000];
		sprintf(str,"./REEF3D_CFD_6DOF/REEF3D_6DOF_mooring_force_%i.dat",line);
		eTout.open(str);
		eTout<<"time \t T"<<endl;	
	}
    
	printtime = 0.0;     
    
    // Initialise breaking
    broken = false;
    curr_time = 0.0;
    breakTension = p->X314 > 0 ? p->X314_T[line]: 0.0;
    breakTime = p->X315 > 0 ? p->X315_t[line]: 0.0;
}


void mooring_Spring::start(lexer *p, fdm *a, ghostcell *pgc)
{
	//- Calculate distance between start and mooring points
    dx = p->X311_xe[line] - p->X311_xs[line];			
    dy = p->X311_ye[line] - p->X311_ys[line];				
    dz = p->X311_ze[line] - p->X311_zs[line];	
	L = sqrt(dx*dx + dy*dy + dz*dz);			
    
    double dL = L - L0;

	//- Calculate tension force in spring
    T = T0;

    if (dL > 0.0)
    {
        T = T0 + k*dL;
    }
    else if (dL < 0.0)
    {
        T = MAX(T0 - k*fabs(dL),0.0);
    }

	//- Calculate reaction forces at mooring points	
    Xme_ = T*fabs(dx/L);
    Yme_ = T*fabs(dy/L);
    Zme_ = T*fabs(dz/L);

	if (dx > 0)	
	{
		Xme_ *= -1.0;
	}
	if (dy > 0)	
	{
		Yme_ *= -1.0;
	}
	if (dz > 0)	
	{
		Zme_ *= -1.0;
	}

	//- Print spring
	print(p);
}


void mooring_Spring::mooringForces
(
	double& Xme, double& Yme, double& Zme
)
{
    // Tension forces if line is not broken
    if (broken == false)
    {
        Xme = Xme_; 
        Yme = Yme_;
        Zme = Zme_;
    }

    // Breakage due to max tension force
    if (breakTension > 0.0 && fabs(T) >= breakTension)
    {
        Xme = 0.0; 
        Yme = 0.0;
        Zme = 0.0;

        broken = true;
    }

    // Breakage due to time limit
    if (breakTime > 0.0 && curr_time >= breakTime)
    {
        Xme = 0.0; 
        Yme = 0.0;
        Zme = 0.0;

        broken = true;
    }
}


