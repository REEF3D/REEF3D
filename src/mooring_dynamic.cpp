/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2026 Tobias Martin

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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"mooring_dynamic.h"
#include"lexer.h"
#include"ghostcell.h"

mooring_dynamic::mooring_dynamic(int number):line(number),beam(number)
{}

mooring_dynamic::~mooring_dynamic(){}

void mooring_dynamic::start(lexer *p, ghostcell *pgc)
{
	// Set mooring time step
	phi_mooring = 0.0;
	double t_new = phi_mooring*p->simtime + (1.0 - phi_mooring)*(p->simtime + p->dt);

	// start() is called at every RK substage of the 6DOF solver, but simtime only
	// advances once per time step: integrate only once per step, otherwise the
	// beam solver is called with t_old == t_new (zero step, division by zero)
	if (t_new - t_mooring <= 1.0e-12*MAX(1.0,fabs(t_new)))
	return;

	t_mooring_n = t_mooring;
	t_mooring = t_new;

	// Update fields
    updateFields(p, pgc);
	
	// Update boundary conditions
    fixPoint << p->X311_xe[line], p->X311_ye[line], p->X311_ze[line]; 

    // Integrate from t_mooring_n to t_mooring
    Integrate(t_mooring_n,t_mooring);

	// Save mooring point
	saveMooringPoint(p);

	// Plot mooring line	
	print(p);
}

void mooring_dynamic::updateFluidVel(lexer *p, ghostcell *pgc, int cmp)
{
}


void mooring_dynamic::updateFields(lexer *p, ghostcell *pgc)
{
    // Get position of points
    getTransPos(c_moor);
	
    // Fluid velocity
	updateFluidVel(p, pgc, 0);
	updateFluidVel(p, pgc, 1);
	updateFluidVel(p, pgc, 2);	
	
    // Fluid acceleration
	for (int i = 0; i < Ne + 1; i++)
	{
        fluid_acc[i][0] = (fluid_vel[i][0] - fluid_vel_n[i][0])/p->dt;
        fluid_acc[i][1] = (fluid_vel[i][1] - fluid_vel_n[i][1])/p->dt;
        fluid_acc[i][2] = (fluid_vel[i][2] - fluid_vel_n[i][2])/p->dt;
	}
}


void mooring_dynamic::mooringForces
(
	double& Xme, double& Yme, double& Zme
)
{
    // Tension forces if line is not broken
    if (broken==false)
    {
        Xme = Xme_; 
        Yme = Yme_;
        Zme = Zme_;
    }

    // Breakage due to max tension force
    if (breakTension > 0.0 && fabs(getTensLoc(Ne)) >= breakTension)
    {
        Xme = 0.0; 
        Yme = 0.0;
        Zme = 0.0;

        broken = true;
    }

    // Breakage due to time limit
    if (breakTime > 0.0 && t_mooring >= breakTime)
    {
        Xme = 0.0; 
        Yme = 0.0;
        Zme = 0.0;

        broken = true;
    }
}


void mooring_dynamic::saveMooringPoint(lexer *p)
{
    // Get position and velocity of points
    getTransPos(c_moor);
    getTransVel(cdot_moor);

    // Save location of line 
    c_moor_n = c_moor;

	// Save acceleration of line 
    cdotdot_moor = (cdot_moor - cdot_moor_n)/p->dt;
    cdot_moor_n = cdot_moor;

	// Save reaction forces at mooring point
    if (p->mpirank==0)
    {
		eTout<<p->simtime<<" \t "<<getTensLoc(Ne)<<endl;
    }
    
    Eigen::Vector3d tension = getTensGlob(Ne);
    Xme_ = -tension(0);
	Yme_ = -tension(1);
	Zme_ = -tension(2);
}
