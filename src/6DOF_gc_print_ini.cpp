/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include<iostream>
#include<fstream>
#include<sys/stat.h>
#include<sys/types.h>
#include"6DOF_gc.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_gc::print_ini_vtp(lexer *p, fdm *a, ghostcell *pgc)
{
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_CFD_6DOF_VTP",0777);
	
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_CFD_6DOF",0777);
	
	
	if(p->P14==0)
	eposout.open("REEF3D_6DOF_E_position.dat");
	if(p->P14==1)
	eposout.open("./REEF3D_CFD_6DOF/REEF3D_6DOF_E_position.dat");
	
	eposout<<"time \t XG \t YG \t ZG \t Phi \t Theta \t Psi"<<endl;
	
	
	
	if(p->P14==0)
	evelout.open("REEF3D_6DOF_E_velocity.dat");
	if(p->P14==1)
	evelout.open("./REEF3D_CFD_6DOF/REEF3D_6DOF_E_velocity.dat");
	
	evelout<<"time \t Ue \t Ve \t We \t Pe \t Qe \t Re"<<endl;
	
	
	if(p->P14==0)
	evelout.open("REEF3D_6DOF_E_force.dat");
	if(p->P14==1)
	eforceout.open("./REEF3D_CFD_6DOF/REEF3D_6DOF_E_force.dat");
	
	eforceout<<"time \t Xe \t Ye \t Ze \t Ke \t Me \t Ne"<<endl;
	
	
	if(p->P14==0)
	evelout.open("REEF3D_6DOF_S_force.dat");
	if(p->P14==1)
	sforceout.open("./REEF3D_CFD_6DOF/REEF3D_6DOF_S_force.dat");
	
	sforceout<<"time \t Xs \t Ys \t Zs \t Ks \t Ms \t Ns"<<endl;
	
}

void sixdof_gc::print_ini_stl(lexer *p, fdm *a, ghostcell *pgc)
{
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_CFD_6DOF_STL",0777);
	
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_CFD_6DOF",0777);
	
	
	if(p->P14==0)
	eposout.open("REEF3D_6DOF_E_position.dat");
	if(p->P14==1)
	eposout.open("./REEF3D_CFD_6DOF/REEF3D_6DOF_E_position.dat");
	
	eposout<<"time \t XG \t YG \t ZG \t Phi \t Theta \t Psi"<<endl;
	
	
	
	if(p->P14==0)
	evelout.open("REEF3D_6DOF_E_velocity.dat");
	if(p->P14==1)
	evelout.open("./REEF3D_CFD_6DOF/REEF3D_6DOF_E_velocity.dat");
	
	evelout<<"time \t Ue \t Ve \t We \t Pe \t Qe \t Re"<<endl;
	
	
	if(p->P14==0)
	evelout.open("REEF3D_6DOF_E_force.dat");
	if(p->P14==1)
	eforceout.open("./REEF3D_CFD_6DOF/REEF3D_6DOF_E_force.dat");
	
	eforceout<<"time \t Xe \t Ye \t Ze \t Ke \t Me \t Ne"<<endl;
	
	
	if(p->P14==0)
	evelout.open("REEF3D_6DOF_S_force.dat");
	if(p->P14==1)
	sforceout.open("./REEF3D_CFD_6DOF/REEF3D_6DOF_S_force.dat");
	
	sforceout<<"time \t Xs \t Ys \t Zs \t Ks \t Ms \t Ns"<<endl;
	
}
