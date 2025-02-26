/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include<sys/stat.h>
#include<iostream>
#include<fstream>
#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_obj::print_ini_vtp(lexer *p, ghostcell *pgc)
{
	if(p->mpirank==0)
    {
        if(p->A10==2)
        {
        mkdir("./REEF3D_SFLOW_6DOF_VTP", 0777);
        mkdir("./REEF3D_SFLOW_6DOF", 0777);
        }
        
        if(p->A10==5)
        {
        mkdir("./REEF3D_NHFLOW_6DOF_VTP", 0777);
        mkdir("./REEF3D_NHFLOW_6DOF", 0777);
        mkdir("./REEF3D_NHFLOW_6DOF_Normals_VTP", 0777);
        }
        
        if(p->A10==6)
        {
        mkdir("./REEF3D_CFD_6DOF_VTP", 0777);
        mkdir("./REEF3D_CFD_6DOF_Normals_VTP", 0777);
        mkdir("./REEF3D_CFD_6DOF", 0777);
        }
    }
	
    ofstream print;
    char str[1000];

    // position
    if(p->A10==2)
    sprintf(str,"./REEF3D_SFLOW_6DOF/REEF3D_6DOF_position_%i.dat",n6DOF);
    
    if(p->A10==5)
    sprintf(str,"./REEF3D_NHFLOW_6DOF/REEF3D_6DOF_position_%i.dat",n6DOF);
    
    if(p->A10==6)
    sprintf(str,"./REEF3D_CFD_6DOF/REEF3D_6DOF_position_%i.dat",n6DOF);
	
    printpos.open(str);
	printpos<<"time \t XG \t YG \t ZG \t Phi \t Theta \t Psi"<<endl;

    
    // velocity
    if(p->A10==2)
    sprintf(str,"./REEF3D_SFLOW_6DOF/REEF3D_6DOF_velocity_%i.dat",n6DOF);
    
    if(p->A10==5)
    sprintf(str,"./REEF3D_NHFLOW_6DOF/REEF3D_6DOF_velocity_%i.dat",n6DOF);
    
    if(p->A10==6)
    sprintf(str,"./REEF3D_CFD_6DOF/REEF3D_6DOF_velocity_%i.dat",n6DOF);
	
    printvel.open(str);
	printvel<<"time \t Ue [m/s] \t Ve [m/s] \t We [m/s] \t Pe [rad/s] \t Qe [rad/s] \t Re [rad/s]"<<endl;

    // force
    if(p->A10==2)
    sprintf(str,"./REEF3D_SFLOW_6DOF/REEF3D_6DOF_forces_%i.dat",n6DOF);
    
    if(p->A10==5)
    sprintf(str,"./REEF3D_NHFLOW_6DOF/REEF3D_6DOF_forces_%i.dat",n6DOF);
    
    if(p->A10==6)
    sprintf(str,"./REEF3D_CFD_6DOF/REEF3D_6DOF_forces_%i.dat",n6DOF);
	
    printforce.open(str);
	printforce<<"time \t Fx \t Fy \t Fz \t Mx \t My \t Mz \t Fx_p \t Fy_p \t Fz_p \t Fx_v \t Fy_v \t Fz_v"<<endl;
    

    curr_time = 0.0;
    
    p->Darray(printtime_wT,p->P35);

    for(int qn=0; qn<p->P35; ++qn)
	printtime_wT[qn]=p->P35_ts[qn];
}

void sixdof_obj::print_ini_stl(lexer *p, ghostcell *pgc)
{
	if(p->mpirank==0)
    {
        if(p->A10==2)
        {
        mkdir("./REEF3D_SFLOW_6DOF_VTP", 0777);
        mkdir("./REEF3D_SFLOW_6DOF", 0777);
        }
        
        if(p->A10==5)
        {
        mkdir("./REEF3D_NHFLOW_6DOF_VTP", 0777);
        mkdir("./REEF3D_NHFLOW_6DOF", 0777);
        }
        
        if(p->A10==6)
        {
        mkdir("./REEF3D_CFD_6DOF_VTP", 0777);
        mkdir("./REEF3D_CFD_6DOF_Normals_VTP", 0777);
        mkdir("./REEF3D_CFD_6DOF", 0777);
        }
    }
	
    ofstream print;
    char str[1000];
    
    if(p->A10==2)
    sprintf(str,"./REEF3D_SFLOW_6DOF/REEF3D_6DOF_position_%i.dat",n6DOF);
    
    if(p->A10==5)
    sprintf(str,"./REEF3D_NHFLOW_6DOF/REEF3D_6DOF_position_%i.dat",n6DOF);
    
    if(p->A10==6)
    sprintf(str,"./REEF3D_CFD_6DOF/REEF3D_6DOF_position_%i.dat",n6DOF);
	
    print.open(str);
	print<<"time \t XG \t YG \t ZG \t Phi \t Theta \t Psi"<<endl;
	print.close();



    if(p->A10==2)
    sprintf(str,"./REEF3D_SFLOW_6DOF/REEF3D_6DOF_velocity_%i.dat",n6DOF);
    
    if(p->A10==5)
    sprintf(str,"./REEF3D_NHFLOW_6DOF/REEF3D_6DOF_velocity_%i.dat",n6DOF);
    
    if(p->A10==6)
    sprintf(str,"./REEF3D_CFD_6DOF/REEF3D_6DOF_velocity_%i.dat",n6DOF);
	
    print.open(str);
	print<<"time \t Ue [m/s] \t Ve [m/s] \t We [m/s] \t Pe [rad/s] \t Qe [rad/s] \t Re [rad/s]"<<endl;
    print.close();
    

    
    if(p->A10==2)
    sprintf(str,"./REEF3D_SFLOW_6DOF/REEF3D_6DOF_forces_%i.dat",n6DOF);
    
    if(p->A10==5)
    sprintf(str,"./REEF3D_NHFLOW_6DOF/REEF3D_6DOF_forces_%i.dat",n6DOF);
    
    if(p->A10==6)
    sprintf(str,"./REEF3D_CFD_6DOF/REEF3D_6DOF_forces_%i.dat",n6DOF);
	
    print.open(str);
	print<<"time \t Fx \t Fy \t Fz \t Mx \t My \t Mz \t Fx_p \t Fy_p \t Fz_p \t Fx_v \t Fy_v \t Fz_v"<<endl;
    print.close();    

    curr_time = 0.0;
    
    p->Darray(printtime_wT,p->P35);

    for(int qn=0; qn<p->P35; ++qn)
	printtime_wT[qn]=p->P35_ts[qn];
}