/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"nhflow_vtp_fsf.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"


void nhflow_vtp_fsf::pvd(lexer *p, fdm_nhf *d, ghostcell* pgc)
{	
    int num=0;
    if(p->P15==1)
    num = printcount;

    if(p->P15==2)
    num = p->count;
    
    
    if(p->count==0)
    {
	sprintf(name,"./REEF3D_NHFLOW_VTP_FSF/REEF3D-NHFLOW-FSF.pvd",num);

	pvdout.open(name);

	pvdout<<"<?xml version=\"1.0\"?>"<<endl;
	pvdout<<"<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">"<<endl;
	}
	
    // timestep="0.000800035" part="0" file="pvdoutput/pvdoutput_0_8.pvtu"/>
    sprintf(name,"./REEF3D_NHFLOW_VTP_FSF/REEF3D-NHFLOW-FSF-%08i.pvtp",num);
    
	pvdout<<"<timestep=\""<<p->simtime<<"\" part=\"0\" file="<<name<<"\"/>"<<endl;


/*
	pvdout<<"</Collection>"<<endl;
	pvdout<<"</VTKFile>"<<endl;

    
	pvdout.close();*/

}

