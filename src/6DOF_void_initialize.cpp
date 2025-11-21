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
--------------------------------------------------------------------*/

#include"6DOF_void.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"net_interface.h"
#include<sys/stat.h>
#include"mooring_void.h"
#include"mooring_barQuasiStatic.h"
#include"mooring_Catenary.h"
#include"mooring_Spring.h"
#include"mooring_dynamic.h"

void sixdof_void::initialize(lexer *p, fdm2D *b, ghostcell *pgc)
{
}

void sixdof_void::initialize(lexer *p, fdm *a, ghostcell *pgc)
{
    //if(p->mpirank==0)
    //cout<<"6DOF: ini_CFD"<<endl;
    
    
	if(p->X310==0)
	{
		pmooring.push_back(new mooring_void());
	}
	else
	{
		pgc->bcast_int(&p->mooring_count,1);	
        
        Xme.resize(p->mooring_count);
		Yme.resize(p->mooring_count);
		Zme.resize(p->mooring_count);
		Kme.resize(p->mooring_count);
		Mme.resize(p->mooring_count);
		Nme.resize(p->mooring_count);

		if(p->mpirank==0)
		{
			mkdir("./REEF3D_CFD_6DOF_Mooring",0777);	
		}		

		pmooring.reserve(p->mooring_count);
			
		for (int i=0; i < p->mooring_count; i++)
		{
			if(p->X310==1)
			{
				pmooring.push_back(new mooring_Catenary(i));
			}	
			else if(p->X310==2)
			{
				pmooring.push_back(new mooring_barQuasiStatic(i)); 
			}	
			else if(p->X310==3)
			{
				pmooring.push_back(new mooring_dynamic(i));
			}
			else if(p->X310==4)
			{
                pmooring.push_back(new mooring_Spring(i));
			}
		
			pmooring[i]->initialize(p,pgc);
		}
	}

    pnetinter->initialize_cfd(p,a,pgc);

    // Ini parameters
	p->ufbi=p->vfbi=p->wfbi=0.0;
	p->pfbi=p->qfbi=p->rfbi=0.0;
    p->xg=p->yg=p->zg=0.0;

    quatRotMat << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
}


void sixdof_void::initialize(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(p->mpirank==0)
    cout<<"6DOF: ini"<<endl;
    
    // Mooring
	if(p->X310==0)
	{
		pmooring.push_back(new mooring_void());
	}
    
	else
	{
		pgc->bcast_int(&p->mooring_count,1);	

		Xme.resize(p->mooring_count);
		Yme.resize(p->mooring_count);
		Zme.resize(p->mooring_count);
		Kme.resize(p->mooring_count);
		Mme.resize(p->mooring_count);
		Nme.resize(p->mooring_count);

		if(p->mpirank==0)
		mkdir("./REEF3D_NHFLOW_6DOF_Mooring",0777);	

		pmooring.reserve(p->mooring_count);
		X311_xen.resize(p->mooring_count,0.0);
		X311_yen.resize(p->mooring_count,0.0);
		X311_zen.resize(p->mooring_count,0.0);
			
		for (int i=0; i < p->mooring_count; i++)
		{
			if(p->X310==1)
			{
				pmooring.push_back(new mooring_Catenary(i));
			}	
			else if(p->X310==2)
			{
				pmooring.push_back(new mooring_barQuasiStatic(i)); 
			}	
			else if(p->X310==3)
			{
                pmooring.push_back(new mooring_dynamic(i));
			}
			else if(p->X310==4)
			{
				pmooring.push_back(new mooring_Spring(i));
			}
		
			X311_xen[i] = p->X311_xe[i] - p->xg;
			X311_yen[i] = p->X311_ye[i] - p->yg;
			X311_zen[i] = p->X311_ze[i] - p->zg;
		
			pmooring[i]->initialize(p,pgc);
		}
	}	
    
    // Net
    pnetinter->initialize_nhflow(p,d,pgc);

    // Ini parameters
	p->ufbi=p->vfbi=p->wfbi=0.0;
	p->pfbi=p->qfbi=p->rfbi=0.0;
    p->xg=p->yg=p->zg=0.0;

    quatRotMat << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
}