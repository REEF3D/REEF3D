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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"6DOF_obj.h"
#include"lexer.h"
#include"momentum.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include"mooring_void.h"
#include"mooring_barQuasiStatic.h"
#include"mooring_Catenary.h"
#include"mooring_Spring.h"
#include"mooring_dynamic.h"
#include"net.h"
#include"net_void.h"
#include"net_barQuasiStatic.h"
#include"net_barDyn.h"
#include"net_sheet.h"

void sixdof_obj::initialize_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, vector<net*>& pnet)
{
    if(p->mpirank==0)
    cout<<"6DOF_df_ini "<<endl;
    
    if(p->mpirank==0)
    mkdir("./REEF3D_NHFLOW_6DOF",0777);
    
    // Initialise folder structure
    if(p->X50==1)
	print_ini_vtp(p,pgc);
    
    if(p->X50==2)
    print_ini_stl(p,pgc);
 
    // Initialise processor boundaries
    ini_parallel(p,pgc);
    
    // Initialise objects
	objects_create(p,pgc);
    
    // Initialise fbvel
	ini_fbvel(p,pgc);
    
    // Level Set for floating body
    ray_cast(p,d,pgc);
	nhflow_reini_RK2(p,d,pgc,d->FB);
    
    // Calculate geometrical properties
    geometry_parameters_nhflow(p,d,pgc);
    
    // Initialise position of bodies
    iniPosition_RBM(p,pgc);
	
	// Recalculate distances
	ray_cast(p,d,pgc);
	nhflow_reini_RK2(p,d,pgc,d->FB);
    
    // Initialise global variables
	update_fbvel(p,pgc);

    // Print initial body 
    if(p->X50==1)
    print_vtp(p,pgc);
    
    if(p->X50==2)
    print_stl(p,pgc);




	// Mooring
	if(p->X310==0)
	{
		pmooring.push_back(new mooring_void());
	}
    
	else
	{
		MPI_Bcast(&p->mooring_count,1,MPI_DOUBLE,0,pgc->mpi_comm);	

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

/*
    // Net
    if (p->X320==0)
    {
        pnet.push_back(new net_void());
    }
    else
    {
		MPI_Bcast(&p->net_count,1,MPI_DOUBLE,0,pgc->mpi_comm);
        
        Xne.resize(p->net_count);
		Yne.resize(p->net_count);
		Zne.resize(p->net_count);
		Kne.resize(p->net_count);
		Mne.resize(p->net_count);
		Nne.resize(p->net_count);

        if(p->mpirank==0)
        {
            mkdir("./REEF3D_CFD_6DOF_Net",0777);	
        }
        else
        {
            p->X320_type = new int[p->net_count];
        }
		
        pnet.reserve(p->net_count);	
  
		for (int ii=0; ii < p->net_count; ii++)
		{
            MPI_Bcast(&p->X320_type[ii],1,MPI_INT,0,pgc->mpi_comm);
			
            if(p->X320_type[ii] > 10)
			{
				pnet.push_back(new net_barDyn(ii,p));
			}
            else if (p->X320_type[ii] < 4)
            {
                pnet.push_back(new net_barQuasiStatic(ii,p));
            }
            else
            {
                 pnet.push_back(new net_sheet(ii,p));
            }
			
            pnet[ii]->initialize(p,a,pgc);
		}
    }*/
    
}


