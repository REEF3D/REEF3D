/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
#include"fdm.h"
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

void sixdof_obj::initialize_cfd(lexer *p, fdm *a, ghostcell *pgc, vector<net*>& pnet)
{
    if(p->mpirank==0)
    cout<<"6DOF_df_ini "<<endl;
    
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
    ray_cast(p,a,pgc);
    //if(p->X188==1)
	reini_RK2(p,a,pgc,a->fb);
    pgc->start4a(p,a->fb,50);
    
    // Calculate geometrical properties
	geometry_parameters(p,a,pgc);
    
    // Initialise position of bodies
    iniPosition_RBM(p,pgc);
	
	// Recalculate distances
	ray_cast(p,a,pgc);
	reini_RK2(p,a,pgc,a->fb);
    pgc->start4a(p,a->fb,50);
    
    // Initialise global variables
	update_fbvel(p,pgc);
   
    // Initialise floating fields
     ULOOP
     a->fbh1(i,j,k) = Hsolidface(p,a,1,0,0);

     VLOOP
     a->fbh2(i,j,k) = Hsolidface(p,a,0,1,0);

     WLOOP
     a->fbh3(i,j,k) = Hsolidface(p,a,0,0,1);

     LOOP
     a->fbh4(i,j,k) = Hsolidface(p,a,0,0,0);

     pgc->start1(p,a->fbh1,10);
     pgc->start2(p,a->fbh2,11);
     pgc->start3(p,a->fbh3,12);
     pgc->start4(p,a->fbh4,40);

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
		mkdir("./REEF3D_CFD_6DOF_Mooring",0777);	

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
    }
    
    // ghostcell update
    pgc->gcdf_update(p,a);
}

void sixdof_obj::ini_parallel(lexer *p, ghostcell *pgc)
{
    p->Darray(xstart, p->mpi_size);
    p->Darray(xend, p->mpi_size);
    p->Darray(ystart, p->mpi_size);
    p->Darray(yend, p->mpi_size);
    p->Darray(zstart, p->mpi_size);
    p->Darray(zend, p->mpi_size);
    
    xstart[p->mpirank] = p->originx;
    ystart[p->mpirank] = p->originy;
    zstart[p->mpirank] = p->originz;
    xend[p->mpirank] = p->endx;
    yend[p->mpirank] = p->endy;
    zend[p->mpirank] = p->endz;
    
    for (int i = 0; i < p->mpi_size; i++)
    {
        MPI_Bcast(&xstart[i],1,MPI_DOUBLE,i,pgc->mpi_comm);
        MPI_Bcast(&xend[i],1,MPI_DOUBLE,i,pgc->mpi_comm);
        MPI_Bcast(&ystart[i],1,MPI_DOUBLE,i,pgc->mpi_comm);
        MPI_Bcast(&yend[i],1,MPI_DOUBLE,i,pgc->mpi_comm);
        MPI_Bcast(&zstart[i],1,MPI_DOUBLE,i,pgc->mpi_comm);
        MPI_Bcast(&zend[i],1,MPI_DOUBLE,i,pgc->mpi_comm);
    }
}    

double sixdof_obj::ramp_vel(lexer *p)
{
    double f=1.0;
    
    if(p->X205==1 && p->X206==1 && p->simtime>=p->X206_ts && p->simtime<p->X206_te)
    {
    f = (p->simtime-p->X206_ts)/(p->X206_te-p->X206_ts);
    }
    
    if(p->X205==2 && p->X206==1 && p->simtime>=p->X206_ts && p->simtime<p->X206_te)
    {
    f = (p->simtime-p->X206_ts)/(p->X206_te-p->X206_ts)-(1.0/PI)*sin(PI*(p->simtime-p->X206_ts)/(p->X206_te-p->X206_ts));
    }
    
    if(p->X206==1 && p->simtime<p->X206_ts)
    f=0.0;
    
    return f;
}

double sixdof_obj::ramp_draft(lexer *p)
{
    double f=1.0;
    
    if(p->X205==1 && p->X207==1 && p->simtime>=p->X207_ts && p->simtime<p->X207_te)
    {
    f = p->simtime/(p->X207_te-p->X207_ts);
    }
    
    if(p->X205==2 && p->X207==1 && p->simtime>=p->X207_ts && p->simtime<p->X207_te)
    {
    f = p->simtime/(p->X207_te-p->X207_ts) - (1.0/PI)*sin(PI*(p->simtime/(p->X207_te-p->X207_ts)));
    }
    
    if(p->X207==1 && p->simtime<p->X207_ts)
    f=0.0;
    
    return f;
}