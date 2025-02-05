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
--------------------------------------------------------------------*/

#include"6DOF_void.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"vrans.h"
#include<sys/stat.h>

#include"mooring_void.h"
#include"mooring_barQuasiStatic.h"
#include"mooring_Catenary.h"
#include"mooring_Spring.h"
#include"mooring_dynamic.h"
#include"net.h"
#include"net_void.h"
#include"net_barDyn.h"
#include"net_barQuasiStatic.h"
#include"net_sheet.h"

void sixdof_void::initialize(lexer *p, fdm2D *b, ghostcell *pgc, vector<net*>& pnet)
{
}

void sixdof_void::initialize(lexer *p, fdm *a, ghostcell *pgc, vector<net*>& pnet)
{
    //if(p->mpirank==0)
    //cout<<"6DOF: ini_CFD"<<endl;
    
    
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

    if (p->X320 == 0)
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
        mkdir("./REEF3D_CFD_6DOF_Net",0777);	

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

    // Ini parameters
	p->ufbi=p->vfbi=p->wfbi=0.0;
	p->pfbi=p->qfbi=p->rfbi=0.0;
    p->xg=p->yg=p->zg=0.0;

    quatRotMat << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
}


void sixdof_void::initialize(lexer *p, fdm_nhf *d, ghostcell *pgc, vector<net*>& pnet)
{
    //if(p->mpirank==0)
    //cout<<"6DOF: ini"<<endl;
    
    /*
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

    if (p->X320 == 0)
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
        mkdir("./REEF3D_CFD_6DOF_Net",0777);	

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

    // Ini parameters
	p->ufbi=p->vfbi=p->wfbi=0.0;
	p->pfbi=p->qfbi=p->rfbi=0.0;
    p->xg=p->yg=p->zg=0.0;

    quatRotMat << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
}