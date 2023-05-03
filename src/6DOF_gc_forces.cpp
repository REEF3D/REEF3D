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

#include"6DOF_gc.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"vrans.h"
#include"mooring.h"
#include"net.h"

void sixdof_gc::fluidForces(lexer *p,fdm* a, ghostcell *pgc)
{
	if (p->X40 == 1)
	{
		forces_lsm(p,a,pgc);
	}
    else if (p->X40 == 2)
    {
        forces_triang_triangulation(p,a,pgc);
		forces_triang_reconstruction(p,a,pgc);
		print_forces_vtp(p,a,pgc);
		forces_triang(p,a,pgc);
		forces_triang_finalize(p,a,pgc);
    }
    else if (p->X40 == 3)
    {
		forces_stl(p,a,pgc);
	}
}


void sixdof_gc::mooringForces(lexer *p, fdm* a, ghostcell *pgc)
{
	for (int i=0; i<p->mooring_count; i++)
	{
		// Update coordinates of end point
		if(p->X13==0)
		{        
			rotation(p->X311_xe[i],p->X311_ye[i],p->X311_ze[i],dphi,dtheta,dpsi);
			
			p->X311_xe[i] += dxg;
			p->X311_ye[i] += dyg;
			p->X311_ze[i] += dzg;
		}
		else
		{
			std::vector<double> point(3,0.0);
				
			point[0] = X311_xen[i];
			point[1] = X311_yen[i];
			point[2] = X311_zen[i];
					
			point = rotation_R(point);
					
			p->X311_xe[i] = point[0] + xg;
			p->X311_ye[i] = point[1] + yg;
			p->X311_ze[i] = point[2] + zg;
		}

		// Advance in time
		pmooring[i]->start(p, a, pgc);
			
		// Get forces
		pmooring[i]->mooringForces(Xme[i],Yme[i],Zme[i]);
			
		// Calculate moments
		Kme[i] = (p->X311_ye[i] - yg)*Zme[i] - (p->X311_ze[i] - zg)*Yme[i];
		Mme[i] = (p->X311_ze[i] - zg)*Xme[i] - (p->X311_xe[i] - xg)*Zme[i];
		Nme[i] = (p->X311_xe[i] - xg)*Yme[i]- (p->X311_ye[i] - yg)*Xme[i];
		
		if( p->mpirank == 0)
		{
			cout<<"Xme"<< i <<" : "<<Xme[i]<<" Yme"<< i <<" : "<<Yme[i]<<" Zme"<< i <<" : "<<Zme[i]
			<<" Kme"<< i <<" : "<<Kme[i]<<" Mme"<< i <<" : "<<Mme[i]<<" Nme"<< i <<" : "<<Nme[i]<<endl;		
		}
		
		// Distribute forces and moments to all processors
		MPI_Bcast(&Xme[i],1,MPI_DOUBLE,0,pgc->mpi_comm);
		MPI_Bcast(&Yme[i],1,MPI_DOUBLE,0,pgc->mpi_comm);
		MPI_Bcast(&Zme[i],1,MPI_DOUBLE,0,pgc->mpi_comm);
		MPI_Bcast(&Kme[i],1,MPI_DOUBLE,0,pgc->mpi_comm);
		MPI_Bcast(&Mme[i],1,MPI_DOUBLE,0,pgc->mpi_comm);
		MPI_Bcast(&Nme[i],1,MPI_DOUBLE,0,pgc->mpi_comm);
	}
}


void sixdof_gc::netForces(lexer *p, fdm* a, ghostcell *pgc, double alpha, vrans *pvrans, vector<net*>& pnet)
{
	for (int ii = 0; ii < p->net_count; ii++)
	{
        // Advance in time	
        pnet[ii]->start(p, a, pgc, alpha, quatRotMat);
        pvrans->start(p, a, pgc, pnet[ii], ii);
        
        // Get forces
        pnet[ii]->netForces(p,Xne[ii],Yne[ii],Zne[ii],Kne[ii],Mne[ii],Nne[ii]);
        
        if( p->mpirank == 0)
        {
            cout<<"Xne"<< ii <<" : "<<Xne[ii]<<" Yne"<< ii <<" : "<<Yne[ii]<<" Zne"<< ii <<" : "<<Zne[ii]
            <<" Kne"<< ii <<" : "<<Kne[ii]<<" Mne"<< ii <<" : "<<Mne[ii]<<" Nne"<< ii <<" : "<<Nne[ii]<<endl;		
        }
		
		// Distribute forces and moments to all processors
		MPI_Bcast(&Xne[ii],1,MPI_DOUBLE,0,pgc->mpi_comm);
		MPI_Bcast(&Yne[ii],1,MPI_DOUBLE,0,pgc->mpi_comm);
		MPI_Bcast(&Zne[ii],1,MPI_DOUBLE,0,pgc->mpi_comm);
		MPI_Bcast(&Kne[ii],1,MPI_DOUBLE,0,pgc->mpi_comm);
		MPI_Bcast(&Mne[ii],1,MPI_DOUBLE,0,pgc->mpi_comm);
		MPI_Bcast(&Nne[ii],1,MPI_DOUBLE,0,pgc->mpi_comm);
	}
}	
