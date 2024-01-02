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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"6DOF_df_object.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"mooring.h"
#include"net.h"
#include"vrans.h"

void sixdof_df_object::externalForces(lexer *p,fdm* a, ghostcell *pgc, double alpha, vrans *pvrans, vector<net*>& pnet)
{
    Xext = Yext = Zext = Kext = Mext = Next = 0.0;

    // Mooring forces
	if (p->X310 > 0)
	{
		mooringForces(p,a,pgc,alpha);
	} 
	
    // Net forces
	if (p->X320 > 0)
	{
		netForces(p,a,pgc,alpha,pvrans,pnet);
	}
}

void sixdof_df_object::mooringForces(lexer *p, fdm* a, ghostcell *pgc, double alpha)
{
	for (int ii=0; ii<p->mooring_count; ii++)
	{
		// Update coordinates of end point
        Eigen::Vector3d point(X311_xen[ii], X311_yen[ii], X311_zen[ii]);
					
        point = R_*point;
					
        p->X311_xe[ii] = point(0) + c_(0);
        p->X311_ye[ii] = point(1) + c_(1);
        p->X311_ze[ii] = point(2) + c_(2);

        // Advance in time
        pmooring[ii]->start(p, a, pgc);
                
        // Get forces
        pmooring[ii]->mooringForces(Xme[ii],Yme[ii],Zme[ii]);
                
        // Calculate moments
        Kme[ii] = (p->X311_ye[ii] - c_(1))*Zme[ii] - (p->X311_ze[ii] - c_(2))*Yme[ii];
        Mme[ii] = (p->X311_ze[ii] - c_(2))*Xme[ii] - (p->X311_xe[ii] - c_(0))*Zme[ii];
        Nme[ii] = (p->X311_xe[ii] - c_(0))*Yme[ii] - (p->X311_ye[ii] - c_(1))*Xme[ii];
            
        /*    
        if( p->mpirank == 0)
        {
            cout<<"Xme"<< ii <<" : "<<Xme[ii]<<" Yme"<< ii <<" : "<<Yme[ii]<<" Zme"<< ii <<" : "<<Zme[ii]
            <<" Kme"<< ii <<" : "<<Kme[ii]<<" Mme"<< ii <<" : "<<Mme[ii]<<" Nme"<< ii <<" : "<<Nme[ii]<<endl;		
        }*/
            
        // Distribute forces and moments to all processors
        MPI_Bcast(&Xme[ii],1,MPI_DOUBLE,0,pgc->mpi_comm);
        MPI_Bcast(&Yme[ii],1,MPI_DOUBLE,0,pgc->mpi_comm);
        MPI_Bcast(&Zme[ii],1,MPI_DOUBLE,0,pgc->mpi_comm);
        MPI_Bcast(&Kme[ii],1,MPI_DOUBLE,0,pgc->mpi_comm);
        MPI_Bcast(&Mme[ii],1,MPI_DOUBLE,0,pgc->mpi_comm);
        MPI_Bcast(&Nme[ii],1,MPI_DOUBLE,0,pgc->mpi_comm);	
        
        // Add to external forces
        Xext += Xme[ii];
        Yext += Yme[ii];
        Zext += Zme[ii];
        
        Kext += Kme[ii];
        Mext += Mme[ii];
        Next += Nme[ii];
    }
}

void sixdof_df_object::netForces(lexer *p, fdm* a, ghostcell *pgc, double alpha, vrans *pvrans, vector<net*>& pnet)
{
    for (int ii = 0; ii < p->net_count; ii++)
    {
        pnet[ii]->start(p, a, pgc, alpha, quatRotMat);
        pvrans->start(p, a, pgc, pnet[ii], ii);
    
        // Forces on rigid body
        pnet[ii]->netForces(p,Xne[ii],Yne[ii],Zne[ii],Kne[ii],Mne[ii],Nne[ii]);
    
        /*
        if( p->mpirank == 0)
        {
            cout<<"Xne"<< ii <<" : "<<Xne[ii]<<" Yne"<< ii <<" : "<<Yne[ii]<<" Zne"<< ii <<" : "<<Zne[ii]
            <<" Kne"<< ii <<" : "<<Kne[ii]<<" Mne"<< ii <<" : "<<Mne[ii]<<" Nne"<< ii <<" : "<<Nne[ii]<<endl;		
        }*/
        
        // Distribute forces and moments to all processors
        MPI_Bcast(&Xne[ii],1,MPI_DOUBLE,0,pgc->mpi_comm);
        MPI_Bcast(&Yne[ii],1,MPI_DOUBLE,0,pgc->mpi_comm);
        MPI_Bcast(&Zne[ii],1,MPI_DOUBLE,0,pgc->mpi_comm);
        MPI_Bcast(&Kne[ii],1,MPI_DOUBLE,0,pgc->mpi_comm);
        MPI_Bcast(&Mne[ii],1,MPI_DOUBLE,0,pgc->mpi_comm);
        MPI_Bcast(&Nne[ii],1,MPI_DOUBLE,0,pgc->mpi_comm);	
        
        // Add to external forces
        Xext += Xne[ii];
        Yext += Yne[ii];
        Zext += Zne[ii];
        Kext += Kne[ii];
        Mext += Mne[ii];
        Next += Nne[ii];
    }
}	

