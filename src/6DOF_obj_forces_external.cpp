/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"mooring.h"

void sixdof_obj::externalForces_cfd(lexer *p, fdm* a, ghostcell *pgc, double alpha)
{
    Xext = Yext = Zext = Kext = Mext = Next = 0.0;

    // Mooring forces
	if (p->X310>0)
	mooringForces(p,pgc,alpha);
	
    // Net forces
	if (p->X320>0)
	netForces_cfd(p,a,pgc,alpha);
}

void sixdof_obj::externalForces_nhflow(lexer *p, fdm_nhf* d, ghostcell *pgc, double alpha)
{
    Xext = Yext = Zext = Kext = Mext = Next = 0.0;

    // Mooring forces
	if (p->X310>0)
	mooringForces(p,pgc,alpha);

	
    // Net forces
	if (p->X320>0)
	netForces_nhflow(p,d,pgc,alpha);
}

void sixdof_obj::mooringForces(lexer *p, ghostcell *pgc, double alpha)
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
        pmooring[ii]->start(p, pgc);
                
        // Get forces
        pmooring[ii]->mooringForces(Xme[ii],Yme[ii],Zme[ii]);
                
        // Calculate moments
        Kme[ii] = (p->X311_ye[ii] - c_(1))*Zme[ii] - (p->X311_ze[ii] - c_(2))*Yme[ii];
        Mme[ii] = (p->X311_ze[ii] - c_(2))*Xme[ii] - (p->X311_xe[ii] - c_(0))*Zme[ii];
        Nme[ii] = (p->X311_xe[ii] - c_(0))*Yme[ii] - (p->X311_ye[ii] - c_(1))*Xme[ii];
            
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

void sixdof_obj::netForces_cfd(lexer *p, fdm* a, ghostcell *pgc, double alpha)
{
    pnetinter->netForces_cfd(p,a,pgc,alpha,Xne,Yne,Zne,Kne,Mne,Nne);
    
    NETLOOP
    {
    // Add to external forces
        Xext += Xne[n];
        Yext += Yne[n];
        Zext += Zne[n];
        Kext += Kne[n];
        Mext += Mne[n];
        Next += Nne[n];
    }
}

void sixdof_obj::netForces_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha)
{
    pnetinter->netForces_nhflow(p,d,pgc,alpha,Xne,Yne,Zne,Kne,Mne,Nne);
    
    NETLOOP
    {
    // Add to external forces
        Xext += Xne[n];
        Yext += Yne[n];
        Zext += Zne[n];
        Kext += Kne[n];
        Mext += Mne[n];
        Next += Nne[n];
    }
}	

