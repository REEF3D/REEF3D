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
Authors: Tobias, Hans Bihs
--------------------------------------------------------------------*/

#include"net_interface.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<sys/stat.h>

#include"net.h"
#include"net_void.h"
#include"net_barDyn.h"
#include"net_barQuasiStatic.h"
#include"net_sheet.h"

void net_ionterface::net_logic(lexer *p, ghostcell *pgc, vector<net*>& pnet)
{

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
        mkdir("./REEF3D_CFD_6DOF_Net",0777);	

        else
        p->X320_type = new int[p->net_count];
		
        pnet.reserve(p->net_count);	
  
		for (int n=0; n < p->net_count; n++)
		{
            MPI_Bcast(&p->X320_type[n],1,MPI_INT,0,pgc->mpi_comm);
			
            if(p->X320_type[n] > 10)
			{
				pnet.push_back(new net_barDyn(n,p));
			}
            else if (p->X320_type[n] < 4)
            {
                pnet.push_back(new net_barQuasiStatic(n,p));
            }
            else
            {
                 pnet.push_back(new net_sheet(n,p));
            }
			
            pnet[n]->initialize(p,a,pgc);
		}
    }
