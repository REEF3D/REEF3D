/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

void net_interface::netForces_cfd(lexer *p, fdm *a, ghostcell *pgc, double alpha, Eigen::Matrix3d quatRotMat, 
                            vector<double> &Xne, vector<double> &Yne, vector<double> &Zne, vector<double> &Kne, 
                            vector<double> &Mne, vector<double> &Nne, bool finalize)
{
    NETLOOP
    {
        pnet[n]->start_cfd(p, a, pgc, alpha, quatRotMat, finalize);
        
        // the forces of the net onto the faces; the ghost cells of all nets after the loop
        dlm_cfd_spread(p, a, n);
        
        // Forces on rigid body
        pnet[n]->netForces(p,Xne[n],Yne[n],Zne[n],Kne[n],Mne[n],Nne[n]);
        
        // Distribute forces and moments to all processor
        pgc->bcast_double(&Xne[n],1);
        pgc->bcast_double(&Yne[n],1);
        pgc->bcast_double(&Zne[n],1);
        pgc->bcast_double(&Kne[n],1);
        pgc->bcast_double(&Mne[n],1);
        pgc->bcast_double(&Nne[n],1);
    }
    
    // one exchange of the net forces for all nets (dlm_cfd: three for each net)
    if(p->net_count>0)
    {
    pgc->start1(p,a->Fext,10); 
    pgc->start2(p,a->Gext,11); 
    pgc->start3(p,a->Hext,12); 
    }
}