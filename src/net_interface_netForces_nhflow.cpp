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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"net_interface.h"
#include"lexer.h"
#include"ghostcell.h"
#include<sys/stat.h>

#include"net.h"
#include"net_void.h"
#include"net_barDyn.h"
#include"net_barQuasiStatic.h"
#include"net_sheet.h"

void net_interface::netForces_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha, Eigen::Matrix3d quatRotMat,
                            vector<double> Xne, vector<double> Yne, vector<double> Zne, vector<double> Kne, vector<double> Mne, vector<double> Nne)
{
    NETLOOP
    {
        pnet[n]->start_nhflow(p, d, pgc, alpha, quatRotMat);
        dlm_nhflow(p, d, pgc, n);
    
        // Forces on rigid body
        pnet[n]->netForces(p,Xne[n],Yne[n],Zne[n],Kne[n],Mne[n],Nne[n]);
        
        // Distribute forces and moments to all processors
        pgc->bcast_double(&Xne[n],1);
        pgc->bcast_double(&Yne[n],1);
        pgc->bcast_double(&Zne[n],1);
        pgc->bcast_double(&Kne[n],1);
        pgc->bcast_double(&Mne[n],1);
        pgc->bcast_double(&Nne[n],1);	
    }
}