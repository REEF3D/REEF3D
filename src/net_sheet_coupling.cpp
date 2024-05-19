/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2024 Tobias Martin

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

#include"net_sheet.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"


void net_sheet::vransCoupling(lexer *p, fdm *a, ghostcell *pgc)
{
    double rho;

    //- Save Lagrangian coordinates and forces
    for (int knotI = 0; knotI < nK; knotI++)
    {
        lagrangePoints[knotI] << x_.row(knotI).transpose();
        
        // Forces
        const Eigen::Vector3d& coordI = lagrangePoints[knotI];
            
        // Density 
        rho = coupledField[knotI][3];
        
        if 
        (
            coordI(0) >= xstart[p->mpirank] && coordI(0) < xend[p->mpirank] &&
            coordI(1) >= ystart[p->mpirank] && coordI(1) < yend[p->mpirank] &&
            coordI(2) >= zstart[p->mpirank] && coordI(2) < zend[p->mpirank] &&
            rho > 900.0
        )
        {
            // Divide by rho since multiplied in vrans
            lagrangeForces[knotI] = forces_knot.row(knotI)/rho; 
        }
        else
        {
            lagrangeForces[knotI] << 0.0, 0.0, 0.0;   
        }
    }   
    
    for (int pI = 0; pI < nK; pI++)
    {
        Eigen::Vector3d& forceI = lagrangeForces[pI];
        forceI << pgc->globalsum(forceI(0)), pgc->globalsum(forceI(1)), pgc->globalsum(forceI(2));
    }

}


