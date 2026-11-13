/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2026 Tobias Martin

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

#include"net_sheet.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include<vector>

void net_sheet::coupling_dlm_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
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
    
    // the force of a point is on the rank that holds it, zero on the others: one reduction for all
    // points (it was three for each point: some 13 000 per stage for the 25 nets of OF1)
    std::vector<double> fsum(3*nK);
    for (int pI = 0; pI < nK; pI++)
    for (int q = 0; q < 3; q++)
    fsum[3*pI+q] = lagrangeForces[pI](q);
    
    pgc->globalsum(fsum.data(),3*nK);
    
    for (int pI = 0; pI < nK; pI++)
    lagrangeForces[pI] << fsum[3*pI], fsum[3*pI+1], fsum[3*pI+2];

}


