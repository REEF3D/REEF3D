/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2025 Tobias Martin

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
Authors: Tobias Martin, Hans Bihs
--------------------------------------------------------------------*/

#include"net_barDyn.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void net_barDyn::updateAcc(lexer *p, ghostcell *pgc)
{
    Vector3d T_knot, x_ij;
    int knotJ, barI;
    double l_ij;

    // Move top as rigid body
    updateTopAcc(p);

    // Calculate acceleration from tension forces
    for (int knotI = 0; knotI < nK; knotI++)
    {
        T_knot << 0.0, 0.0, 0.0;
        
        if (knotI >= nfK[0][0]) // then inner knot
        {
            xdotdot_.row(knotI) *= 0.0; 
            
            for (int k = 1; k < 5; k++)
            {
                barI = nfK[knotI - nfK[0][0]][k];

                if (barI!=-1)
                {
                    // Find bar vector
                    if (Pi[barI]==knotI)
                    {
                        knotJ = Ni[barI];
                    }
                    else
                    {
                        knotJ = Pi[barI];
                    }

                    x_ij = x_.row(knotJ) - x_.row(knotI);

                    l_ij = x_ij.norm();

                    T_knot += T_(barI)*x_ij/l_ij; 
                }
            }

            xdotdot_.row(knotI) = 
                (forces_knot.row(knotI) + T_knot)
                /
                (mass_knot(knotI) + added_mass(knotI));
        }
        else    // rigid motion of top knots
        {
            xdotdot_.row(knotI) = top_xdotdot_.row(knotI);
        }
    } 
}


void net_barDyn::updateTopAcc(lexer *p)
{
    // Update top knot velocities
    for (int i = 0; i < nbK; i++)
    {
        top_xdot_(i,0) = p->ufbi + (x_(i,2) - p->zg)*p->qfbi - (x_(i,1) - p->yg)*p->rfbi;
        top_xdot_(i,1) = p->vfbi + (x_(i,0) - p->xg)*p->rfbi - (x_(i,2) - p->zg)*p->pfbi;
        top_xdot_(i,2) = p->wfbi + (x_(i,1) - p->yg)*p->pfbi - (x_(i,0) - p->xg)*p->qfbi;
    }

    // Calculate top knot acceleration
    top_xdotdot_ = 
          coeffs_(0)*top_xdot_ + coeffs_(1)*xdot_.block(0,0,nbK,3)
        + coeffs_(2)*xdotn_.block(0,0,nbK,3) + coeffs_(3)*xdotnn_.block(0,0,nbK,3);
}



