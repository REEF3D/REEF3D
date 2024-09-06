/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#include"partres.h"
#include"particles_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

    /// @brief Calculate complete intra-particle stress trensor
void partres::particleStressTensor(lexer *p, fdm &a, ghostcell &pgc, particles_obj &PP)
{
        double theta;
        int i,j,k;

        BLOOP
        {
            updateParticleStressTensor(p,a,PP,i,j,k);
        }
        pgc.start4V_par(p,stressTensor,10);
}

    /// @brief Calculate intra-particle stress trensor for cells around (`increment::i`,`increment::j`,`increment::k`)
void partres::particleStressTensorUpdateIJK(lexer *p, fdm &a, particles_obj &PP)
{
        double theta;
        int i,j,k;

        for (int n=-2; n<3; n++)
            for (int m=-2; m<3; m++)
                for (int l=-2; l<3; l++)
                {
                    i=increment::i+n;
                    j=increment::j+m;
                    k=increment::k+l;

                    updateParticleStressTensor(p,a,PP,i,j,k);
                }
}

    /// @brief Calculate intra-particle stress trensor for cell ( \p i , \p j , \p k )
void partres::updateParticleStressTensor(lexer *p, fdm &a, particles_obj &PP, int i, int j, int k)
{
        double theta=theta_s(p,a,PP,i,j,k);
        stressTensor[IJK]=Ps*pow(theta,beta)/max(theta_crit-theta,epsilon*(1.0-theta));
}