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

#include "partres.h"
#include "particles_obj.h"
#include "lexer.h"
#include "fdm.h"
#include "ghostcell.h"
#include "field4a.h"
#include "sediment_fdm.h"
#include "turbulence.h"

    /**
     * @brief Sets up the partres class.
     *
     * This function is responsible for setting up the partres class by calculating
     * the cellSumTopo and columnSum values based on the given parameters.
     *
     * @param p A pointer to the lexer object.
     * @param a A reference to the fdm object.
     * @param diameter The diameter value.
     */

void partres::setup(lexer *p, fdm &a, double &diameter)
{
        BLOOP
        {
            cellSumTopo[IJK] = maxParticlesPerCell(p,a,diameter);
            columnSum[IJ] += cellSumTopo[IJK];
        }
}

void partres::setupState(lexer *p, fdm &a, ghostcell &pgc, particles_obj &PP)
{
        // particlePerCell(p,pgc,PP);
        // PLAINLOOP
        // columnSum[IJ] += cellSumTopo[IJK]+cellSum[IJK];
        // particleStressTensor(p,a,pgc,PP);
}

    /**
     * @brief Determines if a particle should be seeded in the current cell.
     *
     * This function is responsible for determining if a particle should be seeded in the current cell.
     * It does so by checking if the cellSumTopo value is greater than 0 and if the cellSum value is less than the maximum value.
     *
     * @param p A pointer to the lexer object.
     * @param PP A reference to the particles_obj object.
     * @param index The index of the particle.
     * @param max The maximum value.
     * @return True if the particle should be seeded, false otherwise.
     */