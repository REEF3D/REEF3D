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

seedReturn partres::seeding(lexer *p, particles_obj &PP, size_t &index, double max, bool free)
{
        if(free)
        {
            cellSum[IJK] += PP.PackingFactor[index];
        }
        else
        {
            if(cellSumTopo[IJK]>=PP.PackingFactor[index])
                cellSumTopo[IJK] -= PP.PackingFactor[index];
            else if (cellSumTopo[IJK]>0)
            {
                PP.PackingFactor[index] = cellSumTopo[IJK];
                cellSumTopo[IJK] -= PP.PackingFactor[index];
            }
            else
            {
                return seedReturn::REMOVE;
            }
            cellSum[IJK] += PP.PackingFactor[index];
            
        }

        if(cellSum[IJK]>=max)
            return seedReturn::STOP;
        return seedReturn::CONTINUE;
}
