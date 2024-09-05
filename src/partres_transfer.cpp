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
     * @brief Transfers the particle to the current cell.
     *
     * This function is responsible for transferring the particle to the current cell.
     * It does so by adding the ParcelFactor value to the cellSum value.
     *
     * @param p A pointer to the lexer object.
     * @param PP A reference to the particles_obj object.
     * @param index The index of the particle.
     */
void partres::transfer(lexer *p, particles_obj &PP, size_t &index)
{
        cellSum[IJK] += PP.ParcelFactor[index];
    
}

    /**
     * @brief Removes the particle from the current cell.
     *
     * This function is responsible for removing the particle from the current cell.
     * It does so by subtracting the ParcelFactor value from the cellSum value.
     *
     * @param p A pointer to the lexer object.
     * @param PP A reference to the particles_obj object.
     * @param index The index of the particle.
     */
void partres::remove(lexer *p, particles_obj &PP, size_t &index)
{
        cellSum[IJK] -= PP.ParcelFactor[index];
}