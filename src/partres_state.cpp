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


/// @brief Writes the state of the partres class to file.
/// @ToDo Write cellSumTopo 
void partres::writeState(lexer *p, ofstream &result)
{
        float ffn;
        PLAINLOOP
        {
            ffn=cellSum[IJK];
            result.write((char*)&ffn, sizeof (float));
        }
        result.write((char*)&ffn, sizeof (float));  
}

/// Reads the state of the partres class from file.
/// Reconstructs cellSum, columnSum and stressTensor
void partres::readState(lexer *p, ifstream &result)
{
        float ffn;
        PLAINLOOP
        {
            result.read((char*)&ffn, sizeof (float));
            cellSum[IJK]=double(ffn);
        }

}