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

#include "sediment_part.h"

/// @brief SFLOW calculation function
void sediment_part::start_sflow(lexer *p, fdm2D *b, ghostcell *pgc, ioflow* pflow, slice &P, slice &Q)
{}

/// @brief SFLOW initialization function
void sediment_part::ini_sflow(lexer *p, fdm2D *b, ghostcell *pgc)
{}

/// @brief Updates the topography for the SFLOW solver
void sediment_part::update_sflow(lexer *p, fdm2D *b, ghostcell *pgc, ioflow *pflow)
{}