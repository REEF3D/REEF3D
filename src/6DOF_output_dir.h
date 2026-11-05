/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef SIXDOF_OUTPUT_DIR_H_
#define SIXDOF_OUTPUT_DIR_H_

#include"lexer.h"

// Output folder of the 6DOF logs for the active hydrodynamic model:
// "./REEF3D_<model>_6DOF", mooring line geometry goes to the same name with "_Mooring".
inline const char* sixdof_output_dir(lexer *p)
{
    if(p->A10==3)
    return "./REEF3D_FNPF_6DOF";
    
    if(p->A10==5)
    return "./REEF3D_NHFLOW_6DOF";
    
    return "./REEF3D_CFD_6DOF";
}

#endif
