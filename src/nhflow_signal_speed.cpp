/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"nhflow_signal_speed.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm_nhf.h"
#include"patchBC_interface.h"

nhflow_signal_speed::nhflow_signal_speed(lexer* p) 
{

}

nhflow_signal_speed::~nhflow_signal_speed()
{
}

void nhflow_signal_speed::wave_speed_update(lexer* p, ghostcell *pgc, fdm_nhf*, 
                                        double *Us, double *Un, double *Ve, double *Vw, 
                                        slice &Ds,slice &Dn, slice &De, slice &Dw)
{

    
}



