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

#include"lexer.h"

lexer::lexer() : cmu(0.09), position(this), interpolation(this), coordinates(this)
{
    sigT=0.9;
    
	control::ini_default();

    solveriter=0;
	mpirank=0;

	simtime=0.0;
	poissontime=0.0;
	pressval=0;
    alpha=0.0;
    solidread=toporead=porousread=0;
    net_count=0;
    mooring_count=0;
}

lexer::~lexer()
{
}
