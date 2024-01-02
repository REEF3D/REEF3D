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
Author: Hans Bihs
--------------------------------------------------------------------*/

class fdm;
class lexer;
class convection;
class solver;
class ghostcell;
class ioflow;
class reini;
class particle_corr;
class field;

using namespace std;

#ifndef FREESURFACE_H_
#define FREESURFACE_H_

class freesurface
{
public:

	virtual void start(fdm*,lexer*, convection*, solver*, ghostcell*,ioflow*, reini*, particle_corr*,field&)=0;
    virtual void update(lexer*,fdm*,ghostcell*,field&)=0;

};

#endif
