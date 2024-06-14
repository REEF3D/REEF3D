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

#ifndef NHFLOW_PRINTER_H_
#define NHFLOW_PRINTER_H_

class lexer;
class fdm_nhf;
class ghostcell;
class ioflow;

class nhflow_printer
{
public:
	virtual void start(lexer*,fdm_nhf*,ghostcell*,ioflow*)=0;
    virtual void print_vtk(lexer*,fdm_nhf*,ghostcell*)=0;
    virtual void print_stop(lexer*,fdm_nhf*,ghostcell*,ioflow*)=0;
};

#endif
