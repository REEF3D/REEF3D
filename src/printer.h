/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#ifndef PRINTER_H_
#define PRINTER_H_

class lexer;
class fdm;
class fdm_nhf;
class fdm_fnpf;
class ghostcell;
class turbulence;
class nhflow_turbulence;
class heat;
class ioflow;
class expdata;
class concentration;
class multiphase;
class sediment;

class printer
{
public:
    // CFD
    virtual void start(lexer*,fdm*,ghostcell*,turbulence*,heat*,ioflow*,expdata*,concentration*,multiphase*,sediment*){};
    virtual void print_stop(lexer*,fdm*,ghostcell*,turbulence*,heat*,ioflow*,expdata*,concentration*,multiphase*,sediment*){};
    // NHFLOW
    virtual void start(lexer*,fdm_nhf*,ghostcell*,ioflow*,nhflow_turbulence*,sediment*){};
    virtual void print_stop(lexer*,fdm_nhf*,ghostcell*,ioflow*,nhflow_turbulence*,sediment*){};
    // FNPF
    virtual void start(lexer*,fdm_fnpf*,ghostcell*,ioflow*){};
    virtual void print_stop(lexer*,fdm_fnpf*,ghostcell*){};
};

#endif
