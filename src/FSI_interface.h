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

#ifndef FSI_INTERFACE_H_
#define FSI_INTERFACE_H_

class lexer;
class fdm;
class fdm_nhf;
class fdm2D;
class field;
class vrans;
class ghostcell;
class turbulence;

using namespace std;

class FSI_interface
{
public:
	FSI_interface(lexer*,ghostcell*);
	virtual ~FSI_interface();
    
    // Generic
	void FSI_logic(lexer*,ghostcell*);
    void start(lexer*,ghostcell*);
    
	void initialize(lexer*,fdm*,ghostcell*,turbulence*) {};
    void forcing(lexer*,fdm*,ghostcell*,double,field&,field&,field&,field&,field&,field&,bool) {};
    
    
    // CFD
    void isource(lexer*,fdm*,vrans*) {};
    void jsource(lexer*,fdm*,vrans*) {};
    void ksource(lexer*,fdm*,vrans*) {};
    
    // NHFLOW
    void isource(lexer*,fdm_nhf*,vrans*) {};
    void jsource(lexer*,fdm_nhf*,vrans*) {};
    void ksource(lexer*,fdm_nhf*,vrans*) {};
    
    // SFLOW
    void isource(lexer*,fdm2D*,vrans*) {};
    void jsource(lexer*,fdm2D*,vrans*) {};
    void ksource(lexer*,fdm2D*,vrans*) {};
    
private:
};

#endif
