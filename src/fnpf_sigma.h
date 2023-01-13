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

#include"fnpf.h"
#include"increment.h"
#include"slice4.h"

class lexer;
class fdm_fnpf;
class ghostcell;
class fnpf_fsf;

using namespace std;

#ifndef FNPF_SIGMA_H_
#define FNPF_SIGMA_H_

class fnpf_sigma : public increment
{
public:
	fnpf_sigma(lexer*, fdm_fnpf*, ghostcell*);
	virtual ~fnpf_sigma();
    
    virtual void sigma_ini(lexer*, fdm_fnpf*, ghostcell*, fnpf_fsf*, slice&);
    virtual void sigma_update(lexer*, fdm_fnpf*, ghostcell*, fnpf_fsf*, slice&);
        
private:

};

#endif
