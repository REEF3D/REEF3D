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
Authora: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"increment.h"
#include"slice4.h"

class lexer;
class fdm;
class ghostcell;

using namespace std;

#ifndef PARTRES2_H_
#define PARTRES2_H_


class partres2 : public increment
{
public:
        partres2(lexer *);
        ~partres2();
        
        
private:
    const int irand;
	const double drand;
    
    
    
};

#endif