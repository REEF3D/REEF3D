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

#include"ddweno_f_nug.h"
#include"slice4.h"
class lexer;
class slice;
class sliceint;
class ghostcell;

#ifndef NHFLOW_COASTLINE_H_
#define NHFLOW_COASTLINE_H_

using namespace std;

class nhflow_coastline :  public ddweno_f_nug
{
public:
    nhflow_coastline(lexer*);
	virtual ~nhflow_coastline();

   void start(lexer*,ghostcell*,slice&,int*,sliceint&);
   
private:
   void reini(lexer*,ghostcell*,slice&);
   void disc(lexer*,ghostcell*,slice&);
   
   void step(lexer*);
   void time_preproc(lexer*);
   
   slice4 frk1,frk2,L,dt,wet_n;
   
   
   int reiniter,change;
   

};

#endif
