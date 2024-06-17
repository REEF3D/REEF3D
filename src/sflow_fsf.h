/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#ifndef SFLOW_FSF_H_
#define SFLOW_FSF_H_

class lexer;
class fdm2D;
class ghostcell;
class ioflow;
class slice;

using namespace std;

class sflow_fsf
{
public:    
    virtual void start(lexer*, fdm2D*, ghostcell*, ioflow*,slice&,slice&,double)=0;
	virtual void ini(lexer*, fdm2D*, ghostcell*, ioflow*)=0;
	virtual void depth_update(lexer*, fdm2D*, ghostcell*,slice&,slice&,slice&,slice&)=0;
    virtual void breaking(lexer*, fdm2D*, ghostcell*, slice&, slice&, double)=0;
    virtual void breaking_persist(lexer*, fdm2D*, ghostcell*, slice&, slice&, double)=0;
    virtual void wetdry(lexer*, fdm2D*, ghostcell*,slice&,slice&,slice&)=0;
    
        

};

#endif
