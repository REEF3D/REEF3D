/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the B117, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/liceonephases/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"onephase.h"
#include"ddweno_f_nug.h"

using namespace std;

#ifndef ONEPHASE_F_H_
#define ONEPHASE_F_H_

class onephase_f : public onephase, public ddweno_f_nug
{
public:
    onephase_f(lexer*, fdm*, ghostcell*);
	virtual ~onephase_f();
    
	virtual void update(lexer*, fdm*, ghostcell*, ioflow*);
    virtual void ini(lexer*, fdm*, ghostcell*, ioflow*);
    
    virtual void uvel(lexer*, fdm*, ghostcell*, field&);
    virtual void vvel(lexer*, fdm*, ghostcell*, field&);
    virtual void wvel(lexer*, fdm*, ghostcell*, field&);
    
private: 
    void fsf_update(lexer*,fdm*,ghostcell*);
    void disc(lexer *p, fdm* a, field& b);
    int activenum;
    
    double dx, dy, dz, dnorm, sign,deltax;
	double sx,sy,sz,snorm,op;
    double dt;
    double lsv,lsSig;
    double xmin,xplus,ymin,yplus,zmin,zplus;

};

#endif
