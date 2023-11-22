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

#include"onephase_ptf.h"
#include"ddweno_f_nug.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4.h"


using namespace std;

#ifndef ONEPHASE_PTF_F_H_
#define ONEPHASE_PTF_F_H_

class onephase_ptf_f : public onephase_ptf, public ddweno_f_nug
{
public:
    onephase_ptf_f(lexer*, fdm_ptf*, ghostcell*);
	virtual ~onephase_ptf_f();
    
	virtual void update(lexer*, fdm_ptf*, ghostcell*, ioflow*);
    virtual void ini(lexer*, fdm_ptf*, ghostcell*, ioflow*);
    
    field1 uf,urk1;
    field2 vf,vrk1;
    field3 wf,wrk1;
    
    field4 xphi,yphi,zphi;
    
private: 
    void fsf_update_ptf(lexer*,fdm_ptf*,ghostcell*);
    void disc_ptf(lexer*, fdm_ptf* , field&);
    int activenum;
    int gcval_u, gcval_v, gcval_w;
    
    double dx, dy, dz, dnorm, sign,deltax;
	double sx,sy,sz,snorm,op;
    double dt;
    double lsv,lsSig;
    double xmin,xplus,ymin,yplus,zmin,zplus;
    double nx,ny,nz;
    double nx0,ny0,nz0;

};

#endif
