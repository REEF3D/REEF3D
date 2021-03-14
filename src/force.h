/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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

#include"boundarycheck.h"
#include"wave_interface.h"
#include"field1.h"
#include"field2.h"
#include<iostream>
#include<fstream>

class lexer;
class fdm;
class ghostcell;


using namespace std;

#ifndef FORCE_H_
#define FORCE_H_

class force : public boundarycheck, public wave_interface
{

public:
	force(lexer*,fdm*,ghostcell*,int);
	virtual ~force();
	virtual void start(lexer*,fdm*,ghostcell*);
	virtual void ini(lexer*,fdm*,ghostcell*);
	virtual void force_surface(lexer*,fdm*,ghostcell*);

	virtual void velocity(lexer*,fdm*,ghostcell*);
	virtual void cellnodes(lexer*,fdm*,ghostcell*);
	virtual void surfarea(lexer*,fdm*,ghostcell*);
	virtual void cellnodes_gcb(lexer*,fdm*,ghostcell*);
	virtual void surfarea_gcb(lexer*,fdm*,ghostcell*);
	virtual double area(int,int,int);
    virtual void print(lexer*,fdm*,ghostcell*);


private:
	char name[100];
	
	double xs,xe,ys,ye,zs,ze;
	int is,ie,js,je,ks,ke;
	double x_center,y_center,z_center;
	double FDS1,FDS2;
	double Fvert;
	double FD,FL;
	double FDs, FLs;
	double FDs_norm,Fvert_norm,Fvert0;
	double FTs,FTb;
	double FD_Press, FD_Shear;
	double F_morison,F_morison_rect,F_morison_norm;
	double Fmi, Fmd;
	double cd,cm;
	double cD,cL,cDs,cLs;
	double Re,Kc;
	double um,vm,un1,un2,un3,vn1,vn2,vn3,dt1,dt2,dt3,udt,vdt;
	double uinf,vinf,uinf1,vinf1,uinf2,vinf2,uinf_dt, vinf_dt;
	double vt1,vt2;
	double ut1,ut2;
	double udt1,udt2;
	double xm,ym,zm;
	int im,jm,km;
    ofstream fout;
    double cx,cy,cz;
    double px,py,pz;
    int n,q;
    int count,fnum,*fcheck;
    double ***ccpt,*farea, **fn,*fx,*fy,*fz;
    int **fid, *surfnum;
    int pt[16];
    int cl[16];
    int nd[16];
	double epsi;
	double H,phival,fsfmargin;
};

#endif



