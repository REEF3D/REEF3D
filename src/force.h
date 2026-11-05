/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#ifndef FORCE_H_
#define FORCE_H_

#include"fieldint5.h"
#include"field5.h"
#include"increment.h"
#include"vtp3D.h"
#include<iostream>
#include<fstream>

class lexer;
class fdm;
class ghostcell;
class amr_grid_lookup;

using namespace std;

class force :  public increment, private vtp3D
{

public:
	force(lexer*,fdm*,ghostcell*,int,bool=true);    // false: no output files (a patch of the mesh refinement)
	virtual ~force();
	void start(lexer*,fdm*,ghostcell*);
    void ini(lexer*,fdm*,ghostcell*,bool=true);   // false: the surface only (restart from a checkpoint)

    // CFD mesh refinement (cfd_amr_force.cpp): every piece of the solid surface is integrated on the
    // finest grid at its centre.  Level 0 (amr_level0): its pieces under a patch are left out, the
    // patches of this rank add theirs (amr_grid_lookup::solid_forces) before the sums over the ranks.
    // A patch (amr_patch, its level and the box of its interior): its pieces in the box on the finest
    // grid there, summed into F (Fx, Fy, Fz, area) by patch_sum (rank-local)
    void amr_level0(const amr_grid_lookup*);
    void amr_patch(const amr_grid_lookup*, int lev, const double *box);
    void patch_sum(lexer*,fdm*,double *F);

private:
	void triangulation(lexer*, fdm*, ghostcell*, field&);
	void reconstruct(lexer*, fdm*, field&);
	void addpoint(lexer*,fdm*,int,int);
	
	int **tri, **facet, *confac, *numfac,*numpt;
	double **ccpt, **pt, *ls;
	double   dV1,dV2,C1,C2,mi;
	int numtri,numvert, numtri_mem, numvert_mem;
	int count,countM,n,nn,q;
	int ccptcount,facount,check;
	int polygon_sum,polygon_num,vertice_num;
	const double zero,interfac;
    double epsi;
	

	fieldint5 vertice, nodeflag;
    field5 eta;
	
	
    void force_calc(lexer*,fdm*,ghostcell*);
    void force_sum(lexer*,fdm*,double *F);      // this rank's pieces: Fx, Fy, Fz, area; Ax, Ay, Px
    
    const amr_grid_lookup *amr = nullptr;
    int amr_lev = 0;
    double amr_box[6];
    bool amr_own(double,double,double) const;
    void free_arrays();
    
	void print_force(lexer*,fdm*,ghostcell*);
    void print_ini(lexer*,fdm*,ghostcell*);
    void print_vtp(lexer*,fdm*,ghostcell*);
    void pvtp(lexer*,int);

    char name[100];
    int iin,offset[100];
    float ffn;
    int gcval_phi;
    double printtime,printtime2;
	int forceprintcount;
    int gcval_press;
    
    // force
    double Fx,Fy,Fz;
    double A_tot,A;
    
    double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;
    double xc,yc,zc;
    double nx,ny,nz,norm;
    double nxs,nys,nzs;
    double uval,vval,wval,pval,viscosity,density,phival;
    double du,dv,dw;
    double at,bt,ct,st;
    
    ofstream fout;
    
    double xs,xe,ys,ye,zs,ze;
    double xm,ym,zm;
	int is,ie,js,je,ks,ke;
    const int ID;
	

};

#endif


