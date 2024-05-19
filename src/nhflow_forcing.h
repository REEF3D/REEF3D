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

#include"increment.h"

class lexer;
class fdm_nhf;
class ghostcell;
class nhflow_reinidisc_fsf;
#include<vector>

using namespace std;

#ifndef NHFLOW_FORCING_H_
#define NHFLOW_FORCING_H_

class nhflow_forcing : public increment
{
public:
	nhflow_forcing(lexer*);
	virtual ~nhflow_forcing();
    
    void forcing(lexer*, fdm_nhf*, ghostcell*, double, double*, double*, double*, double*, double*, double*);

    void forcing_ini(lexer*, fdm_nhf*, ghostcell*);
    
    void ray_cast(lexer*, fdm_nhf*, ghostcell*);
    void ray_cast_io(lexer*, fdm_nhf*, ghostcell*,int,int);
    void ray_cast_direct(lexer*, fdm_nhf*, ghostcell*,int,int);
    
    void objects_create(lexer*, ghostcell*);
    void objects_allocate(lexer*, ghostcell*);
    
    
    void reini_RK2(lexer*, fdm_nhf*, ghostcell*, double*);
    
private:
    void box(lexer*, ghostcell*, int);
    void cylinder_z(lexer*, ghostcell*, int);
    
    void geometry_refinement(lexer*, ghostcell*);
    void create_triangle(double&,double&,double&,double&,double&,double&,double&,double&,double&,const double&,const double&,const double&);
    
    int *IO,*CR,*CL;
    double *FRK1,*dt,*L;
    
    double **tri_x,**tri_y,**tri_z,**tri_x0,**tri_y0,**tri_z0;
    vector<vector<double> > tri_x_r;
	vector<vector<double> > tri_y_r;
	vector<vector<double> > tri_z_r;
    
    int *tstart,*tend;
    int tricount;
    int entity_count, entity_sum;
    double xs,xe,ys,ye,zs,ze;
    double zmin,zmax;
    
    int reiniter;
    
    const double epsi;
    
    nhflow_reinidisc_fsf *prdisc;


 
};

#endif