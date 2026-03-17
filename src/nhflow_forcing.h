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

#ifndef NHFLOW_FORCING_H_
#define NHFLOW_FORCING_H_

#include"nhflow_geometry.h"
#include"slice4.h"
#include"vtp3D.h"

class lexer;
class fdm_nhf;
class ghostcell;
class slice;
class sixdof;
class mooring;
class fsi;
class nhflow_reinidisc_fsf;

using namespace std;

class nhflow_forcing : public nhflow_geometry, private vtp3D
{
public:
	nhflow_forcing(lexer*, fdm_nhf*, ghostcell*);
	virtual ~nhflow_forcing();
    
    void forcing(lexer*, fdm_nhf*, ghostcell*, sixdof *p6dof, 
                 int, double, double*, double*, double*, slice&, bool);
    
    void solid_forcing(lexer*, fdm_nhf*, ghostcell*, double, double*, double*, double*, slice&);
    void forcing_ini(lexer*, fdm_nhf*, ghostcell*);
        
    void reset(lexer*, fdm_nhf*, ghostcell*);
    
    double Hsolidface(lexer*, fdm_nhf*, int, int, int);
    
    // DLM
    void dlm_forcing(lexer*, fdm_nhf*, ghostcell*, double, double*, double*, double*, slice&);
    void dlm_forcecalc(lexer*, fdm_nhf*, ghostcell*, double, double*, double*, double*, slice&);
    void dlm_forcing_ini(lexer*, ghostcell*);
    double kernel(const double&);
    
private:
    double *FX,*FY,*FZ;
    slice4 fe;
    double starttime;
    
        
    double H,Ht, uf, vf, wf, ef;
    double efc;
	double nx, ny, nz,norm ;
	double phival_sf;
    double dirac;
    
    
    // DLM
    double *EL_L,*EL_dx;
    double **EL_X,**EL_Y,**EL_Z,**EL_V;
    double **EL_FX,**EL_FY,**EL_FZ;
    int **EL_f;
    
    int Ne,Np;
    int ii,jj,kk;
    int gcval_eta;
    double dx,dy,dz;
    double D,dist;
    
    int gcval_u, gcval_v, gcval_w;
    int gcval_uh, gcval_vh, gcval_wh;
};

#endif
