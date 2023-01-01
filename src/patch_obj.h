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

#include"increment.h"

class lexer;
class fdm;
class ghostcell;
class field;

using namespace std;

#ifndef PATCH_OBJ_H_
#define PATCH_OBJ_H_

class patch_obj : public increment
{
public:
	patch_obj(lexer*,int);
	virtual ~patch_obj();
    
    void patch_obj_ini(lexer *p, ghostcell *pgc);
    
    void patch_obj_gcb_generate(lexer *p, ghostcell *pgc);
    
    // Patch DATA 3D
    int ID;
    int IO;
    int gcb_count;
    int **gcb;
    int gcb_flag;
    int gcb_uflag, gcb_pressflag, gcb_phiflag;
    int counter;
    
    /*
    B211=0;        // int patchBC discharge
    B212=0;        // int patchBC pressure BC
    B213=0;        // int patchBC waterlevel
    B214=0;        // int patchBC perpendicular velocity
    B215=0;        // int patchBC velocity components
    B216=0;        // int patchBC horizontal inflow angle
    B217=0;        // int patchBC inflow normals
    */
    
    int Q_flag;
    double Q, Uq;
    
    int velocity_flag;
    double velocity;
    
    int pressure_flag;
    double pressure;
    
    int waterlevel_flag;
    double waterlevel;
    
    int Uio_flag;
    double Uio;
    
    int velcomp_flag;
    double U,V,W;
    
    int flowangle_flag;
    double alpha;
    double sinalpha,cosalpha;
    
    int flownormal_flag;
    double Nx,Ny,Nz;
    
    int pio_flag;
    
    int hydroQ_flag;
    double **hydroQ;
    int hydroQ_count,hydroQ_iter;
    
    int hydroFSF_flag;
    double **hydroFSF;
    int hydroFSF_count,hydroFSF_iter;
    
    // measurement
    double Q0,U0,A0,h0;


};

#endif
