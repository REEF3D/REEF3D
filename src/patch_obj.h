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

#ifndef PATCH_OBJ_H_
#define PATCH_OBJ_H_

class lexer;

using namespace std;

// one patch boundary (all B 440-442 shapes with the same ID)
class patch_obj
{
public:
	patch_obj(lexer*,int);
	virtual ~patch_obj();
    
    // velocity vector (u,v,w) with the component normal to the face (+axis direction) Un,
    // along the patch flow direction (B 415, B 416, B 417; default: normal to the face)
    void velocity(int cs, double Un, double &uval, double &vval, double &wval) const;
    
    int ID;
    int kind;            // PATCH_INLET or PATCH_OUTLET
    int gcb_flag;        // boundary code of the faces: PATCH_INLET(_FSF) / PATCH_OUTLET(_FSF)
    
    // faces: i,j,k,cs and the index in gcb4 (CFD) / gcbsl4 (SFLOW)
    int gcb_count;
    int **gcb;
    
    // inlet: normal velocity from discharge (B 411, B 421) or given (B 414)
    int Q_flag;
    double Q, Uq;
    
    int Uio_flag;
    double Uio;
    
    // inlet: velocity components (B 415)
    int velcomp_flag;
    double U,V,W;
    
    // inlet: flow direction, B 416 horizontal angle from the face normal (rad),
    // B 417 direction vector (unit vector)
    int angle_flag;
    double alpha;
    
    int dir_flag;
    double dirx,diry,dirz;
    
    // outlet: pressure (B 412), SFLOW free stream (B 418)
    int pressure_flag;
    double pressure;
    
    int pio_flag;
    
    // water level (B 413, B 422), inlet or outlet
    int waterlevel_flag;
    double waterlevel;
    
    // hydrographs
    int hydroQ_flag;
    double **hydroQ;
    int hydroQ_count;
    
    int hydroFSF_flag;
    double **hydroFSF;
    int hydroFSF_count;
    
    // measured at the patch: discharge, mean normal velocity, wetted area, water level
    double Q0,U0,A0,h0;
};

#endif
