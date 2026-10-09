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

#ifndef PATCHBC_CORE_H_
#define PATCHBC_CORE_H_

#include"patchBC_interface.h"

using namespace std;

// model independent part of the patch boundary conditions (CFD: patchBC, SFLOW: patchBC_2D):
// patches from the control input, kind and validation, face selection by geometry, hydrographs
class patchBC_core : public patchBC_interface, public increment
{
public:
	patchBC_core() = default;
	virtual ~patchBC_core() = default;
    
protected:
    // creates the patch objects, reads and checks the B 41x/42x settings; aborts on errors
    void patch_setup(lexer*, ghostcell*);
    
    // patch index of a face (cell ii,jj,kk, side cs), -1 if in no patch geometry; twoD: no z test
    int patch_find(lexer*, int ii, int jj, int kk, int cs, bool twoD);
    
    // after the face selection: global face counts, face dependent checks, summary output
    void patch_check(lexer*, ghostcell*);
    
    // allocates patch[qq]->gcb from the selected faces (i,j,k,cs,index)
    void patch_faces(lexer*, int qq, const int *list, int count);
    
    // hydrograph values at the current time
    void patch_hydrograph(lexer*);
    
    // wetted fraction of a side face from the level set of the cell, 0..1
    static double wetfrac(double phival, double dz);
    
    void patch_error(lexer*, ghostcell*, const char*, int ID);
    
private:
    void hydrograph_read(lexer*, ghostcell*, const char*, int ID, double **&, int &);
    double hydrograph_ipol(lexer*, double**, int);
    
    void face_center(lexer*, int ii, int jj, int kk, int cs, bool twoD, double&, double&, double&);
    int patch_index(int ID);
    
    int error_count=0;
};

#endif
