/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef VTR3D_H_
#define VTR3D_H_

#include "vtk3D.h"
#include "increment.h"

class vtr3D : public vtk3D , increment
{
    public:
        vtr3D();
        ~vtr3D();
    public:
        void folder(char*);
        void offset(lexer*, int*, int&);
        void structureWrite(lexer*, fdm*, std::ofstream&);
        void extent(lexer* ,ghostcell*);
        
        void beginning(lexer*, std::ofstream&);
        void beginningParallel(lexer*, std::ofstream&);
        void ending(std::ofstream&, int*, int&);
        void endingParallel(std::ofstream&, char*, int&, int&);
        void fileName(char *name, char *A10, int &num, int &rank){sprintf(name,"./REEF3D_%s_VTR/REEF3D-%s-%08i-%06i.vtr",A10,A10,num,rank);};
        void parallelFileName(char *name, char *A10, int &num){sprintf(name,"./REEF3D_%s_VTR/REEF3D-%s-%08i.pvtr",A10,A10,num);};
    private:
        int *piextent;
        char pname[50];
        char pextent[200];
};

#endif