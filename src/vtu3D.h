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

#ifndef VTU3D_H_
#define VTU3D_H_

#include "vtk3D.h"
#include "increment.h"

class vtu3D : public vtk3D , increment
{
    public:
        vtu3D();
        ~vtu3D();
    public:
        void folder(const char*);
        void offset(lexer*, int*, int&);
        void structureWrite(lexer*, fdm*, std::ofstream&);
        void structureWrite(lexer*, fdm_fnpf*, std::ofstream&);
        void structureWrite(lexer*, fdm_nhf*, std::ofstream&);
        void beginning(lexer*, std::ofstream&);
        void beginningParallel(lexer*, std::ofstream&);
        void ending(std::ofstream&, const int*, int&);
        void endingParallel(std::ofstream&, const char*, const int, const int);
        void fileName(char *name, const char *A10, const int num, const int rank){sprintf(name,"./REEF3D_%s_VTU/REEF3D-%s-%08i-%06i.vtu",A10,A10,num,rank);};
        void parallelFileName(char *name, const char *A10, const int num){sprintf(name,"./REEF3D_%s_VTU/REEF3D-%s-%08i.pvtu",A10,A10,num);};
    private:
        void structureWriteEnd(lexer*, std::ofstream&);
        char pname[50];
};

#endif