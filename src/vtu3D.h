/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#ifndef VTU3D_H_
#define VTU3D_H_

#include "vtk3D.h"
#include "increment.h"

class vtu3D : public vtk3D , increment
{
    public:
        void folder(const char*) override;
        void fileName(char*, const unsigned int, const char*, const int, const int) override;
        void parallelFileName(char*, const unsigned int, const char*, const int) override;

        void offset(lexer*, int*, int&) override;

        void beginning(lexer*, std::ostream&) override;
        void beginningParallel(lexer*, std::ostream&) override;
        void ending(std::ostream&, const int*, int&) override;
        void endingParallel(std::ostream&, const char*, const int, const int) override;

        void structureWrite(lexer*, fdm*, std::vector<char>&, size_t&) override;
        void structureWrite(lexer*, fdm_fnpf*, std::vector<char>&, size_t&) override;
        void structureWrite(lexer*, fdm_nhf*, std::vector<char>&, size_t&) override;
    private:
        void structureWriteEnd(lexer*, std::vector<char>&, size_t&);
        char pname[50];
};

#endif
