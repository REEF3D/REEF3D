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

#ifndef vtk3D_H_
#define vtk3D_H_

#include <fstream>

class lexer;
class fdm;
class ghostcell;

class vtk3D
{
    public:
        virtual void folder()=0;
        virtual void offset(lexer*,int*,int&)=0;
        virtual void structureWrite(lexer*, fdm*, std::ofstream&)=0;
        virtual void extent(lexer*,ghostcell*){};
        virtual void fileName(char*, int&, int&)=0;
        virtual void parallelFileName(char*, int&)=0;
        virtual void beginning(lexer*, std::ofstream&)=0;
        virtual void beginningParallel(lexer*, std::ofstream&)=0;
        virtual void ending(std::ofstream&, int*, int&)=0;
        virtual void endingParallel(std::ofstream&, int&, int&)=0;
    protected:
        void vtkVersion(std::ofstream &result){result<<"version=\"0.1\" byte_order=\"LittleEndian\">"<<std::endl;}; // header_type=\"UInt32\"
        void xmlVersion(std::ofstream &result){result<<"<?xml version=\"1.0\"?>"<<std::endl;};
};

#endif