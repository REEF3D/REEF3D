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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#ifndef VTK3D_H_
#define VTK3D_H_

#include <ostream>
#include <sstream>
#include <vector>
#include <string>
#include <cstring>
#include <iomanip>

class lexer;
class fdm;
class fdm_fnpf;
class fdm_nhf;
class ghostcell;

class vtk3D
{
    public:
        virtual ~vtk3D() = default;

        virtual void folder(const char*) {};
        virtual void fileName(char*, const unsigned int, const char*, const int, const int) {};
        virtual void parallelFileName(char*, const unsigned int, const char*, const int) {};
        virtual void extent(lexer*,ghostcell*) {};

        virtual void offset(lexer*,int*,int&) {};
        virtual void beginning(lexer*, std::ostream&) {};
        virtual void beginning(lexer*, std::ostream&, int, int, int, int, int) {};
        virtual void beginningParallel(lexer*, std::ostream&) {};
        virtual void ending(std::ostream&, const int*, int&) {};
        virtual void ending(std::ostream&) {};
        virtual void endingParallel(std::ostream&, const char*, const int, const int) {};
        virtual void endingParallel(std::ostream&) {};

        virtual void structureWrite(lexer*, fdm*, std::vector<char>&, size_t&) {};
        virtual void structureWrite(lexer*, fdm_fnpf*, std::vector<char>&, size_t&) {};
        virtual void structureWrite(lexer*, fdm_nhf*, std::vector<char>&, size_t&) {};
    protected:
        void xmlVersion(std::ostream &result) {result<<"<?xml version=\"1.0\"?>\n";};
        void vtkVersion(std::ostream &result) {result<<"version=\"1.0\" byte_order=\"LittleEndian\">\n";}; // header_type=\"UInt32\"
        // field data: the simulation time (ParaView reads TimeValue as the time of the
        // file) and, with a run log, the run that wrote it (LAGOON RUN.md). VTK writes
        // ASCII strings as their character codes, each string ending in 0.
        void timeValue(std::ostream &result, const double time, const std::string &run_id=std::string())
        {
            result<<"<FieldData>\n<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\"> "
                  <<std::setprecision(15)<<time<<"\n</DataArray>\n";
            if(!run_id.empty())
            {
                result<<"<DataArray type=\"String\" Name=\"RunId\" NumberOfTuples=\"1\" format=\"ascii\">";
                for(unsigned char c : run_id)
                result<<" "<<int(c);
                result<<" 0\n</DataArray>\n";
            }
            result<<"</FieldData>\n";
        };
        void appendData(std::ostream &result) {result<<"<AppendedData encoding=\"raw\">\n_";};
        void structureWriteEnd(std::vector<char> &buffer, size_t &m) {
            std::stringstream result;
            result<<"\n</AppendedData>\n</VTKFile>"<<std::flush;
            std::memcpy(&buffer[m],result.str().data(),result.str().size());};
    public:
        enum type {none, vtu, vts, vtr};
};

#endif
