/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include "vtp3D.h"

#include "lexer.h"

void vtp3D::beginning(lexer *p, std::ostream &result, int numPoints, int numVerts, int numLines, int numStrips, int numPolys)
{
    xmlVersion(result);
    result<<"<VTKFile type=\"PolyData\" ";
    vtkVersion(result);
    result<<"<PolyData>\n";
    if(p->P16==1)
        timeValue(result,p->simtime);
    result<<"<Piece NumberOfPoints=\""<<numPoints<<"\" NumberOfVerts=\""<<numVerts<<"\" NumberOfLines=\""<<numLines<<"\" NumberOfStrips=\""<<numStrips<<"\" NumberOfPolys=\""<<numPolys<<"\">\n";
}

void vtp3D::beginningParallel(lexer *p, std::ostream &result)
{
    xmlVersion(result);
    result<<"<VTKFile type=\"PPolyData\" ";
    vtkVersion(result);
    result<<"<PPolyData GhostLevel=\"0\">\n";
    if(p->P16==1)
        timeValue(result,p->simtime);
}

void vtp3D::points(std::ostream &result, const int *offset, int &n)
{
    result<<"<Points>\n";
    result<<"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"</Points>\n";
}

void vtp3D::pointsParallel(std::ostream &result)
{
    result<<"<PPoints>\n";
    result<<"<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n";
    result<<"</PPoints>\n";
}

void vtp3D::verts(std::ostream &result, const int *offset, int &n)
{
    result<<"<Verts>\n";
    result<<"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"</Verts>\n";
}

void vtp3D::polys(std::ostream& result, const int* offset, int& n)
{
    result<<"<Polys>\n";
    result<<"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"</Polys>\n";
}

void vtp3D::ending(std::ostream& result)
{
    result<<"</Piece>\n";
    result<<"</PolyData>\n";
    appendData(result);
}

void vtp3D::endingParallel(std::ostream& result)
{
    result<<"</PPolyData>\n";
    result<<"</VTKFile>";
}

void vtp3D::footer(std::ostream& result)
{
    result<<"\n</AppendedData>\n";
    result<<"</VTKFile>"<<flush;
}
