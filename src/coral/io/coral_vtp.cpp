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
Architect: Hans Bihs
--------------------------------------------------------------------*/

#include"coral/io/coral_vtp.h"
#include"coral/core/coral_model.h"
#include"coral/structure/coral_beam.h"
#include"coral/structure/coral_superelement.h"
#include"coral/mooring/coral_mooring.h"
#include<fstream>
#include<iomanip>
#include<vector>

namespace coral
{

bool write_vtp(const model &m, const std::string &file)
{
    std::vector<Vec3> pts;
    std::vector<std::vector<int>> lines;
    std::vector<int> kind;
    for(const auto &b : m.bodies)
        pts.push_back(b.x);
    const int nverts = (int)pts.size();
    for(const auto &l : m.load_list)
    {
        if(const auto *be = dynamic_cast<const beam_element*>(l.get()))
        {
            Vec3 pa, pb;
            be->frame(m, pa, pb);
            const int i0 = (int)pts.size();
            pts.push_back(pa);
            pts.push_back(pb);
            lines.push_back({i0, i0+1});
            kind.push_back(1);
        }
        else if(const auto *se = dynamic_cast<const superelement*>(l.get()))
        {
            std::vector<Vec3> np;
            se->node_positions(m, np);
            const int i0 = (int)pts.size();
            pts.insert(pts.end(), np.begin(), np.end());
            for(const auto &e : se->elements())
            {
                lines.push_back({i0 + e.first, i0 + e.second});
                kind.push_back(3);
            }
        }
        else if(const auto *ms = dynamic_cast<const mooring_system*>(l.get()))
        {
            for(int k=0; k<ms->nlines(); ++k)
            {
                std::vector<int> ln;
                for(const Vec3 &p : ms->line_nodes(k))
                {
                    ln.push_back((int)pts.size());
                    pts.push_back(p);
                }
                lines.push_back(ln);
                kind.push_back(2);
            }
        }
    }
    std::ofstream out(file);
    if(!out)
        return false;
    int nconn = 0;
    for(const auto &l : lines) nconn += (int)l.size();
    out<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n<PolyData>\n";
    out<<"<Piece NumberOfPoints=\""<<pts.size()<<"\" NumberOfVerts=\""<<nverts<<"\" NumberOfLines=\""<<lines.size()
       <<"\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
    out<<std::setprecision(10);
    out<<"<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(const Vec3 &p : pts) out<<p(0)<<" "<<p(1)<<" "<<p(2)<<"\n";
    out<<"</DataArray>\n</Points>\n";
    out<<"<CellData Scalars=\"kind\">\n<DataArray type=\"Int32\" Name=\"kind\" format=\"ascii\">\n";
    for(int i=0; i<nverts; ++i) out<<"0\n";
    for(int k : kind) out<<k<<"\n";
    out<<"</DataArray>\n</CellData>\n";
    out<<"<Verts>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for(int i=0; i<nverts; ++i) out<<i<<"\n";
    out<<"</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for(int i=0; i<nverts; ++i) out<<i+1<<"\n";
    out<<"</DataArray>\n</Verts>\n";
    out<<"<Lines>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for(const auto &l : lines) {for(int i : l) out<<i<<" "; out<<"\n";}
    out<<"</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int off = 0;
    for(const auto &l : lines) {off += (int)l.size(); out<<off<<"\n";}
    (void)nconn;
    out<<"</DataArray>\n</Lines>\n</Piece>\n</PolyData>\n</VTKFile>\n";
    return true;
}

}
