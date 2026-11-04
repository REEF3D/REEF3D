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

#include "fnpf_vtp_bed.h"
#include "fdm_fnpf.h"
#include "ghostcell.h"
#include "lagoon_output.h"
#include "lexer.h"
#include "scour_core.h"
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>

namespace {
lagoon_surface *lagoon_bed = nullptr;  // P 18
}

// the point data of the scour option (P 187) in the order of the file
std::vector<std::pair<const char *, slice *>> fnpf_vtp_bed::scour_fields(scour_core *s)
{
    std::vector<std::pair<const char *, slice *>> f;
    for(int n = 0; n < scour_core::NQ; ++n)
        f.push_back({scour_core::name(n), s->field(n)});
    return f;
}

fnpf_vtp_bed::fnpf_vtp_bed(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if(p->I40 == 0)
        p->printtime = 0.0;

    printcount = 0;

    // Create Folder
    if(p->mpirank == 0)
        mkdir("./REEF3D_FNPF_VTP_BED", 0777);
}

void fnpf_vtp_bed::start(lexer *p, fdm_fnpf *c, ghostcell *pgc, ioflow *pflow)
{
    print2D(p, c, pgc);
}

void fnpf_vtp_bed::print2D(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    int num = 0;
    if(p->P15 == 1)
        num = printcount;
    else if(p->P15 == 2)
        num = p->count;

    if(p->mpirank == 0 && lagoon_output::vtu_files(p))  // P 18 1: the store instead
        pvtp(p, num);

    // offsets
    n = 0;
    offset[n] = 0;
    ++n;

    // Points
    offset[n] = offset[n - 1] + sizeof(float) * p->pointnum2D * 3 + sizeof(int);
    ++n;

    // elevation
    offset[n] = offset[n - 1] + sizeof(float) * p->pointnum2D + sizeof(int);
    ++n;

    // depth
    offset[n] = offset[n - 1] + sizeof(float) * p->pointnum2D + sizeof(int);
    ++n;

    // scour (P 187)
    const int nsc = scour != nullptr ? int(scour_fields(scour).size()) : 0;
    for(int q = 0; q < nsc; ++q)
    {
        offset[n] = offset[n - 1] + sizeof(float) * p->pointnum2D + sizeof(int);
        ++n;
    }

    // Cells
    offset[n] = offset[n - 1] + sizeof(int) * p->polygon_sum * 3 + sizeof(int);
    ++n;
    offset[n] = offset[n - 1] + sizeof(int) * p->polygon_sum + sizeof(int);
    ++n;

    // the piece in memory (lagoon_surface::piece writes the file where it is wanted); every array
    // made in a buffer and written at once: the same bytes as value by value
    sprintf(name, "./REEF3D_FNPF_VTP_BED/REEF3D-FNPF-BED-%08i-%06i.vtp", num, p->mpirank + 1);
    std::ostringstream result(std::ios::out | std::ios::binary);

    vtp3D::beginning(p, result, p->pointnum2D, 0, 0, 0, p->polygon_sum);

    n = 0;
    vtp3D::points(result, offset, n);

    result << "<PointData>\n";
    result << "<DataArray type=\"Float32\" Name=\"elevation\" format=\"appended\" offset=\"" << offset[n] << "\"/>\n";
    ++n;
    result << "<DataArray type=\"Float32\" Name=\"depth\" format=\"appended\" offset=\"" << offset[n] << "\"/>\n";
    ++n;
    for(int q = 0; q < nsc; ++q)
    {
        result << "<DataArray type=\"Float32\" Name=\"" << scour_fields(scour)[q].first << "\" format=\"appended\" offset=\"" << offset[n] << "\"/>\n";
        ++n;
    }
    result << "</PointData>\n";

    vtp3D::polys(result, offset, n);

    vtp3D::ending(result);

    //----------------------------------------------------------------------------

    std::vector<float> fb;
    fb.reserve(size_t(p->pointnum2D) * 3);
    auto put_size = [&](int bytes)
    { result.write((char *)&bytes, sizeof(int)); };
    auto put_floats = [&]()
    { result.write((char*)fb.data(), std::streamsize(sizeof(float)*fb.size())); fb.clear(); };
    // a point array of the slice f, interpolated to the nodes
    auto put_slice = [&](slice &f)
    {
        put_size(sizeof(float) * p->pointnum2D);
        TPSLICELOOP
        fb.push_back(float(p->sl_ipol4(f)));
        put_floats();
    };

    //  XYZ
    put_size(sizeof(float) * p->pointnum2D * 3);
    TPSLICELOOP
    {
        fb.push_back(float(p->Xout(p->XN[IP1], p->YN[JP1])));
        fb.push_back(float(p->Yout(p->XN[IP1], p->YN[JP1])));
        fb.push_back(float(p->sl_ipol4(c->bed)));
    }
    put_floats();

    //  Elevation, depth
    put_slice(c->bed);
    put_slice(c->depth);

    //  Scour (P 187)
    for(int q = 0; q < nsc; ++q)
        put_slice(*scour_fields(scour)[q].second);

    //  Connectivity
    std::vector<int> ib;
    ib.reserve(size_t(p->polygon_sum) * 3);
    put_size(sizeof(int) * (p->polygon_sum) * 3);
    SLICEBASELOOP
    {
        // Triangle 1
        ib.push_back(int(c->nodeval2D(i - 1, j - 1)) - 1);
        ib.push_back(int(c->nodeval2D(i, j - 1)) - 1);
        ib.push_back(int(c->nodeval2D(i, j)) - 1);

        // Triangle 2
        ib.push_back(int(c->nodeval2D(i - 1, j - 1)) - 1);
        ib.push_back(int(c->nodeval2D(i, j)) - 1);
        ib.push_back(int(c->nodeval2D(i - 1, j)) - 1);
    }
    result.write((char *)ib.data(), std::streamsize(sizeof(int) * ib.size()));

    //  Offset of Connectivity
    ib.clear();
    put_size(sizeof(int) * (p->polygon_sum));
    for(n = 0; n < p->polygon_sum; ++n)
        ib.push_back((n + 1) * 3);
    result.write((char *)ib.data(), std::streamsize(sizeof(int) * ib.size()));

    vtp3D::footer(result);

    // the file (P 18 0, 2) and the piece in the LAGOON store (P 18 1: instead of the file, P 18 2: as well)
    std::string piece = std::move(result).str();
    lagoon_surface::piece(p, pgc, lagoon_bed, "FNPF", "bed", "REEF3D_FNPF_VTP_BED", name, num, piece);

    ++printcount;
}
