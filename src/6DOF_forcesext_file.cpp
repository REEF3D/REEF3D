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
Author: Filip Hahs
--------------------------------------------------------------------*/

#include "6DOF_forcesext_file.h"
#include "fdm.h"
#include "ghostcell.h"
#include "lexer.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
class lexer;
class ghostcell;
void sixdof_forcesext_file::read_format(lexer *p, ghostcell *pgc)
{
    char name[100];
    double val, val0, val1;
    double sign, beta, s;
    int count;

    sprintf(name, "6DOF_forces.dat");

    // open file and determine number of columns
    ifstream file(name, ios_base::in);

    if(!file)
    {
        std::cerr << "no '6DOF_forces.dat' file found"
                  << " rank=" << p->mpirank << std::endl;

        std::abort();
    }

    std::string line;
    std::getline(file, line);

    std::istringstream iss(line);

    colnum = 0;

    while(iss >> val)
        ++colnum;

    std::cerr << "rank=" << p->mpirank << " colnum=" << colnum << std::endl;

    // rewind file
    file.clear();
    file.seekg(0);

    // count rows
    count = 0;

    while(true)
    {
        bool complete_row = true;

        for(qn = 0; qn < colnum; ++qn)
        {
            if(!(file >> val))
            {
                complete_row = false;
                break;
            }
        }

        if(!complete_row)
            break;

        ++count;
    }

    ptnum = count;

    file.close();

    // allocate
    p->Darray(data, ptnum, colnum);

    // re-open file
    file.open(name, ios_base::in);

    if(!file)
    {
        std::cout << std::endl
                  << "no '6DOF_forces.dat' file found" << std::endl
                  << std::endl;

        std::abort();
    }

    // read file
    rowcount = 0;

    while(rowcount < ptnum)
    {
        for(qn = 0; qn < colnum; ++qn)
        {
            if(!(file >> data[rowcount][qn]))
            {
                std::abort();
            }
        }

        ++rowcount;
    }

    ts = data[0][0];
    te = data[ptnum - 1][0];
}
void sixdof_forcesext_file::forcesext_trans(lexer *p, ghostcell *pgc)
{
    Xext = 0.0;
    Yext = 0.0;
    Zext = 0.0;
    if(p->simtime < ts || p->simtime > te)
    {
        return;
    }
    int rank = p->mpirank;

    if(!data)
    {
        return;
    }

    if(!data[timecount])
    {
        return;
    }

    while(timecount < ptnum - 1 && p->simtime > data[timecount][0])
    {
        ++timecount;
    }

    if(timecount == 0)
    {
        Xext = data[0][1];
        Yext = data[0][2];
        Zext = data[0][3];
        return;
    }

    const int i0 = timecount - 1;
    const int i1 = timecount;

    const double t0 = data[i0][0];
    const double t1 = data[i1][0];

    const double alpha = (p->simtime - t0) / (t1 - t0);
    Yext = data[i0][1] + alpha * (data[i1][1] - data[i0][1]);

    Xext = data[i0][2] + alpha * (data[i1][2] - data[i0][2]);

    Zext = data[i0][3] + alpha * (data[i1][3] - data[i0][3]);
}
void sixdof_forcesext_file::forcesext_rot(lexer *p, ghostcell *pgc)
{
    Kext = 0.0;
    Mext = 0.0;
    Next = 0.0;
    if(colnum != 7)
        return;

    if(p->simtime < ts || p->simtime > te)
    {
        return;
    }
    int rank = p->mpirank;

    if(!data)
    {
        return;
    }

    if(!data[timecount])
    {
        return;
    }

    while(timecount < ptnum - 1 && p->simtime > data[timecount][0])
    {
        ++timecount;
    }

    if(timecount == 0)
    {
        Kext = data[0][4];
        Mext = data[0][5];
        Next = data[0][6];
        return;
    }

    const int i0 = timecount - 1;
    const int i1 = timecount;

    const double t0 = data[i0][0];
    const double t1 = data[i1][0];

    const double alpha = (p->simtime - t0) / (t1 - t0);
    Kext = data[i0][4] + alpha * (data[i1][4] - data[i0][4]);

    Mext = data[i0][5] + alpha * (data[i1][5] - data[i0][5]);

    Next = data[i0][6] + alpha * (data[i1][6] - data[i0][6]);
}
sixdof_forcesext_file::sixdof_forcesext_file(lexer *p, ghostcell *g)
{
    ini(p, g);
    read_format(p, g);
}

void sixdof_forcesext_file::ini(lexer *p, ghostcell *g) { return; }
