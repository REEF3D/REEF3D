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
Author: Thomas Becker
--------------------------------------------------------------------*/

#ifndef NHFLOW_DEPANDTIME_AVG_VEL_LINEPROBE_H_
#define NHFLOW_DEPANDTIME_AVG_VEL_LINEPROBE_H_

#include"boundarycheck.h"
#include<iostream>
#include<fstream>

class lexer;
class fdm_nhf;
class ghostcell;

using namespace std;

class nhflow_depandtime_avg_vel_lineprobe : public boundarycheck
{
public:
    nhflow_depandtime_avg_vel_lineprobe(lexer*,fdm_nhf*);
    virtual ~nhflow_depandtime_avg_vel_lineprobe();

    void start(lexer*, fdm_nhf*, ghostcell*);

private:
    void ini_location(lexer*, fdm_nhf*);

    char name[220];

    int *flag;
    int **iloc, **jloc, *point_count;
    int max_points;
    int n,q;
    const int probenum;

    double *line_length, *normal_x, *normal_y;
    double *time_accum, *next_print_time;
    double **xpt, **ypt, **dist;
    double **u_timeint, **v_timeint;

    ofstream *pout;
};

#endif
