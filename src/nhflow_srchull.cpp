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
Author: Ahmet Soydan
--------------------------------------------------------------------*/

#include"nhflow_srchull.h"
#include<cstdlib>

namespace nhflow_srchull
{
    std::vector<double> q;
    double x13_mx = 0.0, x13_mz = 0.0, x13_adt = 1.0, x13_dv = 0.0;
    
    bool active()
    {
        static const bool on = getenv("NHF_SRC")!=nullptr;
        return on;
    }
    
    double stream()
    {
        static const double us = getenv("NHF_SRC")!=nullptr ? atof(getenv("NHF_SRC")) : 0.0;
        return us;
    }
    
    int feedback()
    {
        static const int o = getenv("NHF_SRC_FB")!=nullptr ? atoi(getenv("NHF_SRC_FB")) : 0;
        return o;
    }
    
    int opt()
    {
        static const int o = getenv("NHF_SRC_OPT")!=nullptr ? atoi(getenv("NHF_SRC_OPT")) : 0;
        return o;
    }
    
    double k1()
    {
        static const double k = getenv("NHF_SRC_K")!=nullptr ? atof(getenv("NHF_SRC_K")) : 0.0;
        return k;
    }
    
    bool water()
    {
        static const bool w = getenv("NHF_SRC")!=nullptr && getenv("NHF_SRC_WATER")!=nullptr && atoi(getenv("NHF_SRC_WATER"))>0;
        return w;
    }
    
    double qv()
    {
        static const double t = getenv("NHF_SRC_QV")!=nullptr ? atof(getenv("NHF_SRC_QV")) : 0.0;
        return t;
    }
}
