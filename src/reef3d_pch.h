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

// Precompiled header: force-included into every source by the Makefile (-include) and CMake
// (target_precompile_headers). Only headers that rarely change belong here (standard library,
// MPI, Eigen): a change to this file, or to anything it includes, rebuilds every source.
// Parsing these headers is ~80% of the compile time of a typical REEF3D source without it.
// Sources must still include what they use, so that they also build without it (PCH=0).

#ifndef REEF3D_PCH_H_
#define REEF3D_PCH_H_

#include<cmath>
#include<cstdlib>
#include<cstdint>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>
#include<string>
#include<vector>
#include<map>
#include<mutex>
#include<functional>
#include<algorithm>
#include<memory>
#include<mpi.h>
#include<Eigen/Dense>

#endif
