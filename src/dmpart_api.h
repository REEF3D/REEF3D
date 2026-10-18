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

#ifndef DMPART_API_H_
#define DMPART_API_H_

// dmpart: the domain decomposition in REEF3D (M 40). From the global grid file of DIVEMesh
// (DIVEMesh M 50 1, 2) the grid files of nranks ranks, as DIVEMesh would write them for M 10 nranks,
// built in memory (format v2, gridfile_v2.h): rankbuf[r] for rank r.  The decomposition method and
// its keys (M 11-13, 20, 21, 30-35, 41, 45) are those of the DIVEMesh control file.
// false and a message if the file is missing or inconsistent.

#include<vector>
#include<string>

namespace dmpart
{
// weighted: 1 the weighted partition with the load map (M 20 5, M 35 1) whatever the DIVEMesh keys are
// boxes: if given, the subdomains (i0, j0, k0, ni, nj, nk per rank); boxes_only: no grid files built
bool build(const char *fname, int nranks, std::vector<std::vector<char>> &rankbuf, std::string &msg, int weighted=0,
           std::vector<int> *boxes=nullptr, bool boxes_only=false);
}

#endif
