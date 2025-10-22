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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"fnpf_vtp_bed.h"
#include"lexer.h"

void fnpf_vtp_bed::pvtp(lexer *p, int num)
{
    ofstream result;
    sprintf(name,"./REEF3D_FNPF_VTP_BED/REEF3D-FNPF-BED-%08i.pvtp",num);
    result.open(name);

    vtp3D::beginningParallel(p,result);

    vtp3D::pointsParallel(result);

    result<<"<PPointData>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"elevation\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"depth\"/>\n";
    result<<"</PPointData>\n";

    char pname[200];
    for(n=0; n<p->M10; ++n)
    {
        sprintf(pname,"REEF3D-FNPF-BED-%08i-%06i.vtp",num,n+1);
        result<<"<Piece Source=\""<<pname<<"\"/>\n";
    }

    vtp3D::endingParallel(result);

    result.close();
}
