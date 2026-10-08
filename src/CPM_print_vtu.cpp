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
for more details->

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"CPM.h"
#include"lexer.h"
#include"sediment_fdm.h"
#include"ghostcell.h"
#include<fstream>
#include<sstream>
#include<sys/types.h>
#include<cstdio>
#include<cstring>
#include<vector>

void CPM::name_ParaView_parallel_CPM(lexer *p, ofstream &result)
{
    result<<"<PDataArray type=\"Float32\" Name=\"CPM_Ts\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"CPM_Ps\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"CPM_Us\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"CPM_Ws\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"CPM_Kt\"/>\n";
    
    if(p->Q69>0)
    result<<"<PDataArray type=\"Float32\" Name=\"CPM_Pstar\"/>\n";
}

void CPM::name_ParaView_CPM(lexer *p, ostream &result, int *offset, int &n)
{
    result<<"<DataArray type=\"Float32\" Name=\"CPM_Ts\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"CPM_Ps\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"CPM_Us\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"CPM_Ws\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"CPM_Kt\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    
    // seepage (Q 69): piezometric pressure p* = p - rho_f g.x, pore pressure in the bed, fluid above
    if(p->Q69>0)
    {
    result<<"<DataArray type=\"Float32\" Name=\"CPM_Pstar\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    }
}

void CPM::offset_ParaView_CPM(lexer *p, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    
    if(p->Q69>0)
    {
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    }
}

void CPM::print_3D_CPM(lexer* p, ghostcell *pgc, vector<char> &buffer, size_t &m)
{	
	float ffn;
	int iin;
    
    // Ts
    iin=sizeof(float)*p->pointnum;
    memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    
    TPLOOP
    {
        ffn=float(p->ipol4_a(Ts));
        memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);
    }
    
    // Tau
    iin=sizeof(float)*p->pointnum;
    memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    
    TPLOOP
    {
        ffn=float(p->ipol4_a(Tau));
        memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);
    }
    
    // solid velocity u
    iin=sizeof(float)*p->pointnum;
    memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    
    TPLOOP
    {
        ffn=float(p->ipol4_a(Us));
        memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);
    }
    
    // solid velocity w
    iin=sizeof(float)*p->pointnum;
    memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    
    TPLOOP
    {
        ffn=float(p->ipol4_a(Ws));
        memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);
    }
    
    // eddy diffusivity of the parcels (Q 52)
    iin=sizeof(float)*p->pointnum;
    memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    
    TPLOOP
    {
        ffn=float(p->ipol4_a(Kt));
        memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);
    }
    
    // piezometric pressure (Q 69)
    if(p->Q69>0)
    {
    iin=sizeof(float)*p->pointnum;
    memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    
    TPLOOP
    {
        ffn=float(p->ipol4_a(Pse));
        memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);
    }
    }
}

