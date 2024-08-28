/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Alexander Hanke
--------------------------------------------------------------------*/

#include "sedpart.h"

#include "lexer.h"
#include "ghostcell.h"
#include "sediment_fdm.h"

void sedpart::name_pvtu_bedshear(lexer *p, ghostcell *pgc, ofstream &result)
{
    if(p->P79==1)
    {
    result<<"<PDataArray type=\"Float32\" Name=\"ST_tau_eff\"/>"<<endl;
    result<<"<PDataArray type=\"Float32\" Name=\"ST_tau_crit\"/>"<<endl;
    }
    
    if(p->P79==2)
    {
    result<<"<PDataArray type=\"Float32\" Name=\"ST_shearvel_eff\"/>"<<endl;
    result<<"<PDataArray type=\"Float32\" Name=\"ST_shearvel_crit\"/>"<<endl;
    }
    
    if(p->P79==3)
    {
    result<<"<PDataArray type=\"Float32\" Name=\"ST_shields_eff\"/>"<<endl;
    result<<"<PDataArray type=\"Float32\" Name=\"ST_shields_crit\"/>"<<endl;
    }
}

void sedpart::name_vtu_bedshear(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    if(p->P79==1)
    {
    result<<"<DataArray type=\"Float32\" Name=\"ST_tau_eff\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_tau_crit\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    }
    
    if(p->P79==2)
    {
    result<<"<DataArray type=\"Float32\" Name=\"ST_shearvel_eff\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_shearvel_crit\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    }
    
    if(p->P79==3)
    {
    result<<"<DataArray type=\"Float32\" Name=\"ST_shields_eff\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_shields_crit\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    }
}

void sedpart::offset_vtu_bedshear(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
}

void sedpart::print_3D_bedshear(lexer* p, ghostcell *pgc, ofstream &result)
{	
	float ffn;
	int iin;
    double tau_eff,shearvel_eff,shields_eff;
    double tau_crit, shearvel_crit, shields_crit;
    
    if(p->P79==1)
    {
        
	// tau_eff
    pgc->gcsl_start4(p,s.tau_eff,1);
    
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s.tau_eff));
	result.write((char*)&ffn, sizeof (float));
	}
    
    // tau_crit
    pgc->gcsl_start4(p,s.tau_crit,1);
    
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s.tau_crit));
	result.write((char*)&ffn, sizeof (float));
	}
    }
    
    
    if(p->P79==2)
    {
	// shearvel_eff
    pgc->gcsl_start4(p,s.shearvel_eff,1);
    
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s.shearvel_eff));
	result.write((char*)&ffn, sizeof (float));
	}
    
    // shearvel_crit
    pgc->gcsl_start4(p,s.shearvel_crit,1);
    
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s.shearvel_crit));
	result.write((char*)&ffn, sizeof (float));
	}
    }
    
    
    if(p->P79==3)
    {
	// shields_eff
    pgc->gcsl_start4(p,s.shields_eff,1);
    
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s.shields_eff));
	result.write((char*)&ffn, sizeof (float));
	}
    
    // shields_crit
    pgc->gcsl_start4(p,s.shields_crit,1);
    
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s.shields_crit));
	result.write((char*)&ffn, sizeof (float));
	}
    
    }
}

