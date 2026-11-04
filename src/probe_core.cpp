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

#include"probe_core.h"
#include"lexer.h"
#include"ghostcell.h"
#include<iomanip>

void probe_rows(lexer *p, ghostcell *pgc, ofstream *pout, double *val, int num, int ncomp)
{
    if(num*ncomp>0)
    pgc->globalmax(val,num*ncomp);

    if(p->mpirank==0)
    for(int n=0;n<num;++n)
    {
    pout[n]<<setprecision(9)<<p->simtime;

    for(int c=0;c<ncomp;++c)
    pout[n]<<" \t "<<val[n*ncomp+c];

    pout[n]<<endl;
    }
}

void gauge_row(lexer *p, ofstream &out, const double *val, int num, int fileFlushMaxCount)
{
    if(p->mpirank==0)
    {
    out<<setprecision(9)<<p->simtime<<"\t";
    for(int n=0;n<num;++n)
    {
        out<<setprecision(9)<<val[n]<<"\t";
        // flush print to disc limited to prevent data loss for many gauges
        if(n%fileFlushMaxCount==0&&n!=0)
            out<<std::flush;
    }
    out<<endl;
    }
}
