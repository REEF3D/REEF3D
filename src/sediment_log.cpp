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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"sediment_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

void sediment_f::log_ini(lexer *p)
{
    if(p->mpirank==0)
    {
    sedlogout.open("./REEF3D_Log/REEF3D_sedimentlog.dat");

    sedlogout<<"number of cells:  "<<p->cellnumtot<<endl<<endl;
    sedlogout<<"#iteration \t #simtime  \t #dtsed \t| #sedtime \t #sediter \t #slidecells \t| #bedmin \t #bedmax \t|"<<endl;
    }
    
}

void sediment_f::sedimentlog(lexer *p)
{
    if(p->mpirank==0)
    {
    sedlogout<<p->count<<"\t \t \t";
    sedlogout<<setprecision(5)<<p->simtime<<" \t ";
    sedlogout<<setprecision(4)<<p->dtsed<<" \t ";
    sedlogout<<setprecision(4)<<p->sedtime<<" \t ";
    sedlogout<<setprecision(4)<<p->sediter<<" \t ";
    sedlogout<<setprecision(4)<<p->slidecells<<" \t ";

    sedlogout<<setprecision(4)<<p->bedmin<<" \t ";
    sedlogout<<setprecision(4)<<p->bedmax<<" \t "<<endl;
    }
    
}