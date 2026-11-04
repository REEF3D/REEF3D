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
Architect: Hans Bihs
--------------------------------------------------------------------*/

#include"pic_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>
#include<iostream>

pic_f::pic_f(lexer *p, ghostcell *pgc)
{
    static_assert(sizeof(pic_data)==NF*sizeof(double),"pic_data must be plain doubles");

    // the subdomain and the domain
    xloc_s = p->XN[0+marge];
    xloc_e = p->XN[p->knox+marge];
    yloc_s = p->YN[0+marge];
    yloc_e = p->YN[p->knoy+marge];
    zloc_s = p->ZN[0+marge];
    zloc_e = p->ZN[p->knoz+marge];

    xglo_s = pgc->globalmin(xloc_s);
    xglo_e = pgc->globalmax(xloc_e);
    yglo_s = pgc->globalmin(yloc_s);
    yglo_e = pgc->globalmax(yloc_e);
    zglo_s = pgc->globalmin(zloc_s);
    zglo_e = pgc->globalmax(zloc_e);

    npd = p->F83>=1 ? p->F83 : 1;
    npc = npd*npd*(p->j_dir==1 ? npd : 1);

    printtime = 0.0;
    printcount = 0;
    seeded = lost = lost_total = sent = sent_total = 0;

    if(p->mpirank==0)
    {
    mkdir("./REEF3D_CFD_PIC",0777);

    logout.open("./REEF3D_CFD_PIC/REEF3D-CFD-PIC-Log.dat");
    logout<<"t \t particles \t left_domain \t sent_between_ranks"<<endl;

    cout<<"PIC (F 81 "<<p->F81<<"): phase 0, passive particles, "<<npc<<" per cell"<<endl;
    }
}

pic_f::~pic_f()
{
    if(logout.is_open())
    logout.close();
}

void pic_f::ini(lexer *p, fdm *a, ghostcell *pgc)
{
    seed(p,a,pgc);

    for(auto &q : P)
    velocity(p,a,q.x,q.y,q.z,q.u,q.v,q.w);

    // the first output (on a restart, I 45 1, too: the particles are seeded again)
    print(p,pgc);
}

void pic_f::start(lexer *p, fdm *a, ghostcell *pgc)
{
    move(p,a,pgc);
    xchange(p,pgc);

    if(print_due(p))
    print(p,pgc);
}
