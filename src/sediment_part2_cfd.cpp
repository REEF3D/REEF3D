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
Authors: Alexander Hanke, Hans Bihs
--------------------------------------------------------------------*/

#include"sediment_part2.h"
#include"partres.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm.h"
#include"vrans_f.h"
#include"reinitopo.h"
#include"ioflow.h"
#include"bedshear.h"
#include"sediment_fdm.h"
#include"bedslope.h"

void sediment_part2::sediment_algorithm_cfd(lexer* p, fdm* a, ghostcell* pgc, ioflow* pflow,
                                    reinitopo* preto, solver* psolv)
{
    double starttime=pgc->timer();


    if (p->count>=p->Q43)
	{
        // sediment 
        fill_PQ_cfd(p,a,pgc);
        pslope->slope_cds(p,pgc,s);
        pbedshear->taubed(p,a,pgc,s);
        //preduce->start(p,pgc,s);
        pgc->gcsl_start4(p,s->tau_eff,1);
        //pbedshear.taucritbed(p,a,pgc,&s);
        pgc->gcsl_start4(p,s->tau_crit,1);

        //point_source(p,a);
        
        
        //pst->move_RK2(p,*a,*pgc,s,P,*pturb,);
        
        /// topo update
        update_cfd(p,a,pgc,pflow,preto);

	}

    /// print out
	//print_particles(p);
/*
	gparticle_active = pgc->globalisum(PP.size);
    gremoved = pgc->globalisum(removed);
    gxchange = pgc->globalisum(xchanged);

    volume = pst->volume(p,*a,PP);
    volume = pgc->globalsum(volume);

	p->sedsimtime=pgc->timer()-starttime;

    if(p->mpirank==0 && (p->count%p->P12==0))
    	cout<<"Sediment particles: "<<gparticle_active<<" | xch: "<<gxchange<<" rem: "<<gremoved<<" | sim. time: "<<p->sedsimtime<<"\nTotal bed volume change: "<<(volume-volume0)/volume0<<" %, "<<std::setprecision(prec)<<volume-volume0<<" m^3"<<endl;
    debug(p,a,pgc);*/
}
