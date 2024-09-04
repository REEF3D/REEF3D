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

#include "sediment_part.h"
#include "partres.h"

#include "lexer.h"
#include "ghostcell.h"
#include "fdm.h"
#include "vrans_f.h"
#include "reinitopo.h"
#include "ioflow.h"
#include "bedshear.h"
#include "sediment_fdm.h"

/// @brief CFD calculation function
/// @param a fdm object
/// @param pflow IO-flow object
/// @param preto topography reinitialization object
void sediment_part::start_cfd(lexer* p, fdm* a, ghostcell* pgc, ioflow* pflow,
                                    reinitopo* preto, solver* psolv)
{
    double starttime=pgc->timer();
	int xchange=0;
	int removed=0;

    if (p->count>=p->Q43)
	{
        // sediment 
        fill_PQ_cfd(p,a,pgc);
        pslope.slope_cds(p,pgc,&s);
        pbedshear.taubed(p,a,pgc,&s);
        preduce->start(p,pgc,&s);
        pgc->gcsl_start4(p,s.tau_eff,1);
        pbedshear.taucritbed(p,a,pgc,&s);
        pgc->gcsl_start4(p,s.tau_crit,1);
	
        /// runtime seeding
		if(p->Q120==1&&p->count%p->Q121==0)
        posseed_suspended(p,a);
        point_source(p,a);
        
        if(p->Q101>0)
        {
            // topo_influx(p,a);
            // solid_influx(p,a);
            // set_active_topo(p,a);
            // posseed_topo(p,a);
            // seed_topo(p,a);
        }

        /// transport
        if(p->Q101>0)
        pst->erosion(p,*a,PP,s);
        
        pst->move_RK2(p,*a,*pgc,PP,s,*pturb);
        
		xchange=transfer(p,pgc,&PP, pst, maxparticle);
		removed=remove(p,&PP);
        
        if(p->Q101>0)
        pst->deposition(p,*a,PP,s);

        /// topo update
        update_cfd(p,a,pgc,pflow,preto);

        /// cleanup
        if(p->Q20>=0 && p->count%p->Q20==0)
        {
            if(PP.size == 0)
                PP.erase_all();
            PP.optimize();
        }
	}

    /// print out
	print_particles(p);

	gparticle_active = pgc->globalisum(PP.size);
    gremoved = pgc->globalisum(removed);
    gxchange = pgc->globalisum(xchange);

    volume = pst->volume(p,*a,PP);
    volume = pgc->globalsum(volume);

	p->sedsimtime=pgc->timer()-starttime;

    if(p->mpirank==0 && (p->count%p->P12==0))
    	cout<<"Sediment particles: "<<gparticle_active<<" | xch: "<<gxchange<<" rem: "<<gremoved<<" | sim. time: "<<p->sedsimtime<<"\nTotal bed volume change: "<<(volume-volume0)/volume0<<" %, "<<std::setprecision(prec)<<volume-volume0<<" m^3"<<endl;
    debug(p,a,pgc);
}

