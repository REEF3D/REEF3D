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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"6DOF_df.h"
#include"6DOF_df_object.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ddweno_f_nug.h"

sixdof_df::sixdof_df(lexer *p, fdm *a, ghostcell *pgc)
{
    if(p->mpirank==0)
    cout<<"6DOF_df startup ..."<<endl;
    
    number6DOF = 1;
    
    for (int nb = 0; nb < number6DOF; nb++)
    p_df_obj.push_back(new sixdof_df_object(p,a,pgc,nb));
    
    
    alpha[0] = 8.0/15.0;
    alpha[1] = 2.0/15.0;
    alpha[2] = 2.0/6.0;
    
    if(p->N40==3 || p->N40==23 || p->N40==33)
    {
    alpha[0] = 1.0;
    alpha[1] = 0.25;
    alpha[2] = 2.0/3.0;
    }
    
    if(p->N40==2 || p->N40==22)
    {
    alpha[0] = 1.0;
    alpha[1] = 0.5;
    }
    
    gamma[0] = 8.0/15.0;
    gamma[1] = 5.0/12.0;
    gamma[2] = 3.0/4.0;
    
    zeta[0] = 0.0;
    zeta[1] = -17.0/60.0;
    zeta[2] = -5.0/12.0;
}
    
sixdof_df::~sixdof_df()
{
}

void sixdof_df::initialize(lexer *p, fdm *a, ghostcell *pgc, vector<net*>& pnet)
{
    for (int nb = 0; nb < number6DOF; nb++)
    p_df_obj[nb]->initialize(p, a, pgc, pnet);
}
	
void sixdof_df::start(lexer*,fdm*,ghostcell*,double,vrans*,vector<net*>&)
{
}

void sixdof_df::start_forcing(lexer* p, fdm* a, ghostcell* pgc, vrans* pvrans, vector<net*>& pnet, int iter, field &uvel, field &vvel, field &wvel, field &fx, field &fy, field &fz, bool finalise)
{
    // Reset heaviside field
    ULOOP
    a->fbh1(i,j,k) = 0.0;

    VLOOP
    a->fbh2(i,j,k) = 0.0;
    
    WLOOP
    a->fbh3(i,j,k) = 0.0;

    LOOP
    a->fbh4(i,j,k) = 0.0;

    pgc->start1(p,a->fbh1,10);
    pgc->start2(p,a->fbh2,11);
    pgc->start3(p,a->fbh3,12);
    pgc->start4(p,a->fbh4,40);
    
    double starttime,endtime;
 
    for (int nb=0; nb<number6DOF;++nb)
    {
        starttime=pgc->timer();
        
        // Calculate forces
        p_df_obj[nb]->forces_stl(p,a,pgc,alpha[iter],uvel,vvel,wvel);
        
        endtime=pgc->timer();
        
        //if(p->mpirank==0)
        //cout<<"T0: "<<endtime-starttime<<endl;
        
        starttime=pgc->timer();
        
        // Advance body in time
        p_df_obj[nb]->solve_eqmotion(p,a,pgc,iter,pvrans,pnet);
        
        endtime=pgc->timer();
        
        //if(p->mpirank==0)
        //cout<<"T1: "<<endtime-starttime<<endl;
        
        starttime=pgc->timer();

        // Update position and fb level set
        p_df_obj[nb]->transform(p,a,pgc,finalise);  //----> main time consumer
        
        endtime=pgc->timer();
        
        //if(p->mpirank==0)
        //cout<<"T2: "<<endtime-starttime<<endl;
        
        starttime=pgc->timer();

        // Update forcing terms
        p_df_obj[nb]->updateForcing(p,a,pgc,alpha[iter],uvel,vvel,wvel,fx,fy,fz);
        
        endtime=pgc->timer();
        
        //if(p->mpirank==0)
        //cout<<"T3: "<<endtime-starttime<<endl;
        
        starttime=pgc->timer();

        // Save and print
        p_df_obj[nb]->interface(p,true);
        
        endtime=pgc->timer();
        
        //if(p->mpirank==0)
        //cout<<"T4: "<<endtime-starttime<<endl;
        
        starttime=pgc->timer();

        if (finalise == true)
        {
            p_df_obj[nb]->saveTimeStep(p,alpha[iter]);
            
            if(p->X50==1)
            p_df_obj[nb]->print_vtp(p,a,pgc);
            
            if(p->X50==2)
            p_df_obj[nb]->print_stl(p,a,pgc);
            
            p_df_obj[nb]->print_parameter(p, a, pgc);
        }
        endtime=pgc->timer();
        
        //if(p->mpirank==0)
        //cout<<"T5: "<<endtime-starttime<<endl;
        
        starttime=pgc->timer();
    }
    
    // ghostcell update
    pgc->gcdf_update(p,a);
}
