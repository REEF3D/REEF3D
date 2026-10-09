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

#include"iowave.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"vrans_v.h"
#include"vrans_f.h"
#include"patchBC_interface.h"
#include"linear_regression_cont.h"

iowave::iowave(lexer *p, ghostcell *pgc, patchBC_interface *ppBC)  : wave_interface(p,pgc),flowfile_in(p,pgc),epsi(3.0*p->DXM),psi(0.6*p->DXM), 
                                          eta(p),relax1_wg(p),relax1_nb(p),relax2_wg(p),relax2_nb(p),relax4_wg(p),relax4_nb(p),wgflag(p),
                                          vofheight(p),genheight(p)
{
    pBC = ppBC;
    
    // combined beach (B 99 6): in NHFLOW an absorbing edge at x+ (as B 99 3) behind a relaxation
    // beach (length B 96) that relaxes to the low-passed state, so short waves are damped in the
    // zone and long waves pass to the edge; the solver sees B 99 3. Other solvers: B 99 2.
    if(p->B99==6)
    {
        if(p->A10==5)
        {
        beach_lp = true;
        p->B99 = 3;
        }
        else
        {
        if(p->mpirank==0)
        cout<<"iowave: B 99 6 (combined beach) is available for NHFLOW; runs as B 99 2"<<endl;
        p->B99 = 2;
        }
    }
    
    // decomposed precalc (B 89 1) needs the space / time parts of the wave theory, which only
    // the 5th-order Stokes (B 92 5) and the irregular theories (31, 41, 51) provide; with any
    // other wave type it generated no waves at all, so it now falls back to B 89 0. Additional
    // sources (B 500) are decomposed too (wave_field) and need these types as well.
    if(p->B89==1)
    {
        int bad = -1;
        
        for(int q=0; q<wave_nsources(); ++q)
        {
            int sid, stype;
            double srot;
            wave_source_lib(q,sid,stype,srot);
            
            if(!(q==0 && stype==0) && stype!=5 && stype!=31 && stype!=41 && stype!=51)
            bad = stype;
        }
        
        if(bad>=0)
        {
            if(p->mpirank==0)
            cout<<"iowave: decomposed precalc (B 89 1) is available for wave types 5, 31, 41, 51; wave type "<<bad<<" runs with B 89 0"<<endl;
            
            p->B89 = 0;
        }
    }

    if(p->F80==4)
    vofgen = std::make_unique<field4>(p);

    dist1=p->B96_1;
    dist2=p->B96_2;
    
    dist2_fac=1.0;
    
    if(p->B99==1)
    dist2_fac=2.0;
    
    gcval_press=40;
	
	kinval = 0.00001;	
    
    beach_relax=0;
	
	if(p->T10==1 || p->T10==11 || p->T10==21)
    epsval=(pow(0.09,0.75)*pow(kinval,1.5))/(0.5*0.4*p->DXM);

    if(p->T10==2 || p->T10==12 || p->T10==22)
    epsval=(pow(0.09,0.75)*pow(kinval,0.5))/(0.5*0.4*p->DXM);

    if(p->T10==3 || p->T10==13)
    epsval=(pow(0.09,0.75)*pow(kinval,0.5))/(0.5*0.4*p->DXM);	
	
	// ---------------------------------------
    
    if(p->B105==0 && p->B92!=61)
    {
    p->B105_2 = p->global_xmin;
    p->B105_3 = p->global_ymin;
    }
    
    if(p->B105==0 && p->B92==61)
    {
    p->B105_2 = 0.0;
    p->B105_3 = 0.0;
    }
	
	if(p->B106==0)
	{
	p->B106=1;
	p->Darray(p->B106_b,p->B106);
	p->Darray(p->B106_x,p->B106);
	p->Darray(p->B106_y,p->B106);
	
	p->B106_b[0]=p->B105_1;
	p->B106_x[0]=p->xcoormin;
	p->B106_y[0]=p->ycoormin;
	}
	
    if(p->B107>0)
    {
    p->Darray(B1,p->B107,2);
    p->Darray(B2,p->B107,2);
    p->Darray(B3,p->B107,2);
    p->Darray(B4,p->B107,2);
    p->Darray(Bs,p->B107,2);
    p->Darray(Be,p->B107,2);
    
    beach_relax=1;
    }
    
    if(p->B108>0)
    {
    p->Darray(G1,p->B108,2);
    p->Darray(G2,p->B108,2);
    p->Darray(G3,p->B108,2);
    p->Darray(G4,p->B108,2);
    p->Darray(Gs,p->B108,2);
    p->Darray(Ge,p->B108,2);
    }
    
	if(p->B107==0)
	{
	p->B107=1;
    
    p->Darray(B1,p->B107,2);
    p->Darray(B2,p->B107,2);
    p->Darray(B3,p->B107,2);
    p->Darray(B4,p->B107,2);
    p->Darray(Bs,p->B107,2);
    p->Darray(Be,p->B107,2);
    
    p->Darray(p->B107_xs,p->B107);
    p->Darray(p->B107_xe,p->B107);
    p->Darray(p->B107_ys,p->B107);
    p->Darray(p->B107_ye,p->B107);
    p->Darray(p->B107_d,p->B107);

	p->B107_xs[0]=p->xcoormax;
	p->B107_xe[0]=p->xcoormax;
	p->B107_ys[0]=p->ycoormin-10.0*p->DXM;
    p->B107_ye[0]=p->ycoormax+10.0*p->DXM;
    p->B107_d[0]=p->B96_2;
	}
    
    
    if(p->B108==0)
	{
	p->B108=1;
    
    p->Darray(G1,p->B108,2);
    p->Darray(G2,p->B108,2);
    p->Darray(G3,p->B108,2);
    p->Darray(G4,p->B108,2);
    p->Darray(Gs,p->B108,2);
    p->Darray(Ge,p->B108,2);
    
    p->Darray(p->B108_xs,p->B108);
    p->Darray(p->B108_xe,p->B108);
    p->Darray(p->B108_ys,p->B108);
    p->Darray(p->B108_ye,p->B108);
    p->Darray(p->B108_d,p->B108);

	p->B108_xs[0]=p->xcoormin;
	p->B108_xe[0]=p->xcoormin;
	p->B108_ys[0]=p->ycoormin-10.0*p->DXM;
    p->B108_ye[0]=p->ycoormax+10.0*p->DXM;
    p->B108_d[0]=p->B96_1;
	}
    

    distbeach_ini(p);

    distgen_ini(p);

    zones = bc_zone_set::from_legacy(p,pgc);
    bgs.read(p);
    b530_auto(p);
    zones_check(p);
    nhflow_active_beach_edge(p,pgc);
    
    // tidal / current background (NHFLOW)
    if(zones.has_background() || nhf_active_edge)
    {
    bg_on = zones.has_background();
    p->open_xm = zones.open_edge(1)!=nullptr ? 1 : 0;
    p->open_xp = zones.open_edge(2)!=nullptr ? 1 : 0;
    p->open_ym = zones.open_edge(3)!=nullptr ? 1 : 0;
    p->open_yp = zones.open_edge(4)!=nullptr ? 1 : 0;
    }
    
    if(zones.user_beach())
    beach_relax=1;

	
	p->Darray(beta,p->B106);
	p->Darray(tan_beta,p->B106);

	
	alpha = (p->B105_1+90.0)*(PI/180.0);
	gamma = (p->B105_1)*(PI/180.0);
	
	for(n=0;n<p->B106;++n)
	beta[n] = (p->B106_b[n]+90.0)*(PI/180.0);
	
	tan_alpha = tan(alpha);
	
	for(n=0;n<p->B106;++n)
	tan_beta[n] = tan(beta[n]);
	
	gcawa1_count=gcawa2_count=gcawa3_count=gcawa4_count=1;
	p->Iarray(gcawa1, gcawa1_count, 4); 
	p->Iarray(gcawa2, gcawa2_count, 4); 
	p->Iarray(gcawa3, gcawa3_count, 4); 
	p->Iarray(gcawa4, gcawa4_count, 4); 
    p->Darray(Uoutval,gcawa4_count);
    p->Darray(Fioutval, gcawa4_count);
    p->Darray(Fifsfoutval, gcawa4_count);
	
	gcgen1_count=gcgen2_count=gcgen3_count=gcgen4_count=1;
	p->Iarray(gcgen1, gcgen1_count, 4); 
	p->Iarray(gcgen2, gcgen2_count, 4); 
	p->Iarray(gcgen3, gcgen3_count, 4); 
	p->Iarray(gcgen4, gcgen4_count, 4); 
	
	p->Darray(wsfmax,p->knox,p->knoy);
	
	u_switch=1;
	v_switch=1;
	w_switch=1;
	p_switch=1;
	h_switch=1;
    f_switch=0;
	
	if(p->B92==21 || p->B92==22 || p->B92==23)
	{
	u_switch=1;
	v_switch=1;
	w_switch=0;
	p_switch=0;
	h_switch=0;
    f_switch=0;
    
        if(p->B115==1)
        w_switch=1;
	}
    
    if(p->B92==20)
	{
	u_switch=1;
	v_switch=1;
	w_switch=0;
	p_switch=0;
	h_switch=1;
    f_switch=1;
	}
    
    if(p->A10==3)
	{
        u_switch=0;
        v_switch=0;
        w_switch=0;
        p_switch=0;
        h_switch=1;
        f_switch=1;
        
        if(p->B92==21 || p->B92==22 || p->B92==23)
        {
        u_switch=0;
        v_switch=0;
        w_switch=0;
        p_switch=0;
        h_switch=0;
        f_switch=1;
        }
	}
    
    if(p->B99==1 || p->B99==2)
    beach_relax=1;
    
    if(beach_lp)
    {
        beach_relax=1;
        // default Tp / 4: waves much longer than the peak period pass the zone to the edge
        // (|1 - H| = omega tau / sqrt(1 + omega^2 tau^2)), shorter waves are damped; much smaller
        // tau slows the long waves in the zone and reflects them (validation 17)
        beach_tau = p->B526>0.0 ? p->B526 : 0.25*(p->wTp>0.0 ? p->wTp : p->wT);
        
        if(p->mpirank==0)
        cout<<"iowave: combined beach B 99 6: relaxation over "<<p->B96_2<<" m to the state low-passed over "<<beach_tau<<" s, absorbing edge at x+"<<endl;
    }
    
    expinverse = 1.0/(exp(1.0)-1.0);
    
    if(p->mpirank==0)
    timeseries(p,pgc);
    
    linreg = new linear_regression_cont(p);
    
    
    netV=0.0;
    netV_corr_n=0.0;
}

iowave::~iowave()
{
    delete relax4_wg0;
    delete relax4_nb0;
}
