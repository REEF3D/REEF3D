/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

iowave::iowave(lexer *p, ghostcell *pgc, patchBC_interface *ppBC)  : wave_interface(p,pgc),flowfile_in(p,pgc),epsi(3.0*p->DXM),psi(0.6*p->DXM), eta(p),
                                          relax1_wg(p),relax1_nb(p),relax2_wg(p),relax2_nb(p),relax4_wg(p),relax4_nb(p)
{
    pBC = ppBC;
    
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
    
    expinverse = 1.0/(exp(1.0)-1.0);
    
    if(p->mpirank==0)
    timeseries(p,pgc);
}

iowave::~iowave()
{
}
