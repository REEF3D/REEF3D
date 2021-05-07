/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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
--------------------------------------------------------------------*/

#include"fnpf_state.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include<iostream>
#include<fstream>
#include<sys/stat.h>
#include<sys/types.h>

fnpf_state::fnpf_state(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{	
	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_FNPF_State",0777);
	
	printcount=0;
    
    file_version=1;
    
    if(p->P44==1)
    file_version=2;
    
    ini_token=0;
    
    p->Iarray(flag_all,p->M10);
    p->Iarray(is_flag_all,p->M10);
    p->Iarray(ie_flag_all,p->M10);
    p->Iarray(js_flag_all,p->M10);
    p->Iarray(je_flag_all,p->M10);
    p->Iarray(is_global_all,p->M10);
    p->Iarray(ie_global_all,p->M10);
    p->Iarray(js_global_all,p->M10);
    p->Iarray(je_global_all,p->M10);
    
    flag=0;
    is_flag=ie_flag=js_flag=je_flag=0;
    
    is = is_global = 0;
    ie = ie_global = 0;
    
    js = js_global = p->knox;
    je = je_global = p->knoy;
    
    if(p->P43==1)
    {
        
        if(p->P43_xs<p->endx && p->P43_xe>=p->originx && p->P43_ys<p->endy && p->P43_ye >=p->originy)
        {
        flag=1;
        }
        
        // -----
        if(p->P43_xs>=p->originx && p->P43_xs<p->endx)
        {
        is = p->posc_i(p->P43_xs);
        cout<<p->mpirank<<" IS_LOCAL: "<<is<<" ORIG_I_LOCAL: "<<p->origin_i<<" p->P43_xs: "<<p->P43_xs<<endl;
        is_flag=1;
        }
        
        if(p->P43_xe>=p->originx && p->P43_xe<p->endx)
        {
        ie = p->posc_i(p->P43_xe);
        cout<<p->mpirank<<" IE_LOCAL: "<<ie<<" ORIG_I_LOCAL: "<<p->origin_i<<" p->P43_xe: "<<p->P43_xe<<endl;
        ie_flag=1;
        }
        
        if(p->P43_ys>=p->originy && p->P43_ys<p->endy)
        {
        js = p->posc_j(p->P43_ys);
        cout<<p->mpirank<<" JS_LOCAL: "<<js<<" ORIG_J_LOCAL: "<<p->origin_j<<" p->P43_ys: "<<p->P43_ys<<endl;
        js_flag=1;
        }
        
        if(p->P43_ye>=p->originy && p->P43_ye<p->endy)
        {
        je = p->posc_j(p->P43_ye);
         cout<<p->mpirank<<" JE_LOCAL: "<<je<<" ORIG_J_LOCAL: "<<p->origin_j<<" p->P43_ye: "<<p->P43_ye<<endl;
        je_flag=1;
        }
    }
    
    pgc->gather_int(&flag,1,flag_all,1);
    
    // is communication
    if(is_flag==1)
    is_global = is + p->origin_i;    
    
    pgc->gather_int(&is_flag,1,is_flag_all,1);
    pgc->gather_int(&is_global,1,is_global_all,1);
    
    if(p->mpirank==0)
    {
    for(qn=0;qn<p->M10;++qn)
    if(is_flag_all[qn]==1)
    is_global = is_global_all[qn];
    }
    
    pgc->bcast_int(&is_global,1);
    
    cout<<p->mpirank<<" IS_GLOBAL: "<<is_global<<endl;
    
    
    // ie communication
    if(ie_flag==1)
    ie_global = ie + p->origin_i;    
    
    pgc->gather_int(&ie_flag,1,ie_flag_all,1);
    pgc->gather_int(&ie_global,1,ie_global_all,1);
    
    if(p->mpirank==0)
    {
    for(qn=0;qn<p->M10;++qn)
    if(ie_flag_all[qn]==1)
    ie_global = ie_global_all[qn];
    }
    
    pgc->bcast_int(&ie_global,1);
    
    cout<<p->mpirank<<" IE_GLOBAL: "<<ie_global<<endl;
    
    
    // js communication
    if(js_flag==1)
    js_global = js + p->origin_j;    
    
    pgc->gather_int(&js_flag,1,js_flag_all,1);
    pgc->gather_int(&js_global,1,js_global_all,1);
    
    if(p->mpirank==0)
    {
    for(qn=0;qn<p->M10;++qn)
    if(js_flag_all[qn]==1)
    js_global = js_global_all[qn];
    }
    
    pgc->bcast_int(&js_global,1);
    
    cout<<p->mpirank<<" JS_GLOBAL: "<<js_global<<endl;
    
    
    // je communication
    if(je_flag==1)
    je_global = je + p->origin_j;    
    
    pgc->gather_int(&je_flag,1,je_flag_all,1);
    pgc->gather_int(&je_global,1,je_global_all,1);
    
    if(p->mpirank==0)
    {
    for(qn=0;qn<p->M10;++qn)
    if(je_flag_all[qn]==1)
    je_global = je_global_all[qn];
    }
    
    pgc->bcast_int(&je_global,1);
    
    cout<<p->mpirank<<" JE_GLOBAL: "<<je_global<<endl;
    
    
    
    
    
    p->del_Iarray(flag_all,p->M10);
    p->del_Iarray(is_flag_all,p->M10);
    p->del_Iarray(ie_flag_all,p->M10);
    p->del_Iarray(js_flag_all,p->M10);
    p->del_Iarray(je_flag_all,p->M10);
    p->del_Iarray(is_global_all,p->M10);
    p->del_Iarray(ie_global_all,p->M10);
    p->del_Iarray(js_global_all,p->M10);
    p->del_Iarray(je_global_all,p->M10);
}

fnpf_state::~fnpf_state()
{
}

void fnpf_state::write(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    // header file
    if(ini_token==0)
    {
    if(p->mpirank==0)
    ini_mainheader(p,c,pgc);
    
    if(flag==1)
    write_header(p,c,pgc);
    
    ini_token=1;
    }
    
    if(p->mpirank==0)
    write_mainheader(p,c,pgc);
    
    
    // result file
    if(flag==1)
    write_result(p,c,pgc);
    
}