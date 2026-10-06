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
Authors: Tobias Martin, Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"reinidisc_f.h"
#include"nhflow_reinidisc_fsf.h"
#include"6DOF_motionext_fixed.h"
#include"6DOF_motionext_file.h"
#include"6DOF_motionext_CoG.h"
#include"6DOF_motionext_wavemaker.h"
#include"6DOF_motionext_void.h"
#include"net_interface.h"
#include"ship.h"

sixdof_obj::sixdof_obj(lexer *p, ghostcell *pgc, int number) : ddweno_f_nug(p),
                                                                                georay(p),n6DOF(number),
                                                                                epsifb(1.6*p->DXM), epsi(1.6),
                                                                                interfac(1.6),zero(0.0),
                                                                                Mass_fb(rb.mass),
                                                                                p_(rb.p),c_(rb.c),h_(rb.h),dc_(rb.dc),e_(rb.e),
                                                                                R_(rb.R),I_(rb.I),quatRotMat(rb.R),
                                                                                omega_B(rb.omega_B),omega_I(rb.omega_I),
                                                                                phi(rb.phi),theta(rb.theta),psi(rb.psi),
                                                                                Ffb_(rb.F),Mfb_(rb.M),
                                                                                geom(number),amr_hfac(geom.amr_hfac),
                                                                                tri_x(geom.tri_x),tri_y(geom.tri_y),tri_z(geom.tri_z),
                                                                                tri_x0(geom.tri_x0),tri_y0(geom.tri_y0),tri_z0(geom.tri_z0),
                                                                                entity_sum(geom.entity_sum),tstart(geom.tstart),tend(geom.tend),
                                                                                tricount(geom.tricount),entity_count(geom.entity_count)
{
    // rigid-body core: DOF modes (X 11) and linear damping (X 25, X 26)
    rb.dof[0] = p->X11_u;
    rb.dof[1] = p->X11_v;
    rb.dof[2] = p->X11_w;
    rb.dof[3] = p->X11_p;
    rb.dof[4] = p->X11_q;
    rb.dof[5] = p->X11_r;
    rb.twoD = (p->j_dir==0);
    rb.Cdamp_t[0] = p->X26_Cu;
    rb.Cdamp_t[1] = p->X26_Cv;
    rb.Cdamp_t[2] = p->X26_Cw;
    rb.Cdamp_r[0] = p->X25_Cp;
    rb.Cdamp_r[1] = p->X25_Cq;
    rb.Cdamp_r[2] = p->X25_Cr;
    
    // ship module (X 350): hull resistance, damping and propulsion models from ship.dat
    if(p->X350==1)
    pload.push_back(new ship(p,number));

    pnetinter = new net_interface(p,pgc);
    
    triangle_token=0;
    printnormal_count=0;
    
    alpha[0] = 8.0/15.0;
    alpha[1] = 2.0/15.0;
    alpha[2] = 2.0/6.0;
    
    gamma[0] = 8.0/15.0;
    gamma[1] = 5.0/12.0;
    gamma[2] = 3.0/4.0;
    
    zeta[0] = 0.0;
    zeta[1] = -17.0/60.0;
    zeta[2] = -5.0/12.0;
    
    if(((p->N40==3 || p->N40==13 || p->N40==23 || p->N40==33) && p->A10==6) || (p->A510==3 && p->A10==5) || (p->A210==3 && p->A10==2)) 
    {
    alpha[0] = 1.0;
    alpha[1] = 0.25;
    alpha[2] = 2.0/3.0;
    
    gamma[0] = 0.0;
    gamma[1] = 0.0;
    gamma[2] = 0.0;
    
    zeta[0] = 0.0;
    zeta[1] = 0.0;
    zeta[2] = 0.0;
    }
    
    if(((p->N40==2 || p->N40==12 || p->N40==22) && p->A10==6) || (p->A510==2 && p->A10==5) || (p->A210==2 && p->A10==2)) 
    {
    alpha[0] = 1.0;
    alpha[1] = 0.5;
    
    gamma[0] = 0.0;
    gamma[1] = 0.0;
    gamma[2] = 0.0;
    
    zeta[0] = 0.0;
    zeta[1] = 0.0;
    zeta[2] = 0.0;
    }
    
    
    if(p->X210==0 && p->X211==0)
    pmotion = new sixdof_motionext_void(p,pgc);
    
    if((p->X210==1 || p->X211==1) && p->X240==0)
    pmotion = new sixdof_motionext_fixed(p,pgc);
    
    if(p->X240==1)
    pmotion = new sixdof_motionext_file(p,pgc);
    
    if(p->X240==11)
    pmotion = new sixdof_motionext_file_CoG(p,pgc);
    
    if(p->X240==21)
    pmotion = new sixdof_motionext_wavemaker(p,pgc);
    
    Mass_fb =  Rfb = Vfb = 1.0;
    
    // porous floating body (X 16): Darcy-Forchheimer coefficients, same closure as vrans_nhflow_f
    // K = Apor_fb*visc + Bpor_fb*|u - u_fb|  [1/s]
    Xd=Yd=Zd=Kd=Md=Nd=0.0;
    Apor_fb=Bpor_fb=0.0;
    Dpor_t=Dpor_r[0]=Dpor_r[1]=Dpor_r[2]=0.0;
    
    if(p->X16==1)
    {
    Apor_fb = p->X16_alpha*(pow(1.0-p->X16_n,2.0)/pow(p->X16_n,3.0))/pow(p->X16_d50,2.0);
    Bpor_fb = p->X16_beta*(1.0 + 7.5/p->B264)*((1.0-p->X16_n)/pow(p->X16_n,4.0))/p->X16_d50;
    }
    
    
    
    if(p->X10==4)
    {
    p->Darray(uwm,(p->kmax+2)); 
    p->Darray(wwm,(p->kmax+2));   
    }
}

sixdof_obj::~sixdof_obj()
{
    for(size_t ql=0; ql<pload.size(); ++ql)
    delete pload[ql];
}
    
