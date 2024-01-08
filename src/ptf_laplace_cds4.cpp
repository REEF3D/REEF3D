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

#include"ptf_laplace_cds4.h"
#include"lexer.h"
#include"fdm_ptf.h"
#include"solver_ptf.h"
#include"ghostcell.h"

ptf_laplace_cds4::ptf_laplace_cds4() 
{
}

ptf_laplace_cds4::~ptf_laplace_cds4()
{
}

void ptf_laplace_cds4::start(lexer* p, fdm_ptf *e, ghostcell *pgc, solver_ptf *psolv, field &f, slice &Fifsf, slice&eta)
{
    // see p. 1130-1132
    n=0;
    LOOP
	{
	e->M.p[n]  =  1.0/(p->DXP[IP]*p->DXN[IP])*p->x_dir 
                + 1.0/(p->DXP[IM1]*p->DXN[IP])*p->x_dir 
                
                + 1.0/(p->DYP[JP]*p->DYN[JP])*p->y_dir 
                + 1.0/(p->DYP[JM1]*p->DYN[JP])*p->y_dir 
                
                + 1.0/(p->DZP[KP]*p->DZN[KP])*p->z_dir
                + 1.0/(p->DZP[KM1]*p->DZN[KP])*p->z_dir;

   	e->M.n[n] = -1.0/(p->DXP[IP]*p->DXN[IP])*p->x_dir;
	e->M.s[n] = -1.0/(p->DXP[IM1]*p->DXN[IP])*p->x_dir;

	e->M.w[n] = -1.0/(p->DYP[JP]*p->DYN[JP])*p->y_dir;
	e->M.e[n] = -1.0/(p->DYP[JM1]*p->DYN[JP])*p->y_dir;

	e->M.t[n] = -1.0/(p->DZP[KP]*p->DZN[KP])*p->z_dir;
	e->M.b[n] = -1.0/(p->DZP[KM1]*p->DZN[KP])*p->z_dir;

	e->rhsvec.V[n] = 0.0;
	
	++n;
	}
    
    
    n=0;
	LOOP
	{
        if(p->flag4[IJK]>0)
        {
    
            if(p->flag4[Im1JK]<0)
            {
            e->rhsvec.V[n] -= e->M.s[n]*f(i-1,j,k);
            e->M.s[n] = 0.0;
            }
            
            if(p->flag4[Ip1JK]<0)
            {
            e->rhsvec.V[n] -= e->M.n[n]*f(i+1,j,k);
            e->M.n[n] = 0.0;
            }
            
            if(p->flag4[IJm1K]<0)
            {
            e->rhsvec.V[n] -= e->M.e[n]*f(i,j-1,k);
            e->M.e[n] = 0.0;
            }
            
            if(p->flag4[IJp1K]<0)
            {
            e->rhsvec.V[n] -= e->M.w[n]*f(i,j+1,k);
            e->M.w[n] = 0.0;
            }
            
            if(p->flag4[IJKm1]<0)
            {
            e->rhsvec.V[n] -= e->M.b[n]*f(i,j,k-1);
            e->M.b[n] = 0.0;
            }
            
            if(p->flag4[IJKp1]<0)
            {
            e->rhsvec.V[n] -= e->M.t[n]*f(i,j,k+1);
            e->M.t[n] = 0.0;
            }
        }

	++n;
	}
    psolv->start(p,e,pgc,e->Fi,e->rhsvec,5);
    pgc->start4(p,e->Fi,250);
    
    
    // 4th-order
	n=0;
    LOOP
	{
    X0 = -0.5*p->XP[IP2] + 13.0*p->XP[IP1] - 13.0*p->XP[IM1] + 0.5*p->XP[IM2];
    X1 = (-p->XP[IP3] + 27.0*p->XP[IP2] -27.0*p->XP[IP1] + p->XP[IP])*X0;
    X2 = (-p->XP[IP2] + 27.0*p->XP[IP1] -27.0*p->XP[IP] + p->XP[IM1])*X0;
    X3 = (-p->XP[IP1] + 27.0*p->XP[IP] -27.0*p->XP[IM1] + p->XP[IM2])*X0;
    X4 = (-p->XP[IP] + 27.0*p->XP[IM1] -27.0*p->XP[IM2] + p->XP[IM3])*X0;
    
    Y0 = -0.5*p->YP[JP2] + 13.0*p->YP[JP1] - 13.0*p->YP[JM1] + 0.5*p->YP[JM2];
    Y1 = (-p->YP[JP3] + 27.0*p->YP[JP2] -27.0*p->YP[JP1] + p->YP[JP])*Y0;
    Y2 = (-p->YP[JP2] + 27.0*p->YP[JP1] -27.0*p->YP[JP] + p->YP[JM1])*Y0;
    Y3 = (-p->YP[JP1] + 27.0*p->YP[JP] -27.0*p->YP[JM1] + p->YP[JM2])*Y0;
    Y4 = (-p->YP[JP] + 27.0*p->YP[JM1] -27.0*p->YP[JM2] + p->YP[JM3])*Y0;
    
    Z0 = -0.5*p->ZP[KP2] + 13.0*p->ZP[KP1] - 13.0*p->ZP[KM1] + 0.5*p->ZP[KM2];
    Z1 = (-p->ZP[KP3] + 27.0*p->ZP[KP2] -27.0*p->ZP[KP1] + p->ZP[KP])*Z0;
    Z2 = (-p->ZP[KP2] + 27.0*p->ZP[KP1] -27.0*p->ZP[KP] + p->ZP[KM1])*Z0;
    Z3 = (-p->ZP[KP1] + 27.0*p->ZP[KP] -27.0*p->ZP[KM1] + p->ZP[KM2])*Z0;
    Z4 = (-p->ZP[KP] + 27.0*p->ZP[KM1] -27.0*p->ZP[KM2] + p->ZP[KM3])*Z0;
    

	e->M.p[n] = (1.0/X1 + 729.0/X2 + 729.0/X3 + 1.0/X4)*p->x_dir 
             + (1.0/Y1 + 729.0/Y2 + 729.0/Y3 + 1.0/Y4)*p->y_dir 
             + (1.0/Z1 + 729.0/Z2 + 729.0/Z3 + 1.0/Z4)*p->z_dir;
    
    
   	e->M.n[n] = -(27.0/X1 + 729.0/X2 + 27.0/X3)*p->x_dir;
	e->M.s[n] = -(27.0/X2 + 729.0/X3 + 27.0/X4)*p->x_dir;

	e->M.w[n] = -(27.0/Y1 + 729.0/Y2 + 27.0/Y3)*p->y_dir;
	e->M.e[n] = -(27.0/Y2 + 729.0/Y3 + 27.0/Y4)*p->y_dir;

	e->M.t[n] = -(27.0/Z1 + 729.0/Z2 + 27.0/Z3)*p->z_dir;
	e->M.b[n] = -(27.0/Z2 + 729.0/Z3 + 27.0/Z4)*p->z_dir;
    
    
    e->M.nn[n] = (27.0/X1 + 27.0/X2)*p->x_dir; 
    e->M.ss[n] = (27.0/X3 + 27.0/X4)*p->x_dir;
    
    e->M.ww[n] = (27.0/Y1 + 27.0/Y2)*p->y_dir;
    e->M.ee[n] = (27.0/Y3 + 27.0/Y4)*p->y_dir;
     
    e->M.tt[n] = (27.0/Z1 + 27.0/Z2)*p->z_dir; 
    e->M.bb[n] = (27.0/Z3 + 27.0/Z4)*p->z_dir;
     
     
    e->M.nnn[n] = -(1.0/X1)*p->x_dir; 
    e->M.sss[n] = -(1.0/X4)*p->x_dir;
    
    e->M.www[n] = -(1.0/Y1)*p->y_dir;
    e->M.eee[n] = -(1.0/Y4)*p->y_dir;
     
    e->M.ttt[n] = -(1.0/Z1)*p->z_dir; 
    e->M.bbb[n] = -(1.0/Z4)*p->z_dir;
    

//cout<<e->M.p[n]<<" "<<e->M.n[n]<<" "<<e->M.s[n]<<" "<<e->M.t[n]<<" "<<e->M.b[n]<<endl;
	e->rhsvec.V[n] = 0.0;//((f(i+2,j,k)*(-27.0/X1 - 27.0/X2) + f(i-2,j,k)*(-27.0/X3 - 27.0/X4))*p->x_dir
                 // +  (f(i,j+2,k)*(-27.0/Y1 - 27.0/Y2) + f(i,j-2,k)*(-27.0/Y3 - 27.0/Y4))*p->y_dir
                 // +  (f(i,j,k+2)*(-27.0/Z1 - 27.0/Z2) + f(i,j,k-2)*(-27.0/Z3 - 27.0/Z4))*p->z_dir
                  
    //e->rhsvec.V[n] =   (f(i+3,j,k)*(1.0/X1) + f(i-3,j,k)*(1.0/X4))*p->x_dir
      //              +  (f(i,j+3,k)*(1.0/Y1) + f(i,j-3,k)*(1.0/Y4))*p->y_dir
        //            +  (f(i,j,k+3)*(1.0/Z1) + f(i,j,k-3)*(1.0/Z4))*p->z_dir;

	++n;
	}
    
    n=0;
	LOOP
	{
        if(p->flag4[IJK]>0)
        {
    
            if(p->flag4[Im1JK]<0)
            {
            e->rhsvec.V[n] -= e->M.s[n]*f(i-1,j,k);
            e->M.s[n] = 0.0;
            }
            
            if(p->flag4[Ip1JK]<0)
            {
            e->rhsvec.V[n] -= e->M.n[n]*f(i+1,j,k);
            e->M.n[n] = 0.0;
            }
            
            if(p->flag4[IJm1K]<0)
            {
            e->rhsvec.V[n] -= e->M.e[n]*f(i,j-1,k);
            e->M.e[n] = 0.0;
            }
            
            if(p->flag4[IJp1K]<0)
            {
            e->rhsvec.V[n] -= e->M.w[n]*f(i,j+1,k);
            e->M.w[n] = 0.0;
            }
            
            if(p->flag4[IJKm1]<0)
            {
            e->rhsvec.V[n] -= e->M.b[n]*f(i,j,k-1);
            e->M.b[n] = 0.0;
            }
            
            if(p->flag4[IJKp1]<0)
            {
            e->rhsvec.V[n] -= e->M.t[n]*f(i,j,k+1);
            e->M.t[n] = 0.0;
            }
            
            //--
            if(p->flag4[Im2JK]<0)
            {
            e->rhsvec.V[n] -= e->M.ss[n]*f(i-2,j,k);
            e->M.ss[n] = 0.0;
            }
            
            if(p->flag4[Ip2JK]<0)
            {
            e->rhsvec.V[n] -= e->M.nn[n]*f(i+2,j,k);
            e->M.nn[n] = 0.0;
            }
            
            if(p->flag4[IJm2K]<0)
            {
            e->rhsvec.V[n] -= e->M.ee[n]*f(i,j-2,k);
            e->M.ee[n] = 0.0;
            }
            
            if(p->flag4[IJp2K]<0)
            {
            e->rhsvec.V[n] -= e->M.ww[n]*f(i,j+2,k);
            e->M.ww[n] = 0.0;
            }
            
            if(p->flag4[IJKm2]<0)
            {
            e->rhsvec.V[n] -= e->M.bb[n]*f(i,j,k-2);
            e->M.bb[n] = 0.0;
            }
            
            if(p->flag4[IJKp2]<0)
            {
            e->rhsvec.V[n] -= e->M.tt[n]*f(i,j,k+2);
            e->M.tt[n] = 0.0;
            }
            
            
            //--
            
            if(p->flag4[Im3JK]<0)
            {
            e->rhsvec.V[n] -= e->M.sss[n]*f(i-3,j,k);
            e->M.sss[n] = 0.0;
            }
            
            if(p->flag4[Ip3JK]<0)
            {
            e->rhsvec.V[n] -= e->M.nnn[n]*f(i+3,j,k);
            e->M.nnn[n] = 0.0;
            }
            
            if(p->flag4[IJm3K]<0)
            {
            e->rhsvec.V[n] -= e->M.eee[n]*f(i,j-3,k);
            e->M.eee[n] = 0.0;
            }
            
            if(p->flag4[IJp3K]<0)
            {
            e->rhsvec.V[n] -= e->M.www[n]*f(i,j+3,k);
            e->M.www[n] = 0.0;
            }
            
            if(p->flag4[IJKm3]<0)
            {
            e->rhsvec.V[n] -= e->M.bbb[n]*f(i,j,k-3);
            e->M.bbb[n] = 0.0;
            }
            
            if(p->flag4[IJKp3]<0)
            {
            e->rhsvec.V[n] -= e->M.ttt[n]*f(i,j,k+3);
            e->M.ttt[n] = 0.0;
            }
        }

	++n;
	}
    
    
    double starttime=pgc->timer();
    psolv->start(p,e,pgc,e->Fi,e->rhsvec,7);
    double endtime=pgc->timer();
    pgc->start4(p,e->Fi,250);
    
    p->poissoniter=p->solveriter;
    p->poissontime=endtime-starttime;
	if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"Fi_iter: "<<p->solveriter<<"  Fi_time: "<<setprecision(3)<<p->poissontime<<endl;
}




/*
    double ddx,ddy,ddz;
    
    // + (i+2) - 16(i+1) + 30(i) - 16(i-1) + (i-2) / 12dx^2
    
	n=0;
    LOOP
	{
    ddx = -2.0*pow(p->DXP[IP1],2.0) + 8.0*pow(p->DXP[IP],2.0) + 8.0*pow(p->DXP[IM1],2.0) - 2.0*pow(p->DXP[IM2],2.0);
    ddy = -2.0*pow(p->DYP[JP1],2.0) + 8.0*pow(p->DYP[JP],2.0) + 8.0*pow(p->DYP[JM1],2.0) - 2.0*pow(p->DYP[JM2],2.0);
    ddz = -2.0*pow(p->DZP[KP1],2.0) + 8.0*pow(p->DZP[KP],2.0) + 8.0*pow(p->DZP[KM1],2.0) - 2.0*pow(p->DZP[KM2],2.0);
    
    //cout<<ddx<<" "<<ddy<<" "<<ddz<<endl;
	e->M.p[n] =  30.0/ddx*p->x_dir + 30.0/ddy*p->y_dir + 30.0/ddz*p->z_dir;

   	e->M.n[n] = -16.0/ddx*p->x_dir;
	e->M.s[n] = -16.0/ddx*p->x_dir;

	e->M.w[n] = -16.0/ddy*p->y_dir;
	e->M.e[n] = -16.0/ddy*p->y_dir;

	e->M.t[n] = -16.0/ddz*p->z_dir;
	e->M.b[n] = -16.0/ddz*p->z_dir;
    
    e->M.nn[n] = 1.0/ddx*p->x_dir;
	e->M.ss[n] = 1.0/ddx*p->x_dir;

	e->M.ww[n] = 1.0/ddy*p->y_dir;
	e->M.ee[n] = 1.0/ddy*p->y_dir;

	e->M.tt[n] = 1.0/ddz*p->z_dir;
	e->M.bb[n] = 1.0/ddz*p->z_dir;
    
    //cout<<e->M.p[n]<<" "<<e->M.n[n]<<" "<<e->M.s[n]<<" "<<e->M.t[n]<<" "<<e->M.b[n]<<endl;
	
	e->rhsvec.V[n] = 0.0;//(f(i+2,j,k)/ddx*p->x_dir + f(i-2,j,k)/ddx*p->x_dir
                  //+   f(i,j+2,k)/ddy*p->y_dir + f(i,j-2,k)/ddy*p->y_dir
                  //+   f(i,j,k+2)/ddz*p->z_dir + f(i,j,k-2)/ddz*p->z_dir);
	++n;
	}
    
    n=0;
	LOOP
	{
        if(p->flag4[IJK]>0)
        {
    
            if(p->flag4[Im1JK]<0)
            {
            e->rhsvec.V[n] -= e->M.s[n]*f(i-1,j,k);
            e->M.s[n] = 0.0;
            }
            
            if(p->flag4[Ip1JK]<0)
            {
            e->rhsvec.V[n] -= e->M.n[n]*f(i+1,j,k);
            e->M.n[n] = 0.0;
            }
            
            if(p->flag4[IJm1K]<0)
            {
            e->rhsvec.V[n] -= e->M.e[n]*f(i,j-1,k);
            e->M.e[n] = 0.0;
            }
            
            if(p->flag4[IJp1K]<0)
            {
            e->rhsvec.V[n] -= e->M.w[n]*f(i,j+1,k);
            e->M.w[n] = 0.0;
            }
            
            if(p->flag4[IJKm1]<0)
            {
            e->rhsvec.V[n] -= e->M.b[n]*f(i,j,k-1);
            e->M.b[n] = 0.0;
            }
            
            if(p->flag4[IJKp1]<0)
            {
            e->rhsvec.V[n] -= e->M.t[n]*f(i,j,k+1);
            e->M.t[n] = 0.0;
            }
            
            
            //--
            if(p->flag4[Im2JK]<0)
            {
            e->rhsvec.V[n] -= e->M.ss[n]*f(i-2,j,k);
            e->M.ss[n] = 0.0;
            }
            
            if(p->flag4[Ip2JK]<0)
            {
            e->rhsvec.V[n] -= e->M.nn[n]*f(i+2,j,k);
            e->M.nn[n] = 0.0;
            }
            
            if(p->flag4[IJm2K]<0)
            {
            e->rhsvec.V[n] -= e->M.ee[n]*f(i,j-2,k);
            e->M.ee[n] = 0.0;
            }
            
            if(p->flag4[IJp2K]<0)
            {
            e->rhsvec.V[n] -= e->M.ww[n]*f(i,j+2,k);
            e->M.ww[n] = 0.0;
            }
            
            if(p->flag4[IJKm2]<0)
            {
            e->rhsvec.V[n] -= e->M.bb[n]*f(i,j,k-2);
            e->M.bb[n] = 0.0;
            }
            
            if(p->flag4[IJKp2]<0)
            {
            e->rhsvec.V[n] -= e->M.tt[n]*f(i,j,k+2);
            e->M.tt[n] = 0.0;
            }
        }

	++n;
	}
    
    
    
    double starttime=pgc->timer();
    psolv->start(p,a,pgc,e->Fi,e->rhsvec,6);
    double endtime=pgc->timer();
    pgc->start4(p,e->Fi,250);
    
    p->poissoniter=p->solveriter;
    p->poissontime=endtime-starttime;
	if(p->mpirank==0 && innercounter==p->N50-1 && (p->count%p->P12==0))
	cout<<"Fi_iter: "<<p->solveriter<<"  Fi_time: "<<setprecision(3)<<p->poissontime<<endl;
}*/

