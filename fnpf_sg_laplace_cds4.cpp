/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"fnpf_sg_laplace_cds4.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"solver.h"
#include"ghostcell.h"

fnpf_sg_laplace_cds4::fnpf_sg_laplace_cds4() 
{
}

fnpf_sg_laplace_cds4::~fnpf_sg_laplace_cds4()
{
}

void fnpf_sg_laplace_cds4::start(lexer* p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, fnpf_sg_fsfbc *pf, double *f)
{
    // see p. 1130-1132
    
    double sigxyz2;
	n=0;
    
    // 4th-order
	n=0;
    FLOOP
	{
    sigxyz2 = pow(p->sigx[FIJK],2.0) + pow(p->sigy[FIJK],2.0) + pow(p->sigz[FIJK],2.0);
    
    X0 = -0.5*p->XN[IP2] + 13.0*p->XN[IP1] - 13.0*p->XN[IM1] + 0.5*p->XN[IM2];
    X1 = (-p->XN[IP3] + 27.0*p->XN[IP2] -27.0*p->XN[IP1] + p->XN[IP])*X0;
    X2 = (-p->XN[IP2] + 27.0*p->XN[IP1] -27.0*p->XN[IP] + p->XN[IM1])*X0;
    X3 = (-p->XN[IP1] + 27.0*p->XN[IP] -27.0*p->XN[IM1] + p->XN[IM2])*X0;
    X4 = (-p->XN[IP] + 27.0*p->XN[IM1] -27.0*p->XN[IM2] + p->XN[IM3])*X0;
    
    Y0 = -0.5*p->YN[JP2] + 13.0*p->YN[JP1] - 13.0*p->YN[JM1] + 0.5*p->YN[JM2];
    Y1 = (-p->YN[JP3] + 27.0*p->YN[JP2] -27.0*p->YN[JP1] + p->YN[JP])*Y0;
    Y2 = (-p->YN[JP2] + 27.0*p->YN[JP1] -27.0*p->YN[JP] + p->YN[JM1])*Y0;
    Y3 = (-p->YN[JP1] + 27.0*p->YN[JP] -27.0*p->YN[JM1] + p->YN[JM2])*Y0;
    Y4 = (-p->YN[JP] + 27.0*p->YN[JM1] -27.0*p->YN[JM2] + p->YN[JM3])*Y0;
    
    Z0 = -0.5*p->ZN[KP2] + 13.0*p->ZN[KP1] - 13.0*p->ZN[KM1] + 0.5*p->ZN[KM2];
    Z1 = (-p->ZN[KP3] + 27.0*p->ZN[KP2] -27.0*p->ZN[KP1] + p->ZN[KP])*Z0;
    Z2 = (-p->ZN[KP2] + 27.0*p->ZN[KP1] -27.0*p->ZN[KP] + p->ZN[KM1])*Z0;
    Z3 = (-p->ZN[KP1] + 27.0*p->ZN[KP] -27.0*p->ZN[KM1] + p->ZN[KM2])*Z0;
    Z4 = (-p->ZN[KP] + 27.0*p->ZN[KM1] -27.0*p->ZN[KM2] + p->ZN[KM3])*Z0;
    

	c->M.p[n] = (1.0/X1 + 729.0/X2 + 729.0/X3 + 1.0/X4)*p->x_dir 
             + (1.0/Y1 + 729.0/Y2 + 729.0/Y3 + 1.0/Y4)*p->y_dir 
             + sigxyz2*(1.0/Z1 + 729.0/Z2 + 729.0/Z3 + 1.0/Z4)*p->z_dir;
    
    
   	c->M.n[n] = -(27.0/X1 + 729.0/X2 + 27.0/X3)*p->x_dir;
	c->M.s[n] = -(27.0/X2 + 729.0/X3 + 27.0/X4)*p->x_dir;

	c->M.w[n] = -(27.0/Y1 + 729.0/Y2 + 27.0/Y3)*p->y_dir;
	c->M.e[n] = -(27.0/Y2 + 729.0/Y3 + 27.0/Y4)*p->y_dir;

	c->M.t[n] = -sigxyz2*(27.0/Z1 + 729.0/Z2 + 27.0/Z3 + p->sigxx[FIJK]/(p->DZN[KP]+p->DZN[KM1]))*p->z_dir;
	c->M.b[n] = -sigxyz2*(27.0/Z2 + 729.0/Z3 + 27.0/Z4 - p->sigxx[FIJK]/(p->DZN[KP]+p->DZN[KM1]))*p->z_dir;
    
    
    c->M.nn[n] = (27.0/X1 + 27.0/X2)*p->x_dir; 
    c->M.ss[n] = (27.0/X3 + 27.0/X4)*p->x_dir;
    
    c->M.ww[n] = (27.0/Y1 + 27.0/Y2)*p->y_dir;
    c->M.ee[n] = (27.0/Y3 + 27.0/Y4)*p->y_dir;
     
    c->M.tt[n] = sigxyz2*(27.0/Z1 + 27.0/Z2)*p->z_dir; 
    c->M.bb[n] = sigxyz2*(27.0/Z3 + 27.0/Z4)*p->z_dir;
     
     /*
    c->M.nnn[n] = -(1.0/X1)*p->x_dir; 
    c->M.sss[n] = -(1.0/X4)*p->x_dir;
    
    c->M.www[n] = -(1.0/Y1)*p->y_dir;
    c->M.eee[n] = -(1.0/Y4)*p->y_dir;
     
    c->M.ttt[n] = -sigxyz2*(1.0/Z1)*p->z_dir; 
    c->M.bbb[n] = -sigxyz2*(1.0/Z4)*p->z_dir;*/
    

//cout<<c->M.p[n]<<" "<<c->M.n[n]<<" "<<c->M.s[n]<<" "<<c->M.t[n]<<" "<<c->M.b[n]<<endl;
	//c->rhsvec.V[n] = 0.0;//((f(i+2,j,k)*(-27.0/X1 - 27.0/X2) + f(i-2,j,k)*(-27.0/X3 - 27.0/X4))*p->x_dir
                 // +  (f(i,j+2,k)*(-27.0/Y1 - 27.0/Y2) + f(i,j-2,k)*(-27.0/Y3 - 27.0/Y4))*p->y_dir
                 // +  (f(i,j,k+2)*(-27.0/Z1 - 27.0/Z2) + f(i,j,k-2)*(-27.0/Z3 - 27.0/Z4))*p->z_dir
                  
    c->rhsvec.V[n] =   (f[FIp3JK]*(1.0/X1) + f[FIm3JK]*(1.0/X4))*p->x_dir
                    +  (f[FIJp3K]*(1.0/Y1) + f[FIJm3K]*(1.0/Y4))*p->y_dir
                    +  sigxyz2*(f[FIJKp3]*(1.0/Z1) + f[FIJKm3]*(1.0/Z4))*p->z_dir;
        
    c->rhsvec.V[n] += 2.0*p->sigx[IJK]*(f[FIp1JKp1] - f[FIm1JKp1] - f[FIp1JKm1] + f[FIm1JKm1])
                    /((p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]))*p->x_dir
                    
                    +2.0*p->sigy[IJK]*(f[FIJp1Kp1] - f[FIJm1Kp1] - f[FIJp1Km1] + f[FIJm1Km1])
                    /((p->DYN[JP]+p->DYN[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir;

	++n;
	}
    
    n=0;
	FLOOP
	{
        if(p->flag4[IJK]>0)
        {
    
            if(p->flag4[Im1JK]<0)
            {
            c->rhsvec.V[n] -= c->M.s[n]*f[FIJK];
            c->M.s[n] = 0.0;
            }
            
            if(p->flag4[Ip1JK]<0)
            {
            c->rhsvec.V[n] -= c->M.n[n]*f[FIJK];
            c->M.n[n] = 0.0;
            }
            
            if(p->flag4[IJm1K]<0)
            {
            c->rhsvec.V[n] -= c->M.e[n]*f[FIJK];
            c->M.e[n] = 0.0;
            }
            
            if(p->flag4[IJp1K]<0)
            {
            c->rhsvec.V[n] -= c->M.w[n]*f[FIJK];
            c->M.w[n] = 0.0;
            }
            
            if(p->flag4[IJKm1]<0)
            {
            c->rhsvec.V[n] -= c->M.b[n]*f[FIJKm1];
            c->M.b[n] = 0.0;
            }
            
            if(p->flag4[IJKp1]<0)
            {
            c->rhsvec.V[n] -= c->M.t[n]*f[FIJKp1];
            c->M.t[n] = 0.0;
            }
            
            //--
            if(p->flag4[Im2JK]<0)
            {
            c->rhsvec.V[n] -= c->M.ss[n]*f[FIJK];
            c->M.ss[n] = 0.0;
            }
            
            if(p->flag4[Ip2JK]<0)
            {
            c->rhsvec.V[n] -= c->M.nn[n]*f[FIJK];
            c->M.nn[n] = 0.0;
            }
            
            if(p->flag4[IJm2K]<0)
            {
            c->rhsvec.V[n] -= c->M.ee[n]*f[FIJK];
            c->M.ee[n] = 0.0;
            }
            
            if(p->flag4[IJp2K]<0)
            {
            c->rhsvec.V[n] -= c->M.ww[n]*f[FIJK];
            c->M.ww[n] = 0.0;
            }
            
            if(p->flag4[IJKm2]<0)
            {
            c->rhsvec.V[n] -= c->M.bb[n]*f[FIJKm2];
            c->M.bb[n] = 0.0;
            }
            
            if(p->flag4[IJKp2]<0)
            {
            c->rhsvec.V[n] -= c->M.tt[n]*f[FIJKp2];
            c->M.tt[n] = 0.0;
            }
        }
	++n;
	}
    
    
    double starttime=pgc->timer();
    psolv->startF(p,pgc,f,c->rhsvec,c->M,8,250,p->N44);
    double endtime=pgc->timer();
    
    p->poissoniter=p->solveriter;
    p->poissontime=endtime-starttime;
	if(p->mpirank==0 && innercounter==p->N50-1 && (p->count%p->P12==0))
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
	c->M.p[n] =  30.0/ddx*p->x_dir + 30.0/ddy*p->y_dir + 30.0/ddz*p->z_dir;

   	c->M.n[n] = -16.0/ddx*p->x_dir;
	c->M.s[n] = -16.0/ddx*p->x_dir;

	c->M.w[n] = -16.0/ddy*p->y_dir;
	c->M.e[n] = -16.0/ddy*p->y_dir;

	c->M.t[n] = -16.0/ddz*p->z_dir;
	c->M.b[n] = -16.0/ddz*p->z_dir;
    
    c->M.nn[n] = 1.0/ddx*p->x_dir;
	c->M.ss[n] = 1.0/ddx*p->x_dir;

	c->M.ww[n] = 1.0/ddy*p->y_dir;
	c->M.ee[n] = 1.0/ddy*p->y_dir;

	c->M.tt[n] = 1.0/ddz*p->z_dir;
	c->M.bb[n] = 1.0/ddz*p->z_dir;
    
    //cout<<c->M.p[n]<<" "<<c->M.n[n]<<" "<<c->M.s[n]<<" "<<c->M.t[n]<<" "<<c->M.b[n]<<endl;
	
	c->rhsvec.V[n] = 0.0;//(f(i+2,j,k)/ddx*p->x_dir + f(i-2,j,k)/ddx*p->x_dir
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
            c->rhsvec.V[n] -= c->M.s[n]*f(i-1,j,k);
            c->M.s[n] = 0.0;
            }
            
            if(p->flag4[Ip1JK]<0)
            {
            c->rhsvec.V[n] -= c->M.n[n]*f(i+1,j,k);
            c->M.n[n] = 0.0;
            }
            
            if(p->flag4[IJm1K]<0)
            {
            c->rhsvec.V[n] -= c->M.e[n]*f(i,j-1,k);
            c->M.e[n] = 0.0;
            }
            
            if(p->flag4[IJp1K]<0)
            {
            c->rhsvec.V[n] -= c->M.w[n]*f(i,j+1,k);
            c->M.w[n] = 0.0;
            }
            
            if(p->flag4[IJKm1]<0)
            {
            c->rhsvec.V[n] -= c->M.b[n]*f(i,j,k-1);
            c->M.b[n] = 0.0;
            }
            
            if(p->flag4[IJKp1]<0)
            {
            c->rhsvec.V[n] -= c->M.t[n]*f(i,j,k+1);
            c->M.t[n] = 0.0;
            }
            
            
            //--
            if(p->flag4[Im2JK]<0)
            {
            c->rhsvec.V[n] -= c->M.ss[n]*f(i-2,j,k);
            c->M.ss[n] = 0.0;
            }
            
            if(p->flag4[Ip2JK]<0)
            {
            c->rhsvec.V[n] -= c->M.nn[n]*f(i+2,j,k);
            c->M.nn[n] = 0.0;
            }
            
            if(p->flag4[IJm2K]<0)
            {
            c->rhsvec.V[n] -= c->M.ee[n]*f(i,j-2,k);
            c->M.ee[n] = 0.0;
            }
            
            if(p->flag4[IJp2K]<0)
            {
            c->rhsvec.V[n] -= c->M.ww[n]*f(i,j+2,k);
            c->M.ww[n] = 0.0;
            }
            
            if(p->flag4[IJKm2]<0)
            {
            c->rhsvec.V[n] -= c->M.bb[n]*f(i,j,k-2);
            c->M.bb[n] = 0.0;
            }
            
            if(p->flag4[IJKp2]<0)
            {
            c->rhsvec.V[n] -= c->M.tt[n]*f(i,j,k+2);
            c->M.tt[n] = 0.0;
            }
        }

	++n;
	}
    
    
    
    double starttime=pgc->timer();
    psolv->start(p,a,pgc,c->Fi,c->xvec,c->rhsvec,6,250,p->N44);
    double endtime=pgc->timer();
    pgc->start4(p,a->Fi,250);
    
    p->poissoniter=p->solveriter;
    p->poissontime=endtime-starttime;
	if(p->mpirank==0 && innercounter==p->N50-1)
	cout<<"Fi_iter: "<<p->solveriter<<"  Fi_time: "<<setprecision(3)<<p->poissontime<<endl;
}*/

