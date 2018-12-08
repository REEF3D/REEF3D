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

#include"fnpf_sg_laplace_cds4_v2.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"solver.h"
#include"ghostcell.h"

fnpf_sg_laplace_cds4_v2::fnpf_sg_laplace_cds4_v2() 
{
}

fnpf_sg_laplace_cds4_v2::~fnpf_sg_laplace_cds4_v2()
{
}

void fnpf_sg_laplace_cds4_v2::start(lexer* p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, fnpf_sg_fsfbc *pf, double *f)
{
    // see p. 1130-1132
    
    double sigxyz2;
    double ab,abb,denom;
    double fbxm,fbxp,fbym,fbyp;
	n=0;
    
    // 4th-order
	n=0;
    LOOP
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

	c->M.t[n] = -(sigxyz2*(27.0/Z1 + 729.0/Z2 + 27.0/Z3) + p->sigxx[FIJK]/(p->DZN[KP]+p->DZN[KM1]))*p->z_dir;
	c->M.b[n] = -(sigxyz2*(27.0/Z2 + 729.0/Z3 + 27.0/Z4) - p->sigxx[FIJK]/(p->DZN[KP]+p->DZN[KM1]))*p->z_dir;
    
    
    c->M.nn[n] = (27.0/X1 + 27.0/X2)*p->x_dir; 
    c->M.ss[n] = (27.0/X3 + 27.0/X4)*p->x_dir;
    
    c->M.ww[n] = (27.0/Y1 + 27.0/Y2)*p->y_dir;
    c->M.ee[n] = (27.0/Y3 + 27.0/Y4)*p->y_dir;
     
    c->M.tt[n] = sigxyz2*(27.0/Z1 + 27.0/Z2)*p->z_dir; 
    c->M.bb[n] = sigxyz2*(27.0/Z3 + 27.0/Z4)*p->z_dir;
                  
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
	LOOP
	{
        if(p->flag7[FIJK]>0)
        {
    
            if(p->flag7[FIm1JK]<0)
            {
            c->rhsvec.V[n] -= c->M.s[n]*f[FIJK];
            c->M.s[n] = 0.0;
            }
            
            if(p->flag7[FIp1JK]<0)
            {
            c->rhsvec.V[n] -= c->M.n[n]*f[FIJK];
            c->M.n[n] = 0.0;
            }
            
            if(p->flag7[FIJm1K]<0)
            {
            c->rhsvec.V[n] -= c->M.e[n]*f[FIJK];
            c->M.e[n] = 0.0;
            }
            
            if(p->flag7[FIJp1K]<0)
            {
            c->rhsvec.V[n] -= c->M.w[n]*f[FIJK];
            c->M.w[n] = 0.0;
            }
            
            if(p->flag7[FIJKp2]<0 && p->flag7[FIJKp1]>0)
            {
            c->rhsvec.V[n] -= c->M.t[n]*f[FIJKp2];
            c->M.t[n] = 0.0;
            }
            
            
            //--
            if(p->flag7[FIm2JK]<0)
            {
            c->rhsvec.V[n] -= c->M.ss[n]*f[FIJK];
            c->M.ss[n] = 0.0;
            }
            
            if(p->flag7[FIp2JK]<0)
            {
            c->rhsvec.V[n] -= c->M.nn[n]*f[FIJK];
            c->M.nn[n] = 0.0;
            }
            
            if(p->flag7[FIJm2K]<0)
            {
            c->rhsvec.V[n] -= c->M.ee[n]*f[FIJK];
            c->M.ee[n] = 0.0;
            }
            
            if(p->flag7[FIJp2K]<0)
            {
            c->rhsvec.V[n] -= c->M.ww[n]*f[FIJK];
            c->M.ww[n] = 0.0;
            }
            
            
            if(p->flag7[FIJKp3]<0 && p->flag7[FIJKp2]>0)
            {
            c->rhsvec.V[n] -= c->M.tt[n]*f[FIJKp3];
            c->M.tt[n] = 0.0;
            }
            
            
            // KBEDBC
            if(p->flag7[FIJKm1]<0)
            {
            sigxyz2 = pow(p->sigx[FIJK],2.0) + pow(p->sigy[FIJK],2.0) + pow(p->sigz[FIJK],2.0);
            
            Z0 = -0.5*p->ZN[KP2] + 13.0*p->ZN[KP1] - 13.0*p->ZN[KM1] + 0.5*p->ZN[KM2];
            Z1 = (-p->ZN[KP3] + 27.0*p->ZN[KP2] -27.0*p->ZN[KP1] + p->ZN[KP])*Z0;
            Z2 = (-p->ZN[KP2] + 27.0*p->ZN[KP1] -27.0*p->ZN[KP] + p->ZN[KM1])*Z0;
            Z3 = (-p->ZN[KP1] + 27.0*p->ZN[KP] -27.0*p->ZN[KM1] + p->ZN[KM2])*Z0;
            Z4 = (-p->ZN[KP] + 27.0*p->ZN[KM1] -27.0*p->ZN[KM2] + p->ZN[KM3])*Z0;
            

            ab = - (sigxyz2*(27.0/Z2 + 729.0/Z3 + 27.0/Z4) - p->sigxx[FIJK]/(p->DZN[KP]+p->DZN[KM1]))*p->z_dir;
            abb = sigxyz2*(27.0/Z3 + 27.0/Z4)*p->z_dir;
            denom = p->sigz[FIJK] + c->Bx(i,j)*p->sigx[FIJK] + c->By(i,j)*p->sigy[FIJK];
            
            fbxp = f[FIp1JKp1] + (2.0*p->DZP[KP]*(c->Bx(i+1,j)*((c->Fi[FIp2JK]-c->Fi[FIJK])/(p->DXP[IP] + p->DXP[IM1]))
                                                 +c->By(i+1,j)*((c->Fi[FIp1Jp1K]-c->Fi[FIp1Jm1K])/(p->DYP[JP] + p->DYP[JM1]))))
                                /(p->sigz[FIp1JK] + c->Bx(i+1,j)*p->sigx[FIp1JK] + c->By(i+1,j)*p->sigy[FIp1JK]);
            
            fbxm = f[FIm1JKp1] + (2.0*p->DZP[KP]*(c->Bx(i-1,j)*((c->Fi[FIJK]-c->Fi[FIm2JK])/(p->DXP[IP] + p->DXP[IM1]))
                                                 +c->By(i-1,j)*((c->Fi[FIm1Jp1K]-c->Fi[FIm1Jm1K])/(p->DYP[JP] + p->DYP[JM1]))))
                                /(p->sigz[FIm1JK] + c->Bx(i-1,j)*p->sigx[FIm1JK] + c->By(i-1,j)*p->sigy[FIm1JK]);
            
            fbyp = f[FIJp1Kp1] + (2.0*p->DZP[KP]*(c->Bx(i,j+1)*((c->Fi[FIp1Jp1K]-c->Fi[FIm1Jp1K])/(p->DXP[IP] + p->DXP[IM1]))
                                                 +c->By(i,j+1)*((c->Fi[FIJp2K]-c->Fi[FIJK])/(p->DYP[JP] + p->DYP[JM1]))))
                                /(p->sigz[FIJp1K] + c->Bx(i,j+1)*p->sigx[FIJp1K] + c->By(i,j+1)*p->sigy[FIJp1K]);
            
            fbym = f[FIJm1Kp1] + (2.0*p->DZP[KP]*(c->Bx(i,j-1)*((c->Fi[FIp1Jm1K]-c->Fi[FIm1Jm1K])/(p->DXP[IP] + p->DXP[IM1]))
                                                 +c->By(i,j-1)*((c->Fi[FIJK]-c->Fi[FIJm2K])/(p->DYP[JP] + p->DYP[JM1]))))
                                /(p->sigz[FIJm1K] + c->Bx(i,j-1)*p->sigx[FIJm1K] + c->By(i,j-1)*p->sigy[FIJm1K]);
                                
                            
            
            //--
            c->rhsvec.V[n] -=  2.0*p->sigx[FIJK]*(f[FIp1JKp1] - f[FIm1JKp1] - f[FIp1JKm1] + f[FIm1JKm1])
                    /((p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]))*p->x_dir
                    
                    + 2.0*p->sigy[FIJK]*(f[FIJp1Kp1] - f[FIJm1Kp1] - f[FIJp1Km1] + f[FIJm1Km1])
                    /((p->DYN[JP]+p->DYN[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir;
                    
            c->rhsvec.V[n] +=  2.0*p->sigx[FIJK]*(f[FIp1JKp1] - f[FIm1JKp1] - fbxp + fbxm)
                    /((p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]))*p->x_dir
                    
                    + 2.0*p->sigy[FIJK]*(f[FIJp1Kp1] - f[FIJm1Kp1] - fbyp + fbym)
                    /((p->DYN[JP]+p->DYN[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir
                    
                    +  f[FIJKp3]*(1.0/Z4)*p->z_dir;
                    
        
            
            c->M.n[n] += ab*2.0*p->DZN[KP]*c->Bx(i,j)/(denom*(p->DXP[IP] + p->DXP[IM1]));
            c->M.nn[n] += abb*2.0*p->DZN[KP]*c->Bx(i,j)/(denom*(p->DXP[IP] + p->DXP[IM1]));
            
            c->M.s[n] += -ab*2.0*p->DZN[KP]*c->Bx(i,j)/(denom*(p->DXP[IP] + p->DXP[IM1]));
            c->M.ss[n] += -abb*2.0*p->DZN[KP]*c->Bx(i,j)/(denom*(p->DXP[IP] + p->DXP[IM1]));
            
            c->M.e[n] += ab*2.0*p->DZN[KP]*c->By(i,j)/(denom*(p->DYP[JP] + p->DYP[JM1]));
            c->M.ee[n] += abb*2.0*p->DZN[KP]*c->By(i,j)/(denom*(p->DYP[JP] + p->DYP[JM1]));
            
            c->M.w[n] += -ab*2.0*p->DZN[KP]*c->By(i,j)/(denom*(p->DYP[JP] + p->DYP[JM1]));
            c->M.ww[n] += -abb*2.0*p->DZN[KP]*c->By(i,j)/(denom*(p->DYP[JP] + p->DYP[JM1]));
            
            c->M.t[n] += ab;
            c->M.tt[n] += abb;
            
            c->M.b[n] = 0.0;
            c->M.bb[n] = 0.0;
            }
        }
	++n;
	}
    
    
    double starttime=pgc->timer();
    psolv->startF(p,pgc,f,c->rhsvec,c->M,10,250,p->N44);
    double endtime=pgc->timer();
    
    p->poissoniter=p->solveriter;
    p->poissontime=endtime-starttime;
	if(p->mpirank==0 && innercounter==p->N50-1 && (p->count%p->P12==0))
	cout<<"Fi_iter: "<<p->solveriter<<"  Fi_time: "<<setprecision(3)<<p->poissontime<<endl;
}




