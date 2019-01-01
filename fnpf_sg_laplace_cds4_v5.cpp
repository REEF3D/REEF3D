/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#include"fnpf_sg_laplace_cds4_v5.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"solver.h"
#include"ghostcell.h"
#include"fnpf_discrete_weights.h"
#include"fnpf_sg_bed_update.h"

fnpf_sg_laplace_cds4_v5::fnpf_sg_laplace_cds4_v5(lexer *p) 
{
    p->Darray(ckx,p->knox+1+4*marge,5);
    p->Darray(cky,p->knoy+1+4*marge,5);
    p->Darray(ckz,p->knoz+1+4*marge,5);

    fnpf_discrete_weights dw(p);
    
    dw.ck_weights(p, ckx, p->XP, p->knox, 2, 4, 1);
    dw.ck_weights(p, cky, p->YP, p->knoy, 2, 4, 2);
    dw.ck_weights(p, ckz, p->ZN, p->knoz, 2, 4, 3);
    
    bc=2;
    
    pbed = new fnpf_sg_bed_update(p);
}


fnpf_sg_laplace_cds4_v5::~fnpf_sg_laplace_cds4_v5()
{
}

void fnpf_sg_laplace_cds4_v5::start(lexer* p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, fnpf_sg_fsfbc *pf, double *f)
{
    // see p. 1130-1132
    
    double sigxyz2;
    double ab,abb,abbb,denom;
    double fbxm,fbxp,fbym,fbyp;
    double distfac,dist;
    double zdelta;
    int iter; 
	n=0;
    
    p->N46=1;
    
    // 4th-order
    
    iter=0;
    do
    {
    pbed->bedbc_sig(p,c,pgc,c->Fi,pf);
        n=0;
        LOOP
        {
            sigxyz2 = pow(p->sigx[FIJK],2.0) + pow(p->sigy[FIJK],2.0) + pow(p->sigz[FIJK],2.0);
            
            zdelta = -p->ZN[KP2] + 8.0*p->ZN[KP1] - 8.0*p->ZN[KM1] + p->ZN[KM2];
            
            c->M.p[n] = ckx[IP][2]*p->x_dir 
                      + cky[JP][2]*p->y_dir 
                      + sigxyz2*ckz[KP][2]*p->z_dir; 
            
            c->M.n[n] = ckx[IP][3]*p->x_dir; 
            c->M.s[n] = ckx[IP][1]*p->x_dir;

            c->M.w[n] = cky[JP][3]*p->y_dir; 
            c->M.e[n] = cky[JP][1]*p->y_dir; 

            c->M.t[n] = (sigxyz2*ckz[KP][3]  - 8.0*p->sigxx[FIJK]/zdelta)*p->z_dir;
            c->M.b[n] = (sigxyz2*ckz[KP][1]  + 8.0*p->sigxx[FIJK]/zdelta)*p->z_dir;
            
            
            c->M.nn[n] = ckx[IP][4]*p->x_dir;
            c->M.ss[n] = ckx[IP][0]*p->x_dir;
            
            c->M.ww[n] = cky[JP][4]*p->y_dir;
            c->M.ee[n] = cky[JP][0]*p->y_dir;
             
            c->M.tt[n] = (sigxyz2*ckz[KP][4] + p->sigxx[FIJK]/zdelta)*p->z_dir; 
            c->M.bb[n] = (sigxyz2*ckz[KP][0] - p->sigxx[FIJK]/zdelta)*p->z_dir; 
            
            /*
            c->M.nt[n] = -2.0*p->sigx[FIJK]/((p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]))*p->x_dir;
            c->M.nb[n] = 2.0*p->sigx[FIJK]/((p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]))*p->x_dir;
            c->M.st[n] = 2.0*p->sigx[FIJK]/((p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]))*p->x_dir;
            c->M.sb[n] = -2.0*p->sigx[FIJK]/((p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]))*p->x_dir;
            
            c->M.wt[n] = -2.0*p->sigy[IJK]/((p->DYN[JP]+p->DYN[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir;
            c->M.wb[n] = 2.0*p->sigy[IJK]/((p->DYN[JP]+p->DYN[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir;
            c->M.et[n] = 2.0*p->sigy[IJK]/((p->DYN[JP]+p->DYN[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir;
            c->M.eb[n] = -2.0*p->sigy[IJK]/((p->DYN[JP]+p->DYN[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir;*/
           
            c->rhsvec.V[n] = 2.0*p->sigx[IJK]*(f[FIp1JKp1] - f[FIm1JKp1] - f[FIp1JKm1] + f[FIm1JKm1])
                            /((p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]))*p->x_dir
                            
                            +2.0*p->sigy[IJK]*(f[FIJp1Kp1] - f[FIJm1Kp1] - f[FIJp1Km1] + f[FIJm1Km1])
                            /((p->DYN[JP]+p->DYN[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir;
                    
                    // nb, nt, sb, st
                    // wb, wt, eb, et
            ++n;
            }
        
        n=0;
        LOOP
        if(p->flag7[FIJK]>0)
        {
                
                if(p->flag7[FIm1JK]<0)
                {
                c->M.p[n] += c->M.s[n];
                c->M.s[n] = 0.0;
                }
                
                if(p->flag7[FIp1JK]<0)
                {
                c->M.p[n] += c->M.n[n];
                c->M.n[n] = 0.0;
                }
                
                if(p->flag7[FIJm1K]<0)
                {
                c->M.p[n] += c->M.e[n];
                c->M.e[n] = 0.0;
                }
                
                if(p->flag7[FIJp1K]<0)
                {
                c->M.p[n] += c->M.w[n];
                c->M.w[n] = 0.0;
                }
                
                if(p->flag7[FIJKp2]<0)
                {
                c->rhsvec.V[n] -= c->M.t[n]*f[FIJKp1];
                c->M.t[n] = 0.0;
                }
                
                
               //--
                
                if(p->flag7[FIm2JK]<0)
                {
                c->M.p[n] += c->M.ss[n];
                c->M.ss[n] = 0.0;
                }
                
                if(p->flag7[FIp2JK]<0)
                {
                c->M.p[n] += c->M.nn[n];
                c->M.nn[n] = 0.0;
                }
                
                if(p->flag7[FIJm2K]<0)
                {
                c->M.p[n] += c->M.ee[n];
                c->M.ee[n] = 0.0;
                }
                
                if(p->flag7[FIJp2K]<0)
                {
                c->M.p[n] += c->M.ww[n];
                c->M.ww[n] = 0.0;
                }

                if(p->flag7[FIJKp3]<0)
                {
                c->rhsvec.V[n] -= c->M.tt[n]*f[FIJKp3];
                c->M.tt[n] = 0.0;
                }
                
                
        
            // KBEDBC
            if(p->flag7[FIJKm1]<0 || p->flag7[FIJKm2]<0 )
            {
                ab  = c->M.b[n];
                abb = c->M.bb[n];
                
                // ------
                if(p->flag7[FIJKm2]<0 && p->flag7[FIJKm1]>0)
                {
                // bb
                denom = p->sigz[FIJK] + c->Bx(i,j)*p->sigx[FIJK] + c->By(i,j)*p->sigy[FIJK];
                
                double xdelta = (-p->XP[IP2] + 8.0*p->XP[IP1] - 8.0*p->XP[IM1] + p->XP[IM2]);
                double ydelta = (-p->YP[JP2] + 8.0*p->YP[JP1] - 8.0*p->YP[JM1] + p->YP[JM2]);  
                
                dist = 2.0*(p->DZN[KP]);
                
                c->M.n[n] += 8.0*abb*dist*c->Bx(i,j)/(denom*xdelta); 
                c->M.nn[n] += -abb*dist*c->Bx(i,j)/(denom*xdelta);     
                
                c->M.s[n] += -8.0*abb*dist*c->Bx(i,j)/(denom*xdelta); 
                c->M.ss[n] += abb*dist*c->Bx(i,j)/(denom*xdelta); 
                
                c->M.e[n] += -8.0*abb*dist*c->By(i,j)/(denom*ydelta);
                c->M.ee[n] += abb*dist*c->By(i,j)/(denom*ydelta);
                
                c->M.w[n] += 8.0*abb*dist*c->By(i,j)/(denom*ydelta);
                c->M.ww[n] += -abb*dist*c->By(i,j)/(denom*ydelta);
                
                c->M.p[n] += abb;
                c->M.bb[n] = 0.0;
                }


                // ------
                if(p->flag7[FIJKm2]<0 && p->flag7[FIJKm1]<0 && bc==2)
                {
                denom = p->sigz[FIJK] + c->Bx(i,j)*p->sigx[FIJK] + c->By(i,j)*p->sigy[FIJK];
                
                // b
                fbxp = f[FIp1JKp1] + (2.0*p->DZN[KP]*(c->Bx(i+1,j)*((c->Fi[FIp2JK]-c->Fi[FIJK])/(p->DXP[IP] + p->DXP[IM1]))
                                                     +c->By(i+1,j)*((c->Fi[FIp1Jp1K]-c->Fi[FIp1Jm1K])/(p->DYP[JP] + p->DYP[JM1]))))
                                    /(p->sigz[FIp1JK] + c->Bx(i+1,j)*p->sigx[FIp1JK] + c->By(i+1,j)*p->sigy[FIp1JK]);
                
                fbxm = f[FIm1JKp1] + (2.0*p->DZN[KP]*(c->Bx(i-1,j)*((c->Fi[FIJK]-c->Fi[FIm2JK])/(p->DXP[IP] + p->DXP[IM1]))
                                                     +c->By(i-1,j)*((c->Fi[FIm1Jp1K]-c->Fi[FIm1Jm1K])/(p->DYP[JP] + p->DYP[JM1]))))
                                    /(p->sigz[FIm1JK] + c->Bx(i-1,j)*p->sigx[FIm1JK] + c->By(i-1,j)*p->sigy[FIm1JK]);
                
                fbyp = f[FIJp1Kp1] + (2.0*p->DZN[KP]*(c->Bx(i,j+1)*((c->Fi[FIp1Jp1K]-c->Fi[FIm1Jp1K])/(p->DXP[IP] + p->DXP[IM1]))
                                                     +c->By(i,j+1)*((c->Fi[FIJp2K]-c->Fi[FIJK])/(p->DYP[JP] + p->DYP[JM1]))))
                                    /(p->sigz[FIJp1K] + c->Bx(i,j+1)*p->sigx[FIJp1K] + c->By(i,j+1)*p->sigy[FIJp1K]);
                
                fbym = f[FIJm1Kp1] + (2.0*p->DZN[KP]*(c->Bx(i,j-1)*((c->Fi[FIp1Jm1K]-c->Fi[FIm1Jm1K])/(p->DXP[IP] + p->DXP[IM1]))
                                                     +c->By(i,j-1)*((c->Fi[FIJK]-c->Fi[FIJm2K])/(p->DYP[JP] + p->DYP[JM1]))))
                                    /(p->sigz[FIJm1K] + c->Bx(i,j-1)*p->sigx[FIJm1K] + c->By(i,j-1)*p->sigy[FIJm1K]);
                                                     
                
                //--
               /* c->rhsvec.V[n] -= 2.0*p->sigx[IJK]*(f[FIp1JKp1] - f[FIm1JKp1] - f[FIp1JKm1] + f[FIm1JKm1])
                            /((p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]))*p->x_dir
                            
                            +2.0*p->sigy[IJK]*(f[FIJp1Kp1] - f[FIJm1Kp1] - f[FIJp1Km1] + f[FIJm1Km1])
                            /((p->DYN[JP]+p->DYN[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir;
                        
                c->rhsvec.V[n] +=  2.0*p->sigx[FIJK]*(f[FIp1JKp1] - f[FIm1JKp1] - fbxp + fbxm)
                        /((p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]))*p->x_dir
                        
                        + 2.0*p->sigy[FIJK]*(f[FIJp1Kp1] - f[FIJm1Kp1] - fbyp + fbym)
                        /((p->DYN[JP]+p->DYN[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir;*/
                        
        
                double xdelta = (-p->XP[IP2] + 8.0*p->XP[IP1] - 8.0*p->XP[IM1] + p->XP[IM2]);
                double ydelta = (-p->YP[JP2] + 8.0*p->YP[JP1] - 8.0*p->YP[JM1] + p->YP[JM2]);  
                
                dist = 2.0*(p->DZN[KP]);
                
                c->M.n[n] += 8.0*ab*dist*c->Bx(i,j)/(denom*xdelta); 
                c->M.nn[n] += -ab*dist*c->Bx(i,j)/(denom*xdelta);     
                
                c->M.s[n] += -8.0*ab*dist*c->Bx(i,j)/(denom*xdelta); 
                c->M.ss[n] += ab*dist*c->Bx(i,j)/(denom*xdelta); 

       
                 c->M.e[n] += -8.0*ab*dist*c->By(i,j)/(denom*ydelta);
                c->M.ee[n] += ab*dist*c->By(i,j)/(denom*ydelta);
                
                c->M.w[n] += 8.0*ab*dist*c->By(i,j)/(denom*ydelta);
                c->M.ww[n] += -ab*dist*c->By(i,j)/(denom*ydelta);
      
                c->M.t[n] += ab;
                c->M.b[n] = 0.0;
                }
                
                
                
                if(p->flag7[FIJKm2]<0 && p->flag7[FIJKm1]<0 && bc==4)
                {
                denom = p->sigz[FIJK] + c->Bx(i,j)*p->sigx[FIJK] + c->By(i,j)*p->sigy[FIJK];
                
                // b
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
                        /((p->DYN[JP]+p->DYN[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir;
                        
                dist = -(1.0/3.0)*p->ZN[KM1] - (1.0/2.0)*p->ZN[KP] + p->ZN[KP1] - (5.0/3.0)*p->ZN[KP2];
                
                c->M.n[n] += ab*3.0*dist*c->Bx(i,j)/(denom*(p->DXP[IP] + p->DXP[IM1]));     
                c->M.s[n] += -ab*3.0*dist*c->Bx(i,j)/(denom*(p->DXP[IP] + p->DXP[IM1]));    
                c->M.e[n] += ab*3.0*dist*c->By(i,j)/(denom*(p->DYP[JP] + p->DYP[JM1]));
                c->M.w[n] += -ab*3.0*dist*c->By(i,j)/(denom*(p->DYP[JP] + p->DYP[JM1]));  
                c->M.p[n] += -(3.0/2.0)*ab;
                c->M.t[n] += 3.0*ab;
                c->M.tt[n]+= -5.0*ab;
                c->M.b[n] = 0.0;
                }
            }
        ++n;
        }
    ++iter;
    }
    while(iter<20);
    


    double starttime=pgc->timer();
    psolv->startF(p,pgc,f,c->rhsvec,c->M,10,250,p->N44);
    double endtime=pgc->timer();
    
    p->poissoniter=p->solveriter;
    p->poissontime=endtime-starttime;
    
	if(p->mpirank==0 && innercounter==p->N50-1 && (p->count%p->P12==0))
	cout<<"Fi_iter: "<<p->poissoniter<<" Final_residual: "<<p->final_res<<"  Fi_time: "<<setprecision(3)<<p->poissontime<<endl;
}
