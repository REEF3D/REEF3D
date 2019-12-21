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

#include"fnpf_sg_laplace_cds2_v2.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"solver.h"
#include"ghostcell.h"
#include"fnpf_discrete_weights.h"

fnpf_sg_laplace_cds2_v2::fnpf_sg_laplace_cds2_v2(lexer *p) 
{
    p->Darray(ckx,p->knox+1+4*marge,5);
    p->Darray(cky,p->knoy+1+4*marge,5);
    p->Darray(ckz,p->knoz+1+4*marge,5);

    fnpf_discrete_weights dw(p);
    
    dw.ck_weights(p, ckx, p->XP, p->knox, 2, 2, 1);
    dw.ck_weights(p, cky, p->YP, p->knoy, 2, 2, 2);
    dw.ck_weights(p, ckz, p->ZN, p->knoz, 2, 2, 3);
}


fnpf_sg_laplace_cds2_v2::~fnpf_sg_laplace_cds2_v2()
{
}


void fnpf_sg_laplace_cds2_v2::start(lexer* p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, fnpf_sg_fsf *pf, double *f)
{
    // see p. 1130-1132
    double sigxyz2;
    double ab,abb,abbb,denom;
    double fbxm,fbxp,fbym,fbyp;
    double distfac,dist;
    double xdelta,ydelta;    
    
	n=0;
    LOOP
    {
        if(c->wet(i,j)==1)
        {
        sigxyz2 = pow(p->sigx[FIJK],2.0) + pow(p->sigy[FIJK],2.0) + pow(p->sigz[IJ],2.0);
        
        c->M.p[n] = ckx[IP][1]*p->x_dir 
                  + cky[JP][1]*p->y_dir 
                  + sigxyz2*ckz[KP][1]*p->z_dir; 
        
        c->M.n[n] = ckx[IP][2]*p->x_dir; 
        c->M.s[n] = ckx[IP][0]*p->x_dir;

        c->M.w[n] = cky[JP][2]*p->y_dir; 
        c->M.e[n] = cky[JP][0]*p->y_dir; 

        c->M.t[n] = (sigxyz2*ckz[KP][2]  - p->sigxx[FIJK]/(p->DZN[KP]+p->DZN[KM1]))*p->z_dir;
        c->M.b[n] = (sigxyz2*ckz[KP][0]  + p->sigxx[FIJK]/(p->DZN[KP]+p->DZN[KM1]))*p->z_dir;
        
        
       
        c->rhsvec.V[n] =  2.0*p->sigx[FIJK]*(f[FIp1JKp1] - f[FIm1JKp1] - f[FIp1JKm1] + f[FIm1JKm1])
                        /((p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]))*p->x_dir
                        
                        + 2.0*p->sigy[FIJK]*(f[FIJp1Kp1] - f[FIJm1Kp1] - f[FIJp1Km1] + f[FIJm1Km1])
                        /((p->DYN[JP]+p->DYN[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir;
        }    
    ++n;
    }
    
    n=0;
    LOOP
	{
        if(c->wet(i,j)==0)
        {
        c->M.p[n]  =  1.0;


        c->M.n[n] = 0.0;
        c->M.s[n] = 0.0;

        c->M.w[n] = 0.0;
        c->M.e[n] = 0.0;

        c->M.t[n] = 0.0;
        c->M.b[n] = 0.0;
        
        c->rhsvec.V[n] =  0.0;
        }
	++n;
	}
    
    
    n=0;
	LOOP
    if(p->flag7[FIJK]>0)
	{
            if(c->wet(i,j)==1)
            {
        
            // south
            if(p->flag7[FIm1JK]<0 && c->wet(i-1,j)==1 && c->bc(i-1,j)==0)
            {
            c->M.p[n] += c->M.s[n];
            c->M.s[n] = 0.0;
            }
            
            if(c->wet(i-1,j)==0)
            {
            c->M.p[n] += -1.0/(p->DXP[IM1]*p->DXN[IM1])*p->x_dir;
            c->M.s[n] = 0.0;
            }
            
            if(p->flag7[FIm1JK]<0 && c->bc(i-1,j)==1)
            {
            c->rhsvec.V[n] += c->M.s[n]*c->Uin[FIm1JK]*p->DXP[IM1];
            c->M.p[n] += c->M.s[n];
            c->M.s[n] = 0.0;
            }
            
            // north
            if(p->flag7[FIp1JK]<0 && c->wet(i+1,j)==1 && c->bc(i+1,j)==0)
            {
            c->M.p[n] += c->M.n[n];
            c->M.n[n] = 0.0;
            }
            
            if(c->wet(i+1,j)==0)
            {
            c->M.p[n] += c->M.n[n];
            c->M.n[n] = 0.0;
            }
            
            if(p->flag7[FIp1JK]<0 && c->wet(i+1,j)==1 &&  c->bc(i+1,j)==2)
            {
            c->rhsvec.V[n] -= c->M.n[n]*c->Uin[FIp1JK]*p->DXP[IP1];
            c->M.p[n] += c->M.n[n];
            c->M.n[n] = 0.0;
            }

        
            // east
            if(p->flag7[FIJm1K]<0 && c->wet(i,j-1)==1)
            {
            c->M.p[n] += c->M.e[n];
            c->M.e[n] = 0.0;
            }
            
            if(c->wet(i,j-1)==0)
            {
            c->M.p[n] += c->M.e[n];
            c->M.e[n] = 0.0;
            }
            
            
            // west
            if(p->flag7[FIJp1K]<0 && c->wet(i,j+1)==1)
            {
            c->M.p[n] += c->M.w[n];
            c->M.w[n] = 0.0;
            }
            
            if(c->wet(i,j+1)==0)
            {
            c->M.p[n] += c->M.w[n];
            c->M.w[n] = 0.0;
            }
    
            
            // top
            if(p->flag7[FIJKp2]<0 && p->flag7[FIJKp1]>0)
            {
            c->rhsvec.V[n] -= c->M.t[n]*f[FIJKp2];
            c->M.t[n] = 0.0;
            }
 
            // KBEDBC
            if(p->flag7[FIJKm1]<0)
            {
            sigxyz2= pow(p->sigx[FIJK],2.0) + pow(p->sigy[FIJK],2.0) + pow(p->sigz[IJ],2.0);
            
            ab = c->M.b[n];//-(sigxyz2/(p->DZP[KM1]*p->DZN[KM1]) - p->sigxx[FIJK]/(p->DZN[KP]+p->DZN[KM1]));
            
            denom = p->sigz[IJ] + c->Bx(i,j)*p->sigx[FIJK] + c->By(i,j)*p->sigy[FIJK];

                    if(c->wet(i+1,j)==1 && c->wet(i-1,j)==1)
                    {
                    c->M.n[n] +=  ab*2.0*p->DZN[KP]*c->Bx(i,j)/(denom*(p->DXP[IP] + p->DXP[IM1]));
                    c->M.s[n] += -ab*2.0*p->DZN[KP]*c->Bx(i,j)/(denom*(p->DXP[IP] + p->DXP[IM1]));
                    }
                    
                    if(c->wet(i,j-1)==1 && c->wet(i,j+1)==1)
                    {
                    c->M.e[n] +=  ab*2.0*p->DZN[KP]*c->By(i,j)/(denom*(p->DYP[JP] + p->DYP[JM1]));
                    c->M.w[n] += -ab*2.0*p->DZN[KP]*c->By(i,j)/(denom*(p->DYP[JP] + p->DYP[JM1]));
                    }
                
                
                c->M.t[n] += ab;
                c->M.b[n] = 0.0;

            }
            }
	++n;
	}
    
    double starttime=pgc->timer();
    psolv->startF(p,c,pgc,f,c->rhsvec,c->M,8,p->N46,p->N44);
    double endtime=pgc->timer();
    
    p->poissoniter=p->solveriter;
    p->poissontime=endtime-starttime;
    
	if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"Fi_iter: "<<p->poissoniter<<" Final_residual: "<<p->final_res<<"  Fi_time: "<<setprecision(3)<<p->poissontime<<endl;
}