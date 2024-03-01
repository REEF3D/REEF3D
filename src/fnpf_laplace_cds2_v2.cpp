/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"fnpf_laplace_cds2_v2.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"solver_fnpf.h"
#include"hypre_struct_fnpf.h"
#include"hypre_sstruct_fnpf.h"
#include"ghostcell.h"
#include"fnpf_discrete_weights.h"

fnpf_laplace_cds2_v2::fnpf_laplace_cds2_v2(lexer *p, ghostcell *pgc) 
{
    if(p->N10>=11 && p->N10<=19)
    psolv = new hypre_struct_fnpf(p,pgc,p->N10,p->N11);
    
    if(p->N10>=31 && p->N10<=39)
    psolv = new hypre_sstruct_fnpf(p,pgc,p->N10,p->N11);
    
    
    vecsize=p->knox*p->knoy*(p->knoz+1); 
    
    p->Darray(M,vecsize*15);
    p->Darray(x,vecsize);
    p->Darray(rhs,vecsize);

}

fnpf_laplace_cds2_v2::~fnpf_laplace_cds2_v2()
{
}

void fnpf_laplace_cds2_v2::start(lexer* p, fdm_fnpf *c, ghostcell *pgc, solver *psolv_reg, fnpf_fsf *pf, double *f, slice &Fifsf)
{
    if(p->j_dir==0)
    laplace2D(p, c, pgc, psolv_reg, pf, f);
    
    if(p->j_dir==1)
    laplace3D(p, c, pgc, psolv_reg, pf, f);
}

void fnpf_laplace_cds2_v2::laplace2D(lexer* p, fdm_fnpf *c, ghostcell *pgc, solver *psolv_reg, fnpf_fsf *pf, double *f)
{
    double sigxyz2;
    double ab,denom;
    double fbxm,fbxp,fbym,fbyp;
    p->poissoniter=0;
    p->poissontime=0.0;
    
/*
 * {{0,0}, {-1,0}, {1,0},  {0,-1}, {0,1}, {-1,-1},{-1,1},{1,-1},{1,1}};
p  0
s  1
n  2
b  3
t  4
sb 5
st 6
nb 7
nt 8
*/

	n=0;
    KJILOOP
	{
        if(p->flag7[FIJK]>0 && p->wet[IJ]==1)
        {
        sigxyz2 = pow(p->sigx[FIJK],2.0) + pow(p->sigy[FIJK],2.0) + pow(p->sigz[IJ],2.0);
        
        M[n*9]  =     1.0/(p->DXP[IP]*p->DXN[IP])*p->x_dir 
                    + 1.0/(p->DXP[IM1]*p->DXN[IP])*p->x_dir 
                     
                    + (sigxyz2/(p->DZP[KM1]*p->DZN[KP]))*p->z_dir
                    + (sigxyz2/(p->DZP[KM1]*p->DZN[KM1]))*p->z_dir;

        
        M[n*9+1] = -1.0/(p->DXP[IM1]*p->DXN[IP])*p->x_dir;
        M[n*9+2] = -1.0/(p->DXP[IP]*p->DXN[IP])*p->x_dir;
        
        M[n*9+3] = -(sigxyz2/(p->DZP[KM1]*p->DZN[KM1]) - p->sigxx[FIJK]/(p->DZN[KP]+p->DZN[KM1]))*p->z_dir;
        M[n*9+4] = -(sigxyz2/(p->DZP[KM1]*p->DZN[KP])  + p->sigxx[FIJK]/(p->DZN[KP]+p->DZN[KM1]))*p->z_dir;
      
        M[n*9+5]  = -2.0*p->sigx[FIJK]/((p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]))*p->x_dir;
        M[n*9+6]  =  2.0*p->sigx[FIJK]/((p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]))*p->x_dir;
        M[n*9+7]  =  2.0*p->sigx[FIJK]/((p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]))*p->x_dir;
        M[n*9+8]  = -2.0*p->sigx[FIJK]/((p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]))*p->x_dir;

        x[n] = f[FIJK];
        
        rhs[n] =  0.0;
                        
        }
        
        if(p->wet[IJ]==0 || p->flag7[FIJK]<0)
        {
        M[n*9]  =  1.0;
        M[n*9+1] = 0.0;
        M[n*9+2] = 0.0;
        M[n*9+3] = 0.0;
        M[n*9+4] = 0.0;
        M[n*9+5] = 0.0;
        M[n*9+6] = 0.0;
        M[n*9+7] = 0.0;
        M[n*9+8] = 0.0;
        
        x[n] = 0.0;
        rhs[n] =  0.0;
        }
        
	++n;
	}
    
    n=0;
	KJILOOP
	{
            if(p->flag7[FIJK]>0 && p->wet[IJ]==1)
            {
            // south
            if((p->flag7[FIm1JK]<0 || p->wet[Im1J]==0) && c->bc(i-1,j)==0)
            {
            M[n*9] += M[n*9+1];  
            M[n*9+1] = 0.0;        
            }
            
            if(p->flag7[FIm1JK]<0 && c->bc(i-1,j)==1)
            {
            rhs[n] += M[n*9+1]*c->Uin[FIm1JK]*p->DXP[IM1];
            M[n*9] += M[n*9+1];
            M[n*9+1] = 0.0;
            }
            
            // north
            if((p->flag7[FIp1JK]<0 || p->wet[Ip1J]==0) && c->bc(i+1,j)==0)
            {
            M[n*9] += M[n*9+2];
            M[n*9+2] = 0.0;
            }
            
            if(p->flag7[FIp1JK]<0 && c->bc(i+1,j)==2)
            {
            rhs[n] -= M[n*9+2]*c->Uin[FIp1JK]*p->DXP[IP1];
            M[n*9] += M[n*9+2];
            M[n*9+2] = 0.0;
            }

            // top
            if(p->flag7[FIJKp2]<0 && p->flag7[FIJKp1]>0)
            {
            rhs[n] -= M[n*9+4]*f[FIJKp2];
            M[n*9+4] = 0.0;
            }
   
        // diagonal entries
/*
p  0
s  1
n  2
b  3
t  4
sb 5
st 6
nb 7
nt 8
*/
            // st
                // fsfbc
            if(p->flag7[FIm1JKp2]<0 && p->flag7[FIJKp2]<0 && p->flag7[FIJKp1]>0  && c->bc(i-1,j)==0) // fsfbc
            {
            rhs[n] -= M[n*9+6]*f[FIm1JKp2];
            M[n*9+6] = 0.0;
            }
                // wall
            if((p->flag7[FIm1JKp1]<0 && p->flag7[FIJKp2]>0) && c->bc(i-1,j)==0) //
            {
            M[n*9] += M[n*9+6];  
            M[n*9+6] = 0.0;   //cout<<p->mpirank<<" ST i: "<<i<<" k: "<<k<<endl;  
            }
            
                // wavegen
            if(p->flag7[FIm1JKp1]<0 && c->bc(i-1,j)==1)
            {
            rhs[n] += M[n*9+6]*c->Uin[FIm1JKp1]*p->DXP[IM1]/(p->DZN[KP]+p->DZN[KM1]);
            M[n*9+8] += M[n*9+6];
            M[n*9+6] = 0.0;
            }
            
            // nt
                // fsfbc
            if(p->flag7[FIp1JKp2]<0 && p->flag7[FIJKp2]<0 && p->flag7[FIJKp1]>0  && c->bc(i+1,j)==0) 
            {
            rhs[n] -= M[n*9+8]*f[FIp1JKp2];
            M[n*9+8] = 0.0;
            }
            
                // wall
            if(p->flag7[FIp1JKp1]<0 && p->flag7[FIJKp2]>0  && c->bc(i+1,j)==0)
            {
            M[n*9] += M[n*9+8];
            M[n*9+8] = 0.0;
            }
            
                // AWA
            if(p->flag7[FIp1JKp1]<0 && c->bc(i+1,j)==2)
            {
            rhs[n] -= M[n*9+8]*c->Uin[FIp1JKp1]*p->DXP[IP1]/(p->DZN[KP]+p->DZN[KM1]);
            M[n*9] += M[n*9+8];
            M[n*9+8] = 0.0;
            }
            
            // sb -----
                // wall
            if(((p->flag7[FIm1JKm1]<0 && p->flag7[FIJKm1]>0)|| p->wet[Im1J]==0) && c->bc(i-1,j)==0)
            {
            M[n*9] += M[n*9+5];  
            M[n*9+5] = 0.0;        
            }
            
                // wavegen
            if((p->flag7[FIm1JKm1]<0 && p->flag7[FIJKm1]>0) && c->bc(i-1,j)==1)
            {
            rhs[n] += M[n*9+5]*c->Uin[FIm1JKm1]*p->DXP[IM1]/(p->DZN[KP]+p->DZN[KM1]);
            M[n*9+7] += M[n*9+5];
            M[n*9+5] = 0.0;
            //cout<<p->mpirank<<" SB i: "<<i<<" k: "<<k<<endl;   
            }
            
            // nb 
                // wall
            if(((p->flag7[FIp1JKm1]<0 && p->flag7[FIJKm1]>0)|| p->wet[Ip1J]==0) && c->bc(i+1,j)==0)
            {
            M[n*9] += M[n*9+7];  
            M[n*9+7] = 0.0;        
            }
            
                // wavegen
            if((p->flag7[FIp1JKm1]<0 && p->flag7[FIJKm1]>0) && c->bc(i+1,j)==2)
            {
            rhs[n] -= M[n*9+7]*c->Uin[FIp1JKm1]*p->DXP[IP1]/(p->DZN[KP]+p->DZN[KM1]);
            M[n*9] += M[n*9+7];
            M[n*9+7] = 0.0;
            }
                
            // sb KBEDBC
            if(p->flag7[FIm1JKm1]<0 && p->flag7[FIJKm1]<0 && c->bc(i-1,j)==0)
            {
            ab = -2.0*p->sigx[FIJK]/((p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]))*p->x_dir;
            
            denom = p->sigz[IJ] + c->Bx(i-1,j)*p->sigx[FIm1JK];

                    if(p->wet[Ip1J]==1 && p->wet[Im1J]==1)
                    {
                    M[n*9+2] +=  ab*2.0*p->DZN[KP]*c->Bx(i,j)/(denom*(p->DXP[IP] + p->DXP[IM1]));
                    M[n*9+1] += -ab*2.0*p->DZN[KP]*c->Bx(i,j)/(denom*(p->DXP[IP] + p->DXP[IM1]));
                    }
                
                M[n*9+6] += ab;
                M[n*9+5] = 0.0;
            }
            
            // nb KBEDBC
            if(p->flag7[FIp1JKm1]<0 && p->flag7[FIJKm1]<0  && c->bc(i+1,j)==0)
            {
            ab = 2.0*p->sigx[FIJK]/((p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]))*p->x_dir;
            
            denom = p->sigz[IJ] + c->Bx(i,j)*p->sigx[FIJK];

                    if(p->wet[Ip1J]==1 && p->wet[Im1J]==1)
                    {
                    M[n*9+2] +=  ab*2.0*p->DZN[KP]*c->Bx(i,j)/(denom*(p->DXP[IP] + p->DXP[IM1]));
                    M[n*9+1] += -ab*2.0*p->DZN[KP]*c->Bx(i,j)/(denom*(p->DXP[IP] + p->DXP[IM1]));
                    }
                
                M[n*9+8] += ab;
                M[n*9+7] = 0.0;
            }
 
            // KBEDBC
            if(p->flag7[FIJKm1]<0)
            {
            sigxyz2 = pow(p->sigx[FIJK],2.0) + pow(p->sigy[FIJK],2.0) + pow(p->sigz[IJ],2.0);
            
            ab = -(sigxyz2/(p->DZP[KM1]*p->DZN[KM1]) - p->sigxx[FIJK]/(p->DZN[KP]+p->DZN[KM1]));
            
            denom = p->sigz[IJ] + c->Bx(i,j)*p->sigx[FIJK];

                    if(p->wet[Ip1J]==1 && p->wet[Im1J]==1)
                    {
                    M[n*9+2] +=  ab*2.0*p->DZN[KP]*c->Bx(i,j)/(denom*(p->DXP[IP] + p->DXP[IM1]));
                    M[n*9+1] += -ab*2.0*p->DZN[KP]*c->Bx(i,j)/(denom*(p->DXP[IP] + p->DXP[IM1]));
                    }
                    
                M[n*9+4] += ab;
                M[n*9+3] = 0.0;

            }
            }
	++n;
	}

    double starttime=pgc->timer();
    psolv->startF(p,pgc,x,rhs,M,8);
    double endtime=pgc->timer();
    
        n=0;
        KJILOOP
        {
		 FPWDCHECK
        f[FIJK]=x[n];
		
        ++n;
        }
    
    p->poissoniter+=p->solveriter;
    p->poissontime+=endtime-starttime;
    
    
	if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"Fi_iter: "<<p->poissoniter<<" Final_residual: "<<p->final_res<<"  Fi_time: "<<setprecision(3)<<p->poissontime<<endl;
}
/*
 * {{0,0}, {-1,0}, {1,0},  {0,-1}, {0,1}, {-1,-1},{-1,1},{1,-1},{1,1}};
p  0
s  1
n  2
b  3
t  4
sb 5
st 6
nb 7
nt 8

//-----------------

p  0
s  1
n  2
e  3
w  4
b  5
t  6
sb 7
st 8
nb 9
nt 10
eb 11
et 12
wb 13
wt 14
*/

void fnpf_laplace_cds2_v2::laplace3D(lexer* p, fdm_fnpf *c, ghostcell *pgc, solver *psolv_reg, fnpf_fsf *pf, double *f)
{
    double sigxyz2;
    double ab,denom;
    double fbxm,fbxp,fbym,fbyp;
    p->poissoniter=0;
    p->poissontime=0.0;

	n=0;
    KJILOOP
	{
        if(p->wet[IJ]==1 && p->flag7[FIJK]>0)
        {
        sigxyz2 = pow(p->sigx[FIJK],2.0) + pow(p->sigy[FIJK],2.0) + pow(p->sigz[IJ],2.0);
        
        M[n*15]  =  1.0/(p->DXP[IP]*p->DXN[IP])*p->x_dir 
                    + 1.0/(p->DXP[IM1]*p->DXN[IP])*p->x_dir 
                    
                    + 1.0/(p->DYP[JP]*p->DYN[JP])*p->y_dir 
                    + 1.0/(p->DYP[JM1]*p->DYN[JP])*p->y_dir 
                    
                    + (sigxyz2/(p->DZP[KM1]*p->DZN[KP]))*p->z_dir
                    + (sigxyz2/(p->DZP[KM1]*p->DZN[KM1]))*p->z_dir;

        
        M[n*15+1] = -1.0/(p->DXP[IM1]*p->DXN[IP])*p->x_dir;
        M[n*15+2] = -1.0/(p->DXP[IP]*p->DXN[IP])*p->x_dir;
        
        M[n*15+3] = -1.0/(p->DYP[JM1]*p->DYN[JP])*p->y_dir;
        M[n*15+4] = -1.0/(p->DYP[JP]*p->DYN[JP])*p->y_dir;
        
        M[n*15+5] = -(sigxyz2/(p->DZP[KM1]*p->DZN[KM1]) - p->sigxx[FIJK]/(p->DZN[KP]+p->DZN[KM1]))*p->z_dir;
        M[n*15+6] = -(sigxyz2/(p->DZP[KM1]*p->DZN[KP])  + p->sigxx[FIJK]/(p->DZN[KP]+p->DZN[KM1]))*p->z_dir;
        
        M[n*15+7]  = -2.0*p->sigx[FIJK]/((p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]))*p->x_dir;
        M[n*15+8]  =  2.0*p->sigx[FIJK]/((p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]))*p->x_dir;
        M[n*15+9]  =  2.0*p->sigx[FIJK]/((p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]))*p->x_dir;
        M[n*15+10] = -2.0*p->sigx[FIJK]/((p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]))*p->x_dir;
        
        M[n*15+11] = -2.0*p->sigy[FIJK]/((p->DYN[JP]+p->DYN[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir;
        M[n*15+12] =  2.0*p->sigy[FIJK]/((p->DYN[JP]+p->DYN[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir;
        M[n*15+13] =  2.0*p->sigy[FIJK]/((p->DYN[JP]+p->DYN[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir;
        M[n*15+14] = -2.0*p->sigy[FIJK]/((p->DYN[JP]+p->DYN[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir;
        
        x[n] = f[FIJK];
        
        rhs[n] = 0.0;
                        
        }
        
        if(p->wet[IJ]==0 || p->flag7[FIJK]<0)
        {
        M[n*15]  =  1.0;
        M[n*15+1] = 0.0;
        M[n*15+2] = 0.0;
        M[n*15+3] = 0.0;
        M[n*15+4] = 0.0;
        M[n*15+5] = 0.0;
        M[n*15+6] = 0.0;
        M[n*15+7] = 0.0;
        M[n*15+8] = 0.0;
        M[n*15+9] = 0.0;
        M[n*15+10] = 0.0;
        M[n*15+11] = 0.0;
        M[n*15+12] = 0.0;
        M[n*15+13] = 0.0;
        M[n*15+14] = 0.0;
        
        x[n] = 0.0;
        rhs[n] =  0.0;
        }   
        
	++n;
	}
    
    n=0;
	KJILOOP
	{
            if(p->wet[IJ]==1 && p->flag7[FIJK]>0)
            {
            // south
            if((p->flag7[FIm1JK]<0 || p->wet[Im1J]==0) && c->bc(i-1,j)==0)
            {
            M[n*15] += M[n*15+1];  
            M[n*15+1] = 0.0;        
            }
            
            if(p->flag7[FIm1JK]<0 && c->bc(i-1,j)==1)
            {
            rhs[n] += M[n*15+1]*c->Uin[FIm1JK]*p->DXP[IM1];
            M[n*15] += M[n*15+1];
            M[n*15+1] = 0.0;
            }
            
            // north
            if((p->flag7[FIp1JK]<0 || p->wet[Ip1J]==0) && c->bc(i+1,j)==0)
            {
            M[n*15] += M[n*15+2];
            M[n*15+2] = 0.0;
            }
            
            if(p->flag7[FIp1JK]<0 && c->bc(i+1,j)==2)
            {
            rhs[n] -= M[n*15+2]*c->Uin[FIp1JK]*p->DXP[IP1];
            M[n*15] += M[n*15+2];
            M[n*15+2] = 0.0;
            }

            // east
            if(p->flag7[FIJm1K]<0 || p->wet[IJm1]==0)
            {
            M[n*15] += M[n*15+3];
            M[n*15+3] = 0.0;
            }
            
            // west
            if(p->flag7[FIJp1K]<0 || p->wet[IJp1]==0)
            {
            M[n*15] += M[n*15+4];
            M[n*15+4] = 0.0;
            }
            
            // top
            if(p->flag7[FIJKp2]<0 && p->flag7[FIJKp1]>0)
            {
            rhs[n] -= M[n*15+6]*f[FIJKp2];
            M[n*15+6] = 0.0;
            }
            
    
         
            // st
            if(p->flag7[FIm1JKp2]<0 && p->flag7[FIm1JKp1]>0)
            {
            rhs[n] -= M[n*15+8]*f[FIm1JKp2];
            M[n*15+8] = 0.0;
            }
            
            // nt
            if(p->flag7[FIp1JKp2]<0 && p->flag7[FIp1JKp1]>0)
            {
            rhs[n] -= M[n*15+10]*f[FIp1JKp2];
            M[n*15+10] = 0.0;
            }
            
            // sb KBEDBC
            if(p->flag7[FIm1JKm1]<0)
            {
            ab = -2.0*p->sigx[FIJK]/((p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]))*p->x_dir;
            
           denom = p->sigz[Im1J] + c->Bx(i-1,j)*p->sigx[FIm1JK] + c->By(i-1,j)*p->sigy[FIm1JK];

                    if(p->wet[Ip1J]==1 && p->wet[Im1J]==1)
                    {
                    M[n*15+2] +=  ab*2.0*p->DZN[KP]*c->Bx(i,j)/(denom*(p->DXP[IP] + p->DXP[IM1]));
                    M[n*15+1] += -ab*2.0*p->DZN[KP]*c->Bx(i,j)/(denom*(p->DXP[IP] + p->DXP[IM1]));
                    }
                    
                    if(p->wet[IJm1]==1 && p->wet[IJp1]==1)
                    {
                    M[n*15+4] +=  ab*2.0*p->DZN[KP]*c->By(i,j)/(denom*(p->DYP[JP] + p->DYP[JM1]));
                    M[n*15+3] += -ab*2.0*p->DZN[KP]*c->By(i,j)/(denom*(p->DYP[JP] + p->DYP[JM1]));
                    }
                
                M[n*15+8] += ab;
                M[n*15+7] = 0.0;

            }
            
            // nb KBEDBC
            if(p->flag7[FIp1JKm1]<0)
            {
            ab = 2.0*p->sigx[FIJK]/((p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]))*p->x_dir;
            
            denom = p->sigz[Ip1J] + c->Bx(i+1,j)*p->sigx[FIp1JK] + c->By(i+1,j)*p->sigy[FIp1JK];

                    if(p->wet[Ip1J]==1 && p->wet[Im1J]==1)
                    {
                    M[n*15+2] +=  ab*2.0*p->DZN[KP]*c->Bx(i,j)/(denom*(p->DXP[IP] + p->DXP[IM1]));
                    M[n*15+1] += -ab*2.0*p->DZN[KP]*c->Bx(i,j)/(denom*(p->DXP[IP] + p->DXP[IM1]));
                    }
                    
                    if(p->wet[IJm1]==1 && p->wet[IJp1]==1)
                    {
                    M[n*15+4] +=  ab*2.0*p->DZN[KP]*c->By(i,j)/(denom*(p->DYP[JP] + p->DYP[JM1]));
                    M[n*15+3] += -ab*2.0*p->DZN[KP]*c->By(i,j)/(denom*(p->DYP[JP] + p->DYP[JM1]));
                    }
                
                M[n*15+10] += ab;
                M[n*15+9] = 0.0;

            }
            
            
  /*
 * {{0,0}, {-1,0}, {1,0},  {0,-1}, {0,1}, {-1,-1},{-1,1},{1,-1},{1,1}};
p  0
s  1
n  2
b  3
t  4
sb 5
st 6
nb 7
nt 8

//-----------------

p  0
s  1
n  2
e  3
w  4
b  5
t  6
sb 7
st 8
nb 9
nt 10
eb 11
et 12
wb 13
wt 14
*/          
            
            // et
            if(p->flag7[FIJm1Kp2]<0 && p->flag7[FIJm1Kp1]>0)
            {
            rhs[n] -= M[n*15+12]*f[FIJp1Kp2];
            M[n*15+12] = 0.0;
            }
            
            // wt
            if(p->flag7[FIJp1Kp2]<0 && p->flag7[FIJp1Kp1]>0)
            {
            rhs[n] -= M[n*15+14]*f[FIJp1Kp2];
            M[n*15+14] = 0.0;
            }
            
            // eb KBEDBC
            if(p->flag7[FIJm1Km1]<0)
            {
            ab = -2.0*p->sigy[FIJK]/((p->DYN[JP]+p->DYN[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir;
            
            denom = p->sigz[IJm1] + c->Bx(i,j-1)*p->sigx[FIJm1K] + c->By(i,j-1)*p->sigy[FIJm1K];

                    if(p->wet[Ip1J]==1 && p->wet[Im1J]==1)
                    {
                    M[n*15+2] +=  ab*2.0*p->DZN[KP]*c->Bx(i,j)/(denom*(p->DXP[IP] + p->DXP[IM1]));
                    M[n*15+1] += -ab*2.0*p->DZN[KP]*c->Bx(i,j)/(denom*(p->DXP[IP] + p->DXP[IM1]));
                    }
                    
                    if(p->wet[IJm1]==1 && p->wet[IJp1]==1)
                    {
                    M[n*15+4] +=  ab*2.0*p->DZN[KP]*c->By(i,j)/(denom*(p->DYP[JP] + p->DYP[JM1]));
                    M[n*15+3] += -ab*2.0*p->DZN[KP]*c->By(i,j)/(denom*(p->DYP[JP] + p->DYP[JM1]));
                    }
                
                M[n*15+12] += ab;
                M[n*15+11] = 0.0;

            }
            
            // wb KBEDBC
            if(p->flag7[FIJp1Km1]<0)
            {
            ab = 2.0*p->sigy[FIJK]/((p->DYN[JP]+p->DYN[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir;
            
            denom = p->sigz[IJp1] + c->Bx(i,j+1)*p->sigx[FIJp1K] + c->By(i,j+1)*p->sigy[FIJp1K];

                    if(p->wet[Ip1J]==1 && p->wet[Im1J]==1)
                    {
                    M[n*15+2] +=  ab*2.0*p->DZN[KP]*c->Bx(i,j)/(denom*(p->DXP[IP] + p->DXP[IM1]));
                    M[n*15+1] += -ab*2.0*p->DZN[KP]*c->Bx(i,j)/(denom*(p->DXP[IP] + p->DXP[IM1]));
                    }
                    
                    if(p->wet[IJm1]==1 && p->wet[IJp1]==1)
                    {
                    M[n*15+4] +=  ab*2.0*p->DZN[KP]*c->By(i,j)/(denom*(p->DYP[JP] + p->DYP[JM1]));
                    M[n*15+3] += -ab*2.0*p->DZN[KP]*c->By(i,j)/(denom*(p->DYP[JP] + p->DYP[JM1]));
                    }
                
                M[n*15+14] += ab;
                M[n*15+13] = 0.0;

            }
            
            
            
            
            
 
            // KBEDBC
            if(p->flag7[FIJKm1]<0)
            {
            sigxyz2= pow(p->sigx[FIJK],2.0) + pow(p->sigy[FIJK],2.0) + pow(p->sigz[IJ],2.0);
            
            ab = -(sigxyz2/(p->DZP[KM1]*p->DZN[KM1]) - p->sigxx[FIJK]/(p->DZN[KP]+p->DZN[KM1]));
            
            denom = p->sigz[IJ] + c->Bx(i,j)*p->sigx[FIJK] + c->By(i,j)*p->sigy[FIJK];

                    if(p->wet[Ip1J]==1 && p->wet[Im1J]==1)
                    {
                    M[n*15+2] +=  ab*2.0*p->DZN[KP]*c->Bx(i,j)/(denom*(p->DXP[IP] + p->DXP[IM1]));
                    M[n*15+1] += -ab*2.0*p->DZN[KP]*c->Bx(i,j)/(denom*(p->DXP[IP] + p->DXP[IM1]));
                    }
                    
                    if(p->wet[IJm1]==1 && p->wet[IJp1]==1)
                    {
                    M[n*15+4] +=  ab*2.0*p->DZN[KP]*c->By(i,j)/(denom*(p->DYP[JP] + p->DYP[JM1]));
                    M[n*15+3] += -ab*2.0*p->DZN[KP]*c->By(i,j)/(denom*(p->DYP[JP] + p->DYP[JM1]));
                    }
                
                
                M[n*15+6] += ab;
                M[n*15+5] = 0.0;

            }
            }
	++n;
	}

    double starttime=pgc->timer();
    psolv->startF(p,pgc,x,rhs,M,8);
    double endtime=pgc->timer();
    
        n=0;
        KJILOOP
        {
		 FPWDCHECK
        f[FIJK]=x[n];
		
        ++n;
        }
    
    p->poissoniter+=p->solveriter;
    p->poissontime+=endtime-starttime;
    
    
	if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"Fi_iter: "<<p->poissoniter<<" Final_residual: "<<p->final_res<<"  Fi_time: "<<setprecision(3)<<p->poissontime<<endl;
}
