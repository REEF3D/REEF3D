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
Author: Fabian Knoblauch
--------------------------------------------------------------------*/

#include"VOF_PLIC.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"


#include"convection.h"
#include"solver.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"fluid_update_vof.h"
#include"heat.h"
#include"hires.h"
#include"weno_hj.h"
#include"hric.h"

void VOF_PLIC::reconstructPlane_alt(fdm* a, lexer* p, field& voffield)
{
    double V0, n_1, n_2, n_3, V,r0, r, n_a, n_b, n_c;
    V0=voffield(i,j,k);
    //get the normal vector
    switch(p->F88)
    {
        case 0:
                break;
                
        case 1:
                break;
        case 2: 
                calcNormalWeymouth(a,p,voffield);
                break;
        case 3:
                calcNormalWang(a,p);
                break;
        case 4:
                calcNormalFO(a,p,voffield);
                break;
        case 5:
                calcNormalLS(a,p,voffield);
                break;
        case 6:
                calcNormalWENO(a,p,voffield);
                break;
        case 7:
                calcNormalPhi(a,p);
                break;
        case 8:
                calcNormalMassCentre(a,p,voffield);
                break;
        case 9:
                calcNormalELVIRA2D(a,p,voffield);
                break;
        case 10:
                calcNormalMYC2D(a,p,voffield);
                break;
        case 11:
                calcNormalMYC2D_V2(a,p,voffield);
                break;
        case 12:
                calcNormalMYC2D_V3(a,p,voffield);
                break;
        case 13:
                calcNormalMYC2D_V4(a,p,voffield);
                break;
        
    }
    //normalise normal vector (to be sure)
    
    if(p->F88 != 9)
    {
    double vecsum = sqrt(nx(i,j,k)*nx(i,j,k)+ny(i,j,k)*ny(i,j,k)+nz(i,j,k)*nz(i,j,k));
    nx(i,j,k)=nx(i,j,k)/vecsum;
    ny(i,j,k)=ny(i,j,k)/vecsum;
    nz(i,j,k)=nz(i,j,k)/vecsum;
  //  cout<<"nx:"<<nx(i,j,k)<<" ny:"<<ny(i,j,k)<<" nz:"<<nz(i,j,k)<<" i:"<<i<<" j:"<<j<<" k:"<<k<<endl;
    
    //scale with cellsize
    n_a=fabs(nx(i,j,k))*p->DXN[IP];
    n_b=fabs(ny(i,j,k))*p->DYN[JP];
    n_c=fabs(nz(i,j,k))*p->DZN[KP];
    
    //normalise scaled vector
    double vecsum2 = sqrt(n_a*n_a+n_b*n_b+n_c*n_c);
    n_a=n_a/vecsum2;
    n_b=n_b/vecsum2;
    n_c=n_c/vecsum2;
    
    //sort vector components
    if(n_b>=n_a)
    {
        n_1=n_a;
        n_2=n_b;
    }
    else
    {
        n_1=n_b;
        n_2=n_a;
    }
    
    if(n_c>=n_2)
    {
        n_3=n_c;
    }
    else if(n_c<n_1)
    {
        n_3=n_2;
        n_2=n_1;
        n_1=n_c;
    }
    else
    {
        n_3=n_2;
        n_2=n_c;
    }
    
    //reduced symmetry transformation
    V=0.5-fabs(V0-0.5);
    if(n_1+n_2<=2.0*V*n_3)   //case 5
    {
        r=V*n_3+0.5*(n_1+n_2);
    }
    else if((3.0*n_2*(2.0*V*n_3+n_1-n_2)<=n_1*n_1) && (n_1*n_1<=6.0*V*n_2*n_3)) //case 2
    {
        r=0.5*n_1+sqrt(2*V*n_2*n_3-1.0/12.0*n_1*n_1);
    }
    
    else if(6.0*V*n_2*n_3<n_1*n_1) //case 1
    {
        r=cbrt(6.0*V*n_1*n_2*n_3);
    }
    
    else //cases 3 und 4
    {   
        double x3, y32, u3,f3;
        x3=81.0*n_1*n_2*(n_1+n_2-2.0*n_3*V);
        y32=fdim(23328.0*n_1*n_1*n_1*n_2*n_2*n_2,x3*x3);
        u3=cbrt(x3*x3+y32);
        f3=n_1+n_2-(7.5595264*n_1*n_2+0.26456684*u3)*1.0/sqrt(u3)*sin(0.5235988-1.0/3.0*atan(sqrt(y32)/x3));
        if(f3<=n_3) //case 3
        {
            r=f3;
        }
        else //case 4
        {   
            double t4,x4,y42,u4,f4;
            t4=9.0*(n_1+n_2+n_3)*(n_1+n_2+n_3)-18.0;
            x4=max(n_1*n_2*n_3*(324.0-648.0*V),1E-35);
            y42=fdim(4.0*t4*t4*t4,x4*x4);
            u4=cbrt(x4*x4+y42);
            f4=0.5*(n_1+n_2+n_3)-(0.20998684*t4+0.13228342*u4)*1.0/sqrt(u4)*sin(0.5235988-1.0/3.0*atan(sqrt(y42)/x4));
            
            r=f4;
        }
    }
    //reverse reduced symmetry
    if(V0-0.5>=0.0)
    {
        r0=0.5*(n_1+n_2+n_3)-r;
    }
    else
    {
        r0=-(0.5*(n_1+n_2+n_3)-r);
    }
    
    alpha(i,j,k)=r0*sqrt(nx(i,j,k)*nx(i,j,k)*p->DXN[IP]*p->DXN[IP]+ny(i,j,k)*ny(i,j,k)*p->DYN[JP]*p->DYN[JP]+nz(i,j,k)*nz(i,j,k)*p->DZN[KP]*p->DZN[KP]);
   // cout<<"alpha: "<<alpha(i,j,k)<<" nx:"<<nx(i,j,k)<<" ny:"<<ny(i,j,k)<<" nz:"<<nz(i,j,k)<<endl;
   }
}
    
double VOF_PLIC::calcAlphaFromInput(fdm* a, lexer* p, double n_x, double n_y, double n_z, double d_x, double d_y, double d_z, double V0)
{
    double vecsum = sqrt(n_x*n_x+n_y*n_y+n_z*n_z);
    double n_a,n_b,n_c,n_1,n_2,n_3,V,r0,ret,r;
    n_x=n_x/vecsum;
    n_y=n_y/vecsum;
    n_z=n_z/vecsum;
  //  cout<<"nx:"<<nx(i,j,k)<<" ny:"<<ny(i,j,k)<<" nz:"<<nz(i,j,k)<<" i:"<<i<<" j:"<<j<<" k:"<<k<<endl;
    
    //scale with cellsize
    n_a=fabs(n_x)*d_x;
    n_b=fabs(n_y)*d_y;
    n_c=fabs(n_z)*d_z;
    
    //normalise scaled vector
    double vecsum2 = sqrt(n_a*n_a+n_b*n_b+n_c*n_c);
    n_a=n_a/vecsum2;
    n_b=n_b/vecsum2;
    n_c=n_c/vecsum2;
    
    //sort vector components
    if(n_b>=n_a)
    {
        n_1=n_a;
        n_2=n_b;
    }
    else
    {
        n_1=n_b;
        n_2=n_a;
    }
    
    if(n_c>=n_2)
    {
        n_3=n_c;
    }
    else if(n_c<n_1)
    {
        n_3=n_2;
        n_2=n_1;
        n_1=n_c;
    }
    else
    {
        n_3=n_2;
        n_2=n_c;
    }
    
    //reduced symmetry transformation
    V=0.5-fabs(V0-0.5);
    if(n_1+n_2<=2.0*V*n_3)   //case 5
    {
        r=V*n_3+0.5*(n_1+n_2);
    }
    else if((3.0*n_2*(2.0*V*n_3+n_1-n_2)<=n_1*n_1) && (n_1*n_1<=6.0*V*n_2*n_3)) //case 2
    {
        r=0.5*n_1+sqrt(2*V*n_2*n_3-1.0/12.0*n_1*n_1);
    }
    
    else if(6.0*V*n_2*n_3<n_1*n_1) //case 1
    {
        r=cbrt(6.0*V*n_1*n_2*n_3);
    }
    
    else //cases 3 und 4
    {   
        double x3, y32, u3,f3;
        x3=81.0*n_1*n_2*(n_1+n_2-2.0*n_3*V);
        y32=fdim(23328.0*n_1*n_1*n_1*n_2*n_2*n_2,x3*x3);
        u3=cbrt(x3*x3+y32);
        f3=n_1+n_2-(7.5595264*n_1*n_2+0.26456684*u3)*1.0/sqrt(u3)*sin(0.5235988-1.0/3.0*atan(sqrt(y32)/x3));
        if(f3<=n_3) //case 3
        {
            r=f3;
        }
        else //case 4
        {   
            double t4,x4,y42,u4,f4;
            t4=9.0*(n_1+n_2+n_3)*(n_1+n_2+n_3)-18.0;
            x4=max(n_1*n_2*n_3*(324.0-648.0*V),1E-35);
            y42=fdim(4.0*t4*t4*t4,x4*x4);
            u4=cbrt(x4*x4+y42);
            f4=0.5*(n_1+n_2+n_3)-(0.20998684*t4+0.13228342*u4)*1.0/sqrt(u4)*sin(0.5235988-1.0/3.0*atan(sqrt(y42)/x4));
            
            r=f4;
        }
    }
    //reverse reduced symmetry
    if(V0-0.5>=0.0)
    {
        r0=0.5*(n_1+n_2+n_3)-r;
    }
    else
    {
        r0=-(0.5*(n_1+n_2+n_3)-r);
    }
    
    ret=r0*sqrt(n_x*n_x*d_x*d_x+n_y*n_y*d_y*d_y+n_z*n_z*d_z*d_z);
    
    return ret;
}

double VOF_PLIC::return_alpha_reconstructPlane_alt(fdm* a, lexer* p, field& voffield, int i ,int j , int k)
{
    double V0, n_1, n_2, n_3, V,r0, r, n_a, n_b, n_c;
    V0=voffield(i,j,k);
    //get the normal vector
    switch(p->F88)
    {
        case 0:
                break;
                
        case 1:
                break;
        case 2: 
                calcNormalWeymouth(a,p,voffield);
                break;
        case 3:
                calcNormalWang(a,p);
                break;
        case 4:
                calcNormalFO(a,p,voffield);
                break;
        case 5:
                calcNormalLS(a,p,voffield);
                break;
        case 6:
                calcNormalWENO(a,p,voffield);
                break;
        case 7:
                calcNormalPhi(a,p);
                break;
        case 8:
                calcNormalMassCentre(a,p,voffield);
                break;
        case 9:
                calcNormalELVIRA2D(a,p,voffield);
                break;
        case 10:
                calcNormalMYC2D(a,p,voffield);
                break;
    }
    //normalise normal vector (to be sure)
    
    if(p->F88 != 9)
    {
    double vecsum = sqrt(nx(i,j,k)*nx(i,j,k)+ny(i,j,k)*ny(i,j,k)+nz(i,j,k)*nz(i,j,k));
    nx(i,j,k)=nx(i,j,k)/vecsum;
    ny(i,j,k)=ny(i,j,k)/vecsum;
    nz(i,j,k)=nz(i,j,k)/vecsum;
  //  cout<<"nx:"<<nx(i,j,k)<<" ny:"<<ny(i,j,k)<<" nz:"<<nz(i,j,k)<<" i:"<<i<<" j:"<<j<<" k:"<<k<<endl;
    
    //scale with cellsize
    n_a=fabs(nx(i,j,k))*p->DXN[IP];
    n_b=fabs(ny(i,j,k))*p->DYN[JP];
    n_c=fabs(nz(i,j,k))*p->DZN[KP];
    
    //normalise scaled vector
    double vecsum2 = sqrt(n_a*n_a+n_b*n_b+n_c*n_c);
    n_a=n_a/vecsum2;
    n_b=n_b/vecsum2;
    n_c=n_c/vecsum2;
    
    //sort vector components
    if(n_b>=n_a)
    {
        n_1=n_a;
        n_2=n_b;
    }
    else
    {
        n_1=n_b;
        n_2=n_a;
    }
    
    if(n_c>=n_2)
    {
        n_3=n_c;
    }
    else if(n_c<n_1)
    {
        n_3=n_2;
        n_2=n_1;
        n_1=n_c;
    }
    else
    {
        n_3=n_2;
        n_2=n_c;
    }
    
    //reduced symmetry transformation
    V=0.5-fabs(V0-0.5);
    if(n_1+n_2<=2.0*V*n_3)   //case 5
    {
        r=V*n_3+0.5*(n_1+n_2);
    }
    else if((3.0*n_2*(2.0*V*n_3+n_1-n_2)<=n_1*n_1) && (n_1*n_1<=6.0*V*n_2*n_3)) //case 2
    {
        r=0.5*n_1+sqrt(2*V*n_2*n_3-1.0/12.0*n_1*n_1);
    }
    
    else if(6.0*V*n_2*n_3<n_1*n_1) //case 1
    {
        r=cbrt(6.0*V*n_1*n_2*n_3);
    }
    
    else //cases 3 und 4
    {   
        double x3, y32, u3,f3;
        x3=81.0*n_1*n_2*(n_1+n_2-2.0*n_3*V);
        y32=fdim(23328.0*n_1*n_1*n_1*n_2*n_2*n_2,x3*x3);
        u3=cbrt(x3*x3+y32);
        f3=n_1+n_2-(7.5595264*n_1*n_2+0.26456684*u3)*1.0/sqrt(u3)*sin(0.5235988-1.0/3.0*atan(sqrt(y32)/x3));
        if(f3<=n_3) //case 3
        {
            r=f3;
        }
        else //case 4
        {   
            double t4,x4,y42,u4,f4;
            t4=9.0*(n_1+n_2+n_3)*(n_1+n_2+n_3)-18.0;
            x4=max(n_1*n_2*n_3*(324.0-648.0*V),1E-35);
            y42=fdim(4.0*t4*t4*t4,x4*x4);
            u4=cbrt(x4*x4+y42);
            f4=0.5*(n_1+n_2+n_3)-(0.20998684*t4+0.13228342*u4)*1.0/sqrt(u4)*sin(0.5235988-1.0/3.0*atan(sqrt(y42)/x4));
            
            r=f4;
        }
    }
    //reverse reduced symmetry
    if(V0-0.5>=0.0)
    {
        r0=0.5*(n_1+n_2+n_3)-r;
    }
    else
    {
        r0=-(0.5*(n_1+n_2+n_3)-r);
    }
    
    alpha(i,j,k)=r0*sqrt(nx(i,j,k)*nx(i,j,k)*p->DXN[IP]*p->DXN[IP]+ny(i,j,k)*ny(i,j,k)*p->DYN[JP]*p->DYN[JP]+nz(i,j,k)*nz(i,j,k)*p->DZN[KP]*p->DZN[KP]);
   // cout<<"alpha: "<<alpha(i,j,k)<<" nx:"<<nx(i,j,k)<<" ny:"<<ny(i,j,k)<<" nz:"<<nz(i,j,k)<<endl;
   }
   
   return alpha(i,j,k);
}