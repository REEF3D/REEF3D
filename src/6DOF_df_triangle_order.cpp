/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_df_object.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_df_object::triangle_order(lexer *p, fdm *a, ghostcell *pgc)
{
    double x0,x1,x2,y0,y1,y2,z0,z1,z2;
	double xc,yc,zc;
	double at,bt,ct,st;
	double nx,ny,nz,norm;
    double n0,n1,n2;
    double fbval;
    
    for (int n = 0; n < tricount; ++n)
    {
        x0 = tri_x[n][0];
        y0 = tri_y[n][0];
        z0 = tri_z[n][0];
        
        x1 = tri_x[n][1];
        y1 = tri_y[n][1];
        z1 = tri_z[n][1];
        
        x2 = tri_x[n][2];
        y2 = tri_y[n][2];
        z2 = tri_z[n][2]; 
        
        nx = (y1 - y0)*(z2 - z0) - (y2 - y0)*(z1 - z0);
        ny = (x2 - x0)*(z1 - z0) - (x1 - x0)*(z2 - z0); 
        nz = (x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0);

        norm = sqrt(nx*nx + ny*ny + nz*nz);
			
        nx /= norm > 1.0e-20 ? norm : 1.0e20;
        ny /= norm > 1.0e-20 ? norm : 1.0e20;
        nz /= norm > 1.0e-20 ? norm : 1.0e20;
        
        
        // Center of triangle
		xc = (x0 + x1 + x2)/3.0;
		yc = (y0 + y1 + y2)/3.0;
		zc = (z0 + z1 + z2)/3.0;
        
        i = p->posc_i(xc);
        j = p->posc_j(yc);
        k = p->posc_k(zc);
        
        xc += nx*p->DXN[IP];
        yc += ny*p->DYN[JP];
        zc += nz*p->DZN[KP];
        
        fbval = p->ccipol4_a(a->fb,xc,yc,zc);
        
        
        // Normal vector sign
            n0 = (a->fb(i+1,j,k) - a->fb(i-1,j,k))/(2.0*p->DXN[IP]);
            n1 = (a->fb(i,j+1,k) - a->fb(i,j-1,k))/(2.0*p->DYN[JP]);
            n2 = (a->fb(i,j,k+1) - a->fb(i,j,k-1))/(2.0*p->DZN[KP]);
    
            norm = sqrt(n0*n0 + n1*n1 + n2*n2);
            
            n0 /= norm>1.0e-20?norm:1.0e20;
            n1 /= norm>1.0e-20?norm:1.0e20;
            n2 /= norm>1.0e-20?norm:1.0e20;
            
            
        
        if(fbval<0.0 && (xc >= p->originx && xc < p->endx &&
                         yc >= p->originy && yc < p->endy &&
                         zc >= p->originz && zc < p->endz))
        {
        cout<<"TRIANGLE SWITCH "<<fbval<<" | nx: "<<nx<<" ny: "<<ny<<" nz: "<<nz<<" || n0: "<<n0<<" n1: "<<n1<<" n2: "<<n2;
        
        cout<<" | xc: "<<xc<<" yc: "<<yc<<" zc: "<<zc<<endl;
        
        if(triangle_token==1)
        cout<<"TRIANGLE SWITCH !!!!!!!!"<<endl;
        
        
        // --- test
        tri_x[n][1] = x2;
        tri_y[n][1] = y2;
        tri_z[n][1] = z2;
        
        tri_x[n][2] = x1;
        tri_y[n][2] = y1;
        tri_z[n][2] = z1;
        
        
        x0 = tri_x[n][0];
        y0 = tri_y[n][0];
        z0 = tri_z[n][0];
        
        x1 = tri_x[n][1];
        y1 = tri_y[n][1];
        z1 = tri_z[n][1];
        
        x2 = tri_x[n][2];
        y2 = tri_y[n][2];
        z2 = tri_z[n][2]; 
        
        nx = (y1 - y0)*(z2 - z0) - (y2 - y0)*(z1 - z0);
        ny = (x2 - x0)*(z1 - z0) - (x1 - x0)*(z2 - z0); 
        nz = (x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0);

        norm = sqrt(nx*nx + ny*ny + nz*nz);
			
        nx /= norm > 1.0e-20 ? norm : 1.0e20;
        ny /= norm > 1.0e-20 ? norm : 1.0e20;
        nz /= norm > 1.0e-20 ? norm : 1.0e20;
        
        // Center of triangle
		xc = (x0 + x1 + x2)/3.0;
		yc = (y0 + y1 + y2)/3.0;
		zc = (z0 + z1 + z2)/3.0;
        
        i = p->posc_i(xc);
        j = p->posc_j(yc);
        k = p->posc_k(zc);
        
        xc += nx*p->DXN[IP];
        yc += ny*p->DYN[JP];
        zc += nz*p->DZN[KP];
        
        fbval = p->ccipol4_a(a->fb,xc,yc,zc);
        
        //cout<<" ||| nx: "<<nx<<" ny: "<<ny<<" nz: "<<nz<<" | fbval: "<<fbval<<endl;
        }
    }
    ++triangle_token;
    
    
}