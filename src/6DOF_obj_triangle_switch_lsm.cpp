/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"6DOF_obj.h"
#include"lexer.h"
#include"ghostcell.h"

void sixdof_obj::triangle_switch_lsm(lexer *p, ghostcell *pgc)
{
    /*
    double x0,x1,x2,y0,y1,y2,z0,z1,z2;
	double xc,yc,zc;
	double at,bt,ct,st;
	double nx,ny,nz,norm;
    double n0,n1,n2;
    double fbval;
    
    // switch list allocate
    p->Iarray(tri_switch,tricount);
    p->Iarray(tri_switch_id,tricount);
    p->Iarray(tricount_local_list,p->M10+1);
    p->Iarray(tricount_local_displ,p->M10+1);
    
    for (int n=0;n<tricount;++n)
    {
    tri_switch[n]=0;
    tri_switch_id[n]=0;
    }
    
    tricount_local=0;
    
    for(n=0;n<p->M10;++n)
    {
    tricount_local_list[n]=0;
    tricount_local_displ[n]=0;
    }
    
    
    // start loop part 1: find local triangles
    for (int n=0;n<tricount;++n)
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
        
        // switch list
        if( xc >= p->originx && xc <= p->endx &&
            ((yc >= p->originy && yc <= p->endy) || p->j_dir==0) &&
            zc >= p->originz && zc <= p->endz)
            {
            ++tricount_local;
            }
    }
            
    
    // Exchange 1: switch list count
    pgc->allgather_int(&tricount_local,1,tricount_local_list,1);
    

        
    p->Iarray(tri_switch_local,tricount_local);
    p->Iarray(tri_switch_local_id,tricount_local);
    
    tricount_local_displ[0]=0;
    
    for(q=1;q<p->M10;++q)
    {
    tricount_local_displ[q] = tricount_local_displ[q-1] + tricount_local_list[q-1];
    
    //cout<<p->mpirank<<" tricount_local: "<<tricount_local<<" tricount_local_list[q]: "<<tricount_local_list[q]<<" tricount_local_list[q]: "<<tricount_local_displ[q]<<endl;
    }

        
    // start loop part 2: determine triangle switch
    count=0;
    for (int n=0;n<tricount;++n)
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
        
        //if((p->j_dir==1 || fabs(ny)<0.001))
        //cout<<"yc: "<<yc<<" j: "<<j<<" fbval: "<<fbval<<endl;
   
        if( xc >= p->originx && xc <= p->endx &&
            ((yc >= p->originy && yc <= p->endy) || p->j_dir==0) &&
            zc >= p->originz && zc <= p->endz)
        {
            i = p->posc_i(xc);
            j = p->posc_j(yc);
            k = p->posc_k(zc);
            
            xc += 0.1*nx*p->DXN[IP];
            yc += 0.1*ny*p->DYN[JP];
            zc += 0.1*nz*p->DZN[KP];
            
            fbval = p->ccipol4_a(a->fb,xc,yc,zc);
            
            // Normal vector sign
            n0 = (a->fb(i+1,j,k) - a->fb(i-1,j,k))/(2.0*p->DXN[IP]);
            n1 = (a->fb(i,j+1,k) - a->fb(i,j-1,k))/(2.0*p->DYN[JP]);
            n2 = (a->fb(i,j,k+1) - a->fb(i,j,k-1))/(2.0*p->DZN[KP]);

            norm = sqrt(n0*n0 + n1*n1 + n2*n2);
            
             n0 /= norm > 1.0e-20 ? norm : 1.0e20;
			n1 /= norm > 1.0e-20 ? norm : 1.0e20;
			n2 /= norm > 1.0e-20 ? norm : 1.0e20;

            if((n0>0.0 && nx<0.0 || n0<0.0 && nx>0.0) || (n2>0.0 && nz<0.0 || n2<0.0 && nz>0.0))
            fbval=-1.0;
        
        
            if(fbval>=0.0)
            {
            tri_switch_local[count] = 0;
            
            }
            
            if(fbval<0.0)
            {
            tri_switch_local[count] = 1;

            }
            
            tri_switch_local_id[count]=n;

            ++count;
        }
    }
        
        // Exchange 2:  switch list
        pgc->allgatherv_int(tri_switch_local_id, tricount_local_list[p->mpirank], tri_switch_id, tricount_local_list, tricount_local_displ);
        
        pgc->allgatherv_int(tri_switch_local, tricount_local_list[p->mpirank], tri_switch, tricount_local_list, tricount_local_displ);
        
        
        // start loop part 3: global switch
        tricount_switch_total=0;
        for (int n=0;n<tricount;++n)
        if(tri_switch[tri_switch_id[n]]==1)
        {   
        cout<<n<<" tri_switch_id[n]: "<<tri_switch_id[n]<<endl;
        x1 = tri_x[tri_switch_id[n]][1];
        y1 = tri_y[tri_switch_id[n]][1];
        z1 = tri_z[tri_switch_id[n]][1];
        
        x2 = tri_x[tri_switch_id[n]][2];
        y2 = tri_y[tri_switch_id[n]][2];
        z2 = tri_z[tri_switch_id[n]][2];
        
        
        tri_x[tri_switch_id[n]][1] = x2;
        tri_y[tri_switch_id[n]][1] = y2;
        tri_z[tri_switch_id[n]][1] = z2;
        
        tri_x[tri_switch_id[n]][2] = x1;
        tri_y[tri_switch_id[n]][2] = y1;
        tri_z[tri_switch_id[n]][2] = z1;
        
        ++tricount_switch_total;
        }
    
    ++triangle_token;
    
    if(p->mpirank==0)
    cout<<"6DOF STL triangle switch count: "<<tricount_switch_total<<endl;
    
    // free allocated arrays
    p->del_Iarray(tri_switch,tricount);
    p->del_Iarray(tri_switch_id,tricount);
    p->del_Iarray(tricount_local_list,p->M10+1);
    p->del_Iarray(tricount_local_displ,p->M10+1);
    p->del_Iarray(tri_switch_local,tricount_local);
    p->del_Iarray(tri_switch_local_id,tricount_local);*/
}