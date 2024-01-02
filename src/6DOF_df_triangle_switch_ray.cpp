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

#include"6DOF_df_object.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"fieldint.h"

void sixdof_df_object::triangle_switch_ray(lexer *p, fdm *a, ghostcell *pgc)
{
	double ys,ye,zs,ze;
	double Px,Py,Pz;
	double Qx,Qy,Qz;
    double Sx,Sy,Sz;
    double Tx,Ty,Tz;
	double Rx,Ry,Rz;
	double Ax,Ay,Az;
	double Bx,By,Bz;
	double Cx,Cy,Cz;
	double PQx,PQy,PQz;
	double PAx,PAy,PAz;
	double PBx,PBy,PBz;
	double PCx,PCy,PCz;
	double Mx,My,Mz;
	double u,v,w;
	double denom;	
	double psi = 1.0e-8*p->DXM;
    
    double x0,x1,x2,y0,y1,y2,z0,z1,z2;
	double xc,yc,zc;
	double at,bt,ct,st;
	double nx,ny,nz,norm;
    int tricount_local_max,sum;
    int cutnum;
    int domdir;
    
    // allocate
    p->Iarray(tricount_local_list,p->M10+1);
    p->Iarray(tricount_local_displ,p->M10+1);
    
    p->Iarray(tri_switch,tricount);
    
    for (n=0;n<tricount;++n)
    tri_switch[n]=0;
    
    // divide triangles to local processors
    tricount_local = int(tricount/p->M10);
    
    sum = 0;
    for(q=0; q<p->M10-1; ++q)
    {
    tricount_local_list[q]=tricount_local;
    sum += tricount_local;
    }
    
    tricount_local_list[p->M10-1] = tricount - sum;

    
    // displacement
    tricount_local_displ[0]=0;
    
    for(q=1;q<p->M10+1;++q)
    tricount_local_displ[q] = tricount_local_displ[q-1] + tricount_local_list[q-1];
    
    tricount_local_max = 0;
    for(q=0; q<p->M10; ++q)
    tricount_local_max = MAX(tricount_local_list[q],tricount_local_max);
    
    p->Iarray(tri_switch_local,tricount_local_max);
    
    /*
    cout<<p->mpirank<<" tricount_local_displ[p->mpirank+1] - tricount_local_displ[p->mpirank] "<<tricount_local_displ[p->mpirank+1] - tricount_local_displ[p->mpirank]<<endl;
    
    if(p->mpirank==1)
    for(q=0;q<p->M10+1;++q)
    cout<<"tricount_local_displ[q]: "<<tricount_local_displ[q]<<endl;
    */

    // ray cast
	for(q=tricount_local_displ[p->mpirank];q<tricount_local_displ[p->mpirank+1];++q)
	{
        //cout<<p->mpirank<<"q: "<<q<<endl;
        
        // triangle points
        x0 = tri_x[q][0];
        y0 = tri_y[q][0];
        z0 = tri_z[q][0];
        
        x1 = tri_x[q][1];
        y1 = tri_y[q][1];
        z1 = tri_z[q][1];
        
        x2 = tri_x[q][2];
        y2 = tri_y[q][2];
        z2 = tri_z[q][2]; 
        
        // normals
        nx = (y1 - y0)*(z2 - z0) - (y2 - y0)*(z1 - z0);
        ny = (x2 - x0)*(z1 - z0) - (x1 - x0)*(z2 - z0); 
        nz = (x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0);

        norm = sqrt(nx*nx + ny*ny + nz*nz);
			
        nx /= norm>1.0e-20?norm:1.0e20;
        ny /= norm>1.0e-20?norm:1.0e20;
        nz /= norm>1.0e-20?norm:1.0e20;
        
        domdir = 0;
        
        if(fabs(nx)>=fabs(ny))
        {
        domdir=1;
        
            if(fabs(nx)<fabs(nz))
            domdir=3;
        }
        
        if(fabs(nx)<fabs(ny))
        {
        domdir=2;
        
            if(fabs(ny)<fabs(nz))
            domdir=3;
        }
        
        // Center of triangle
		xc = (x0 + x1 + x2)/3.0;
		yc = (y0 + y1 + y2)/3.0;
		zc = (z0 + z1 + z2)/3.0;
        
        
            cutnum=0;
            // ray cast loop
            for(n=0; n<tricount; ++n)
            if(n!=q)
            {
            Ax = tri_x[n][0];
            Ay = tri_y[n][0];
            Az = tri_z[n][0];
            
            Bx = tri_x[n][1];
            By = tri_y[n][1];
            Bz = tri_z[n][1];
            
            Cx = tri_x[n][2];
            Cy = tri_y[n][2];
            Cz = tri_z[n][2]; 
            
            
            xs = MIN3(Ax,Bx,Cx);
            xe = MAX3(Ax,Bx,Cx);
    
            ys = MIN3(Ay,By,Cy);
            ye = MAX3(Ay,By,Cy);
            
            zs = MIN3(Az,Bz,Cz);
            ze = MAX3(Az,Bz,Cz);
            
            
            if((yc>=ys-psi && yc<=ye+psi && zc>=zs-psi && zc<=ze+psi && domdir==1)
            || (xc>=xs-psi && xc<=xe+psi && zc>=zs-psi && zc<=ze+psi && domdir==2 && p->y_dir==1) 
            || (xc>=xs-psi && xc<=xe+psi && yc>=ys-psi && yc<=ye+psi && domdir==3)) 
            {
            //cout<<"TRIANGLE SWITCH "<<domdir<<" | nx: "<<nx<<" ny: "<<ny<<" nz: "<<nz<<endl;
        
    
                Px = xc;
                Py = yc;
                Pz = zc;
                
                Qx = xc;
                Qy = yc;
                Qz = zc;
                
                // domninate direction 1
                if(domdir==1 && nx>0.0)
                {
                Px = p->global_xmin-10.0*p->DXM;
                Qx = p->global_xmax+10.0*p->DXM;
                
                Sx = xc;
                Tx = p->global_xmax+10.0*p->DXM;
                }
                
                if(domdir==1 && nx<0.0)
                {
                Px = p->global_xmin-10.0*p->DXM;
                Qx = p->global_xmax+10.0*p->DXM;
                
                Sx = p->global_xmin-10.0*p->DXM;
                Tx = xc;
                }
                
                // domninate direction 2
                if(domdir==2 && ny>0.0)
                {
                Py = p->global_ymin-1000.0*p->DXM;
                Qy = p->global_ymax+1000.0*p->DXM;
                
                Sy = yc;
                Ty = p->global_ymax+1000.0*p->DXM;
                }
                
                if(domdir==2 && ny<0.0)
                {
                Py = p->global_ymin-1000.0*p->DXM;
                Qy = p->global_ymax+1000.0*p->DXM;
                
                Sy = p->global_ymin-1000.0*p->DXM;
                Ty = yc;
                }
                
                // domninate direction 3
                if(domdir==3 && nz>0.0)
                {
                Pz = p->global_zmin-1000.0*p->DXM;
                Qz = p->global_zmax+1000.0*p->DXM;
                
                Sz = zc;
                Tz = p->global_zmax+1000.0*p->DXM;
                }
                
                if(domdir==3 && nz<0.0)
                {
                Pz = p->global_zmin-1000.0*p->DXM;
                Qz = p->global_zmax+1000.0*p->DXM;
                
                Sz = p->global_zmin-1000.0*p->DXM;
                Tz = zc;
                }
                
                
                PQx = Qx-Px;
                PQy = Qy-Py;
                PQz = Qz-Pz;
                
                PAx = Ax-Px;
                PAy = Ay-Py;
                PAz = Az-Pz;
                
                PBx = Bx-Px;
                PBy = By-Py;
                PBz = Bz-Pz;
                
                PCx = Cx-Px;
                PCy = Cy-Py;
                PCz = Cz-Pz;
                
                // uvw
                Mx = PQy*Pz - PQz*Py;
                My = PQz*Px - PQx*Pz;
                Mz = PQx*Py - PQy*Px;

                
                u = PQx*(Cy*Bz - Cz*By) + PQy*(Cz*Bx - Cx*Bz) + PQz*(Cx*By - Cy*Bx)
                  + Mx*(Cx-Bx) + My*(Cy-By) + Mz*(Cz-Bz);
                  
                v = PQx*(Ay*Cz - Az*Cy) + PQy*(Az*Cx - Ax*Cz) + PQz*(Ax*Cy - Ay*Cx)
                  + Mx*(Ax-Cx) + My*(Ay-Cy) + Mz*(Az-Cz);
                  
                w = PQx*(By*Az - Bz*Ay) + PQy*(Bz*Ax - Bx*Az) + PQz*(Bx*Ay - By*Ax)
                  + Mx*(Bx-Ax) + My*(By-Ay) + Mz*(Bz-Az);
                
                
                int check=1;
                if(fabs(u)<1.0e-15 && fabs(v)<1.0e-15 && fabs(w)<1.0e-15)
                check = 0;
                
                    if((u>=0.0 && v>=0.0 && w>=0.0) || (u<0.0 && v<0.0 && w<0.0) && check==1)
                    {
                    denom = 1.0/(u+v+w);
                    u *= denom;
                    v *= denom;
                    w *= denom;
                    
                        if(domdir==1)
                        {
                        Rx = u*Ax + v*Bx + w*Cx;
                        
                        if(Rx>=Sx && Rx<=Tx)
                        ++cutnum;
                        }
                        
                        if(domdir==2 && p->j_dir==1)
                        {
                        Ry = u*Ay + v*By + w*Cy;
                        
                        if(Ry>=Sy && Ry<=Ty)
                        ++cutnum;
                        }
                        
                        if(domdir==3)
                        {
                        Rz = u*Az + v*Bz + w*Cz;
                        
                        if(Rz>=Sz && Rz<=Tz)
                        ++cutnum;
                        }
                    }
                
        
            }
                    
            
        } // ray cast loop end
        
        if(cutnum%2==0 || cutnum==0) // outside poiting -> no switch
        {
        tri_switch_local[q-tricount_local_displ[p->mpirank]] = 0;
        //cout<<"cutnum1: "<<cutnum<<endl;
        }
            
                
        if(cutnum%2!=0 && cutnum>0) // inside poiting -> switch
        {
        tri_switch_local[q-tricount_local_displ[p->mpirank]] = 1;
        //cout<<p->mpirank<<" cutnum2: "<<cutnum<<endl;
        }
	}
    
    count=0;
    for (q=tricount_local_displ[p->mpirank];q<tricount_local_displ[p->mpirank+1];++q)
    if(tri_switch_local[q-tricount_local_displ[p->mpirank]]==1)
    ++count;
    
    //cout<<p->mpirank<<"  LOCAL SWITCH: "<<count<<"  | local trinum: "<<tricount_local_list[p->mpirank]<<endl;
    
    pgc->allgatherv_int(tri_switch_local, tricount_local_list[p->mpirank], tri_switch, tricount_local_list, tricount_local_displ);
        
    // start loop part 3: global switch
        tricount_switch_total=0;
        for (int n=0;n<tricount;++n)
        if(tri_switch[n]==1)
        {   
        x1 = tri_x[n][1];
        y1 = tri_y[n][1];
        z1 = tri_z[n][1];
        
        x2 = tri_x[n][2];
        y2 = tri_y[n][2];
        z2 = tri_z[n][2];
        
        
        tri_x[n][1] = x2;
        tri_y[n][1] = y2;
        tri_z[n][1] = z2;
        
        tri_x[n][2] = x1;
        tri_y[n][2] = y1;
        tri_z[n][2] = z1;

        ++tricount_switch_total;
        }
        
        if(p->mpirank==0)
        cout<<"6DOF STL triangle switch count: "<<tricount_switch_total<<endl;
        
        // free allocated arrays
    p->del_Iarray(tri_switch,tricount);
    p->del_Iarray(tricount_local_list,p->M10+1);
    p->del_Iarray(tricount_local_displ,p->M10+1);
    p->del_Iarray(tri_switch_local,tricount_local);	
}