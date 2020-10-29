/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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
along with this nuogram; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/
/*
#include"6DOF_f.h"
#include"net_QuasiStatic.h"
#include"net_void.h"
#include"mooring_void.h"
#include"mooring_DGSEM.h"
#include"mooring_QuasiStatic.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"reinidisc.h"
#include"reinidisc_f.h"
#include"reinidisc_fsf.h"
#include"momentum_FSI.h"

#include <Eigen/Dense>

void sixdof_f::elastic_beam_ini
(
	lexer *p,
	fdm *a, 
	ghostcell *pgc
)
{
    //- Input
    
    N = 10;                // Number of elements

    double D = 0.1;           // Outer diameter
    double th = 0.00297;      // Thickness of pipe

    double E = 0.8e9;         // Elasticity modulus (0.8e9 for HDPE, 210.0e9 for steel)
    double nu = 0.45;         // Poisson ratio (0.45 for HDPE, 0.3 for steel)
    
    double rho = 940.0;      // Density (940 for HDPE, 7850 for steel)


    double max_z = -10000.0;
    double min_z = 10000.0;

	for(n=0; n<tricount; ++n)
    {
        for(q=0; q<3; ++q)
        {
            if (tri_z[n][q] > max_z) max_z = tri_z[n][q];
            if (tri_z[n][q] < min_z) min_z = tri_z[n][q];
        }
    }

//---------------------------------------------

    //- Geometrical calculations
    
    Eigen::Vector3d xs(xg, yg, min_z);  // Start coordinates
    Eigen::Vector3d xe(xg, yg, max_z);  // End coordinates

    double L = xe(2) - xs(2);       // Vertical cylinder assumed
    double Le = L/N;                // Element length
    
    ndof = (N + 1)*6;           // Number of degrees of freedom    

    zcoord.resize(N+1);    // z-coordinates of nodes
    for (int i = 0; i < N+1; i++)
    {
        zcoord(i) = xs(2) + i*Le;   // coord of element n is between z(n) and z(n+1)
    };

    double Di = D - 2.0*th;         // Inner diameter
    double Ar = PI*(pow(0.5*D, 2.0) - pow(0.5*Di, 2.0));    // Area of cross section
    
    double G = 0.5*E/(1.0 + nu);    // shear modulus in [Pa]

    // Torsion constant, 2nd-order moments of area, polar moment 
    double It = 2.0/3.0*PI*0.5*(D + Di)/2.0*th*th*th;           // for thin-walled cylinder 
    double Iz = PI/4.0*(pow(D/2.0, 4.0) - pow(Di/2.0, 4.0));       
    double Iy = PI/4.0*(pow(D/2.0, 4.0) - pow(Di/2.0, 4.0));
    double Ip = PI/2.0*(pow(D/2.0, 4.0) - pow(Di/2.0, 4.0));  
    
    // shear parameters (= 0 for Bernoulli beam)
    double Asy = Ar;                                // -> How to calculate this?
    double Asz = Ar;

    double Py = 0.0; // 12.0*E*Iz/(G*Asy*Le*Le);
    double Pz = 0.0; // 12.0*E*Iy/(G*Asz*Le*Le);

    //- Ini matrices
    invA.resize(ndof, ndof);
    B.resize(ndof);
    
    F_el.resize(N);
    x.resize(ndof);
    xdot.resize(ndof); 
    xdotdot.resize(ndof); 
    x_new.resize(ndof); 
    xdot_new.resize(ndof); 
    xdotdot_new.resize(ndof);



    
    //- Local mass and stiffness matrix

    Eigen::MatrixXd Ke11(6, 6);
    Eigen::MatrixXd Ke21(6, 6);
    Eigen::MatrixXd Ke22(6, 6);
    Eigen::MatrixXd Ke(12, 12);

    Ke11(0,0) = E*Ar*Le*Le;
    Ke11(1,1) = 12.0*E*Iz/(1.0 + Py);
    Ke11(1,5) = 6.0*E*Iz*Le/(1.0 + Py);
    Ke11(2,2) = 12.0*E*Iy/(1.0 + Pz);
    Ke11(2,4) = -6.0*E*Iy*Le/(1.0 + Pz);
    Ke11(3,3) = G*It*Le*Le;
    Ke11(4,2) = -6.0*E*Iy*Le/(1.0 + Pz);
    Ke11(4,4) = E*Iy*Le*Le*(4.0 + Pz)/(1.0 + Pz);
    Ke11(5,1) = 6.0*E*Iz*Le/(1.0 + Py);
    Ke11(5,5) = E*Iz*Le*Le*(4.0 + Py)/(1.0 + Py);

    Ke22(0,0) = E*Ar*Le*Le;
    Ke22(1,1) = 12.0*E*Iz/(1.0 + Py);
    Ke22(1,5) = -6.0*E*Iz*Le/(1.0 + Py);
    Ke22(2,2) = 12.0*E*Iy/(1.0 + Pz);
    Ke22(2,4) = 6.0*E*Iy*Le/(1.0 + Pz);
    Ke22(3,3) = G*It*Le*Le;
    Ke22(4,2) = 6.0*E*Iy*Le/(1.0 + Pz);
    Ke22(4,4) = E*Iy*Le*Le*(4.0 + Pz)/(1.0 + Pz);
    Ke22(5,1) = -6.0*E*Iz*Le/(1.0 + Py);
    Ke22(5,5) = E*Iz*Le*Le*(4.0 + Py)/(1.0 + Py);
    
    Ke21(0,0) = -E*Ar*Le*Le; 
    Ke21(1,1) = -12.0*E*Iz/(1.0 + Py);
    Ke21(1,5) = -6.0*E*Iz*Le/(1.0 + Py);
    Ke21(2,2) = -12.0*E*Iy/(1.0 + Pz);
    Ke21(2,4) = 6.0*E*Iy*Le/(1.0 + Pz);
    Ke21(3,3) = -G*It*Le*Le;
    Ke21(4,2) = -6.0*E*Iy*Le/(1.0 + Pz);
    Ke21(4,4) = E*Iy*Le*Le*(2.0 - Pz)/(1.0 + Pz);
    Ke21(5,1) = 6.0*E*Iz*Le/(1.0 + Py);
    Ke21(5,5) = E*Iz*Le*Le*(2.0 - Py)/(1.0 + Py);

    Ke << Ke11, Ke21.transpose(), 
          Ke21, Ke22;

    Ke *= 1.0/(Le*Le*Le);


    Eigen::MatrixXd Me11(6, 6);
    Eigen::MatrixXd Me21(6, 6);
    Eigen::MatrixXd Me22(6, 6);
    Eigen::MatrixXd Me(12, 12);
    
    Me11(0,0) = 1.0/3.0;
    Me11(1,1) = 13.0/35.0 + 6.0/5.0*Iz/(Ar*Le*Le);
    Me11(1,5) = 11.0/210.0*Le + 1.0/10.0*Iz/(Ar*Le);
    Me11(2,2) = 13.0/35.0 + 6.0/5.0*Iy/(Ar*Le*Le);
    Me11(2,4) = -11.0/210.0*Le - 1.0/10.0*Iy/(Ar*Le);
    Me11(3,3) = 1.0/3.0*Ip/Ar;
    Me11(4,2) = -11.0/210.0*Le - 1.0/10.0*Iy/(Ar*Le);
    Me11(4,4) = Le*Le/105.0 + 2.0/15.0*Iy/Ar;
    Me11(5,1) = 11.0/210.0*Le + 1.0/10.0*Iz/(Le*Ar);
    Me11(5,5) = 1.0/105.0*Le*Le + 2.0/15.0*Iz/Ar;
    
    Me22(0,0) = 1.0/3.0;
    Me22(1,1) = 13.0/35.0 + 6.0/5.0*Iz/(Ar*Le*Le);
    Me22(1,5) = -11.0/210.0*Le - 1.0/10.0*Iz/(Ar*Le);
    Me22(2,2) = 13.0/35.0 + 6.0/5.0*Iy/(Ar*Le*Le);
    Me22(2,4) = 11.0/210.0*Le + 1.0/10.0*Iy/(Ar*Le);
    Me22(3,3) = 1.0/3.0*Ip/Ar;
    Me22(4,2) = 11.0/210.0*Le + 1.0/10.0*Iy/(Ar*Le);
    Me22(4,4) = Le*Le/105.0 + 2.0/15.0*Iy/Ar;
    Me22(5,1) = -11.0/210.0*Le - 1.0/10.0*Iz/(Le*Ar);
    Me22(5,5) = 1.0/105.0*Le*Le + 2.0/15.0*Iz/Ar;

    Me21(0,0) = 1.0/6.0;
    Me21(1,1) = 9.0/70.0 - 6.0/5.0*Iz/(Ar*Le*Le);
    Me21(1,5) = 13.0/420.0*Le - 1.0/10.0*Iz/(Ar*Le);
    Me21(2,2) = 9.0/70.0 - 6.0/5.0*Iy/(Ar*Le*Le);
    Me21(2,4) = -13.0/420.0*Le + 1.0/10.0*Iy/(Ar*Le);
    Me21(3,3) = 1.0/6.0*Ip/Ar;
    Me21(4,2) = 13.0/420.0*Le - 1.0/10.0*Iy/(Ar*Le);
    Me21(4,4) = -Le*Le/140.0 - 2.0/30.0*Iy/Ar;
    Me21(5,1) = -13.0/420.0*Le + 1.0/10.0*Iz/(Ar*Le);
    Me21(5,5) = -Le*Le/140.0 - 2.0/30.0*Iz/Ar;

    Me << Me11, Me21.transpose(), 
          Me21, Me22;

    Me *= rho*Ar*Le;
    
    //- Global mass and stiffness matrices
    
    M = Eigen::MatrixXd::Zero(ndof, ndof);
    K = Eigen::MatrixXd::Zero(ndof, ndof);
    
    for (int i = 1; i <= N; i++)
    {
        int ind_j = 0;
        for (int j = 6*(i-1); j < 6*(i + 1); j++)                             
        {
            int ind_k = 0;
            for (int k = 6*(i-1); k < 6*(i + 1); k++)
            {
                K(j, k) += Ke(ind_j, ind_k); 
                M(j, k) += Me(ind_j, ind_k);

                ind_k++;
            }

            ind_j++;
        }
    } 


    // Damping matrix
    C = Eigen::MatrixXd::Zero(ndof, ndof);
    double alpha1 = 0.0;
    double alpha2 = 0.0;
    C = alpha1*M + alpha2*K;


    // Boundary condition, left end Dirichlet: strike rows and colums, and set diagonal to 1
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < ndof; j++)
        {
            K(i,j) = 0.0;
            K(j,i) = 0.0;
            M(i,j) = 0.0;
            M(j,i) = 0.0;
            C(i,j) = 0.0;
            C(j,i) = 0.0;
        }
        
        K(i,i) = 1.0;
        M(i,i) = 1.0;
        C(i,i) = 1.0;
    }

    // Initial conditions
    x = Eigen::VectorXd::Zero(ndof);
    xdot = Eigen::VectorXd::Zero(ndof);
    xdotdot = Eigen::VectorXd::Zero(ndof);
    F = Eigen::VectorXd::Zero(ndof);
    xdotdot = M.inverse()*(F - C*xdot - K*x);
    
    elastic_beam_forces(p, a, pgc);
    elastic_beam_print(p);
}


void sixdof_f::elastic_beam_start
(
	lexer *p,
	fdm *a, 
	ghostcell *pgc
)
{
    //- Global force vector
    F = Eigen::VectorXd::Zero(ndof);
    Eigen::VectorXd Fe(12);
   
    elastic_beam_forces(p, a, pgc);
    
    for (int i = 1; i <= N; i++)
    {
        Fe << 0.0, 0.0, 0.5*F_knot(i-1), 0.0, 0.0, 0.0, 
              0.0, 0.0, 0.5*F_knot(i),   0.0, 0.0, 0.0;
        
        // Assemble global matrix
        int ind_j = 0;
        for (int j = 6*(i-1); j < 6*(i + 1); j++)                             
        {
            F(j) += Fe(ind_j);
            
            ind_j++;
        }
    } 

    //- Run-time process
    double gamma = 0.5; 
    double beta = 0.25;
    
    invA = ((1.0/(beta*p->dt*p->dt))*M + (gamma/(beta*p->dt))*C + K).inverse();


    B = F 
        + M*((1.0/(beta*p->dt*p->dt))*x + (1.0/(beta*p->dt))*xdot + (1.0/(2.0*beta)-1.0)*xdotdot)
        + C*((gamma/(beta*p->dt))*x + (gamma/beta-1.0)*xdot + (gamma/beta-2.0)*(p->dt/2.0)*xdotdot);

    x_new = invA*B;
        
    xdotdot_new = (1.0/(beta*p->dt*p->dt))*(x_new - x) - (1.0/(beta*p->dt))*xdot - ((1.0/(2.0*beta)) - 1.0)*xdotdot;

    xdot_new = xdot + (1.0 - gamma)*p->dt*xdotdot + gamma*p->dt*xdotdot_new;
        
    x = x_new;
    xdot = xdot_new;
    xdotdot = xdotdot_new;


    //- Print deformation
    print_E_deformation(p,a,pgc);
    elastic_beam_print(p);


    //- Velocity boundary condition to fluid
    pgc->setfbvel1(p,zcoord,xdot);
}


void sixdof_f::elastic_beam_forces
(
	lexer *p,
	fdm *a, 
	ghostcell *pgc
)
{
	double x0,x1,x2,y0,y1,y2,z0,z1,z2;
	double xc,yc,zc;
	double at,bt,ct,st;
	double nx,ny,nz,norm;
	double A_triang,A;
	double p_int,rho_int,nu_int,u_int,v_int,w_int;
	double du,dv,dw;
	double xlocvel,ylocvel,zlocvel;
	double Fx,Fy,Fz;

    pgc->start4(p,a->press,40);
	pgc->start4(p,a->press,401);
    pgc->start4(p,a->press,401);

    pgc->gcfb_update_extra_gcb(p,a,a->press);
    
    pgc->dgcpol(p,a->press,p->dgc4,p->dgc4_count,14);
    a->press.ggcpol(p);
    
	
    A=0.0;
            
    F_el = Eigen::VectorXd::Zero(N);

    for (int n = 0; n < tricount; ++n)
    {     
		// Vertices of triangle
	
		x0 = tri_x[n][0];
		x1 = tri_x[n][1];
		x2 = tri_x[n][2];
		
		y0 = tri_y[n][0];
		y1 = tri_y[n][1];
		y2 = tri_y[n][2];
		
		z0 = tri_z[n][0];
		z1 = tri_z[n][1];
		z2 = tri_z[n][2];  
           
		   
		// Center of triangle
		
		xc = (x0 + x1 + x2)/3.0;
		yc = (y0 + y1 + y2)/3.0;
		zc = (z0 + z1 + z2)/3.0;
            
		if 
		(
			xc >= p->originx && xc < p->endx &&
			yc >= p->originy && yc < p->endy &&
			zc >= p->originz && zc < p->endz
		)
		{
			// Area of triangle using Heron's formula

			at = sqrt(pow(x1-x0,2.0) + pow(y1-y0,2.0) + pow(z1-z0,2.0));
			bt = sqrt(pow(x1-x2,2.0) + pow(y1-y2,2.0) + pow(z1-z2,2.0));
			ct = sqrt(pow(x2-x0,2.0) + pow(y2-y0,2.0) + pow(z2-z0,2.0));
				
			st = 0.5*(at+bt+ct);
				
			A_triang = sqrt(MAX(0.0,st*(st-at)*(st-bt)*(st-ct)));
				

			// Normal vectors (always pointing outwards)      
				
			nx = (y1 - y0) * (z2 - z0) - (y2 - y0) * (z1 - z0);
			ny = (x2 - x0) * (z1 - z0) - (x1 - x0) * (z2 - z0); 
			nz = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);

			norm = sqrt(nx*nx + ny*ny + nz*nz);
			
			nx /= norm > 1.0e-20 ? norm : 1.0e20;
			ny /= norm > 1.0e-20 ? norm : 1.0e20;
			nz /= norm > 1.0e-20 ? norm : 1.0e20;
			
			   
			// Interpolate values at centre of triangle

			p_int = p->ccipol4_a(a->press,xc,yc,zc);
			nu_int = p->ccipol4_a(a->visc,xc,yc,zc);
			rho_int = p->ccipol4_a(a->ro,xc,yc,zc);
			
			
			// Interpolate velocities in fluid near centre of triangle
			
			i = p->posc_i(xc);
			j = p->posc_j(yc);
			k = p->posc_k(zc);
	 
			xlocvel = xc + nx*p->DXP[IP];
			ylocvel = yc + ny*p->DYP[JP];
			zlocvel = zc + nz*p->DZP[KP];
			
			u_int = p->ccipol1_a(a->u,xlocvel,ylocvel,zlocvel);
			v_int = p->ccipol2_a(a->v,xlocvel,ylocvel,zlocvel);
			w_int = p->ccipol3_a(a->w,xlocvel,ylocvel,zlocvel);
				
			du = u_int/p->DXP[IP];
			dv = v_int/p->DYP[JP];
			dw = w_int/p->DZP[KP];
				

			// Calculate forces on triangle
				
			Fx = -nx*p_int*A_triang + rho_int*nu_int*A_triang*(du*ny+du*nz);
			Fy = -ny*p_int*A_triang + rho_int*nu_int*A_triang*(dv*nx+dv*nz);                    
            Fz = -nz*p_int*A_triang + rho_int*nu_int*A_triang*(dw*nx+dw*ny);
			

			// Add forces to element
		
            for (int el = 0; el < N; el++)
            {
                if (zc <= zcoord(el+1))
                {
                    F_el(el) += Fx;
                    A += A_triang;
                    break;
                }
            }
		}
	}		
 
	// Communication with other processors
	
    A = pgc->globalsum(A);
	
    for (int el = 0; el < N; el++)
    {
        F_el(el) = pgc->globalsum(F_el(el));
    }

    F_knot = Eigen::VectorXd::Zero(N + 1);

    for (int el = 1; el < N - 1; el++)
    {
        F_knot(el) = 0.5*F_el(el - 1) + 0.5*F_el(el);
    }
    F_knot(0) = 0.0; // 0.5*F_el(0);
    F_knot(N) = 0.5*F_el(N - 1);


	// Print results
	
	if (p->mpirank == 0)
	{
		cout<<"Beam forces = "<<F_knot.transpose()<<endl; 
	}
}


void sixdof_f::elastic_beam_print(lexer *p)
{
	int num=0;
	
	if(p->P15==1)
    num = p->printcount_sixdof;

    if(p->P15==2)
    num = p->count;
	
	if(num<0)
	num=0;
    
    if
	(
		p->mpirank==0 && (((p->count%p->P20==0) && p->P30<0.0)  
		|| (p->simtime>printtime && p->P30>0.0)   
		|| p->count==0)
	)
	{

        if(p->P14==1)
		{
			if(num<10)
			sprintf(name,"./REEF3D_6DOF_Elastic/REEF3D-Elastic-00000%d.vtk",num);

			if(num<100&&num>9)
			sprintf(name,"./REEF3D_6DOF_Elastic/REEF3D-Elastic-0000%d.vtk",num);

			if(num<1000&&num>99)
			sprintf(name,"./REEF3D_6DOF_Elastic/REEF3D-Elastic-000%d.vtk",num);

			if(num<10000&&num>999)
			sprintf(name,"./REEF3D_6DOF_Elastic/REEF3D-Elastic-00%d.vtk",num);

			if(num<100000&&num>9999)
			sprintf(name,"./REEF3D_6DOF_Elastic/REEF3D-Elastic-0%d.vtk",num);

			if(num>99999)
			sprintf(name,"./REEF3D_6DOF_Elastic/REEF3D-Elastic-%d.vtk",num);
		}	

		ofstream result;
		result.open(name, ios::binary);
		
		result << "# vtk DataFile Version 2.0" << endl;
		result << "Elastic beam"<< endl;
        result << "ASCII \nDATASET UNSTRUCTURED_GRID" << endl;
		result << "POINTS " << N + 1 << " float" <<endl;
		
		for(int n = 0; n < N + 1; ++n)
		{
			result<<p->X133_xc + x(2 + n*6)<<" "<<p->X133_yc<<" "<<zcoord(n)<<endl;
		}
		
		result << "\nCELLS " << N << " " << N*3 <<endl;	
		
		for(int n = 0; n < N; ++n)
		{
			result<<"2 "<< n << " " << n+1 << endl;
		}
		
		result << "\nCELL_TYPES " << N << endl;	
		
		for(int n = 0; n < N; ++n)
		{
			result<<"3"<<endl;
		}	

		result<<"\nPOINT_DATA " << N + 1 <<endl;
		
        result<<"SCALARS Displacement float 1 \nLOOKUP_TABLE default"<<endl;
        
        for (int n = 0; n < N + 1; n++)
        {
             result<<x(2 + n*6)<<endl;
        }
        
        result<<"SCALARS Rotation float 1 \nLOOKUP_TABLE default"<<endl;
        
        for (int n = 0; n < N + 1; n++)
        {
             result<<x(4 + n*6)<<endl;
        }
        
        result<<"SCALARS Velocity-x float 1 \nLOOKUP_TABLE default"<<endl;
        
        for (int n = 0; n < N + 1; n++)
        {
             result<<xdot(2 + n*6)<<endl;
        }
        
        result<<"SCALARS Acceleration-x float 1 \nLOOKUP_TABLE default"<<endl;
        
        for (int n = 0; n < N + 1; n++)
        {
             result<<xdotdot(2 + n*6)<<endl;
        }

        result<<"SCALARS Forces float 1 \nLOOKUP_TABLE default"<<endl;
        result<<F_knot<<endl;
		
        result.close();
	}
}
*/


