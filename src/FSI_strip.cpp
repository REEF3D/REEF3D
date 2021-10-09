/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2021 Tobias Martin

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
--------------------------------------------------------------------*/

#include"FSI_strip.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

fsi_strip::fsi_strip(int num):nstrip(num),beam(num)
{}
    
fsi_strip::~fsi_strip(){}

void fsi_strip::initialize(lexer *p, fdm *a, ghostcell *pgc)
{
	// Initialise parameter
    double L = p->Z11_l[nstrip];     // Length of strip [m]
    double W = p->Z11_w[nstrip];     // Width of strip [m]
    double T = p->Z11_t[nstrip];     // Thickness of strip [m]
	rho_s = p->Z11_rho[nstrip];      // Density of material [kg/m3]
	double E = p->Z11_e[nstrip];   	 // Young modulus [N/m^2]
	double Ix = p->Z11_ix[nstrip];   // X-moment of area [m^4]
	double Iy = p->Z11_iy[nstrip];   // Y-moment of area [m^4]
	double Iz = p->Z11_iz[nstrip];   // Z-moment of area [m^4]
	double Nu = p->Z11_nu[nstrip];   // Poisson ratio [-]
	Ne = p->Z11_n[nstrip];           // Number of elements

    gravity_vec << a->gi, a->gj, a->gk;
    rho_f = p->W1;
    A_el = W*T;

    // Initialise beam
    iniBeam(Ne, E, A_el, rho_s, L, E/(2.0*(1.0 + Nu)), Ix, Iy, Iz);

    // Initialise material
    iniMaterial();

    // Initialise damping and compression effects
    iniDamping(0,0,0,0,0,0,true);

    // Meshing
    Eigen::Matrix3Xd ini_coord = Eigen::Matrix3Xd::Zero(3,Ne+1); 
    for (int n = 0; n < Ne+1; n++)
    {
        ini_coord.col(n) << T/2.0, W/2.0, L/Ne*n;
    }
    meshBeam(ini_coord.row(0), ini_coord.row(1), ini_coord.row(2));

    // Initialise solver
    iniSolver();
	
	// Initialise communication 
	ini_parallel(p, a, pgc);

    // Initialise cell size
    get_cellsize(p, a, pgc);

    // Initialise Lagrangian fields
    lagrangePoints.resize(Ne);  
    lagrangeVel.resize(Ne); 
    lagrangeVelCoup.resize(Ne); 
    lagrangeForceCoup.resize(Ne); 
    lagrangeArea.resize(Ne);    
    Xil.resize(Ne); 
    Xil_0.resize(Ne);
 
    l_el = L/Ne;
    int nl = ceil(l_el/dx_body);
    int nw = ceil(W/dx_body);
    double dl = l_el/nl;
    double dw = W/nw;

    for (int n = 0; n < Ne; n++)
    {
        lagrangePoints[n] = Eigen::Matrix3Xd::Zero(3,nl*nw);   
        lagrangeVel[n] = Eigen::MatrixXd::Zero(3,nl*nw);   
        lagrangeVelCoup[n] = Eigen::MatrixXd::Zero(3,nl*nw);   
        lagrangeForceCoup[n] = Eigen::MatrixXd::Zero(3,nl*nw);   
        lagrangeArea[n] = Eigen::VectorXd::Zero(nl*nw);   
        Xil[n] = Eigen::Matrix3Xd::Zero(3,nl*nw);   
        Xil_0[n] = Eigen::Matrix3Xd::Zero(3,nl*nw);   
        
        double l_0 = n*l_el;
        Eigen::Vector3d cg; cg << T/2.0, W/2.0, 0.5*l_el + l_0; 

        int ind = 0;
        for (int ii = 0; ii < nl; ii++)
        {
            for (int jj = 0; jj < nw; jj++)
            {
                lagrangePoints[n].col(ind) << T/2.0, 0.5*dw+dw*jj, l_0 + 0.5*dl + dl*ii;
                lagrangeArea[n](ind) = dw*dl;
                Xil_0[n].col(ind) << lagrangePoints[n].col(ind) - cg;
                ind++;
            }
        }
    }

    // Initialise field
    getTransPos(x_el);
    getTransVel(xdot_el);
    omega_el = Eigen::Matrix3Xd::Zero(3,Ne+2);   
    getAngVel(omega_el);
    
    t_strip = 0.0;
    t_strip_n = 0.0;

    F_el = Eigen::Matrix3Xd::Zero(3,Ne+2);   
    P_el = Eigen::Matrix3Xd::Zero(3,Ne+2);   
    P_el_n = Eigen::Matrix3Xd::Zero(3,Ne+2);   
    M_el = Eigen::Matrix3Xd::Zero(3,Ne+2);   
    I_el = Eigen::Matrix3Xd::Zero(3,Ne+2);   
    I_el_n = Eigen::Matrix3Xd::Zero(3,Ne+2);   
    /*
    c_moor_n = Matrix3Xd::Zero(3,Ne+1); 
    cdot_moor = Matrix3Xd::Zero(3,Ne+1); 
    cdot_moor_n = Matrix3Xd::Zero(3,Ne+1);
    cdotdot_moor = Matrix3Xd::Zero(3,Ne+1);

    getTransPos(c_moor); c_moor_n = c_moor;

    vector<double> three(3, 0);
	fluid_vel.resize(Ne+1, three);
	fluid_vel_n.resize(Ne+1, three);
	fluid_acc.resize(Ne+1, three);
*/
}

void fsi_strip::ini_parallel(lexer *p, fdm *a, ghostcell *pgc)
{
	p->Darray(xstart, p->mpi_size);
	p->Darray(xend, p->mpi_size);
	p->Darray(ystart, p->mpi_size);
	p->Darray(yend, p->mpi_size);
	p->Darray(zstart, p->mpi_size);
	p->Darray(zend, p->mpi_size);
	
	xstart[p->mpirank] = p->originx;
	ystart[p->mpirank] = p->originy;
	zstart[p->mpirank] = p->originz;
	xend[p->mpirank] = p->endx;
	yend[p->mpirank] = p->endy;
	zend[p->mpirank] = p->endz;
	
	for (int i = 0; i < p->mpi_size; i++)
	{
		MPI_Bcast(&xstart[i],1,MPI_DOUBLE,i,pgc->mpi_comm);
		MPI_Bcast(&xend[i],1,MPI_DOUBLE,i,pgc->mpi_comm);
		MPI_Bcast(&ystart[i],1,MPI_DOUBLE,i,pgc->mpi_comm);
		MPI_Bcast(&yend[i],1,MPI_DOUBLE,i,pgc->mpi_comm);
		MPI_Bcast(&zstart[i],1,MPI_DOUBLE,i,pgc->mpi_comm);
		MPI_Bcast(&zend[i],1,MPI_DOUBLE,i,pgc->mpi_comm);
	}
}
	

void fsi_strip::interpolate_vel(lexer* p, fdm* a, ghostcell* pgc, field& uvel, field& vvel, field& wvel)
{
    int ii, jj, kk;
    double dx, dy, dz, dist, D;

    for (int eI = 0; eI < Ne; eI++)
    {
        lagrangeVel[eI] = Eigen::MatrixXd::Zero(3,lagrangePoints[eI].cols());   
    
        for (int pI = 0; pI < lagrangePoints[eI].cols(); pI++)
        {
            const Eigen::Vector3d& coordI = lagrangePoints[eI].col(pI);

            if 
            (
                coordI(0) >= xstart[p->mpirank] && coordI(0) < xend[p->mpirank] &&
                coordI(1) >= ystart[p->mpirank] && coordI(1) < yend[p->mpirank] &&
                coordI(2) >= zstart[p->mpirank] && coordI(2) < zend[p->mpirank]
            )
            {
                ii = p->posc_i(coordI(0));
                jj = p->posc_j(coordI(1));
                kk = p->posc_k(coordI(2));
                
                dx = p->DXN[ii + marge];
                dy = p->DYN[jj + marge];
                dz = p->DZN[kk + marge];

                for (int i_it = ii - 2; i_it <= ii + 2; i_it++)
                {
                    for (int j_it = jj - 2; j_it <= jj + 2; j_it++)
                    {
                        for (int k_it = kk - 2; k_it <= kk + 2; k_it++)
                        {
                            dist = (p->XN[i_it + 1 + marge] - coordI(0))/dx;
                            D = kernel_roma(dist);
                            dist = (p->YP[j_it + marge] - coordI(1))/dy;
                            D *= kernel_roma(dist);
                            dist = (p->ZP[k_it + marge] - coordI(2))/dz;
                            D *= kernel_roma(dist);
                            
                            lagrangeVel[eI](0,pI) += uvel(i_it,j_it,k_it)*D;

                            dist = (p->XP[i_it + marge] - coordI(0))/dx;
                            D = kernel_roma(dist);
                            dist = (p->YN[j_it + 1 + marge] - coordI(1))/dy;
                            D *= kernel_roma(dist);
                            dist = (p->ZP[k_it + marge] - coordI(2))/dz;
                            D *= kernel_roma(dist);
                                
                            lagrangeVel[eI](1,pI) += vvel(i_it,j_it,k_it)*D;
                            
                            dist = (p->XP[i_it + marge] - coordI(0))/dx;
                            D = kernel_roma(dist);
                            dist = (p->YP[j_it + marge] - coordI(1))/dy;
                            D *= kernel_roma(dist);
                            dist = (p->ZN[k_it + 1 + marge] - coordI(2))/dz;
                            D *= kernel_roma(dist);
                             
                            lagrangeVel[eI](2,pI) += wvel(i_it,j_it,k_it)*D;
                        }
                    }
                }
            }
        }
        
        for (int pI = 0; pI < lagrangePoints[eI].cols(); pI++)
        {
            lagrangeVel[eI].col(pI) << pgc->globalsum(lagrangeVel[eI](0,pI)), pgc->globalsum(lagrangeVel[eI](1,pI)), pgc->globalsum(lagrangeVel[eI](2,pI));
        }
    }
/*
    char name[100];
    sprintf(name,"./REEF3D-CFD-Lagrange.csv");
    ofstream result;
    result.open(name, ios::binary);

    for (int eI = 0; eI < lagrangePoints.size(); eI++)
    {
        for (int pI = 0; pI < lagrangePoints[eI].cols(); pI++)
        {
            const Eigen::Vector3d& coordI = lagrangePoints[eI].col(pI);
            const Eigen::Vector3d& velI = lagrangeVel[eI].col(pI);
        
            result<<coordI(0)<<","<<coordI(1)<<","<<coordI(2)<<","<<velI(0)<<","<<velI(1)<<","<<velI(2)<<endl;
        }
    }
    result.close();
*/
}
    
double fsi_strip::kernel_roma(const double& dist)
{
    double D = 0.0;

    if (fabs(dist) <= 0.5)
    {
        D = 1.0/3.0*(1.0 + sqrt(-3*dist*dist + 1));
    }
    else if (fabs(dist) <= 1.5)
    {    
        D = 1.0/6.0*(5.0 - 3.0*fabs(dist) - sqrt(-3*(1 - fabs(dist))*(1 - fabs(dist)) + 1));
    }
    
    return D;
}

void fsi_strip::coupling_vel()
{
    getTransVel(xdot_el);
    getAngVel(omega_el);

    for (int eI = 0; eI < Ne; eI++)
    {
        for (int pI = 0; pI < lagrangePoints[eI].cols(); pI++)
        {
            lagrangeVelCoup[eI].col(pI) = (xdot_el.col(eI+1) + xdot_el.col(eI))/2.0 + omega_el.col(eI+1).cross(Xil[eI].col(pI));
        }
    }
}

void fsi_strip::coupling_force(lexer *p, double alpha)
{
    for (int eI = 0; eI < Ne; eI++)
    {
        lagrangeForceCoup[eI] = (lagrangeVelCoup[eI] - lagrangeVel[eI])/(alpha*p->dt);
    }
}

void fsi_strip::distribute_forces(lexer *p, fdm *a, ghostcell *pgc, field1& fx, field2& fy, field3& fz)
{
    int ii, jj, kk;
    double dx, dy, dz, dV, dist, D;

    for (int eI = 0; eI < Ne; eI++)
    {
        for (int pI = 0; pI < lagrangePoints[eI].cols(); pI++)
        {
            const Eigen::Vector3d& coordI = lagrangePoints[eI].col(pI);
            const Eigen::Vector3d& forceI = lagrangeForceCoup[eI].col(pI);
            const double& areaI = lagrangeArea[eI](pI);

            ii = p->posc_i(coordI(0));
            jj = p->posc_j(coordI(1));
            kk = p->posc_k(coordI(2));

            dx = p->DXN[ii + marge];
            dy = p->DYN[jj + marge];
            dz = p->DZN[kk + marge];

            for (int i_it = ii - 2; i_it <= ii + 2; i_it++)
            {
                for (int j_it = jj - 2; j_it <= jj + 2; j_it++)
                {
                    for (int k_it = kk - 2; k_it <= kk + 2; k_it++)
                    {
                        dV = areaI*dx_body;

                        dist = (p->XN[i_it + 1 + marge] - coordI(0))/dx;
                        D = kernel_roma(dist);
                        dist = (p->YP[j_it + marge] - coordI(1))/dy;
                        D *= kernel_roma(dist);
                        dist = (p->ZP[k_it + marge] - coordI(2))/dz;
                        D *= kernel_roma(dist);
                        
                        fx(i_it,j_it,k_it) += forceI(0)*D*dV/(dx*dy*dz);
      

                        dist = (p->XP[i_it + marge] - coordI(0))/dx;
                        D = kernel_roma(dist);
                        dist = (p->YN[j_it + 1 + marge] - coordI(1))/dy;
                        D *= kernel_roma(dist);
                        dist = (p->ZP[k_it + marge] - coordI(2))/dz;
                        D *= kernel_roma(dist);
                        
                        fy(i_it,j_it,k_it) += forceI(1)*D*dV/(dx*dy*dz);
                

                        dist = (p->XP[i_it + marge] - coordI(0))/dx;
                        D = kernel_roma(dist);
                        dist = (p->YP[j_it + marge] - coordI(1))/dy;
                        D *= kernel_roma(dist);
                        dist = (p->ZN[k_it + 1 + marge] - coordI(2))/dz;
                        D *= kernel_roma(dist);
                        
                        fz(i_it,j_it,k_it) += forceI(2)*D*dV/(dx*dy*dz);
                    }
                }
            }
        }
    }
    
    pgc->start1(p,fx,10);
    pgc->start2(p,fy,11);
    pgc->start3(p,fz,12);

    ULOOP
    a->test(i,j,k) = fx(i,j,k); 
    
    pgc->start4(p,a->test,10);
}

void fsi_strip::get_cellsize(lexer *p, fdm *a, ghostcell *pgc)
{
    Eigen::Vector3d coordI;
    coordI << p->Z11_t[nstrip]/2.0, p->Z11_w[nstrip]/2.0, p->Z11_l[nstrip]/2.0;
    
    if 
    (
        coordI(0) >= xstart[p->mpirank] && coordI(0) < xend[p->mpirank] &&
        coordI(1) >= ystart[p->mpirank] && coordI(1) < yend[p->mpirank] &&
        coordI(2) >= zstart[p->mpirank] && coordI(2) < zend[p->mpirank]
    )
    {
        int ii = p->posc_i(coordI(0));
        int jj = p->posc_j(coordI(1));
        int kk = p->posc_k(coordI(2));
        
        dx_body = p->DXN[ii + marge];
    }
    else
    {
         dx_body = 0.0;
    }

    dx_body = pgc->globalsum(dx_body);
}

void fsi_strip::update_points()
{
    for (int eI = 0; eI < lagrangePoints.size(); eI++)
    {
        for (int pI = 0; pI < lagrangePoints[eI].cols(); pI++)
        {
            Xil[eI].col(pI) = Xil_0[eI].col(pI); //rotVec(Xil_0[eI].col(pI),eI+1);
            lagrangePoints[eI].col(pI) = (x_el.col(eI+1) + x_el.col(eI))/2.0 + Xil[eI].col(pI);
        }
    }
cout<<"no transformation!"<<endl;
}


void fsi_strip::start(lexer *p, fdm *a, ghostcell *pgc, double alpha)
{
	// Set mooring time step
	double phi_strip = 0.0;
	t_strip_n = t_strip;
	t_strip = phi_strip*p->simtime + (1.0 - phi_strip)*(p->simtime + alpha*p->dt);

    // Integrate from t_mooring_n to t_mooring
    Integrate(t_strip_n,t_strip);
};

void fsi_strip::setFieldBC(Matrix3Xd& c_, Matrix3Xd& cdot_, Matrix4Xd& q_, Matrix4Xd& q0_, Matrix4Xd& qdot_, Matrix3Xd& f_, Matrix4Xd& m0_, Matrix3Xd& rhs_cdot_, double time , int ind)
{
    if (ind == 0){}
    else if (ind == 1)
    {
        // BC: Free translatory end with vanishing forces
        f_.col(Ne+1) = -f_.col(Ne); // correct?
    }
    else if (ind == 2)
    {
        // BC: Free rotatory end with vanishing moments 
        m0_.col(Ne) = Eigen::Vector4d::Zero(4); q_.col(Ne+1) = q_.col(Ne); q0_.col(Ne+1) = q0_.col(Ne);
    }
    else if (ind == 3)
    {
        // BC: Fixed translatory end
        rhs_cdot_.col(0) = Eigen::Vector3d::Zero(3);
    }
}


void fsi_strip::setConstantLoads(Matrix3Xd& Fext_, Matrix4Xd& Mext_, const Matrix3Xd& c_, const Matrix3Xd& cdot_, const Matrix4Xd& q_, const Matrix4Xd& qdot_)
{
}


void fsi_strip::setVariableLoads(Matrix3Xd& Fext_, Matrix4Xd& Mext_, const Matrix3Xd& c_, const Matrix3Xd& cdot_, const Matrix4Xd& q_, const Matrix4Xd& qdot_, const double time)
{
    double dm_el, m_el;
    Eigen::Vector3d P_el_star,I_el_star,s0;
    Eigen::Matrix3d J0,Xil_0_skew;

    for (int eI = 1; eI < Ne+1; eI++)
    {
        m_el = 0.0;   
        P_el.col(eI) << 0.0, 0.0, 0.0;
        P_el_star << 0.0, 0.0, 0.0;
        I_el_star << 0.0, 0.0, 0.0;
        s0 << 0.0, 0.0, 0.0;
        J0 << Eigen::Matrix3d::Zero();
        for (int pI = 0; pI < lagrangePoints[eI-1].cols(); pI++)
        {
            // Mass of element
            dm_el = rho_f*dx_body*lagrangeArea[eI-1](pI);
            m_el += dm_el;

            // Preliminary linear momentum
            P_el_star += dm_el*lagrangeVel[eI-1].col(pI);
    
            // Static moment
            s0 += dm_el*Xil_0[eI-1].col(pI);
            
            // Preliminary angular momentum
            I_el_star += dm_el*Xil[eI-1].col(pI).cross(lagrangeVel[eI-1].col(pI));

            // Quaternionic tensor of inertia
            Xil_0_skew << 0, -Xil_0[eI-1](2,pI), Xil_0[eI-1](1,pI), Xil_0[eI-1](2,pI), 0, -Xil_0[eI-1](0,pI), -Xil_0[eI-1](1,pI), Xil_0[eI-1](0,pI), 0;
            Xil_0_skew = Xil_0_skew.transpose()*Xil_0_skew;
            J0 += dm_el*Xil_0_skew;
        }

        // Determine linear momentum
        P_el.col(eI) = m_el*(cdot_.col(eI-1)+cdot_.col(eI))/2.0 + omega_el.col(eI).cross(rotVec(s0,eI));;

        // Determine coupling force
        F_el.col(eI) = -(P_el.col(eI) - P_el_n.col(eI))/(time - t_strip_n) + (P_el_n.col(eI) - P_el_star)/(t_strip - t_strip_n);

        // Determine angular momentum
        Eigen::Vector3d omega_0;
        I_el.col(eI) = rotVec(s0,eI).cross((cdot_.col(eI-1)+cdot_.col(eI))/2.0) + rotVec(J0*omega_0,eI);

        // Determine coupling moment
        M_el.col(eI) = -(I_el.col(eI) - I_el_n.col(eI))/(time - t_strip_n) + (I_el_n.col(eI) - I_el_star)/(t_strip - t_strip_n);
    }
    
    // Assign external forces
    for (int eI = 0; eI < Ne+1; eI++)
    {
        Fext_.col(eI) = (F_el.col(eI) + F_el.col(eI+1))/(2.0*rho_s*A_el*l_el) + (1.0 - rho_f/rho_s)*gravity_vec;
    }
    
    // Assign external moments
    for (int eI = 0; eI < Ne+2; eI++)
    {
        Mext_.col(eI) << 0.0, M_el.col(eI)/l_el;
    }
}
