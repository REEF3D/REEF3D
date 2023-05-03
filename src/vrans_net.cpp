/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2023 Tobias Martin

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

#include"vrans_net.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"net.h"

vrans_net::vrans_net(lexer *p, ghostcell *pgc) 
: Fx_net(p), Fy_net(p), Fz_net(p),kernel_x(p), kernel_y(p), kernel_z(p)
{
    //initialize(p,a,pgc);
}

vrans_net::~vrans_net(){}

void vrans_net::start(lexer *p, fdm *a, ghostcell *pgc, net *&ppnet, int nNet)
{ 
    // Distribute net forces on surrounding cells
    int ii, jj, kk;
    double dist, D, dx, dy, dz;
    double test = 0.0;

    const EigenMat& lagrangePoints = ppnet->getLagrangePoints();
    const EigenMat& lagrangeForces = ppnet->getLagrangeForces();

    ULOOP
    {
        kernel_x(i,j,k) = 0.0;
    }
    VLOOP
    {
        kernel_y(i,j,k) = 0.0;
    }
    WLOOP
    {
        kernel_z(i,j,k) = 0.0;
    }
    
    if (nNet == 0)
    {
        ULOOP
        {
            Fx_net(i,j,k) = 0.0;
        }        
        VLOOP
        {
            Fy_net(i,j,k) = 0.0;
        }        
        WLOOP
        {
            Fz_net(i,j,k) = 0.0;
        }        
    }
    

    for (int pI = 0; pI < lagrangePoints.size(); pI++)
    {
        const Eigen::Vector3d& forcesI = lagrangeForces[pI];
        
        if (forcesI.norm() != 0.0)
        {
            const Eigen::Vector3d& coordI = lagrangePoints[pI];

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
                        dist = sqrt(pow(p->XN[i_it + 1 + marge] - coordI(0), 2.0))/dx;
                        D = kernel_peskin(dist);
                        dist = sqrt(pow(p->YP[j_it + marge] - coordI(1), 2.0))/dy;
                        D *= kernel_peskin(dist);
                        dist = sqrt(pow(p->ZP[k_it + marge] - coordI(2), 2.0))/dz;
                        D *= kernel_peskin(dist);
                        
                        test = forcesI(0)*D/(dx*dy*dz);
                        Fx_net(i_it,j_it,k_it) += test;

                        dist = sqrt(pow(p->XP[i_it + marge] - coordI(0), 2.0))/dx;
                        D = kernel_peskin(dist);
                        dist = sqrt(pow(p->YN[j_it + 1 + marge] - coordI(1), 2.0))/dy;
                        D *= kernel_peskin(dist);
                        dist = sqrt(pow(p->ZP[k_it + marge] - coordI(2), 2.0))/dz;
                        D *= kernel_peskin(dist);
                            
                        test = forcesI(1)*D/(dx*dy*dz);
                        Fy_net(i_it,j_it,k_it) += test;
                        
                        dist = sqrt(pow(p->XP[i_it + marge] - coordI(0), 2.0))/dx;
                        D = kernel_peskin(dist);
                        dist = sqrt(pow(p->YP[j_it + marge] - coordI(1), 2.0))/dy;
                        D *= kernel_peskin(dist);
                        dist = sqrt(pow(p->ZN[k_it + 1 + marge] - coordI(2), 2.0))/dz;
                        D *= kernel_peskin(dist);
                         
                        test = forcesI(2)*D/(dx*dy*dz);
                        Fz_net(i_it,j_it,k_it) += test;

                        if (p->mpirank == 1)
                        {
                        //cout<<"in"<<endl;
                        //    cout<<i_it<<" "<<j_it<<" "<<k_it<<" "<<test<<" "<<Fx_net(i_it,j_it,k_it)<<" "<<Fy_net(i_it,j_it,k_it)<<" "<<Fz_net(i_it,j_it,k_it)<<endl;
                        //cout<<"out"<<endl;
                        
                        }
                    }
                }
            }
        }
    } 
    
    pgc->start1(p,Fx_net,10); 
    pgc->start2(p,Fy_net,11); 
    pgc->start3(p,Fz_net,12); 


    //distributeNetForces_x(p, a, pgc, ppnet, nNet);
    //distributeNetForces_y(p, a, pgc, ppnet, nNet);
    //distributeNetForces_z(p, a, pgc, ppnet, nNet);  

    // Distribute collar forces 
    // if (p->X10 == 10) distributeCollarForces(p, a, pgc, ppnet, nNet);
}


double vrans_net::kernel_peskin(const double& dist)
{
    double D = 0.0;

    if (dist <= 2.0)
    {
        D = 0.25*(1.0 + cos(PI*dist/2.0));
    }
    
    return D;
}


void vrans_net::distributeNetForces_x(lexer *p, fdm *a, ghostcell *pgc, net *&ppnet, int nNet)
{ 
    int ii, jj, kk;
    double dist, D, dx, dy, dz;

    const EigenMat& lagrangePoints = ppnet->getLagrangePoints();
    const EigenMat& lagrangeForces = ppnet->getLagrangeForces();
    
    ULOOP
    {
        kernel_x(i,j,k) = 0.0;
    }
    
    if (nNet == 0)
    {
        ULOOP
        {
            Fx_net(i,j,k) = 0.0;
        }        
    }


    for (int pI = 0; pI < lagrangePoints.size(); pI++)
    {
        const Eigen::Vector3d& forcesI = lagrangeForces[pI];
        const Eigen::Vector3d& coordI = lagrangePoints[pI];
         
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
                    dist = sqrt(pow(p->XN[i_it + 1 + marge] - coordI(0), 2.0))/dx;
                    D = kernel_peskin(dist);
                    
                    dist = sqrt(pow(p->YP[j_it + marge] - coordI(1), 2.0))/dy;
                    D *= kernel_peskin(dist);

                    dist = sqrt(pow(p->ZP[k_it + marge] - coordI(2), 2.0))/dz;
                    D *= kernel_peskin(dist);
                    
                   Fx_net(i_it,j_it,k_it) += forcesI(0)*D/(dx*dy*dz);
                }
            }
        }
    } 
    
    pgc->start1(p,Fx_net,10); 
}


void vrans_net::distributeNetForces_y(lexer *p, fdm *a, ghostcell *pgc, net *&ppnet, int nNet)
{ 
    int ii, jj, kk;
    double dist, D, dx, dy, dz;

    const EigenMat& lagrangePoints = ppnet->getLagrangePoints();
    const EigenMat& lagrangeForces = ppnet->getLagrangeForces();
    
    VLOOP
    {
        kernel_y(i,j,k) = 0.0;
    }
    
    if (nNet == 0)
    {
        VLOOP
        {
            Fy_net(i,j,k) = 0.0;
        }        
    }
    
    for (int pI = 0; pI < lagrangePoints.size(); pI++)
    {
        const Eigen::Vector3d& forcesI = lagrangeForces[pI];
        const Eigen::Vector3d& coordI = lagrangePoints[pI];
         
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
                    dist = sqrt(pow(p->XP[i_it + marge] - coordI(0), 2.0))/dx;
                    D = kernel_peskin(dist);
                            
                    dist = sqrt(pow(p->YN[j_it + 1 + marge] - coordI(1), 2.0))/dy;
                    D *= kernel_peskin(dist);

                    dist = sqrt(pow(p->ZP[k_it + marge] - coordI(2), 2.0))/dz;
                    D *= kernel_peskin(dist);
                        
                    Fy_net(i_it,j_it,k_it) += forcesI(1)*D/(dx*dy*dz);
                }
            }
        }
    } 

    pgc->start2(p,Fy_net,11); 
}


void vrans_net::distributeNetForces_z(lexer *p, fdm *a, ghostcell *pgc, net *&ppnet, int nNet)
{ 
    int ii, jj, kk;
    double dist, D, dx, dy, dz;

    const EigenMat& lagrangePoints = ppnet->getLagrangePoints();
    const EigenMat& lagrangeForces = ppnet->getLagrangeForces();
    
    WLOOP
    {
        kernel_z(i,j,k) = 0.0;
    }
    
    if (nNet == 0)
    {
        WLOOP
        {
            Fz_net(i,j,k) = 0.0;
        }        
    }
    
    for (int pI = 0; pI < lagrangePoints.size(); pI++)
    {
        const Eigen::Vector3d& forcesI = lagrangeForces[pI];
        const Eigen::Vector3d& coordI = lagrangePoints[pI];
         
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
                    dist = sqrt(pow(p->XP[i_it + marge] - coordI(0), 2.0))/dx;
                    D = kernel_peskin(dist);
                    
                    dist = sqrt(pow(p->YP[j_it + marge] - coordI(1), 2.0))/dy;
                    D *= kernel_peskin(dist);

                    dist = sqrt(pow(p->ZN[k_it + 1 + marge] - coordI(2), 2.0))/dz;
                    D *= kernel_peskin(dist);
                     
                    Fz_net(i_it,j_it,k_it) += forcesI(2)*D/(dx*dy*dz);
                }
            }
        }
    } 

    pgc->start3(p,Fz_net,12); 
}


void vrans_net::distributeCollarForces(lexer *p, fdm *a, ghostcell *pgc, net *&ppnet, int nNet)
{ 
    int ii, jj, kk;
    double dist, D, dx, dy, dz;
    
    const EigenMat& collarPoints = ppnet->getCollarPoints();
    Eigen::MatrixXd U_l = Eigen::MatrixXd::Zero(collarPoints.size(),3);
    
    for (int pI = 0; pI < collarPoints.size(); pI++)
    {
        const Eigen::Vector3d& coordI = collarPoints[pI];
         
        ii = p->posc_i(coordI(0));
        jj = p->posc_j(coordI(1));
        kk = p->posc_k(coordI(2));
        
        dx = p->DXN[ii + marge];
        dy = p->DYN[jj + marge];
        dz = p->DZN[kk + marge];

        // x-direction
        for (int i_it = ii - 2; i_it <= ii + 2; i_it++)
        {
            for (int j_it = jj - 2; j_it <= jj + 2; j_it++)
            {
                for (int k_it = kk - 2; k_it <= kk + 2; k_it++)
                {
                    dist = sqrt(pow(p->XN[i_it + 1 + marge] - coordI(0), 2.0))/dx;
                    D = kernel_peskin(dist);
                    
                    dist = sqrt(pow(p->YP[j_it + marge] - coordI(1), 2.0))/dy;
                    D *= kernel_peskin(dist);

                    dist = sqrt(pow(p->ZP[k_it + marge] - coordI(2), 2.0))/dz;
                    D *= kernel_peskin(dist);
                    
                    U_l(pI,0) += a->u(i_it, j_it, k_it)*D;
                }
            }
        }
        
        // y-direction
        for (int i_it = ii - 2; i_it <= ii + 2; i_it++)
        {
            for (int j_it = jj - 2; j_it <= jj + 2; j_it++)
            {
                for (int k_it = kk - 2; k_it <= kk + 2; k_it++)
                {
                    dist = sqrt(pow(p->XP[i_it + marge] - coordI(0), 2.0))/dx;
                    D = kernel_peskin(dist);
                            
                    dist = sqrt(pow(p->YN[j_it + 1 + marge] - coordI(1), 2.0))/dy;
                    D *= kernel_peskin(dist);

                    dist = sqrt(pow(p->ZP[k_it + marge] - coordI(2), 2.0))/dz;
                    D *= kernel_peskin(dist);
                        
                    U_l(pI,1) += a->v(i_it, j_it, k_it)*D;
                }
            }
        }
        
        // z-direction
        for (int i_it = ii - 2; i_it <= ii + 2; i_it++)
        {
            for (int j_it = jj - 2; j_it <= jj + 2; j_it++)
            {
                for (int k_it = kk - 2; k_it <= kk + 2; k_it++)
                {
                    dist = sqrt(pow(p->XP[i_it + marge] - coordI(0), 2.0))/dx;
                    D = kernel_peskin(dist);
                    
                    dist = sqrt(pow(p->YP[j_it + marge] - coordI(1), 2.0))/dy;
                    D *= kernel_peskin(dist);

                    dist = sqrt(pow(p->ZN[k_it + 1 + marge] - coordI(2), 2.0))/dz;
                    D *= kernel_peskin(dist);
                     
                    U_l(pI,2) += a->w(i_it, j_it, k_it)*D;
                }
            }
        }
    }

    const EigenMat& collarVel = ppnet->getCollarVel();
    Eigen::MatrixXd forcing_l = Eigen::MatrixXd::Zero(collarPoints.size(),3); 
  
    for (int pI = 0; pI < collarPoints.size(); pI++)
    {
        forcing_l(pI,0) = (collarVel[pI](0) - U_l(pI,0))/p->dt;
        forcing_l(pI,1) = (collarVel[pI](1) - U_l(pI,1))/p->dt;
        forcing_l(pI,2) = (collarVel[pI](2) - U_l(pI,2))/p->dt;
    }

    
    double delta_net = p->X322_D[nNet]*sin(PI/p->X321_nd[nNet]);
    double dV;

    for (int pI = 0; pI < collarPoints.size(); pI++)
    {
        const Eigen::Vector3d& coordI = collarPoints[pI];
        const Eigen::Vector3d& forcesI = forcing_l.row(pI);
         
        ii = p->posc_i(coordI(0));
        jj = p->posc_j(coordI(1));
        kk = p->posc_k(coordI(2));
        
        dx = p->DXN[ii + marge];
        dy = p->DYN[jj + marge];
        dz = p->DZN[kk + marge];

        dV = delta_net/pow(dx*dy*dz, 1.0/3.0);

        // x-direction
        for (int i_it = ii - 2; i_it <= ii + 2; i_it++)
        {
            for (int j_it = jj - 2; j_it <= jj + 2; j_it++)
            {
                for (int k_it = kk - 2; k_it <= kk + 2; k_it++)
                {
                    dist = sqrt(pow(p->XN[i_it + 1 + marge] - coordI(0), 2.0))/dx;
                    D = kernel_peskin(dist);
                    
                    dist = sqrt(pow(p->YP[j_it + marge] - coordI(1), 2.0))/dy;
                    D *= kernel_peskin(dist);

                    dist = sqrt(pow(p->ZP[k_it + marge] - coordI(2), 2.0))/dz;
                    D *= kernel_peskin(dist);
                    
                    Fx_net(i_it,j_it,k_it) -= forcesI(0)*D*dV;
                }
            }
        }
        
        // y-direction
        for (int i_it = ii - 2; i_it <= ii + 2; i_it++)
        {
            for (int j_it = jj - 2; j_it <= jj + 2; j_it++)
            {
                for (int k_it = kk - 2; k_it <= kk + 2; k_it++)
                {
                    dist = sqrt(pow(p->XP[i_it + marge] - coordI(0), 2.0))/dx;
                    D = kernel_peskin(dist);
                            
                    dist = sqrt(pow(p->YN[j_it + 1 + marge] - coordI(1), 2.0))/dy;
                    D *= kernel_peskin(dist);

                    dist = sqrt(pow(p->ZP[k_it + marge] - coordI(2), 2.0))/dz;
                    D *= kernel_peskin(dist);
                        
                    Fy_net(i_it,j_it,k_it) -= forcesI(1)*D*dV;
                }
            }
        }
        
        // z-direction
        for (int i_it = ii - 2; i_it <= ii + 2; i_it++)
        {
            for (int j_it = jj - 2; j_it <= jj + 2; j_it++)
            {
                for (int k_it = kk - 2; k_it <= kk + 2; k_it++)
                {
                    dist = sqrt(pow(p->XP[i_it + marge] - coordI(0), 2.0))/dx;
                    D = kernel_peskin(dist);
                    
                    dist = sqrt(pow(p->YP[j_it + marge] - coordI(1), 2.0))/dy;
                    D *= kernel_peskin(dist);

                    dist = sqrt(pow(p->ZN[k_it + 1 + marge] - coordI(2), 2.0))/dz;
                    D *= kernel_peskin(dist);
                     
                    Fz_net(i_it,j_it,k_it) -= forcesI(2)*D*dV;
                }
            }
        }
    }
    
    pgc->start1(p, Fx_net, 10); 
    pgc->start2(p, Fy_net, 11); 
    pgc->start3(p, Fz_net, 12); 
}

void vrans_net::veltimesave(lexer *p, fdm *a, ghostcell *pgc)
{    
}

