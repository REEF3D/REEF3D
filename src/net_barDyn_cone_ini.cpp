/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2024 Tobias Martin

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

#include<sys/stat.h>
#include"net_barDyn.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void net_barDyn::cone_ini(lexer *p, fdm *a, ghostcell *pgc)
{
    D = p->X322_D[nNet];        // Width of top of cone
    L = p->X322_L[nNet];        // Height of cone
    nd = p->X321_nd[nNet];      // Number of meshes over width  
    nl = p->X321_nl[nNet];      // Number of meshes over height    
    
    l_c = p->X321_lambda[nNet];  // Length of twine
    d_c = p->X321_d[nNet];       // Diameter of twine
    rho_c = p->X321_rho[nNet];   // Density of material
    
    origin_x = p->X322_x0[nNet];  
    origin_y = p->X322_y0[nNet];  
    origin_z = p->X322_z0[nNet];  
    phi = p->X322_phi[nNet]*PI/180;  
    theta = p->X322_theta[nNet]*PI/180;  
    psi = p->X322_psi[nNet]*PI/180;  

    kappa = 0.01;                   // Elasticity constant
    C1_ = 1160.0;    
    C2_ = 37300.0;  
   
    sinker_m = p->X323_m;    // Sinker mass in air [kg]
    sinker_d = p->X323_d;    // Sinker diameter [m]
    sinker_l = p->X323_l;    // Sinker length [m]

    knot_d = p->X321_dk[nNet];  // Knot diameter
    

 /*---------------------------------------*/

    //- Initialise values

    ad = D*sin(PI/nd);              // Length of bars around top cylinder
    al = L/nl;                      // Length of bars along height

    nf = 2*nd*nl; 
    niK = nd*nl;                    // Number of inner knots
    nbK = nd;                       // Number of boundary knots = number of boundary bars    
    nK = niK + nbK;                 // Total number of knots
     
    //- Initialise fields
    
    p->Darray(coupledField, nK, 4);	// fluid coupling matrix (velocity 1,2,3 + phi 4)
    p->Darray(coupledFieldn, nK, 4);
    
    p->Darray(l0, nf);			    // initial bar length

    p->Iarray(Pb,nbK);              // boundary owner knots
    p->Iarray(Nb,nbK);              // boundary neighbour knots
    p->Iarray(Pi,nf);               // inner owner knots
    p->Iarray(Ni,nf);               // inner neighbour knots
    p->Darray(K,nK,3);              // knot coordinates

    p->Iarray(nfK, niK, 5);		          // Bars per Knot -> first entry is knot ID, then up to four bars
    p->Iarray(nfbK, nbK);		          // Bars per Knot -> first entry is knot ID, then up to four bars
    meshID.resize(nd*nl,vector<int>(4));  // List of knots in each mesh
 
    x0_ = MatrixXd::Zero(nK,3); 
    x_ = MatrixXd::Zero(nK,3);
    xn_ = MatrixXd::Zero(nK,3);
    xnn_ = MatrixXd::Zero(nK,3);
    xdot_ = MatrixXd::Zero(nK,3);
    xdotn_ = MatrixXd::Zero(nK,3);
    xdotnn_ = MatrixXd::Zero(nK,3);
    xdotdot_ = MatrixXd::Zero(nK,3);
    
    top_xdot_ = MatrixXd::Zero(nbK,3);
    top_xdotdot_ = MatrixXd::Zero(nbK,3);

    mass_knot = VectorXd::Zero(nK);
    weight_knot = VectorXd::Zero(nK);
    sinker_knot = VectorXd::Zero(nK);
    added_mass = VectorXd::Zero(nK);
    forces_knot = MatrixXd::Zero(nK, 3);

    A_ = MatrixXd::Zero(nf,nf);   // System matrix
    B_ = VectorXd::Zero(nf);      // Rhs
    T_ = VectorXd::Zero(nf);      // Tension forces
    T_old = VectorXd::Zero(nf);   
    T_backup = VectorXd::Zero(nf);   
 
    sinker_m /= nd;    // Sinker mass per bottom knot
    sinker_l /= nd;    // Sinker length per bottom knot






    //- Get standard net coordinates in K
    
    int index = 0;
    for (int i = 0; i <= nl; i++)
    {
        for (int j = 0; j < nd; j++)
        {
            K[index][0] = j;
            K[index][1] = i;
            K[index][2] = 0.0;

            index++;
        }
    }


    // Initialise inner owner and neighbour lists
    vector<double> coord_owner(3,0.0);
    vector<double> coord_neigh(3,0.0);
    
    index = 0;

    // vertical twines
    for (int j = 0; j < nl; j++)
    {
        for (int i = 0; i < nd; i++)
        {
            coord_owner[0] = i;
            coord_owner[1] = j;
            coord_owner[2] = 0.0;
            
            coord_neigh[0] = i;
            coord_neigh[1] = j+1;
            coord_neigh[2] = 0.0;
           
            for (int k = 0; k < nK; k++)
            {
                if (K[k][0] == coord_owner[0] && K[k][1] == coord_owner[1] && K[k][2] == coord_owner[2])
                {
                    Pi[index] = k;
                }
                else if (K[k][0] == coord_neigh[0] && K[k][1] == coord_neigh[1] && K[k][2] == coord_neigh[2])
                {
                    Ni[index] = k;
                }
            }
            
            index++;
        }
    }

    // horizontal twines
    for (int j = 1; j <= nl; j++)
    {
        for (int i = 0; i < nd; i++)
        {
            coord_owner[0] = i;
            coord_owner[1] = j;
            coord_owner[2] = 0.0;
            
            coord_neigh[0] = i+1;
            coord_neigh[1] = j;
            coord_neigh[2] = 0.0;
            
            // Close cylinder
            if (i == nd-1) coord_neigh[0] = 0;              

            for (int k = 0; k < nK; k++)
            {
                if (K[k][0] == coord_owner[0] && K[k][1] == coord_owner[1] && K[k][2] == coord_owner[2])
                {
                    Pi[index] = k;
                }
                else if (K[k][0] == coord_neigh[0] && K[k][1] == coord_neigh[1] && K[k][2] == coord_neigh[2])
                {
                    Ni[index] = k;
                }
            }
            
            index++;
        }
    }    
 
           
    // Initialise boundary owner and neighbour lists - left is always owner
    index = 0;
    
    for (int i = 0; i < nd; i++)
    {
        coord_owner[0] = i;
        coord_owner[1] = 0.0;
        coord_owner[2] = 0.0;
            
        coord_neigh[0] = i+1;
        coord_neigh[1] = 0.0;
        coord_neigh[2] = 0.0;  
    
        // Close cylinder
        if (i == nd-1) coord_neigh[0] = 0.0;      
   
        for (int k = 0; k < nK; k++)
        {
            if (K[k][0] == coord_owner[0] && K[k][1] == coord_owner[1] && K[k][2] == coord_owner[2])
            {
                Pb[index] = k;
            }
            else if (K[k][0] == coord_neigh[0] && K[k][1] == coord_neigh[1] && K[k][2] == coord_neigh[2])
            {
                Nb[index] = k;
            }
        }
            
        index++;
    }
   
 
   // List of knots in each mesh
    index = 0;
    
    vector<double> v1(3,0.0), v2(3,0.0), v3(3,0.0), v4(3,0.0);

    for (int j = 0; j < nl; j++)
    {
        for (int i = 0; i < nd; i++)
        {
            v1[0] = i;
            v1[1] = j;
            v2[0] = i;
            v2[1] = j+1;
            v3[0] = i+1;
            v3[1] = j+1;
            v4[0] = i+1;
            v4[1] = j;
            
            // Close cylinder
            if (i == nd-1)
            {
                v3[0] = 0.0;
                v4[0] = 0.0;
            }
            
            for (int s = 0; s < nf; s++)
            {
                // left half, bar always pointing downwards
                if 
                (
                    v1[0] == K[Pi[s]][0] && v1[1] == K[Pi[s]][1] && v1[2] == K[Pi[s]][2]
                    && v2[0] == K[Ni[s]][0] && v2[1] == K[Ni[s]][1] && v2[2] == K[Ni[s]][2]
                )
                {
                    meshID[index][0] = Pi[s];
                    meshID[index][1] = Ni[s];
                }
                
                // right half, bar always pointing downwards
                if 
                (
                    v4[0] == K[Pi[s]][0] && v4[1] == K[Pi[s]][1] && v4[2] == K[Pi[s]][2]
                    && v3[0] == K[Ni[s]][0] && v3[1] == K[Ni[s]][1] && v3[2] == K[Ni[s]][2]
                )
                {
                    meshID[index][2] = Pi[s];
                    meshID[index][3] = Ni[s];
                }
            }
 
            index++;
        }
    } 



    // List of bars per knot nfK
    double tmp;
    bool bK;
    int indexK = 0;
    int indexbK = 0;
    
    for (int i = 0; i < niK; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            nfK[i][j] = -1;
        }
    }
    
    for (int i = 0; i < nK; i++)
    {
        bK = false;
		
        // Check whether it is a boundary knot
        for (int k = 0; k < nbK; k++)
        {
            if (i == Pb[k] || i == Nb[k])
            {
                bK = true;
                
                // Search for adjoint vertical bar
                for (int j = 0; j < nf; j++)
                {
                    if (Pi[j] == i || Ni[j] == i)
                    {
                        nfbK[indexbK] = j;
                        indexbK++;
                    }
                }
                
                break;
            }
        }

        if (bK == false) // then a inner knot
        {
            nfK[indexK][0] = i;

            // Search adjointed bars
            index = 1;
            for (int j = 0; j < nf; j++)
            {
                if (Pi[j] == i || Ni[j] == i)
                {
                    nfK[indexK][index] = j;
                    index++;
                }
            }
            
            // Switch bottom edge knots to be consistent with side edge knots
            if (K[i][1] == nl)
            {
                tmp = nfK[indexK][3];
                nfK[indexK][3] = nfK[indexK][1];
                nfK[indexK][1] = tmp;
            }         
            
            indexK++;
        }
    }
    
    // Stretch net
    vector<int> transK(nK,0);
    double dz;

    for (int i = 0; i < nd; i++)
    {
        for (int j = 0; j <= nl; j++)
        {
            dz = (nl - j)*al;

            coord_owner[0] = (D/2.0*dz/L + 0.1)*cos(2*PI/nd*i);
            coord_owner[1] = (D/2.0*dz/L + 0.1)*sin(2*PI/nd*i);
            coord_owner[2] = dz;
        
            for (int k = 0; k < nK; k++)
            {
                if (K[k][0] == i && K[k][1] == j && K[k][2] == 0.0 && transK[k] == 0)
                { 
                    K[k][0] = coord_owner[0];
                    K[k][1] = coord_owner[1];
                    K[k][2] = coord_owner[2];
                    
                    transK[k] = 1;
                }
            } 
        }    
    }
   

    // Initialise length of bars
    double mag;
    Vector3d fi;
    
    for (int i = 0; i < nf; i++)
    {
        fi(0) = K[Ni[i]][0] - K[Pi[i]][0];
        fi(1) = K[Ni[i]][1] - K[Pi[i]][1];
        fi(2) = K[Ni[i]][2] - K[Pi[i]][2];
        
        mag = fi.norm();
        
        fi /= (mag + 1e-10);

        l0[i] = mag;
    }
    

    // Shift and rotate K
    for (int j = 0; j < nK; j++)
    {
        double a = K[j][0];
        double b = K[j][1];
        double c = K[j][2];

        x_(j,0) = a*(cos(psi)*cos(theta)) + b*(cos(theta)*sin(psi)) - c*sin(theta);
        x_(j,1) = a*(cos(psi)*sin(phi)*sin(theta)-cos(phi)*sin(psi)) + b*(cos(phi)*cos(psi)+sin(phi)*sin(psi)*sin(theta)) + c*(cos(theta)*sin(phi));
        x_(j,2) = a*(sin(phi)*sin(psi)+cos(phi)*cos(psi)*sin(theta)) + b*(cos(phi)*sin(psi)*sin(theta)-cos(psi)*sin(phi)) + c*(cos(phi)*cos(theta));    

        x_(j,0) += origin_x;
        x_(j,1) += origin_y;
        x_(j,2) += origin_z;
    } 
    
    x0_ = x_;


    // Calculate mesh areas
    vector<double> meshArea;
    meshArea.resize(nd*nl,0.0); 
    Vector3d p1, p2, p3, p4;
    double at,bt,ct,st,x0,x1,x2,y0,y1,y2,z0,z1,z2;

    for (int id = 0; id < meshID.size(); id++)
    {
        p1 << x_(meshID[id][0], 0), x_(meshID[id][0], 1), x_(meshID[id][0], 2);
        p2 << x_(meshID[id][1], 0), x_(meshID[id][1], 1), x_(meshID[id][1], 2);
        p3 << x_(meshID[id][2], 0), x_(meshID[id][2], 1), x_(meshID[id][2], 2);
        p4 << x_(meshID[id][3], 0), x_(meshID[id][3], 1), x_(meshID[id][3], 2);

        // Triangle 1
        x0 = p1(0); x1 = p2(0); x2 = p3(0);
        y0 = p1(1); y1 = p2(1); y2 = p3(1);
        z0 = p1(2); z1 = p2(2); z2 = p3(2);
        at = sqrt(pow(x1 - x0, 2.0) + pow(y1 - y0, 2.0) + pow(z1 - z0, 2.0));
        bt = sqrt(pow(x1 - x2, 2.0) + pow(y1 - y2, 2.0) + pow(z1 - z2, 2.0));
        ct = sqrt(pow(x2 - x0, 2.0) + pow(y2 - y0, 2.0) + pow(z2 - z0, 2.0));
        st = 0.5*(at+bt+ct);
            
        meshArea[id] += sqrt(MAX(0.0,st*(st-at)*(st-bt)*(st-ct)));

        // Triangle 2
        x0 = p4(0); x1 = p2(0); x2 = p3(0);
        y0 = p4(1); y1 = p2(1); y2 = p3(1);
        z0 = p4(2); z1 = p2(2); z2 = p3(2);
        at = sqrt(pow(x1 - x0, 2.0) + pow(y1 - y0, 2.0) + pow(z1 - z0, 2.0));
        bt = sqrt(pow(x1 - x2, 2.0) + pow(y1 - y2, 2.0) + pow(z1 - z2, 2.0));
        ct = sqrt(pow(x2 - x0, 2.0) + pow(y2 - y0, 2.0) + pow(z2 - z0, 2.0));
        st = 0.5*(at+bt+ct);
            
        meshArea[id] += sqrt(MAX(0.0,st*(st-at)*(st-bt)*(st-ct)));
    }

    
    // Mass lumping at each inner knot
    double A_panel, A_solid, l_solid;
    index = 0;
    
    for (int i = nfK[0][0]; i < nK; i++)
    {
        A_panel = 0.0;
        A_solid = 0.0;
        l_solid = 0.0;
        
        // Lookup adjoint meshes and sum up area
        for (int id = 0; id < meshID.size(); id++)
        {
            for (int k = 0; k < 4; k++)
            {
                if (meshID[id][k] == nfK[index][0])
                {
                    A_panel += 0.25*meshArea[id];
                }
            }
        }
        
        A_solid = p->X321_Sn[nNet]*A_panel;
        l_solid = A_solid/d_c;    
        
        // Mass (in air)
        mass_knot(i) = rho_c*PI/4.0*d_c*d_c*l_solid; 
        
        // Weight (in water)
        weight_knot(i) = p->W1*PI/4.0*d_c*d_c*l_solid; 

        // Added mass assuming ca = 1.0
        added_mass(i) = p->W1*PI/4.0*d_c*d_c*1.0*l_solid; 
  
        // Add sinker to bottom row
        if (i >= nK - nd)
        {
            mass_knot(i) += sinker_m;
           
            weight_knot(i) += p->W1*PI/4.0*sinker_d*sinker_d*sinker_l;

            added_mass(i) += p->W1*PI/4.0*sinker_d*sinker_d*sinker_l*1.0;
        }
        
        index++;
    }
  

    // Initialise probe points
    Eigen::Vector3d ppI;
    double dist, dist_new;
    
    probeKnot.resize(p->X324);

    for (int pp = 0; pp < p->X324; pp++)
    {
        ppI << p->X324_x[pp], p->X324_y[pp], p->X324_z[pp];
        
        dist = 1e10;

        for (int kI = 0; kI < nK; kI++)
        {
            dist_new = (x_.row(kI).transpose() - ppI).norm();

            if (dist_new < dist) 
            {
                dist = dist_new;
                probeKnot(pp) = kI;
            }
        }
        
        // Initialise print
        if(p->mpirank==0 && p->P14==1)
        {
            char str[1000];
            sprintf(str,"./REEF3D_CFD_6DOF/REEF3D_6DOF_Net_%i_Point_Probe_%i.dat",nNet,pp+1);
            ofstream header_out;
            header_out.open(str);
            header_out<<"Knot point probe located near "<<ppI.transpose()<<endl;
            header_out<<"time [s] \t x [m] \t y [m] \t z [m]"<<endl;
            header_out.close();
        }		
    }




    // Initialise force print
    Tne = 0.0;
    if(p->mpirank==0 && p->P14==1)
    {
        char str[1000];
        sprintf(str,"./REEF3D_CFD_6DOF/REEF3D_6DOF_Net_Forces_%i.dat",nNet);
        ofstream header_out;
        header_out.open(str);
        header_out<<"time [s] \t Ttop [N] \t Fx [N] \t Fy [N] \t Fz [N]"<<endl;
        header_out.close();
    }		
    printtime = 0.0;

    // Initialise communication 
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
