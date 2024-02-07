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
#include"net_barQuasiStatic.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void net_barQuasiStatic::wall_ini(lexer *p, fdm *a, ghostcell *pgc)
{
    D = p->X322_D[nNet];        // Width of wall
    L = p->X322_L[nNet];        // Height of wall
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

/*---------------------------------------*/


    // Initialise values
    ad = D/nd;                          // Length of bars along width
    al = L/nl;                          // Length of bars along height
    
    nf = (2*nd+1) * nl;
    niK = (nd+1)*nl;                    // Number of inner knots
    nbK = nd;                           // Number of boundary bars
    nK = niK + nbK + 1;                 // Total number of knots

    
    
    // Initialise fields
    p->Darray(coupledField, nK, 4);		// fluid coupling matrix (velocity 1,2,3 + phi 4)
    p->Darray(v_t, nf, 3);		// tangential direction
    p->Darray(v_n, nf, 3);		// normal direction
    
    p->Darray(l0, nf);			// initial bar length
    p->Darray(l, nf);			// bar length
    
    fi = MatrixXd::Zero(nf,3);  // inner bar matrix      
    p->Darray(fb,nd,3);         // boundary bar matrix
    
    p->Iarray(Pb,nd);           // boundary bar owner knots
    p->Iarray(Nb,nd);           // boundary bar neighbour knots
    p->Iarray(Pi,nf);           // inner bar owner knots
    p->Iarray(Ni,nf);           // inner bar neighbour knots
    p->Darray(K,nK,3);          // knot coordinates
    p->Darray(K_, nK, 3);       // deformed knot coordinates
      
    A = MatrixXd::Zero(nf,nf);  // Tension force matrix
    B = MatrixXd::Zero(nf,3);   // External force matrix
    Bh = MatrixXd::Zero(nf,3);  // External force matrix with hydrodynamics
    
    p->Iarray(nfK, niK, 5);                // Bars per inner Knot -> first entry is knot ID, then up to four bars (no bar = -1)
    meshID.resize(nd*nl,vector<int>(4));   // List of knots in each mesh


    // Get standard net coordinates in K
    int index = 0;
    for (int i = 0; i <= nl; i++)
    {
        for (int j = 0; j <= nd; j++)
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
        for (int i = 0; i <= nd; i++)
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

    // List of bars per knot nfK
    double tmp;
    bool bK;
    int indexK = 0;
    
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
        for (int k = 0; k < nd; k++)
        {
            if (i == Pb[k] || i == Nb[k])
            {
                bK = true;
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
            if ((K[i][1] == nl) && (K[i][0] != 0.0) && (K[i][0] != nd))
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

    for (int i = 0; i <= nd; i++)
    {
        coord_owner[0] = 0.0;
        coord_owner[1] = ad*i;
      
        for (int j = 0; j <= nl; j++)
        { 
            coord_owner[2] = (nl-j)*al;
        
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

    // Shift and rotate K
    for (int j = 0; j < nK; j++)
    {
        double a = K[j][0];
        double b = K[j][1];
        double c = K[j][2];

        K[j][0] = a*(cos(psi)*cos(theta)) + b*(cos(theta)*sin(psi)) - c*sin(theta);
        K[j][1] = a*(cos(psi)*sin(phi)*sin(theta)-cos(phi)*sin(psi)) + b*(cos(phi)*cos(psi)+sin(phi)*sin(psi)*sin(theta)) + c*(cos(theta)*sin(phi));
        K[j][2] = a*(sin(phi)*sin(psi)+cos(phi)*cos(psi)*sin(theta)) + b*(cos(phi)*sin(psi)*sin(theta)-cos(psi)*sin(phi)) + c*(cos(phi)*cos(theta));    

        K[j][0] += origin_x;
        K[j][1] += origin_y - D/2.0;
        K[j][2] += origin_z;		
    } 
 
 
 
 
    // Initialise fi
    double mag;
    
    for (int i = 0; i < nf; i++)
    {
        fi(i,0) = K[Ni[i]][0] - K[Pi[i]][0];
        fi(i,1) = K[Ni[i]][1] - K[Pi[i]][1];
        fi(i,2) = K[Ni[i]][2] - K[Pi[i]][2];
        
        mag = fi.row(i).norm();
       
        fi(i,0) /= mag;
        fi(i,1) /= mag;
        fi(i,2) /= mag;
        
        if (fabs(fi(i,2)) > 1e-5) // vertical bar
        {
            l0[i] = al;
            l[i] = al;
        }
        else // horizontal bar
        {
            l0[i] = ad;
            l[i] = ad;
        }
    }
    
    // Initialise fb
    for (int i = 0; i < nd; i++)
    {
        fb[i][0] = K[Nb[i]][0] - K[Pb[i]][0];
        fb[i][1] = K[Nb[i]][1] - K[Pb[i]][1];
        fb[i][2] = K[Nb[i]][2] - K[Pb[i]][2];
    }   
 

    // Initialise system of equations

    // Force equilibrium at each inner knot
    
    Fg = 9.81*(rho_c - p->W1)/rho_c*d_c*p->X321_Sn[nNet];  // Gravity force per unit area  
    
    double As;
    double F_ini = 4.0*l_c*Fg;
    int nBars = 0;
 
    index = 0;
            
    for (int i = 0; i < niK; i++)
    {
        nBars = 0;
        
        for (int j = 0; j < nf; j++)
        {
            if (Pi[j] == nfK[i][0])
            {
                A(index,j) = F_ini;
                nBars++;
            }
            else if (Ni[j] == nfK[i][0])
            {
                A(index,j) = -F_ini;
                nBars++;
            }
        }
        
        // Calculate total length of bars adjoint to knot nfK[i]
        if (nBars == 2) 
        {
            As = 0.25*l_c*l_c;
        }
        else if (nBars == 3)
        {
            As = 0.5*l_c*l_c;
        }  
        else
        {
            As = l_c*l_c;
        }

        B(index,2) = Fg*As;
        
        index++;
    }
 
   
    // Sinker 
    int index_ = index - 1;
    
    
    for (int i = 0; i < nd+1; i++)
    {
        B(index_,2) += 0.0645*9.81/(nd+1);
        
        index_--;
    }

    
    
    
    
  
    // Geometrical constraints for boundary meshes, positiv clockwise
    int addOwner, addNeigh;
    
    for (int i = 0; i < nd; i++)
    {
        // Vertical bars
        for (int j = 0; j < nf; j++)
        {
            if (Pi[j] == Pb[i])
            {
                A(index,j) = -al;
                addOwner = Ni[j];
            }
            else if (Ni[j] == Pb[i])
            {
                A(index,j) = al;
                addOwner = Pi[j];               
            }
            else if (Pi[j] == Nb[i])
            {
                A(index,j) = al;
                addNeigh = Ni[j];               
            }
            else if (Ni[j] == Nb[i])
            {
                A(index,j) = -al;
                addNeigh = Pi[j];              
            }
        }  

        // Horizontal bar from addOwner, addNeigh
        for (int j = 0; j < nf; j++)
        {
            if (Pi[j] == addOwner && Ni[j] == addNeigh)
            {
                A(index,j) = -ad;
            }
            else if (Pi[j] == addNeigh && Ni[j] == addOwner)
            {
                A(index,j) = al;              
            }
        }

        B(index,0) = fb[i][0];
        B(index,1) = fb[i][1];
        B(index,2) = fb[i][2];
        
        index++;
    }

    // Geometrical constraints for inner meshes, positiv clockwise (first nd meshes are boundary meshes)
    for (int i = nd; i < meshID.size(); i++)
    {
        for (int j = 0; j < nf; j++)
        {
            // vertical bars
            if (Pi[j] == meshID[i][0] && Ni[j] == meshID[i][1])
            {
                A(index,j) = -al;
            }
            else if (Pi[j] == meshID[i][2] && Ni[j] == meshID[i][3])
            {
                A(index,j) = al;
            }  
          
            // horizontal bars, bars always positiv to the right
            if (Pi[j] == meshID[i][0] && Ni[j] == meshID[i][2])
            {
                A(index,j) = ad;
            }
            else if (Pi[j] == meshID[i][1] && Ni[j] == meshID[i][3])
            {
                A(index,j) = -ad;
            }    
        }
        
        // B[index][] = [0 0 0];
        
        index++;
    }


    // Initialise Bh
    Bh = B;


    // Initialise print
    if(p->mpirank==0 && p->P14==1)
    {
        char str[1000];
        sprintf(str,"./REEF3D_CFD_6DOF/REEF3D_6DOF_Net_Forces_%i.dat",nNet);
        eTout.open(str);
        eTout<<"time [s] \t Tmax [N] \t Fx [N] \t Fy [N] \t Fz [N]"<<endl;
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
