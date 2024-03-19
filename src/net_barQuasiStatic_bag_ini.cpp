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

void net_barQuasiStatic::bag_ini(lexer *p, fdm *a, ghostcell *pgc)
{
    // Input data
 
	// Net wall in current
    L = 1.0;        // length of net in y-direction
    D = 1.0;        // length of net in x-direction b

    nd = 8;          // number of meshes in x-direction  n
    nl = 8;          // number of meshes in y-direction  m 

    EA = 1e7;       // Elasticity
    al = p->X321_lambda[nNet];         // Length of twine
    d_c = p->X321_d[nNet];       // Diameter of twine
    rho_c = 1300;    // density material
    w = 0.08;       // weight per meter in air

    origin_x = 1.5;
    origin_y = 0.5;
    origin_z = 0.2+1.0;
    phi = 0.0;
    theta = 0*PI/180;
    psi = 0.0;
    
/*
	// Net wall in waves
    L = 1.0;        // length of net in y-direction
    D = 1.0;        // length of net in x-direction 

    nd = 8;          // number of meshes in x-direction  
    nl = 8;          // number of meshes in y-direction   

    al = 0.13;     	// length of mesh in [m]
    EA = 1e9;       // Elasticity
    d_c = 0.002;    // diameter
    rho_c = 1100;    // density material
    w = 0.007;       // weight per meter in air

    origin_x = 1.0;
    origin_y = 0.75;
    origin_z = 1.2;
    phi = 0.0;
    theta = 0*PI/180;
    psi = 0.0;   
 */   

 
    if (2*nd*al < (D / sin(atan(D/L))))
    {
        cout<<"wrong input data"<<endl;
    }


    // Initialise values
    
    beta = D/nd;                    // length of boundary bars in x-direction 
    gamma = L/nl;                   // length of boundary bars in y-direction 
    nf = 4*nd*nl;                   // number of bars
    niK = nl*nd+(nl-1)*(nd-1);      // number of inner knots
    nbK = 2*(nd+1)+2*(nl-1);        // number of boundary knots = number of boundary bars
    nK = niK + nbK;                 // number of knots

    Fg = 9.81*(rho_c - p->W1)/rho_c*w;

    
    // Initialise fields
    p->Darray(v_t, nf, 3);		// tangential direction
    p->Darray(v_n, nf, 3);		// normal direction
    
    p->Darray(l0, nf);			// initial bar length
    p->Darray(l, nf);			// bar length
    
    fi = MatrixXd::Zero(nf,3);  // inner bar matrix      
    p->Darray(fb,nbK,3);         // boundary bar matrix
    
    p->Iarray(Pb,nbK);           // boundary owner knots
    p->Iarray(Nb,nbK);           // boundary neighbour knots
    p->Iarray(Pi,nf);           // inner owner knots
    p->Iarray(Ni,nf);           // inner neighbour knots
    p->Darray(K,nK,3);          // knot coordinates
    p->Darray(K_, nK, 3);

    p->Iarray(nfK, niK, 4);		// Bars per Knot
 

    // Get standard net coordinates in K
    genericNet();

    // Initialise inner owner and neighbour lists
    iniInnerKnots();

    // Initialise boundary owner and neighbour lists
    iniBoundaryKnots();   												

    // Stretch net
    stretch();

    // Initialise fi
    for (int i = 0; i < nf; i++)
    {
        l0[i] = al;
        l[i] = al;
        
        fi(i,0) = K[Ni[i]][0] - K[Pi[i]][0];
        fi(i,1) = K[Ni[i]][1] - K[Pi[i]][1];
        fi(i,2) = K[Ni[i]][2] - K[Pi[i]][2];
        
        double mag = fi.row(i).norm();
  
        fi(i,0) /= mag;
        fi(i,1) /= mag;
        fi(i,2) /= mag;
    }
 
    // Initialise system of equations
    iniLSE(p);

    // Initialise print
    if(p->mpirank==0)
    {
        char str[1000];
        sprintf(str,"./REEF3D_CFD_6DOF_Net/REEF3D_6DOF_net_max_force_%i.dat",nNet);
        eTout.open(str);
        eTout<<"time \t Tmax \t Fx \t Fy \t Fz"<<endl;
    }		
    printtime = 0.0;

    // Initialise communication 
    ini_parallel(p, a, pgc);
}


void net_barQuasiStatic::ini_parallel(lexer *p, fdm *a, ghostcell *pgc)
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


void net_barQuasiStatic::genericNet()
{
    int p = 0;
    int index = 0;
    
    for (int i = 0; i < 2*nl+1; i++)
    {
        if (i%2!=0)   // odd index
        {
            for (int q = 0; q < nd; q++)
            {
                K[index][0] = q;
                K[index][1] = p;
                index++;
            }
        
            p++;             
        }
        else            // even index
        {  
            for (int q = 1; q <= nd+1; q++)
            {
                K[index][0] = q - 1;
                K[index][1] = p;
                index++;
            }
        
            p++;     
        }
    }
}


void net_barQuasiStatic::iniInnerKnots()
{
    double xCoord, yCoord;
    int k, Kij;
    bool newK;
    
    for (int w = 0; w < nf; w++)
    {
        Pi[w] = -1;
        Ni[w] = -1;
    }    

    
    // vertical direction
    k = 0;
    
    for (int i = 0; i < K[nK-1][0]; i++)							
    {
        for (int j = 0; j < nK; j++)
        {
            if (K[j][0]==i)
            {
                newK = true;
                
                for (int w = 0; w < nf; w++)
                {
                    if (Pi[w]==j || K[j][1]==K[nK-1][1])
                    {
                        newK = false;
                    }
                }
            
                if (newK==true)
                {
                    Pi[k] = j;
                     k++;
                }
            }
            
        }
    }
    
   
    // horizontal direction
    k = 0;

    for (int i = 0; i < K[nK-1][0]; i++)
    {
        for (int j = 0; j < nK; j++)
        {
            if (K[j][0]==i)
            {
                newK = true;
                
                for (int w = 0; w < nf; w++)
                {
                    if (Ni[w]==j || K[j][1]==0)
                    {
                        newK = false;
                    }
                }
            
                if (newK==true)
                {
                    Ni[k] = j;
                    k++;
                }
            }
            
        }
        
    }
    
    
    // diagonal direction
    
    for (int i = 1; i <= K[nK-1][0]; i++)
    {
        xCoord = i;
        
        for (int j = 0; j < K[nK-1][1]; j++)
        {
            yCoord = j;

            for (int w = 0; w < nK; w++)
            {
                if ((K[w][0]==i) && (K[w][1]==j))
                {
                    Kij = w;
                }
            }

            if (j%2==0)   // even index
            {
                Pi[k] = Kij;
                
                for (int w = 0; w < nK; w++)
                {
                    if ((K[w][0]==i-1) && (K[w][1]==j+1))
                    {
                        Ni[k] = w;
                    }
                }

                k++;
            }
            else            // odd index
            {
                for (int w = 0; w < nK; w++)
                {
                    if ((K[w][0]==i-1) && (K[w][1]==j))
                    {
                        Pi[k] = w;
                    }

                    if ((K[w][0]==i) && (K[w][1]==j+1))
                    {
                        Ni[k] = w;
                    }
                }
                
                k++;
            }
        }
    }
}


void net_barQuasiStatic::iniBoundaryKnots()
{
   int index;
    
    for (int i = 0; i < nd; i++)
    {
        fb[i][0] = beta;
        fb[i][1] = 0.0;
        fb[i][2] = 0.0;
        
        fb[i + nd + nl][0] = beta;
        fb[i + nd + nl][1] = 0.0;
        fb[i + nd + nl][2] = 0.0;
        
        Pb[i] = i;
        Pb[i + nd + nl] = nK - i - 2;
        
        Nb[i] = i + 1;
        Nb[i + nd + nl] = nK - i - 1;
    }
 
    index = 0;
    
    for (int i = nd ; i < nd + nl; i++)
    {
        fb[i][0] = 0.0;
        fb[i][1] = gamma;
        fb[i][2] = 0.0;
        
        fb[i + nd + nl][0] = 0.0;
        fb[i + nd + nl][1] = gamma;
        fb[i + nd + nl][2] = 0.0;
        
        Pb[i] = nd + index*(2*nd + 1);
        Pb[i + nd + nl] = nK - (nd + 1) - (index + 1)*(2*nd + 1);
     
        Nb[i] = nd + (index + 1)*(2*nd + 1);
        Nb[i + nd + nl] = nK - (nd + 1) - index*(2*nd + 1);

        index++;
    }
}


void net_barQuasiStatic::stretch()
{
    // Stretch according to gamma and beta
    for (int j = 0; j < nK; j++)
    {
        K[j][2] = 0.0;
        
        if (j%2!=0)   // odd index   
        {
            for (int k = 0; k < nK; k++)
            {
                if (K[k][1]==j)
                {
                    K[k][0] = K[k][0]*beta + beta/2.0;
                    K[k][1] *= gamma/2.0;
                }
            }
        }
        else        // even index
        {
            for (int k = 0; k < nK; k++)
            {
                if (K[k][1]==j)
                {
                    K[k][0] *= beta;
                    K[k][1] *= gamma/2.0; 
                }
            }
        }
    }
    
    // Shift and rotate K_
    for (int j = 0; j < nK; j++)
    {
        double a = K[j][0];
        double b = K[j][1];
        double c = K[j][2];

        K[j][0] = a*(cos(psi)*cos(theta)) + b*(cos(theta)*sin(psi)) - c*sin(theta);
        K[j][1] = a*(cos(psi)*sin(phi)*sin(theta)-cos(phi)*sin(psi)) + b*(cos(phi)*cos(psi)+sin(phi)*sin(psi)*sin(theta)) + c*(cos(theta)*sin(phi));
        K[j][2] = a*(sin(phi)*sin(psi)+cos(phi)*cos(psi)*sin(theta)) + b*(cos(phi)*sin(psi)*sin(theta)-cos(psi)*sin(phi)) + c*(cos(phi)*cos(theta));    

        K[j][0] += origin_x;
        K[j][1] += origin_y;
        K[j][2] += origin_z;
    }
}


void net_barQuasiStatic::iniLSE(lexer *p)
{
    int index = 0; 
    int numk;
    double Fini;
    bool bk;
    
 
    // Static force equilibria
    Fini = 1.0;
    
    for (int i = 0; i < nK; i++)
    {
        numk = 0;
        bk = false;
        
        for (int j = 0; j < 2*nd+2*nl; j++)
        {
            if (i==Pb[j] || i==Nb[j])
            {
                bk = true;
                break;
            }
        }
        
        if (bk==false)
        {
            for (int j = 0; j < nf; j++)
            {
                if (Pi[j]==i)
                {
                    A(index,j) = Fini;
        
                    nfK[index][numk] = j;
                    numk++;
                }
                else if (Ni[j]==i)
                {
                    A(index,j) = -Fini;
                    
                    nfK[index][numk] = j;
                    numk++;
                }
            }
            
            B(index,0) = 0.0;
            B(index,1) = 0.0;
            B(index,2) = 2.0*Fg*2*al;
  
//if (fabs(K[i][1]-origin_y) < 0.1 && fabs(K[i][2]-1.66) < 0.1) cout<<i<<endl; B[index][0] -= 0.1;  
      
            index++;
        }
    }


    // Geometrical boundary equations
    double dotProd;
    
    for (int i = 0; i < (2*nd+2*nl); i++)
    {
        for (int j = 0; j < nf; j++)
        {
            dotProd = fi(j,0)*fb[i][0] + fi(j,1)*fb[i][1] + fi(j,2)*fb[i][2];
            
            if (Pi[j]==Pb[i] && dotProd > 0)
            {
                A(index,j) = al;
            }
            else if (Pi[j]==Nb[i] && dotProd < 0)
            {
                A(index,j) = -al;
            }
            else if (Ni[j]==Pb[i] && dotProd < 0)
            {
                A(index,j) = -al;
            }
            else if (Ni[j]==Nb[i] && dotProd > 0)
            {
                A(index,j) = al;
            }
  
        }
        B(index,0) = fb[i][0]*(cos(psi)*cos(theta)) + fb[i][1]*(cos(theta)*sin(psi)) - fb[i][2]*sin(theta);
        B(index,1) = fb[i][0]*(cos(psi)*sin(phi)*sin(theta)-cos(phi)*sin(psi)) + fb[i][1]*(cos(phi)*cos(psi)+sin(phi)*sin(psi)*sin(theta)) + fb[i][2]*(cos(theta)*sin(phi));
        B(index,2) = fb[i][0]*(sin(phi)*sin(psi)+cos(phi)*cos(psi)*sin(theta)) + fb[i][1]*(cos(phi)*sin(psi)*sin(theta)-cos(psi)*sin(phi)) + fb[i][2]*(cos(phi)*cos(theta)); 
        
        Bh.row(index) = B.row(index);
     /*   Bh(index,0) = B(index,0);
        Bh(index,1) = B(index,1);
        Bh(index,2) = B(index,2);*/
        
        index++;
    }


    // Geometrical mesh equations
   
    int *knotID;
    p->Iarray(knotID,5);
    
    int currRow = 0;
    int leftK = nd + 1;
  
    while (index < nf)
    {
        if (currRow%2==0)   // even index 
        {
            for (int i = 0; i < nd - 1; i++)
            {
                if (index < nf)
                {
                    knotID[0] = leftK; knotID[1] = leftK + nd + 1;
                    knotID[2] = leftK + 1; knotID[3] = leftK - nd;
                     knotID[4] = leftK;
              
                    for (int j = 0; j < nf; j++)
                    {
                        for (int k = 0; k < 4; k++)
                        {
                            if (Pi[j]==knotID[k] && Ni[j]==knotID[k + 1])
                            {
                                A(index,j) = al;
                            }
                            else if (Ni[j]==knotID[k] && Pi[j]==knotID[k + 1])
                            {
                                A(index,j) = -al;
                            }
                        }
                    }
                    leftK++;
                    
                    index++;
                }
            }
            leftK++;
            currRow++;
        }
        else        // odd index
        {
            for (int i = 0; i < nd; i++)
            {
                if (index < nf)
                {
                    knotID[0] = leftK; knotID[1] = leftK + nd + 1;
                    knotID[2] = leftK + 1; knotID[3] = leftK - nd;
                    knotID[4] = leftK;
                    
                    for (int j = 0; j < nf; j++)
                    {
                        for (int k = 0; k < 4; k++)
                        {
                            if (Pi[j]==knotID[k] && Ni[j]==knotID[k + 1])
                            {
                                A(index,j) = al;
                            }
                            else if (Ni[j]==knotID[k] && Pi[j]==knotID[k + 1])
                            {
                                A(index,j) = -al;
                            }
                        }
                    }
                    leftK++;
                    
                    index++;
                }
            }
            leftK++;
            currRow++;
        }
    }

    p->del_Iarray(knotID,5);
}

