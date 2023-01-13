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
--------------------------------------------------------------------*/

#include"net_barDyn.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"reinidisc_fsf.h"	
#include"vrans.h"


void net_barDyn::fillNonLinSystem(lexer *p, fdm *a, ghostcell *pgc)
{
    bool foundOwner, foundNeigh;
    int ownerKnot, neighKnot, otherKnot;
    Vector3d X_ij, x_ij, x_ijn, x_ijnn, x_jk, x_ik, b_ij, b_ji;
    Vector3d V_ij, v_ij, v_ijn, v_ijnn;
    Vector3d contrib, sum_neigh, sum_own, bracket, A_ij, F_i, F_j;
    double l, l0_ij, mass_i, mass_j;
    Eigen::ArrayXi ownerBars(5);
    Eigen::ArrayXi neighBars(5);

    A_ *= 0.0;
    
    for (int barI = 0; barI < nf; barI++)
    {
        ownerKnot = Pi[barI];
        neighKnot = Ni[barI];        
        
        // Masses lumped at knots
        mass_i = mass_knot(ownerKnot) + added_mass(ownerKnot);
        mass_j = mass_knot(neighKnot) + added_mass(neighKnot);

        // Forces lumped at knots
        F_i = forces_knot.row(ownerKnot);
        F_j = forces_knot.row(neighKnot);
        
        // Bar difference vector
        x_ij = x_.row(neighKnot) - x_.row(ownerKnot);
        x_ijn = xn_.row(neighKnot) - xn_.row(ownerKnot);
        x_ijnn = xnn_.row(neighKnot) - xnn_.row(ownerKnot);

        X_ij = -coeffs_(0)*
            (   
                coeffs_(1)*x_ij + coeffs_(2)*x_ijn + coeffs_(3)*x_ijnn
            );        
        
        // Velocity difference vector
        v_ij = xdot_.row(neighKnot) - xdot_.row(ownerKnot);
        v_ijn = xdotn_.row(neighKnot) - xdotn_.row(ownerKnot);
        v_ijnn = xdotnn_.row(neighKnot) - xdotnn_.row(ownerKnot);

        V_ij = -1.0*
            (   
                coeffs_(1)*v_ij + coeffs_(2)*v_ijn + coeffs_(3)*v_ijnn
            );        
 
        // Acceleration vector contribution A_ij
        A_ij << 0.0, 0.0, 0.0;
        if (mass_i > 0.0) A_ij -= F_i/mass_i;   
        if (mass_j > 0.0) A_ij += F_j/mass_j;
        
 
        // Length of original bar vector
        l0_ij = (x0_.row(neighKnot) - x0_.row(ownerKnot)).norm();

        
        // Find bars attached to the two knots
        foundOwner = false;
        foundNeigh = false;
        
        ownerBars << -1, -1, -1, -1, -1;
        neighBars << -1, -1, -1, -1, -1;
        
        for (int k = 0; k < niK; k++)
        {
            if (nfK[k][0] == ownerKnot)
            {
                ownerBars << nfK[k][0], nfK[k][1], nfK[k][2], nfK[k][3], nfK[k][4];
                foundOwner = true;
            }
            else if (nfK[k][0] == neighKnot)
            {
                neighBars << nfK[k][0], nfK[k][1], nfK[k][2], nfK[k][3], nfK[k][4];
                foundNeigh = true;
            }
        }

        // Add coefficents of neighbour knot bars
        sum_neigh << 0.0, 0.0, 0.0;
        if (foundNeigh == true)    
        {
            for (int j = 1; j < 5; j++)
            {
                int& bar_jk = neighBars(j);
                
                if (bar_jk != -1)
                {
                    // Bar vector pointing from neighKnot to otherKnot
                    if (Pi[bar_jk] == neighKnot)
                    {
                        otherKnot = Ni[bar_jk];
                    }
                    else
                    {
                        otherKnot = Pi[bar_jk];
                    }
                    
                    x_jk = x_.row(otherKnot) - x_.row(neighKnot);

                    // Length of bar vector
                    l = x_jk.norm();

                    // Add coefficient
                    if (l > 0.0 && mass_j > 0.0)
                    {                        
                        sum_neigh += x_jk/l*T_(bar_jk)/mass_j;
                    }
                }
            }  
        }
    
        // Add coefficents of owner knot bars
        sum_own << 0.0, 0.0, 0.0;
        if (foundOwner == true)    
        {
            for (int i = 1; i < 5; i++)
            {
                int& bar_ik = ownerBars(i);
                
                if (bar_ik != -1)
                {
                    // Bar vector pointing from neighKnot to otherKnot
                    if (Pi[bar_ik] == ownerKnot)
                    {
                        otherKnot = Ni[bar_ik];
                    }
                    else
                    {
                        otherKnot = Pi[bar_ik];
                    }
                    
                    x_ik = x_.row(otherKnot) - x_.row(ownerKnot);

                    // Length of bar vector
                    l = x_ik.norm();

                    // Add coefficient
                    if (l > 0.0 && mass_i > 0.0)
                    {
                        sum_own += x_ik/l*T_(bar_ik)/mass_i;
                    }
                }
            }  
        } 


        // Bracket term
        bracket = sum_neigh - sum_own + A_ij + V_ij + X_ij; 


        // Main diagonal coefficient
        b_ij = (x_.row(neighKnot) - x_.row(ownerKnot)).normalized();
        b_ji = (x_.row(ownerKnot) - x_.row(neighKnot)).normalized();
        
        if (mass_j > 0.0 && mass_i > 0.0)
        {
            contrib = b_ji/mass_j - b_ij/mass_i;
        }
        else if (mass_i > 0.0)
        {
            contrib = - b_ij/mass_i;
        }  
        else
        {
            contrib = b_ji/mass_j;
        } 
        
        A_(barI,barI) = 
            2.0*bracket.dot(contrib) - pow(coeffs_(0), 4.0)*pow(l0_ij, 2.0)/C2_
            *
            (
                (-C1_ + 2.0*C2_ + sqrt(C1_*C1_ + 4.0*C2_*T_(barI)))
                /
                (sqrt(C1_*C1_ + 4.0*C2_*T_(barI)))
            );
    
 
        // Owner bar coefficients
        if (foundOwner == true)    
        {
            for (int i = 1; i < 5; i++)
            {
                int& bar_ik = ownerBars(i);
                
                if (bar_ik != -1 && bar_ik != barI)
                {
                    // Bar vector pointing from neighKnot to otherKnot
                    if (Pi[bar_ik] == ownerKnot)
                    {
                        otherKnot = Ni[bar_ik];
                    }
                    else
                    {
                        otherKnot = Pi[bar_ik];
                    }
                    
                    x_ik = (x_.row(otherKnot) - x_.row(ownerKnot)).normalized();

                    // Add coefficient
                    if (mass_i > 0.0)
                    { 
                        A_(barI,bar_ik) = -2.0/mass_i*x_ik.dot(bracket);
                    }
                }
            }  
        } 
        
        // Neighbour bar coefficients
        if (foundNeigh == true)    
        {
            for (int j = 1; j < 5; j++)
            {
                int& bar_jk = neighBars(j);
                
                if (bar_jk != -1 && bar_jk != barI)
                {
                    // Bar vector pointing from neighKnot to otherKnot
                    if (Pi[bar_jk] == neighKnot)
                    {
                        otherKnot = Ni[bar_jk];
                    }
                    else
                    {
                        otherKnot = Pi[bar_jk];
                    }
                    
                    x_jk = (x_.row(otherKnot) - x_.row(neighKnot)).normalized();

                    // Add coefficient
                    if (mass_j > 0.0)
                    {                        
                        A_(barI,bar_jk) = 2.0/mass_j*x_jk.dot(bracket);
                    }
                }
            }  
        }        
  
    }
}


void net_barDyn::fillNonLinRhs(lexer *p, fdm *a, ghostcell *pgc)
{
    bool foundOwner, foundNeigh;
    int ownerKnot, neighKnot, otherKnot;
    Vector3d X_ij, x_ij, x_ijn, x_ijnn, x_jk, x_ik;
    Vector3d V_ij, v_ij, v_ijn, v_ijnn;
    Vector3d FbarI, sum_neigh, sum_own, A_ij, F_i, F_j;
    double l, l0_ij, mass_i, mass_j;
    Eigen::ArrayXi ownerBars(5);
    Eigen::ArrayXi neighBars(5);

    B_ *= 0.0;

    for (int barI = 0; barI < nf; barI++)
    {
        ownerKnot = Pi[barI];
        neighKnot = Ni[barI];         
        
        // Masses lumped at knots
        mass_i = mass_knot(ownerKnot) + added_mass(ownerKnot);
        mass_j = mass_knot(neighKnot) + added_mass(neighKnot);     
        
         // Forces lumped at knots
        F_i = forces_knot.row(ownerKnot);
        F_j = forces_knot.row(neighKnot);         

        // Bar difference vector
        x_ij = x_.row(neighKnot) - x_.row(ownerKnot);
        x_ijn = xn_.row(neighKnot) - xn_.row(ownerKnot);
        x_ijnn = xnn_.row(neighKnot) - xnn_.row(ownerKnot);

        X_ij = -coeffs_(0)*
            (   
                coeffs_(1)*x_ij + coeffs_(2)*x_ijn + coeffs_(3)*x_ijnn
            );        
        
        // Velocity difference vector
        v_ij = xdot_.row(neighKnot) - xdot_.row(ownerKnot);
        v_ijn = xdotn_.row(neighKnot) - xdotn_.row(ownerKnot);
        v_ijnn = xdotnn_.row(neighKnot) - xdotnn_.row(ownerKnot);

        V_ij = -1.0*
            (   
                coeffs_(1)*v_ij + coeffs_(2)*v_ijn + coeffs_(3)*v_ijnn
            );        
 
        // Acceleration vector contribution A_ij
        A_ij << 0.0, 0.0, 0.0;
        if (mass_i > 0.0) 
        {
            A_ij -= F_i/mass_i;   
        }   
        else
        {
            A_ij -= xdotdot_.row(ownerKnot);
        }

        if (mass_j > 0.0) 
        {
            A_ij += F_j/mass_j;
        }
        else
        {
            A_ij += xdotdot_.row(neighKnot); 
        }
        
        // Length of original bar vector
        l0_ij = (x0_.row(neighKnot) - x0_.row(ownerKnot)).norm();     
 
 
        // Find bars attached to the two knots
        foundOwner = false;
        foundNeigh = false;
        
        ownerBars << -1, -1, -1, -1, -1;
        neighBars << -1, -1, -1, -1, -1;
        
        for (int k = 0; k < niK; k++)
        {
            if (nfK[k][0] == ownerKnot)
            {
                ownerBars << nfK[k][0], nfK[k][1], nfK[k][2], nfK[k][3], nfK[k][4];
                foundOwner = true;
            }
            else if (nfK[k][0] == neighKnot)
            {
                neighBars << nfK[k][0], nfK[k][1], nfK[k][2], nfK[k][3], nfK[k][4];
                foundNeigh = true;
            }
        } 
  
         // Add coefficents of neighbour knot bars
        sum_neigh << 0.0, 0.0, 0.0;
        if (foundNeigh == true)    
        {
            for (int j = 1; j < 5; j++)
            {
                int& bar_jk = neighBars(j);
                
                if (bar_jk != -1)
                {
                    // Bar vector pointing from neighKnot to otherKnot
                    if (Pi[bar_jk] == neighKnot)
                    {
                        otherKnot = Ni[bar_jk];
                    }
                    else
                    {
                        otherKnot = Pi[bar_jk];
                    }
                    
                    x_jk = x_.row(otherKnot) - x_.row(neighKnot);

                    // Length of bar vector
                    l = x_jk.norm();

                    // Add coefficient
                    if (l > 0.0 && mass_j > 0.0)
                    {                        
                        sum_neigh += x_jk/l*T_(bar_jk)/mass_j;
                    }
                }
            }  
        }
    
        // Add coefficents of owner knot bars
        sum_own << 0.0, 0.0, 0.0;
        if (foundOwner == true)    
        {
            for (int i = 1; i < 5; i++)
            {
                int& bar_ik = ownerBars(i);
                
                if (bar_ik != -1)
                {
                    // Bar vector pointing from neighKnot to otherKnot
                    if (Pi[bar_ik] == ownerKnot)
                    {
                        otherKnot = Ni[bar_ik];
                    }
                    else
                    {
                        otherKnot = Pi[bar_ik];
                    }
                    
                    x_ik = x_.row(otherKnot) - x_.row(ownerKnot);

                    // Length of bar vector
                    l = x_ik.norm();

                    // Add coefficient
                    if (l > 0.0 && mass_i > 0.0)
                    {
                        sum_own += x_ik/l*T_(bar_ik)/mass_i;
                    }
                }
            }  
        } 
  
        // Sum up terms for B_(barI)
        FbarI = sum_neigh - sum_own + A_ij + V_ij + X_ij;
        
        B_(barI) = 
            FbarI.dot(FbarI) 
            - pow(coeffs_(0), 4.0)*l0_ij*l0_ij/(4.0*C2_*C2_)
            *pow((-C1_ + 2.0*C2_ + sqrt(C1_*C1_ + 4.0*C2_*T_(barI))), 2.0);     
    }    
}
