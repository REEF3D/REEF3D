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

#include"net_barDyn.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"reinidisc_fsf.h"	
#include"vrans.h"


void net_barDyn::fillLinSystem(lexer *p, fdm *a, ghostcell *pgc)
{
    bool foundOwner, foundNeigh;
    int ownerKnot, neighKnot, otherKnot;
    Vector3d x_ij, x_ijn, x_ijnn, x_jk, x_ik, x0_ij;
    double l, mass_i, mass_j;
    
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
        
        // Bar vector
        x_ij = x_.row(neighKnot) - x_.row(ownerKnot);
        x_ijn = xn_.row(neighKnot) - xn_.row(ownerKnot);
        x_ijnn = xnn_.row(neighKnot) - xnn_.row(ownerKnot);

        x_ij = 
            1.0/coeffs_(0)*
            (   
                -coeffs_(1)*x_ij - coeffs_(2)*x_ijn - coeffs_(3)*x_ijnn
            );        
        
        
        // Add main diagonal
        x0_ij = x0_.row(neighKnot) - x0_.row(ownerKnot);
        l = x0_ij.norm();
        
        A_(barI,barI) = -coeffs_(0)*l*l*kappa;     

        
        // Find bars attached to the two knots
        foundOwner = false;
        foundNeigh = false;
        
        ownerBars << -1, -1, -1, -1, -1;
        neighBars << -1, -1, -1, -1, -1;
        
        for (int k = 0; k < niK; k++)
        {
            if (nfK[k][0]==ownerKnot)
            {
                ownerBars << nfK[k][0], nfK[k][1], nfK[k][2], nfK[k][3], nfK[k][4];
                foundOwner = true;
            }
            else if (nfK[k][0]==neighKnot)
            {
                neighBars << nfK[k][0], nfK[k][1], nfK[k][2], nfK[k][3], nfK[k][4];
                foundNeigh = true;
            }
        }

        // Add coefficents of neighbour knot bars
        if (foundNeigh==true)    
        {
            for (int j = 1; j < 5; j++)
            {
                int& bar_jk = neighBars(j);
                
                if (bar_jk!=-1)
                {
                    // Bar vector pointing from neighKnot to otherKnot
                    if (Pi[bar_jk]==neighKnot)
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
                        A_(barI,bar_jk) += x_ij.dot(x_jk/l)/mass_j;
                    }
                }
            }  
        }
    
        // Add coefficents of owner knot bars
        if (foundOwner==true)    
        {
            for (int i = 1; i < 5; i++)
            {
                int& bar_ik = ownerBars(i);
                
                if (bar_ik!=-1)
                {
                    // Bar vector pointing from neighKnot to otherKnot
                    if (Pi[bar_ik]==ownerKnot)
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
                        A_(barI,bar_ik) -= x_ij.dot(x_ik/l)/mass_i;
                    }
                }
            }  
        }     
    }
}


void net_barDyn::fillLinRhs(lexer *p, fdm *a, ghostcell *pgc)
{
    int ownerKnot, neighKnot;
    Vector3d x_ij, x_ijn, x_ijnn, v_ij, v_ijn, v_ijnn, F_i, F_j;
    double l0_ij, mass_i, mass_j;

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

        // Velocity difference vector
        v_ij = xdot_.row(neighKnot) - xdot_.row(ownerKnot);
        v_ijn = xdotn_.row(neighKnot) - xdotn_.row(ownerKnot);
        v_ijnn = xdotnn_.row(neighKnot) - xdotnn_.row(ownerKnot);

        v_ij = 
            1.0/coeffs_(0)*
            (   
                -coeffs_(1)*v_ij - coeffs_(2)*v_ijn - coeffs_(3)*v_ijnn
            ); 
 
     
        // Bar vector
        x_ij = x_.row(neighKnot) - x_.row(ownerKnot);
        x_ijn = xn_.row(neighKnot) - xn_.row(ownerKnot);
        x_ijnn = xnn_.row(neighKnot) - xnn_.row(ownerKnot);

        x_ij = 
            1.0/coeffs_(0)*
            (   
                -coeffs_(1)*x_ij - coeffs_(2)*x_ijn - coeffs_(3)*x_ijnn
            );  

        
        // Length of original bar vector
        l0_ij = (x0_.row(neighKnot) - x0_.row(ownerKnot)).norm();
     
  
        // Fill rhs
        B_(barI) = 
            coeffs_(0)*coeffs_(0)/2.0*
            (
                l0_ij*l0_ij - x_ij.dot(x_ij) - 2.0/coeffs_(0)*v_ij.dot(x_ij)
            ); 
        
        if (mass_j > 0.0)
        {          
            B_(barI) -= F_j.dot(x_ij)/mass_j;
        }
        else
        {
            B_(barI) -= xdotdot_.row(ownerKnot).dot(x_ij);
        }
        
        if (mass_i > 0.0)
        {
            B_(barI) += F_i.dot(x_ij)/mass_i;
        }       
        else
        {
            B_(barI) += xdotdot_.row(neighKnot).dot(x_ij);
        }
    }    
}
