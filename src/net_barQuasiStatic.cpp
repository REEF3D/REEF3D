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

#include"net_barQuasiStatic.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"reinidisc_fsf.h"	
#include"vrans.h"

net_barQuasiStatic::net_barQuasiStatic(int number, lexer *p):nNet(number),f_(p),dt(p),frk1(p),frk2(p),L_(p), cutl(p), cutr(p){}

net_barQuasiStatic::~net_barQuasiStatic(){}


void net_barQuasiStatic::initialize(lexer *p, fdm *a, ghostcell *pgc)
{
    prdisc = new reinidisc_fsf(p);    
   
    //- Initialise net model
    if (p->X320_type[nNet] == 1)
    {
        bag_ini(p,a,pgc);
        
        buildNet_bag(p);
    }
    else if (p->X320_type[nNet] == 2)   
    {
        cyl_ini(p,a,pgc);
        
        buildNet_cyl(p); 
    }
    else if (p->X320_type[nNet] == 3)   
    {
        wall_ini(p,a,pgc);
        
        buildNet_wall(p);    
    }  
    
    //- Update porous zone
    vransCoupling(p,a,pgc);
}


void net_barQuasiStatic::start
(
	lexer *p, 
	fdm *a, 
	ghostcell *pgc,
    double alpha,
    Eigen::Matrix3d quatRotMat 
)
{
    double starttime1=pgc->timer();     
    
    double norm, error;
    int iter = 1;
   
    //- Get velocities at knots
    
    updateField(p, a, pgc, 0);
    updateField(p, a, pgc, 1);	
    updateField(p, a, pgc, 2);
    
    //- Get density at knots
    
    updateField(p, a, pgc, 3);
    

	//- Solving the system of equations
   
	while(iter < 200)
	{
        //- Fill right-hand side Bh with gravity and hydrodynamic forces
        if (p->X320_type[nNet] == 1)
        {
            fillRhs_bag(p);
        }
        else
        {
            // fillRhs_Morison(p);  
            fillRhs_Screen(p);
        }
		
        //- Solve the system A * fi = Bh
        // fi = A.lu().solve(Bh); 

		//-  Correct system such that length of normal vectors equals one
		error = 0.0;
		for (int j = 0; j < nf; j++)
		{
            norm = fi.row(j).norm();

            fi(j,0) /= norm;
            fi(j,1) /= norm;
            fi(j,2) /= norm;
        
			for (int k = 0; k < niK; k++) 
			{
				A(k,j) *= norm;
			}
            
            error = max(error,fabs(norm-1.0));
		}	

        //- Check convergence
		if (error < 1e-2)
        {
            break;
        }
        else
        {
            iter++;
        }

		//- Correct length of bars
        updateLength();
	}
   
    if (p->mpirank == 0)
    {
        cout<<"Number of iterations = "<<iter<<setprecision(5)<<" with error = "<<error<<endl;
    }
    
    
	//- Build and save net	
  
	print(p);	
 
    
    //- Update porous zone and coefficients
    
    vransCoupling(p,a,pgc);
    

    double endtime1 = pgc->timer()-starttime1; 
    if (p->mpirank == 0) cout<<"Net time: "<<endtime1<<endl;    
}


void net_barQuasiStatic::fillRhs_Morison(lexer *p)
{
    double vn_mag, cn, ct;
    double vel_x, vel_y, vel_z, vt_fi;

    int index = 0;

    //- Calculate force directions and coefficients using Morison formula
    
    for (int j = 0; j < nf; j++)
    {
        vel_x = (coupledField[Pi[j]][0] + coupledField[Ni[j]][0])/2.0;
        vel_y = (coupledField[Pi[j]][1] + coupledField[Ni[j]][1])/2.0;
        vel_z = (coupledField[Pi[j]][2] + coupledField[Ni[j]][2])/2.0;    

        // tangential velocity v_t = scalarProd(v,fi)*fi
        vt_fi = vel_x*fi(j,0) + vel_y*fi(j,1) + vel_z*fi(j,2);

        v_t[j][0] = vt_fi*fi(j,0);
        v_t[j][1] = vt_fi*fi(j,1);
        v_t[j][2] = vt_fi*fi(j,2);

        // normal velocity v_n = v - v_t
        v_n[j][0] = vel_x - v_t[j][0];
        v_n[j][1] = vel_y - v_t[j][1];
        v_n[j][2] = vel_z - v_t[j][2];   
    }	

    
    for (int i = 0; i < niK; i++)
    {
        // Assign gravity force to knot
        Bh.row(i) = B.row(i);
        
        // Assign hydrodynamic forces to knot from each adjoint bar
        for (int k = 1; k < 5; k++)
        {
            int& nfKik = nfK[i][k];
            
            if (nfKik > -1)
            {
                vn_mag = 
                    sqrt
                    (
                        v_n[nfKik][0]*v_n[nfKik][0] 
                      + v_n[nfKik][1]*v_n[nfKik][1] 
                      + v_n[nfKik][2]*v_n[nfKik][2]
                    );
                        
                morisonForceCoeff(cn,ct,vn_mag);
                      
                Bh(index,0) -= 
                    (
                        p->W1/2.0*d_c*l[nfKik]/2.0*cn*vn_mag*v_n[nfKik][0] 
                      + ct*v_t[nfKik][0]
                    );
       /*        Bh(index,1) -= 
                    (
                        p->W1/2.0*d_c*l[nfKik]/2.0*cn*vn_mag*v_n[nfKik][1] 
                      + ct*v_t[nfKik][1]
                    );
                Bh(index,2) -= 
                    (
                        p->W1/2.0*d_c*l[nfKik]/2.0*cn*vn_mag*v_n[nfKik][2] 
                      + ct*v_t[nfKik][2]
                    ); */
            }   
        }

        index++;
    }
}


void net_barQuasiStatic::fillRhs_Screen(lexer *p)
{
    // Screen force model of Kristiansen (2012)
    // Assign hydrodynamic forces to knot from each adjoint screen
    
    Fx = 0.0;
    Fy = 0.0;
    Fz = 0.0;
    
    int nBars;
    int index = 0;
    
    Vector3d v_rel(1,3);
    Vector3d n_v(1,3);
    double v_mag, rho;
    
    for (int i = 0; i < niK; i++)
    {
        int*& barsiKI = nfK[i];
        int& kI = barsiKI[0];
        
        // Assign gravity force to knot
        Bh.row(i) = B.row(i);

        // Count number of screens
        nBars = 4;
        
        for (int k = 1; k < 5; k++)
        {
            if (barsiKI[k] == -1) nBars--;
        }
        
        // Access velocity at knot
        v_rel << coupledField[kI][0], coupledField[kI][1], coupledField[kI][2];

        v_mag = v_rel.norm();

        n_v = v_rel/(v_mag + 1e-10);
        
        // Access density at knot
        rho = coupledField[kI][3];
 
        // Calculate hydrodynamic forces
        if (nBars == 2)         // Corner screens
        {
            Bh.row(index) -= screenForce(p,rho,v_rel,n_v,v_mag,barsiKI[1],barsiKI[2]);
        }      
        else if (nBars == 3)    // Edge screens
        {
            Bh.row(index) -= screenForce(p,rho,v_rel,n_v,v_mag,barsiKI[1],barsiKI[3]);
            Bh.row(index) -= screenForce(p,rho,v_rel,n_v,v_mag,barsiKI[2],barsiKI[3]);
        }
        else                    // Inner screens
        {
            Bh.row(index) -= screenForce(p,rho,v_rel,n_v,v_mag,barsiKI[1],barsiKI[4]);
            Bh.row(index) -= screenForce(p,rho,v_rel,n_v,v_mag,barsiKI[4],barsiKI[2]);
            Bh.row(index) -= screenForce(p,rho,v_rel,n_v,v_mag,barsiKI[2],barsiKI[3]);
            Bh.row(index) -= screenForce(p,rho,v_rel,n_v,v_mag,barsiKI[1],barsiKI[3]);            
        }

        index++;        
    }
}


void net_barQuasiStatic::fillRhs_bag(lexer *p)
{
    double vn_mag, cn, ct;
    
	int index = 0;
    bool bk;
		
    for (int j = 0; j < nK; j++)
    {
        bk = false;
		
        for (int k = 0; k < 2*nd+2*nl; k++)
        {
            if (j == Pb[k] || j == Nb[k])
            {
                bk = true;
                break;
            }
        }
			
        if (bk == false)
        {
            Bh.row(index) = B.row(index);

            for (int k = 0; k < 4; k++)
            {
                int nfKik = nfK[index][k]; 	

                vn_mag = sqrt(v_n[nfKik][0]*v_n[nfKik][0] + v_n[nfKik][1]*v_n[nfKik][1] + v_n[nfKik][2]*v_n[nfKik][2]);
                    
                morisonForceCoeff(cn,ct,vn_mag);
					
                Bh(index,0) -= (p->W1/2.0*d_c*l[nfKik]/2.0*cn*vn_mag*v_n[nfKik][0] + ct*v_t[nfKik][0]);
                Bh(index,1) -= (p->W1/2.0*d_c*l[nfKik]/2.0*cn*vn_mag*v_n[nfKik][1] + ct*v_t[nfKik][1]);
                Bh(index,2) -= (p->W1/2.0*d_c*l[nfKik]/2.0*cn*vn_mag*v_n[nfKik][2] + ct*v_t[nfKik][2]);
            }
							
            index++;
        }
    } 
}
