/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"benchmark_TaylorGreen.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<sys/stat.h>

benchmark_TaylorGreen::benchmark_TaylorGreen(lexer *p, fdm *a) : gradient(p)
{
/*    double L = 1.0;
    double U = 0.024256;
    
    double x,y,z;
*/
	double x,y,z,b,r,x_0,z_0;
    double L = 1.0;
    double U_inf = 1.0;
    double U_a = 0.8;
    b = 1.0 / sqrt(log(2.0));
//    x_0 = -18.75 * L;
    x_0 = 0.0 * L;
    z_0 = 0.0;
	
	
    ULOOP
    {
/*        x = p->pos1_x();
        y = p->pos1_y();
        z = p->pos1_z();
        a->u(i,j,k) = U*(sin(x/L)*cos(y/L)*cos(z/L));
//        a->u(i,j,k) = 1.0;
*/
        x = p->pos1_x();
        y = p->pos1_y();
        z = p->pos1_z();
		r = sqrt(pow(x-x_0,2) + pow(z-z_0,2));
        a->u(i,j,k) = U_inf + (U_a * exp((1.0 - pow(r/b,2.0))/2.0) * ((z - z_0) / b));

    }   
    
    VLOOP
    {
/*        x = p->pos2_x();
        y = p->pos2_y();
        z = p->pos2_z();
        a->v(i,j,k) = -U*(cos(x/L)*sin(y/L)*cos(z/L));
//        a->v(i,j,k) = 1.0;
*/
        a->v(i,j,k) = 0.0;
    }
    
    WLOOP
    {
//        a->w(i,j,k) = 0.0;
        x = p->pos3_x();
        y = p->pos3_y();
        z = p->pos3_z();
		r = sqrt(pow(x-x_0,2) + pow(z-z_0,2));
        a->w(i,j,k) = (-U_a * exp((1.0 - pow(r/b,2.0))/2.0) * ((x - x_0) / b));
    }

    if(p->mpirank==0)
    {
        mkdir("./REEF3D_CFD_TaylorGreen_Diss", 0777);
	
        ofstream print;
        print.open("./REEF3D_CFD_TaylorGreen_Diss/REEF3D_CFD_TG_Diss.dat");
        print<<"time \t epsilon "<<endl;
        print.close();


	print.open("./REEF3D_CFD_TaylorGreen_Diss/REEF3D_CFD_TG_Energy.dat");
        print<<"time \t energy "<<endl;
        print.close();

    }

    	

    if(p->mpirank==0)
    {
        mkdir("./REEF3D_CFD_Data",0777);
        mkdir("./REEF3D_CFD_position",0777);
        mkdir("./REEF3D_CFD_Growth_Ratio", 0777);
        mkdir("./REEF3D_CFD_RMS_error",0777);
    }

}

benchmark_TaylorGreen::~benchmark_TaylorGreen()
{
}

void benchmark_TaylorGreen::start(lexer* p, fdm *a, ghostcell *pgc, convection *pconvec )
{
//Data file	
	char name[250];
        int numberi = p->count;
	
	ofstream printdata, printpos;
//	double x,y,z,xp,yp,zp;
	double X1NM1P,growfactor1;

//Constansts for Analytical Solution of Isentropic Vortex
    double x,y,z,r,w_analytical,error;
    double L = 1.0;
    double U_inf = 1.0;
    double U_a = 0.8;
    double b = 1.0 / sqrt(log(2.0));
    double x_0 = 0.0 * L;
    double z_0 = 0.0;
    double total_error = 0.0;
    double nn = 0.0;
	
//Variables for Taylor Green test case
/*    double vx, vy, vz, uu, vv, ww, vmag2, velmag2, vol;

    double vVolAvg = 0.0;
    double velVolAvg = 0.0;
    double volTot = 0.0;

    // Local calculation of volume averaged vorticity
    LOOP
    {
        if(p->j_dir==1)
        {
            vx = pvdz(p,a) - pwdy(p,a); 
            vy = pudz(p,a) - pwdx(p,a);
            vz = pvdx(p,a) - pudy(p,a);
            uu = 0.5*(a->u(i,j,k) + a->u(i-1,j,k));
	    vv = 0.5*(a->v(i,j,k) + a->v(i,j-1,k));
	    ww = 0.5*(a->w(i,j,k) + a->w(i,j,k-1));
        }
        else
        {
            vx = 0.0;
            vy = pudz(p,a) - pwdx(p,a);
            vz = 0.0;
        }
       
        vmag2 = vx*vx + vy*vy + vz*vz;

	velmag2 = uu*uu + vv*vv + ww*ww;
        
        vol = p->DXN[IP]*p->DYN[JP]*p->DZN[KP];

        vVolAvg += vmag2*vol;

	velVolAvg += velmag2*vol;

        volTot += vol; 
    }
*/


//RMS_error Calculation of Isentropic Vortex after 50s, time step=0.01
    if(numberi == 5001.0)
    {
    WLOOP
    {
	//Analytical Solution
        x = p->pos3_x();
        y = p->pos3_y();
        z = p->pos3_z();
	r = sqrt(pow(x-x_0,2) + pow(z-z_0,2));
        w_analytical = (-U_a * exp((1.0 - pow(r/b,2.0))/2.0) * ((x - x_0) / b));

	error = pow((a->w(i,j,k) - w_analytical),2.0);
	total_error += error;

	nn++;
        
    }
    }


    // Sum up vorticity and volume and send to rank 0 
//    double dissipation = 0.0;
//    double energy = 0.0;
//    double volume = 0.0;
    double RMS_error = 0.0;
    double Nt = 0.0;
//    MPI_Reduce(&vVolAvg, &dissipation, 1, MPI_DOUBLE, MPI_SUM, 0, pgc->mpi_comm);
//    MPI_Reduce(&velVolAvg, &energy, 1, MPI_DOUBLE, MPI_SUM, 0, pgc->mpi_comm);
//    MPI_Reduce(&volTot, &volume, 1, MPI_DOUBLE, MPI_SUM, 0, pgc->mpi_comm);
    MPI_Reduce(&total_error, &RMS_error, 1, MPI_DOUBLE, MPI_SUM, 0, pgc->mpi_comm);
    MPI_Reduce(&nn, &Nt, 1, MPI_DOUBLE, MPI_SUM, 0, pgc->mpi_comm);

    if(p->mpirank == 0)
    {
//        dissipation = p->W2*dissipation/volume;

//	energy = 0.5*energy/volume;

	RMS_error = sqrt(RMS_error/Nt);
        
        ofstream print;
/*        print.open("./REEF3D_CFD_TaylorGreen_Diss/REEF3D_CFD_TG_Diss.dat", ofstream::app);
        print<<setprecision(15)<<p->simtime<<" \t "<<dissipation<<endl;
        print.close();


        print.open("./REEF3D_CFD_TaylorGreen_Diss/REEF3D_CFD_TG_Energy.dat", ofstream::app);
        print<<setprecision(15)<<p->simtime<<" \t "<<energy<<endl;
        print.close();
*/
	if(numberi == 5001.0 || numberi == 10001.0)
        print.open("./REEF3D_CFD_RMS_error/REEF3D_CFD_RMSerror.dat", ofstream::app);
        print<<setprecision(15)<<p->simtime<<" \t "<<RMS_error<<endl;
	print<<setprecision(15)<<p->simtime<<" \t "<<Nt<<endl;
        print.close();
    }



    //W_Velocity_X_Line
//    if(p->mpirank==0)
//    {
            //ofstream printdata;
			//double x,y,z;
			// open file
/*
			if(numberi == 5001.0 || numberi == 10001.0)
			{
			
			if(numberi<10)
			sprintf(name,"./REEF3D_CFD_Data/REEF3D-CFD-WX-%i-00000%i.dat",numberi,p->mpirank);

			if(numberi<100&&numberi>9)
			sprintf(name,"./REEF3D_CFD_Data/REEF3D-CFD-WX-%i-0000%i.dat",numberi,p->mpirank);

			if(numberi<1000&&numberi>99)
			sprintf(name,"./REEF3D_CFD_Data/REEF3D-CFD-WX-%i-000%i.dat",numberi,p->mpirank);

			if(numberi<10000&&numberi>999)
			sprintf(name,"./REEF3D_CFD_Data/REEF3D-CFD-WX-%i-00%i.dat",numberi,p->mpirank);

			if(numberi<100000&&numberi>9999)
			sprintf(name,"./REEF3D_CFD_Data/REEF3D-CFD-WX-%i-0%i.dat",numberi,p->mpirank);

			if(numberi>99999)
			sprintf(name,"./REEF3D_CFD_Data/REEF3D-CFD-WX-%i-%i.dat",numberi,p->mpirank);
			
			
			printdata.open(name);

			printdata<<"Data ID:  "<<1<<endl<<endl;
//			lineout[n]<<"x_start     x_end     y_start     y_end     z_start     z_end"<<endl;
			
//			lineout[n]<<p->P62_xs[n]<<"\t "<<p->P62_xe[n]<<"\t "<<p->P62_ys[n]<<"\t "<<p->P62_ye[n]<<"\t "<<p->P62_zs[n]<<"\t "<<p->P62_ze[n]<<endl;

			printdata<<endl<<endl;
			
            printdata<<"X \t Y \t Z \t W \t"<<endl;
//			printdata<<"X \t Y \t Z \t U \t V \t W \t P"<<endl;
			
			printdata<<endl<<endl;

			for(i=0; i<p->knox; ++i)
			for(j=0; j<p->knoy; ++j)
			for(k=0; k<p->knoz-p->wlast; ++k)
			if (j==0 && k==p->knoz/2.0)
			{
			
			
//			x = i;
//			y = j;
//			z = k;

			x = p->pos3_x();
        		y = p->pos3_y();
        		z = p->pos3_z();
			
            printdata<<x<<" \t "<<y<<" \t "<<z<<" \t "<<a->w(i,j,k)<<" \t "<<endl;
//          printdata<<p->pos1_x()<<" \t "<<p->pos1_y()<<" \t "<<p->pos1_z()<<" \t "<<a->u(i,j,k)<<" \t "<<V<<" \t "<<W;
//			printdata<<" \t "<<P<<endl;
			}
			
			printdata.close();
			
			
			
    }

    //Growth_Ratio
//    if(p->mpirank==0)
//    {
            //ofstream printdata;
			//double x,y,z;
			// open file

			if(numberi == 1.0 || numberi == 2.0)
			{
			
			if(numberi<10)
			sprintf(name,"./REEF3D_CFD_Growth_Ratio/REEF3D_CFD_Gr_Ra-%i-00000%i.dat",numberi,p->mpirank);

			if(numberi<100&&numberi>9)
			sprintf(name,"./REEF3D_CFD_Growth_Ratio/REEF3D_CFD_Gr_Ra-%i-0000%i.dat",numberi,p->mpirank);

			if(numberi<1000&&numberi>99)
			sprintf(name,"./REEF3D_CFD_Growth_Ratio/REEF3D_CFD_Gr_Ra-%i-000%i.dat",numberi,p->mpirank);

			if(numberi<10000&&numberi>999)
			sprintf(name,"./REEF3D_CFD_Growth_Ratio/REEF3D_CFD_Gr_Ra-%i-00%i.dat",numberi,p->mpirank);

			if(numberi<100000&&numberi>9999)
			sprintf(name,"./REEF3D_CFD_Growth_Ratio/REEF3D_CFD_Gr_Ra-%i-0%i.dat",numberi,p->mpirank);

			if(numberi>99999)
			sprintf(name,"./REEF3D_CFD_Growth_Ratio/REEF3D_CFD_Gr_Ra-%i-%i.dat",numberi,p->mpirank);
			
			
			printdata.open(name);

			printdata<<"Data ID:  "<<1<<endl<<endl;
//			lineout[n]<<"x_start     x_end     y_start     y_end     z_start     z_end"<<endl;
			
//			lineout[n]<<p->P62_xs[n]<<"\t "<<p->P62_xe[n]<<"\t "<<p->P62_ys[n]<<"\t "<<p->P62_ye[n]<<"\t "<<p->P62_zs[n]<<"\t "<<p->P62_ze[n]<<endl;

			printdata<<endl<<endl;
			
            printdata<<"position \t GrowthRatio \t"<<endl;
//			printdata<<"X \t Y \t Z \t U \t V \t W \t P"<<endl;
			
			printdata<<endl<<endl;

			for(i=1;i<p->knox-1;++i)
			for(j=0;j<p->knoy;++j)
			for(k=0;k<p->knoz;++k)
			if (j==0 && k==p->knoz/2.0)
			{
			
			
//			x = i;
//			y = j;
//			z = k;

//			x = p->pos3_x();
//       		y = p->pos3_y();
//        		z = p->pos3_z();

			X1NM1P = p->XP[IP];
        if (p->XP[IP] < 0)
        {
            growfactor1 = 1.0/(p->DXP[IP]/p->DXP[IM1]);
        }
        else
        {
            growfactor1 = (p->DXP[IP]/p->DXP[IM1]);
        }
        
        printdata<<X1NM1P<<" \t "<<growfactor1<<" \t "<<endl;

        
			
//            printdata<<x<<" \t "<<y<<" \t "<<z<<" \t "<<a->w(i,j,k)<<" \t "<<endl;
//          printdata<<p->pos1_x()<<" \t "<<p->pos1_y()<<" \t "<<p->pos1_z()<<" \t "<<a->u(i,j,k)<<" \t "<<V<<" \t "<<W;
//			printdata<<" \t "<<P<<endl;
			}
			
			printdata.close();
			
			
			
    }




*/


}


/*
benchmark_TaylorGreen::benchmark_TaylorGreen(lexer *p, fdm *a) : gradient(p)
{
    double L = 1.0;
    double U = 1600.0;
    
    double x,y,z;
    
    ULOOP
    {
        x = p->pos1_x();
        y = p->pos1_y();
        z = p->pos1_z();
        a->u(i,j,k) = U*(sin(x/L)*cos(y/L)*cos(z/L));
    }   
    
    VLOOP
    {
        x = p->pos2_x();
        y = p->pos2_y();
        z = p->pos2_z();
        a->v(i,j,k) = -U*(cos(x/L)*sin(y/L)*cos(z/L));
    }
    
    WLOOP
    a->w(i,j,k) = 0.0;
    
    if(p->mpirank==0)
    {
        mkdir("./REEF3D_CFD_TaylorGreen_Diss", 0777);
	
        ofstream print;
        print.open("./REEF3D_CFD_TaylorGreen_Diss/REEF3D_CFD_TG_Diss.dat");
        print<<"time \t epsilon "<<endl;
        print.close();
    }
}

benchmark_TaylorGreen::~benchmark_TaylorGreen()
{
}

void benchmark_TaylorGreen::start(lexer* p, fdm *a, ghostcell *pgc, convection *pconvec )
{
    double vx, vy, vz, vmag2, vol;

    double vVolAvg = 0.0;
    double volTot = 0.0;

    // Local calculation of volume averaged vorticity
    LOOP
    {
        if(p->j_dir==1)
        {
            vx = pvdz(p,a) - pwdy(p,a); 
            vy = pudz(p,a) - pwdx(p,a);
            vz = pvdx(p,a) - pudy(p,a); 
        }
        else
        {
            vx = 0.0;
            vy = pudz(p,a) - pwdx(p,a);
            vz = 0.0;
        }
       
        vmag2 = vx*vx + vy*vy + vz*vz;
        
        vol = p->DXN[IP]*p->DYN[JP]*p->DZN[KP];

        vVolAvg += vmag2*vol;
        volTot += vol; 
    }

    // Sum up vorticity and volume and send to rank 0 
    double dissipation = 0.0;
    double volume = 0.0;
    MPI_Reduce(&vVolAvg, &dissipation, 1, MPI_DOUBLE, MPI_SUM, 0, pgc->mpi_comm);
    MPI_Reduce(&volTot, &volume, 1, MPI_DOUBLE, MPI_SUM, 0, pgc->mpi_comm);

    if(p->mpirank == 0)
    {
        dissipation = p->W2*dissipation/volume;
        
        ofstream print;
        print.open("./REEF3D_CFD_TaylorGreen_Diss/REEF3D_CFD_TG_Diss.dat", ofstream::app);
        print<<p->simtime<<" \t "<<dissipation<<endl;
        print.close();
    }
}*/
