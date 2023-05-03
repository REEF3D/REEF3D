/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include<iomanip>
#include"sflow_print_wsfline_y.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"wave_interface.h"
#include<sys/stat.h>
#include<sys/types.h>

sflow_print_wsfline_y::sflow_print_wsfline_y(lexer *p, fdm2D *b, ghostcell *pgc)
{	
	p->Iarray(iloc,p->P56);

    maxknoy=pgc->globalimax(p->knoy);
    sumknoy=pgc->globalisum(maxknoy);
	
    p->Darray(yloc,p->P56+2,maxknoy);
    p->Darray(wsf,p->P56+2,maxknoy);
    p->Iarray(flag,p->P56+2,maxknoy);
	p->Iarray(wsfpoints,p->P56+2);
	

    p->Darray(yloc_all,p->P56+2,sumknoy);
    p->Darray(wsf_all,p->P56+2,sumknoy);
	p->Iarray(flag_all,p->P56+2,sumknoy);
	p->Iarray(rowflag,sumknoy);

    for(q=0;q<p->P56;++q)
    for(n=0;n<maxknoy;++n)
    {
    yloc[q][n]=0.0;
    wsf[q][n]=0.0;
    }

    for(q=0;q<p->P56;++q)
    for(n=0;n<sumknoy;++n)
    {
    yloc_all[q][n]=0.0;
    wsf_all[q][n]=0.0;
	flag_all[q][n]=0;
	rowflag[n]=0;
    }

    ini_location(p,b,pgc);
	
	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_SFLOW_WSFLINE_Y",0777);
}

sflow_print_wsfline_y::~sflow_print_wsfline_y()
{
    wsfout.close();
}

void sflow_print_wsfline_y::start(lexer *p, fdm2D *b, ghostcell *pgc, ioflow *pflow, slice &f)
{
	
    char name[250];
    double zval=0.0;
    int num,check;
	
    num = p->count;

    if(p->mpirank==0)
    {
		// open file
		if(p->P14==0)
		{
		if(num<10)
		sprintf(name,"REEF3D-SFLOW-wsfline_y-00000%i.dat",num);

		if(num<100&&num>9)
		sprintf(name,"REEF3D-SFLOW-wsfline_y-0000%i.dat",num);

		if(num<1000&&num>99)
		sprintf(name,"REEF3D-SFLOW-wsfline_y-000%i.dat",num);

		if(num<10000&&num>999)
		sprintf(name,"REEF3D-SFLOW-wsfline_y-00%i.dat",num);

		if(num<100000&&num>9999)
		sprintf(name,"REEF3D-SFLOW-wsfline_y-0%i.dat",num);

		if(num>99999)
		sprintf(name,"REEF3D-SFLOW-wsfline_y-%i.dat",num);
		}
		
		if(p->P14==1)
		{
		if(num<10)
		sprintf(name,"./REEF3D_SFLOW_WSFLINE_Y/REEF3D-SFLOW-wsfline_y-00000%i.dat",num);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_SFLOW_WSFLINE_Y/REEF3D-SFLOW-wsfline_y-0000%i.dat",num);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_SFLOW_WSFLINE_Y/REEF3D-SFLOW-wsfline_y-000%i.dat",num);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_SFLOW_WSFLINE_Y/REEF3D-SFLOW-wsfline_y-00%i.dat",num);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_SFLOW_WSFLINE_Y/REEF3D-SFLOW-wsfline_y-0%i.dat",num);

		if(num>99999)
		sprintf(name,"./REEF3D_SFLOW_WSFLINE_Y/REEF3D-SFLOW-wsfline_y-%i.dat",num);
		}
		
		wsfout.open(name);

		wsfout<<"simtime:  "<<p->simtime<<endl;
		wsfout<<"number of wsf-lines:  "<<p->P56<<endl<<endl;
		wsfout<<"line_No     x_coord"<<endl;
		for(q=0;q<p->P56;++q)
		wsfout<<q+1<<"\t "<<p->P56_x[q]<<endl;

		if(p->P53==1)
		wsfout<<q+1<<"\t "<<" Wave Theory "<<endl;

		wsfout<<endl<<endl;

		
		for(q=0;q<p->P56;++q)
		{
		wsfout<<"Y "<<q+1;
		wsfout<<"\t P "<<q+1<<" \t \t ";
		if(p->P53==1)
		wsfout<<"\t \t W "<<q+1;
		}

		wsfout<<endl<<endl;
    }

    //-------------------

    for(q=0;q<p->P56;++q)
    for(n=0;n<maxknoy;++n)
    {
    yloc[q][n]=1.0e20;
    wsf[q][n]=-1.0e20;
    }

    for(q=0;q<p->P56;++q)
    {
        JLOOP
        if(flag[q][j]>0)
        {
        i=iloc[q];


            wsf[q][j]=f(i,j)+p->phimean;
            yloc[q][j]=p->pos_y();
        }
    }
	
	
	for(q=0;q<p->P56;++q)
    wsfpoints[q]=sumknoy;
	
    // gather
    for(q=0;q<p->P56;++q)
    {
    pgc->gather_double(yloc[q],maxknoy,yloc_all[q],maxknoy);
    pgc->gather_double(wsf[q],maxknoy,wsf_all[q],maxknoy);
	pgc->gather_int(flag[q],maxknoy,flag_all[q],maxknoy);

		
        if(p->mpirank==0)
        {
        sort(yloc_all[q], wsf_all[q], flag_all[q], 0, wsfpoints[q]-1);
        remove_multientry(p,yloc_all[q], wsf_all[q], flag_all[q], wsfpoints[q]); 
        }
		
    }
	
    // write to file
    if(p->mpirank==0)
    {
		for(n=0;n<sumknoy;++n)
		rowflag[n]=0;
		
		for(n=0;n<sumknoy;++n)
        {
			check=0;
		    for(q=0;q<p->P56;++q)
			if(flag_all[q][n]>0 && yloc_all[q][n]<1.0e20)
			check=1;
			
			if(check==1)
			rowflag[n]=1;
		}
        

        for(n=0;n<sumknoy;++n)
        {
			check=0;
		    for(q=0;q<p->P56;++q)
			{
				if(flag_all[q][n]>0 && yloc_all[q][n]<1.0e20)
				{
				wsfout<<setprecision(5)<<yloc_all[q][n]<<" \t ";
				wsfout<<setprecision(5)<<wsf_all[q][n]<<" \t  ";
				
				
					if(p->P53==1)
					wsfout<<pflow->wave_fsf(p,pgc,yloc_all[q][n])<<" \t  ";
					
				check=1;
				}
				
				if((flag_all[q][n]<0 || yloc_all[q][n]>=1.0e20) && rowflag[n]==1)
				{
					wsfout<<setprecision(5)<<" \t ";
					wsfout<<setprecision(5)<<" \t ";
					
				}
			}

            
			if(check==1)
            wsfout<<endl;
        }

    wsfout.close();
    }
}

void sflow_print_wsfline_y::ini_location(lexer *p, fdm2D *b, ghostcell *pgc)
{
    int check,count;
    
    
    for(q=0;q<p->P56;++q)
    {
        count=0;
        JLOOP
        {
        iloc[q]=p->posc_i(p->P56_x[q]);

        if(iloc[q]>=0 && iloc[q]<p->knox)
        flag[q][count]=1;
        ++count;
        }
    }
}

void sflow_print_wsfline_y::sort(double *a, double *b, int *c, int left, int right)
{

  if (left < right)
  {

    double pivot = a[right];
    int l = left;
    int r = right;

    do {
      while (a[l] < pivot) l++;

      while (a[r] > pivot) r--;

      if (l <= r) {
          double swap = a[l];
          double swapd = b[l];
		  int swapc = c[l];

          a[l] = a[r];
          a[r] = swap;

          b[l] = b[r];
          b[r] = swapd;
		  
		  c[l] = c[r];
          c[r] = swapc;

          l++;
          r--;
      }
    } while (l <= r);

    sort(a,b,c, left, r);
    sort(a,b,c, l, right);
  }
}

void sflow_print_wsfline_y::remove_multientry(lexer *p, double* b, double* c, int *d, int& num)
{
    int oldnum=num;
    double xval=-1.12e23;

    int count=0;

    double *f,*g;
	int *h;
	
	p->Darray(f,num);
	p->Darray(g,num);
	p->Iarray(h,num);

    for(n=0;n<num;++n)
    g[n]=-1.12e22;


    for(n=0;n<oldnum;++n)
    {
        if(xval<=b[n]+0.001*p->DXM && xval>=b[n]-0.001*p->DXM && count>0)
        g[count-1]=MAX(g[count-1],c[n]);

        if(xval>b[n]+0.001*p->DXM || xval<b[n]-0.001*p->DXM)
        {
        f[count]=b[n];
        g[count]=c[n];
		h[count]=d[n];
        ++count;
        }

    xval=b[n];
    }

    for(n=0;n<count;++n)
    {
    b[n]=f[n];
    c[n]=g[n];
	d[n]=h[n];
    }

    
    p->del_Darray(f,num);
	p->del_Darray(g,num);
	p->del_Iarray(h,num);
	
	num=count;

}


int sflow_print_wsfline_y::conv(double a)
{

int b,c;
double d,diff;

c= int( a);
d=double(c);
diff=a-d;

b=c;

if(diff>0.5)
b=c+1;

if(diff<=-0.5)
b=c-1;

return b;

}
