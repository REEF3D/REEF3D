/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"

void ghostcell::gcfb_velupdate(lexer *p, fdm *a, int **cellmem, int cellcount,double xdelt, double ydelt, double zdelt, int dir)
{
    int nn,aa,bb,cc;

    for(nn=0;nn<cellcount;++nn)
	{
		// new fluid cell
		if(cellmem[nn][3]==1) 
		{
		i=cellmem[nn][0];
		j=cellmem[nn][1];
		k=cellmem[nn][2];
		
		if(dir==1)
		GC1LOOP
		if(p->gcb1[n][4]==41)
		if(i==p->gcb1[n][0] && j==p->gcb1[n][1] && k==p->gcb1[n][2])
		p->gcb1[n][4]=42;
		
		if(dir==2)
		GC2LOOP
		if(p->gcb2[n][4]==41)
		if(i==p->gcb2[n][0] && j==p->gcb2[n][1] && k==p->gcb2[n][2])
		p->gcb2[n][4]=42;
		
		if(dir==3)
		GC3LOOP
		if(p->gcb3[n][4]==41)
		if(i==p->gcb3[n][0] && j==p->gcb3[n][1] && k==p->gcb3[n][2])
		p->gcb3[n][4]=42;
		}
        
	
		// new solid cell
		if(cellmem[nn][3]==2)
		{
		i=cellmem[nn][0];
		j=cellmem[nn][1];
		k=cellmem[nn][2];
		
		//cout<<"$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
		
		if(dir==1)
		GC1LOOP
		if(p->gcb1[n][4]==41)
		{
				
			aa=bb=cc=0;
			if(p->gcb1[n][3]==1)
			aa=-1;
			
			if(p->gcb1[n][3]==4)
			aa=1;
			
			if(p->gcb1[n][3]==3)
			bb=-1;
			
			if(p->gcb1[n][3]==2)
			bb=1;
			
			if(p->gcb1[n][3]==5)
			cc=-1;
			
			if(p->gcb1[n][3]==6)
			cc=1;
			
			
			if(i-aa==p->gcb1[n][0] && j-bb==p->gcb1[n][1] && k-cc==p->gcb1[n][2])
			{
			p->gcb1[n][4]=43;
			//cout<<"++++++++++++++++++++++++"<<endl;
			}
		}
		
		
		if(dir==2)
		GC2LOOP
		if(p->gcb2[n][4]==41)
		if(i==p->gcb2[n][0] && j==p->gcb2[n][1] && k==p->gcb2[n][2])
		p->gcb2[n][4]=43;
		
		if(dir==3)
		GC3LOOP
		if(p->gcb3[n][4]==41)
		{
				
			aa=bb=cc=0;
			if(p->gcb3[n][3]==1)
			aa=-1;
			
			if(p->gcb3[n][3]==4)
			aa=1;
			
			if(p->gcb3[n][3]==3)
			bb=-1;
			
			if(p->gcb3[n][3]==2)
			bb=1;
			
			if(p->gcb3[n][3]==5)
			cc=-1;
			
			if(p->gcb3[n][3]==6)
			cc=1;
			
			
			if(i-aa==p->gcb3[n][0] && j-bb==p->gcb3[n][1] && k-cc==p->gcb3[n][2])
			{
			p->gcb3[n][4]=43;
			//cout<<"++++++++++++++++++++++++"<<endl;
			}
		}

		}
	}
}

void ghostcell::gcfb_scalarupdate(lexer *p, fdm *a, int **cellmem, int cellcount, field &f)
{
    int nn;
	double nx,ny,nz,norm;
	double posx,posy,posz;
	double locx,locy,locz;
	double fbval,fval;
	double ui,vi,wi,un,vn,wn;


    for(nn=0;nn<cellcount;++nn)
	{
		if(cellmem[nn][3]==1)
		{
		i=cellmem[nn][0];
		j=cellmem[nn][1];
		k=cellmem[nn][2];
		
		fbval = a->fb(i,j,k);
		
		nx=(a->fb(i+1,j,k)-a->fb(i-1,j,k))/(2.0*p->DXM);
		ny=(a->fb(i,j+1,k)-a->fb(i,j-1,k))/(2.0*p->DXM);
		nz=(a->fb(i,j,k+1)-a->fb(i,j,k-1))/(2.0*p->DXM);

		norm=sqrt(nx*nx + ny*ny + nz*nz);

		nx/=norm;
		ny/=norm;
		nz/=norm;

		posx = p->pos_x() + nx*(fabs(fbval) + p->DXM);
		posy = p->pos_y() + ny*(fabs(fbval) + p->DXM);
		posz = p->pos_z() + nz*(fabs(fbval) + p->DXM);
		
		
		fval = p->ccipol4_a(f,posx,posy,posz);
		
		ui = p->ufbi + (p->pos_z()-p->zg)*p->qfbi - (p->pos_y()-p->yg)*p->rfbi;
		vi = p->vfbi + (p->pos_x()-p->xg)*p->rfbi - (p->pos_z()-p->zg)*p->pfbi;
		wi = p->wfbi + (p->pos_y()-p->yg)*p->pfbi - (p->pos_x()-p->xg)*p->qfbi;
		
		un = p->ufbn + (p->pos_z()-p->zg)*p->qfbn - (p->pos_y()-p->yg)*p->rfbn;
		vn = p->vfbn + (p->pos_x()-p->xg)*p->rfbn - (p->pos_z()-p->zg)*p->pfbn;
		wn = p->wfbn + (p->pos_y()-p->yg)*p->pfbn - (p->pos_x()-p->xg)*p->qfbn;

		f(i,j,k)=fval;//nx*(ui-un)/p->dt + ny*(vi-vn)/p->dt + nz*(wi-wn)/p->dt;
		}
		
		/*
		if(cellmem[nn][3]==2)
		{
		i=cellmem[nn][0];
		j=cellmem[nn][1];
		k=cellmem[nn][2];
		
		fbval = a->fb(i,j,k);
		
		nx=(a->fb(i+1,j,k)-a->fb(i-1,j,k))/(2.0*p->DXM);
		ny=(a->fb(i,j+1,k)-a->fb(i,j-1,k))/(2.0*p->DXM);
		nz=(a->fb(i,j,k+1)-a->fb(i,j,k-1))/(2.0*p->DXM);

		norm=sqrt(nx*nx + ny*ny + nz*nz);

		nx/=norm;
		ny/=norm;
		nz/=norm;

		posx = p->pos_x() + nx*(fabs(fbval) + p->DXM);
		posy = p->pos_y() + ny*(fabs(fbval) + p->DXM);
		posz = p->pos_z() + nz*(fabs(fbval) + p->DXM);
		
		fval = ccipol4(a,f,posx,posy,posz);

		f(i,j,k)=fval;
		}*/
		
	}

}

/*
void ghostcell::gcfb_velupdate(lexer *p, fdm *a, int **cellmem, int cellcount,double xdelt, double ydelt, double zdelt, int dir)
{
    int nn;

    for(nn=0;nn<cellcount;++nn)
	{
		if(cellmem[nn][3]==1)
		{
		i=cellmem[nn][0];
		j=cellmem[nn][1];
		k=cellmem[nn][2];

		if(dir==1)
		a->u(i,j,k)=1.0*(p->ufbi + (p->pos_z()-p->zg)*p->qfbi - (p->pos_y()-p->yg)*p->rfbi);

		if(dir==2)
		a->v(i,j,k)=1.0*(p->vfbi + (p->pos_x()-p->xg)*p->rfbi - (p->pos_z()-p->zg)*p->pfbi);

		if(dir==3)
		a->w(i,j,k)=1.0*(p->wfbi + (p->pos_y()-p->yg)*p->pfbi - (p->pos_x()-p->xg)*p->qfbi);
		}
        

		
		if(cellmem[nn][3]==2)
		{
		i=cellmem[nn][0];
		j=cellmem[nn][1];
		k=cellmem[nn][2];
	
		//if(p->mpirank==0)
		//cout<<"CELLMEM: "<<i<<" "<<j<<" "<<k<<endl;

			if(dir==1)
			{
			if(p->flag1[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]>0) 
			a->u(i-1,j,k)=(p->ufbi + (p->pos_z()-p->zg)*p->qfbi - (p->pos_y()-p->yg)*p->rfbi);
			
			if(p->flag1[(i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]>0) 
			a->u(i+1,j,k)=(p->ufbi + (p->pos_z()-p->zg)*p->qfbi - (p->pos_y()-p->yg)*p->rfbi);
			}

			if(dir==2)
			{
			if(p->flag2[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin]>0) 
			a->v(i,j-1,k)=(p->vfbi + (p->pos_x()-p->xg)*p->rfbi - (p->pos_z()-p->zg)*p->pfbi);
			
			if(p->flag2[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin]>0) 
			a->v(i,j+1,k)=(p->vfbi + (p->pos_x()-p->xg)*p->rfbi - (p->pos_z()-p->zg)*p->pfbi);			
			}

			if(dir==3)
			{
			if(p->flag3[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1]>0) 
			a->w(i,j,k-1)=1.0*(p->wfbi + (p->pos_y()-p->yg)*p->pfbi - (p->pos_x()-p->xg)*p->qfbi);
			
			if(p->flag3[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+1]>0) 
			a->w(i,j,k+1)=1.0*(p->wfbi + (p->pos_y()-p->yg)*p->pfbi - (p->pos_x()-p->xg)*p->qfbi);
			}
		}
	}
}

void ghostcell::gcfb_scalarupdate(lexer *p, fdm *a, int **cellmem, int cellcount, field &f)
{
    int nn;
	double nx,ny,nz,norm;
	double posx,posy,posz;
	double locx,locy,locz;
	double fbval,fval;
	double ui,vi,wi,un,vn,wn;


    for(nn=0;nn<cellcount;++nn)
	{
		if(cellmem[nn][3]==1)
		{
		i=cellmem[nn][0];
		j=cellmem[nn][1];
		k=cellmem[nn][2];
		
		fbval = a->fb(i,j,k);
		
		nx=(a->fb(i+1,j,k)-a->fb(i-1,j,k))/(2.0*p->DXM);
		ny=(a->fb(i,j+1,k)-a->fb(i,j-1,k))/(2.0*p->DXM);
		nz=(a->fb(i,j,k+1)-a->fb(i,j,k-1))/(2.0*p->DXM);

		norm=sqrt(nx*nx + ny*ny + nz*nz);

		nx/=norm;
		ny/=norm;
		nz/=norm;

		posx = p->pos_x() + nx*(fabs(fbval) + p->DXM);
		posy = p->pos_y() + ny*(fabs(fbval) + p->DXM);
		posz = p->pos_z() + nz*(fabs(fbval) + p->DXM);
		
		
		fval = ccipol4_a(a,f,posx,posy,posz);
		
		ui = p->ufbi + (p->pos_z()-p->zg)*p->qfbi - (p->pos_y()-p->yg)*p->rfbi;
		vi = p->vfbi + (p->pos_x()-p->xg)*p->rfbi - (p->pos_z()-p->zg)*p->pfbi;
		wi = p->wfbi + (p->pos_y()-p->yg)*p->pfbi - (p->pos_x()-p->xg)*p->qfbi;
		
		un = p->ufbn + (p->pos_z()-p->zg)*p->qfbn - (p->pos_y()-p->yg)*p->rfbn;
		vn = p->vfbn + (p->pos_x()-p->xg)*p->rfbn - (p->pos_z()-p->zg)*p->pfbn;
		wn = p->wfbn + (p->pos_y()-p->yg)*p->pfbn - (p->pos_x()-p->xg)*p->qfbn;

		f(i,j,k)=fval;//nx*(ui-un)/p->dt + ny*(vi-vn)/p->dt + nz*(wi-wn)/p->dt;
		}
		
		
		if(cellmem[nn][3]==2)
		{
		i=cellmem[nn][0];
		j=cellmem[nn][1];
		k=cellmem[nn][2];
		
		fbval = a->fb(i,j,k);
		
		nx=(a->fb(i+1,j,k)-a->fb(i-1,j,k))/(2.0*p->DXM);
		ny=(a->fb(i,j+1,k)-a->fb(i,j-1,k))/(2.0*p->DXM);
		nz=(a->fb(i,j,k+1)-a->fb(i,j,k-1))/(2.0*p->DXM);

		norm=sqrt(nx*nx + ny*ny + nz*nz);

		nx/=norm;
		ny/=norm;
		nz/=norm;

		posx = p->pos_x() + nx*(fabs(fbval) + p->DXM);
		posy = p->pos_y() + ny*(fabs(fbval) + p->DXM);
		posz = p->pos_z() + nz*(fabs(fbval) + p->DXM);
		
		fval = ccipol4(a,f,posx,posy,posz);

		f(i,j,k)=fval;
		}
		
	}

}*/


