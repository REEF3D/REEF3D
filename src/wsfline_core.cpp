/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"wsfline_core.h"
#include"lexer.h"
#include"ghostcell.h"
#include<iomanip>

wsfline_core::wsfline_core(lexer *p, ghostcell *pgc, int nl, int nloc) : nline(nl)
{
    maxn=pgc->globalimax(nloc);
    sumn=pgc->globalisum(maxn);

    p->Darray(loc,nline+1,maxn);
    p->Darray(wsf,nline+1,maxn);
    p->Iarray(flag,nline+1,maxn);
	p->Iarray(wsfpoints,nline+1);

    p->Darray(loc_all,nline+1,sumn);
    p->Darray(wsf_all,nline+1,sumn);
	p->Iarray(flag_all,nline+1,sumn);
	p->Iarray(rowflag,sumn);
}

wsfline_core::~wsfline_core()
{
}

void wsfline_core::reset()
{
    for(int q=0;q<nline;++q)
    for(int n=0;n<maxn;++n)
    {
    loc[q][n]=1.0e20;
    wsf[q][n]=-1.0e20;
    }
}

void wsfline_core::write_rows(lexer *p, ghostcell *pgc, ofstream &wsfout, int prec, const function<double(double)> &coord,
                              const function<double(double)> &theory, const char *fill)
{
    int check;

	for(int q=0;q<nline;++q)
    wsfpoints[q]=sumn;

    // gather
    for(int q=0;q<nline;++q)
    {
    pgc->gather_double(loc[q],maxn,loc_all[q],maxn);
    pgc->gather_double(wsf[q],maxn,wsf_all[q],maxn);
	pgc->gather_int(flag[q],maxn,flag_all[q],maxn);

        if(p->mpirank==0)
        {
        sort(loc_all[q], wsf_all[q], flag_all[q], 0, wsfpoints[q]-1);
        remove_multientry(p,loc_all[q], wsf_all[q], flag_all[q], wsfpoints[q]);
        }
    }

    // write to file
    if(p->mpirank==0)
    {
		for(int n=0;n<sumn;++n)
		rowflag[n]=0;

		for(int n=0;n<sumn;++n)
        {
			check=0;
		    for(int q=0;q<nline;++q)
			if(flag_all[q][n]>0 && loc_all[q][n]<1.0e20)
			check=1;

			if(check==1)
			rowflag[n]=1;
		}

        for(int n=0;n<sumn;++n)
        {
			check=0;
		    for(int q=0;q<nline;++q)
			{
				if(flag_all[q][n]>0 && loc_all[q][n]<1.0e20)
				{
				wsfout<<setprecision(prec)<<coord(loc_all[q][n])<<" \t ";
				wsfout<<setprecision(prec)<<wsf_all[q][n]<<" \t  ";

					if(theory)
					wsfout<<theory(loc_all[q][n])<<" \t  ";

				check=1;
				}

				if((flag_all[q][n]<0 || loc_all[q][n]>=1.0e20) && rowflag[n]==1)
				{
					wsfout<<setprecision(5)<<fill;
					wsfout<<setprecision(5)<<fill;
				}
			}

			if(check==1)
            wsfout<<endl;
        }
    }
}

void wsfline_core::sort(double *a, double *b, int *c, int left, int right)
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

void wsfline_core::remove_multientry(lexer *p, double* b, double* c, int *d, int& num)
{
    int oldnum=num;
    double xval=-1.12e23;

    int count=0;

    double *f,*g;
	int *h;

	p->Darray(f,num);
	p->Darray(g,num);
	p->Iarray(h,num);

    for(int n=0;n<num;++n)
    g[n]=-1.12e22;

    for(int n=0;n<oldnum;++n)
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

    for(int n=0;n<count;++n)
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
