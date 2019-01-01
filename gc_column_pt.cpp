/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"
#include"fieldint6.h"

void ghostcell::column_pt1_update(lexer* p, fdm* a)
{
    fieldint1 cval1(p);
    cvalupdate1(p,a,cval1);	
	column_pt1(p,a,cval1);
    cval_gcb1(p,a,cval1);
    cval_gcpara1(p,a,cval1);
}

void ghostcell::column_pt2_update(lexer* p, fdm* a)
{
    fieldint2 cval2(p);
    cvalupdate2(p,a,cval2);	
	column_pt2(p,a,cval2);
    cval_gcb2(p,a,cval2);
    cval_gcpara2(p,a,cval2);
}

void ghostcell::column_pt3_update(lexer* p, fdm* a)
{
    fieldint3 cval3(p);
    cvalupdate3(p,a,cval3);	
	column_pt3(p,a,cval3);
    cval_gcb3(p,a,cval3);
    cval_gcpara3(p,a,cval3);
}

void ghostcell::column_pt4_update(lexer* p, fdm* a)
{
    fieldint4 cval4(p);
    cvalupdate4(p,a,cval4);	
	column_pt4(p,a,cval4);
    cval_gcb4(p,a,cval4);
    cval_gcpara4(p,a,cval4);
}

void ghostcell::column_pt4a_update(lexer* p, fdm* a)
{
    fieldint4a cval4a(p);
    cvalupdate4a(p,a,cval4a);	
	column_pt4a(p,a,cval4a);
    cval_gcb4a(p,a,cval4a);
    cval_gcpara4a(p,a,cval4a);
}

void ghostcell::column_pt6_update(lexer* p, fdm* a)
{
    fieldint6 cval6(p);
    cvalupdate6(p,a,cval6);	
	column_pt6(p,a,cval6);
    cval_gcb6(p,a,cval6);
    cval_gcpara6(p,a,cval6);
}

void ghostcell::column_pt1(lexer* p, fdm* a, fieldint &cval1)
{
	n=0;
    ULOOP
	{	
	a->C1.p[n] = cval1(i,j,k);

	a->C1.n[n] = cval1(i+1,j,k);
	a->C1.s[n] = cval1(i-1,j,k);
	
	a->C1.w[n] = cval1(i,j+1,k);
	a->C1.e[n] = cval1(i,j-1,k);
	
	a->C1.t[n] = cval1(i,j,k+1);
	a->C1.b[n] = cval1(i,j,k-1);
	++n;
	}

	
	GGC1LOOP
    {
    i=p->gcb1[g][0];
    j=p->gcb1[g][1];
    k=p->gcb1[g][2];

        if(p->gcb1[g][3]==1)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C1.n[n] = cval1(i-q,j,k);
			a->C1.s[n] = cval1(i-q-2,j,k);
			}

        ++n;
        }

        if(p->gcb1[g][3]==2)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{	
			a->C1.e[n] = cval1(i,j+q,k);
			a->C1.w[n] = cval1(i,j+q+2,k);
			}
			
        ++n;
        }

        if(p->gcb1[g][3]==3)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{	
			a->C1.e[n] = cval1(i,j-q-2,k);
			a->C1.w[n] = cval1(i,j-q,k);
			}

        ++n;
        }

        if(p->gcb1[g][3]==4)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
            a->C1.s[n] = cval1(i+q,j,k);
			 a->C1.n[n] = cval1(i+q+2,j,k);
			
			}
		++n;
        }

        if(p->gcb1[g][3]==5)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C1.b[n] = cval1(i,j,k-q-2);
			a->C1.t[n] = cval1(i,j,k-q);
			}

		++n;
        }

        if(p->gcb1[g][3]==6)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C1.b[n] = cval1(i,j,k+q);
			a->C1.t[n] = cval1(i,j,k+q+2);
			}
		
		++n;
        }
    }


	for(g=0;g<p->gcpara1_count;++g)
    {
    i=p->gcpara1[g][0];
    j=p->gcpara1[g][1];
    k=p->gcpara1[g][2];
        
        if(p->gcpara1[g][3]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C1.n[n] = cval1(i-q,j,k); 
			a->C1.s[n] = cval1(i-q-2,j,k); 
			}
		++n;
		}
	}
	
	for(g=0;g<p->gcpara2_count;++g)
    {
    i=p->gcpara2[g][0];
    j=p->gcpara2[g][1];
    k=p->gcpara2[g][2];
        
        if(p->gcpara2[g][3]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C1.e[n] = cval1(i,j+q,k);
			a->C1.w[n] = cval1(i,j+q+2,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara3_count;++g)
    {
    i=p->gcpara3[g][0];
    j=p->gcpara3[g][1];
    k=p->gcpara3[g][2];
        
        if(p->gcpara3[g][3]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C1.e[n] = cval1(i,j-q-2,k);
			a->C1.w[n] = cval1(i,j-q,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara4_count;++g)
    {
    i=p->gcpara4[g][0];
    j=p->gcpara4[g][1];
    k=p->gcpara4[g][2];
        
        if(p->gcpara4[g][3]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C1.n[n] = cval1(i+q+2,j,k);
			a->C1.s[n] = cval1(i+q,j,k); 
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara5_count;++g)
    {
    i=p->gcpara5[g][0];
    j=p->gcpara5[g][1];
    k=p->gcpara5[g][2];
        
        if(p->gcpara5[g][3]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C1.b[n] = cval1(i,j,k-q-2);
			a->C1.t[n] = cval1(i,j,k-q);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara6_count;++g)
    {
    i=p->gcpara6[g][0];
    j=p->gcpara6[g][1];
    k=p->gcpara6[g][2];
        
        if(p->gcpara6[g][3]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C1.b[n] = cval1(i,j,k+q);
			a->C1.t[n] = cval1(i,j,k+q+2);
			}
		
		++n;
		}
	}
    
    //cout<<p->mpirank<<"Cpt1_n: "<<n<<"  p->veclength: "<<p->veclength<<endl;
}

void ghostcell::column_pt2(lexer* p, fdm* a, fieldint &cval2)
{
	n=0;
    VLOOP
	{	
	a->C2.p[n] = cval2(i,j,k);

	a->C2.n[n] = cval2(i+1,j,k);
	a->C2.s[n] = cval2(i-1,j,k);
	
	a->C2.w[n] = cval2(i,j+1,k);
	a->C2.e[n] = cval2(i,j-1,k);
	
	a->C2.t[n] = cval2(i,j,k+1);
	a->C2.b[n] = cval2(i,j,k-1);
	++n;
	}
	
	//int cellnum1=n;
	GGC2LOOP
    {
    i=p->gcb2[g][0];
    j=p->gcb2[g][1];
    k=p->gcb2[g][2];

        if(p->gcb2[g][3]==1)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C2.s[n] = cval2(i-q-2,j,k);
			a->C2.n[n] = cval2(i-q,j,k);
			}

        ++n;
        }

        if(p->gcb2[g][3]==2)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C2.e[n] = cval2(i,j+q,k);
			a->C2.w[n] = cval2(i,j+q+2,k);
			}
			
        ++n;
        }

        if(p->gcb2[g][3]==3)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C2.e[n] = cval2(i,j-q-2,k);
			a->C2.w[n] = cval2(i,j-q,k);
			}

        ++n;
        }

        if(p->gcb2[g][3]==4)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C2.s[n] = cval2(i+q,j,k);
			a->C2.n[n] = cval2(i+q+2,j,k);
			}

		++n;
        }

        if(p->gcb2[g][3]==5)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C2.b[n] = cval2(i,j,k-q-2);
			a->C2.t[n] = cval2(i,j,k-q);
			}

		++n;
        }

        if(p->gcb2[g][3]==6)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C2.b[n] = cval2(i,j,k+q);
			a->C2.t[n] = cval2(i,j,k+q+2);
			}
		
		++n;
        }
    }
	

	for(g=0;g<p->gcpara1_count;++g)
    {
    i=p->gcpara1[g][0];
    j=p->gcpara1[g][1];
    k=p->gcpara1[g][2];
        
        if(p->gcpara1[g][4]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C2.s[n] = cval2(i-q-2,j,k);
			a->C2.n[n] = cval2(i-q,j,k);
			}
		++n;
		}
	}
	
	for(g=0;g<p->gcpara2_count;++g)
    {
    i=p->gcpara2[g][0];
    j=p->gcpara2[g][1];
    k=p->gcpara2[g][2];
        
        if(p->gcpara2[g][4]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C2.e[n] = cval2(i,j+q,k);
			a->C2.w[n] = cval2(i,j+q+2,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara3_count;++g)
    {
    i=p->gcpara3[g][0];
    j=p->gcpara3[g][1];
    k=p->gcpara3[g][2];
        
        if(p->gcpara3[g][4]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C2.e[n] = cval2(i,j-q-2,k);
			a->C2.w[n] = cval2(i,j-q,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara4_count;++g)
    {
    i=p->gcpara4[g][0];
    j=p->gcpara4[g][1];
    k=p->gcpara4[g][2];
	
        if(p->gcpara4[g][4]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C2.s[n] = cval2(i+q,j,k);
			a->C2.n[n] = cval2(i+q+2,j,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara5_count;++g)
    {
    i=p->gcpara5[g][0];
    j=p->gcpara5[g][1];
    k=p->gcpara5[g][2];
        
        if(p->gcpara5[g][4]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C2.b[n] = cval2(i,j,k-q-2);
			a->C2.t[n] = cval2(i,j,k-q);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara6_count;++g)
    {
    i=p->gcpara6[g][0];
    j=p->gcpara6[g][1];
    k=p->gcpara6[g][2];
        
        if(p->gcpara6[g][4]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{	
			a->C2.b[n] = cval2(i,j,k+q);
			a->C2.t[n] = cval2(i,j,k+q+2);
			}
		
		++n;
		}
	}
	
}

void ghostcell::column_pt3(lexer* p, fdm* a, fieldint &cval3)
{
	n=0;
    WLOOP
	{	
	a->C3.p[n] = cval3(i,j,k);

	a->C3.n[n] = cval3(i+1,j,k);
	a->C3.s[n] = cval3(i-1,j,k);
	
	a->C3.w[n] = cval3(i,j+1,k);
	a->C3.e[n] = cval3(i,j-1,k);
	
	a->C3.t[n] = cval3(i,j,k+1);
	a->C3.b[n] = cval3(i,j,k-1);
	++n;
	}
	
	GGC3LOOP
    {
    i=p->gcb3[g][0];
    j=p->gcb3[g][1];
    k=p->gcb3[g][2];

        if(p->gcb3[g][3]==1)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C3.n[n] = cval3(i-q,j,k);
			a->C3.s[n] = cval3(i-q-2,j,k);
			}

        ++n;
        }

        if(p->gcb3[g][3]==2)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C3.e[n] = cval3(i,j+q,k);
			a->C3.w[n] = cval3(i,j+q+2,k);
			}
			
        ++n;
        }

        if(p->gcb3[g][3]==3)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C3.e[n] = cval3(i,j-q-2,k);
			a->C3.w[n] = cval3(i,j-q,k);
			}

        ++n;
        }

        if(p->gcb3[g][3]==4)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C3.n[n] = cval3(i+q+2,j,k);
			a->C3.s[n] = cval3(i+q,j,k);
			}

		++n;
        }

        if(p->gcb3[g][3]==5)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C3.t[n] = cval3(i,j,k-q);
			a->C3.b[n] = cval3(i,j,k-q-2);
			}

		++n;
        }

        if(p->gcb3[g][3]==6)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C3.t[n] = cval3(i,j,k+q+2);
			a->C3.b[n] = cval3(i,j,k+q);
			}
		
		++n;
        }
    }
	
	for(g=0;g<p->gcpara1_count;++g)
    {
    i=p->gcpara1[g][0];
    j=p->gcpara1[g][1];
    k=p->gcpara1[g][2];
        
        if(p->gcpara1[g][5]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C3.n[n] = cval3(i-q,j,k);
			a->C3.s[n] = cval3(i-q-2,j,k);
			}
 
		++n;
		}
	}
	
	for(g=0;g<p->gcpara2_count;++g)
    {
    i=p->gcpara2[g][0];
    j=p->gcpara2[g][1];
    k=p->gcpara2[g][2];
        
        if(p->gcpara2[g][5]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C3.e[n] = cval3(i,j+q,k);
			a->C3.w[n] = cval3(i,j+q+2,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara3_count;++g)
    {
    i=p->gcpara3[g][0];
    j=p->gcpara3[g][1];
    k=p->gcpara3[g][2];
        
        if(p->gcpara3[g][5]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C3.e[n] = cval3(i,j-q-2,k);
			a->C3.w[n] = cval3(i,j-q,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara4_count;++g)
    {
    i=p->gcpara4[g][0];
    j=p->gcpara4[g][1];
    k=p->gcpara4[g][2];
        
        if(p->gcpara4[g][5]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C3.n[n] = cval3(i+q+2,j,k);
			a->C3.s[n] = cval3(i+q,j,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara5_count;++g)
    {
    i=p->gcpara5[g][0];
    j=p->gcpara5[g][1];
    k=p->gcpara5[g][2];
        
        if(p->gcpara5[g][5]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C3.t[n] = cval3(i,j,k-q);
			a->C3.b[n] = cval3(i,j,k-q-2);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara6_count;++g)
    {
    i=p->gcpara6[g][0];
    j=p->gcpara6[g][1];
    k=p->gcpara6[g][2];
        
        if(p->gcpara6[g][5]==1)
		for(q=0;q<margin;++q)   
		{
			if(q<margin-1)
			{
			a->C3.t[n] = cval3(i,j,k+q+2);
			a->C3.b[n] = cval3(i,j,k+q);
			}
		
		++n;
		}
	}
}

void ghostcell::column_pt4(lexer* p, fdm* a, fieldint &cval4)
{
	n=0;
    LOOP
	{	
	a->C4.p[n] = cval4(i,j,k);
	
	a->C4.n[n] = cval4(i+1,j,k);
	a->C4.s[n] = cval4(i-1,j,k);
	
	a->C4.w[n] = cval4(i,j+1,k);
	a->C4.e[n] = cval4(i,j-1,k);
	
	a->C4.t[n] = cval4(i,j,k+1);
	a->C4.b[n] = cval4(i,j,k-1);
	++n;
	}
	
	//int cellnum1=n;
	GGC4LOOP
    {
    i=p->gcb4[g][0];
    j=p->gcb4[g][1];
    k=p->gcb4[g][2];

        if(p->gcb4[g][3]==1)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C4.n[n] = cval4(i-q,j,k);
			a->C4.s[n] = cval4(i-q-2,j,k);
			}

        ++n;
        }

        if(p->gcb4[g][3]==2)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C4.e[n] = cval4(i,j+q,k);
			a->C4.w[n] = cval4(i,j+q+2,k);
			}
			
        ++n;
        }

        if(p->gcb4[g][3]==3)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C4.e[n] = cval4(i,j-q-2,k); 
			a->C4.w[n] = cval4(i,j-q,k); 
			}

        ++n;
        }

        if(p->gcb4[g][3]==4)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C4.n[n] = cval4(i+q+2,j,k);
			a->C4.s[n] = cval4(i+q,j,k);
			}

		++n;
        }

        if(p->gcb4[g][3]==5)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C4.b[n] = cval4(i,j,k-q-2);
			a->C4.t[n] = cval4(i,j,k-q);
			}

		++n;
        }

        if(p->gcb4[g][3]==6)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C4.b[n] = cval4(i,j,k+q);
			a->C4.t[n] = cval4(i,j,k+q+2);
			}
		
		++n;
        }
    }
	
	//cout<<p->mpirank<<" gcbcount4: "<<p->gcb4_count*3<<" N4: "<<n-cellnum1<<endl;
	
	for(g=0;g<p->gcpara1_count;++g)
    {
    i=p->gcpara1[g][0];
    j=p->gcpara1[g][1];
    k=p->gcpara1[g][2];
        
        if(p->gcpara1[g][6]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C4.n[n] = cval4(i-q,j,k);
			a->C4.s[n] = cval4(i-q-2,j,k);
			}
 
		++n;
		}
	}
	
	for(g=0;g<p->gcpara2_count;++g)
    {
    i=p->gcpara2[g][0];
    j=p->gcpara2[g][1];
    k=p->gcpara2[g][2];
        
        if(p->gcpara2[g][6]==1)
		for(q=0;q<margin;++q)   
		{
			if(q<margin-1)
			{
			a->C4.e[n] = cval4(i,j+q,k);
			a->C4.w[n] = cval4(i,j+q+2,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara3_count;++g)
    {
    i=p->gcpara3[g][0];
    j=p->gcpara3[g][1];
    k=p->gcpara3[g][2];
        
        if(p->gcpara3[g][6]==1)
		for(q=0;q<margin;++q)   
		{
			if(q<margin-1)
			{
			a->C4.e[n] = cval4(i,j-q-2,k);
			a->C4.w[n] = cval4(i,j-q,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara4_count;++g)
    {
    i=p->gcpara4[g][0];
    j=p->gcpara4[g][1];
    k=p->gcpara4[g][2];
        
        if(p->gcpara4[g][6]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C4.n[n] = cval4(i+q+2,j,k);
			a->C4.s[n] = cval4(i+q,j,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara5_count;++g)
    {
    i=p->gcpara5[g][0];
    j=p->gcpara5[g][1];
    k=p->gcpara5[g][2];
        
        if(p->gcpara5[g][6]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C4.b[n] = cval4(i,j,k-q-2);
			a->C4.t[n] = cval4(i,j,k-q);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara6_count;++g)
    {
    i=p->gcpara6[g][0];
    j=p->gcpara6[g][1];
    k=p->gcpara6[g][2];
        
        if(p->gcpara6[g][6]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C4.b[n] = cval4(i,j,k+q);
			a->C4.t[n] = cval4(i,j,k+q+2);
			}
		
		++n;
		}
	}
}

void ghostcell::column_pt4a(lexer* p, fdm* a, fieldint &cval4a)
{
	n=0;
    ALOOP
	{	
	a->C4a.p[n] = cval4a(i,j,k);
	
	a->C4a.n[n] = cval4a(i+1,j,k);
	a->C4a.s[n] = cval4a(i-1,j,k);
	
	a->C4a.w[n] = cval4a(i,j+1,k);
	a->C4a.e[n] = cval4a(i,j-1,k);
	
	a->C4a.t[n] = cval4a(i,j,k+1);
	a->C4a.b[n] = cval4a(i,j,k-1);
	++n;
	}
	
	//int cellnum1=n;
	GGC4ALOOP
    {
    i=p->gcb4a[g][0];
    j=p->gcb4a[g][1];
    k=p->gcb4a[g][2];

        if(p->gcb4a[g][3]==1)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C4a.n[n] = cval4a(i-q,j,k);
			a->C4a.s[n] = cval4a(i-q-2,j,k);
			}

        ++n;
        }

        if(p->gcb4a[g][3]==2)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C4a.e[n] = cval4a(i,j+q,k);
			a->C4a.w[n] = cval4a(i,j+q+2,k);
			}
			
        ++n;
        }

        if(p->gcb4a[g][3]==3)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C4a.e[n] = cval4a(i,j-q-2,k);
			a->C4a.w[n] = cval4a(i,j-q,k);
			}

        ++n;
        }

        if(p->gcb4a[g][3]==4)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C4a.n[n] = cval4a(i+q+2,j,k);
			a->C4a.s[n] = cval4a(i+q,j,k);
			}

		++n;
        }

        if(p->gcb4a[g][3]==5)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C4a.b[n] = cval4a(i,j,k-q-2);
			a->C4a.t[n] = cval4a(i,j,k-q);
			}

		++n;
        }

        if(p->gcb4a[g][3]==6)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C4a.b[n] = cval4a(i,j,k+q);
			a->C4a.t[n] = cval4a(i,j,k+q+2);
			}
		
		++n;
        }
    }
	
	
	for(g=0;g<p->gcpara1_count;++g)
    {
    i=p->gcpara1[g][0];
    j=p->gcpara1[g][1];
    k=p->gcpara1[g][2];
	
        if(p->gcpara1[g][7]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C4a.n[n] = cval4a(i-q,j,k);
			a->C4a.s[n] = cval4a(i-q-2,j,k);
			}
 
		++n;
		}
	}
	
	for(g=0;g<p->gcpara2_count;++g)
    {
    i=p->gcpara2[g][0];
    j=p->gcpara2[g][1];
    k=p->gcpara2[g][2];
        
        if(p->gcpara2[g][7]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C4a.e[n] = cval4a(i,j+q,k);
			a->C4a.w[n] = cval4a(i,j+q+2,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara3_count;++g)
    {
    i=p->gcpara3[g][0];
    j=p->gcpara3[g][1];
    k=p->gcpara3[g][2];
        
        if(p->gcpara3[g][7]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C4a.e[n] = cval4a(i,j-q-2,k);
			a->C4a.w[n] = cval4a(i,j-q,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara4_count;++g)
    {
    i=p->gcpara4[g][0];
    j=p->gcpara4[g][1];
    k=p->gcpara4[g][2];
        
        if(p->gcpara4[g][7]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C4a.n[n] = cval4a(i+q+2,j,k);
			a->C4a.s[n] = cval4a(i+q,j,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara5_count;++g)
    {
    i=p->gcpara5[g][0];
    j=p->gcpara5[g][1];
    k=p->gcpara5[g][2];
        
        if(p->gcpara5[g][7]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C4a.b[n] = cval4a(i,j,k-q-2);
			a->C4a.t[n] = cval4a(i,j,k-q);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara6_count;++g)
    {
    i=p->gcpara6[g][0];
    j=p->gcpara6[g][1];
    k=p->gcpara6[g][2];
        
        if(p->gcpara6[g][7]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C4a.b[n] = cval4a(i,j,k+q);
			a->C4a.t[n] = cval4a(i,j,k+q+2);
			}
		
		++n;
		}
	}
}

void ghostcell::column_pt6(lexer* p, fdm* a, fieldint &cval6)
{
	n=0;
    BASELOOP
	{	
	a->C6.p[n] = cval6(i,j,k);
	
	a->C6.n[n] = cval6(i+1,j,k);
	a->C6.s[n] = cval6(i-1,j,k);
	
	a->C6.w[n] = cval6(i,j+1,k);
	a->C6.e[n] = cval6(i,j-1,k);
	
	a->C6.t[n] = cval6(i,j,k+1);
	a->C6.b[n] = cval6(i,j,k-1);
	++n;
	}
  
	//int cellnum1=n;
	GGC6LOOP
    {
    i=p->gcb4[g][0];
    j=p->gcb4[g][1];
    k=p->gcb4[g][2];

        if(p->gcb4[g][3]==1)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C6.n[n] = cval6(i-q,j,k);
			a->C6.s[n] = cval6(i-q-2,j,k);
			}

        ++n;
        }

        if(p->gcb4[g][3]==2)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C6.e[n] = cval6(i,j+q,k);
			a->C6.w[n] = cval6(i,j+q+2,k);
			}
			
        ++n;
        }

        if(p->gcb4[g][3]==3)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C6.e[n] = cval6(i,j-q-2,k); 
			a->C6.w[n] = cval6(i,j-q,k); 
			}

        ++n;
        }

        if(p->gcb4[g][3]==4)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C6.n[n] = cval6(i+q+2,j,k);
			a->C6.s[n] = cval6(i+q,j,k);
			}

		++n;
        }

        if(p->gcb4[g][3]==5)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C6.b[n] = cval6(i,j,k-q-2);
			a->C6.t[n] = cval6(i,j,k-q);
			}

		++n;
        }

        if(p->gcb4[g][3]==6)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C6.b[n] = cval6(i,j,k+q);
			a->C6.t[n] = cval6(i,j,k+q+2);
			}
		
		++n;
        }
    }
	


	for(g=0;g<p->gcpara1_count;++g)
    {
    i=p->gcpara1[g][0];
    j=p->gcpara1[g][1];
    k=p->gcpara1[g][2];
        
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C6.n[n] = cval6(i-q,j,k);
			a->C6.s[n] = cval6(i-q-2,j,k);
			}
 
		++n;
		}
	}
	
	for(g=0;g<p->gcpara2_count;++g)
    {
    i=p->gcpara2[g][0];
    j=p->gcpara2[g][1];
    k=p->gcpara2[g][2];
        
		for(q=0;q<margin;++q)   
		{
			if(q<margin-1)
			{
			a->C6.e[n] = cval6(i,j+q,k);
			a->C6.w[n] = cval6(i,j+q+2,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara3_count;++g)
    {
    i=p->gcpara3[g][0];
    j=p->gcpara3[g][1];
    k=p->gcpara3[g][2];
        
		for(q=0;q<margin;++q)   
		{
			if(q<margin-1)
			{
			a->C6.e[n] = cval6(i,j-q-2,k);
			a->C6.w[n] = cval6(i,j-q,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara4_count;++g)
    {
    i=p->gcpara4[g][0];
    j=p->gcpara4[g][1];
    k=p->gcpara4[g][2];
        
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C6.n[n] = cval6(i+q+2,j,k);
			a->C6.s[n] = cval6(i+q,j,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara5_count;++g)
    {
    i=p->gcpara5[g][0];
    j=p->gcpara5[g][1];
    k=p->gcpara5[g][2];
        
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C6.b[n] = cval6(i,j,k-q-2);
			a->C6.t[n] = cval6(i,j,k-q);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara6_count;++g)
    {
    i=p->gcpara6[g][0];
    j=p->gcpara6[g][1];
    k=p->gcpara6[g][2];
        
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C6.b[n] = cval6(i,j,k+q);
			a->C6.t[n] = cval6(i,j,k+q+2);
			}
		
		++n;
		}
	}
}