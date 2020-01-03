/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

#include"fillvec.h"
#include"lexer.h"
#include"fdm.h"

fillvec::fillvec()
{
	margin=3;
}

fillvec::~fillvec()
{
}

void fillvec::fillxvec1(lexer* p, fdm* a, field& f)
{
	count=0;
    ULOOP
    {
    a->xvec.V[count]=f(i,j,k);
    ++count;
    }

    GC1LOOP
    {
    i=p->gcb1[n][0];
    j=p->gcb1[n][1];
    k=p->gcb1[n][2];

        if(p->gcb1[n][3]==1)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i-1-q,j,k);
        ++count;
        }

        if(p->gcb1[n][3]==2)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j+1+q,k);
        ++count;
        }

        if(p->gcb1[n][3]==3)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j-1-q,k);
        ++count;
        }

        if(p->gcb1[n][3]==4)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i+1+q,j,k);
        ++count;
        }

        if(p->gcb1[n][3]==5)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j,k-1-q);
        ++count;
        }

        if(p->gcb1[n][3]==6)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j,k+1+q);
        ++count;
        }
    }
	
	// parallel boundary
    for(q=0;q<p->gcpara1_count;++q)
    {
    i=p->gcpara1[q][0];
    j=p->gcpara1[q][1];
    k=p->gcpara1[q][2];
    
    n=p->gcpara1[q][9];
        
        if(p->gcpara1[q][3]==1)
        {
        a->xvec.V[Im1_J_K_1]=f(i-1,j,k);
        a->xvec.V[Im2_J_K_1]=f(i-2,j,k);
        a->xvec.V[Im3_J_K_1]=f(i-3,j,k);
        }
    }

	for(q=0;q<p->gcpara2_count;++q)
	{
    i=p->gcpara2[q][0];
    j=p->gcpara2[q][1];
    k=p->gcpara2[q][2];
        
    n=p->gcpara2[q][9];
        
        if(p->gcpara2[q][3]==1)
        {
        a->xvec.V[I_Jp1_K_1]=f(i,j+1,k);
        a->xvec.V[I_Jp2_K_1]=f(i,j+2,k);
        a->xvec.V[I_Jp3_K_1]=f(i,j+3,k);
        }
	}	
	
	for(q=0;q<p->gcpara3_count;++q)
	{
    i=p->gcpara3[q][0];
    j=p->gcpara3[q][1];
    k=p->gcpara3[q][2];
        
    n=p->gcpara3[q][9];
        
        if(p->gcpara3[q][3]==1)
        {
        a->xvec.V[I_Jm1_K_1]=f(i,j,k);
        a->xvec.V[I_Jm2_K_1]=f(i,j,k);
        a->xvec.V[I_Jm3_K_1]=f(i,j,k);
        }
	}

	for(q=0;q<p->gcpara4_count;++q)
	{
    i=p->gcpara4[q][0];
    j=p->gcpara4[q][1];
    k=p->gcpara4[q][2];
        
    n=p->gcpara4[q][9];
        
        if(p->gcpara4[q][3]==1)
        {
        a->xvec.V[Ip1_J_K_1]=f(i+1,j,k);
        a->xvec.V[Ip2_J_K_1]=f(i+2,j,k);
        a->xvec.V[Ip3_J_K_1]=f(i+3,j,k);
        }
	}

    for(q=0;q<p->gcpara5_count;++q)
    {
    i=p->gcpara5[q][0];
    j=p->gcpara5[q][1];
    k=p->gcpara5[q][2];
        
    n=p->gcpara5[q][9];
        
        if(p->gcpara5[q][3]==1)
        {
        a->xvec.V[I_J_Km1_1]=f(i,j,k-1);
        a->xvec.V[I_J_Km2_1]=f(i,j,k-2);
        a->xvec.V[I_J_Km3_1]=f(i,j,k-3);
        }
    }

	for(q=0;q<p->gcpara6_count;++q)
	{
    i=p->gcpara6[q][0];
    j=p->gcpara6[q][1];
    k=p->gcpara6[q][2];
        
    n=p->gcpara6[q][9];
        
        if(p->gcpara6[q][3]==1)
        {
        a->xvec.V[I_J_Kp1_1]=f(i,j,k+1);
        a->xvec.V[I_J_Kp2_1]=f(i,j,k+2);
        a->xvec.V[I_J_Kp3_1]=f(i,j,k+3);
        }
	}
}

void fillvec::fillxvec2(lexer* p, fdm* a, field& f)
{
	count=0;
    VLOOP
    {
    a->xvec.V[count]=f(i,j,k);
    ++count;
    }

    GC2LOOP
    {
    i=p->gcb2[n][0];
    j=p->gcb2[n][1];
    k=p->gcb2[n][2];

        if(p->gcb2[n][3]==1)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i-1-q,j,k);
        ++count;
        }

        if(p->gcb2[n][3]==2)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j+1+q,k);
        ++count;
        }

        if(p->gcb2[n][3]==3)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j-1-q,k);
        ++count;
        }

        if(p->gcb2[n][3]==4)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i+1+q,j,k);
        ++count;
        }

        if(p->gcb2[n][3]==5)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j,k-1-q);
        ++count;
        }

        if(p->gcb2[n][3]==6)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j,k+1+q);
        ++count;
        }
    }
	
	// parallel boundary
    for(q=0;q<p->gcpara1_count;++q)
    {
    i=p->gcpara1[q][0];
    j=p->gcpara1[q][1];
    k=p->gcpara1[q][2];
        
    n=p->gcpara1[q][9];
        
        if(p->gcpara1[q][3]==1)
        {
        a->xvec.V[Im1_J_K_2]=f(i-1,j,k);
        a->xvec.V[Im2_J_K_2]=f(i-2,j,k);
        a->xvec.V[Im3_J_K_2]=f(i-3,j,k);
        }
    }

	for(q=0;q<p->gcpara2_count;++q)
	{
    i=p->gcpara2[q][0];
    j=p->gcpara2[q][1];
    k=p->gcpara2[q][2];
        
    n=p->gcpara2[q][10];
        
        if(p->gcpara2[q][3]==1)
        {
        a->xvec.V[I_Jp1_K_2]=f(i,j+1,k);
        a->xvec.V[I_Jp2_K_2]=f(i,j+2,k);
        a->xvec.V[I_Jp3_K_2]=f(i,j+3,k);
        }
	}	
	
	for(q=0;q<p->gcpara3_count;++q)
	{
    i=p->gcpara3[q][0];
    j=p->gcpara3[q][1];
    k=p->gcpara3[q][2];
        
    n=p->gcpara3[q][10];
        
        if(p->gcpara3[q][3]==1)
        {
        a->xvec.V[I_Jm1_K_2]=f(i,j-1,k);
        a->xvec.V[I_Jm2_K_2]=f(i,j-2,k);
        a->xvec.V[I_Jm3_K_2]=f(i,j-3,k);
        }
	}

	for(q=0;q<p->gcpara4_count;++q)
	{
    i=p->gcpara4[q][0];
    j=p->gcpara4[q][1];
    k=p->gcpara4[q][2];
        
    n=p->gcpara4[q][10];
        
        if(p->gcpara4[q][3]==1)
        {
        a->xvec.V[I_J_K_2]=f(i,j,k);
        a->xvec.V[I_J_K_2]=f(i,j,k);
        a->xvec.V[I_J_K_2]=f(i,j,k);
        }
	}

    for(q=0;q<p->gcpara5_count;++q)
    {
    i=p->gcpara5[q][0];
    j=p->gcpara5[q][1];
    k=p->gcpara5[q][2];
        
    n=p->gcpara5[q][10];
        
        if(p->gcpara5[q][3]==1)
        {
        a->xvec.V[I_J_Km1_2]=f(i,j,k-1);
        a->xvec.V[I_J_Km2_2]=f(i,j,k-2);
        a->xvec.V[I_J_Km3_2]=f(i,j,k-3);
        }
    }

	for(q=0;q<p->gcpara6_count;++q)
	{
    i=p->gcpara6[q][0];
    j=p->gcpara6[q][1];
    k=p->gcpara6[q][2];
        
    n=p->gcpara6[q][10];
        
        if(p->gcpara6[q][3]==1)
        {
        a->xvec.V[I_J_Kp1_2]=f(i,j,k+1);
        a->xvec.V[I_J_Kp2_2]=f(i,j,k+2);
        a->xvec.V[I_J_Kp3_2]=f(i,j,k+3);
        }
	}
}

void fillvec::fillxvec3(lexer* p, fdm* a, field& f)
{
	count=0;
    WLOOP
    {
    a->xvec.V[count]=f(i,j,k);
    ++count;
    }

    GC3LOOP
    {
    i=p->gcb3[n][0];
    j=p->gcb3[n][1];
    k=p->gcb3[n][2];

        if(p->gcb3[n][3]==1)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i-1-q,j,k);
        ++count;
        }

        if(p->gcb3[n][3]==2)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j+1+q,k);
        ++count;
        }

        if(p->gcb3[n][3]==3)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j-1-q,k);
        ++count;
        }

        if(p->gcb3[n][3]==4)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i+1+q,j,k);
        ++count;
        }

        if(p->gcb3[n][3]==5)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j,k-1-q);
        ++count;
        }

        if(p->gcb3[n][3]==6)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j,k+1+q);
        ++count;
        }
    }
	
	// parallel boundary
    for(q=0;q<p->gcpara1_count;++q)
    {
    i=p->gcpara1[q][0];
    j=p->gcpara1[q][1];
    k=p->gcpara1[q][2];
        
    n=p->gcpara1[q][11];
        
        if(p->gcpara1[q][3]==1)
        {
        a->xvec.V[Im1_J_K_3]=f(i-1,j,k);
        a->xvec.V[Im2_J_K_3]=f(i-2,j,k);
        a->xvec.V[Im3_J_K_3]=f(i-3,j,k);
        }
    }

	for(q=0;q<p->gcpara2_count;++q)
	{
    i=p->gcpara2[q][0];
    j=p->gcpara2[q][1];
    k=p->gcpara2[q][2];
        
    n=p->gcpara2[q][11];
        
        if(p->gcpara2[q][3]==1)
        {
        a->xvec.V[I_Jp1_K_3]=f(i,j+1,k);
        a->xvec.V[I_Jp2_K_3]=f(i,j+2,k);
        a->xvec.V[I_Jp3_K_3]=f(i,j+3,k);
        }
	}	
	
	for(q=0;q<p->gcpara3_count;++q)
	{
    i=p->gcpara3[q][0];
    j=p->gcpara3[q][1];
    k=p->gcpara3[q][2];
        
    n=p->gcpara3[q][11];
        
        if(p->gcpara3[q][3]==1)
        {
        a->xvec.V[I_Jp1_K_3]=f(i,j+1,k);
        a->xvec.V[I_Jp2_K_3]=f(i,j+2,k);
        a->xvec.V[I_Jp3_K_3]=f(i,j+3,k);
        }
	}

	for(q=0;q<p->gcpara4_count;++q)
	{
    i=p->gcpara4[q][0];
    j=p->gcpara4[q][1];
    k=p->gcpara4[q][2];
        
    n=p->gcpara4[q][11];
        
        if(p->gcpara4[q][3]==1)
        {
        a->xvec.V[Ip1_J_K_3]=f(i+1,j,k);
        a->xvec.V[Ip2_J_K_3]=f(i+2,j,k);
        a->xvec.V[Ip3_J_K_3]=f(i+3,j,k);
        }
	}

    for(q=0;q<p->gcpara5_count;++q)
    {
    i=p->gcpara5[q][0];
    j=p->gcpara5[q][1];
    k=p->gcpara5[q][2];
        
    n=p->gcpara5[q][11];
        
        if(p->gcpara5[q][3]==1)
        {
        a->xvec.V[I_J_Km1_3]=f(i,j,k-1);
        a->xvec.V[I_J_Km2_3]=f(i,j,k-2);
        a->xvec.V[I_J_Km3_3]=f(i,j,k-3);
        }
    }

	for(q=0;q<p->gcpara6_count;++q)
	{
    i=p->gcpara6[q][0];
    j=p->gcpara6[q][1];
    k=p->gcpara6[q][2];
        
    n=p->gcpara6[q][11];
        
        if(p->gcpara6[q][3]==1)
        {
        a->xvec.V[I_J_Kp1_3]=f(i,j,k+1);
        a->xvec.V[I_J_Kp2_3]=f(i,j,k+2);
        a->xvec.V[I_J_Kp3_3]=f(i,j,k+3);
        }
	}
}

void fillvec::fillxvec4(lexer* p, fdm* a, field& f)
{	
	count=0;
    LOOP
    {
    a->xvec.V[count]=f(i,j,k);
    ++count;
    }

    GC4LOOP
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];

        if(p->gcb4[n][3]==1)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i-1-q,j,k);
        ++count;
        }

        if(p->gcb4[n][3]==2)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j+1+q,k);
        ++count;
        }

        if(p->gcb4[n][3]==3)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j-1-q,k);
        ++count;
        }

        if(p->gcb4[n][3]==4)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i+1+q,j,k);
        ++count;
        }

        if(p->gcb4[n][3]==5)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j,k-1-q);
        ++count;
        }

        if(p->gcb4[n][3]==6)
        for(q=0;q<margin;++q)
        {
        a->xvec.V[count]=f(i,j,k+1+q);
        ++count;
        }
    }
	
	// parallel boundary
    for(q=0;q<p->gcpara1_count;++q)
    {
    i=p->gcpara1[q][0];
    j=p->gcpara1[q][1];
    k=p->gcpara1[q][2];
        
    n=p->gcpara1[q][12];
        
        if(p->gcpara1[q][3]==1)
        {
        a->xvec.V[Im1_J_K_4]=f(i-1,j,k);
        a->xvec.V[Im2_J_K_4]=f(i-2,j,k);
        a->xvec.V[Im3_J_K_4]=f(i-3,j,k);
        }
    }

	for(q=0;q<p->gcpara2_count;++q)
	{
    i=p->gcpara2[q][0];
    j=p->gcpara2[q][1];
    k=p->gcpara2[q][2];
        
    n=p->gcpara2[q][12];
        
        if(p->gcpara2[q][3]==1)
        {
        a->xvec.V[I_Jp1_K_4]=f(i,j+1,k);
        a->xvec.V[I_Jp2_K_4]=f(i,j+2,k);
        a->xvec.V[I_Jp3_K_4]=f(i,j+3,k);
        }
	}	
	
	for(q=0;q<p->gcpara3_count;++q)
	{
    i=p->gcpara3[q][0];
    j=p->gcpara3[q][1];
    k=p->gcpara3[q][2];
        
    n=p->gcpara3[q][12];
        
        if(p->gcpara3[q][3]==1)
        {
        a->xvec.V[I_Jm1_K_4]=f(i,j-1,k);
        a->xvec.V[I_Jm2_K_4]=f(i,j-2,k);
        a->xvec.V[I_Jm3_K_4]=f(i,j-3,k);
        }
	}

	for(q=0;q<p->gcpara4_count;++q)
	{
    i=p->gcpara4[q][0];
    j=p->gcpara4[q][1];
    k=p->gcpara4[q][2];
    
    n=p->gcpara4[q][12];
        
        if(p->gcpara4[q][3]==1)
        for(n=0;n<margin;++n)
        {
        a->xvec.V[count]=f(i+n+1,j,k);
        ++count;
        }
	}

    for(q=0;q<p->gcpara5_count;++q)
    {
    i=p->gcpara5[q][0];
    j=p->gcpara5[q][1];
    k=p->gcpara5[q][2];
        
    n=p->gcpara5[q][12];
        
        if(p->gcpara5[q][3]==1)
        {
        a->xvec.V[I_J_Km1_4]=f(i,j,k-1);
        a->xvec.V[I_J_Km2_4]=f(i,j,k-2);
        a->xvec.V[I_J_Km3_4]=f(i,j,k-3);
        }
    }

	for(q=0;q<p->gcpara6_count;++q)
	{
    i=p->gcpara6[q][0];
    j=p->gcpara6[q][1];
    k=p->gcpara6[q][2];
        
    n=p->gcpara6[q][12];
        
        if(p->gcpara6[q][3]==1)
        {
        a->xvec.V[I_J_Kp1_4]=f(i,j,k+1);
        a->xvec.V[I_J_Kp2_4]=f(i,j,k+2);
        a->xvec.V[I_J_Kp3_4]=f(i,j,k+3);
        }
	}

}

void fillvec::fillxfield1(lexer *p, fdm *a, field &xfield)
{
        count=0;
        ULOOP
        {
        xfield(i,j,k)=a->xvec.V[count];
        ++count;
        }
}

void fillvec::fillxfield2(lexer *p, fdm *a, field &xfield)
{
        count=0;
        VLOOP
        {
        xfield(i,j,k)=a->xvec.V[count];
        ++count;
		}
}

void fillvec::fillxfield3(lexer *p, fdm *a, field &xfield)
{	
        count=0;
        WLOOP
        {
        xfield(i,j,k)=a->xvec.V[count];
        ++count;
        }
}

void fillvec::fillxfield4(lexer *p, fdm *a, field &xfield)
{
        count=0;
        LOOP
        {
        xfield(i,j,k)=a->xvec.V[count];
        ++count;
        }
}
