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

#include"ghostcell.h"
#include"lexer.h"
#include"field.h"

void ghostcell::dgcpol(lexer* p,field& f,int **dgc, int dgc_count, int gcv)
{/*
    double val0,val1,val2,val3,val4,val5,val6;
    int aa,bb,cc;
    int acheck,bcheck,ccheck;
    int a1,a2,b1,b2,c1,c2;

    val0=val1=val2=val3=val4=val5=val6=0.0;


    for(n=0;n<dgc_count;++n)
    {
    i=dgc[n][0];
    j=dgc[n][1];
    k=dgc[n][2];

    acheck=bcheck=ccheck=0;

    // pre-check
    for(q=0;q<dgc[n][3];++q)
    {
        if(dgc[n][4+q]==1)
        ++acheck;

        if(dgc[n][4+q]==4)
        ++acheck;

        if(dgc[n][4+q]==2)
        ++bcheck;

        if(dgc[n][4+q]==3)
        ++bcheck;

        if(dgc[n][4+q]==5)
        ++ccheck;

        if(dgc[n][4+q]==6)
        ++ccheck;
    }

// -------------------------------
    if(dgc[n][3]==2 && acheck!=2 && bcheck!=2 && ccheck!=2)
    {
        aa=bb=cc=0;
        val1=val2=val3=val4=val5=val6=0.0;

        for(q=0;q<2;++q)
        {
            if(dgc[n][4+q]==1)
            {
            val1=f(i-1,j,k);
            aa=-1;
            }

            if(dgc[n][4+q]==2)
            {
            val2=f(i,j+1,k);
            bb=1;
            }

            if(dgc[n][4+q]==3)
            {
            val3=f(i,j-1,k);
            bb=-1;
            }

            if(dgc[n][4+q]==4)
            {
            val4=f(i+1,j,k);
            aa=1;
            }

            if(dgc[n][4+q]==5)
            {
            val5=f(i,j,k-1);
            cc=-1;
            }

            if(dgc[n][4+q]==6)
            {
            val6=f(i,j,k+1);
            cc=1;
            }
        }
        f(i+aa,j+bb,k+cc)=0.5*(val1+val2+val3+val4+val5+val6);

    }

// -------------------------------
    if(dgc[n][3]==3 && acheck==1 && bcheck==1 && ccheck==1)
    {

         aa=bb=cc=0;
         val1=val2=val3=val4=val5=val6=0.0;
        //1
        for(q=0;q<3;++q)
        {
            if(dgc[n][4+q]==1)
            {
            val1=f(i-1,j,k);
            aa=-1;
            }

            if(dgc[n][4+q]==2)
            {
            val2=f(i,j+1,k);
            bb=1;
            }

            if(dgc[n][4+q]==3)
            {
            val3=f(i,j-1,k);
            bb=-1;
            }

            if(dgc[n][4+q]==4)
            {
            val4=f(i+1,j,k);
            aa=1;
            }

            if(dgc[n][4+q]==5)
            {
            val5=f(i,j,k-1);
            cc=-1;
            }

            if(dgc[n][4+q]==6)
            {
            val6=f(i,j,k+1);
            cc=1;
            }
        }

        f(i+aa,j+bb,k+cc)=0.333*(val1+val2+val3+val4+val5+val6);

        if(aa<0)
        {
            if(bb<0)
            {
                if(cc<0)
                {
                 f(i+aa,j+bb,k)=0.5*(val3+val1);
                 f(i+aa,j,k+cc)=0.5*(val1+val5);
                 f(i,j+bb,k+cc)=0.5*(val3+val5);
                }

                if(cc>0)
                {
                 f(i+aa,j+bb,k)=0.5*(val3+val1);
                 f(i+aa,j,k+cc)=0.5*(val1+val6);
                 f(i,j+bb,k+cc)=0.5*(val3+val6);
                }
            }

            if(bb>0)
            {
                if(cc<0)
                {
                 f(i+aa,j+bb,k)=0.5*(val2+val1);
                 f(i+aa,j,k+cc)=0.5*(val1+val5);
                 f(i,j+bb,k+cc)=0.5*(val2+val5);

                }

                if(cc>0)
                {
                 f(i+aa,j+bb,k)=0.5*(val2+val1);
                 f(i+aa,j,k+cc)=0.5*(val1+val6);
                 f(i,j+bb,k+cc)=0.5*(val2+val6);
                }
            }
        }


        if(aa>0)
        {
            if(bb<0)
            {
                if(cc<0)
                {
                 f(i+aa,j+bb,k)=0.5*(val3+val4);
                 f(i+aa,j,k+cc)=0.5*(val4+val5);
                 f(i,j+bb,k+cc)=0.5*(val3+val5);
                }

                if(cc>0)
                {
                 f(i+aa,j+bb,k)=0.5*(val3+val4);
                 f(i+aa,j,k+cc)=0.5*(val4+val6);
                 f(i,j+bb,k+cc)=0.5*(val3+val6);
                }
            }

            if(bb>0)
            {
                if(cc<0)
                {
                 f(i+aa,j+bb,k)=0.5*(val2+val4);
                 f(i+aa,j,k+cc)=0.5*(val4+val5);
                 f(i,j+bb,k+cc)=0.5*(val2+val5);

                }

                if(cc>0)
                {
                 f(i+aa,j+bb,k)=0.5*(val2+val4);
                 f(i+aa,j,k+cc)=0.5*(val4+val6);
                 f(i,j+bb,k+cc)=0.5*(val2+val6);
                }
            }
        }

    }

// -------------------------------
    if(dgc[n][3]==3 && (acheck==2 || bcheck==2 || ccheck==2))
    {
        a1=a2=b1=b2=c1=c2=0;
        val1=val2=val3=val4=val5=val6=0.0;

        //1
        for(q=0;q<3;++q)
        {
            if(dgc[n][4+q]==1)
            {
            val1=f(i-1,j,k);
            a1=-1;
            }

            if(dgc[n][4+q]==2)
            {
            val2=f(i,j+1,k);
            b2=1;
            }

            if(dgc[n][4+q]==3)
            {
            val3=f(i,j-1,k);
            b1=-1;
            }

            if(dgc[n][4+q]==4)
            {
            val4=f(i+1,j,k);
            a2=1;
            }

            if(dgc[n][4+q]==5)
            {
            val5=f(i,j,k-1);
            c1=-1;
            }

            if(dgc[n][4+q]==6)
            {
            val6=f(i,j,k+1);
            c2=1;
            }
        }


        if(acheck==2)
        {
            if(b1!=0)
            {
            f(i+a1,j+b1,k)=0.5*(val1+val3);
            f(i+a2,j+b1,k)=0.5*(val4+val3);
            }

            if(b2!=0)
            {
            f(i+a1,j+b2,k)=0.5*(val1+val2);
            f(i+a2,j+b2,k)=0.5*(val4+val2);
            }


            if(c1!=0)
            {
            f(i+a1,j,k+c1)=0.5*(val1+val5);
            f(i+a2,j,k+c1)=0.5*(val4+val5);
            }

            if(c2!=0)
            {
            f(i+a1,j,k+c2)=0.5*(val1+val6);
            f(i+a2,j,k+c2)=0.5*(val4+val6);
            }
        }


        if(bcheck==2)
        {
            if(a1!=0)
            {
            f(i+a1,j+b1,k)=0.5*(val3+val1);
            f(i+a1,j+b2,k)=0.5*(val2+val1);
            }

            if(a2!=0)
            {
            f(i+a2,j+b1,k)=0.5*(val3+val4);
            f(i+a2,j+b2,k)=0.5*(val2+val4);
            }


            if(c1!=0)
            {
            f(i,j+b1,k+c1)=0.5*(val3+val5);
            f(i,j+b2,k+c1)=0.5*(val2+val5);
            }

            if(c2!=0)
            {
            f(i,j+b1,k+c2)=0.5*(val3+val6);
            f(i,j+b2,k+c2)=0.5*(val2+val6);
            }
        }


        if(ccheck==2)
        {
            if(a1!=0)
            {
            f(i+a1,j,k+c1)=0.5*(val5+val1);
            f(i+a1,j,k+c2)=0.5*(val6+val1);
            }

            if(a2!=0)
            {
            f(i+a2,j,k+c1)=0.5*(val5+val4);
            f(i+a2,j,k+c2)=0.5*(val6+val4);
            }


            if(b1!=0)
            {
            f(i,j+b1,k+c1)=0.5*(val5+val3);
            f(i,j+b1,k+c2)=0.5*(val6+val3);
            }

            if(b2!=0)
            {
            f(i,j+b2,k+c1)=0.5*(val5+val2);
            f(i,j+b2,k+c2)=0.5*(val6+val2);
            }
        }
    }

// -------------------------------
    if(dgc[n][3]==4 && (acheck==2 || bcheck==2 || ccheck==2))
    {
        a1=a2=b1=b2=c1=c2=0;
        val1=val2=val3=val4=val5=val6=0.0;

        val1 = f(i,j,k);

        //1
        for(q=0;q<4;++q)
        {
            if(dgc[n][4+q]==1)
            a1=-1;

            if(dgc[n][4+q]==2)
            b2=1;

            if(dgc[n][4+q]==3)
            b1=-1;

            if(dgc[n][4+q]==4)
            a2=1;

            if(dgc[n][4+q]==5)
            c1=-1;

            if(dgc[n][4+q]==6)
            c2=1;
        }


        if(acheck==2)
        {
            if(b1!=0 && c1!=0) // 2 u 6
            {
            f(i,j+b1,k+c1)=val1; // 1

            f(i+a1,j+b1,k)=val1; // 3
            f(i+a2,j+b1,k)=val1; // 2

            f(i+a1,j+b1,k+c1)=val1; // 5
            f(i+a2,j+b1,k+c1)=val1; // 4

            f(i+a1,j,k+c1)=val1; // 7
            f(i+a2,j,k+c1)=val1; // 6
            }

            if(b1!=0 && c2!=0) // 2 u 6
            {
            f(i,j+b1,k+c2)=val1; // 1

            f(i+a1,j+b1,k)=val1; // 3
            f(i+a2,j+b1,k)=val1; // 2

            f(i+a1,j+b1,k+c2)=val1; // 5
            f(i+a2,j+b1,k+c2)=val1; // 4

            f(i+a1,j,k+c2)=val1; // 7
            f(i+a2,j,k+c2)=val1; // 6
            }

            if(b2!=0 && c1!=0) // 2 u 6
            {
            f(i,j+b2,k+c1)=val1; // 1

            f(i+a1,j+b2,k)=val1; // 3
            f(i+a2,j+b2,k)=val1; // 2

            f(i+a1,j+b2,k+c1)=val1; // 5
            f(i+a2,j+b2,k+c1)=val1; // 4

            f(i+a1,j,k+c1)=val1; // 7
            f(i+a2,j,k+c1)=val1; // 6
            }

            if(b2!=0 && c2!=0) // 2 u 6
            {
            f(i,j+b2,k+c2)=val1; // 1

            f(i+a1,j+b2,k)=val1; // 3
            f(i+a2,j+b2,k)=val1; // 2

            f(i+a1,j+b2,k+c2)=val1; // 5
            f(i+a2,j+b2,k+c2)=val1; // 4

            f(i+a1,j,k+c2)=val1; // 7
            f(i+a2,j,k+c2)=val1; // 6
            }
        }

        if(bcheck==2)
        {
            if(a1!=0 && c1!=0) // 2 u 6
            {
            f(i+a1,j,k+c1)=val1; // 1

            f(i+a1,j+b1,k)=val1; // 3
            f(i+a1,j+b2,k)=val1; // 2

            f(i+a1,j+b1,k+c1)=val1; // 5
            f(i+a1,j+b2,k+c1)=val1; // 4

            f(i,j+b1,k+c1)=val1; // 7
            f(i,j+b2,k+c1)=val1; // 6
            }

            if(a1!=0 && c2!=0) // 2 u 6
            {
            f(i+a1,j,k+c2)=val1; // 1

            f(i+a1,j+b1,k)=val1; // 3
            f(i+a1,j+b2,k)=val1; // 2

            f(i+a1,j+b1,k+c2)=val1; // 5
            f(i+a1,j+b2,k+c2)=val1; // 4

            f(i,j+b1,k+c2)=val1; // 7
            f(i,j+b2,k+c2)=val1; // 6
            }

            if(a2!=0 && c1!=0) // 2 u 6
            {
            f(i+a2,j,k+c1)=val1; // 1

            f(i+a2,j+b1,k)=val1; // 3
            f(i+a2,j+b2,k)=val1; // 2

            f(i+a2,j+b1,k+c1)=val1; // 5
            f(i+a2,j+b2,k+c1)=val1; // 4

            f(i,j+b1,k+c1)=val1; // 7
            f(i,j+b2,k+c1)=val1; // 6
            }

            if(a2!=0 && c2!=0) // 2 u 6
            {
            f(i+a2,j,k+c2)=val1; // 1

            f(i+a2,j+b1,k)=val1; // 3
            f(i+a2,j+b2,k)=val1; // 2

            f(i+a2,j+b1,k+c2)=val1; // 5
            f(i+a2,j+b2,k+c2)=val1; // 4

            f(i,j+b1,k+c2)=val1; // 7
            f(i,j+b2,k+c2)=val1; // 6
            }
        }


        if(ccheck==2)
        {
            if(a1!=0 && b1!=0) // 2 u 6
            {
            f(i+a1,j+b1,k)=val1; // 1

            f(i+a1,j,k+c1)=val1; // 3
            f(i+a1,j,k+c2)=val1; // 2

            f(i+a1,j+b1,k+c1)=val1; // 5
            f(i+a1,j+b1,k+c2)=val1; // 4

            f(i,j+b1,k+c1)=val1; // 7
            f(i,j+b1,k+c2)=val1; // 6
            }

            if(a1!=0 && b2!=0) // 2 u 6
            {
            f(i+a1,j+b2,k)=val1; // 1

            f(i+a1,j,k+c1)=val1; // 3
            f(i+a1,j,k+c2)=val1; // 2

            f(i+a1,j+b2,k+c1)=val1; // 5
            f(i+a1,j+b2,k+c2)=val1; // 4

            f(i,j+b2,k+c1)=val1; // 7
            f(i,j+b2,k+c2)=val1; // 6
            }

            if(a2!=0 && b1!=0) // 2 u 6
            {
            f(i+a2,j+b1,k)=val1; // 1

            f(i+a2,j,k+c1)=val1; // 3
            f(i+a2,j,k+c2)=val1; // 2

            f(i+a2,j+b1,k+c1)=val1; // 5
            f(i+a2,j+b1,k+c2)=val1; // 4

            f(i,j+b1,k+c1)=val1; // 7
            f(i,j+b1,k+c2)=val1; // 6
            }

            if(a2!=0 && b2!=0) // 2 u 6
            {
            f(i+a2,j+b2,k)=val1; // 1

            f(i+a2,j,k+c1)=val1; // 3
            f(i+a2,j,k+c2)=val1; // 2

            f(i+a2,j+b2,k+c1)=val1; // 5
            f(i+a2,j+b2,k+c2)=val1; // 4

            f(i,j+b2,k+c1)=val1; // 7
            f(i,j+b2,k+c2)=val1; // 6
            }
        }

    }


    }
*/
}
