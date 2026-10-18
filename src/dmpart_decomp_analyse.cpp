/*--------------------------------------------------------------------
DIVEMesh
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D (port of the DIVEMesh decomposition, dmpart).

DIVEMesh is free software; you can redistribute it and/or modify it
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

#include "dmpart_decomp.h"

namespace dmpart
{

void decomp::partition_analyse(lexer *p, dive *a)
{
    for(i=0;i<a->knox;++i)
    xcross[i]=0;

    for(j=0;j<a->knoy;++j)
    ycross[j]=0;

    for(k=0;k<a->knoz;++k)
    zcross[k]=0;

    LOOP
    if(a->flag(i,j,k)>0 && a->solid(i,j,k)>0)
    {
        xcross[i] += 1;
        ycross[j] += 1;
        zcross[k] += 1;
    }

    xcross_m=0;
    ycross_m=0;
    zcross_m=0;

    for(i=0;i<a->knox;++i)
    xcross_m+=xcross[i];

    for(j=0;j<a->knoy;++j)
    ycross_m+=ycross[j];

    for(k=0;k<a->knoz;++k)
    zcross_m+=zcross[k];

    xcross_m=xcross_m/a->knox;
    ycross_m=ycross_m/a->knoy;
    zcross_m=zcross_m/a->knoz;

    ddout<<"xcross_m: "<<xcross_m<<endl;
    ddout<<"ycross_m: "<<ycross_m<<endl;
    ddout<<"zcross_m: "<<zcross_m<<endl;

    xvar=yvar=zvar=0.0;

    for(i=0;i<a->knox;++i)
    xvar += pow(xcross[i]-xcross_m,2.0);

    for(j=0;j<a->knoy;++j)
    yvar += pow(ycross[j]-ycross_m,2.0);

    for(k=0;k<a->knoz;++k)
    zvar += pow(zcross[k]-zcross_m,2.0);

    xvar=sqrt(xvar/a->knox);
    yvar=sqrt(yvar/a->knoy);
    zvar=sqrt(zvar/a->knoz);

    xvar_per = xvar/xcross_m;
    yvar_per = yvar/ycross_m;
    zvar_per = zvar/zcross_m;

    ddout<<"xvar: "<<xvar<<endl;
    ddout<<"yvar: "<<yvar<<endl;
    ddout<<"zvar: "<<zvar<<endl;

    ddout<<"xvar_per: "<<xvar_per<<endl;
    ddout<<"yvar_per: "<<yvar_per<<endl;
    ddout<<"zvar_per: "<<zvar_per<<endl;

    cellnum=0;
    LOOP
    if(a->flag(i,j,k)>0 && a->solid(i,j,k)>0)
    ++cellnum;

    paverage=int(cellnum/p->M10);

    ddout<<"average per process: "<<paverage<<endl;

    // maindir
    if(a->mz>=a->mx && a->mz>=a->my)
    maindir=3;

    if(a->my>=a->mx && a->my>=a->mz)
    maindir=2;

    if(a->mx>=a->my && a->mx>=a->mz)
    maindir=1;

    cout<<"xvar: "<<xvar<<" yvar: "<<yvar<<" zvar: "<<zvar<<endl;

    // turn off partition directions
    if((xvar_per>=p->M41 && yvar_per<p->M41 && zvar_per<p->M41)
    || (xvar_per<p->M41 && yvar_per>=p->M41 && zvar_per<p->M41)
    || (xvar_per<p->M41 && yvar_per<p->M41 && zvar_per>=p->M41))
    {
        if(xvar>=p->M45*yvar && xvar>=p->M45*zvar)
        if((p->M12==1 && p->knoy>p->M10) || (p->M13==1 && p->knoz>p->M10))
        p->M11=0;

        if(yvar>p->M45*xvar && yvar>=p->M45*zvar)
        if((p->M11==1 && p->knox>p->M10) || (p->M13==1 && p->knoz>p->M10))
        p->M12=0;

        if(zvar>p->M45*xvar && zvar>p->M45*yvar)
        if((p->M11==1 && p->knox>p->M10) || (p->M12==1 && p->knoy>p->M10))
        p->M13=0;
    }

    //---

    if(xvar_per>=p->M41 && yvar_per>=p->M41 && zvar_per<p->M41)
    {
        if(xvar>yvar)
        if((p->M12==1 && p->knoy>p->M10) || (p->M13==1 && p->knoz>p->M10))
        p->M11=0;


        if(xvar<yvar)
        if((p->M11==1 && p->knox>p->M10) || (p->M13==1 && p->knoz>p->M10))
        p->M12=0;

        if(xvar==yvar)
        {
            if(maindir==1)
            if((p->M11==1 && p->knox>p->M10) || (p->M13==1 && p->knoz>p->M10))
            p->M12=0;

            if(maindir==2)
            if((p->M12==1 && p->knoy>p->M10) || (p->M13==1 && p->knoz>p->M10))
            p->M11=0;
        }
    }

    if(xvar_per>=p->M41 && yvar_per<p->M41 && zvar_per>=p->M41)
    {
        if(xvar>zvar)
        if((p->M12==1 && p->knoy>p->M10) || (p->M13==1 && p->knoz>p->M10))
        p->M11=0;

        if(xvar<zvar)
        if((p->M11==1 && p->knox>p->M10) || (p->M12==1 && p->knoy>p->M10))
        p->M13=0;

        if(xvar==zvar)
        {
            if(maindir==1)
            if((p->M11==1 && p->knox>p->M10) || (p->M12==1 && p->knoy>p->M10))
            p->M13=0;

            if(maindir==3)
            if((p->M12==1 && p->knoy>p->M10) || (p->M13==1 && p->knoz>p->M10))
            p->M11=0;
        }
    }

    if(xvar_per<p->M41 && yvar_per>=p->M41 && zvar_per>=p->M41)
    {
        if(yvar>zvar)
        if((p->M11==1 && p->knox>p->M10) || (p->M13==1 && p->knoz>p->M10))
        p->M12=0;

        if(yvar<zvar)
        if((p->M11==1 && p->knox>p->M10) || (p->M12==1 && p->knoy>p->M10))
        p->M13=0;

        if(yvar==zvar)
        {
            if(maindir==2)
            if((p->M11==1 && p->knox>p->M10) || (p->M12==1 && p->knoy>p->M10))
            p->M13=0;

            if(maindir==3)
            if((p->M11==1 && p->knox>p->M10) || (p->M13==1 && p->knoz>p->M10))
            p->M12=0;
        }

    }

    int check=0;
    if(xvar_per>=p->M41 && yvar_per>=p->M41 && zvar_per>=p->M41)
    {
        if(xcross_m<=ycross_m && xcross_m<=zcross_m)
        if(p->M11==1 && p->knox>p->M10)
        {
            p->M12=0;
            p->M13=0;
            check=1;
        }

        if(ycross_m<=xcross_m && ycross_m<=zcross_m)
        if(p->M12==1 && p->knoy>p->M10)
        {
            p->M11=0;
            p->M13=0;
            check=1;
        }

        if(zcross_m<=xcross_m && zcross_m<=ycross_m)
        if(p->M13==1 && p->knoz>p->M10)
        {
            p->M11=0;
            p->M12=0;
            check=1;
        }

        if(check==0)
        {
            if(a->knox>=a->knoy && a->knox>=a->knoz)
            if(p->M11==1 && p->knox>p->M10)
            {
                p->M12=0;
                p->M13=0;
            }

            if(a->knoy>a->knox && a->knoy>=a->knoz)
            if(p->M12==1 && p->knoy>p->M10)
            {
                p->M11=0;
                p->M13=0;
            }

            if(a->knoz>a->knox && a->knoz>a->knoy)
            if(p->M13==1 && p->knoz>p->M10)
            {
                p->M11=0;
                p->M12=0;
            }
        }
    }

    // maindir change
    if(maindir==1 && p->M11==0)
    {
        if(a->my>=a->mz && p->M12==1)
        {
            cout<<"MX: "<<a->mx<<" MY: "<<a->my<<endl;
            maindir=2;
        }

        if(a->mz>=a->my && p->M13==1)
        maindir=3;

        if(p->M12==0)
        maindir=3;

        if(p->M13==0)
        maindir=2;
    }

    if(maindir==2 && p->M12==0)
    {
        if(a->mx>=a->mz && p->M11==1)
        maindir=1;

        if(a->mz>=a->mx && p->M13==1)
        maindir=3;

        if(p->M11==0)
        maindir=3;

        if(p->M13==0)
        maindir=1;
    }

    if(maindir==3 && p->M13==0)
    {
        if(a->mx>=a->my && p->M11==1)
        maindir=1;

        if(a->my>=a->mx && p->M12==1)
        maindir=2;

        if(p->M11==0)
        maindir=2;

        if(p->M12==0)
        maindir=1;
    }

    ddout<<"maindir_final: "<<maindir<<endl;
    ddout<<"M11: "<<p->M11<<endl;
    ddout<<"M12: "<<p->M12<<endl;
    ddout<<"M13: "<<p->M13<<endl;
}
}
