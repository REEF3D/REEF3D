/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"nhflow_u_profile.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

nhflow_u_profile::nhflow_u_profile(lexer *p, fdm_nhf *d) : probenum(p->P67)
{
    p->Iarray(iloc,probenum);
    p->Iarray(jloc,probenum);
    p->Iarray(flag,probenum);
    
    // Create Folder
    if(p->mpirank==0)
    mkdir("./REEF3D_NHFLOW_U-Profile",0777);
    
    pout = new ofstream[probenum];
    
    ini_location(p,d);
}

nhflow_u_profile::~nhflow_u_profile()
{
    for(n=0;n<probenum;++n)
    pout[n].close();
    
    delete [] pout;
}

void nhflow_u_profile::start(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // prepare
    if(probenum>0)
    {
        
        // open file
        for(n=0;n<probenum;++n)
        if(flag[n]==1)
        {
        sprintf(name,"./REEF3D_NHFLOW_U-Profile/REEF3D-NHFLOW-U-Profile-%i-%i.dat",n+1,p->count);
        
        pout[n].open(name);
        
        //cout<<pout[n].is_open()<<" "<<n+1<<endl;
        
        

        pout[n]<<"U-Profile ID:  "<<n<<endl<<endl;
        pout[n]<<"x_coord     y_coord"<<endl;
        pout[n]<<n+1<<"\t "<<p->P67_x[n]<<"\t "<<p->P67_y[n]<<"\t "<<endl;
        pout[n]<<"simtime:  "<<p->simtime<<endl<<endl;

        pout[n]<<endl<<endl;
        
        pout[n]<<"z \t U  "<<endl;
        }
    
    
    // print 
    double xp,yp;
    
    for(n=0;n<probenum;++n)
    if(flag[n]==1)
    for(k=0;k<p->knoz;++k)
    {
    
        if(flag[n]>0)
        {
        xp=p->P67_x[n];
        yp=p->P67_y[n];
        
        i = iloc[n];
        j = jloc[n];
    
        uval = p->ccipol4V(d->U, d->WL, d->bed, xp, yp, p->ZSP[IJK]);
        }
 
    pout[n]<<setprecision(9)<<p->ZSP[IJK]<<" \t "<<uval<<endl;
    } 
    
    for(n=0;n<probenum;++n)
    if(flag[n]==1)
    pout[n].close();
    
    }
}

void nhflow_u_profile::ini_location(lexer *p, fdm_nhf *d)
{
    int check;
    int ii,jj;

    for(n=0;n<probenum;++n)
    {
    check=0;
    
    iloc[n]=p->posc_i(p->P67_x[n]);
    
    if(p->j_dir==0)
    {
    jloc[n]=0;
    j=0;
    p->P67_y[n] = 0.5*p->YP[JP];
    }
    
    if(p->j_dir==1)
    jloc[n]=p->posc_j(p->P67_y[n]);
    
    

    check=0;
    
    ii=iloc[n];
    jj=jloc[n];

    if(ii>=0 && ii<p->knox)
    if(jj>=0 && jj<p->knoy)
    check=1;
    
    cout<<p->mpirank<<" PROBE check: "<<check<<" i: "<<iloc[n]<<" j: "<<jloc[n]<<" ZSN: "<<p->ZSN[10+marge]<<endl;
    if(check==1)
    flag[n]=1;
    }
}
