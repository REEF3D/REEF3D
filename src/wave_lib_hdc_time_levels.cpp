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
Architect: Hans Bihs
--------------------------------------------------------------------*/

#include"wave_lib_hdc.h"
#include"lexer.h"
#include"ghostcell.h"
#include<iostream>
#include<cstdlib>

/*--------------------------------------------------------------------
The time levels of the HDC input around t = wavetime + I 241:
q1 is the last level with simtime <= t, q2 the next one. E1.. hold the
data of q1 and E2.. that of q2; L1 and L2 track which levels they hold,
so each level is read exactly once and in order, for single files
(P 45 1 of the source) and continuous files (P 45 2) alike, and also when
the target time step spans several source levels.
--------------------------------------------------------------------*/

void wave_lib_hdc::time_levels(lexer *p, ghostcell *pgc, bool cfd)
{
    const double t = p->wavetime+p->I241;
    
    if(startup==0)
    {
        q1 = diter;
        q2 = diter+1;
        L1 = L2 = -1;
        qs = diter;
        
        if(file_conti==2)
        {
        filename_continuous(p,pgc);
        result.open(name, ios::binary);
        }
        
        startup=1;
    }
    
    // find q1
    while(q1+1-diter<numiter && simtime[q1+1-diter]<=t)
    ++q1;
        
    // find q2
    while(q2-diter<numiter && simtime[q2-diter]<t)
    ++q2;
    
    if(q2>=numiter+diter)
    endseries=1;
        
    q1=MIN(q1,numiter+diter-1);
    q2=MIN(q2,numiter+diter-1);
    
    if(q1==q2 && q2<numiter+diter-1)
    ++q2;
    
    if(q1==q2)
    endseries=1;
    
    if(endseries==1)
    return;
    
    // data of q1: from E2.. if q2 has become q1, otherwise read
    if(q1!=L1)
    {
        if(q1==L2)
        {
            if(cfd)
            fill_conti_cfd(p,pgc);
            
            if(!cfd)
            fill_conti_fnpf(p,pgc);
        }
        
        if(q1!=L2)
        read_level(p,pgc,cfd,1,q1);
        
        L1=q1;
    }
    
    // data of q2
    if(q2!=L2)
    {
        read_level(p,pgc,cfd,2,q2);
        L2=q2;
    }

    deltaT = simtime[q2-diter]-simtime[q1-diter];
    deltaT = deltaT>0.0?deltaT:1.0e20;
    
    t1 = (simtime[q2-diter]-t)/deltaT;
    t2 = (t-simtime[q1-diter])/deltaT;
    
    if(p->mpirank==0)
    cout<<"HDC  q1: "<<q1<<" q2: "<<q2<<" t1: "<<t1<<" t2: "<<t2<<" deltaT: "<<deltaT<<" simtime[q1]: "<<simtime[q1-diter]<<" simtime[q2]: "<<simtime[q2-diter]<<endl;
}

void wave_lib_hdc::read_level(lexer *p, ghostcell *pgc, bool cfd, int slot, int q)
{
    // continuous file: skip the records before q (they are never needed again)
    if(file_conti==2)
    {
        if(q<qs)
        {
        cout<<endl<<"!!! HDC: time level "<<q<<" is behind the continuous input file (next record "<<qs<<") !!!"<<endl<<endl;
        exit(1);
        }
        
        while(qs<q)
        {
            if(cfd)
            read_result_cfd(p,pgc,E2,U2,V2,W2,qs);
            
            if(!cfd)
            read_result_fnpf(p,pgc,E2,F2,qs);
            
            ++qs;
        }
        
        ++qs;
    }
    
    if(cfd && slot==1)
    read_result_cfd(p,pgc,E1,U1,V1,W1,q);
    
    if(cfd && slot==2)
    read_result_cfd(p,pgc,E2,U2,V2,W2,q);
    
    if(!cfd && slot==1)
    read_result_fnpf(p,pgc,E1,F1,q);
    
    if(!cfd && slot==2)
    read_result_fnpf(p,pgc,E2,F2,q);
}
