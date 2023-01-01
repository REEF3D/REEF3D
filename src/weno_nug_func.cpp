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

#include"weno_nug_func.h"
#include"lexer.h"
#include"fdm.h"
#include"flux_face_CDS2.h"
#include"flux_face_CDS2_vrans.h"
#include"flux_face_FOU.h"
#include"flux_face_FOU_vrans.h"
#include"flux_face_QOU.h"

weno_nug_func::weno_nug_func(lexer* p):epsilon(0.0),psi(1.0e-6)
{
    ini(p);
    
    pp=p;
}

weno_nug_func::~weno_nug_func()
{
}

void weno_nug_func::ini(lexer* p)
{
    if(iniflag==0)
    {
    p->Darray(qfx,p->knox+8,2,6,2);
    p->Darray(qfy,p->knoy+8,2,6,2);
    p->Darray(qfz,p->knoz+8,2,6,2);
    
    p->Darray(cfx,p->knox+8,2,6);
    p->Darray(cfy,p->knoy+8,2,6);
    p->Darray(cfz,p->knoz+8,2,6);
    
    p->Darray(isfx,p->knox+8,2,6,3);
    p->Darray(isfy,p->knoy+8,2,6,3);
    p->Darray(isfz,p->knoz+8,2,6,3);
    
    precalc_qf(p);
    precalc_cf(p);
    precalc_isf(p);
    
    /*
    if(p->mpirank==0)
    {
        IBLOOP
        {
            
        cout<<"QFX_0: "<<qfx[IP][0][0][0]<<" "<<qfx[IP][0][0][1]<<" | "<<qfx[IP][0][1][0]<<" "<<qfx[IP][0][1][1]<<" | "<<qfx[IP][0][2][0]<<" "<<qfx[IP][0][2][1]<<" | "<<qfx[IP][0][3][0]<<" "<<qfx[IP][0][3][1]<<" | "
            <<qfx[IP][0][4][0]<<" "<<qfx[IP][0][4][1]<<" | "<<qfx[IP][0][5][0]<<" "<<qfx[IP][0][5][1]<<endl;

        cout<<"QFX_1: "<<qfx[IP][1][0][0]<<" "<<qfx[IP][1][0][1]<<" | "<<qfx[IP][1][1][0]<<" "<<qfx[IP][0][1][1]<<" | "<<qfx[IP][1][2][0]<<" "<<qfx[IP][1][2][1]<<" | "<<qfx[IP][1][3][0]<<" "<<qfx[IP][1][3][1]<<" | "
            <<qfx[IP][1][4][0]<<" "<<qfx[IP][1][4][1]<<" | "<<qfx[IP][1][5][0]<<" "<<qfx[IP][1][5][1]<<endl;
        }
        
        cout<<endl<<endl<<endl;
        
        
        JBLOOP
        {

        cout<<"QFY_0: "<<qfy[JP][0][0][0]<<" "<<qfy[JP][0][0][1]<<" | "<<qfy[JP][0][1][0]<<" "<<qfy[JP][0][1][1]<<" | "<<qfy[JP][0][2][0]<<" "<<qfy[JP][0][2][1]<<" | "<<qfy[JP][0][3][0]<<" "<<qfy[JP][0][3][1]<<" | "
            <<qfy[JP][0][4][0]<<" "<<qfy[JP][0][4][1]<<" | "<<qfy[JP][0][5][0]<<" "<<qfy[JP][0][5][1]<<endl;

        cout<<"QFY_1: "<<qfy[JP][1][0][0]<<" "<<qfy[JP][1][0][1]<<" | "<<qfy[JP][1][1][0]<<" "<<qfy[JP][0][1][1]<<" | "<<qfy[JP][1][2][0]<<" "<<qfy[JP][1][2][1]<<" | "<<qfy[JP][1][3][0]<<" "<<qfy[JP][1][3][1]<<" | "
            <<qfy[JP][1][4][0]<<" "<<qfy[JP][1][4][1]<<" | "<<qfy[JP][1][5][0]<<" "<<qfy[JP][1][5][1]<<endl;
        }
                         
        
        
        
        cout<<endl<<endl<<endl;
        
        KBLOOP
        {

        cout<<"QFZ_0: "<<qfz[KP][0][0][0]<<" "<<qfz[KP][0][0][1]<<" | "<<qfz[KP][0][1][0]<<" "<<qfz[KP][0][1][1]<<" | "<<qfz[KP][0][2][0]<<" "<<qfz[KP][0][2][1]<<" | "<<qfz[KP][0][3][0]<<" "<<qfz[KP][0][3][1]<<" | "
            <<qfz[KP][0][4][0]<<" "<<qfz[KP][0][4][1]<<" | "<<qfz[KP][0][5][0]<<" "<<qfz[KP][0][5][1]<<endl;

        cout<<"QFZ_1: "<<qfz[KP][1][0][0]<<" "<<qfz[KP][1][0][1]<<" | "<<qfz[KP][1][1][0]<<" "<<qfz[KP][0][1][1]<<" | "<<qfz[KP][1][2][0]<<" "<<qfz[KP][1][2][1]<<" | "<<qfz[KP][1][3][0]<<" "<<qfz[KP][1][3][1]<<" | "
            <<qfz[KP][1][4][0]<<" "<<qfz[KP][1][4][1]<<" | "<<qfz[KP][1][5][0]<<" "<<qfz[KP][1][5][1]<<endl;
        }
                         
        
        
        
        cout<<endl<<endl<<endl;
        
        IBLOOP
        {

        cout<<"IFX_0: "<<isfx[IP][0][0][0]<<" "<<isfx[IP][0][0][1]<<" "<<isfx[IP][0][0][2]<<" | "<<isfx[IP][0][1][0]<<" "<<isfx[IP][0][1][1]<<" "<<isfx[IP][0][1][2]<<" | "<<isfx[IP][0][2][0]<<" "<<isfx[IP][0][2][1]<<" "<<isfx[IP][0][2][2]
        <<" | "<<isfx[IP][0][3][0]<<" "<<isfx[IP][0][3][1]<<" "<<isfx[IP][0][3][2]<<" | "<<isfx[IP][0][4][0]<<" "<<isfx[IP][0][4][1]<<" "<<isfx[IP][0][4][2]<<" | "<<isfx[IP][0][5][0]<<" "<<isfx[IP][0][5][1]<<" "<<isfx[IP][0][5][2]<<endl;
            
        cout<<"IFX_1: "<<isfx[IP][1][0][0]<<" "<<isfx[IP][1][0][1]<<isfx[IP][1][0][2]<<" | "<<isfx[IP][1][1][0]<<" "<<isfx[IP][0][1][1]<<" "<<isfx[IP][0][1][2]<<" | "<<isfx[IP][1][2][0]<<" "<<isfx[IP][1][2][1]<<" "<<isfx[IP][1][2][2]
        <<" | "<<isfx[IP][1][3][0]<<" "<<isfx[IP][1][3][1]<<" "<<isfx[IP][1][3][2]<<" | "<<isfx[IP][1][4][0]<<" "<<isfx[IP][1][4][1]<<" "<<isfx[IP][1][4][2]<<" | "<<isfx[IP][1][5][0]<<" "<<isfx[IP][1][5][1]<<" "<<isfx[IP][1][5][2]<<endl;
        }
        
        cout<<endl<<endl<<endl;
        
        JBLOOP
        {
        cout<<"IFY_0: "<<isfy[JP][0][0][0]<<" "<<isfy[JP][0][0][1]<<" "<<isfy[JP][0][0][2]<<" | "<<isfy[JP][0][1][0]<<" "<<isfy[JP][0][1][1]<<" "<<isfy[JP][0][1][2]<<" | "<<isfy[JP][0][2][0]<<" "<<isfy[JP][0][2][1]<<" "<<isfy[JP][0][2][2]
        <<" | "<<isfy[JP][0][3][0]<<" "<<isfy[JP][0][3][1]<<" "<<isfy[JP][0][3][2]<<" | "<<isfy[JP][0][4][0]<<" "<<isfy[JP][0][4][1]<<" "<<isfy[JP][0][4][2]<<" | "<<isfy[JP][0][5][0]<<" "<<isfy[JP][0][5][1]<<" "<<isfy[JP][0][5][2]<<endl;

        cout<<"IFY_1: "<<isfy[JP][1][0][0]<<" "<<isfy[JP][1][0][1]<<isfy[JP][1][0][2]<<" | "<<isfy[JP][1][1][0]<<" "<<isfy[JP][0][1][1]<<" "<<isfy[JP][0][1][2]<<" | "<<isfy[JP][1][2][0]<<" "<<isfy[JP][1][2][1]<<" "<<isfy[JP][1][2][2]
        <<" | "<<isfy[JP][1][3][0]<<" "<<isfy[JP][1][3][1]<<" "<<isfy[JP][1][3][2]<<" | "<<isfy[JP][1][4][0]<<" "<<isfy[JP][1][4][1]<<" "<<isfy[JP][1][4][2]<<" | "<<isfy[JP][1][5][0]<<" "<<isfy[JP][1][5][1]<<" "<<isfy[JP][1][5][2]<<endl;
        }
        
        cout<<endl<<endl<<endl;
        
        KBLOOP
        {

        cout<<"IFZ_0: "<<isfz[KP][0][0][0]<<" "<<isfz[KP][0][0][1]<<" "<<isfz[KP][0][0][2]<<" | "<<isfz[KP][0][1][0]<<" "<<isfz[KP][0][1][1]<<" "<<isfz[KP][0][1][2]<<" | "<<isfz[KP][0][2][0]<<" "<<isfz[KP][0][2][1]<<" "<<isfz[KP][0][2][2]
        <<" | "<<isfz[KP][0][3][0]<<" "<<isfz[KP][0][3][1]<<" "<<isfz[KP][0][3][2]<<" | "<<isfz[KP][0][4][0]<<" "<<isfz[KP][0][4][1]<<" "<<isfz[KP][0][4][2]<<" | "<<isfz[KP][0][5][0]<<" "<<isfz[KP][0][5][1]<<" "<<isfz[KP][0][5][2]<<endl;

        cout<<"IFZ_1: "<<isfz[KP][1][0][0]<<" "<<isfz[KP][1][0][1]<<" "<<isfz[KP][1][0][2]<<" | "<<isfz[KP][1][1][0]<<" "<<isfz[KP][0][1][1]<<" "<<isfz[KP][0][1][2]<<" | "<<isfz[KP][1][2][0]<<" "<<isfz[KP][1][2][1]<<" "<<isfz[KP][1][2][2]
        <<" | "<<isfz[KP][1][3][0]<<" "<<isfz[KP][1][3][1]<<" "<<isfz[KP][1][3][2]<<" | "<<isfz[KP][1][4][0]<<" "<<isfz[KP][1][4][1]<<" "<<isfz[KP][1][4][2]<<" | "<<isfz[KP][1][5][0]<<" "<<isfz[KP][1][5][1]<<" "<<isfz[KP][1][5][2]<<endl;
        }
        
        cout<<endl<<endl<<endl;
        
        
        IBLOOP
        {


        cout<<"CFX_0: "<<cfx[IP][0][0]<<" | "<<cfx[IP][0][1]<<" | "<<cfx[IP][0][2]<<" | "<<cfx[IP][0][3]<<" | "
            <<cfx[IP][0][4]<<" | "<<cfx[IP][0][5]<<endl;

        cout<<"CFX_1: "<<cfx[IP][1][0]<<" | "<<cfx[IP][1][1]<<" | "<<cfx[IP][1][2]<<" | "<<cfx[IP][1][3]<<" | "
            <<cfx[IP][1][4]<<" | "<<cfx[IP][1][5]<<endl;
        }
        
        cout<<endl<<endl<<endl;
        
        
        JBLOOP
        {
        cout<<"CFY_0: "<<cfy[JP][0][0]<<" | "<<cfy[JP][0][1]<<" | "<<cfy[JP][0][2]<<" | "<<cfy[JP][0][3]<<" | "
            <<cfy[JP][0][4]<<" | "<<cfy[JP][0][5]<<endl;

        cout<<"CFY_1: "<<cfy[JP][1][0]<<" | "<<cfy[JP][1][1]<<" | "<<cfy[JP][1][2]<<" | "<<cfy[JP][1][3]<<" | "
            <<cfy[JP][1][4]<<" | "<<cfy[JP][1][5]<<endl;
        }
        
        cout<<endl<<endl<<endl;
        
        
        KBLOOP
        {
        cout<<"CFZ_0: "<<cfz[KP][0][0]<<" | "<<cfz[KP][0][1]<<" | "<<cfz[KP][0][2]<<" | "<<cfz[KP][0][3]<<" | "
            <<cfz[KP][0][4]<<" | "<<cfz[KP][0][5]<<endl;
            
        cout<<"CFZ_1: "<<cfz[KP][1][0]<<" | "<<cfz[KP][1][1]<<" | "<<cfz[KP][1][2]<<" | "<<cfz[KP][1][3]<<" | "
            <<cfz[KP][1][4]<<" | "<<cfz[KP][1][5]<<endl;
        }
    }*/
                      
                      
                      
    iniflag=1;    
    }              
                      
}


double ****weno_nug_func::qfx,****weno_nug_func::qfy,****weno_nug_func::qfz;
double ***weno_nug_func::cfx,***weno_nug_func::cfy,***weno_nug_func::cfz;
double ****weno_nug_func::isfx,****weno_nug_func::isfy,****weno_nug_func::isfz;
int weno_nug_func::iniflag(0);
