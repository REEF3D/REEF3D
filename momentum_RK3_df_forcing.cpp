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

#include"momentum_RK3_df.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"6DOF_df.h"
#include"density_f.h"

void momentum_RK3_df::forcing(lexer* p, fdm* a, ghostcell* pgc, sixdof_df* p6dof_df, field& uvel, field& vvel, field& wvel, field& uveln, field& vveln, field& wveln, double alpha, vrans* pvrans, vector<net*>& pnet)
{
    //- Calculate rigid body velocity u_fb
    Eigen::Matrix<double, 6, 1> u_fb;

    // Quaternion Eulerian approach
    bool conv = true;
    forces(p ,a, pgc, p6dof_df, uveln, vveln, wveln, conv, alpha);
    p6dof_df->corrector(p,a,pgc,alpha,pvrans,pnet);
    p6dof_df->updateFSI(p,a,pgc,conv);
    
    for (int it = 0; it < 1; it++)
    {

        if (p->knoy == 1)
        {
            u_fb << p->ufbi, 0.0, p->wfbi, 0.0, p->qfbi, 0.0;
        }
        else
        {
            u_fb << p->ufbi, p->vfbi, p->wfbi, p->pfbi, p->qfbi, p->rfbi;
        }


        //- Calculate rigid body forcing field

        ULOOP
        {
            uf(i,j,k) = u_fb(0) + u_fb(4)*(p->pos1_z() - p->zg) - u_fb(5)*(p->pos1_y() - p->yg);
        }
        VLOOP
        {
            vf(i,j,k) = u_fb(1) + u_fb(5)*(p->pos2_x() - p->xg) - u_fb(3)*(p->pos2_z() - p->zg);
        }
        WLOOP
        {
            wf(i,j,k) = u_fb(2) + u_fb(3)*(p->pos3_y() - p->yg) - u_fb(4)*(p->pos3_x() - p->xg);
        }

        pgc->start1(p,uf,gcval_u);  
        pgc->start2(p,vf,gcval_v);  
        pgc->start3(p,wf,gcval_w);

        //- Calculate forcing field
        ULOOP
        {
           fx(i,j,k) = a->fbh1(i,j,k)*(uf(i,j,k) - uvel(i,j,k))/(alpha*p->dt);
        }
        VLOOP
        {
           fy(i,j,k) = a->fbh2(i,j,k)*(vf(i,j,k) - vvel(i,j,k))/(alpha*p->dt);
        }   
        
        WLOOP
        {
           fz(i,j,k) = a->fbh3(i,j,k)*(wf(i,j,k) - wvel(i,j,k))/(alpha*p->dt);
        }
      
        pgc->start1(p,fx,gcval_u);
        pgc->start2(p,fy,gcval_v);
        pgc->start3(p,fz,gcval_w);           
   
        //- Correct rigid body forces
        double fz_old = Zfb;

        if (fabs(Zfb - fz_old)/fabs(fz_old) < 0.01) break;
    }

    p6dof_df->saveTimeStep(p,alpha);
    p6dof_df->interface(p,true);
}


void momentum_RK3_df::forces(lexer* p, fdm* a, ghostcell* pgc, sixdof_df* p6dof,field& uveln, field& vveln, field& wveln, bool converged, double alpha)
{
    double fn, dm;
    
    Xfb = 0.0;
    Yfb = 0.0;
    Zfb = 0.0;
    Kfb = 0.0;
    Mfb = 0.0;
    Nfb = 0.0;
    
    if (p->count > 1)
    {
        ULOOP
        {
            // Forcing term
            fn = -fx(i,j,k);

            // Mass
            dm = a->fbh1(i,j,k)*pdensity->roface(p,a,1,0,0)*p->DXP[IP]*p->DYN[JP]*p->DZN[KP];

            // Forces and moments on body
            Xfb += fn*dm + dm*(p->ufbi - uveln(i,j,k))/(p->dt);
           
            Mfb += fn*dm*(p->pos1_z() - p->zg);
            Nfb -= fn*dm*(p->pos1_y() - p->yg);
        }
            
        VLOOP
        {
            // Forcing term
            fn = -fy(i,j,k);

            // Mass
            dm = a->fbh2(i,j,k)*pdensity->roface(p,a,0,1,0)*p->DXN[IP]*p->DYP[JP]*p->DZN[KP];
            
            // Forces and moments on body
            Yfb += fn*dm + dm*(p->vfbi - vveln(i,j,k))/(p->dt);
            
            Kfb -= fn*dm*(p->pos2_z() - p->zg);
            Nfb += fn*dm*(p->pos2_x() - p->xg);
        }

        WLOOP
        {
            // Forcing term
            fn = -fz(i,j,k);
            a->test(i,j,k) = fn;
            
            // Mass
            dm = a->fbh3(i,j,k)*pdensity->roface(p,a,0,0,1)*p->DXN[IP]*p->DYN[JP]*p->DZP[KP];
                
            // Forces and moments on body
            Zfb += fn*dm + dm*(p->wfbi - wveln(i,j,k))/(p->dt);

            Kfb += fn*dm*(p->pos3_y() - p->yg); 
            Mfb -= fn*dm*(p->pos3_x() - p->xg); 
        }
    }

    Xfb = pgc->globalsum(Xfb);
    Yfb = pgc->globalsum(Yfb);
    Zfb = pgc->globalsum(Zfb);

    Kfb = pgc->globalsum(Kfb);
    Mfb = pgc->globalsum(Mfb);
    Nfb = pgc->globalsum(Nfb);

    if (p->mpirank == 0) 
    {
       if (alpha == 1.0)
       {
        ofstream eposout;
        eposout.open("./REEF3D_6DOF/REEF3D_6DOF_forces.dat",std::ios_base::app);
        eposout<<p->simtime<<" \t "<<Xfb<<" \t "<<Yfb<<" \t "<<Zfb<<" \t "<<Kfb
        <<" \t "<<Mfb<<" \t "<<Nfb<<" \t "<<cd<<" \t "<<cq<<" \t "<<cl<<endl;   
        }
    }

    p6dof->forces_stl(p,a,pgc,alpha,Xfb,Yfb,Zfb,Kfb,Mfb,Nfb);
}
