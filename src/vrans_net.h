/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2024 Tobias Martin

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"vrans.h"
#include"increment.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4.h"
#include <Eigen/Dense>
#include<vector>

using namespace std;

#ifndef VRANS_NET_H_
#define VRANS_NET_H_

class vrans_net : public vrans, public increment
{
public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    
    typedef vector<Eigen::Vector3d> EigenMat;
	typedef Eigen::Matrix<double, 3, 3> Matrix3d;
    typedef vector<vector<double> > MatrixVd;    
    
	vrans_net(lexer*, ghostcell*);
	virtual ~vrans_net();

	virtual void initialize(lexer*, fdm*, ghostcell*);	
	virtual void start(lexer*, fdm*, ghostcell*, net*&, int);
    virtual void sed_update(lexer*, fdm*, ghostcell*){};
    
	virtual void u_source(lexer*, fdm*);
	virtual void v_source(lexer*, fdm*);
	virtual void w_source(lexer*, fdm*);
    
    virtual void ke_source(lexer*, fdm*, field&){};
    virtual void kw_source(lexer*, fdm*, field&){};
    virtual void eps_source(lexer*, fdm*, field&, field&){};
    virtual void omega_source(lexer*, fdm*, field&, field&){};
    
    virtual void eddyv_func(lexer*, fdm*){};
    
    virtual void veltimesave(lexer*,fdm*,ghostcell*);
   
private:

    double kernel_peskin(const double&);
    
    void distributeNetForces_x(lexer*,fdm*,ghostcell*, net*&, int);
    void distributeNetForces_y(lexer*,fdm*,ghostcell*, net*&, int);
    void distributeNetForces_z(lexer*,fdm*,ghostcell*, net*&, int);    
    void distributeCollarForces(lexer*,fdm*,ghostcell*, net*&, int);
    
    double kernel_radius;

    double *xstart, *xend, *ystart, *yend, *zstart, *zend;
    
    field1 Fx_net, kernel_x;
    field2 Fy_net, kernel_y;
    field3 Fz_net, kernel_z;
    
	int count;
};

#endif
