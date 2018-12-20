/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"fnpf_discrete_weights.h"
#include"lexer.h"
#include"field.h"
#include"slice.h"
#include"vec.h"

fnpf_discrete_weights::fnpf_discrete_weights(lexer* p)
{
}

fnpf_discrete_weights::~fnpf_discrete_weights()
{
}

void fnpf_discrete_weights::ck_weights(lexer *p, double *ck, int kval, int nval, int numpt)
{
   /*
  subroutine apply_fd(nin, maxorder, xdata, ydata, xtgt, out)
    integer, intent(in)    :: nin, maxorder
    real(dp), intent(in)    :: xdata(0:), ydata(0:), xtgt
    real(dp), intent(out) :: out(0:)

    integer :: j
    real(dp) :: c(0:nin-1, 0:maxorder)

    call populate_weights(xtgt, xdata, nin-1, maxorder, c)
    forall(j=0:maxorder) out(j) = sum(c(:, j)*ydata)
  end subroutine


  subroutine populate_weights(z, x, nd, m, c)
    !
    !  Input Parameters
    !    z            -  location where approximations are to be
    !                    accurate
    !    x(0:nd)      -  grid point locations, found in x(0:n)
    !    nd           -  dimension of x- and c-arrays in calling
    !                    program x(0:nd) and c(0:nd, 0:m), respectively
    !    m            -  highest derivative for which weights are
    !                    sought
    !
    !  Output Parameter
    !    c(0:nd,0:m)  -  weights at grid locations x(0:n) for
    !                    derivatives of order 0:m, found in c(0:nd, 0:m)
    !
    real(dp), intent(in)    :: z
    integer,  intent(in)    :: nd, m
    real(dp), intent(in)    :: x(0:nd)
    real(dp), intent(out)   :: c(0:nd, 0:m)

    real(dp) :: c1, c2, c3, c4, c5
    integer  :: i, j, k, mn*/

    // ini
    
    
    
    /*c1 = 1.0;
    c4 = x[0] - z;
    c = 0.0;
    c[0][0] = 1.0;
    

    for(i=0;i<nd; ++i)
    {
    
      mn = MIN(i,m)
      c2 = 1.0;
      c5 = c4;
      c4 = x[i] - z;
      

      for(j=0; j<i-1; ++j)
      {
        c3 = x[i] - x[j];
        c2 = c2*c3;
        
        if (j==i-1) 
        for(k=mn;k>1;k--)
        {
        c[i][k] = c1*(k*c[i-1][k-1] - c5*c[i-1][k]))/c2
 
        c(i, 0) = -c1*c5*c[i-1][0]/c2;
        }
        
        for(k=mn;k>1;k--)
          c(j, k) = (c4*c(j, k) - k*c(j, k-1))/c3
   
        c(j, 0) = c4*c(j, 0)/c3
      }
      c1 = c2;
      
    }

end module
}




n=length(x); c=zeros(m+1,n); c1=1; c4=x(1)-z; c(1,1)=1;
for i=2:n
   mn=min(i,m+1); c2=1; c5=c4; c4=x(i)-z;
   for j=1:i-1
      c3=x(i)-x(j);  c2=c2*c3;
      if j==i-1 
         c(2:mn,i)=c1*((1:mn-1)'.*c(1:mn-1,i-1)-c5*c(2:mn,i-1))/c2;
         c(1,i)=-c1*c5*c(1,i-1)/c2;
      end
      c(2:mn,j)=(c4*c(2:mn,j)-(1:mn-1)'.*c(1:mn-1,j))/c3;
      c(1,j)=c4*c(1,j)/c3;
   end
   c1=c2;
end


*/

}


