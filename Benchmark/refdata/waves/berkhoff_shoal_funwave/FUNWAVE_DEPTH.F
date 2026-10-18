!This program is used to generate the bathymetry for the Berkhoff(1982) shoal experiment

           Program shoal
           real(8) x,y,x1,y1,x0,y0
           real(8) dx,dy,theta
           real(8) hs,r2
           real(8) h(1000,1000)
           integer mx, ny, i,j          

           dx=0.1
           dy=0.2
                 
           mx=300
           ny=101
           x0=0.0
           y0=0.0          
           theta=4*atan(1.0)/180.0*20.0                    !theta =20 degree 

           h0 = .45
           h02 = 0.07

           ytmp  = 5.84-(h0-h02)/.02
           do i=1,mx
               do j=1,ny

                  x=-10.0+(j-1)*dy
                  y=+13.0-(i-1)*dx

                  x1=(x-x0)*cos(theta)-(y-y0)*sin(theta)               
                  y1=(y-y0)*cos(theta)+(x-x0)*sin(theta)

                  xp=(x-x0)*cos(theta)-(y-y0)*sin(theta)               
                  yp=(y-y0)*cos(theta)+(x-x0)*sin(theta)

c-------------longwen------------------------
                  if(y1.lt.5.84)then
                     hs=0.45-0.02*(5.84-y1)
                  else
                     hs=0.45
                  endif

                  r2=x1**2/16.0+y1**2/9.0                !r square

                  if (r2.lt.1.0)then
                     h(i,j)=hs+0.3
     &                     -0.5*sqrt(1-x1**2/25.0-y1**2/3.75/3.75)
                  else
                     h(i,j)=hs
                  endif
c-----------------

c------------Qin Chen---------------------
               ctmp1 = (xp/4.)**2+(yp/3.)**2
               ctmp2 = (xp/5.)**2+(yp/3.75)**2

               if (yp.ge.5.84) then
                  h(i,j) = h0
               elseif (yp.le.ytmp) then 
                  h(i,j) = h02
               elseif (ctmp1.gt.1.) then
                  h(i,j) = h0-0.02*(5.84-yp)
               else
                  h(i,j) = .75-.02*(5.84-yp)-0.5*sqrt(1.-ctmp2)
               endif
c-------------

           enddo
           enddo
          
           open(1,file='depth.dat')  !This is input bathymetry for funwave2d
               do j=1,ny
                  write(1,100)(h(i,j),i=1,mx)
               enddo
100            format(10000(f10.5,1x))     
           close(1)

           open(2,file='shoal.dat')   !This is for plotting the bathymetry by matlab
               do j=1,ny
                  write(2,200)(-h(i,j),i=1,mx)
               enddo
200            format(10000(f10.5,1x))
           close(2)
   
           end

