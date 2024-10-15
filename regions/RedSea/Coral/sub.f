	  module sub

	  use param

      contains

      SUBROUTINE derivs (x,y,dydx)

        integer,parameter :: nr=2
        real :: dydx(nr),y(nr),x

        if (x .le. dt) then 
        
        dH=JHGm
	    dS=JSGm

	    else

		dS=(JSG-JST)!*S
		dH=(JHG-JHT)!*H

		end if

        dydx(1)=dH
        dydx(2)=dS


        return
      END SUBROUTINE derivs


      SUBROUTINE rk4(y,dydx,nr,x,h,yout,derivs)
          implicit none
          integer :: i,nr,nmax
          real :: h,x,dydx(nr),y(nr),yout(nr)
          external derivs
          parameter (nmax=50)
          real ::  h6,hh,xh,dym(nmax),dyt(nmax),yt(nmax)
          hh = h*0.5d0
          h6 = h/6d0
          xh = x+hh
        do i=1,nr
          yt(i) = y(i) + hh*dydx(i)
        end do
          call derivs(xh,yt,dyt)
        do i=1,nr
          yt(i) = y(i) + hh*dyt(i)
        end do
          call derivs(xh,yt,dym)
        do i=1,nr
          yt(i) = y(i) + h*dym(i)
          dym(i) = dyt(i) + dym(i)
        end do
          call derivs(x+h,yt,dyt)
        do i=1,nr
          yout(i) = y(i) + h6*(dydx(i)+dyt(i)+2d0*dym(i))
        end do
        return
      END SUBROUTINE rk4


      Subroutine synth(ysynth,x,y,m)

      implicit none
      real ysynth,m,x,y

	  ysynth=1 / ((1 / m) + (1 / x) + (1 / y) - (1 / (x + y)))

	  end Subroutine synth


	end module sub



