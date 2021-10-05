      program convert
c
      integer ndim,ne,np
      parameter(ndim=1000000)
      real*8 x(ndim),y(ndim),slam(ndim),sfea(ndim),slam0,sfea0,
     $     b(ndim),pi
      pi=atan2(0.0,-1.0)
      write(*,*) 'enter slam0,sfea0 '
      read(*,*) slam0,sfea0
      slam0=slam0*pi/180
      sfea0=sfea0*pi/180
c
      open(14,file='fort.14')
      read(14,*)
      read(14,*) ne,np
      do i=1,np
         read(14,*) j,x(i),y(i),b(i)
         call invcp(x(i),y(i),slam(i),sfea(i),slam0,sfea0)
      enddo
      close(14)
      open(14,file='fort.14.latlon')
      write(14,*)
      write(14,*) ne,np
      do i=1,np
         write(14,*) i,slam(i)*180/pi-360,sfea(i)*180/pi,b(i)
      enddo
      close(14)
      stop
      end
c
      SUBROUTINE INVCP(XXCP,YYCP,RLAMBDA,PHI,RLAMBDA0,PHI0)
      IMPLICIT NONE
      REAL*8 XXCP,YYCP,RLAMBDA,PHI,RLAMBDA0,PHI0,R
      R=6378206.4d0
      RLAMBDA=RLAMBDA0+XXCP/(R*COS(PHI0))
      PHI=YYCP/R
      RETURN
      END

