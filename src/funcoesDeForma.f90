!=================================================================================
!
!         programa de elementos finitos em fortran 90, 
!         baseado em: The Finite Element Method, Hughes, T. J. R., (2003)
!         + implementacoes de Abimael Loula
!
!         novo projeto: 
!         Eduardo Garcia,        bidu@lncc.br [4]
!         Tuane Lopes,           tuane@lncc.br [5]
!
!         LNCC/MCT
!         Petropolis, 02.2014
!=================================================================================
module mFuncoesDeForma

   use mGlobaisEscalares,only: zero,one,two,three,four,five,six,pt1667,pt25,pt5

   implicit none

   public :: oneshl, oneshg, shap2m, shlqrt, shgqrt, shlt, shgq
   public :: shg3d, shlq, shg2q, shl2q, shlq3d

   contains

   !=======================================================================       
   subroutine shgSurface_points(xl, det, shl, shg, nen, npoints, nsd)
      ! program to calculate global derivatives of shape functions and
      ! jacobian determinants at the integration points

      !     xl(j,i)    = global coordinates of nodes
      !     det(l)     = jacobian determinant
      !     shl(1,i,l) = local ("xi") derivative of shape function
      !     shl(2,i,l) = local ("eta") derivative of shape function
      !     shl(3,i,l) = local shape function
      !     shg(1,i,l) = x-derivative of shape function
      !     shg(2,i,l) = y-derivative of shape function
      !     shg(3,i,l) = shl(3,i,l)
      !              l = integration-point number
      !           nint = number of integration points

      ! program to calculate global derivatives of shape functions
      ! and jacobian determinants for a codimensional one element at the set of points:
      ! surface element if nsd==3 or line element of nsd == 2
      !
      !        xl(i,j) = global coordinates of nodes
      !     shl(1,j,l) = local ("xi") derivative of shape function j at a point l (nsd==2)
      !     shl(2,j,l) = local ("eta") derivative of shape function j at a point l (nsd==3)
      !   shl(nsd,j,l) = value of shape function j at a point l
      !     shg(1,j,l) = global x-derivative of shape function j at a point l
      !     shg(2,j,l) = global y-derivative of shape function j at a point l
      !     shg(3,j,l) = global z-derivative of shape function j at a point l
      ! shg(nsd+1,j,l) = shl(nsd,i,l)
      !         det(l) = square root of the determinant of the g matrix at a point l
      !        npoints = number of points
      !              l = point number

      integer :: l
      integer, intent(in) :: npoints, nen, nsd
      real*8, intent(in)  :: xl(nsd, nen)
      real*8, intent(in)  :: shl(nsd,  nen, npoints)
      real*8, intent(out) :: shg(nsd+1,nen, npoints)
      real*8, intent(out) :: det(npoints)

      do l = 1, npoints
         call shgSurface_point(xl, det(l), shl(:,:,l), shg(:,:,l), nen, nsd)
      end do
   end subroutine shgSurface_points

   !=======================================================================
   subroutine shgSurface_point(xl, det, shl, shg, nen, nsd)
      ! program to calculate global derivatives of shape functions 
      ! and jacobian determinants for a codimensional one element at a given point:
      ! surface element if nsd==3 or line element of nsd == 2
      !
      !      xl(i,j) = global coordinates of nodes
      !     shl(1,j) = local ("xi") derivative of shape function j at a given point (nsd==2)
      !     shl(2,j) = local ("eta") derivative of shape function j at a given point (nsd==3)
      !   shl(nsd,j) = value of shape function j at a given point
      !     shg(1,j) = global x-derivative of shape function j at a given point
      !     shg(2,j) = global y-derivative of shape function j at a given point
      !     shg(3,j) = global z-derivative of shape function j at a given point
      ! shg(nsd+1,j) = shl(nsd,i)
      !          det = square root of the determinant of the g matrix
      !
      !       d(i,j) = local j-derivative of the i-coordinate of a given point
      !            g = first fundamental form (d^T d)
      !         detG = determinant of matrix g
      !        g_inv = inverse of matrix g

      integer, intent(in) :: nen, nsd
      real*8, intent(in)  :: xl(nsd,nen), shl(nsd,nen)
      real*8, intent(out) :: shg(nsd+1,nen), det

      real*8 :: d(nsd,nsd-1), g(nsd-1,nsd-1), g_inv(nsd-1,nsd-1), detG
      integer :: i,j,no

      d=0.d0
      g=0.d0
      g_inv=0.d0
      detG=0.d0

      if(nen == 1) return

      do i=1,nsd
         do j=1,nsd-1
            d(i,j) = dot_product(xl(i,:),shl(j,:))
         end do
      end do

      do i=1,nsd-1
         do j=1,nsd-1
            g(i,j) = dot_product(d(:,i),d(:,j))
         end do
      end do

      if (nsd==2) then
         detG = g(1,1)
         g_inv = 1d0/detG
      else if (nsd==3) then
         detG = g(1,1)*g(2,2) - g(1,2)*g(2,1)
         g_inv(1,1) = g(2,2)
         g_inv(2,2) = g(1,1)
         g_inv(1,2) = -g(1,2)
         g_inv(2,1) = -g(2,1)
         g_inv(:,:) = g_inv(:,:)/detG
      else
         print*, "Erro na subrotina shgSurface_point: nsd deve ser 2 ou 3"
      end if

      if (detG <= 0d0) then
         print*, "Erro na subrotina shgSurface_point: determinante negativo"
         stop
      endif
      det = sqrt(detG)

      do no=1,nen
         do i=1,nsd
            shg(i,no) = 0d0
            do j=1,nsd-1
               shg(i,no) = shg(i,no) + d(i,j)*dot_product(g_inv(j,:),shl(1:nsd-1,no))
            end do
         end do
         shg(nsd+1,no) = shl(nsd,no)
      end do
   end subroutine shgSurface_point

   !=======================================================================
   subroutine oneshl(shl,w,npint,nen)
      ! .... program to calculate integration-rule weights, shape functions 
      ! and local derivatives for a two, three or four  node, 
      ! one-dimensional element

      !                 r = local element coordinate ("xi")
      !        shl(1,i,l) = local ("xi") derivative of shape function
      !        shl(2,i,l) = shape function
      !              w(l) = integration-rule weight
      !                 i = local node number 
      !                 l = integration-point number
      !              npint = number of integration points, eq. 1, 2, 3, 4, 6, 7 ,8

      implicit none

      !.... remove above card for single precision operation

      integer*4 ::  npint, nen, i, j, l, k
      real*8 :: shl(3,nen,npint),w(npint),ra(10),xa(10)
      real*8 :: r, aa, bb, aax, daj
      real*8, parameter ::  five9=0.5555555555555555d0
      real*8, parameter :: eight9=0.8888888888888888d0

      ! write(*,'(a)') '   em oneshl'
      ! write(*,*) ' em oneshl 1, shl  ', shl
      ! write(*,*) ' em oneshl 1,  w  ', w
      ! write(*,*) ' em oneshl 1, npint, nen  ', npint,nen

      if (npint.eq.1) then
         w(1)  = two
         ra(1) = zero
      endif

      if (nen.eq.1) xa(1) = zero

      if (npint.eq.2) then
         w(1) = one
         w(2) = one
         ra(1)=-.577350269189626
         ra(2)= .577350269189625
      endif

      if (nen.eq.2) then
         xa(1) = -one
         xa(2) =  one
      endif

      if (npint.eq.3) then
         w(1) = five9
         w(2) = eight9
         w(3) = five9
         ra(1)=-.774596669241483
         ra(2)= zero
         ra(3)= .774596669241483
      endif

      if(nen.eq.3) then
         xa(1)= -one
         xa(2)= one
         xa(3)= zero
      endif

      if (npint.eq.4) then
         w(1) = .347854845137454
         w(2) = .652145154862546
         w(3) = .652145154862546
         w(4) = .347854845137454
         ra(1)=-.861136311594053
         ra(2)=-.339981043584856
         ra(3)= .339981043584856
         ra(4)= .861136311594053
      endif

      if (nen.eq.4) then
         xa(1) = -one
         xa(2) = one
         xa(3) = -.333333333333333
         xa(4) =  .333333333333333
      endif

      if(npint.eq.5) then
         w(1) = .236926885056189
         w(2) = .478628670499366
         w(3) = .568888888888888
         w(4) = .478628670499366
         w(5) = .236926885056189
         ra(1)=-.906179845938664
         ra(2)=-.538469310105683
         ra(3)= zero
         ra(4)= .538469310105683
         ra(5)= .906179845938664
      endif
       
      if(nen.eq.5) then
         xa(1)= -one 
         xa(2)=  one 
         xa(3)= -pt5
         xa(4)= zero
         xa(5)= pt5
      endif
      
      if(npint.eq.6) then
         w(1) = .171324492397170
         w(2) = .360761573048139
         w(3) = .467913934572691
         w(4) = .467913934572691
         w(5) = .360761573048139
         w(6) = .171324492397170

         ra(1)=-.932469514203152
         ra(2)=-.661209386466265
         ra(3)=-.238619186083197
         ra(4)= .238619186083197
         ra(5)= .661209386466365
         ra(6)= .932469514203152
      endif

      if(nen.eq.6) then
         xa(1) = -one
         xa(2) =  one
         xa(3) = -.600000000000000
         xa(4) = -.200000000000000
         xa(5) =  .200000000000000
         xa(6) =  .600000000000000
      endif

      if(npint.eq.7) then
         w(1) = .129484966168870
         w(2) = .279705391489277 
         w(3) = .381830050505119
         w(4) = .417959183673469
         w(5) = .381830050505119
         w(6) = .279705391489277
         w(7) = .129484966168870

         ra(1)=-.949107912342759
         ra(2)=-.741531185599394
         ra(3)=-.405845151377397
         ra(4)= zero
         ra(5)= .405845151377397
         ra(6)= .741531185599394
         ra(7)= .949107912342759
      endif

      if(nen.eq.7) then
         xa(1) = -one
         xa(2) =  one
         xa(3) = -.666666666666666
         xa(4) = -.333333333333333
         xa(5) = zero
         xa(6) =  .333333333333333
         xa(7) =  .666666666666666  
      endif

      if(npint.eq.8) then
         w(1) = .101228536290376
         w(2) = .222381034453374
         w(3) = .313706645877887
         w(4) = .362683783378362
         w(5) = .362683783378362
         w(6) = .313706645877887
         w(7) = .222381034453374
         w(8) = .101228536290376

         ra(1)=-.960289856497536
         ra(2)=-.796666477413627
         ra(3)=-.525532409916329
         ra(4)=-.183434642495650
         ra(5)= .183434642495650
         ra(6)= .525532409916329
         ra(7)= .796666477413627
         ra(8)= .960289856497536
      endif
      if(nen.eq.8) then
         xa(1) = -one
         xa(2) =  one
         xa(3) = -0.71428571428571
         xa(4) = -0.42857142857143
         xa(5) = -0.14285714285714
         xa(6) = 0.14285714285714
         xa(7) = 0.42857142857143
         xa(8) = 0.71428571428571
      endif
      
      do 100 l = 1, npint
         r = ra(l)
         if(nen.eq.1) then
            shl(1,1,l) = zero
            shl(2,1,l) = one
            go to 100
         endif

         do 50 i = 1, nen
            aa = one
            bb = one
            aax = zero
            do 40 j =1, nen
               daj = one
               if (i .ne. j)then
                  aa = aa * (     r - xa(j) )
                  bb = bb * ( xa(i) - xa(j) )
                  do 30 k = 1, nen
                     if(k .ne. i .and. k .ne. j) daj = daj * ( r - xa(k))
                  30     continue
                  aax =aax + daj
               endif
            40    continue
            shl(3,i,l) = aa/bb
            shl(1,i,l) = aax/bb
            ! shl(2,i,l) = aax/bb
         50  continue
      100  continue
      
      return
   end subroutine
   
   !=======================================================================
   subroutine oneshg(xl,det,shl,shg,nen,npint,nesd,ns,nel,neg)
      ! .... program to calculate global derivatives of shape functions 
      ! and jacobian determinants for the bi-dimensional,
      ! elastic beam element

      !           xl(j,l) = global coordinates of integration points
      !        shl(1,i,l) = local ("xi") derivative of shape function
      !        shl(2,i,l) = shape function
      !        shg(1,i,l) = global ("arc-length") derivative of shape ftn
      !        shg(2,i,l) = shl(2,i,l)
      !            det(l) = euclidean length 
      !                 i = local node number 
      !                 j = global coordinate number
      !                 l = integration-point number
      !              npint = number of integration points

      use mLeituraEscrita, only: iecho

      implicit none

      ! .... remove above card for single precision operation
      integer*4 ::  npint, nen, nesd, ns, nel, neg
      real*8 :: xl(nesd,*),det(*),shl(3,nen,*),shg(3,nen,*)
      integer*4:: i, j, l
      real*8 :: x1, x2, dsqrt

      ! write(*,'(a)') '   em oneshg'
      do 400 l=1,npint
         det(l)=zero
         x1=0.d0
         x2=0.d0
         do 100 j=1,nen
            x1=x1+shl(1,j,l)*xl(1,j)
            x2=x2+shl(1,j,l)*xl(2,j)
         100   continue
         det(l)=dsqrt(x1*x1+x2*x2)

         if (det(l).le.zero) then
            write(iecho,1000) ns,nel,neg
            stop
         endif

         do 300 i=1,nen
            shg(1,i,l)=shl(1,i,l)/det(l)
            shg(3,i,l)=shl(3,i,l)
         300   continue
      400 continue
   
      1000 format(///,' oneshg - non-positive determinant in side ',i10,/, &
               ' in element ',i10,5x,' in element group  ',i10)
      
      return
   end subroutine

   !=======================================================================
   subroutine shap2m(s,xl,det,sh,nen,nesd)
      
      implicit none

      real*8  :: s
      integer*4:: nen, nesd
      double precision :: xl(nesd,*),sh(2,*)

      real*8  :: corr, corrh, det
      integer*4:: l
     
      ! shape function
      sh(2,1)=(1.d00-s)/2.d00
      sh(2,2)=(1.d00+s)/2.d00

      sh(1,1)=-.5d00
      sh(1,2)=.5d00

      ! 3 node correction
      if(nen.eq.3) then
         corr=1.d00-s*s
         corrh=corr/2
         sh(2,1)=sh(2,1)-corrh
         sh(2,2)=sh(2,2)-corrh
         sh(2,3)=corr

         corr=-2.d00*s
         corrh=-s
         sh(1,1)=sh(1,1)-corrh
         sh(1,2)=sh(1,2)-corrh
         sh(1,3)=corr
      end if

      det=0.d00
      do 100 l=1,nen
         det=det+xl(1,l)*sh(1,l)
      100   continue

      ! global derivatives
      do 200 l=1,nen
         sh(1,l)=sh(1,l)/det
      200   continue
      return
   end subroutine
   
   !=======================================================================      
   subroutine shlqrt(numLadosElem,npint,w,shlrt) 
      ! .... program to calculate integration-rule weights and shape functions
      ! for a raviart-thomas quadrilateral element

      !               s,t = local element coordinates ("xi", "eta", resp.)       
      !        shl(1,i,l) = local shape function ("xi")
      !        shl(2,i,l) = local shape function ("eta")
      !        shl(3,i,l) = divergence
      !              w(l) = integration-rule weight
      !                 i = local edge number
      !                 l = integration point number
      !              npint = number of integration points, eq. 1 or 4

      implicit none

      integer*4:: numLadosElem,npint
      real(8), dimension(*)        :: w
      real(8), dimension(3,numLadosElem,*) :: shlrt

      integer*4:: i,l
      real(8), dimension(25) :: ra,sa
      real(8), parameter :: zero=0.d0,pt5=0.5d0
      real(8), parameter :: r1=0.d0,w1=2.d0
      real(8), parameter :: r2=0.577350269189626d0,w2=1.d0      
      real(8), parameter :: r3a=0.774596669241483d00
      real(8), parameter :: w3a=0.555555555555556d00
      real(8), parameter :: r3b=0.d00
      real(8), parameter :: w3b=0.888888888888889d00
      real(8) :: r,s,w3ab,w3asq

      if (npint.eq.1) then
         w(1)=w1*w1
         ra(1)=r1
         sa(1)=r1
      end if

      if(npint.eq.4) then
         
         do i=1,4
            w(i)=w2
         end do
                
         ra(1)=r2
         sa(1)=r2
         ra(2)=-r2
         sa(2)=r2
         ra(3)=-r2
         sa(3)=-r2
         ra(4)=r2
         sa(4)=-r2     
      end if ! npint.eq.4
     
      if(npint.eq.9) then
      
         w3asq=w3a*w3a 
         do  i=1,4 
            w(i)=w3asq
         end do   
     
         ra(1)=r3a 
         sa(1)=r3a 
         ra(2)=-r3a 
         sa(2)=r3a 
         ra(3)=-r3a 
         sa(3)=-r3a 
         ra(4)=r3a 
         sa(4)=-r3a
     
         w3ab=w3a*w3b 
         do i=5,8 
            w(i)=w3ab 
         end do 
           
         ra(5)=r3a 
         sa(5)=r3b 
         ra(6)=r3b 
         sa(6)=r3a 
         ra(7)=-r3a 
         sa(7)=r3b 
         ra(8)=r3b 
         sa(8)=-r3a 

         w(9)=w3b*w3b 
         ra(9)=r3b 
         sa(9)=r3b
      end if ! npint.eq.9

      !  Numeracao local dos lados
      !   
      !            3
      !         _______
      !        |       |
      !      4 |       | 2
      !        |_______|
      !   
      !            1

      do l=1,npint
         if(numLadosElem.eq.4) then
            r=ra(l)
            s=sa(l)
           
            shlrt(1,1,l)=zero
            shlrt(2,1,l)=(s-1.d0)*pt5
            shlrt(3,1,l)=pt5
           
            shlrt(1,2,l)=(r+1.d0)*pt5
            shlrt(2,2,l)=zero
            shlrt(3,2,l)=pt5
           
            shlrt(1,3,l)=zero
            shlrt(2,3,l)=(s+1.d0)*pt5
            shlrt(3,3,l)=pt5

            shlrt(1,4,l)=(r-1.d0)*pt5
            shlrt(2,4,l)=zero
            shlrt(3,4,l)=pt5
         else
            write(*,*) 'numLadosElem.ne.4 em shlqrt'
            stop
         end if ! numLadosElem.eq.4
      end do

      return
   end subroutine

   !======================================================================
   subroutine shgqrt(numLadosElem,npint,hx,hy,shlrt,shgrt)
      !.... program to calculate shape functions
      ! for a raviart-thomas quadrilateral element

      !              s,t = local element coordinates ("xi", "eta", resp.)       
      !        shg(1,i,l) = global shape function ("xi")
      !        shg(2,i,l) = global shape function ("eta")
      !        shg(3,i,l) = divergence
      !              w(l) = integration-rule weight
      !                 i = local edge number
      !                 l = integration point number
      !              npint = number of integration points, eq. 1 or 4

      implicit none

      integer*4:: numLadosElem,npint      
      real(8) :: hx,hy
      real(8), dimension(3,numLadosElem,*) :: shlrt,shgrt

      integer*4:: i,j,l
      
      do l=1,npint
         do j=1,numLadosElem
            do i=1,2
               shgrt(i,j,l)=shlrt(i,j,l)
            end do     
            shgrt(3,j,l)=shlrt(3,j,l)/hx+shlrt(3,j,l)/hy
         end do
      end do

      return
   end subroutine
   
   !======================================================================
   subroutine shlt(shl,w,npint,nen)
      !.... program to calculate integration-rule weights, shape functions
      ! and local derivatives for a triangular element

      !        c1, c2, c3 = local element coordinates ("l1", "l2", "l3".)
      !        shl(j,i,l) = local ("j") derivative of shape function
      !        shl(3,i,l) = local  shape function
      !              w(l) = integration-rule weight
      !                 i = local node number
      !                 l = integration point number
      !              npint = number of integration points, eq. 1 or 4

      ! use mGlobaisEscalares
      implicit none

      !.... remove above card for single precision operation

      integer*4:: npint, nen
      real*8 :: shl(3,nen,*),w(*)

      real*8 :: c1, c2, c3, cl1(16),cl2(16),cl3(16)
      real*8 :: r1,w1,r2,w2,r3a,w3a,r3b1,w3b,r3b2,r7a,w7a,r7b,w7b,r7c,r7d,w7d,&
               r7e,r13a,w13a,r13b,w13b,r13c,w13c,r13d,w13d,r13e,r13f,r13g,r13h
      integer*4:: i, l

      data r1/0.33333333333333333333d00/,w1/1.d00/,&
           r2/0.5d00                   /,w2/0.3333333333333333333d00/,&
           r3a/0.3333333333333333333d00/,w3a/-0.5625d00/,&
           r3b1/0.6d00                 /,w3b/0.520833333333333333d00/,&
           r3b2/0.2d00                 /,&
           r7a/0.3333333333333333333d00/,w7a/0.225d00/,&
           r7b/0.0597158717d00         /,w7b/0.1323941527d00/,&
           r7c/0.4701420641d00         /,&
           r7d/0.7974269853d00         /,w7d/0.1259391805d00/,&
           r7e/0.1012865073d00         /,&
           r13a/0.3333333333333333333d00/,w13a/-0.1495700444676d0/,&
           r13b/0.479308067841923d0         /,w13b/0.175615257433204d0/,&
           r13c/0.260345966079038d0         /,w13c/0.053347235608839d0/&
           r13d/0.869739794195568d0         /,w13d/0.077113760890257d0/,&
           r13e/0.065130102902216d0        /,&
           r13f/0.638444188568809d0         /,&
           r13g/0.312865496004875d0         /&
           r13h/0.048690315425316d0         /
           
      if (npint.eq.1) then
         w(1)=w1/two
         cl1(1)=r1
         cl2(1)=r1
         cl3(1)=one-r1-r1
      end if

      if(npint.eq.3) then
         do 10 i=1,3
            w(i)=w2/two
         10       continue
         cl1(1)=r2
         cl2(1)=r2
         cl3(1)=zero
         cl1(2)=zero
         cl2(2)=r2
         cl3(2)=r2
         cl1(3)=r2
         cl2(3)=zero
         cl3(3)=r2
      end if

      if(npint.eq.4) then
         w(1)= w3a/two
         do 20 i=2,4
            w(i)=w3b/two
         20       continue
         cl1(1)=r3a
         cl2(1)=r3a
         cl3(1)=one - r3a - r3a
         cl1(2)=r3b1
         cl2(2)=r3b2
         cl3(2)=r3b2
         cl1(3)=r3b2
         cl2(3)=r3b1
         cl3(3)=r3b2
         cl1(4)=r3b2
         cl2(4)=r3b2
         cl3(4)=r3b1
      end if

      if(npint.eq.7) then
         w(1)= w7a/two
         do 30 i=2,4
            w(i)=w7b/two
         30       continue
         do 40 i=5,7
            w(i)=w7d/two
         40       continue
         cl1(1)=r7a
         cl2(1)=r7a
         cl3(1)=r7a
         do 50 i=2,4
            cl1(i)=r7c
            cl2(i)=r7c
            cl3(i)=r7c
         50 continue
         cl1(2)=r7b
         cl2(3)=r7b
         cl3(4)=r7b
         do 60 i=5,7
            cl1(i)=r7e
            cl2(i)=r7e
            cl3(i)=r7e
         60 continue
         cl1(5)=r7d
         cl2(6)=r7d
         cl3(7)=r7d
      end if

      if (npint.eq.13) then
         w(1)= w13a/two
         do 31 i=2,4
            w(i)=w13b/two
         31       continue
         do 41 i=5,7
            w(i)=w13c/two
         41       continue
         do 51 i=8,13
            w(i)=w13d/two
         51       continue
         cl1(1)=r13a
         cl2(1)=r13a
         cl3(1)=r13a
         do 52 i=2,4
            cl1(i)=r13c
            cl2(i)=r13c
            cl3(i)=r13c
         52 continue

         cl1(2)=r13b
         cl2(3)=r13b
         cl3(4)=r13b
         do 53 i=5,7
            cl1(i)=r13e
            cl2(i)=r13e
            cl3(i)=r13e
         53 continue

         cl1(5)=r13d
         cl2(6)=r13d
         cl3(7)=r13d

         cl1(8)=r13f
         cl2(8)=r13g
         cl3(8)=r13h
         cl1(9)=r13h
         cl2(9)=r13f
         cl3(9)=r13g
         cl1(10)=r13g
         cl2(10)=r13h
         cl3(10)=r13f

         cl1(11)=r13f
         cl2(11)=r13h
         cl3(11)=r13g
         cl1(12)=r13g
         cl2(12)=r13f
         cl3(12)=r13h
         cl1(13)=r13h
         cl2(13)=r13g
         cl3(13)=r13f
      endif

      do 200 l=1,npint

         c1 = cl1(l)
         c2 = cl2(l)
         c3 = cl3(l)
         shl(1,1,l)= one
         shl(2,1,l)= zero
         shl(3,1,l)= c1
         shl(1,2,l)= zero
         shl(2,2,l)= one
         shl(3,2,l)= c2
         shl(1,3,l)=-one
         shl(2,3,l)=-one
         shl(3,3,l)= c3
         if(nen.eq.6) then
            shl(1,4,l)= four * c2
            shl(2,4,l)= four * c1
            shl(3,4,l)= four * c1 * c2
            shl(1,5,l)=-four * c2
            shl(2,5,l)= four * (c3 - c2)
            shl(3,5,l)= four * c2 * c3
            shl(1,6,l)= four * (c3 - c1)
            shl(2,6,l)=-four * c1
            shl(3,6,l)= four * c3 * c1

            do 70 i=1,3
               shl(i,1,l)=shl(i,1,l)-pt5*(shl(i,4,l)+shl(i,6,l))
               shl(i,2,l)=shl(i,2,l)-pt5*(shl(i,4,l)+shl(i,5,l))
               shl(i,3,l)=shl(i,3,l)-pt5*(shl(i,5,l)+shl(i,6,l))
            70             continue
         end if

         if (nen.eq.10) then
            shl(1,1,l)= pt5*((three*c1-one)*(three*c1-two)+c1*three*(three*c1-two) &
                     +c1*three*(three*c1-one))
            shl(2,1,l)= zero
            shl(3,1,l)= pt5*c1*(three*c1-one)*(three*c1-two)

            shl(1,2,l)= zero
            shl(2,2,l)= pt5*((three*c2-one)*(three*c2-two)+c2*three*(three*c2-two) &
                     +c2*three*(three*c2-one))
            shl(3,2,l)= pt5*c2*(three*c2-one)*(three*c2-two)

            shl(1,3,l)= -pt5*((three*c3-one)*(three*c3-two)+c3*three*(three*c3-two) &
                     +c3*three*(three*c3-one))
            shl(2,3,l)=-pt5*((three*c3-one)*(three*c3-two)+c3*three*(three*c3-two) &
                     +c3*three*(three*c3-one))
            shl(3,3,l)=  pt5*c3*(three*c3-one)*(three*c3-two)

            shl(1,4,l) = (9.0d0/two)*(c2*(three*c1-one)+c1*c2*three)
            shl(2,4,l) = (9.0d0/two)*c1*(three*c1-one)
            shl(3,4,l) = (9.0d0/two)*c1*c2*(three*c1-one)

            shl(1,5,l) = (9.0d0/two)*c2*(three*c2-one)
            shl(2,5,l) = (9.0d0/two)*(c1*(three*c2-one)+c1*c2*three)
            shl(3,5,l) = (9.0d0/two)*c1*c2*(three*c2-one)

            shl(1,6,l) =- (9.0d0/two)*c2*(three*c2-one)
            shl(2,6,l) =  (9.0d0/two)*(c3*(three*c2-one)+c2*c3*three-c2*(three*c2-one))
            shl(3,6,l) = (9.0d0/two)*c3*c2*(three*c2-one)

            shl(1,7,l) = (9.0d0/two)*(-c2*(three*c3-one)-c2*c3*three)
            shl(2,7,l) =  (9.0d0/two)*(c3*(three*c3-one)-c2*c3*three -c2*(three*c3-one))
            shl(3,7,l) = (9.0d0/two)*c2*c3*(three*c3-one)

            shl(1,9,l)= (9.0d0/two)*(-c1*(three*c1-one)+c1*c3*three +c3*(three*c1-one))
            shl(2,9,l) =  (9.0d0/two)*(-c1*(three*c1-one))
            shl(3,9,l) = (9.0d0/two)*c3*c1*(three*c1-one)


            shl(1,8,l) = (9.0d0/two)*(-c1*(three*c3-one)-c1*c3*three +c3*(three*c3-one))
            shl(2,8,l) =  (9.0d0/two)*(-c1*c3*three -c1*(three*c3-one))
            shl(3,8,l) = (9.0d0/two)*c1*c3*(three*c3-one)

            shl(1,10,l) = 27.0d0*(c2*c3-c1*c2) 
            shl(2,10,l) = 27.0d0*(c1*c3-c1*c2)
            shl(3,10,l) = 27.0d0*c1*c2*c3
         endif
      200      continue
      
      return
   end subroutine

   
   !======================================================================
   subroutine shgq(xl,det,shl,shg,npint,nel,quad,nen)  
      ! .... program to calculate global derivatives of shape functions and       
      ! jacobian determinants for a  quadrilateral element

      !        xl(j,i)    = global coordinates
      !        det(l)     = jacobian determinant
      !        shl(1,i,l) = local ("xi") derivative of shape function
      !        shl(2,i,l) = local ("eta") derivative of shape function
      !        shl(3,i,l) = local  shape function
      !        shg(1,i,l) = x-derivative of shape function
      !        shg(2,i,l) = y-derivative of shape function
      !        shg(3,i,l) = shl(3,i,l)
      !        xs(i,j)    = jacobian matrix
      !                 i = local node number or global coordinate number
      !                 j = global coordinate number
      !                 l = integration-point number
      !              npint = number of integration points, eq. 1 or 4

      use mLeituraEscrita, only: iecho
      implicit none

      ! .... remove above card for single precision operation
      integer*4:: npint, nel, nen
      real*8  :: xl(2,*), det(*), shl(3,nen,npint), shg(3,nen,npint)
      logical :: quad

      real*8 :: xs(2,2), temp
      integer*4 :: me, i, l, j
      integer*4, parameter :: dois = 2, tres = 3 

      ! call move(shg,shl,3*nen*npint) !copia shl para shg
      shg=shl

      do 700 l=1,npint   
         if (.not.quad) then
            do 100 i=1,3
               shg(i,3,l) = shl(i,3,l) + shl(i,4,l)
               shg(i,4,l) = zero
            100    continue
         endif

         do 300 j=1,2
            do 200 i=1,2
               xs(i,j) = rowdot(shg(i,1,l),xl(j,1),tres,dois,nen)
            200 continue
         300 continue

         det(l) = xs(1,1)*xs(2,2)-xs(1,2)*xs(2,1)
         ! print*,"shg", shg(1,1:4,l)
         ! print*,"xl", xl(1,1:4)
         ! print*, "l=", l, det(l), xs(1,1),xs(2,2),xs(1,2),xs(2,1)
      
         if (det(l).le.zero) then
            write(iecho,1000) nel
            write(*,1000) nel
            stop
         endif

         do 500 j=1,2
            do 400 i=1,2
               xs(i,j) = xs(i,j)/det(l) 
            400 continue
         500 continue

         do 600 i=1,nen
            temp = xs(2,2)*shg(1,i,l) - xs(1,2)*shg(2,i,l)
            shg(2,i,l) = - xs(2,1)*shg(1,i,l) + xs(1,1)*shg(2,i,l)
            shg(1,i,l) = temp
         600 continue

      700 continue

      1000 format(///,'non-positive determinant in element number  ',i10,&
               ' in element group  ',i10)

      return
   end subroutine

   !======================================================================
   subroutine shg3d(xl,det,shl,shg,npint,nel,nen)
      ! .... program to calculate global derivatives of shape functions and
      ! jacobian determinants for a  quadrilateral element

      !        xl(j,i)    = global coordinates
      !        det(l)     = jacobian determinant
      !        shl(1,i,l) = local ("xi") derivative of shape function
      !        shl(2,i,l) = local ("eta") derivative of shape function
      !        shl(3,i,l) = local ("zeta") derivative of shape function
      !        shl(4,i,l) = local  shape function
      !        shg(1,i,l) = x-derivative of shape function
      !        shg(2,i,l) = y-derivative of shape function
      !        shg(3,i,l) = z-derivative of shape funciton
      !        xs(i,j)    = jacobian matrix
      !                 i = local node number or global coordinate number
      !                 j = global coordinate number
      !                 l = integration-point number
      !              npint = number of integration points, eq. 1 or 4

      use mLeituraEscrita, only: iecho

      implicit none

      !.... remove above card for single precision operation
      integer*4 npint, nel, neg, nen
      real*8 xl(3,*), det(*), shl(4,nen,npint), shg(4,nen,npint)

      real*8 xs(3,3)
      integer*4 i, j, l
      real*8 cof11, cof12, cof13, cof21, cof22, cof23, cof31, cof32, cof33
      real*8 temp1, temp2, temp3
      integer*4 me
      integer*4, parameter :: tres = 3, quatro = 4 


      ! common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      ! common /iounit/ iin,iecho,isaid,ipres,iconc

      ! call move(shg,shl,4*nen*npint)
      shg=shl

      do 700 l=1,npint
         
         do 300 j=1,3
            do 200 i=1,3
               xs(i,j) = zero
               xs(i,j) = rowdot(shg(i,1,l),xl(j,1),quatro,tres,nen)
               ! print*, "rowdot", xs(i,j)
            200 continue
         300 continue

         !.. definition of the cofactors
         !.. (recall the definition of matrix inverse : A^-1 = (cof)^T / det A)

         cof11 = xs(2,2)*xs(3,3) - xs(3,2)*xs(2,3)
         cof12 = xs(2,3)*xs(3,1) - xs(2,1)*xs(3,3)
         cof13 = xs(2,1)*xs(3,2) - xs(2,2)*xs(3,1)

         cof21 = xs(3,2)*xs(1,3) - xs(3,3)*xs(1,2)
         cof22 = xs(3,3)*xs(1,1) - xs(1,3)*xs(3,1)
         cof23 = xs(3,1)*xs(1,2) - xs(3,2)*xs(1,1)

         cof31 = xs(1,2)*xs(2,3) - xs(1,3)*xs(2,2)
         cof32 = xs(1,3)*xs(2,1) - xs(1,1)*xs(2,3)
         cof33 = xs(1,1)*xs(2,2) - xs(2,1)*xs(1,2)

         ! print*,  ' em shg3d', cof11, cof12 ,cof13
         ! print*,  ' em shg3d',cof21, cof22, cof23
         ! print*,  ' em shg3d', cof31, cof32, cof33

         det(l) = xs(1,1)*cof11 + xs(1,2)*cof12 + xs(1,3)*cof13
         ! print*, "det(l)=",det(l)

         if (det(l).le.zero) then
            write(iecho,1000) me,nel,neg
            write(*,1000) me,nel,neg
            stop ' em shg3d'
         endif

         do 600 i=1,nen
            temp1=shg(1,i,l)
            temp2=shg(2,i,l)
            temp3=shg(3,i,l)

            shg(3,i,l) =(temp1*cof13+temp2*cof23+temp3*cof33)/det(l)

            shg(2,i,l) =(temp1*cof12+temp2*cof22+temp3*cof32)/det(l)

            shg(1,i,l) =(temp1*cof11+temp2*cof21+temp3*cof31)/det(l)

         600 continue

      700 continue

      1000 format(///,'shg3d, me=',i10,/ &
               'non-positive determinant in element number  ',i10,/,&
               ' in element group  ',i10)

      return
   end subroutine

   !======================================================================
   subroutine shlq(shl,w,npint,nen)
      ! .... program to calculate integration-rule weights, shape functions
      ! and local derivatives for a four-node quadrilateral element

      !               s,t = local element coordinates ("xi", "eta", resp.)       
      !        shl(1,i,l) = local ("xi") derivative of shape function
      !        shl(2,i,l) = local ("eta") derivative of shape function
      !        shl(3,i,l) = local  shape function
      !              w(l) = integration-rule weight
      !                 i = local node number
      !                 l = integration point number
      !              npint = number of integration points, eq. 1 or 4

      implicit none

      !.... remove above card for single precision operation

      integer*4 npint, nen
      real*8 :: shl(3,nen,*),w(*),ra(40),sa(40)

      real*8 r1,w1,r2,w2,r3a,w3a,r3b,w3b,r4a,w4a,r4b,w4b,r5a,w5a,r5b,w5b,w5c
      integer*4 i, l,k
      real*8 w3asq,w3ab,w4asq,w4bsq,w4ab,w5asq,w5bsq,w5ab
      real*8 r, s, onepr,onemr,oneps,onems,onemrsq,onemssq,onep3r,onem3r, onep3s,onem3s
      real*8 f1, f2, f3, f4, f1x, f2x, f3x, f4x
      real*8 g1, g2, g3, g4, g1x, g2x, g3x, g4x

      data r1/0.d00/,w1/2.d00/,&
           r2/0.577350269189626d00/,w2/1.d00/,&
           r3a/0.774596669241483d00/,w3a/0.555555555555556d00/,&
           r3b/0.d00/,w3b/0.888888888888889d00/,&
           r4a/0.861136311594053d00/,w4a/0.347854845137454d00/,&
           r4b/0.339981043584856d00/,w4b/0.652145154862546d00/,&
           r5a/0.906179845938664d00/,w5a/0.236926885056189d00/,&
           r5b/0.538469310105683d00/,w5b/0.478628670499366d00/,&
                                     w5c/0.568888888888888d00/

      if (npint.eq.1) then
         w(1)=w1*w1
         ra(1)=r1
         sa(1)=r1
      end if

      if(npint.eq.4) then
         do 1111 i=1,4
            w(i)=w2
         1111        continue     
         ra(1)=r2
         sa(1)=r2
         ra(2)=-r2
         sa(2)=r2
         ra(3)=-r2
         sa(3)=-r2
         ra(4)=r2
         sa(4)=-r2
      end if

      if(npint.eq.9) then
         w3asq=w3a*w3a
         do 2222 i=1,4
            w(i)=w3asq
         2222        continue     
         ra(1)=r3a
         sa(1)=r3a
         ra(2)=-r3a
         sa(2)=r3a
         ra(3)=-r3a
         sa(3)=-r3a
         ra(4)=r3a
         sa(4)=-r3a

         w3ab=w3a*w3b
         do 3333 i=5,8
            w(i)=w3ab
         3333        continue         
         ra(5)=r3a
         sa(5)=r3b
         ra(6)=r3b
         sa(6)=r3a
         ra(7)=-r3a
         sa(7)=r3b
         ra(8)=r3b
         sa(8)=-r3a

         w(9)=w3b*w3b
         ra(9)=r3b
         sa(9)=r3b
      end if

      if(npint.eq.16) then
         w4asq=w4a*w4a
         do 4444 i=1,4
            w(i)=w4asq
         4444        continue     
         ra(1)=r4a
         sa(1)=r4a
         ra(2)=-r4a
         sa(2)=r4a
         ra(3)=-r4a
         sa(3)=-r4a
         ra(4)=r4a
         sa(4)=-r4a

         w4bsq=w4b*w4b
         do 5555 i=5,8
            w(i)=w4bsq
         5555       continue    
         ra(5)=r4b
         sa(5)=r4b
         ra(6)=-r4b
         sa(6)=r4b
         ra(7)=-r4b
         sa(7)=-r4b
         ra(8)=r4b
         sa(8)=-r4b

         w4ab=w4a*w4b
         do 6666 i=9,16
            w(i)=w4ab
         6666        continue          
         ra(9)=r4b
         sa(9)=r4a
         ra(10)=-r4b
         sa(10)=r4a
         ra(11)=-r4a
         sa(11)=r4b
         ra(12)=-r4a
         sa(12)=-r4b
         ra(13)=-r4b
         sa(13)=-r4a
         ra(14)=r4b
         sa(14)=-r4a
         ra(15)=r4a
         sa(15)=-r4b
         ra(16)=r4a
         sa(16)=r4b
      end if

      if(npint.eq.25) then
         w5asq=w5a*w5a
         do 5544 i=1,4
            w(i  )=w5asq
         5544        continue     
         ra(1)=r5a
         sa(1)=r5a
         ra(2)=-r5a
         sa(2)=r5a
         ra(3)=-r5a
         sa(3)=-r5a
         ra(4)=r5a
         sa(4)=-r5a

         w5bsq=w5b*w5b
         do 5566 i=5,8
            w(i)=w5bsq
         5566       continue    
         ra(5)=r5b
         sa(5)=r5b
         ra(6)=-r5b
         sa(6)=r5b
         ra(7)=-r5b
         sa(7)=-r5b
         ra(8)=r5b
         sa(8)=-r5b

         w5ab=w5a*w5b
         do 6677 i=9,16
            w(i)=w5ab
         6677        continue          
         ra(9)=r5b
         sa(9)=r5a
         ra(10)=-r5b
         sa(10)=r5a
         ra(11)=-r5a
         sa(11)=r5b
         ra(12)=-r5a
         sa(12)=-r5b
         ra(13)=-r5b
         sa(13)=-r5a
         ra(14)=r5b
         sa(14)=-r5a
         ra(15)=r5a
         sa(15)=-r5b
         ra(16)=r5a
         sa(16)=r5b

         do 6688 i=17,20
            w(i)=w5c*w5a
         6688       continue
         ra(17)=zero
         sa(17)=r5a
         ra(18)=-r5a
         sa(18)=zero
         ra(19)=zero
         sa(19)=-r5a
         ra(20)=r5a
         sa(20)=zero

         do 6699 i=21,24
            w(i)=w5c*w5b
         6699       continue
         ra(21)=zero
         sa(21)=r5b
         ra(22)=-r5b
         sa(22)=zero
         ra(23)=zero
         sa(23)=-r5b
         ra(24)=r5b
         sa(24)=zero

         w(25)=w5c*w5c
         ra(25)=zero
         sa(25)=zero
      end if

      do 200 l=1,npint

         r=ra(l)
         s=sa(l)
         onepr=one+r
         onemr=one-r
         oneps=one+s
         onems=one-s
         shl(1,1,l)=oneps*pt25
         shl(2,1,l)=onepr*pt25
         shl(3,1,l)=onepr*oneps*pt25
         shl(1,2,l)=-oneps*pt25
         shl(2,2,l)=onemr*pt25
         shl(3,2,l)=onemr*oneps*pt25
         shl(1,3,l)=-onems*pt25
         shl(2,3,l)=-onemr*pt25
         shl(3,3,l)=onemr*onems*pt25
         shl(1,4,l)=onems*pt25
         shl(2,4,l)=-onepr*pt25
         shl(3,4,l)=onepr*onems*pt25
         if(nen.eq.9) then
            onemrsq=one-r*r
            onemssq=one-s*s
            shl(1,5,l)=-r*oneps
            shl(2,5,l)=onemrsq*pt5
            shl(3,5,l)=onemrsq*oneps*pt5
            shl(1,6,l)=-onemssq*pt5
            shl(2,6,l)=-s*onemr
            shl(3,6,l)=onemssq*onemr*pt5
            shl(1,7,l)=-r*onems
            shl(2,7,l)=-onemrsq*pt5
            shl(3,7,l)=onemrsq*onems*pt5
            shl(1,8,l)=onemssq*pt5
            shl(2,8,l)=-s*onepr
            shl(3,8,l)=onemssq*onepr*pt5
            shl(1,9,l)=-two*r*onemssq
            shl(2,9,l)=-two*s*onemrsq
            shl(3,9,l)=onemrsq*onemssq

            do 7777 k=5,8
               do 8888 i=1,3
                  shl(i,k,l)=shl(i,k,l)-pt5*shl(i,9,l)
               8888                    continue  
            7777              continue 

            do 9999 i=1,3
               shl(i,1,l)=shl(i,1,l)-pt5*(shl(i,5,l)+shl(i,8,l))-pt25*shl(i,9,l)      
               shl(i,2,l)=shl(i,2,l)-pt5*(shl(i,6,l)+shl(i,5,l))-pt25*shl(i,9,l)
               shl(i,3,l)=shl(i,3,l)-pt5*(shl(i,7,l)+shl(i,6,l))-pt25*shl(i,9,l)
               shl(i,4,l)=shl(i,4,l)-pt5*(shl(i,8,l)+shl(i,7,l))-pt25*shl(i,9,l)
            9999               continue   
         end if
         if (nen.eq.16) then
            onemrsq=one-r*r
            onemssq=one-s*s
            onep3r=one+three*r
            onem3r=one-three*r
            onep3s=one+three*s
            onem3s=one-three*s
            f1=-1.d00/16.d00*(9.d00*r*r-1.d00)*(r-1.d00)
            f2=9.d00/16.d00*(1.d00-r*r)*onem3r
            f3=9.d00/16.d00*(1.d00-r*r)*onep3r
            f4=1.d00/16.d00*(9.d00*r*r-1.d00)*(r+1.d00)

            f1x=-1.d00/16.d00*(18.d00*r)*(r-1.d00)-1.d00/16.d00*(9.d00*r*r-1.d00)
            f2x=9.d00/16.d00*(-2.d00*r)*onem3r-9.d00/16.d00*(1.d00-r*r)*3.d00
            f3x=9.d00/16.d00*(-2.d00*r)*onep3r-9.d00/16.d00*(r*r-1.d00)*3.d00
            f4x=1.d00/16.d00*(18.d00*r)*(r+1.d00)+1.d00/16.d00*(9.d00*r*r-1.d00)

            g1=-1.d00/16.d00*(9.d00*s*s-1.d00)*(s-1.d00)
            g2=9.d00/16.d00*(1.d00-s*s)*onem3s
            g3=9.d00/16.d00*(1.d00-s*s)*onep3s
            g4=1.d00/16.d00*(9.d00*s*s-1.d00)*(s+1.d00)

            g1x=-1.d00/16.d00*(18.d00*s)*(s-1.d00)-1.d00/16.d00*(9.d00*s*s-1.d00)
            g2x=9.d00/16.d00*(-2.d00*s)*onem3s-9.d00/16.d00*(1.d00-s*s)*3.d00
            g3x=9.d00/16.d00*(-2.d00*s)*onep3s-9.d00/16.d00*(s*s-1.d00)*3.d00
            g4x=1.d00/16.d00*(18.d00*s)*(s+1.d00)+1.d00/16.d00*(9.d00*s*s-1.d00)

            shl(3,1,l)=f1*g1
            shl(3,2,l)=f4*g1
            shl(3,3,l)=f4*g4
            shl(3,4,l)=f1*g4
	     
            shl(3,5,l)=f2*g1
            shl(3,6,l)=f3*g1
            shl(3,7,l)=f4*g2
            shl(3,8,l)=f4*g3
            shl(3,9,l)=f3*g4
            shl(3,10,l)=f2*g4
            shl(3,11,l)=f1*g3
            shl(3,12,l)=f1*g2
            shl(3,13,l)=f2*g2
            shl(3,14,l)=f3*g2
            shl(3,15,l)=f3*g3
            shl(3,16,l)=f2*g3

            shl(1,1,l)=f1x*g1
            shl(1,2,l)=f4x*g1
            shl(1,3,l)=f4x*g4
            shl(1,4,l)=f1x*g4
	     
            shl(1,5,l)=f2x*g1
            shl(1,6,l)=f3x*g1
            shl(1,7,l)=f4x*g2
            shl(1,8,l)=f4x*g3
            shl(1,9,l)=f3x*g4
            shl(1,10,l)=f2x*g4
            shl(1,11,l)=f1x*g3
            shl(1,12,l)=f1x*g2
            shl(1,13,l)=f2x*g2
            shl(1,14,l)=f3x*g2
            shl(1,15,l)=f3x*g3
            shl(1,16,l)=f2x*g3

            shl(2,1,l)=f1*g1x
            shl(2,2,l)=f4*g1x
            shl(2,3,l)=f4*g4x
            shl(2,4,l)=f1*g4x
	     
            shl(2,5,l)=f2*g1x
            shl(2,6,l)=f3*g1x
            shl(2,7,l)=f4*g2x
            shl(2,8,l)=f4*g3x
            shl(2,9,l)=f3*g4x
            shl(2,10,l)=f2*g4x
            shl(2,11,l)=f1*g3x
            shl(2,12,l)=f1*g2x
            shl(2,13,l)=f2*g2x
            shl(2,14,l)=f3*g2x
            shl(2,15,l)=f3*g3x
            shl(2,16,l)=f2*g3x

         end if
      200 continue
      
      return
   end subroutine

   !======================================================================
   subroutine shg2q(xl,shl,shl2,shg2,npint,nel,neg,quad,nen) 
      ! .... program to calculate global second derivatives of shape functions
      ! for a  quadrilateral element

      !        xl(j,i)    = global coordinates
      !        shl(1,i,l) = local ("xi") derivative of shape function
      !        shl(2,i,l) = local ("eta") derivative of shape function
      !        shl(3,i,l) = local  shape function
      !        shl2(1,i,l) = local second ("xi") derivative of shape function
      !        shl2(2,i,l) = local second ("eta") derivative of shape function
      !        shl2(3,i,l) = local second ("xi*eta") derivative of shape function
      !        shg2(1,i,l) = second x-derivative of shape function
      !        shg2(2,i,l) = second y-derivative of shape function
      !        shg2(3,i,l) = second xy-derivative of shape function
      !        xs(i,j)    = jacobian matrix
      !                 i = local node number or global coordinate number
      !                 j = global coordinate number
      !                 l = integration-point number
      !              npint = number of integration points, eq. 1 or 4

      use mLeituraEscrita, only: iecho

      implicit none

      !.... remove above card for single precision operation

      integer*4 npint, nel, neg, nen
      real*8 :: xl(2,*),shl(3,nen,*),shl2(3,nen,*),shg2(3,nen,*)
      logical quad

      real*8 :: xs(2,2),t1(3,2),t2(3,3),c1(3,2), det
      real*8 :: dxirx, dxiry, detarx, detary
      integer*4:: i,l,j,k
      integer*4, parameter :: dois = 2, tres = 3


      do 960 l=1,npint

         if (.not.quad) then
            do 100 i=1,3
               shl(i,3,l) = shl(i,3,l) + shl(i,4,l)
               shl(i,4,l) = zero
               shl2(i,3,l) = shl2(i,3,l) + shl2(i,4,l)
               shl2(i,4,l) = zero
            100    continue
         endif

         do 300 j=1,2
            do 200 i=1,2
               xs(i,j) = rowdot(shl(i,1,l),xl(j,1),tres,dois,nen)
            200 continue
         300 continue

         det = xs(1,1)*xs(2,2)-xs(1,2)*xs(2,1)
         if (det.le.zero) then
            write(iecho,1000) nel,neg
            write(*,1000) nel,neg
            stop
         endif

         do 500 j=1,2
            do 400 i=1,2
               xs(i,j) = xs(i,j)/det
            400 continue
         500 continue

         ! jacobian elements
         dxirx=xs(2,2)
         dxiry=xs(1,2)
         detarx=xs(2,1)
         detary=xs(1,1)
         xs(1,1)=dxirx
         xs(1,2)=detarx
         xs(2,1)=dxiry
         xs(2,2)=detary

         ! calculation of global second order derivatives
     
         ! {d2global}=[t1]{dlocal}+[t2]{d2local}
         ! [t1]=-[t2][c1][xs]

         ! .... form t2
         t2(1,1)=xs(1,1)**2
         t2(1,2)=xs(1,2)**2
         t2(1,3)=two*xs(1,1)*xs(1,2)
         t2(2,1)=xs(2,1)**2
         t2(2,2)=xs(2,2)**2
         t2(2,3)=two*xs(2,1)*xs(2,2)
         t2(3,1)=xs(1,1)*xs(2,1)
         t2(3,2)=xs(1,2)*xs(2,2)
         t2(3,3)=xs(1,1)*xs(2,2)+xs(1,2)*xs(2,1)

         ! .... form c1
         do i=1,3
            do j=1,2
               c1(i,j) = rowdot(shl2(i,1,l),xl(j,1),tres,dois,nen)
            enddo
         enddo

         ! .... form t1
         do i=1,3
            do j=1,2
               do k=1,3
                  t1(i,j)=t1(i,j)-t2(i,k)*(c1(k,1)*xs(1,j)+c1(k,2)*xs(2,j))
               enddo
            enddo
         enddo

         ! .... transformation from natural coor. to global coor.
         do j=1,nen
            do i=1,3
               do k=1,2
                  shg2(i,j,l)=shg2(i,j,l)+t1(i,k)*shl(k,j,l)
               enddo
               do k=1,3
                  shg2(i,j,l)=shg2(i,j,l)+t2(i,k)*shl2(k,j,l)
               enddo
            enddo
         enddo
      960 continue

      1000 format(///,'non-positive determinant in element number  ',i10,&
               ' in element group  ',i10)

      return
   end subroutine

   !======================================================================
   subroutine shl2q(shl2,npint,nen)
      ! .... program to calculate local second  derivatives 
      ! for a four-node quadrilateral element

      !               s,t = local element coordinates ("xi", "eta", resp.)       
      !        shl2(1,i,l) = local second ("xi") derivative of shape function
      !        shl2(2,i,l) = local second ("eta") derivative of shape function
      !        shl2(3,i,l) = local second ("xi*eta") derivative of shape function
      !                 i = local node number
      !                 l = integration point number
      !              npint = number of integration points, eq. 1 or 4
      
      implicit none

      ! .... remove above card for single precision operation

      integer*4:: npint, nen
      real*8  :: shl2(3,nen,*)

      real*8  :: ra(16),sa(16)
      real*8  :: r1,r2,r3a,r3b,r4a,r4b
      integer*4:: l,k,i
      real*8  :: r,s,onepr,onemr,onems,oneps,onemssq,onemrsq,onep3r,onem3r,onep3s,onem3s

      data r1/0.d00/,&
           r2/0.577350269189626d00/,&
           r3a/0.774596669241483d00/,&
           r3b/0.d00/,&
           r4a/0.861136311594053d00/,&
           r4b/0.339981043584856d00/

      if (npint.eq.1) then
         ra(1)=r1
         sa(1)=r1
      end if

      if(npint.eq.4) then
         ra(1)=r2
         sa(1)=r2
         ra(2)=-r2
         sa(2)=r2
         ra(3)=-r2
         sa(3)=-r2
         ra(4)=r2
         sa(4)=-r2
      end if

      if(npint.eq.9) then
         ra(1)=r3a
         sa(1)=r3a
         ra(2)=-r3a
         sa(2)=r3a
         ra(3)=-r3a
         sa(3)=-r3a
         ra(4)=r3a
         sa(4)=-r3a

         ra(5)=r3a
         sa(5)=r3b
         ra(6)=r3b
         sa(6)=r3a
         ra(7)=-r3a
         sa(7)=r3b
         ra(8)=r3b
         sa(8)=-r3a

         ra(9)=r3b
         sa(9)=r3b
      end if

      if(npint.eq.16) then
         ra(1)=r4a
         sa(1)=r4a
         ra(2)=-r4a
         sa(2)=r4a
         ra(3)=-r4a
         sa(3)=-r4a
         ra(4)=r4a
         sa(4)=-r4a

         ra(5)=r4b
         sa(5)=r4b
         ra(6)=-r4b
         sa(6)=r4b
         ra(7)=-r4b
         sa(7)=-r4b
         ra(8)=r4b
         sa(8)=-r4b

         ra(9)=r4b
         sa(9)=r4a
         ra(10)=-r4b
         sa(10)=r4a
         ra(11)=-r4a
         sa(11)=r4b
         ra(12)=-r4a
         sa(12)=-r4b
         ra(13)=-r4b
         sa(13)=-r4a
         ra(14)=r4b
         sa(14)=-r4a
         ra(15)=r4a
         sa(15)=-r4b
         ra(16)=r4a
         sa(16)=r4b
      end if

      do 200 l=1,npint
         r=ra(l)
         s=sa(l)
         shl2(1,1,l)=zero
         shl2(2,1,l)=zero
         shl2(3,1,l)=pt25
         shl2(1,2,l)=zero
         shl2(2,2,l)=zero
         shl2(3,2,l)=-pt25
         shl2(1,3,l)=zero
         shl2(2,3,l)=zero
         shl2(3,3,l)=pt25
         shl2(1,4,l)=zero
         shl2(2,4,l)=zero
         shl2(3,4,l)=-pt25
         if(nen.eq.9) then
            onepr=one+r
            onemr=one-r
            oneps=one+s
            onems=one-s
            shl2(1,5,l)=-oneps
            shl2(2,5,l)=zero
            shl2(3,5,l)=-r
            shl2(1,6,l)=zero
            shl2(2,6,l)=-onemr
            shl2(3,6,l)=s
            shl2(1,7,l)=-onems
            shl2(2,7,l)=zero
            shl2(3,7,l)=r
            shl2(1,8,l)=zero
            shl2(2,8,l)=-onepr
            shl2(3,8,l)=-s
            shl2(1,9,l)=-two*onemssq
            shl2(2,9,l)=-two*onemrsq
            shl2(3,9,l)=four*r*s

            do 1111 k=5,8
               do 2222 i=1,3
                  shl2(i,k,l)=shl2(i,k,l)-pt5*shl2(i,9,l)
               2222                    continue
            1111              continue

            do 3333 i=1,3
               shl2(i,1,l)=shl2(i,1,l)-pt5*(shl2(i,5,l)+shl2(i,8,l))-pt25*shl2(i,9,l)         
               shl2(i,2,l)=shl2(i,2,l)-pt5*(shl2(i,6,l)+shl2(i,5,l))-pt25*shl2(i,9,l)
               shl2(i,3,l)=shl2(i,3,l)-pt5*(shl2(i,7,l)+shl2(i,6,l))-pt25*shl2(i,9,l)
               shl2(i,4,l)=shl2(i,4,l)-pt5*(shl2(i,8,l)+shl2(i,7,l))-pt25*shl2(i,9,l)
   
            3333              continue
         end if
         if (nen.eq.16) then
            onemrsq=one-r*r
            onemssq=one-s*s
            onep3r=one+three*r
            onem3r=one-three*r
            onep3s=one+three*s
            onem3s=one-three*s

            ! not inplemented
         end if
      200 continue

      return
   end subroutine


   ! 3D
   !======================================================================
   subroutine shlq3d (shl,w,npint,nen)
      ! .... program to calculate integration-rule weights, shape functions
      ! and local derivatives for a eight-node tridimensional element

      !              s,t,u = local element coordinates ("xi", "eta"," resp.)       
      !         shl(1,i,l) = local ("xi") derivative of shape function
      !        shl(2,i,l) = local ("eta") derivative of shape function
      !        shl(3,i,l) = local (" zeta") derivative of shape function
      !        shl(4,i,l) = local  shape function
      !              w(l) = integration-rule weight
      !                 i = local node number
      !                 l = integration point number
      !              npint = number of integration points, eq. 1 or 4

      implicit none

      !.... remove above card for single precision operation

      integer*4:: npint, nen
      real*8 :: shl(4,nen,*),w(*)

      real*8 :: ra(27),sa(27),ta(27)
      real*8 :: r2,r3a,w3a,r3b,w3b,w1
      integer :: i,j
      real*8  :: eight

      data  r2/0.577350269189626d00/, &
               &r3a/0.774596669241483d00/,w3a/0.555555555555556d00/,&
               &r3b/0.d00/,w3b/0.888888888888889d00/

      eight = five+three

      ! inserido por tuane
      if (npint.eq.1) then
         ra(1) = zero
         sa(1) = zero
         ta(1) = zero
         w1=2.d0
         w(1) = w1*w1*w1

      endif

      if (npint.eq.6) then

         ra(1) = one
         sa(1) = zero
         ta(1) = zero

         ra(2) = -one
         sa(2) = zero
         ta(2) = zero

         ra(3) = zero
         sa(3) = one
         ta(3) = zero

         ra(4) = zero
         sa(4) = -one
         ta(4) = zero

         ra(5) = zero
         sa(5) = zero
         ta(5) = one

         ra(6) = zero
         sa(6) = zero
         ta(6) = -one

         do 10 i=1,6
            w(i) = four/three

         10   continue
      endif

      if (npint.eq.8) then

         ! ra(1) =-r2
         ! sa(1) =-r2
         ! ta(1) =-r2

         ! ra(2) = r2
         ! sa(2) =-r2
         ! ta(2) =-r2

         ! ra(3) = r2
         ! sa(3) = r2
         ! ta(3) =-r2

         ! ra(4) =-r2
         ! sa(4) = r2
         ! ta(4) =-r2

         ! ra(5) = -r2
         ! sa(5) = -r2
         ! ta(5) = -r2

         ! ra(6) = r2
         ! sa(6) = -r2
         ! ta(6) = r2
   

         ! ra(7) = r2
         ! sa(7) = r2
         ! ta(7) = r2

         ! ra(8) = -r2
         ! sa(8) = r2
         ! ta(8) = r2

         ! original
         ra(1) = -r2
         sa(1) = -r2
         ta(1) = -r2

         ra(2) = r2
         sa(2) =-r2
         ta(2) =-r2

         ra(3) =-r2
         sa(3) = r2
         ta(3) =-r2

         ra(4) = r2
         sa(4) = r2
         ta(4) =-r2

         ra(5) =-r2
         sa(5) =-r2
         ta(5) = r2

         ra(6) = r2
         sa(6) =-r2
         ta(6) = r2
   
         ra(7) =-r2
         sa(7) = r2
         ta(7) = r2

         ra(8) = r2
         sa(8) = r2
         ta(8) = r2

         do 15 i=1,8
            w(i) = one
         15   continue
      endif

      if (npint.eq.27) then
         ra(1) = -r3a
         sa(1) = -r3a
         ta(1) = -r3a
         w(1) = w3a*w3a*w3a

         ra(2) =  r3b
         sa(2) = -r3a
         ta(2) = -r3a
         w(2) = w3b*w3a*w3a

         ra(3) = r3a
         sa(3) = -r3a
         ta(3) = -r3a
         w(3) = w(1)

         ra(4) = -r3a
         sa(4) = r3b
         ta(4) = -r3a
         w(4) = w(2)

         ra(5) = r3b
         sa(5) = r3b
         ta(5) = -r3a
         w(5) = w3b*w3b*w3a

         ra(6) = r3a
         sa(6) = r3b
         ta(6) = -r3a
         w(6) = w(2)

         ra(7) = -r3a
         sa(7) = r3a
         ta(7) = -r3a
         w(7) = w(1)

         ra(8) = r3b
         sa(8) = r3a
         ta(8) = -r3a
         w(8) = w(2)

         ra(9) = r3a
         sa(9) = r3a
         ta(9) = -r3a
         w(9) = w(1)

         ra(10) = -r3a
         sa(10) = -r3a
         ta(10) =  r3b
         w(10) = w(2)

         ra(11) =  r3b
         sa(11) = -r3a
         ta(11) = - r3b
         w(11) =w3b*w3a*w3b

         ra(12) = r3a
         sa(12) = -r3a
         ta(12) = r3b
         w(12) = w(2)

         ra(13) = -r3a
         sa(13) = r3b
         ta(13) = r3b
         w(13) = w3a*w3b*w3b

         ra(14) = r3b
         sa(14) = r3b
         ta(14) = r3b
         w(14) = w3b*w3b*w3b

         ra(15) = r3a
         sa(15) = r3b
         ta(15) = r3b
         w(15) = w(13)

         ra(16) = -r3a
         sa(16) = r3a
         ta(16) = r3b
         w(16) = w(2)

         ra(17) = r3b
         sa(17) = r3a
         ta(17) = r3b
         w(17) = w(13)

         ra(18) = r3a
         sa(18) = r3a
         ta(18) = r3b
         w(18) = w(2)

         ra(19) = -r3a
         sa(19) = -r3a
         ta(19) =  r3a
         w(19) = w(1)

         ra(20) =  r3b
         sa(20) = -r3a
         ta(20) = r3a
         w(20) = w(2)

         ra(21) = r3a
         sa(21) = -r3a
         ta(21) = r3a
         w(21) = w(1)

         ra(22) = -r3a
         sa(22) = r3b
         ta(22) = r3a
         w(22) = w(2)

         ra(23) = r3b
         sa(23) = r3b
         ta(23) = r3a
         w(23) = w(13)

         ra(24) = r3a
         sa(24) = r3b
         ta(24) = r3a
         w(24) = w(2)

         ra(25) = -r3a
         sa(25) = r3a
         ta(25) = r3a
         w(25) = w(1)

         ra(26) = r3b
         sa(26) = r3a
         ta(26) = r3a
         w(26) = w(2)

         ra(27) = r3a
         sa(27) = r3a
         ta(27) = r3a
         w(27) = w(1)
      endif

      if (nen.eq.8) then
         do 20 j=1,npint
            shl(1,1,j) = -(one - sa(j))*(one - ta(j))/eight
            shl(2,1,j) = -(one - ra(j))*(one - ta(j))/eight
            shl(3,1,j) = -(one - ra(j))*(one - sa(j))/eight
            shl(4,1,j) =  (one - ra(j))*(one - sa(j))*(one - ta(j))/eight

            shl(1,2,j) =  (one - sa(j))*(one - ta(j))/eight           
            shl(2,2,j) = -(one + ra(j))*(one - ta(j))/eight           
            shl(3,2,j) = -(one + ra(j))*(one - sa(j))/eight
            shl(4,2,j) =  (one + ra(j))*(one - sa(j))*(one - ta(j))/eight

            shl(1,3,j) =  (one + sa(j))*(one - ta(j))/eight           
            shl(2,3,j) =  (one + ra(j))*(one - ta(j))/eight
            shl(3,3,j) = -(one + ra(j))*(one + sa(j))/eight
            shl(4,3,j) =  (one + ra(j))*(one + sa(j))*(one - ta(j))/eight           

            shl(1,4,j) = -(one + sa(j))*(one - ta(j))/eight           
            shl(2,4,j) =  (one - ra(j))*(one - ta(j))/eight
            shl(3,4,j) = -(one - ra(j))*(one + sa(j))/eight
            shl(4,4,j) =  (one - ra(j))*(one + sa(j))*(one - ta(j))/eight

            shl(1,5,j) = -(one - sa(j))*(one + ta(j))/eight           
            shl(2,5,j) = -(one - ra(j))*(one + ta(j))/eight
            shl(3,5,j) =  (one - ra(j))*(one - sa(j))/eight
            shl(4,5,j) =  (one - ra(j))*(one - sa(j))*(one + ta(j))/eight
           
            shl(1,6,j) =  (one - sa(j))*(one + ta(j))/eight           
            shl(2,6,j) = -(one + ra(j))*(one + ta(j))/eight
            shl(3,6,j) =  (one + ra(j))*(one - sa(j))/eight           
            shl(4,6,j) =  (one + ra(j))*(one - sa(j))*(one + ta(j))/eight
           
            shl(1,7,j) =  (one + sa(j))*(one + ta(j))/eight
            shl(2,7,j) =  (one + ra(j))*(one + ta(j))/eight           
            shl(3,7,j) =  (one + ra(j))*(one + sa(j))/eight           
            shl(4,7,j) =  (one + ra(j))*(one + sa(j))*(one + ta(j))/eight
           
            shl(1,8,j) = -(one + sa(j))*(one + ta(j))/eight           
            shl(2,8,j) =  (one - ra(j))*(one + ta(j))/eight
            shl(3,8,j) =  (one - ra(j))*(one + sa(j))/eight
            shl(4,8,j) =  (one - ra(j))*(one + sa(j))*(one + ta(j))/eight
      
            ! shl(4,1,j) = (one - ra(j))*(one - sa(j))*(one - ta(j))/eight
            ! shl(4,2,j) = (one + ra(j))*(one - sa(j))*(one - ta(j))/eight
            ! shl(4,3,j) = (one + ra(j))*(one + sa(j))*(one - ta(j))/eight
            ! shl(4,4,j) = (one - ra(j))*(one + sa(j))*(one - ta(j))/eight
            ! shl(4,5,j) = (one - ra(j))*(one - sa(j))*(one + ta(j))/eight
            ! shl(4,6,j) = (one + ra(j))*(one - sa(j))*(one + ta(j))/eight
            ! shl(4,7,j) = (one + ra(j))*(one + sa(j))*(one + ta(j))/eight
            ! shl(4,8,j) = (one - ra(j))*(one + sa(j))*(one + ta(j))/eight

            ! shl(3,1,j) = -(one - ra(j))*(one - sa(j))/eight
            ! shl(3,2,j) = -(one + ra(j))*(one - sa(j))/eight
            ! shl(3,3,j) = -(one + ra(j))*(one + sa(j))/eight
            ! shl(3,4,j) = -(one - ra(j))*(one + sa(j))/eight
            ! shl(3,5,j) = (one - ra(j))*(one - sa(j))/eight
            ! shl(3,6,j) = (one + ra(j))*(one - sa(j))/eight
            ! shl(3,7,j) = (one + ra(j))*(one + sa(j))/eight
            ! shl(3,8,j) = (one - ra(j))*(one + sa(j))/eight

            ! shl(2,1,j) = -(one - ra(j))*(one - ta(j))/eight
            ! shl(2,2,j) = -(one + ra(j))*(one - ta(j))/eight
            ! shl(2,3,j) = (one + ra(j))*(one - ta(j))/eight
            ! shl(2,4,j) = (one - ra(j))*(one - ta(j))/eight
            ! shl(2,5,j) = -(one - ra(j))*(one + ta(j))/eight
            ! shl(2,6,j) = -(one + ra(j))*(one + ta(j))/eight
            ! shl(2,7,j) = (one + ra(j))*(one + ta(j))/eight
            ! shl(2,8,j) = (one - ra(j))*(one + ta(j))/eight

            ! shl(1,1,j) = -(one - sa(j))*(one - ta(j))/eight
            ! shl(1,2,j) = (one - sa(j))*(one - ta(j))/eight
            ! shl(1,3,j) = (one + sa(j))*(one - ta(j))/eight
            ! shl(1,4,j) = -(one + sa(j))*(one - ta(j))/eight
            ! shl(1,5,j) = -(one - sa(j))*(one + ta(j))/eight
            ! shl(1,6,j) = (one - sa(j))*(one + ta(j))/eight
            ! shl(1,7,j) = (one + sa(j))*(one + ta(j))/eight
            ! shl(1,8,j) = -(one + sa(j))*(one + ta(j))/eight

         20   continue
      else
         if (nen.eq.27) then
            do 30 j=1,npint
               shl(4,1,j) = (ra(j)-one)*(sa(j)-one)*(ta(j)-one)&
                        *ra(j)*sa(j)*ta(j)/eight
               shl(4,2,j) = (ra(j)+one)*(sa(j)-one)*(ta(j)-one)&
                        *ra(j)*sa(j)*ta(j)/eight
               shl(4,3,j) = (ra(j)+one)*(sa(j)+one)*(ta(j)-one)&
                        *ra(j)*sa(j)*ta(j)/eight
               shl(4,4,j) = (ra(j)-one)*(sa(j)+one)*(ta(j)-one)&
                        *ra(j)*sa(j)*ta(j)/eight
               shl(4,5,j) = (ra(j)-one)*(sa(j)-one)*(ta(j)+one)&
                        *ra(j)*sa(j)*ta(j)/eight
               shl(4,6,j) = (ra(j)+one)*(sa(j)-one)*(ta(j)+one)&
                        *ra(j)*sa(j)*ta(j)/eight
               shl(4,7,j) = (ra(j)+one)*(sa(j)+one)*(ta(j)+one)&
                        *ra(j)*sa(j)*ta(j)/eight
               shl(4,8,j) = (ra(j)-one)*(sa(j)+one)*(ta(j)+one)&
                        *ra(j)*sa(j)*ta(j)/eight
               shl(4,9,j) = (one - ra(j)**2)*(sa(j)-one)*(ta(j)-one)&
                        *sa(j)*ta(j)/four
               shl(4,10,j) = (ra(j)+one)*(one - sa(j)**2)*(ta(j)-one)&
                        *ra(j)*ta(j)/four
               shl(4,11,j) = (one - ra(j)**2)*(sa(j)+one)*(ta(j)-one)&
                        *sa(j)*ta(j)/four
               shl(4,12,j) = (ra(j)-one)*(one - sa(j)**2)*(ta(j)-one)&
                        *ra(j)*ta(j)/four
               shl(4,13,j) = (one - ra(j)**2)*(sa(j)-one)*(ta(j)+one)&
                        *sa(j)*ta(j)/four
               shl(4,14,j) =  (ra(j)+one)*(one - sa(j)**2)*(ta(j)+one)&
                        *ra(j)*ta(j)/four
               shl(4,15,j) =(one - ra(j)**2)*(sa(j)+one)*(ta(j)+one)&
                        *sa(j)*ta(j)/four
               shl(4,16,j) = (ra(j)-one)*(one - sa(j)**2)*(ta(j)+one)&
                        *ra(j)*ta(j)/four
               shl(4,17,j) =  (ra(j)-one)*(sa(j)-one)*(one-ta(j)**2)&
                        *ra(j)*sa(j)/four
               shl(4,18,j) =  (ra(j)+one)*(sa(j)-one)*(one-ta(j)**2)&
                        *ra(j)*sa(j)/four
               shl(4,19,j) =  (ra(j)+one)*(sa(j)+one)*(one-ta(j)**2)&
                        *ra(j)*sa(j)/four
               shl(4,20,j) =(ra(j)-one)*(sa(j)+one)*(one-ta(j)**2)&
                        *ra(j)*sa(j)/four
               shl(4,21,j) = (one-ra(j)**2)*(one-sa(j)**2)*(ta(j)-one)&
                        *ta(j)/two
               shl(4,22,j) =  (one-ra(j)**2)*(one-sa(j)**2)*(ta(j)+one)&
                        *ta(j)/two
               shl(4,23,j) =  (one-ra(j)**2)*(sa(j)-one)*(one-ta(j)**2)&
                        *sa(j)/two
               shl(4,24,j) =  (one-ra(j)**2)*(sa(j)+one)*(one-ta(j)**2)&
                        *sa(j)/two
               shl(4,25,j) = (ra(j)-one)* (one-sa(j)**2)*(one-ta(j)**2)&
                        *ra(j)/two
               shl(4,26,j) =  (ra(j)+one)* (one-sa(j)**2)*(one-ta(j)**2)&
                        *ra(j)/two
               shl(4,27,j) =(one - ra(j)**2)*(one - sa(j)**2)*(one - ta(j)**2)

               shl(3,1,j) = ((ra(j)-one)*(sa(j)-one)*ra(j)*sa(j)*ta(j)&
                        +(ra(j)-one)*(sa(j)-one)*(ta(j)-one)*ra(j)*sa(j))/eight
               shl(3,2,j) = ((ra(j)+one)*(sa(j)-one)*ra(j)*sa(j)*ta(j)&
                        + (ra(j)+one)*(sa(j)-one)*(ta(j)-one)*ra(j)*sa(j))/eight
               shl(3,3,j) = ((ra(j)+one)*(sa(j)+one)*ra(j)*sa(j)*ta(j)&
                        +(ra(j)+one)*(sa(j)+one)*(ta(j)-one)*ra(j)*sa(j))/eight
               shl(3,4,j) = ((ra(j)-one)*(sa(j)+one)*ra(j)*sa(j)*ta(j)&
                        +(ra(j)-one)*(sa(j)+one)*(ta(j)-one)*ra(j)*sa(j))/eight
               shl(3,5,j) = ((ra(j)-one)*(sa(j)-one)*ra(j)*sa(j)*ta(j)&
                        + (ra(j)-one)*(sa(j)-one)*(ta(j)+one)*ra(j)*sa(j))/eight
               shl(3,6,j) = ((ra(j)+one)*(sa(j)-one)*ra(j)*sa(j)*ta(j)&
                        +(ra(j)+one)*(sa(j)-one)*(ta(j)+one)*ra(j)*sa(j))/eight
               shl(3,7,j) = ((ra(j)+one)*(sa(j)+one)*ra(j)*sa(j)*ta(j)&
                        +(ra(j)+one)*(sa(j)+one)*(ta(j)+one)*ra(j)*sa(j))/eight
               shl(3,8,j) = ((ra(j)-one)*(sa(j)+one)*ra(j)*sa(j)*ta(j)&
                        +(ra(j)-one)*(sa(j)+one)*(ta(j)+one)*ra(j)*sa(j))/eight
               shl(3,9,j) =( (one - ra(j)**2)*(sa(j)-one)*sa(j)*ta(j)&
                        + (one - ra(j)**2)*(sa(j)-one)*(ta(j)-one)*sa(j))/four
               shl(3,10,j) = ((ra(j)+one)*(one - sa(j)**2)*ra(j)*ta(j)&
                        + (ra(j)+one)*(one - sa(j)**2)*(ta(j)-one)*ra(j))/four
               shl(3,11,j) =( (one - ra(j)**2)*(sa(j)+one)*sa(j)*ta(j)&
                        + (one - ra(j)**2)*(sa(j)+one)*(ta(j)-one)*sa(j))/four
               shl(3,12,j) = ((ra(j)-one)*(one - sa(j)**2)*ra(j)*ta(j)&
                        +(ra(j)-one)*(one - sa(j)**2)*(ta(j)-one)*ra(j))/four
               shl(3,13,j) = ((one - ra(j)**2)*(sa(j)-one)*sa(j)*ta(j)&
                        + (one - ra(j)**2)*(sa(j)-one)*(ta(j)+one)*sa(j))/four
               shl(3,14,j) = ( (ra(j)+one)*(one - sa(j)**2)*ra(j)*ta(j)&
                        + (ra(j)+one)*(one - sa(j)**2)*(ta(j)+one)*ra(j))/four
               shl(3,15,j) =((one - ra(j)**2)*(sa(j)+one)*sa(j)*ta(j)&
                        +(one - ra(j)**2)*(sa(j)+one)*(ta(j)+one)*sa(j))/four
               shl(3,16,j) = ((ra(j)-one)*(one - sa(j)**2)*ra(j)*ta(j)&
                        +(ra(j)-one)*(one - sa(j)**2)*(ta(j)+one)*ra(j))/four
               shl(3,17,j) =  (ra(j)-one)*(sa(j)-one)*(-two*ta(j))&
                        *ra(j)*sa(j)/four
               shl(3,18,j) =  (ra(j)+one)*(sa(j)-one)*(-two*ta(j))&
                        *ra(j)*sa(j)/four
               shl(3,19,j) =  (ra(j)+one)*(sa(j)+one)*(-two*ta(j))&
                        *ra(j)*sa(j)/four
               shl(3,20,j) =(ra(j)-one)*(sa(j)+one)*(-two*ta(j))&
                        *ra(j)*sa(j)/four
               shl(3,21,j) = ((one-ra(j)**2)*(one-sa(j)**2)*ta(j)&
                        +(one-ra(j)**2)*(one-sa(j)**2)*(ta(j)-one))/two
               shl(3,22,j) = ( (one-ra(j)**2)*(one-sa(j)**2)*ta(j)&
                        + (one-ra(j)**2)*(one-sa(j)**2)*(ta(j)+one))/two
               shl(3,23,j) =  (one-ra(j)**2)*(sa(j)-one)*(-two*ta(j))&
                        *sa(j)/two
               shl(3,24,j) =  (one-ra(j)**2)*(sa(j)+one)*(-two*ta(j))&
                        *sa(j)/two
               shl(3,25,j) = (ra(j)-one)* (one-sa(j)**2)*(-two*ta(j))&
                        *ra(j)/two
               shl(3,26,j) =  (ra(j)+one)* (one-sa(j)**2)*(-two*ta(j))&
                        *ra(j)/two
               shl(3,27,j) =(one - ra(j)**2)*(one - sa(j)**2)*(-two*ta(j))

               shl(2,1,j) = ((ra(j)-one)*(ta(j)-one)*ra(j)*sa(j)*ta(j)&
                        +(ra(j)-one)*(sa(j)-one)*(ta(j)-one)*ra(j)*ta(j))/eight
               shl(2,2,j) = ((ra(j)+one)*(ta(j)-one)*ra(j)*sa(j)*ta(j)&
                        + (ra(j)+one)*(sa(j)-one)*(ta(j)-one)*ra(j)*ta(j))/eight
               shl(2,3,j) =( (ra(j)+one)*(ta(j)-one)*ra(j)*ta(j)*sa(j)&
                        +(ra(j)+one)*(sa(j)+one)*(ta(j)-one)*ra(j)*ta(j))/eight
               shl(2,4,j) = ((ra(j)-one)*(ta(j)-one)*ra(j)*sa(j)*ta(j)&
                        + (ra(j)-one)*(sa(j)+one)*(ta(j)-one)*ra(j)*ta(j))/eight
               shl(2,5,j) = ((ra(j)-one)*(ta(j)+one)*ra(j)*sa(j)*ta(j)&
                        + (ra(j)-one)*(sa(j)-one)*(ta(j)+one)*ra(j)*ta(j))/eight
               shl(2,6,j) = ((ra(j)+one)*(ta(j)+one)*ra(j)*sa(j)*ta(j)&
                        +(ra(j)+one)*(sa(j)-one)*(ta(j)+one)*ra(j)*ta(j))/eight
               shl(2,7,j) = ((ra(j)+one)*(ta(j)+one)*ra(j)*sa(j)*ta(j)&
                        +(ra(j)+one)*(sa(j)+one)*(ta(j)+one)*ra(j)*ta(j))/eight
               shl(2,8,j) = ((ra(j)-one)*(ta(j)+one)*ra(j)*sa(j)*ta(j)&
                        + (ra(j)-one)*(sa(j)+one)*(ta(j)+one)*ra(j)*ta(j))/eight
               shl(2,9,j) = ((one - ra(j)**2)*(ta(j)-one)*sa(j)*ta(j)&
                        + (one - ra(j)**2)*(sa(j)-one)*(ta(j)-one)*ta(j))/four
               shl(2,10,j) = (ra(j)+one)*(- two*sa(j))*(ta(j)-one)&
                        *ra(j)*ta(j)/four
               shl(2,11,j) = ((one - ra(j)**2)*(ta(j)-one)*sa(j)*ta(j)&
                        + (one - ra(j)**2)*(sa(j)+one)*(ta(j)-one)*ta(j))/four
               shl(2,12,j) = (ra(j)-one)*( -two*sa(j))*(ta(j)-one)&
                        *ra(j)*ta(j)/four
               shl(2,13,j) =( (one - ra(j)**2)*(ta(j)+one)*sa(j)*ta(j)&
                        + (one - ra(j)**2)*(sa(j)-one)*(ta(j)+one)*ta(j))/four
               shl(2,14,j) =  (ra(j)+one)*( -two*sa(j))*(ta(j)+one)&
                        *ra(j)*ta(j)/four
               shl(2,15,j) =((one - ra(j)**2)*(ta(j)+one)*sa(j)*ta(j)&
                        +(one - ra(j)**2)*(sa(j)+one)*(ta(j)+one)*ta(j))/four
               shl(2,16,j) = (ra(j)-one)*( - two*sa(j))*(ta(j)+one)&
                        *ra(j)*ta(j)/four
               shl(2,17,j) =  ((ra(j)-one)*(one-ta(j)**2)*ra(j)*sa(j)&
                        +(ra(j)-one)*(sa(j)-one)*(one-ta(j)**2)*ra(j))/four
               shl(2,18,j) =  ((ra(j)+one)*(one-ta(j)**2)*ra(j)*sa(j)&
                        + (ra(j)+one)*(sa(j)-one)*(one-ta(j)**2)*ra(j))/four
               shl(2,19,j) =  ((ra(j)+one)*(one-ta(j)**2)*ra(j)*sa(j)&
                        + (ra(j)+one)*(sa(j)+one)*(one-ta(j)**2)*ra(j))/four
               shl(2,20,j) =((ra(j)-one)*(one-ta(j)**2)*ra(j)*sa(j)&
                        +(ra(j)-one)*(sa(j)+one)*(one-ta(j)**2)*ra(j))/four
               shl(2,21,j) = (one-ra(j)**2)*(-two*sa(j))*(ta(j)-one)&
                        *ta(j)/two
               shl(2,22,j) =  (one-ra(j)**2)*(-two*sa(j))*(ta(j)+one)&
                        *ta(j)/two
               shl(2,23,j) =  ((one-ra(j)**2)*(one-ta(j)**2)*sa(j)&
                        + (one-ra(j)**2)*(sa(j)-one)*(one-ta(j)**2))/two
               shl(2,24,j) = ( (one-ra(j)**2)*(one-ta(j)**2)*sa(j)&
                        +  (one-ra(j)**2)*(sa(j)+one)*(one-ta(j)**2))/two
               shl(2,25,j) = (ra(j)-one)*(-two*sa(j))*(one-ta(j)**2)&
                        *ra(j)/two
               shl(2,26,j) =  (ra(j)+one)*(-two*sa(j))*(one-ta(j)**2)&
                        *ra(j)/two
               shl(2,27,j) =(one - ra(j)**2)*( - two*sa(j))&
                        *(one - ta(j)**2)

               shl(1,1,j) = ( (sa(j)-one)*(ta(j)-one)*ra(j)*sa(j)*ta(j)&
                        + (ra(j)-one)*(sa(j)-one)*(ta(j)-one)*sa(j)*ta(j))/eight
               shl(1,2,j) = ( (sa(j)-one)*(ta(j)-one)*ra(j)*sa(j)*ta(j)&
                        + (ra(j)+one)*(sa(j)-one)*(ta(j)-one)*sa(j)*ta(j))/eight
               shl(1,3,j) = ( (sa(j)+one)*(ta(j)-one)*ra(j)*sa(j)*ta(j)&
                        + (ra(j)+one)*(sa(j)+one)*(ta(j)-one)*sa(j)*ta(j))/eight
               shl(1,4,j) = ( (sa(j)+one)*(ta(j)-one)*ra(j)*sa(j)*ta(j)&
                        + (ra(j)-one)*(sa(j)+one)*(ta(j)-one)*sa(j)*ta(j))/eight
               shl(1,5,j) = ( (sa(j)-one)*(ta(j)+one)*ra(j)*sa(j)*ta(j)&
                        + (ra(j)-one)*(sa(j)-one)*(ta(j)+one)*sa(j)*ta(j))/eight
               shl(1,6,j) = ( (sa(j)-one)*(ta(j)+one)*ra(j)*sa(j)*ta(j)&
                        + (ra(j)+one)*(sa(j)-one)*(ta(j)+one)*sa(j)*ta(j))/eight
               shl(1,7,j) = ( (sa(j)+one)*(ta(j)+one)*ra(j)*sa(j)*ta(j)&
                        + (ra(j)+one)*(sa(j)+one)*(ta(j)+one)*sa(j)*ta(j))/eight
               shl(1,8,j) = ( (sa(j)+one)*(ta(j)+one)*ra(j)*sa(j)*ta(j)&
                        + (ra(j)-one)*(sa(j)+one)*(ta(j)+one)*sa(j)*ta(j))/eight
               shl(1,9,j) = (- two*ra(j))*(sa(j)-one)*(ta(j)-one)&
                        *sa(j)*ta(j)/four
               shl(1,10,j) = ( (one - sa(j)**2)*(ta(j)-one)*ra(j)*ta(j)&
                        + (ra(j)+one)*(one - sa(j)**2)*(ta(j)-one)*ta(j))/four
               shl(1,11,j) = (-two*ra(j))*(sa(j)+one)*(ta(j)-one)&
                        *sa(j)*ta(j)/four
               shl(1,12,j) = ( (one - sa(j)**2)*(ta(j)-one)*ra(j)*ta(j)&
                        + (ra(j)-one)*(one - sa(j)**2)*(ta(j)-one)*ta(j))/four
               shl(1,13,j) = ( - two*ra(j))*(sa(j)-one)*(ta(j)+one)&
                        *sa(j)*ta(j)/four
               shl(1,14,j) = ( (one - sa(j)**2)*(ta(j)+one)*ra(j)*ta(j)&
                        + (ra(j)+one)*(one - sa(j)**2)*(ta(j)+one)*ta(j))/four
               shl(1,15,j) =(- two*ra(j))*(sa(j)+one)*(ta(j)+one)&
                        *sa(j)*ta(j)/four
               shl(1,16,j) = ( (one - sa(j)**2)*(ta(j)+one)*ra(j)*ta(j)&
                        + (ra(j)-one)*(one - sa(j)**2)*(ta(j)+one)*ta(j))/four
               shl(1,17,j) =  ((sa(j)-one)*(one-ta(j)**2)*ra(j)*sa(j)&
                        + (ra(j)-one)*(sa(j)-one)*(one-ta(j)**2)*sa(j))/four
               shl(1,18,j) =  ((sa(j)-one)*(one-ta(j)**2)*ra(j)*sa(j)&
                        + (ra(j)+one)*(sa(j)-one)*(one-ta(j)**2)*sa(j))/four
               shl(1,19,j) =  ((sa(j)+one)*(one-ta(j)**2)*ra(j)*sa(j)&
                        + (ra(j)+one)*(sa(j)+one)*(one-ta(j)**2)*sa(j))/four
               shl(1,20,j) =  ((sa(j)+one)*(one-ta(j)**2)*ra(j)*sa(j)&
                        + (ra(j)-one)*(sa(j)+one)*(one-ta(j)**2)*sa(j))/four
               shl(1,21,j) = (-two*ra(j))*(one-sa(j)**2)*(ta(j)-one)&
                        *ta(j)/two
               shl(1,22,j) =  (-two*ra(j))*(one-sa(j)**2)*(ta(j)+one)&
                        *ta(j)/two
               shl(1,23,j) =  (-two*ra(j))*(sa(j)-one)*(one-ta(j)**2)&
                        *sa(j)/two
               shl(1,24,j) =  (-two*ra(j))*(sa(j)+one)*(one-ta(j)**2)&
                        *sa(j)/two
               shl(1,25,j) = ( (one-sa(j)**2)*(one-ta(j)**2)*ra(j)&
                        +(ra(j)-one)* (one-sa(j)**2)*(one-ta(j)**2))/two
               shl(1,26,j) =  ( (one-sa(j)**2)*(one-ta(j)**2)*ra(j)&
                        + (ra(j)+one)* (one-sa(j)**2)*(one-ta(j)**2))/two
               shl(1,27,j) =(-two*ra(j))*(one - sa(j)**2)*(one - ta(j)**2)
            30    continue
         endif
      endif

      return
   end subroutine

   !======================================================================
   function rowdot(a,b,ma,mb,n)
      ! .... program to compute the dot product of vectors stored row-wise

      implicit none

      !.... remove above card for single precision operation

      integer*4:: ma, mb, n
      real*8  :: a(ma,*),b(mb,*)

      real*8  :: rowdot
      integer*4:: i

      rowdot = 0.0d0

      do i=1,n
         rowdot = rowdot + a(1,i)*b(1,i)
      enddo
      
     return
   end function

   !======================================================================
   subroutine shlten(shlen, nen)
      ! Calcula funcoes de interpolacao e suas derivadas locais para
      ! elementos triangulares nos nos do elemento

      ! - ----------------------------------------------- - -------------------- -
        
      !            s,t = local element coordinates ("xi", "eta", resp.)
      !   shlen(1,i,l) = local ("xi") derivative of shape function i in node l
      !   shlen(2,i,l) = local ("eta") derivative of shape function i in node l
      !   shlen(3,i,l) = local shape function i in node l
      !            nen = number of element nodes
      !
      implicit none

      integer :: l, nen
      real*8 :: r, s
      real*8 :: shlen(3, nen, nen), ra(nen), sa(nen)

      if (nen == 3) then
         ra(1) = one
         sa(1) = zero
         ra(2) = zero
         sa(2) = one
         ra(3) = zero
         sa(3) = zero
      end if

      do l = 1, nen
         r = ra(l)
         s = sa(l)
         shlen(1, 1, l) = one
         shlen(2, 1, l) = zero
         shlen(3, 1, l) = r
         shlen(1, 2, l) = zero
         shlen(2, 2, l) = one
         shlen(3, 2, l) = s
         shlen(1, 3, l) = -one
         shlen(2, 3, l) = -one
         shlen(3, 3, l) = 1 - r - s
      end do
      return
   end subroutine

   !======================================================================
   subroutine shlqen(shlen, nen) !alterei esta rotina... corrigi a numeração
      ! Calcula funcoes de interpolacao e suas derivadas locais para
      ! elementos quadrangulares nos nos do elemento
        
      ! - ----------------------------------------------- - -------------------- -
        
      !            s,t = local element coordinates ("xi", "eta", resp.)
      !   shlen(1,i,l) = local ("xi") derivative of shape function i in node l
      !   shlen(2,i,l) = local ("eta") derivative of shape function i in node l
      !   shlen(3,i,l) = local shape function i in node l
      !            nen = number of element nodes
        
      implicit none

      integer :: l, nen
      real*8 :: onemr, onems, onepr, oneps
      real*8 :: r, s
      real*8 :: shlen(3, nen, nen), ra(nen), sa(nen)

      if (nen == 4) then
         ra(1) = -one
         sa(1) = -one
         ra(2) =  one
         sa(2) = -one
         ra(3) =  one
         sa(3) =  one
         ra(4) = -one
         sa(4) =  one
      end if

      do l = 1, nen
         r = ra(l)
         s = sa(l)
         onepr = one + r
         onemr = one - r
         oneps = one + s
         onems = one - s

         shlen(1, 1, l) = -onems * pt25
         shlen(2, 1, l) = -onemr * pt25
         shlen(3, 1, l) = onemr * onems * pt25
            
         shlen(1, 2, l) = onems * pt25
         shlen(2, 2, l) = -onepr * pt25
         shlen(3, 2, l) = onepr * onems * pt25

         shlen(1, 3, l) = oneps * pt25
         shlen(2, 3, l) = onepr * pt25
         shlen(3, 3, l) = onepr * oneps * pt25
            
         shlen(1, 4, l) = -oneps * pt25
         shlen(2, 4, l) = onemr * pt25
         shlen(3, 4, l) = onemr * oneps * pt25
      end do
      return
   end subroutine

   !======================================================================
   subroutine shlqen3D(shl, nen)
      ! Calcula funcoes de interpolacao e suas derivadas locais para
      ! elementos quadrangulares nos nos do elemento
        
      ! - ----------------------------------------------- - -------------------- -
        
      !            s,t = local element coordinates ("xi", "eta", resp.)
      !   shlen(1,i,l) = local ("xi") derivative of shape function i in node l
      !   shlen(2,i,l) = local ("eta") derivative of shape function i in node l
      !   shlen(3,i,l) = local shape function i in node l
      !            nen = number of element nodes
        
      implicit none

      integer :: j, nen
      real*8 :: onemr, onems, onepr, oneps, onept, onemt
      real*8 :: r, s, t, eight
      real*8 :: shl(4, nen, nen), ra(nen), sa(nen), ta(nen)
        
      eight=8.d0
        
      if (nen == 8) then
         ra(1) =-one
         sa(1) =-one
         ta(1) =-one

         ra(2) = one
         sa(2) =-one
         ta(2) =-one

         ra(4) =-one
         sa(4) = one
         ta(4) =-one

         ra(3) = one
         sa(3) = one
         ta(3) =-one

         ra(5) =-one
         sa(5) =-one
         ta(5) = one

         ra(6) = one
         sa(6) =-one
         ta(6) = one
   
         ra(8) =-one
         sa(8) = one
         ta(8) = one

         ra(7) = one
         sa(7) = one
         ta(7) = one

      end if
        
      do j=1,nen
         shl(1,1,j) = -(one - sa(j))*(one - ta(j))/eight
         shl(2,1,j) = -(one - ra(j))*(one - ta(j))/eight
         shl(3,1,j) = -(one - ra(j))*(one - sa(j))/eight
         shl(4,1,j) =  (one - ra(j))*(one - sa(j))*(one - ta(j))/eight

         shl(1,2,j) =  (one - sa(j))*(one - ta(j))/eight           
         shl(2,2,j) = -(one + ra(j))*(one - ta(j))/eight           
         shl(3,2,j) = -(one + ra(j))*(one - sa(j))/eight
         shl(4,2,j) =  (one + ra(j))*(one - sa(j))*(one - ta(j))/eight

         shl(1,3,j) =  (one + sa(j))*(one - ta(j))/eight           
         shl(2,3,j) =  (one + ra(j))*(one - ta(j))/eight
         shl(3,3,j) = -(one + ra(j))*(one + sa(j))/eight
         shl(4,3,j) =  (one + ra(j))*(one + sa(j))*(one - ta(j))/eight           

         shl(1,4,j) = -(one + sa(j))*(one - ta(j))/eight           
         shl(2,4,j) =  (one - ra(j))*(one - ta(j))/eight
         shl(3,4,j) = -(one - ra(j))*(one + sa(j))/eight
         shl(4,4,j) =  (one - ra(j))*(one + sa(j))*(one - ta(j))/eight

         shl(1,5,j) = -(one - sa(j))*(one + ta(j))/eight           
         shl(2,5,j) = -(one - ra(j))*(one + ta(j))/eight
         shl(3,5,j) =  (one - ra(j))*(one - sa(j))/eight
         shl(4,5,j) =  (one - ra(j))*(one - sa(j))*(one + ta(j))/eight
           
         shl(1,6,j) =  (one - sa(j))*(one + ta(j))/eight           
         shl(2,6,j) = -(one + ra(j))*(one + ta(j))/eight
         shl(3,6,j) =  (one + ra(j))*(one - sa(j))/eight           
         shl(4,6,j) =  (one + ra(j))*(one - sa(j))*(one + ta(j))/eight
           
         shl(1,7,j) =  (one + sa(j))*(one + ta(j))/eight
         shl(2,7,j) =  (one + ra(j))*(one + ta(j))/eight           
         shl(3,7,j) =  (one + ra(j))*(one + sa(j))/eight           
         shl(4,7,j) =  (one + ra(j))*(one + sa(j))*(one + ta(j))/eight
           
         shl(1,8,j) = -(one + sa(j))*(one + ta(j))/eight           
         shl(2,8,j) =  (one - ra(j))*(one + ta(j))/eight
         shl(3,8,j) =  (one - ra(j))*(one + sa(j))/eight
         shl(4,8,j) =  (one - ra(j))*(one + sa(j))*(one + ta(j))/eight
      end do

      ! if (nen == 8) then
      !    ra(1) =-one
      !    sa(1) =-one
      !    ta(1) =-one
      
      !    ra(2) = one
      !    sa(2) =-one
      !    ta(2) =-one

      !    ra(3) =-one
      !    sa(3) = one
      !    ta(3) =-one

      !    ra(4) = one
      !    sa(4) = one
      !    ta(4) =-one

      !    ra(5) =-one
      !    sa(5) =-one
      !    ta(5) = one

      !    ra(6) = one
      !    sa(6) =-one
      !    ta(6) = one
   
      !    ra(7) =-one
      !    sa(7) = one
      !    ta(7) = one

      !    ra(8) = one
      !    sa(8) = one
      !    ta(8) = one
      ! end if

      ! do j=1,nen
      !    shl(4,1,j) = (one - ra(j))*(one - sa(j))*(one - ta(j))/eight
      !    shl(4,2,j) = (one + ra(j))*(one - sa(j))*(one - ta(j))/eight
      !    shl(4,3,j) = (one + ra(j))*(one + sa(j))*(one - ta(j))/eight
      !    shl(4,4,j) = (one - ra(j))*(one + sa(j))*(one - ta(j))/eight
      !    shl(4,5,j) = (one - ra(j))*(one - sa(j))*(one + ta(j))/eight
      !    shl(4,6,j) = (one + ra(j))*(one - sa(j))*(one + ta(j))/eight
      !    shl(4,7,j) = (one + ra(j))*(one + sa(j))*(one + ta(j))/eight
      !    shl(4,8,j) = (one - ra(j))*(one + sa(j))*(one + ta(j))/eight

      !    shl(3,1,j) = -(one - ra(j))*(one - sa(j))/eight
      !    shl(3,2,j) = -(one + ra(j))*(one - sa(j))/eight
      !    shl(3,3,j) = -(one + ra(j))*(one + sa(j))/eight
      !    shl(3,4,j) = -(one - ra(j))*(one + sa(j))/eight
      !    shl(3,5,j) =  (one - ra(j))*(one - sa(j))/eight
      !    shl(3,6,j) =  (one + ra(j))*(one - sa(j))/eight
      !    shl(3,7,j) =  (one + ra(j))*(one + sa(j))/eight
      !    shl(3,8,j) =  (one - ra(j))*(one + sa(j))/eight

      !    shl(2,1,j) = -(one - ra(j))*(one - ta(j))/eight
      !    shl(2,2,j) = -(one + ra(j))*(one - ta(j))/eight
      !    shl(2,3,j) =  (one + ra(j))*(one - ta(j))/eight
      !    shl(2,4,j) =  (one - ra(j))*(one - ta(j))/eight
      !    shl(2,5,j) = -(one - ra(j))*(one + ta(j))/eight
      !    shl(2,6,j) = -(one + ra(j))*(one + ta(j))/eight
      !    shl(2,7,j) =  (one + ra(j))*(one + ta(j))/eight
      !    shl(2,8,j) =  (one - ra(j))*(one + ta(j))/eight

      !    shl(1,1,j) = -(one - sa(j))*(one - ta(j))/eight
      !    shl(1,2,j) =  (one - sa(j))*(one - ta(j))/eight
      !    shl(1,3,j) =  (one + sa(j))*(one - ta(j))/eight
      !    shl(1,4,j) = -(one + sa(j))*(one - ta(j))/eight
      !    shl(1,5,j) = -(one - sa(j))*(one + ta(j))/eight
      !    shl(1,6,j) =  (one - sa(j))*(one + ta(j))/eight
      !    shl(1,7,j) =  (one + sa(j))*(one + ta(j))/eight
      !    shl(1,8,j) = -(one + sa(j))*(one + ta(j))/eight
      ! end do
      
      return
   end subroutine
   !======================================================================
   subroutine gerarPtosIntegracaoLaterais2D
      use mGlobaisEscalares, only: npint, npintFratura
      use mMalha, only: nsd
      select case(npint)
      case(1)
         npintFratura = 1
      case(3)
         npintFratura = 2
      case(4)
         if(nsd==2) npintFratura = 2
         if(nsd==3) npintFratura = 3
      case(7)
         npintFratura = 3
      case(8)
         npintFratura = 4
      case(9)
         npintFratura = 3
      case(16)
         npintFratura = 4
      case(25)
         npintFratura = 5
      case(36)
         npintFratura = 6
      case(49)
         npintFratura = 7
      case(64)
         npintFratura = 8
      end select
   end subroutine gerarPtosIntegracaoLaterais2D

   !======================================================================
   subroutine shl1D_elemNodes(shlen, nen)
      ! Calcula funcoes de interpolacao e suas derivadas locais avaliadas nos nos.
      ! Implementacao para elementos unidimensionais
      
      ! - ----------------------------------------------- - -------------------- -
      
      !            r = local element coordinate ("xi")
      !   shlen(1,i,l) = local ("xi") derivative of shape function i in node l
      !   shlen(2,i,l) = local shape function i in node l
      !            nen = number of element nodes
      !            nsd = number of the space dimension

      real*8, intent(out) :: shlen(2, nen, nen)
      integer, intent(in)     :: nen

      real*8 :: r(nen)
      integer :: l

      call masterNodeCoords1D(r, nen)

      do l = 1, nen
         call shl1D_point(r(l), nen, shlen(:,:,l))
      end do
   end subroutine shl1D_elemNodes

   !======================================================================
   subroutine shl1D_intPoints(shl, w, nint, nen)
      ! Calculate integration-rule weights, shape functions
      ! and local derivatives for an one-dimensional element
      
      !              r = local element coordinate ("xi")
      !     shl(1,i,l) = local ("xi") derivative of shape function
      !     shl(2,i,l) = local  shape function
      !           w(l) = integration-rule weight
      !              i = local node number
      !              l = integration point number
      !           nint = number of integration points
      
      integer, intent(in) :: nint, nen
      real*8, intent(out) :: shl(2, nen, nint), w(nint)
      integer :: l
      real*8 :: r(nint)

      !Atribuicao dos pesos e posicoes da quadratura
      call masterIntCoords1D(r, w, nint)

      !Atribuicao dos valores das shape functions
      do l = 1, nint
         call shl1D_point(r(l), nen, shl(:,:,l))
      end do
   end subroutine shl1D_intPoints

   !======================================================================
   subroutine shl1D_point(r, nen, shl)
      ! Shape functions and its first derivative evaluated at a given point r
      real *8, intent(in)  :: r
      integer, intent(in)      :: nen
      real *8, intent(out) :: shl(2,nen)
      integer                 :: i, j, k
      real*8              :: aa, aax, bb, daj, rElemNodes(nen)

      call masterNodeCoords1D(rElemNodes, nen)

      if (nen == 1) then
         shl(1,1) = zero
         shl(2,1) = one
         return
      end if

      do i = 1, nen
         aa = one
         bb = one
         aax = zero
         do j = 1, nen
            daj = one
            if (i /= j)then
               aa = aa * ( r - rElemNodes(j))
               bb = bb * ( rElemNodes(i) - rElemNodes(j))
               do k = 1, nen
                  if (k /= i .and. k /= j) daj = daj * ( r - rElemNodes(k))
               end do
               aax = aax + daj
            endif
         end do
         shl(2,i) = aa/bb
         shl(1,i) = aax/bb
      end do
   end subroutine shl1D_point

   !======================================================================
   subroutine shg1D_point(xl, det, shl, shg, nen, nsd)
      ! program to calculate global derivatives of shape functions
      ! and jacobian determinants for the one-dimensional element

      !     xl(i,j) = global coordinates of nodes
      !  shl(1,j,l) = local ("xi") derivative of shape function j on integration point l
      !  shl(2,j,l) = value of shape function j on integration point l
      !  shg(1,j,l) = global ("arc-length") derivative of shape function j on integration point l
      !  shg(2,j,l) = shl(2,i,l)
      !      det(l) = euclidean length
      !           i = global coordinate number
      !           j = local node number
      !           l = integration-point number
      !        nint = number of integration points

      integer, intent(in)     :: nen, nsd
      real*8, intent(in)  :: xl(nsd,nen), shl(2,nen)
      real*8, intent(out) :: shg(2,nen), det
      integer :: sinal
      real*8 :: projecao
      real*8  :: diff(nsd), length, sinTheta, cosTheta

      integer                 :: i

      if(nen == 1) return

      det = det1D(xl, shl, nen, nsd)

      if (det <= zero) then
         write(*, "(a)") "Error: shg1D_point found a non-positive determinant"
         stop
      end if

      diff(1) = xl(1,2) - xl(1,1)
      diff(2) = xl(2,2) - xl(2,1)
      length = sqrt(diff(1)*diff(1) + diff(2)*diff(2))
      sinTheta = abs(diff(2))/length
      if (diff(2) > zero) then
         cosTheta = diff(1)/length
      else if (diff(2) < zero) then
         cosTheta = -diff(1)/length
      else
         cosTheta = -abs(diff(1))/length
      end if
        
      projecao=diff(1)*cosTheta+diff(2)*sinTheta
      if(projecao>0)  then
         sinal=1
      else
         sinal=-1
      endif
      
      do i=1,nen
         shg(1,i)=sinal*shl(1,i)/det
         shg(2,i)=shl(2,i)
      end do
   end subroutine shg1D_point

   !======================================================================
   subroutine shg1D_points(xl, det, shl, shg, nen, npoints, nsd)
      ! program to calculate global derivatives of shape functions and
      ! jacobian determinants at the integration points
      
      !     xl(j,i)    = global coordinates
      !     det(l)     = jacobian determinant
      !     shl(1,i,l) = local ("xi") derivative of shape function
      !     shl(2,i,l) = local ("eta") derivative of shape function
      !     shl(3,i,l) = local  shape function
      !     shg(1,i,l) = x-derivative of shape function
      !     shg(2,i,l) = y-derivative of shape function
      !     shg(3,i,l) = shl(3,i,l)
      !              l = integration-point number
      !           nint = number of integration points
      
      integer :: l
      integer, intent(in) :: npoints, nen, nsd
      real*8, intent(in) :: xl(nsd, nen)
      real*8, intent(in) :: shl(2, nen, npoints)
      real*8, intent(out) :: shg(2, nen, npoints)
      real*8, intent(out) :: det(npoints)

      do l = 1, npoints
         call shg1D_point(xl, det(l), shl(:,:,l), shg(:,:,l), nen, nsd)
      end do
   end subroutine shg1D_points
   
   !======================================================================
   subroutine shgCodimOneSurface_point(xl, det, shl, shg, nen, nsd)
      ! program to calculate global derivatives of shape functions
      ! and jacobian determinants for a codimensional one surface element at a given point:
      ! surface element if nsd==3 or curve element of nsd == 2
      
      !      xl(i,j) = global coordinates of nodes
      !     shl(1,j) = local ("xi") derivative of shape function j at a given point (nsd==2)
      !     shl(2,j) = local ("eta") derivative of shape function j at a given point (nsd==3)
      !   shl(nsd,j) = value of shape function j at a given point
      !     shg(1,j) = global x-derivative of shape function j at a given point
      !     shg(2,j) = global y-derivative of shape function j at a given point
      !     shg(3,j) = global z-derivative of shape function j at a given point
      ! shg(nsd+1,j) = shl(nsd,i)
      !          det = square root of the determinant of the g matrix
      
      !       d(i,j) = local j-derivative of the i-coordinate of a given point
      !            g = first fundamental form (d^T d)
      !         detG = determinant of matrix g
      !        g_inv = inverse of matrix g

      integer, intent(in) :: nen, nsd
      real*8, intent(in)  :: xl(nsd,nen), shl(nsd,nen)
      real*8, intent(out) :: shg(nsd+1,nen), det

      real*8 :: d(nsd,nsd-1), g(nsd-1,nsd-1), g_inv(nsd-1,nsd-1), detG
      integer :: i,j,no

      if(nen == 1) return

      do i=1,nsd
         do j=1,nsd-1
            d(i,j) = dot_product(xl(i,:),shl(j,:))
         end do
      end do

      do i=1,nsd-1
         do j=1,nsd-1
            g(i,j) = dot_product(d(:,i),d(:,j))
         end do
      end do

      if (nsd==2) then
         detG = g(1,1)
         g_inv = 1d0/detG
      else if (nsd==3) then
         detG = g(1,1)*g(2,2) - g(1,2)*g(2,1)
         g_inv(1,1) = g(2,2)
         g_inv(2,2) = g(1,1)
         g_inv(1,2) = -g(1,2)
         g_inv(2,1) = -g(2,1)
         g_inv(:,:) = g_inv(:,:)/detG
      else
         print*, "Erro na subrotina shgSurface_point: nsd deve ser 2 ou 3"
         stop
      end if

      if (detG <= 0d0) then
         print*, "Erro na subrotina shgSurface_point: determinante nulo ou negativo"
         stop
      end if
      det = sqrt(detG)

      do no=1,nen
         do i=1,nsd
            shg(i,no) = 0d0
            do j=1,nsd-1
               shg(i,no) = shg(i,no) + d(i,j)*dot_product(g_inv(j,:),shl(1:nsd-1,no))
            end do
         end do
         shg(nsd+1,no) = shl(nsd,no)
      end do
   end subroutine shgCodimOneSurface_point

   !======================================================================
   subroutine shgCodimOneSurface_points(xl, det, shl, shg, nen, npoints, nsd)
      ! program to calculate global derivatives of shape functions
      ! and jacobian determinants for a codimensional one surface element at the set of points:
      ! surface element if nsd==3 or line element of nsd == 2
      
      !        xl(i,j) = global coordinates of nodes
      !     shl(1,j,l) = local ("xi") derivative of shape function j at a point l (nsd==2)
      !     shl(2,j,l) = local ("eta") derivative of shape function j at a point l (nsd==3)
      !   shl(nsd,j,l) = value of shape function j at a point l
      !     shg(1,j,l) = global x-derivative of shape function j at a point l
      !     shg(2,j,l) = global y-derivative of shape function j at a point l
      !     shg(3,j,l) = global z-derivative of shape function j at a point l
      ! shg(nsd+1,j,l) = shl(nsd,i,l)
      !         det(l) = square root of the determinant of the g matrix at a point l
      !        npoints = number of points
      
      !              l = point number

      integer :: l
      integer, intent(in) :: npoints, nen, nsd
      real*8, intent(in) :: xl(nsd, nen)
      real*8, intent(in) :: shl(nsd, nen, npoints)
      real*8, intent(out) :: shg(nsd+1, nen, npoints)
      real*8, intent(out) :: det(npoints)

      do l = 1, npoints
         call shgCodimOneSurface_point(xl, det(l), shl(:,:,l), shg(:,:,l), nen, nsd)
      end do
   end subroutine shgCodimOneSurface_points

   !======================================================================
   subroutine masterNodeCoords1D(r, nen)
      ! Set the master element coordinates of the nodes
      integer, intent(in)     :: nen
      real*8, intent(out) :: r(nen)

      select case(nen)
      case (1)
         r(1) = zero
      case (2)
         r(1) = -one
         r(2) = one
      case (3)
         r(1)= -one
         r(2)= one
         r(3)= zero
      case (4)
         r(1) = -one
         r(2) = one
         r(3) = -.333333333333333d0
         r(4) =  .333333333333333d0
      case (5)
         r(1)= -one
         r(2)=  one
         r(3)= -pt5
         r(4)= zero
         r(5)= pt5
      case (6)
         r(1) = -one
         r(2) =  one
         r(3) = -.6d0
         r(4) = -.2d0
         r(5) =  .2d0
         r(6) =  .6d0
      case (7)
         r(1) = -one
         r(2) =  one
         r(3) = -.666666666666666d0
         r(4) = -.333333333333333d0
         r(5) = zero
         r(6) =  .333333333333333d0
         r(7) =  .666666666666666d0
      case (8)
         r(1) = -one
         r(2) =  one
         r(3) = -0.71428571428571d0
         r(4) = -0.42857142857143d0
         r(5) = -0.14285714285714d0
         r(6) = 0.14285714285714d0
         r(7) = 0.42857142857143d0
         r(8) = 0.71428571428571d0
      case default
         write(*, "(a, i2, a)") 'Error: subrotina masterNodeCoords1D não está implementada para elementos 1D com ', &
                  nen, ' nós!'
         stop
      end select

   end subroutine masterNodeCoords1D

   !======================================================================
   subroutine masterIntCoords1D(r, w, nint)
      ! Set the master element coordinates of the integration points and its respective weights
      integer, intent(in)     :: nint
      real*8, intent(out) :: r(nint), w(nint)
      real*8, parameter   :: five9=0.5555555555555555d0
      real*8, parameter   :: eight9=0.8888888888888888d0

      !Atribuicao dos pesos e posicoes da quadratura
      select case(nint)
      case(1)
         w(1) = two
         r(1) = zero
      case (2)
         w(1) = one
         w(2) = one
         r(1)=-.577350269189626d0
         r(2)= .577350269189625d0
      case (3)
         w(1) = five9
         w(2) = five9
         w(3) = eight9
         r(1)=-.774596669241483d0
         r(2)= .774596669241483d0
         r(3)= zero
      case (4)
         w(1) = .347854845137454d0
         w(2) = .347854845137454d0
         w(3) = .652145154862546d0
         w(4) = .652145154862546d0
         r(1)=-.861136311594053d0
         r(2)= .861136311594053d0
         r(3)=-.339981043584856d0
         r(4)= .339981043584856d0
      case (5)
         w(1) = .236926885056189d0
         w(2) = .236926885056189d0
         w(3) = .478628670499366d0
         w(4) = .478628670499366d0
         w(5) = .568888888888888d0
         r(1)=-.906179845938664d0
         r(2)= .906179845938664d0
         r(3)=-.538469310105683d0
         r(4)= .538469310105683d0
         r(5)= zero
      case (6)
         w(1) = .171324492397170d0
         w(2) = .171324492397170d0
         w(3) = .360761573048139d0
         w(4) = .360761573048139d0
         w(5) = .467913934572691d0
         w(6) = .467913934572691d0
         r(1)=-.932469514203152d0
         r(2)= .932469514203152d0
         r(3)=-.661209386466265d0
         r(4)= .661209386466365d0
         r(5)=-.238619186083197d0
         r(6)= .238619186083197d0
      case (7)
         w(1) = .129484966168870d0
         w(2) = .129484966168870d0
         w(3) = .279705391489277d0
         w(4) = .279705391489277d0
         w(5) = .381830050505119d0
         w(6) = .381830050505119d0
         w(7) = .417959183673469d0
         r(1)=-.949107912342759d0
         r(2)= .949107912342759d0
         r(3)=-.741531185599394d0
         r(4)= .741531185599394d0
         r(5)=-.405845151377397d0
         r(6)= .405845151377397d0
         r(7)= zero
      case (8)
         w(1) = .101228536290376d0
         w(2) = .101228536290376d0
         w(3) = .222381034453374d0
         w(4) = .222381034453374d0
         w(5) = .313706645877887d0
         w(6) = .313706645877887d0
         w(7) = .362683783378362d0
         w(8) = .362683783378362d0
         r(1)=-.960289856497536d0
         r(2)= .960289856497536d0
         r(3)=-.796666477413627d0
         r(4)= .796666477413627d0
         r(5)=-.525532409916329d0
         r(6)= .525532409916329d0
         r(7)=-.183434642495650d0
         r(8)= .183434642495650d0
      case default
         write(*, "(a, i2, a)") 'Error: subrotina masterIntCoords1D não está implementada para elementos 1D com ', &
                  nint, ' pontos de integração!'
         stop
      end select
   end subroutine masterIntCoords1D

   !======================================================================
   function det1D(xl, shl, nen, nsd)
      real*8 :: det1D
      integer, intent(in)     :: nen, nsd
      real*8, intent(in)  :: xl(nsd,nen), shl(2,nen)

      integer    :: i
      real*8 :: x(nsd)

      x(:) = zero
      do i=1,nen
         x(:) = x(:) + shl(1,i)*xl(:,i)
      end do
      det1D = dsqrt(dot_product(x, x))
   end function det1D

   !======================================================================
   subroutine shltPoint3D(r,s,t,nen,shl)
      ! Shape functions and its first derivatives evaluated at a given point of coords (r,s)

      use mGlobaisEscalares, only: one, zero

      real*8, intent(in)  :: r,s,t
      integer, intent(in) :: nen
      real*8, intent(out) :: shl(4,nen)

      real*8 :: u

      u = 1d0-r-s-t

      shl(1,1) = one
      shl(2,1) = zero
      shl(3,1) = zero
      shl(4,1) = r

      shl(1,2) = zero
      shl(2,2) = zero
      shl(3,2) = one
      shl(4,2) = t

      shl(1,3) = zero
      shl(2,3) = one
      shl(3,3) = zero
      shl(4,3) = s

      shl(1,4) = -one
      shl(2,4) = -one
      shl(3,4) = -one
      shl(4,4) = u
   end subroutine shltPoint3D

   !======================================================================
   subroutine shlt3D(shl,w,npint,nen)
      ! .... program to calculate integration-rule weights, shape functions
      ! and local derivatives for a triangular element
              
      !        r, s, t = local element coordinates
      !        shl(j,i,l) = local ("j") derivative of shape function
      !        shl(3,i,l) = local  shape function
      !              w(l) = integration-rule weight
      !              npint = number of integration points, eq. 1 or 4

      ! use mGlobaisEscalares

      integer*4, intent(in) :: npint, nen
      real*8, intent(out) :: shl(4,nen,npint),w(npint)

      integer :: i
      real*8 :: r(npint), s(npint), t(npint)

      call masterIntCoords3D_tetra(r, s, t, w, npint)
      do i=1,npint
         call shltPoint3D(r(i), s(i), t(i), nen, shl(:,:,i))
      end do

   end subroutine shlt3D

   !======================================================================
   subroutine shlten3D(shl,nen)
      ! .... program to calculate integration-rule weights, shape functions
      ! and local derivatives for a triangular element
              
      !        c1, c2, c3 = local element coordinates ("l1", "l2", "l3".)
      !        shl(j,i,l) = local ("j") derivative of shape function
      !        shl(3,i,l) = local  shape function
      !                 i = local node number

      ! use mGlobaisEscalares

      integer*4, intent(in) :: nen
      real*8, intent(out) :: shl(4,nen,nen)

      integer :: i
      real*8 :: r(nen), s(nen), t(nen)

      call masterNodeCoords3D_tetra(r, s, t, nen)
      do i=1,nen
         call shltPoint3D(r(i), s(i), t(i), nen, shl(:,:,i))
      end do

   end subroutine shlten3D

   !======================================================================
   subroutine masterIntCoords3D_tetra(r, s, t, w, nint)
      ! Set the master element coordinates of the integration points and its respective weights

      use mGlobaisEscalares, only: one, zero, pt25, pt45, pt5, pt8, pt1667, six

      integer, intent(in)     :: nint
      real*8, intent(out) :: r(nint), s(nint), t(nint), w(nint)

      !Atribuicao dos pesos e posicoes da quadratura
      select case(nint)
      case(1)
         w(1) = one/six

         r(1) = pt25
         s(1) = pt25
         t(1) = pt25

      case (4)
         w(1) = pt25/six
         w(2) = pt25/six
         w(3) = pt25/six
         w(4) = pt25/six

         r(1) = 0.58541020d0
         s(1) = 0.13819660d0
         t(1) = 0.13819660d0

         r(2) = 0.13819660d0
         s(2) = 0.58541020d0
         t(2) = 0.13819660d0

         r(3) = 0.13819660d0
         s(3) = 0.13819660d0
         t(3) = 0.58541020d0

         r(4) = 0.13819660d0
         s(4) = 0.13819660d0
         t(4) = 0.13819660d0

      case (5)
         w(1) = pt45/six
         w(2) = pt45/six
         w(3) = pt45/six
         w(4) = pt45/six
         w(5) = -pt8/six

         r(1) = pt5
         s(1) = pt1667
         t(1) = pt1667

         r(2) = pt1667
         s(2) = pt5
         t(2) = pt1667

         r(3) = pt1667
         s(3) = pt1667
         r(3) = pt5

         r(4) = pt1667
         s(4) = pt1667
         t(4) = pt1667

         r(5) = pt25
         s(5) = pt25
         t(5) = pt25

      case default
         write(*, "(a, i2, a)") 'Error: subrotina masterIntCoords3D_tetra não está implementada para tetraedros com ', &
                  nint, ' pontos de integração!'
         stop
      end select
   end subroutine masterIntCoords3D_tetra

   !======================================================================
   subroutine masterNodeCoords3D_tetra(r, s, t, nen)
      ! Set the master element coordinates of the integration points and its respective weights

      use mGlobaisEscalares, only: one, zero

      integer, intent(in) :: nen
      real*8, intent(out) :: r(nen), s(nen), t(nen)

      !Atribuicao dos pesos e posicoes da quadratura
      select case(nen)
      case (4)
         r(1) = one
         s(1) = zero
         t(1) = zero

         r(2) = zero
         s(2) = one
         t(2) = zero

         r(3) = zero
         s(3) = zero
         t(3) = one

         r(4) = zero
         s(4) = zero
         t(4) = zero

      case default
         write(*, "(a, i2, a)") 'Error: subrotina masterNodeCoords3D_tetra nao esta implementada para tetraedros com ', &
                  nen, ' nos!'
         stop
      end select
   end subroutine masterNodeCoords3D_tetra

   !======================================================================
   subroutine shg3D_Josue(xl,det,shl,shg,npint,nel,nen)
      ! .... program to calculate global derivatives of shape functions and
      ! jacobian determinants for a  quadrilateral element
         
      !        xl(j,i)    = global coordinates
      !        det(l)     = jacobian determinant
      !        shl(1,i,l) = local ("xi") derivative of shape function
      !        shl(2,i,l) = local ("eta") derivative of shape function
      !        shl(3,i,l) = local ("zeta") derivative of shape function
      !        shl(4,i,l) = local  shape function
      !        shg(1,i,l) = x-derivative of shape function
      !        shg(2,i,l) = y-derivative of shape function
      !        shg(3,i,l) = z-derivative of shape funciton
      !        dx_dr(i,j) = jacobian matrix (derivative of x_i with respect to r_j )
      !                 i = global coordinate number
      !                 j = local coordinate number
      !                 l = integration-point number
      !              npint = number of integration points, eq. 1 or 4

      use mLeituraEscrita, only: iecho   

      implicit none
         
      !.... remove above card for single precision operation
      integer*4, intent(in) :: npint, nel, nen
      real*8, intent(in) :: xl(3,nen), shl(4,nen,npint)
      real*8, intent(out) :: det(npint), shg(4,nen,npint)
         
      real*8 :: dx_dr(3,3)
      integer*4 i, j, l
      real*8 cof11, cof12, cof13, cof21, cof22, cof23, cof31, cof32, cof33
      real*8 temp1, temp2, temp3
      shg=shl

      do l=1,npint
         do i=1,3
            do j=1,3
               dx_dr(i,j) = dot_product(shl(j,:,l),xl(i,:))
            end do
         end do

         ! .. definition of the cofactors
         ! .. (recall the definition of matrix inverse : A^-1 = (cof)^T / det A) => (dx_dr)^-1 = (cof)^T / det

         cof11 = dx_dr(2,2)*dx_dr(3,3) - dx_dr(2,3)*dx_dr(3,2)
         cof12 = dx_dr(2,3)*dx_dr(3,1) - dx_dr(2,1)*dx_dr(3,3)
         cof13 = dx_dr(2,1)*dx_dr(3,2) - dx_dr(2,2)*dx_dr(3,1)

         cof21 = dx_dr(3,2)*dx_dr(1,3) - dx_dr(3,3)*dx_dr(1,2)
         cof22 = dx_dr(3,3)*dx_dr(1,1) - dx_dr(3,1)*dx_dr(1,3)
         cof23 = dx_dr(3,1)*dx_dr(1,2) - dx_dr(3,2)*dx_dr(1,1)

         cof31 = dx_dr(1,2)*dx_dr(2,3) - dx_dr(1,3)*dx_dr(2,2)
         cof32 = dx_dr(1,3)*dx_dr(2,1) - dx_dr(1,1)*dx_dr(2,3)
         cof33 = dx_dr(1,1)*dx_dr(2,2) - dx_dr(1,2)*dx_dr(2,1)

         det(l) = dx_dr(1,1)*cof11 + dx_dr(1,2)*cof12 + dx_dr(1,3)*cof13

         if (det(l) < zero) then
            write(iecho,"(a,i10,a,i3)") "Error in shg3D: non-positive determinant in element number ", nel, &
                     " and integration point ", l
            write(*,"(a,i10,a,i3)") "Error in shg3D: non-positive determinant in element number ", nel, " and integration point ", l
            stop ' em shg3d'
         endif

         do i=1,nen

            temp1 = shl(1,i,l)
            temp2 = shl(2,i,l)
            temp3 = shl(3,i,l)

            shg(1,i,l) = (temp1*cof11 + temp2*cof12 + temp3*cof13)/det(l)
            shg(2,i,l) = (temp1*cof21 + temp2*cof22 + temp3*cof23)/det(l)
            shg(3,i,l) = (temp1*cof31 + temp2*cof32 + temp3*cof33)/det(l)
         end do
      end do

   end subroutine shg3D_Josue

end module
