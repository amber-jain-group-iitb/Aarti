      module mod1_tully
      implicit none        
      integer hop
      !potential        
      integer :: i,j,k,l,m,p
      double precision,allocatable::V(:,:),d(:,:),der_V(:,:),E(:,:),F(:)
      double precision ::x,z
      double precision::t,prob_hop,time
      
      !setup
      integer :: n,nstate,ntraj,momentum,newstate,ini_mom
      double precision :: mass
      
      !equation
      double precision :: r_dot,rnd
      complex*16,allocatable :: c_dot(:),cstate(:)
      complex*16,allocatable :: a(:,:)
      double precision,allocatable :: b(:,:)
      double precision :: h
      complex*16,parameter :: iota = (0.d0,1.d0)
      
      double precision, parameter :: xmin=-10d0,xmax=10d0,tmin=0d0

      contains

!-------------------------------------------------------------------------------------------------------------------!              
      subroutine main
      implicit none  
      double precision::prob1_trans,prob2_trans

      !j=0
      !k=0
      call setup
      time=0
      do momentum=8,8         
      j=0
      k=0
      p=0
      m=0
      do l=1,1!ntraj
      call initialize
      do while(x>xmin .and. x<xmax)
      call potential
      call eigenvalues
      call force
      !write(3,*)x,b(1,2),prob_hop,rnd,nstate
      call hoppingprobability
      call rk4
      call velocityverlet 
      time=time+h  
      enddo
      if (nstate==1) then
      if (x<0) then  
      j=j+1
      else 
      k=k+1        
      endif
      endif
      if (nstate==2) then
      if (x<0) then  
      p=p+1
      else 
      m=m+1        
      endif
      endif

      enddo
      prob1_trans=float(j)/float(ntraj)
      prob2_trans=float(k)/float(ntraj)
      write(1,*)momentum,prob1_trans,prob2_trans,float(p)/float(ntraj),float(m)/float(ntraj)
      enddo
      !if (x>=10) then
      !j=j+1        
      !else if (x<=-10) then
      !k=k+1        
      !endif
      !enddo
      !enddo
      !write(*,*)j,k   
      end subroutine main
!--------------------------------------------------------------------------------------------------------------------!
    
      subroutine setup
      character :: ch      
      open(10,file='input1.inp')
      read(10,*)ntraj
      read(10,*)nstate
      read(10,*)newstate
      read(10,*)n
      read(10,*)mass
      read(10,*)h
      read(10,*)ch
      close(10)
       
      if(ch.ne.'y') then
      write(6,*) "problem in reading input file"
      stop
      endif

      allocate(c_dot(nstate),cstate(nstate),a(nstate,nstate))
      allocate(b(nstate,nstate),V(nstate,nstate),D(nstate,nstate))
      allocate(der_V(nstate,nstate),E(nstate,nstate),F(nstate))

      end subroutine setup        
!---------------------------------------------------------------------------------------!
      
      subroutine potential        
      !double precision, parameter :: a=0.01d0,b=1.6d0,c=0.005d0,d=1d0
      double precision :: dx
      
!---------------------------------------------------------------------------------------!
      !problem and its derivative
!---------------------------------------------------------------------------------------!
      
      if (x<0) then
      v(1,1)=-0.01d0*(1-exp(1.6d0*x))
      der_v(1,1)=0.01d0*1.6d0*exp(1.6d0*x)
      else if (x>0) then
      v(1,1)=0.01d0*(1-exp(-1.6d0*x))
      der_v(1,1)=0.01d0*1.6d0*exp(-1.6d0*x)
      endif
      V(2,2)=-V(1,1)
      V(1,2)=0.005d0*exp(-1d0*x**2)
      V(2,1)=V(1,2)
      der_V(2,2)=-der_V(1,1)
      der_V(1,2)=-2*0.005d0*1d0*x*exp(-1d0*x**2)
      der_V(2,1)=der_V(1,2)

!-----------------------------------------------------------------------------------------!
     !coupling
!-----------------------------------------------------------------------------------------!
     z=(V(1,1)-V(2,2))/(V(1,2)*2d0)
      
      d(1,2)=1d0/(2d0*(1+z**2))*(1d0/(2d0*(V(1,2)**2)))*&
      &(V(1,2)*(der_V(1,1)-der_V(2,2))-der_V(1,2)*(V(1,1)-V(2,2)))
      !d(1,2)=d12/50
      
      d(2,1)=-d(1,2)       
      d(1,1)=0
      d(2,2)=0
      end subroutine potential
!------------------------------------------------------------------------------------------!

       Subroutine eigenvalues
       implicit none
      
       E(2,2)=(V(1,1)+V(2,2)+sqrt((V(1,1)+V(2,2))**2-(4*(V(1,1)*V(2,2)&
       &  -(V(1,2))**2))))/2.d0
     
      E(1,1)=(V(1,1)+V(2,2)-sqrt((V(1,1)+V(2,2))**2-(4*(V(1,1)*V(2,2)&
         &  -(V(1,2))**2))))/2.d0

      
      E(1,2)=0
      E(2,1)=0

      end subroutine

!-----------------------------------------------------------------------------------------!      
      Subroutine force
      implicit none        
       !force
      F(1)=(1/sqrt((V(1,1)**2)+(V(1,2)**2)))*((V(1,1)*der_V(1,1))&
           &   +(V(1,2)*der_V(1,2)))
      F(2)=-(1/sqrt((V(1,1)**2)+(V(1,2)**2)))*((V(1,1)*der_V(1,1))&
           &   +(V(1,2)*der_V(1,2)))
      
      end subroutine
!-----------------------------------------------------------------------------------------!
      subroutine equation
       implicit none
       
      call eigenvalues        
              
      c_dot(1)=((E(1,1)/iota))*cstate(1)+&
      & (-(r_dot*d(1,2)))*cstate(2)

      c_dot(2)=((E(2,2)/iota))&
         & *cstate(2)+(-(r_dot*d(2,1)))*cstate(1)

       !return 

      end subroutine equation

!------------------------------------------------------------------------------------------!

      subroutine rk4
      implicit none
      complex*16 :: x1,x2,x3,x4,cstate2 
      complex*16 :: y1,y2,y3,y4,cstate1 
      
      call equation
      
      cstate1=cstate(1)
      cstate2=cstate(2)

      x1=h*c_dot(1)
      y1=h*c_dot(2)
   
      cstate(1)=cstate1+x1/2.0d0
      cstate(2)=cstate2+y1/2.0d0
      
      call equation

      x2=h*c_dot(1)
      y2=h*c_dot(2)

      cstate(1)=cstate1+x2/2.0d0
      cstate(2)=cstate2+y2/2.0d0
      
      x2=x1+2.0d0*x2
      y2=y1+2.0d0*y2

      call equation

      x3=h*c_dot(1)
      y3=h*c_dot(2)
      cstate(1)=cstate1+x3
      cstate(2)=cstate2+y3
      
      x3=x2+2.0d0*x3
      y3=y2+2.0d0*y3

      call equation

      x4=(x3+h*c_dot(1))/6
      y4=(y3+h*c_dot(2))/6

      cstate(1)=cstate1+x4
      cstate(2)=cstate2+y4
      
      end subroutine rk4

!-----------------------------------------------------------------------------------------------------------------!

      subroutine initialize
      implicit none
      x=-9.9999
      r_dot=momentum/mass
      cstate(1)=1
      cstate(2)=0
      nstate=1
      end subroutine initialize

!-----------------------------------------------------------------------------------------------------------------!

      subroutine hoppingprobability
      implicit none
      integer :: i1,j1

      !do i1=1,nstate
      !do j1=1,nstate
      !a(i1,j1)=cstate(i1)*conjg(cstate(j1))
      !b(i1,j1)=2d0*aimag(conjg(a(i1,j1))*E(i1,j1))-&
      !    & 2d0*real(conjg(a(i1,j1))*r_dot*d(i1,j1))
      !enddo
      !enddo 
      
      a(1,2)=cstate(1)*conjg(cstate(2))
      a(2,1)=cstate(2)*conjg(cstate(1))
      a(1,1)=cstate(1)*conjg(cstate(1))
      a(2,2)=cstate(2)*conjg(cstate(2))
      b(1,2)=2d0*aimag(conjg(a(1,2))*E(1,2))&
          & -2d0*real((a(1,2)))*r_dot*d(1,2)
      b(2,1)=2d0*aimag(conjg(a(2,1))*E(2,1))&
           &-2d0*real((a(2,1)))*r_dot*d(2,1)


      if (nstate==1)    then
      if ((0.5d0*mass*(r_dot)**2)>(E(2,2)-E(1,1))) then
     prob_hop=(h*b(2,1))/real(a(1,1))
      call random_number(rnd)
      !write(1,*)x,prob_hop,rnd,nstate
      if (prob_hop>rnd)    then 
      newstate=2
      hop=1
      r_dot=sqrt(2*(E(1,1)-E(2,2)+(0.5*mass*r_dot**2))/mass)
      else
      hop=0 
      nstate=nstate
      endif
      
      else
      nstate=nstate       
      endif      

      else     
      if (nstate==2)    then
      prob_hop=(h*b(1,2))/real(a(2,2))
      call random_number(rnd)
      !write(2,*)x,prob_hop,rnd,nstate
      if (prob_hop>rnd)    then 
      newstate=1
      hop=1
      if (r_dot>0) then
      r_dot=sqrt(2*(E(2,2)-E(1,1)+0.5*mass*r_dot**2)/mass)
      else
      r_dot=-sqrt(2*(E(2,2)-E(1,1)+0.5*mass*r_dot**2)/mass)
      endif        
      else
      hop=0        
      nstate=nstate
      endif
      endif
      endif
        
        if (hop==1)  then 
      nstate=newstate
      endif

      end subroutine 
!!-------------------------------------------------------------------------------------------------!!
   
      subroutine velocityverlet
      implicit none

      i=nstate
      x=x+r_dot*h+(1/(2d0*mass))*F(i)*(h**2)
      r_dot=r_dot+(1/mass)*f(i)*h
      
      write(2,*)x,0.5*mass*r_dot**2+E(i,i) 
      end subroutine
!!-------------------------------------------------------------------------------------------------!! 

     end


