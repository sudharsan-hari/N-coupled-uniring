program main
    implicit none
    integer i,j,count1,k, rlz
    integer, parameter :: num=100,N=3*num,steps=40000000
    real(kind=8):: x(N),fx(N),x12,t=0,h=0.01,m,mean,s,sd,te,te1,y12,x12avg,i1,i2,i3,i4,ic1,ic2,ic3,ic4,ss = 1.0D0
    real(kind=8), dimension (steps) :: x21

!Do rlz = 1,1,1
    !print*, rlz
    !CALL RANDOM_SEED()
    !CALL RANDOM_NUMBER(i1)
    !CALL RANDOM_NUMBER(i2)
    !CALL RANDOM_NUMBER(i3)
    !CALL RANDOM_NUMBER(i4)
	!ic1 = 2.0D0*(2*ss*i1-ss)
	!ic2 = 10.0D0*(2*ss*i2-ss)
	!ic3 = 2.0D0*(2*ss*i3-ss)
	!ic4 = 10.0D0*(2*ss*i4-ss)

    !print*, ic1, ic2, ic3, ic4
	
    do i=1,N,6
        x(i)=1
        x(i+1)=0.5
        x(i+2)=1
        x(i+3)=0.5
        x(i+4)=1
        x(i+5)=0.5
    end do

    open(12, file='spt.txt')
    open(13, file='avg.txt')
    open(14, file='snap.txt')
    open(15, file='osc1.txt')
    open(16, file='osc2.txt')
   
    do i=1,steps,1
        if (mod(i,100000)==0) then
            print*,i/100
        end if
        x12=0
        y12=0
        call ODE_RK4(x,t,h,N)
        t=t+h
        do j=1,N,3
            x12=x12+x(j)
            y12=y12+x(j+1)
        end do
		!if (i >= 3640000 .and. i <=3642000)	then
		! do j=1,N,3
		!  write(12,'(2G17.8,1x,i4)') x(j), t, (j/3)+1
		!  if (t .gt. 356265.00D0 .and. t .lt. 356266.01D0) write (14,'(3G17.8)') x(j), (j/3)+1, t	
		! enddo 
		!endif  
        x21(i)=x12/(num)
        x12avg=x12/(num)
        y12=y12/(num)
        if (i >= 39000000 .and. i <=40000000) then
          write(13,*) t,x21(i)
          write(15,*) t,x(1)
          write(16,*) t,x(4)
		  !do j=1,N,3
		    
		  !end do
        end if     
    end do
    !m=0
    !do k=100000000,steps,1
    !    m=m+x21(k)
    !end do
    !mean=m/steps
    !print*,"Mean=",mean
    !s=0
    !do k=100000000,steps,1
    !    s=s+(x21(k)-mean)**2
    !end do

    !sd=sqrt(s/steps)
    !print*,"sd=",sd
    !te=mean+6*sd
    !print*, "TE=",te
    !te1=mean-6*sd
    !print*,"TE1=",te1
    !count1=0
    !do k=100000000,steps,1
    !    if (x21(k)>te) then
    !        count1=count1+1
    !    end if
    !end do
    !WRITE(13,*) rlz, count1 
!ENDDO    

end program main

subroutine ODE_RK4(x,t,h,n)
    implicit none
    integer i,n
    real(kind=8):: k1(n),k2(n),k3(n),k4(n),k5(n),k6(n),temp(n),fx(n),x(n),t,h

    call rhs(x,t,fx,n)
    do i=1,n
        k1(i)=h*fx(i)
        temp(i)=x(i)+((0.25)*(k1(i)))
    end do

    call rhs(temp,t+((0.25)*h),fx,N)
    do i=1,n
        k2(i)=h*fx(i)
        temp(i)=x(i)+((0.09375)*k1(i))+((0.28125)*k2(i))
    end do

    call rhs(temp,t+((0.375)*h),fx,N)
    do i=1,n
        k3(i)=h*fx(i)
        temp(i)=x(i)+((0.8793809741)*k1(i))-((3.2771961766)*k2(i))+((3.3208921256)*k3(i))
    end do

    call rhs(temp,t+((0.9230769231)*h),fx,N)
    do i=1,n
        k4(i)=h*fx(i)
        temp(i)=x(i)+((2.0324074074)*k1(i))+((-8)*k2(i))+((7.1734892788)*k3(i))-((0.2058966862)*k4(i))
    end do

    call rhs(temp,t+h,fx,N)
    do i=1,n
        k5(i)=h*fx(i)
        temp(i)=x(i)-((0.2962962963)*k1(i))+(2*k2(i))-((1.3816764133)*k3(i))+((0.4529727096)*k4(i))-((0.275)*k5(i))
    end do

    call rhs(temp,t+((0.5)*h),fx,N)
    do i=1,n
        k6(i)=h*fx(i)
        x(i)=x(i)+(0.1157407407)*k1(i)+(0.5489278752)*k3(i)+(0.5357229944)*k4(i)-(0.2)*k5(i)
    end do
end subroutine ODE_RK4

subroutine rhs(x,t,f,N)
    integer k,N
    real(kind=8):: x(N),t,f(N)
    real(kind=8):: a=1,b=3,c=1,d=5,xr=-1.6,r=0.01,s=5,I=4,vs=2,ga,l1=0.0

        do k=1,N,3
            f(k)=x(k+1)+(b*x(k)**2)-(a*(x(k)**3))-x(k+2)+I-(l1*(x(k)-vs)*ga(k,x,n))
            f(k+1)=c-(d*x(k)**2)-x(k+1)
            f(k+2)=r*(s*(x(k)-xr)-x(k+2))
        end do
    

   
end subroutine

function ga(k,x,N)
    implicit none
    integer k, n,i,r
    real(kind=8):: ag1,x(N),lam=10,theta=-0.25,ga,ag2,ra=1
    r=3*ra
    ag1=0
    ag2=0
    do i=k+3,k+r,3
        if (i<=N) then
            ag1=ag1+(1/(1+exp(-lam*(x(i)-theta))))
        end if
        if (i>N) then
            ag2=ag2+1/(1+exp(-lam*(x(i-N)-theta)))
        end if
    end do
    ga=(ag1+ag2)/ra
end function

