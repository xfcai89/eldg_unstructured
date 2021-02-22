    !*************************************************************************
    ! Example 1:
    !    u_t + u_x + u_y = 0, x\in[-pi,pi], y\in[-pi, pi]
    !    u(x,y,0) = sin(x+y)
    !
    !    u(x,y,t) = sin(x+y-2t).
    !    types of meshes:
    !       regular mesh (to January 15th, 2020, we programmed this.)
    ! Example 2:
    !    rotation with 3 shapes
    ! Example 3:
    !    swirling deformation flow with a cosine bell.
    ! Example 4: 
    !    rotation with u(x,y,0) = exp(-3x^2-3y^2)
    !*************************************************************************
    module module_prob
    implicit none
    integer,public :: ng_tri

    integer,public :: iexample
    real,public :: xleft,xright,ybottom,ytop
    real,public :: tfinal,cfl
    integer,public :: mrk ! for choosing a rk.
    real,public :: rk_d(0:3),trk_stage 
    integer,public :: nghost
    integer,public :: n_moment ! degree of freedom in an element
    real,public :: h

    !***************************************
    integer,public :: nk
    ! gauss points for numerical flux of RKDG.
    integer,public :: ng_line ! number of gauss point for numerical flux
    real,public :: gau_line(1:3,1:2)
    !#######################################
    integer,public :: localrk ! used in characteristic tracing
    real,parameter,private :: pi=4.*atan(1.)


    contains
    subroutine setup_prob
    implicit none

    !******************************************************
    ! method setting
    nk = 2
    ng_line = nk+1
    n_moment = (nk+2)*(nk+1)/2

    localrk = 5

    ng_tri = 9
    !******************************************************
    ! example setting
    iexample = 4
    ! iexample = 1 is linear transport with  sin( x+y  )
    ! iexample = 2 is rotation with three shapes

    tfinal = 2.*pi

    h = 0.1

    !******************************************************
    if(nk==0)then
        mrk = 1
        cfl = 0.8
    elseif(nk==1)then
        cfl = 0.3
        mrk = 2
    elseif(nk==2)then
        cfl = 0.15
        mrk =3
    endif
    nghost = 3
    if(iexample==1)then
        xleft = -pi
        xright = pi
        ybottom = -pi
        ytop = pi
    elseif(iexample==2)then
        xleft = -pi
        xright = pi
        ybottom = -pi
        ytop = pi
    elseif(iexample==3)then
        xleft = -pi
        xright = pi
        ybottom = -pi
        ytop = pi
    elseif(iexample==4)then
        xleft = -pi
        xright = pi
        ybottom = -pi
        ytop = pi
    endif

    !***************************************
    ! gauss points for numerical flux of RKDG.
    if(ng_line==1)then
        gau_line(1,1) = 0.
        gau_line(1,2) = 1.
    elseif(ng_line==2)then
        gau_line(1,1) = -sqrt(3.)/6.
        gau_line(1,2) = 0.5

        gau_line(2,1) = sqrt(3.)/6.
        gau_line(2,2) = 0.5
    elseif(ng_line==3)then
        gau_line(1,1) = -sqrt(15.)/10.
        gau_line(1,2) = 5./18.

        gau_line(2,1) = 0.
        gau_line(2,2) = 4./9.

        gau_line(3,1) = sqrt(15.)/10.
        gau_line(3,2) = 5./18.
    else
        write(911,*)  'gauss'
    endif
    
    if(mrk==1)then
        rk_d(0) = 0.
        rk_d(1) = 1.
    elseif(mrk==2)then
        rk_d(0) = 0.
        rk_d(1) = 1.
        rk_d(2) = 1.
    elseif(mrk==3)then
        rk_d(0) = 0.
        rk_d(1) = 1.
        rk_d(2) = 0.5
        rk_d(3) = 1.
    endif    

    end subroutine setup_prob
    !*************************************************************************
    real function fun_init(coor)
    implicit none
    real,intent(in) :: coor(1:2)

    real :: x,y
    real :: rb,rc,rs
    !real,parameter :: pi=4.*atan(1.)
    real :: pp

    real :: rb0

    x = coor(1); y = coor(2);

    if(iexample==1)then
        fun_init = sin( x+y  )
    elseif(iexample==2)then

        pp = pi*2.

        rb = sqrt( (x+0.25*pp)**2 + (y)**2  )
        rc = sqrt( (x)**2 + (y+0.25*pp)**2 )
        rs = sqrt( (x)**2 + (y-0.25*pp)**2 )
        if(rb < 0.2*pp )then
            fun_init = 0.25*( 1. + cos(pi* rb/(0.2*pp) ) )

        elseif(rc <0.15*pp)then
            fun_init = 1 - rc/(0.15*pp)
        elseif(rs < 0.15*pp)then
            fun_init = 1.

        else
            fun_init = 0.
        endif

        if(x>-0.025*pp .and. x<0.025*pp .and. y>0.1*pp .and. y<0.35*pp)then
            fun_init = 0.
        endif

    elseif(iexample==3)then
        rb = sqrt( (x-0.15/0.5*pi)**2 + (y)**2  )
        rb0 = 0.3*pi
        if(rb < rb0)then
            fun_init =  rb0 * cos( pi*rb/(2.*rb0) )**6
        else
            fun_init = 0.
        endif
    elseif(iexample==4)then
        fun_init = exp( -3.*x**2 - 3.*y**2 )
    endif

    return
    end function fun_init
    !*************************************************************************
    ! begin of flux functions
    real function fun_f(a,coor ,t)
    implicit none
    real,intent(in) :: a,coor(1:2),t
    real :: x,y
    
    x = coor(1); y = coor(2)

    if( iexample==1 )then
        fun_f = a
    elseif( iexample == 2 )then
        fun_f = -y*a
    elseif( iexample == 3 )then
        fun_f = -cos(x/2)**2 * sin(y) * cos(pi*t/tfinal )*pi *a
    elseif( iexample == 4 )then
        fun_f = -y*a
    endif

    end function fun_f
    !*************************************************************************
    real function fun_fp(a,coor ,t)
    implicit none
    real,intent(in) :: a,coor(1:2),t
    real :: x,y
    x = coor(1); y = coor(2);

    if(iexample==1)then
        fun_fp = 1.
    elseif( iexample == 2 )then
        fun_fp = -y
    elseif( iexample == 3 )then
        fun_fp = -cos(x/2)**2 * sin(y) * cos(pi*t/tfinal )*pi
    elseif( iexample == 4 )then
        fun_fp = -y
    endif

    end function fun_fp
    !*************************************************************************
    real function fun_g(a,coor ,t)
    implicit none
    real,intent(in) :: a,coor(1:2),t
    
    real :: x,y
    
    x = coor(1); y = coor(2);

    if(iexample==1)then
        fun_g = a
    elseif( iexample == 2 )then
        fun_g = x*a
    elseif( iexample == 3 )then
        fun_g = sin(x)*cos(y/2.)**2 * cos(pi*t/tfinal )*pi*a 
    elseif( iexample == 4 )then
        fun_g = x*a
    endif

    end function fun_g
    !*************************************************************************
    real function fun_gp(a,coor ,t)
    implicit none
    real,intent(in) :: a,coor(1:2),t
    real :: x,y
    x = coor(1); y = coor(2);

    if(iexample==1)then
        fun_gp = 1.
    elseif( iexample == 2 )then
        fun_gp = x
    elseif( iexample == 3 )then
        fun_gp = sin(x)*cos(y/2.)**2 * cos(pi*t/tfinal )*pi
    elseif( iexample == 4 )then
        fun_gp = x
    endif

    end function fun_gp
    ! end of flux functions
    !*************************************************************************
    !*************************************************************************
    real function fun_exact( coor,t )
    implicit none

    real,intent(in) :: coor(1:2),t

    real :: x,y

    real :: rb,rb0
    !real,parameter :: pi = 4.*atan(1.)

    x = coor(1); y = coor(2);


    if(iexample==1)then
        fun_exact = sin( x+y-2*t  )
    elseif(iexample==2)then

    elseif(iexample==3)then
        rb = sqrt( (x-0.15/0.5*pi)**2 + (y)**2  )
        rb0 = 0.3*pi
        if(rb < rb0)then
            fun_exact =  rb0 * cos( pi*rb/(2.*rb0) )**6
        else
            fun_exact = 0.
        endif
    elseif( iexample==4 )then
        fun_exact = exp( -3.*x**2 -3.*y**2 )
    endif

    return
    end function fun_exact

    !*********************************************************************
    ! the velocity field of passive linea transport problem
    real function  velx( coor,t )
    implicit none

    real,intent(in) :: coor(1:2),t

    real :: x,y

    x = coor(1); y = coor(2)


    if( iexample == 1 )then
        velx = 1.
    elseif( iexample == 2 )then
        velx = -y
    elseif( iexample == 3 )then
        velx= -cos(x/2)**2 * sin(y) * cos(pi*t/tfinal )*pi
    elseif( iexample == 4 )then
        velx = -y 
    endif

    end function velx
    !*********************************************************************
    real function  vely(coor,t )
    implicit none

    real,intent(in) :: coor(1:2),t

    real :: x,y

    x = coor(1); y = coor(2);

    if( iexample == 1 )then
        vely = 1.
    elseif( iexample == 2 )then
        vely = x
    elseif( iexample == 3 )then
        vely = sin(x)*cos(y/2.)**2 * cos(pi*t/tfinal )*pi
    elseif( iexample == 4 )then
        vely = x
    endif

    end function vely
    !***********************************

    end module module_prob