    !*************************************************************************
    module module_gauss_tri
    ! the number of gauss quadrature points was set in module_prob
    use module_prob,only : ng_tri ! the number of gauss quadrature points!
    implicit none
    
    real,allocatable :: w(:),sx(:),sy(:)
    !!************************************************
    !! 6 quadrature points that is exact for polynomials
    !! of total degree less than 4
    !! which can be found in page of book
    !! Discontinuous Galerkin Methods for solving Elliptic and Parabolic
    !! Equations: Theory and Implementation.
    !REAL, PARAMETER, DIMENSION(1:6) :: w = (/ &
    !    0.11169079483901,&
    !    0.11169079483901,&
    !    0.11169079483901,&
    !    0.05497587182766,&
    !    0.05497587182766,&
    !    0.05497587182766/)
    !
    !REAL, PARAMETER, DIMENSION(1:6) :: sx = (/ &
    !    0.445948490915965,&
    !    0.108103018168070,&
    !    0.445948490915965,&
    !    0.091576213509771,&
    !    0.816847572980459,&
    !    0.091576213509771/)
    !
    !REAL, PARAMETER, DIMENSION(1:6) :: sy = (/ &
    !    0.445948490915965,&
    !    0.445948490915965,&
    !    0.108103018168070,&
    !    0.091576213509771,&
    !    0.091576213509771,&
    !    0.816847572980459/)
    contains
    !*************************************************************************
    ! quadrature rules for triangles.
    !
    ! quadrature rules for the triangle are generally defined over T1,
    ! the unit triangle whose vertices are (0,0), (1,0), (0,1).
    !
    ! from http://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html
    !
    ! References:
    ! 1. Jarle Berntsen, Terje Espelid,
    ! Algorithm 706,
    ! DCUTRI: an algorithm for adaptive cubature over a collection of triangles,
    ! ACM Transactions on Mathematical Software,
    ! Volume 18, Number 3, September 1992, pages 329-342.
    !
    ! 2. Elise deDoncker, Ian Robinson,
    ! Algorithm 612: Integration over a Triangle Using Nonlinear Extrapolation,
    ! ACM Transactions on Mathematical Software,
    ! Volume 10, Number 1, March 1984, pages 17-22.
    !
    ! 3. Dirk Laurie,
    ! Algorithm 584, CUBTRI, Automatic Cubature Over a Triangle,
    ! ACM Transactions on Mathematical Software,
    ! Volume 8, Number 2, 1982, pages 210-218.
    !
    ! 4. James Lyness, Dennis Jespersen,
    ! Moderate Degree Symmetric Quadrature Rules for the Triangle,
    ! Journal of the Institute of Mathematics and its Applications,
    ! Volume 15, Number 1, February 1975, pages 19-32.
    !
    ! 5. Hans Rudolf Schwarz,
    ! Finite Element Methods,
    ! Academic Press, 1988,
    ! ISBN: 0126330107,
    ! LC: TA347.F5.S3313.
    !
    ! 6. Gilbert Strang, George Fix,
    ! An Analysis of the Finite Element Method,
    ! Cambridge, 1973,
    ! ISBN: 096140888X,
    ! LC: TA335.S77.
    !
    ! 7. Arthur Stroud,
    ! Approximate Calculation of Multiple Integrals,
    ! Prentice Hall, 1971,
    ! ISBN: 0130438936,
    ! LC: QA311.S85.
    !
    ! 8. Olgierd Zienkiewicz,
    ! The Finite Element Method,
    ! Sixth Edition,
    ! Butterworth-Heinemann, 2005,
    ! ISBN: 0750663200,
    ! LC: TA640.2.Z54
    ! 9. David Dunavant,
    ! High Degree Efficient Symmetrical Gaussian Quadrature Rules for the Triangle,
    ! International Journal for Numerical Methods in Engineering,
    ! Volume 21, 1985, pages 1129-1148.
    !*************************************************************************
    subroutine quadrature_rules_tri
    implicit none


    if( ng_tri==1 )then
        w = (/  0.5 /)
        sx = (/  1./3. /)
        sy = (/  1./3. /)
    elseif( ng_tri==6 )then
        !************************************************
        ! 6 quadrature points that is exact for polynomials
        ! of total degree less than 4
        ! which can be found in page of book
        ! Discontinuous Galerkin Methods for solving Elliptic and Parabolic
        ! Equations: Theory and Implementation.
        !
        ! David Dunavant,
        !High Degree Efficient Symmetrical Gaussian Quadrature Rules for the Triangle,
        !International Journal for Numerical Methods in Engineering,
        !Volume 21, 1985, pages 1129-1148.
        w = (/ &
            0.11169079483901,&
            0.11169079483901,&
            0.11169079483901,&
            0.05497587182766,&
            0.05497587182766,&
            0.05497587182766/)

        sx = (/ &
            0.445948490915965,&
            0.108103018168070,&
            0.445948490915965,&
            0.091576213509771,&
            0.816847572980459,&
            0.091576213509771/)

        sy = (/ &
            0.445948490915965,&
            0.445948490915965,&
            0.108103018168070,&
            0.091576213509771,&
            0.091576213509771,&
            0.816847572980459/)
    elseif( ng_tri==9 )then

        w = (/ &
            0.205950504760887, &
            0.205950504760887, &
            0.205950504760887, &
            0.063691414286223, &
            0.063691414286223, &
            0.063691414286223, &
            0.063691414286223, &
            0.063691414286223, &
            0.063691414286223 /)
        w = 0.5*w

        sx = (/ &
            0.124949503233232, &
            0.437525248383384, &
            0.437525248383384, &
            0.797112651860071, &
            0.797112651860071, &
            0.165409927389841, &
            0.165409927389841, &
            0.037477420750088, &
            0.037477420750088 /)

        sy = (/ &
            0.437525248383384, &
            0.124949503233232, &
            0.437525248383384, &
            0.165409927389841, &
            0.037477420750088, &
            0.797112651860071, &
            0.037477420750088, &
            0.797112651860071,&
            0.165409927389841         /)
    endif
    end subroutine quadrature_rules_tri
    !*************************************************************************
    subroutine allocate_gau_para
    implicit none


    allocate( w(ng_tri) )
    allocate( sx(ng_tri) )
    allocate( sy(ng_tri) )

    end subroutine allocate_gau_para
    !*************************************************************************
    subroutine deallocate_gau_para
    implicit none

    deallocate( w )
    deallocate( sx )
    deallocate( sy )

    end subroutine deallocate_gau_para
    !*************************************************************************
    end module module_gauss_tri