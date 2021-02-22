    module module_polynomials_tri

    implicit none


    contains
    !*************************************************************************
    real function ortho_basis_tri(k,coor,center,area,a )
    ! test functions of dg
    implicit none
    integer,intent(in) :: k
    real,intent(in) :: coor(1:2),center(1:2),area
    real,intent(in) :: a(1:5,1:5)
    real :: x,y,x0,y0

    x = coor(1); y = coor(2);
    x0 = center(1); y0 = center(2);

    if(k .eq. 1)then
        ortho_basis_tri = 1
    elseif(k .eq. 2) then
        ortho_basis_tri = (x-x0)/sqrt( area )
    elseif(k .eq. 3)then
        ortho_basis_tri = a(2,1)*( x-x0 )/sqrt(area)+( y-y0 )/sqrt(area)
    elseif(k .eq. 4)then
        ortho_basis_tri = (x-x0)**2/area+a(3,1)*(x-x0)/sqrt(area)  &
            +a(3,2)*(y-y0)/sqrt(area)+a(3,3)
    elseif(k .eq. 5)then
        ortho_basis_tri = a(4,1)*(x-x0)**2/area+(x-x0)*(y-y0)/area &
            +a(4,2)*(x-x0)/sqrt(area)+a(4,3)*(y-y0)/sqrt(area) + a(4,4)
    elseif(k .eq. 6)then
        ortho_basis_tri = a(5,1)*(x-x0)**2/area &
            + a(5,2)*(x-x0)*(y-y0)/area + (y-y0)**2/area &
            + a(5,3)*(x-x0)/sqrt(area) + a(5,4)*(y-y0)/sqrt(area)+a(5,5)
    else
        write(911,*) 'wrong in ortho_basis_tri'
    endif

    end function ortho_basis_tri
    !*************************************************************************
    !*************************************************************************
    real function poly_tri( n_moment,a, coor , center ,area,oa  )
    implicit none
    integer,intent(in) :: n_moment
    real,intent(in) :: a(1:n_moment)
    real,intent(in) :: coor(1:2)
    real,intent(in) :: center(1:2)
    real,intent(in) :: area,oa(1:5,1:5)

    real :: x,y,x0,y0

    x = coor(1); y = coor(2);
    x0 = center(1); y0 = center(2);


    if( n_moment ==1 )then
        poly_tri = a(1)
    elseif( n_moment == 3 )then
        poly_tri=a(1)+a(2)*( x -x0)/sqrt(area) &
            +a(3)* (  oa(2,1)*(x-x0)/sqrt(area)+(y-y0)/sqrt(area)  )
    elseif( n_moment == 6 )then
        poly_tri = a(1)+a(2)*(x-x0)/sqrt(area) &
            +a(3)* (  oa(2,1)*(x-x0)/sqrt(area)+(y-y0)/sqrt(area)  ) &
            +a(4)*( (x-x0)**2/area+oa(3,1)*(x-x0)/sqrt(area)  &
            +oa(3,2)*(y-y0)/sqrt(area)+oa(3,3) ) &
            +a(5)*(  oa(4,1)*(x-x0)**2/area+(x-x0)*(y-y0)/area &
            +oa(4,2)*(x-x0)/sqrt(area)+oa(4,3)*(y-y0)/sqrt(area) + oa(4,4)) &
            +a(6)*(   oa(5,1)*(x-x0)**2/area &
            + oa(5,2)*(x-x0)*(y-y0)/area + (y-y0)**2/area &
            + oa(5,3)*(x-x0)/sqrt(area) + oa(5,4)*(y-y0)/sqrt(area)+oa(5,5) )
    endif

    end function poly_tri
    !*************************************************************************
    !*************************************************************************
    real function dxh_ortho_basis_tri(k,xh,yh,center,area,a,x1,y1,x2,y2,x3,y3)
    ! test functions of dg
    implicit none
    integer,intent(in) :: k
    real,intent(in) :: xh,yh,center(1:2),area
    real,intent(in) :: a(1:5,1:5)
    real,intent(in) :: x1,y1,x2,y2,x3,y3

    real :: x2m1,x3m1,y2m1,y3m1

    real :: x0,y0

    x0 = center(1);   y0 = center(2);

    x2m1 = x2 - x1
    x3m1 = x3 - x1
    y2m1 = y2 - y1
    y3m1 = y3 - y1

    if(k .eq. 1)then
        dxh_ortho_basis_tri = 0.
    elseif(k .eq. 2) then
        dxh_ortho_basis_tri = x2m1/sqrt(area)
    elseif(k .eq. 3)then
        dxh_ortho_basis_tri = a(2,1)*x2m1/sqrt(area) + y2m1/sqrt(area)
    elseif(k .eq. 4)then
        dxh_ortho_basis_tri = 2*(x2m1*xh+x3m1*yh+x1-x0)*x2m1/area &
            + a(3,1)*x2m1/sqrt(area) &
            + a(3,2)*y2m1/sqrt(area)
    elseif(k .eq. 5)then
        dxh_ortho_basis_tri = a(4,1)*2.*(x2m1*xh+x3m1*yh+x1-x0)*x2m1/area &
            +x2m1*(y2m1*xh+y3m1*yh+y1-y0)/area &
            +(x2m1*xh+x3m1*yh+x1-x0)*y2m1/area &
            + a(4,2)*x2m1/sqrt(area)  &
            + a(4,3)*y2m1/sqrt(area)

    elseif(k .eq. 6)then

        dxh_ortho_basis_tri = a(5,1)*2.*( x2m1*xh+x3m1*yh+x1-x0 )*x2m1/area &
            + a(5,2)/area*( x2m1*(y2m1*xh+y3m1*yh+y1-y0) &
            +( x2m1*xh+x3m1*yh+x1-x0 )*y2m1  ) &
            + 2.*( y2m1*xh+y3m1*yh+y1-y0 )*y2m1/area  &
            + a(5,3)*x2m1/sqrt(area)  &
            + a(5,4)*y2m1/sqrt(area)
    else
        write(911,*) 'wrong in dxh_ortho_basis_tri'
    endif

    end function dxh_ortho_basis_tri
    !*************************************************************************
    real function dyh_ortho_basis_tri(k,xh,yh,center,area,a,x1,y1,x2,y2,x3,y3)
    ! test functions of dg
    implicit none
    integer,intent(in) :: k
    real,intent(in) :: xh,yh,center(1:2),area
    real,intent(in) :: a(1:5,1:5)
    real,intent(in) :: x1,y1,x2,y2,x3,y3
    real :: x0,y0

    real :: x2m1,x3m1,y2m1,y3m1

    x0 = center(1); y0 = center(2);

    x2m1 = x2 - x1
    x3m1 = x3 - x1
    y2m1 = y2 - y1
    y3m1 = y3 - y1

    if(k .eq. 1)then
        dyh_ortho_basis_tri = 0.
    elseif(k .eq. 2) then

        dyh_ortho_basis_tri = x3m1/sqrt(area)
    elseif(k .eq. 3)then

        dyh_ortho_basis_tri = a(2,1)*x3m1/sqrt(area) + y3m1/sqrt(area)
    elseif(k .eq. 4)then

        dyh_ortho_basis_tri = 2.*(x2m1*xh+x3m1*yh+x1-x0)*x3m1/area &
            + a(3,1)*x3m1/sqrt(area) &
            + a(3,2)*y3m1/sqrt(area)
    elseif(k .eq. 5)then


        dyh_ortho_basis_tri = a(4,1)*2.*(x2m1*xh+x3m1*yh+x1-x0)*x3m1/area &
            +x3m1*(y2m1*xh+y3m1*yh+y1-y0)/area &
            +(x2m1*xh+x3m1*yh+x1-x0)*y3m1/area &
            + a(4,2)*x3m1/sqrt(area)  &
            + a(4,3)*y3m1/sqrt(area)

    elseif(k .eq. 6)then

        dyh_ortho_basis_tri = a(5,1)*2.*( x2m1*xh+x3m1*yh+x1-x0 )*x3m1/area &
            + a(5,2)/area*( x3m1*(y2m1*xh+y3m1*yh+y1-y0) &
            +( x2m1*xh+x3m1*yh+x1-x0 )*y3m1  ) &
            + 2.*( y2m1*xh+y3m1*yh+y1-y0 )*y3m1/area  &
            + a(5,3)*x3m1/sqrt(area)  &
            + a(5,4)*y3m1/sqrt(area)
    else
        write(911,*) 'wrong in dyh_ortho_basis_tri'
    endif

    end function dyh_ortho_basis_tri
    !*************************************************************************
    end module module_polynomials_tri