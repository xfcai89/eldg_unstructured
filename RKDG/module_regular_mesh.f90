    !*************************************************************************
    ! please set up the parameters regarding the rectangular compuational
    ! domain
    !             [xleft,xright]\times[ybottom,ytop]
    ! in subroutine parameters_mesh
    !-------------------------------------------------------------------------
    ! Note that the ghost cells
    ! is formed by the following precedence,
    ! first x direction;
    ! and then y direction.
    !*************************************************************************
    module module_regular_mesh
    use module_prob,only : gau_line,ng_line,n_moment,nghost
    use module_prob,only : xleft,xright,ybottom,ytop ! computational domain
  
    
    use module_regular_mesh_data
  
    implicit none
 
    integer :: nx,ny
 
    integer,public :: n_node   ! number of nodes
    integer,public :: n_element ! number of elements.
 
    integer :: n_side ! number of sides.
    integer :: n_node_all,n_element_all,n_side_all
    real :: dx,dy

    real,public :: radius !
    real,public :: speed
 
 

    contains
    !
 
    !*************************************************************************
    subroutine setup_trimesh(kkkk )
    use module_prob, only : iexample
    implicit none
    integer,intent(in) :: kkkk
    real,parameter :: pi= 4.*atan(1.)


    nx = 10*2**(kkkk-1)
    ny = 10*2**(kkkk-1)  

    !nx = 4
    !ny = 4
    ! for generating nx*ny*2 elements.
    !nghost = 1
    !nghost = 3


    !***********************************
    dx = (xright-xleft)/nx
    dy = (ytop-ybottom)/ny

    n_node = (nx+1)*(ny+1)
    n_element = 2*nx*ny
    n_side = nx*(ny+1) + ny*(nx+1) + nx*ny

    n_node_all = ( nx+1+2*nghost )*( ny+1+2*nghost )
    n_element_all = 2.*( nx+2*nghost )*( ny+2*nghost )
    n_side_all = (nx+2*nghost)*(ny+2*nghost+1) &
        + (ny+2*nghost)*(nx+2*nghost+1) + (nx+2*nghost)*(ny+2*nghost )

    end subroutine setup_trimesh
    !*************************************************************************
    subroutine allocate_trimesh
    implicit none


    allocate( node(1:n_node_all ) )
    allocate(element(1:n_element_all) )
    !allocate( side(1:n_side) )
    allocate( side(1:n_side_all) )
    !--------------------------------------
    allocate( side_lower(1:nx,1:nghost) )
    allocate( side_upper(1:nx,1:nghost) )
    allocate( side_left(1-nghost:ny+nghost,1:nghost) )
    allocate( side_right(1-nghost:ny+nghost,1:nghost) )

    end subroutine allocate_trimesh
    !*************************************************************************
    subroutine deallocate_trimesh
    implicit none

    deallocate( node )
    deallocate( element )
    deallocate( side )
    !--------------------------------------
    deallocate( side_lower )
    deallocate( side_upper )
    deallocate( side_left )
    deallocate( side_right )

    end subroutine deallocate_trimesh
    !*************************************************************************
    subroutine get_nodes_elements
    implicit none
    integer :: i,j,ip,iele,iside

    integer :: n_node_s,n_node_t
    integer :: n_element_s,n_element_t

    !allocate( point(1:n_node ) )

    ! Get nodes in the computational domain
    do j = 1 ,ny+1
        do i = 1,nx+1
            ip = (j-1)*(nx+1) + i
            node(ip)%coor(1) = xleft + (i-1)*dx
            node(ip)%coor(2) = ybottom + (j-1)*dy
        enddo
    enddo

    !*************************************************************************
    !allocate(element(1:n_element) )


    iele = 0
    ! get elements, anti-clockwise
    do j = 1,ny
        do i = 1,nx
            iele = iele + 1
            element(iele)%node(1) = (j-1)*(nx+1) + i
            element(iele)%node(2) = (j-1)*(nx+1) + i + 1
            element(iele)%node(3) = (j)*(nx+1) + i

            element(iele)%mark = 1
            iele = iele + 1
            !
            element(iele)%node(1) = (j-1)*(nx+1) + i + 1
            element(iele)%node(2) = (j)*(nx+1) + i + 1
            element(iele)%node(3) = (j)*(nx+1) + i

            element(iele)%mark = 1
        enddo
    enddo

    !call output_tri(n_node,n_element)

    !*************************************************
    !*************************************************
    ! for forming ghost cells
    !allocate( side_lower(1:nx,1:nghost) )
    !allocate( side_upper(1:nx,1:nghost) )
    do j = 1,nghost
        do i = 1 ,nx
            side_lower(i,j)%node(1) = (j-1)*(nx+1) + i
            side_lower(i,j)%node(2) = (j-1)*(nx+1) + i +1

            side_lower(i,j)%ie(1) = ( i + (j-1)*nx )*2-1
            side_lower(i,j)%ie(2) = ( i + (j-1)*nx )*2
            !*******************
            side_upper(i,j)%node(1) =  (ny+1-j)*(nx+1) + i
            side_upper(i,j)%node(2) =  (ny+1-j)*(nx+1) + i + 1

            side_upper(i,j)%ie(1) = ( i + (ny-j)*nx )*2 - 1
            side_upper(i,j)%ie(2) = ( i + (ny-j)*nx )*2
        enddo
    enddo

    !allocate( side_left(1-nghost:ny+nghost,1:nghost) )
    !allocate( side_right(1-nghost:ny+nghost,1:nghost) )

    do i = 1,nghost
        do j = 1,ny

            side_left(j,i)%node(1) = (j-1)*(nx+1) + i
            side_left(j,i)%node(2) = (j)*(nx+1) + i
            side_left(j,i)%ie(1) = (i+(j-1)*nx)*2 -1
            side_left(j,i)%ie(2) = (i+(j-1)*nx)*2

            !************************************************
            side_right(j,i)%node(1) = (j-1)*(nx+1) + (nx+2-i)
            side_right(j,i)%node(2) = (j)*(nx+1) + (nx+2-i)
            side_right(j,i)%ie(1) = (nx+1-i+ (j-1)*nx)*2 -1
            side_right(j,i)%ie(2) = (nx+1-i+ (j-1)*nx)*2
        enddo
    enddo

    !**********************************************************************
    n_node_s = n_node
    n_element_s = n_element

    call add_upper_ghost(n_node_s,n_node_t,dy,n_element_s,n_element_t)
    !call output_tri(n_node_t,n_element_t)
    !**********************************************************************
    n_node_s = n_node_t
    n_element_s = n_element_t
    call add_lower_ghost(n_node_s,n_node_t,dy,n_element_s,n_element_t)
    !call output_tri(n_node_t,n_element_t)
    !**********************************************************************
    n_node_s = n_node_t
    n_element_s = n_element_t
    call add_left_ghost(n_node_s,n_node_t,dx,n_element_s,n_element_t)
    !call output_tri(n_node_t,n_element_t)
    !**********************************************************************
    n_node_s = n_node_t
    n_element_s = n_element_t
    call add_right_ghost(n_node_s,n_node_t,dx,n_element_s,n_element_t)
    !call output_tri(n_node_t,n_element_t)



    end subroutine get_nodes_elements
    !*************************************************************************
    subroutine add_upper_ghost(n_node_s,n_node_t,dy,n_element_s,n_element_t)
    implicit none
    integer,intent(in) :: n_node_s
    integer,intent(out) :: n_node_t
    real,intent(in) :: dy
    integer,intent(in) :: n_element_s
    integer,intent(out) :: n_element_t

    integer :: j,i,ip,nnt,iele,ii
    type(type_side_b), allocatable :: side_lower_copy(:,:)

    integer :: jj,nn

    allocate( side_lower_copy(1:nx,1:nghost+1) )

    do i = 1,nx
        side_lower_copy(i,1)%node(1) = side_upper(i,1)%node(1)
        side_lower_copy(i,1)%node(2) = side_upper(i,1)%node(2)
        !side_lower_copy(i,1)%ie(1) = side_lower(i,1)%ie(1)
        !side_lower_copy(i,1)%ie(2) = side_lower(i,1)%ie(2)
    enddo

    ip = n_node_s
    iele = n_element_s
    do j = 1,nghost
        do i = 1 ,nx
            !********
            ip = ip +1
            nnt = side_lower_copy(i,j)%node(1)

            node(ip)%coor(1) =  node(nnt)%coor(1)
            node(ip)%coor(2) =  node(nnt)%coor(2) + dy
            !*********
        enddo
        ip = ip +1
        nnt = side_lower_copy(nx,j)%node(2)
        node(ip)%coor(1) = node(nnt)%coor(1)
        node(ip)%coor(2) = node(nnt)%coor(2) + dy
        !---------------------------------------------------
        n_node_t = ip
        ii = n_node_t - (nx+1)
        !if( j<nghost )then
        do i = 1,nx
            side_lower_copy(i,j+1)%node(1) =  ii +i
            side_lower_copy(i,j+1)%node(2) =  ii +i+1
        enddo
        !endif
        !------------------------------------------------------
        !add new elements
        do i = 1,nx
            iele = iele + 1
            element(iele)%node(1) = side_lower_copy(i,j )%node(1)
            element(iele)%node(2) = side_lower_copy(i,j )%node(2)
            element(iele)%node(3) = ii+i

            element(iele)%mark = 0

            element(iele)%i_periodic = side_lower(i,j)%ie(1)

            !##############################################################
            if(i<=nghost)then
                jj = ny + j
                side_left(jj,i)%ie(1) = iele
            elseif( i>nx-nghost )then
                nn = nx+1-i
                jj = ny + j
                side_right(jj,nn)%ie(1) =iele
            else

            endif
            !##############################################################
            iele = iele + 1

            !
            element(iele)%node(1) = side_lower_copy(i,j )%node(2)
            element(iele)%node(2) = ii +i+1
            element(iele)%node(3) = ii + i

            element(iele)%mark = 0

            element(iele)%i_periodic = side_lower(i,j)%ie(2)

            !##############################################################
            if(i<=nghost)then
                jj = ny + j
                side_left(jj,i)%ie(2) = iele
            elseif( i>nx-nghost )then
                nn = nx+1-i
                jj = ny + j
                side_right(jj,nn)%ie(2) =iele
            else
            endif
            !##############################################################
        enddo
    enddo

    n_node_t = n_node_s + (nx+1) *nghost
    n_element_t = n_element_s + nx*nghost*2


    !*******************************************************
    ! update side_left and side_right
    do i = 1,nghost
        do j = ny+1,ny+nghost
            jj = j-ny
            side_left(j,i)%node(1) = side_lower_copy(i,jj)%node(1)

            side_left(j,i)%node(2) = side_lower_copy(i,jj+1)%node(1)


            !**********************************************************
            side_right(j,i)%node(1) = side_lower_copy(nx+1-i,jj)%node(2)

            side_right(j,i)%node(2) = side_lower_copy(nx+1-i,jj+1)%node(2)

        enddo
    enddo


    deallocate(side_lower_copy)

    end subroutine add_upper_ghost
    !*************************************************************************
    subroutine add_lower_ghost(n_node_s,n_node_t,dy,n_element_s,n_element_t)
    implicit none
    integer,intent(in) :: n_node_s
    integer,intent(out) :: n_node_t
    real,intent(in) :: dy
    integer,intent(in) :: n_element_s
    integer,intent(out) :: n_element_t

    integer :: i,j,ip,iele,nnt,ii
    type(type_side_b), allocatable :: side_upper_copy(:,:)

    integer :: jj,nn

    allocate( side_upper_copy(1:nx,1:nghost+1) )

    do i = 1,nx
        side_upper_copy(i,1)%node(1) = side_lower(i,1)%node(1)
        side_upper_copy(i,1)%node(2) = side_lower(i,1)%node(2)
        !side_upper_copy(i,1)%ie(1) = side_upper(i,1)%ie(1)
        !side_upper_copy(i,1)%ie(2) = side_upper(i,1)%ie(2)
    enddo

    ip = n_node_s
    iele = n_element_s
    do j = 1,nghost
        do i = 1 ,nx
            !********
            ip = ip +1
            nnt = side_upper_copy(i,j)%node(1)

            node(ip)%coor(1) =  node(nnt)%coor(1)
            node(ip)%coor(2) =  node(nnt)%coor(2) - dy
            !*********
        enddo
        ip = ip +1
        nnt = side_upper_copy(nx,j)%node(2)
        node(ip)%coor(1) = node(nnt)%coor(1)
        node(ip)%coor(2) = node(nnt)%coor(2) - dy
        !---------------------------------------------------
        n_node_t = ip
        ii = n_node_t - (nx+1)
        !if( j<nghost )then
        do i = 1,nx
            side_upper_copy(i,j+1)%node(1) =  ii +i
            side_upper_copy(i,j+1)%node(2) =  ii +i+1
        enddo
        !endif
        !------------------------------------------------------
        !add new elements
        do i = 1,nx
            iele = iele + 1
            element(iele)%node(1) = side_upper_copy(i,j )%node(1)
            element(iele)%node(2) = ii + i
            element(iele)%node(3) = ii+i +1

            element(iele)%mark = 0

            element(iele)%i_periodic = side_upper(i,j )%ie(1)

            !##############################################################
            if(i<=nghost)then
                jj = 1-j
                side_left(jj,i)%ie(1) = iele
            elseif( i>nx-nghost )then
                jj = 1-j
                nn = nx+1-i
                side_right(jj,nn)%ie(1) =iele
            else

            endif
            !##############################################################
            iele = iele + 1
            !
            element(iele)%node(1) = side_upper_copy(i,j )%node(1)
            element(iele)%node(2) = ii +i+1
            element(iele)%node(3) =   side_upper_copy(i,j )%node(2)

            element(iele)%mark = 0

            element(iele)%i_periodic = side_upper(i,j )%ie(2)

            !##############################################################
            if(i<=nghost)then
                jj = 1-j
                side_left(jj,i)%ie(2) = iele
            elseif( i>nx-nghost )then
                nn = nx+1-i
                jj = 1-j
                side_right(jj,nn)%ie(2) =iele
            else
            endif
            !##############################################################
        enddo
    enddo

    n_node_t = n_node_s + (nx+1) *nghost
    n_element_t = n_element_s + nx*nghost*2

    !*******************************************************
    ! update side_left and side_right
    do i = 1,nghost
        do j =  1-nghost,1-1
            jj = 1-j

            side_left(j,i)%node(1) = side_upper_copy(i,jj+1)%node(1)

            side_left(j,i)%node(2) = side_upper_copy(i,jj)%node(1)
            !**********************************************************
            side_right(j,i)%node(1) = side_upper_copy(nx+1-i,jj+1)%node(2)

            side_right(j,i)%node(2) = side_upper_copy(nx+1-i,jj)%node(2)
        enddo
    enddo

    deallocate( side_upper_copy )

    end subroutine add_lower_ghost
    !*************************************************************************
    subroutine add_left_ghost(n_node_s,n_node_t,dx,n_element_s,n_element_t)
    implicit none
    integer,intent(in) :: n_node_s
    integer,intent(out) :: n_node_t
    real,intent(in) :: dx
    integer,intent(in) :: n_element_s
    integer,intent(out) :: n_element_t

    integer :: i,j,ip,iele,nnt
    type(type_side_b), allocatable :: side_right_copy(:,:)
    integer :: ii


    allocate( side_right_copy(1-nghost:ny+nghost,1:nghost+1 ) )
    do j = 1-nghost,ny+nghost
        side_right_copy(j,1)%node(1) = side_left(j,1)%node(1)
        side_right_copy(j,1)%node(2) = side_left(j,1)%node(2)
    enddo


    ip = n_node_s
    iele = n_element_s
    do i = 1,nghost
        do j = 1-nghost,ny+nghost
            !********
            ip = ip +1
            nnt = side_right_copy(j,i)%node(1)

            node(ip)%coor(1) =  node(nnt)%coor(1) - dx
            node(ip)%coor(2) =  node(nnt)%coor(2)
            !*********
        enddo
        ip = ip +1
        nnt = side_right_copy(ny+nghost,i)%node(2)
        node(ip)%coor(1) = node(nnt)%coor(1) - dx
        node(ip)%coor(2) = node(nnt)%coor(2)
        !---------------------------------------------------
        n_node_t = ip
        ii = n_node_t - (ny+nghost+1)
        !if( j<nghost )then
        do j= 1-nghost,ny+nghost
            side_right_copy(j,i+1)%node(1) =  ii +j
            side_right_copy(j,i+1)%node(2) =  ii +j+1
        enddo
        !endif
        !------------------------------------------------------
        !add new elements
        do j = 1-nghost,ny+nghost
            iele = iele + 1
            element(iele)%node(1) = ii + j
            element(iele)%node(2) = side_right_copy(j,i )%node(1)
            element(iele)%node(3) = ii+j+1

            element(iele)%mark = 0

            element(iele)%i_periodic = side_right(j,i)%ie(1)

            !****************************************************
            iele = iele + 1
            !
            element(iele)%node(1) = ii + j + 1
            element(iele)%node(2) = side_right_copy(j,i )%node(1)
            element(iele)%node(3) = side_right_copy(j,i )%node(2)

            element(iele)%mark = 0

            element(iele)%i_periodic = side_right(j,i )%ie(2)

        enddo
    enddo

    n_node_t = n_node_s + (ny+2*nghost+1) *nghost
    n_element_t = n_element_s +(ny+2*nghost ) *nghost*2

    deallocate( side_right_copy )

    end subroutine add_left_ghost
    !*************************************************************************
    subroutine add_right_ghost(n_node_s,n_node_t,dx,n_element_s,n_element_t)
    implicit none
    integer,intent(in) :: n_node_s
    integer,intent(out) :: n_node_t
    real,intent(in) :: dx
    integer,intent(in) :: n_element_s
    integer,intent(out) :: n_element_t

    integer :: i,j,ip,iele,nnt
    type(type_side_b), allocatable :: side_left_copy(:,:)
    integer :: ii

    allocate( side_left_copy(1-nghost:ny+nghost,1:nghost+1 ) )
    do j = 1-nghost,ny+nghost
        side_left_copy(j,1)%node(1) = side_right(j,1)%node(1)
        side_left_copy(j,1)%node(2) = side_right(j,1)%node(2)
    enddo

    ip = n_node_s
    iele = n_element_s
    do i = 1,nghost
        do j = 1-nghost,ny+nghost
            !********
            ip = ip +1
            nnt = side_left_copy(j,i)%node(1)

            node(ip)%coor(1) =  node(nnt)%coor(1) + dx
            node(ip)%coor(2) =  node(nnt)%coor(2)
            !*********
        enddo
        ip = ip +1
        nnt = side_left_copy(ny+nghost,i)%node(2)
        node(ip)%coor(1) = node(nnt)%coor(1) + dx
        node(ip)%coor(2) = node(nnt)%coor(2)
        !---------------------------------------------------
        n_node_t = ip
        ii = n_node_t - (ny+nghost+1)
        !if( j<nghost )then
        do j= 1-nghost,ny+nghost
            side_left_copy(j,i+1)%node(1) =  ii +j
            side_left_copy(j,i+1)%node(2) =  ii +j+1
        enddo
        !endif
        !------------------------------------------------------
        !add new elements
        do j = 1-nghost,ny+nghost
            iele = iele + 1
            element(iele)%node(1) = side_left_copy(j,i )%node(2)
            element(iele)%node(2) = side_left_copy(j,i )%node(1)
            element(iele)%node(3) = ii+j

            element(iele)%mark = 0

            element(iele)%i_periodic = side_left(j,i)%ie(1)

            !*********************************************************
            iele = iele + 1
            !
            element(iele)%node(1) =  side_left_copy(j,i )%node(2)
            element(iele)%node(2) = ii+j
            element(iele)%node(3) = ii+j+1

            element(iele)%mark = 0

            element(iele)%i_periodic = side_left(j,i)%ie(2)

        enddo
    enddo

    n_node_t = n_node_s + (ny+2*nghost+1) *nghost
    n_element_t = n_element_s +(ny+2*nghost ) *nghost*2
    deallocate( side_left_copy )

    end subroutine add_right_ghost
    !*************************************************************************
    !*************************************************************************
    !*************************************************************************
    !*************************************************************************
    subroutine get_tri_info(n_ele)
    use module_prob,only : ng_tri
    use module_gauss_tri, only : sx,sy,w

    use module_gau_elimi
    implicit none
    integer,intent(in) :: n_ele
    integer :: iele
    integer :: n1,n2,n3
    real :: x1,y1,x2,y2,x3,y3,det
    real :: x2m1,x3m1,y2m1,y3m1
    integer :: ig
    real :: x0,y0

    real :: xh,yh
    real :: xmx0_map,ymy0_map

    real :: xy,xs,ys,int_xy,int_ys,int_xs

    real ::     xc,    yc,    xsv3
    real :: int_xc,int_yc,int_xsv3

    real :: aa(1:5,1:5),b(1:5),xx(1:5)

    real :: xcy,        xsy,    xf,    xys,    xsys,    xyc
    real :: int_xcy,int_xsy,int_xf,int_xys,int_xsys,int_xyc
    real :: a21
    real :: a31,a32,a33
    real :: a41,a42,a43,a44
    real :: a51,a52,a53,a54,a55

    real :: ar,v4v4
    real :: v5v5,v6v6

    do iele = 1,n_ele

        n1 = element( iele )%node(1)
        n2 = element( iele )%node(2)
        n3 = element( iele )%node(3)

        x1 = node(n1)%coor(1)
        y1 = node(n1)%coor(2)

        x2 = node(n2)%coor(1)
        y2 = node(n2)%coor(2)

        x3 = node(n3)%coor(1)
        y3 = node(n3)%coor(2)

        det = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)

        element(iele)%area = 0.5*det

        ! afine mapping
        x2m1 = x2 - x1
        x3m1 = x3 - x1
        y2m1 = y2 - y1
        y3m1 = y3 - y1

        !**************************************
        ! mapping
        !   | x |   =  | x2-x1    x3-x1 | | hx |  +  | x1 |
        !   | y |      | y2-y1    y3-y1 | | hy |     | y1 |

        do ig = 1, ng_tri
            element(iele)%gauss(ig,1) = x2m1*sx(ig) + x3m1*sy(ig) + x1
            element(iele)%gauss(ig,2) = y2m1*sx(ig) + y3m1*sy(ig) + y1
        enddo


        !**************************************************************
        ! get bary center of a triangle!
        element(iele)%bary(1) = (x1+x2+x3)/3.
        element(iele)%bary(2) = (y1+y2+y3)/3.


        !
        element(iele)%jac_inv(1,1) = (y3-y1)/det
        element(iele)%jac_inv(1,2) = -(x3-x1)/det
        element(iele)%jac_inv(2,1) = -(y2-y1)/det
        element(iele)%jac_inv(2,2) = (x2-x1)/det

        !**************************************************************
        ! get the coefficient of the local orthogonal basis
        ! p1
        ! we adopt a local orthogonal basis over a target element Delta0:
        ! v0 = 1,
        ! v1 = ( x-x0 )/sqrt(area),
        ! v2 = a21*( x-x0 )/sqrt(area)+( y-y0 )/sqrt(area).
        !
        ! Let
        !           X = ( x-x0 )/sqrt(area)
        !           Y = ( y-y0 )/sqrt(area)
        !
        !  xy = [ ( x-x0 )/sqrt(area) ] * [ ( y-y0 )/sqrt(area) ]
        !  xs = [ ( x-x0 )/sqrt(area) ] * [ ( x-x0 )/sqrt(area) ]
        !  ys = [ ( y-y0 )/sqrt(area) ] * [ ( y-y0 )/sqrt(area) ]
        ! Then
        ! int_xy = int_Delta0 xy dxdy
        ! int_xs = int_Delta0 xs dxdy
        ! int_ys = int_Delta0 ys dxdy
        !
        ! a21 = - int_xy/int_xs
        !
        ! Mass matrix
        ! diag( w1,w2,w3 )

        !                  w1 = area
        !                  w2 = int_xs
        !                  w3 = a21^2*int_xs  + 2a21*int_xy + int_ys
        !-------------------------------------------------------------------
        ! Numerical quadratures for
        ! int_xy = int_Delta0 xy dxdy
        ! int_xs = int_Delta0 xs dxdy
        ! int_ys = int_Delta0 ys dxdy
        !
        ! map to xh yh
        !  x = x2m1 * xh + x3m1 * yh + x1
        !  y = y2m1 * xh + y3m1 * yh + y1
        ! where x2m1 = x2-x1; x3m1 = x3-x1; y2y1 = y2-y1; y3m1 = y3 - y1;
        ! by maple
        !   a = area
        !
        ! xy = [  x2m1*xh^2*y2m1+x2m1*xh*y3m1*yh+x3m1*xh*y2m1*yh  &
        !  +x3m1*y3m1*yh^2-x0*xh*y2m1-x0*y3m1*yh+x1*xh*y2m1+x1*y3m1*yh &
        !   -x2m1*xh*y0+x2m1*xh*y1-x3m1*y0*yh+x3m1*y1*yh     &
        !   +x0*y0-x0*y1-x1*y0+x1*y1 ]/area
        !
        ! xs = [ x2m1^2*xh^2 + 2*x2m1*x3m1*xh*yh+x3m1^2*yh^2-2*x0*x2m1*xh &
        !  -2*x0*x3m1*yh+2*x1*x2m1*xh+2*x1*x2m1*xh + 2*x1*x3m1*yh &
        ! (x0-x1)^2 ]/area
        !
        ! ys = [ y2m1^2*xh^2 + 2*y2m1*y3m1*xh*yh+y3m1^2*yh^2-2*y0*y2m1*xh &
        ! +2*y1*y2m1*xh+ 2*y0*y3m1*yh + 2*y1*y3m1*yh +(y0-y1)^2 ]/area

        x0 = element(iele)%bary(1)
        y0 = element(iele)%bary(2)

        if(n_moment>=1)then
            element(iele)%orth_d(1) = element(iele)%area
        endif


        !int_xy_t
        if(n_moment >=3)then
            xy = 0.; xs = 0.; ys = 0.
            do ig = 1, ng_tri
                xh = sx(ig)
                yh = sy(ig)

                xmx0_map = x2m1*xh+x3m1*yh+x1-x0
                ymy0_map = y2m1*xh+y3m1*yh+y1-y0

                xy = xy + xmx0_map * ymy0_map*w(ig)/element(iele)%area

                xs = xs + xmx0_map**2 *w(ig)/element(iele)%area

                ys = ys + ymy0_map**2 * w(ig)/element(iele)%area

            enddo

            int_xy = xy * 2.*element(iele)%area
            int_xs = xs * 2.*element(iele)%area
            int_ys = ys * 2.*element(iele)%area


            element(iele)%a(2,1) = -int_xy/int_xs


            element(iele)%orth_d(2) = int_xs
            element(iele)%orth_d(3) = element(iele)%a(2,1)**2*int_xs &
                + 2*element(iele)%a(2,1)*int_xy + int_ys
        endif
        !*************************************************************************
        !*************************************************************************
        if(n_moment>=6)then
            ! solve(A,b,xx,N)
            xc = 0.
            xsv3 = 0.

            a21 = element(iele)%a(2,1)
            do ig = 1, ng_tri
                xh = sx(ig)
                yh = sy(ig)

                xmx0_map = x2m1*xh+x3m1*yh+x1-x0
                ymy0_map = y2m1*xh+y3m1*yh+y1-y0

                xc = xc + xmx0_map**3 * w(ig)/element(iele)%area**(1.5)

                xsv3=xsv3+xmx0_map**2 *( a21*xmx0_map+ymy0_map )* w(ig) &
                    /element(iele)%area**(1.5)

            enddo
            int_xc = xc*2.*element(iele)%area
            int_xsv3 = xsv3*2.*element(iele)%area

            aa(1,1) = 0.
            aa(1,2) = 0.
            aa(1,3) = element(iele)%area
            aa(2,1) = int_xs
            aa(2,2) = int_xy
            aa(2,3) = 0.
            aa(3,1) = a21*int_xs+int_xy
            aa(3,2) = a21*int_xy+int_ys
            aa(3,3) = 0.

            b(1) = -int_xs
            b(2) = -int_xc
            b(3) = -int_xsv3

            call solve(aa(1:3,1:3),b(1:3),xx(1:3),3)
            element(iele)%a(3,1:3) = xx(1:3)


            !v4v4 = v4
            v4v4 = 0.
            a31 = element(iele)%a(3,1)
            a32 = element(iele)%a(3,2)
            a33 = element(iele)%a(3,3)
            ar = element(iele)%area
            do ig = 1,ng_tri

                xh = sx(ig)
                yh = sy(ig)

                xmx0_map = x2m1*xh+x3m1*yh+x1-x0
                ymy0_map = y2m1*xh+y3m1*yh+y1-y0

                v4v4=v4v4+v4v4_map(xmx0_map,ymy0_map,a31,a32,a33,ar)*w(ig)

            enddo
            element(iele)%orth_d(4) = v4v4 * 2. * element(iele)%area
            !*****************************************************************
            !******************************************************************

            xsy = 0.
            xcy = 0.
            xys = 0.
            xf = 0.
            do ig = 1, ng_tri
                xh = sx(ig)
                yh = sy(ig)

                xmx0_map = x2m1*xh+x3m1*yh+x1-x0
                ymy0_map = y2m1*xh+y3m1*yh+y1-y0

                xsy = xsy + xmx0_map**2 * ymy0_map *w(ig) &
                    /element(iele)%area**(1.5)

                xcy = xcy + xmx0_map**3 * ymy0_map *w(ig) &
                    /element(iele)%area**2

                xys = xys + xmx0_map * ymy0_map**2 *w(ig) &
                    /element(iele)%area**(1.5)

                xf = xf + xmx0_map**4  *w(ig) &
                    /element(iele)%area**2
            enddo
            int_xsy = xsy * 2. * element(iele)%area
            int_xcy = xcy * 2. * element(iele)%area
            int_xys = xys * 2. * element(iele)%area
            int_xf = xf * 2. * element(iele)%area


            ! assemble matrix
            aa(1,1)=int_xs;aa(1,2)=0.;aa(1,3)=0.;aa(1,4)=element(iele)%area;
            aa(2,1) = int_xc;  aa(2,2) = int_xs; aa(2,3) = int_xy;
            aa(2,4) = 0.;
            aa(3,1) = a21*int_xc+int_xsy;
            aa(3,2) = a21*int_xs+int_xy;
            aa(3,3) = a21*int_xy+int_ys;
            aa(3,4) = 0.;
            aa(4,1) = int_xf+a31*int_xc+a32*int_xsy+a33*int_xs;
            aa(4,2) = int_xc+a31*int_xs+a32*int_xy;
            aa(4,3) = int_xsy+a31*int_xy+a32*int_ys;
            aa(4,4) = int_xs+element(iele)%area;

            b(1) = - int_xy
            b(2) = - int_xsy
            b(3) = - (a21*int_xsy + int_xys )
            b(4) = - ( int_xcy+a31*int_xsy+a32*int_xys+a33*int_xy )

            call solve(aa(1:4,1:4),b(1:4),xx(1:4),4)
            element(iele)%a(4,1:4) = xx(1:4)


            v5v5 = 0.
            a41 = element(iele)%a(4,1)
            a42 = element(iele)%a(4,2)
            a43 = element(iele)%a(4,3)
            a44 = element(iele)%a(4,4)
            !ar = element(iele)%area
            do ig = 1,ng_tri

                xh = sx(ig)
                yh = sy(ig)

                xmx0_map = x2m1*xh+x3m1*yh+x1-x0
                ymy0_map = y2m1*xh+y3m1*yh+y1-y0

                v5v5=v5v5+v5v5_map(xmx0_map,ymy0_map,a41,a42,a43,a44,ar)*w(ig)
            enddo
            element(iele)%orth_d(5) = v5v5 * 2. * element(iele)%area
            !******************************************************************
            !******************************************************************

            yc = 0.
            xsys = 0.
            xyc = 0.
            do ig =1,ng_tri

                xh = sx(ig)
                yh = sy(ig)
                xmx0_map = x2m1*xh+x3m1*yh+x1-x0
                ymy0_map = y2m1*xh+y3m1*yh+y1-y0

                yc = yc + ymy0_map**3*w(ig)/element(iele)%area**(1.5)

                xsys = xsys +xmx0_map**2*ymy0_map**2 *w(ig) &
                    /element(iele)%area**2

                xyc = xyc + xmx0_map *ymy0_map**3 *w(ig) &
                    /element(iele)%area**2
            enddo
            int_yc = yc * 2. *  element(iele)%area
            int_xsys = xsys * 2. *  element(iele)%area
            int_xyc = xyc * 2. *  element(iele)%area

            aa(1,1) = int_xs;
            aa(1,2) = int_xy
            aa(1,3) = 0.
            aa(1,4) = 0.
            aa(1,5) = element(iele)%area

            aa(2,1) = int_xc;
            aa(2,2) = int_xsy;
            aa(2,3) = int_xs;
            aa(2,4) = int_xy;
            aa(2,5) = 0.;

            aa(3,1) = a21*int_xc+int_xsy;
            aa(3,2) = a21*int_xsy+int_xys;
            aa(3,3) = a21*int_xs+int_xy;
            aa(3,4) = a21*int_xy+int_ys;
            aa(3,5) = 0.

            aa(4,1) = int_xf+a31*int_xc+a32*int_xsy+a33*int_xs  ;
            aa(4,2) = int_xcy+a31*int_xsy+a32*int_xys+a33*int_xy;
            aa(4,3) = int_xc+a31*int_xs+a32*int_xy;
            aa(4,4) = int_xsy+a31*int_xy+a32*int_ys;
            aa(4,5) = int_xs+a33*element(iele)%area

            aa(5,1) = a41*int_xf+int_xcy+a42*int_xc+a43*int_xsy+a44*int_xs;
            aa(5,2) = a41*int_xcy+int_xsys+a42*int_xsy+a43*int_xys+a44*int_xy;
            aa(5,3) = a41*int_xc+int_xsy+a42*int_xs+a43*int_xy;
            aa(5,4) = a41*int_xsy+int_xys+a42*int_xy+a43*int_ys;
            aa(5,5) = a41*int_xs+int_xy+a44*element(iele)%area

            b(1) = -int_ys
            b(2) = -int_xys
            b(3) = -(a21*int_xys+int_yc)
            b(4) = -(int_xsys+a31*int_xys+a32*int_yc+a33*int_ys)
            b(5) = -( a41*int_xsys+int_xyc+a42*int_xys+a43*int_yc+a44*int_ys )

            call solve(aa(1:5,1:5),b(1:5),xx(1:5),5)
            element(iele)%a(5,1:5) = xx(1:5)


            v6v6 = 0.
            a51 = element(iele)%a(5,1)
            a52 = element(iele)%a(5,2)
            a53 = element(iele)%a(5,3)
            a54 = element(iele)%a(5,4)
            a55 = element(iele)%a(5,5)

            do ig = 1,ng_tri

                xh = sx(ig)
                yh = sy(ig)

                xmx0_map = x2m1*xh+x3m1*yh+x1-x0
                ymy0_map = y2m1*xh+y3m1*yh+y1-y0

                v6v6=v6v6+ &
                    v6v6_map(xmx0_map,ymy0_map,a51,a52,a53,a54,a55,ar)*w(ig)
            enddo
            element(iele)%orth_d(6) = v6v6 * 2. * element(iele)%area
        endif

    enddo


    end subroutine get_tri_info
    !*************************************************************************
    real function v4v4_map( xmx0_map,ymy0_map,a31,a32,a33,area )
    implicit none
    real,intent(in) :: xmx0_map,ymy0_map,a31,a32,a33,area

    real :: v

    v = xmx0_map**2/area+a31*xmx0_map/sqrt(area)+a32*ymy0_map/sqrt(area)+a33
    v4v4_map = v**2

    end function v4v4_map
    !*************************************************************************
    real function v5v5_map( xmx0_map,ymy0_map,a41,a42,a43,a44,area )
    implicit none
    real,intent(in) :: xmx0_map,ymy0_map,a41,a42,a43,a44,area

    v5v5_map = ( a41*xmx0_map**2/area+xmx0_map*ymy0_map/area  &
        +a42*xmx0_map/sqrt(area)+a43*ymy0_map/sqrt(area)+a44  )**2

    end function v5v5_map
    !*************************************************************************
    real function v6v6_map(xmx0_map,ymy0_map,a51,a52,a53,a54,a55,area)
    implicit none
    real,intent(in) :: xmx0_map,ymy0_map,a51,a52,a53,a54,a55,area

    real :: v

    v =a51*xmx0_map**2/area+a52*xmx0_map*ymy0_map/area &
        +ymy0_map**2/area+a53*xmx0_map/sqrt(area)+a54*ymy0_map/sqrt(area)+a55

    v6v6_map = v**2

    end function v6v6_map
    !*************************************************************************
    !*************************************************************************
    subroutine node_connectivity_to_elements
    implicit none
    integer :: iele
    integer :: n(1:3)
    integer :: i

    do i = 1,n_node_all
        node(i)%ne = 0
    enddo

    do iele = 1,n_element_all
        do i = 1,3
            n(i)= element(iele)%node(i)
            node(i)%ne = node(i)%ne + 1
            node(i)%ielement(  node(i)%ne )=iele
        enddo

    enddo

    end subroutine node_connectivity_to_elements
    !*************************************************************************
    subroutine get_regular_side3
    implicit none
    integer :: iside,i,j

    !*************************************************************************
    ! get sides
    ! by the following way
    !
    ! step 1. get sides along x direction.
    ! step 2. get sides along y direction.
    ! step 3. get symmetric sides.

    !n_side = nx*(ny+1) + ny*(nx+1) + nx*ny

    !allocate( side(1:n_side) )
    iside = 0
    ! step 1. nx*(ny+1)
    do j = 1,ny+1
        do i = 1,nx
            iside = iside +1
            side(iside)%node(1) = (j-1)*(nx+1) + i
            side(iside)%node(2) = (j-1)*(nx+1) + i +1
        enddo
    enddo

    ! step 2. ny*(nx+1)
    do i = 1 , nx+1
        do j = 1 , ny
            iside = iside + 1
            side(iside)%node(1) = i + (j-1)*(nx+1)
            side(iside)%node(2) = i + j*(nx+1)
        enddo
    enddo
    ! step 3. nx*ny
    do i = 1,nx
        do j = 1,ny
            iside = iside + 1
            side(iside)%node(1) = i + j*(nx+1)
            side(iside)%node(2) = i+1 +(j-1)*(nx+1)
        enddo
    enddo

    end subroutine get_regular_side3
    !*************************************************************************
    subroutine get_sides(ngs,gau_f)
    implicit none
    integer,intent(in) :: ngs
    real,intent(in) :: gau_f(1:3,1:2)

    integer :: iele,iside
    integer :: i1,i2
    integer :: n(1:3)
    integer :: ii,ns,it

    integer :: n1,n2
    real :: x1,y1,x2,y2,r
    integer :: ig

    do ii = 1, n_node_all
        node( ii )%nside = 0
    enddo
    !**************************************
    iside = 0
    do iele = 1 , n_element
        do i1 = 1,3
            if(i1==3)then
                i2 = 1
            else
                i2 = i1 + 1
            endif
            n( i1 ) = element(iele)%node(i1)
            n( i2 ) = element(iele)%node(i2)
            ! check n(i1) n(i2)
            it = 1
            do ii = 1,node( n(i1) )%nside
                if( n(i2) == node( n(i1) )%is2(ii) )then
                    !
                    it = 0

                    ns = node( n(i1) )%sside(ii)
                    side( ns )%be(2) = iele
                    side( ns )%ibe(2) = i1
                    element(iele)%side(i1)%n = ns
                    element(iele)%side(i1)%clock = -1.
                    side( ns )%mark = 2
                endif
            enddo
            if(it==1)then
                iside = iside + 1
                side(iside)%node(1) = n(i1)
                side(iside)%node(2) = n(i2)

                side( iside )%mark = 1
                side(iside)%be(1) = iele
                side( iside )%ibe(1) = i1
                element(iele)%side(i1)%n = iside
                element(iele)%side(i1)%clock = 1.
                !******************************
                node( n(i1) )%nside = node( n(i1) )%nside + 1
                ns = node( n(i1) )%nside
                node( n(i1) )%is2( ns ) = n(i2)

                node( n(i1) )%sside(ns) = iside
                !@
                node( n(i2) )%nside = node( n(i2) )%nside + 1
                ns = node( n(i2) )%nside
                node( n(i2) )%is2( ns ) = n(i1)

                node( n(i2) )%sside(ns) = iside
            endif

            !******************
        enddo

    enddo !iele


    !!!!!!
    do iele = n_element+1 , n_element_all
        do i1 = 1,3
            if(i1==3)then
                i2 = 1
            else
                i2 = i1 + 1
            endif
            n( i1 ) = element(iele)%node(i1)
            n( i2 ) = element(iele)%node(i2)

            ! check n(i1) n(i2)
            it = 1
            do ii = 1,node( n(i1) )%nside
                if( n(i2) == node( n(i1) )%is2(ii) )then
                    !
                    it = 0

                    ns = node( n(i1) )%sside(ii)
                    side( ns )%be(2) = iele
                    side( ns )%ibe(2) = i1
                    element(iele)%side(i1)%n = ns
                    element(iele)%side(i1)%clock = -1.
                    !side( ns )%mark = 2
                endif
            enddo
            if(it==1)then
                iside = iside + 1
                side(iside)%node(1) = n(i1)
                side(iside)%node(2) = n(i2)

                !side( iside )%mark = 1
                side(iside)%be(1) = iele
                side(iside)%ibe(1) = i1
                element(iele)%side(i1)%n = iside
                element(iele)%side(i1)%clock = 1.
                !******************************
                node( n(i1) )%nside = node( n(i1) )%nside + 1

                ns = node( n(i1) )%nside
                node( n(i1) )%is2( ns ) = n(i2)

                node( n(i1) )%sside(ns) = iside
                !@
                node( n(i2) )%nside = node( n(i2) )%nside + 1

                ns = node( n(i2) )%nside
                node( n(i2) )%is2( ns ) = n(i1)

                node( n(i2) )%sside(ns) = iside
            endif

            !******************
        enddo
    enddo !iele

    do ii = 1,n_side
        n1 = side(ii)%node(1)
        n2 = side(ii)%node(2)

        x1 = node(n1)%coor(1)
        y1 = node(n1)%coor(2)

        x2 = node(n2)%coor(1)
        y2 = node(n2)%coor(2)

        r = sqrt( (x2-x1)**2 + (y2-y1)**2  )

        side(ii)%v_normal(1) = (y2-y1)/r
        side(ii)%v_normal(2) = -(x2-x1)/r
        side(ii)%length = r


        do ig = 1,ngs
            side(ii)%gau(ig,1) = (x1+x2)/2.+(x2-x1)*gau_f(ig,1)
            side(ii)%gau(ig,2) = (y1+y2)/2.+(y2-y1)*gau_f(ig,1)
        enddo
    enddo

    !************************************************************
    ! determine the neighbor of elements,
    ! it is demanded in the numerical flux.
    do iele = 1,n_element
        do ii = 1,3

            if( element(iele)%side(ii)%clock >0. )then
                ns = element(iele)%side(ii)%n
                element(iele)%neighbor(ii) = side(ns)%be(2)
            else
                ns = element(iele)%side(ii)%n
                element(iele)%neighbor(ii) = side(ns)%be(1)
            endif
        enddo
    enddo

    end subroutine get_sides
    !*************************************************************************
    subroutine get_radius_inscribed_circle
    implicit none
    integer :: iele
    integer :: n1,n2,n3
    real :: x1,y1,x2,y2,x3,y3
    real :: d1,d2,d3
    real :: rtmp

    radius = 20200330.

    do iele = 1,n_element
        n1 = element( iele )%node(1)
        n2 = element( iele )%node(2)
        n3 = element( iele )%node(3)

        x1 = node(n1)%coor(1)
        y1 = node(n1)%coor(2)

        x2 = node(n2)%coor(1)
        y2 = node(n2)%coor(2)

        x3 = node(n3)%coor(1)
        y3 = node(n3)%coor(2)

        d1 = sqrt( (x1-x2)**2+(y1-y2)**2 )
        d2 = sqrt( (x2-x3)**2+(y2-y3)**2 )
        d3 = sqrt( (x3-x1)**2+(y3-y1)**2 )

        rtmp = 2.*( element(iele)%area )/(d1+d2+d3)
 
        radius = min(radius,rtmp)

    enddo



    end subroutine get_radius_inscribed_circle
    !*************************************************************************
    subroutine get_max_speed(t)
    use module_prob, only : velx,vely
    implicit none
    real,intent(in) :: t
    real :: v(1:2)
    integer :: i
    real :: gt(1:3,1:2)
    real :: vx,vy,st
    integer :: ig

    speed =0.
    do i  = 1,n_side
        v(1:2) = side(i)%v_normal(1:2)
        gt(1:ng_line,1:2) = side(i)%gau(1:ng_line,1:2)

        do ig =1,ng_line
            vx = velx(gt(ig,1:2) ,t)
            vy = vely( gt(ig,1:2) ,t )
            st = vx*v(1) + vy*v(2)
            speed=max(speed,abs(st) )
        enddo
    enddo

    end subroutine get_max_speed

    end module module_regular_mesh