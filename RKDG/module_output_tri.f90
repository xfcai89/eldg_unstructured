    module module_output_tri
    use module_regular_mesh_data, only : node,element
    implicit none

    contains
    !*************************************************************************
    subroutine output_tri_mesh(n_node,n_element)

    implicit none
    integer,intent(in) :: n_node,n_element
    integer :: i

    open(11,file='mesh.plt')
    write(11,*) 'title="Hermes data file"'
    write(11,*) 'variables="x","y","den"'
    write(11,*) 'zone n=', n_node, ',e=',n_element,',f=fepoint,et=triangle'

    do i = 1, n_node
        write(11,101) node(i)%coor(1),node(i)%coor(2),0.
    enddo

    do i = 1,n_element
        write(11,*) element(i)%node(1),element(i)%node(2),element(i)%node(3)
    enddo

101 FORMAT(1X,3F14.4)

    close(11)

    end subroutine output_tri_mesh
    !*************************************************************************
    subroutine output_tri(n_node,n_element)

    implicit none
    integer,intent(in) :: n_node,n_element
    integer :: i

    open(11,file='solution.plt')
    write(11,*) 'title="Hermes data file"'
    write(11,*) 'variables="x","y","den"'
    write(11,*) 'zone n=', n_node, ',e=',n_element,',f=fepoint,et=triangle'


    do i = 1, n_node

        write(11,101) node(i)%coor(1),node(i)%coor(2),node(i)%sol

    enddo

    do i = 1,n_element

        write(11,*) element(i)%node(1),element(i)%node(2),element(i)%node(3)
    enddo

101 FORMAT(1X,3F20.4)

    close(11)

    end subroutine output_tri
    !*************************************************************************
    subroutine solution_nodes(io,n_node,n_element)
    use module_polynomials_tri, only : poly_tri
    use module_regular_mesh, only : n_moment
    implicit none

    integer,intent(in) :: io,n_node,n_element

    integer :: iele,i
    integer :: no(1:3)
    real :: pt(1:2)
    real :: ar,oa(1:5,1:5)
    real :: pc(1:2)

    do i = 1,n_node
        node(i)%ivalue = 0
    enddo


    do iele = 1,n_element
        ! update node value
        do i = 1,3
            no(i) = element(iele)%node(i)
            if( node( no(i) )%ivalue==0 )then
                node( no(i) )%ivalue=1
                pt(1:2) = node( no(i) )%coor(1:2)

                ar = element(iele)%area
                oa(1:5,1:5) = element(iele)%a(1:5,1:5)
                pc(1:2) = element(iele)%bary(1:2)
   
                node( no(i) )%sol = poly_tri( n_moment, &
                    element(iele)%umodal(1:n_moment,io),&
                    pt(1:2),pc(1:2),ar,oa(1:5,1:5))

            elseif( node(no(i))%ivalue==1 )then
                ! do nothing.
            endif
        enddo !i
    enddo

    end subroutine solution_nodes
    !*************************************************************************

    end module module_output_tri