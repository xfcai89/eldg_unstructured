    !*************************************************************************
    subroutine boundary_peri(io)
    implicit none

    integer,intent(in) :: io

    integer :: iele,ip,i,no(3)

    do iele = n_element+1,n_element_all
        ip = element(iele)%i_periodic

        element(iele)%umodal(1:n_moment,io)=element(ip)%umodal(1:n_moment,io)
        ! update node value
        ! the comment part is for testing the boundary conditions
        !do i = 1,3
        !    no(i) = element(iele)%node(i)
        !    if( node( no(i) )%ivalue==0 )then
        !        node( no(i) )%ivalue=1
        !        node( no(i) )%sol = element(iele)%umodal(1)
        !    elseif( node(no(i))%ivalue==1 )then
        !        ! do nothing.
        !    endif
        !enddo
    enddo


    !call output_tri(n_node_all,n_element_all)

    end subroutine boundary_peri