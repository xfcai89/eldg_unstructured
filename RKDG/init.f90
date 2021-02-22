    !*************************************************************************
    subroutine init
    use module_polynomials_tri, only : ortho_basis_tri

    !use module_output_tri
    !use module_trouble_limiter
    implicit none
    ! ngauss w are from module_gauss_tri
    ! local variables
    integer :: iele,k,ig
    real :: utmp,pg(1:2),pc(1:2),area,oa(1:5,1:5)
    real :: test


    do iele = 1 , n_element
        !****************************************
        ! for calculating the polynomial
        pc(1:2) = element(iele)%bary(1:2)
        area = element(iele)%area
        oa(1:5,1:5) = element(iele)%a(1:5,1:5)
        !****************************************
        do k = 1 , n_moment
            utmp = 0.
            do ig = 1, ng_tri
                pg(1:2) = element(iele)%gauss(ig,1:2)

                test = ortho_basis_tri(k,pg(1:2),pc(1:2),area,oa(1:5,1:5))
                utmp = utmp + fun_init( pg(1:2) )*test*w(ig)*2.*area
            enddo
            element(iele)%umodal(k,0) = utmp/element(iele)%orth_d(k)

        enddo
    enddo

    !call solution_nodes(0,n_node,n_element)
    !call output_tri(n_node,n_element)
    !
    !call trouble(0)
    !call wenolimiter(0)
    !
    !call solution_nodes(0,n_node,n_element)
    !call output_tri(n_node,n_element)


    end subroutine init