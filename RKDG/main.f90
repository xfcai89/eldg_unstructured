    !*************************************************************************
    ! Runge-Kutta Disscontinuous Galerkin method on unstructured mesh
    ! for solving
    !   u_t + f(u)_x + g(u)_y =0.
    !                     Authors:
    !                                  Xiaofeng Cai
    !                              University of Delaware
    !                                March 29th, 2020
    !*************************************************************************
    program main
    use module_gauss_tri    ! please set up the gaussian quadrature.
    use module_regular_mesh

    use module_globals
    !*************************************************************************
    use module_prob
    implicit none


    call setup_prob

    !****************************
    !gaussian quadrature
    call allocate_gau_para
    call quadrature_rules_tri
    !*****************************

    do kkkk = 1,5
        !generating a periodic mesh including ghost cells.
        call setup_trimesh(kkkk )
        call allocate_trimesh

        call allocate_dg

        allocate( error_distri(1:n_element) )

        call get_nodes_elements
        call get_sides(ng_line,gau_line(1:3,1:2 ))


        call get_tri_info(n_element_all)

        call get_radius_inscribed_circle



        ! intial condition
        call init
        ! periodic boundary conditions in both x and y directions.
        !call solution_nodes(0)
        !call output_tri(n_node,n_element)
        !pause

        nt = 0
        time = 0.

        do while(time<tfinal  )
            call get_max_speed(time)

            call setdt
            do io = 0,mrk-1

                call boundary_peri(io)
                trk_stage  = time  + rk_d(io)*dt
                call term_flux(io,trk_stage )
                call term_volume(io,trk_stage  )

                call runge_kutta(io)

            enddo

            ! update solution and time
            call update_solution
            time = time + dt

            !call order_dg

            !pause
            !goto 911
        enddo
        !911 continue
        call order_dg

        !call solution_nodes(0)
        !call output_tri(n_node,n_element)

        call deallocate_trimesh
        call deallocate_dg


        deallocate( error_distri )
    enddo


    !****************************
    !gaussian quadrature
    call deallocate_gau_para
    !*****************************
    contains


    include "init.f90"
    include "boundary_peri.f90"
    include "order_dg.f90"

    include "setdt.f90"

    include "term_flux.f90"
    include "term_volume.f90"

    include "allocate_dg.f90"


    include "runge_kutta.f90"

    include "update_solution.f90"

    end program main