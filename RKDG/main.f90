    !*************************************************************************
    !
    !  PROGRAM: main
    !
    !  PURPOSE:
    !
    !            By Xiaofeng Cai, University of Delaware, 04/18/2020
    !*************************************************************************

    program main
    use module_globals
    use module_unstructured_mesh_data
    use module_unstructured_mesh

    use module_gauss_tri
    use module_read_gmsh

    use module_prob
    !*************************************************************************
    use module_output_tri
    implicit none
    character ( len = 255 ) gmsh_filename

    ! Body of main program

    call setup_prob

    do kkkk = 1,5
        ! read the mesh from a gmsh file
        ! real the number of nodes, n_node
        ! real the number of element, n_element
        !**************************
        if(kkkk==1)then
            gmsh_filename = 'circle_ele160_h0d8'
        elseif( kkkk==2 )then
            gmsh_filename = 'circle_ele522_h0d4'
        elseif(kkkk==3)then
            gmsh_filename = 'circle_ele1884_h0d2'
        elseif( kkkk==4 )then
            gmsh_filename = 'circle_ele7432_h0d1'
        elseif( kkkk==5 )then
            gmsh_filename = 'circle_ele28996_h0d05'
        endif

        call gmsh_nodes_elements(gmsh_filename,n_node,n_element)
        !**************************
        norder(kkkk)=n_element
        !****************************
        !gaussian quadrature
        call allocate_gau_para
        call quadrature_rules_tri( ng_tri,w(1:ng_tri),sx(1:ng_tri),sy(1:ng_tri) )
        !****************************

        call allocate_dg

        call get_sides_no_ghost
        call get_sides_vector_gauss(ng_line,gau_line(1:3,1:2 ))

        call get_tri_info(n_element)

        call get_radius_inscribed_circle

        ! intial condition
        call init

        nt = 0
        time = 0.

        do while( time<tfinal )
            call get_max_speed(time)

            call setdt
            do io = 0,mrk-1
                trk_stage  = time  + rk_d(io)*dt

                call term_flux_no_bc(io,trk_stage )
                call term_volume(io,trk_stage  )

                call runge_kutta(io)

            enddo !io

            ! update solution and time
            call update_solution
            time = time + dt
        enddo
        call solution_nodes(0,n_node,n_element)
        call output_tri(n_node,n_element)

        call order_dg_circle

        call deallocate_dg

        ! deallocating the gmsh
        !**************************
        call deallocate_gmsh
        !**************************

        !****************************
        !gaussian quadrature
        call deallocate_gau_para
        !****************************

    enddo

    contains

    include "allocate_dg.f90"

    include "init.f90"
    include "setdt.f90"

    include "term_flux_no_bc.f90"
    include "term_volume.f90"
    include "runge_kutta.f90"

    include "update_solution.f90"

    include "order_dg_circle.f90"

    end program main