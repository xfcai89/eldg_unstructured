    !*************************************************************************
    subroutine term_flux(io,tt)
    use module_polynomials_tri
    implicit none
    integer,intent(in) :: io
    real,intent(in) :: tt

    integer :: iele,ii,ig
    real :: fi
    real :: xx,yy,ww,uu 
    real :: v1,v2
    integer :: ns

    real :: pc(1:2),ar,oa(1:5,1:5)
    real :: co(1:2),cc,test
    real :: ff,gg

    integer :: nm

    real :: fnum,smax

    integer :: nei
    real :: x0n,y0n,arn,oan(1:5,1:5)
    integer :: iin

    real :: f1,f2,g1,g2

    real :: pi
    integer :: n1,n2

    real :: fins,gins,fout,gout

    pi = 4.*atan(1.)
 

    ! assemble all solution at gauss point on each edge.
    ! below we organize them in the anticlock-wise.

    do iele = 1, n_element_all

        pc(1:2) = element(iele )%bary(1:2)
 
        ar = element(iele )%area
        oa(1:5,1:5) = element(iele )%a(1:5,1:5)


        !for each side of the element(iele)
        do ii = 1,3

            ns = element(iele)%side(ii)%n
            cc = element(iele)%side(ii)%clock

            !v1 = side(ns)%v_normal(1)*cc
            !v2 = side(ns)%v_normal(2)*cc


            do ig = 1,ng_line

                co(1:2) = side( ns )%gau(ig,1:2)
             
                uu = poly_tri( n_moment, &
                    element(iele)%umodal(1:n_moment,io),co(1:2),pc(1:2),ar,oa(1:5,1:5)  )


                !ff = fun_f(uu,co(1:2),tt)
                !gg = fun_g(uu,co(1:2),tt)

                !element(iele)%fg( ig,ii )= ff*v1 + gg*v2

                !element(iele)%sol( ig,ii ) = uu
                if(cc>0.)then
                    element(iele)%sol( ig,ii ) = uu
                else
                    element(iele)%sol( ng_line+1-ig,ii ) = uu
                endif

                !print *,iele,'side',ii,ig,cc
                !print *,co(1:2)
                !pause

            enddo
        enddo
    enddo !iele

    do iele = 1,n_element

        pc(1:2) = element(iele)%bary(1:2)
     
        ar = element(iele)%area
        oa(1:5,1:5) = element(iele)%a(1:5,1:5)
        !*******************
        do nm = 1,n_moment

            flux(iele,nm) = 0.

            do ii = 1,3

                ns = element(iele)%side(ii)%n

                !----------------------\
                ! for debugging
                n1 = side(ns)%node(1)
                n2 = side(ns)%node(2)
                !----------------------/

                cc =element(iele)%side(ii)%clock

                v1 = side(ns)%v_normal(1)*cc
                v2 = side(ns)%v_normal(2)*cc
                !*****************************************
                ! find which side of the neighbor element.
                ! iin side
                if( element(iele)%side(ii)%clock>0. )then
                    iin = side( ns )%ibe(2)
                else
                    iin = side( ns )%ibe(1)
                endif


                !*****************************************
                nei = element(iele)%neighbor(ii)

                !--
                !x0n = element(nei)%bary(1)
                !y0n = element(nei)%bary(2)
                !arn = element(nei)%area
                !oan(1:5,1:5) = element(nei)%a(1:5,1:5)

                !--
                !nk = 1
                do ig = 1,ng_line
                    if(cc>0.)then
                        co(1:2) = side( ns )%gau(ig,1:2)
 
                    else
                        co(1:2) = side( ns )%gau(ng_line+1-ig,1:2)
 
                    endif
 

                    ! Lax-Friedrichs flux
                    f1 = fun_fp( element(iele)%sol( ig,ii ),co(1:2),tt )
                    g1 = fun_gp( element(iele)%sol( ig,ii ),co(1:2),tt )
                    f2 = fun_fp( element(nei)%sol(ng_line+1-ig,iin),co(1:2),tt )
                    g2 = fun_gp( element(nei)%sol(ng_line+1-ig,iin),co(1:2),tt )

                    smax = max( abs(f1*v1+g1*v2),abs(f2*v1+g2*v2) )


                    fins = fun_f( element(iele)%sol( ig,ii ),co(1:2),tt )
                    gins = fun_g( element(iele)%sol( ig,ii ),co(1:2),tt )

                    fout = fun_f( element(nei)%sol(ng_line+1-ig,iin),co(1:2),tt )
                    gout = fun_g( element(nei)%sol(ng_line+1-ig,iin),co(1:2),tt )

                    fnum = 0.5*( fins*v1+gins*v2 + fout*v1+gout*v2 &
                        +smax*( element(iele)%sol( ig,ii ) - element(nei)%sol(ng_line+1-ig,iin)  ) )

                    test = ortho_basis_tri(nm,co(1:2),pc(1:2),ar ,oa(1:5,1:5))
                    flux(iele,nm) = flux(iele,nm) + fnum*test*side(ns)%length*gau_line(ig,2)

  

                enddo ! ig
            enddo !ii

        enddo !nm

    enddo  !iele

    end subroutine term_flux