    !*************************************************************************
    subroutine runge_kutta(io)
    implicit none
    integer,intent(in) :: io

    integer :: iele,m

    type(type_element),pointer :: p

    if( mrk==1 )then

        do iele = 1,n_element
            do m = 1 , n_moment
                p =>element(iele)
                p%umodal(m,1) = p%umodal(m,0)+dt*(volume(iele,m)-flux(iele,m))/p%orth_d(m)

            enddo
        enddo

    elseif( mrk==2 )then

        if( io==0 )then
            do iele = 1,n_element
                do m = 1 , n_moment
                    p =>element(iele)
                    p%umodal(m,1) = p%umodal(m,0)+dt*(volume(iele,m)-flux(iele,m))/p%orth_d(m)

                enddo
            enddo
        elseif( io==1 )then
            do iele = 1,n_element
                do m = 1 , n_moment
                    p =>element(iele)
                    p%umodal(m,2) = 0.5*p%umodal(m,0)&
                    +0.5*(p%umodal(m,1)+dt*(volume(iele,m)-flux(iele,m))/p%orth_d(m) )

                enddo
            enddo
        endif

    elseif(mrk==3)then

        if( io==0 )then
            do iele = 1,n_element
                do m = 1 , n_moment
                    p =>element(iele)
                    p%umodal(m,1) = p%umodal(m,0)+dt*(volume(iele,m)-flux(iele,m))/p%orth_d(m)

                enddo
            enddo
        elseif( io==1 )then
            do iele = 1,n_element
                do m = 1 , n_moment
                    p =>element(iele)
                    p%umodal(m,2) = 0.75*p%umodal(m,0)&
                    +0.25*(p%umodal(m,1)+dt*(volume(iele,m)-flux(iele,m))/p%orth_d(m) )

                enddo
            enddo
        elseif(io==2)then
            do iele = 1,n_element
                do m = 1 , n_moment
                    p =>element(iele)
                    p%umodal(m,3) = 1./3.*p%umodal(m,0)&
                    +2./3.*(p%umodal(m,2)+dt*(volume(iele,m)-flux(iele,m))/p%orth_d(m) )

                enddo
            enddo            
        endif
    endif

    end subroutine runge_kutta