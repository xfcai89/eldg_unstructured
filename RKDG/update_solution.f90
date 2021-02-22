    subroutine update_solution
    implicit none
    
    integer :: iele,nm

    do iele = 1,n_element
        do nm = 1,n_moment
            element(iele)%umodal(nm,0) = element(iele)%umodal(nm,mrk)
        enddo
    enddo

    end subroutine update_solution