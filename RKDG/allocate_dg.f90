    subroutine allocate_dg
    implicit none
    
    !*******************************************
    ! for RKDG flux and volume 
    allocate( flux(1:n_element,1:n_moment) )
    allocate( volume(1:n_element,1:n_moment) )
    ! end for RKDG flux and volume
    !*******************************************

    end subroutine allocate_dg
    !*************************************************************************
    !*************************************************************************
    subroutine deallocate_dg
    implicit none

    !*******************************************
    ! for RKDG flux and volume 
    deallocate( flux )
    deallocate( volume )
    ! end for RKDG flux and volume
    !*******************************************

    end subroutine deallocate_dg