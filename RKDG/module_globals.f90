    module module_globals

    integer  :: kkkk
    real :: eri,er1,er2
 
    real :: time,  dt
    integer :: nt
    real,allocatable :: flux(:,:) ! flux term
    real,allocatable :: volume(:,:)

    integer :: io

    ! for error distribution
    real,allocatable :: error_distri(:)
    end module module_globals