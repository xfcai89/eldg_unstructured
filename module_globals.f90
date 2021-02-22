    module module_globals

    real :: time, dt
    integer :: nt

    integer :: io
    
    real,allocatable :: flux(:,:) ! flux term
    real,allocatable :: volume(:,:)

    integer  :: kkkk,norder(1:6)
    real :: eri,er1,er2
    
    end module module_globals