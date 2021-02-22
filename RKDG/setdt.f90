    subroutine setdt
    implicit none

    dt = cfl*radius/speed
   

    if(time+dt>tfinal) dt = tfinal- time


    !time = time +dt

    print *,time,dt
   

    end subroutine setdt