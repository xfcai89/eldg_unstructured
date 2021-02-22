    subroutine setdt
    implicit none

    dt = cfl*radius/speed
    dt = cfl*radius/3.14

    if(time+dt>tfinal) dt = tfinal- time


    !time = time +dt

    print *,time,dt
   

    end subroutine setdt