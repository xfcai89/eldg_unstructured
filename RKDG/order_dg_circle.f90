    !*************************************************************************
    subroutine order_dg_circle
    use module_polynomials_tri, only : polyt

    implicit none
    integer :: iele,ig
    real :: error1,error2,errori,rr1,rr2,rri

    real :: s
    real :: pg(1:2),pc(1:2),ar,oa(1:5,1:5)

    real :: etmp
    real :: umod(1:n_moment)
    
    real,parameter :: pi = 4.*atan(1.)
    real :: tmp


    errori = 0.

    error1 = 0.
    error2 = 0.
    do iele = 1,n_element
        pc(1:2) = element(iele)%bary(1:2)
        ar = element(iele)%area
        oa(1:5,1:5) = element(iele)%a(1:5,1:5)
        umod(1:n_moment)=element(iele)%umodal(1:n_moment,0)
        do ig = 1, ng_tri

            pg(1:2) = element(iele)%gauss(ig,1:2)

            s=polyt(n_moment,umod(1:n_moment),pg(1:2),pc(1:2),ar,oa(1:5,1:5))


            errori = max(errori,  s-fun_exact( pg(1:2),tfinal ) )
 

            !************************
            etmp =  abs(s-fun_exact( pg(1:2),tfinal ) )
            error1 = error1 +etmp*w(ig)*2.*ar
            error2 = error2 + etmp**2*w(ig)*2.*ar
        enddo


    enddo
    error1=error1/( pi**3 )
    error2=sqrt(error2/( pi**3  ))
    !********************************
    ! error distribution

    !********************************

    if(kkkk.eq.1) write(123,103) n_element,error1,error2,errori
    write(*,*) error1,error2,errori,'n',n_element

    if(kkkk.gt.1) then
        tmp = sqrt(  real(n_element)/real(norder(kkkk-1)  )     )
        rr1=log(er1/error1)/log( tmp )
        rr2=log(er2/error2)/log( tmp )
        !rri=log(eri/errori)/log( real(kkkk)/real(kkkk-1) )
        rri=log(eri/errori)/log( tmp )
        write(123,102) n_element,error1,rr1,error2, rr2,errori, rri
        write(*,*) n_element,rr1,rr2,rri
    endif

    close(111)


111 format(4(1x,f12.4))
102 format( i6,1x,3('&',1x, es12.2e2,1x,'&' 1x,f8.2 ,1x),'&',i8,'&',i8,'\\',1x,'\hline')
103 format( i6,1x,3('&',1x,es12.2E2,1x,'&',1x) ,'&',i8,'&',i8,'\\',1x,'\hline')

    eri = errori
    er1 = error1
    er2 = error2

    !call solution_nodes(0)
    !call output_tri(n_node,n_element)

    end subroutine order_dg_circle
    !*************************************************************************