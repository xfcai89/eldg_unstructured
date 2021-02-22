    !*************************************************************************
    subroutine order_dg
    use module_polynomials_tri, only : poly_tri

    implicit none
    integer :: iele,ig
    real :: error1,error2,errori,rr1,rr2,rri

    real :: sol
    real :: pg(1:2),pc(1:2),ar,oa(1:5,1:5)
    
    real :: etmp


    errori = 0.
    
    error1 = 0.
    error2 = 0.
    do iele = 1,n_element
        
        do ig = 1, ng_tri
 
            pg(1:2) = element(iele)%gauss(ig,1:2)
            
            pc(1:2) = element(iele)%bary(1:2)
            ar = element(iele)%area
            oa(1:5,1:5) = element(iele)%a(1:5,1:5)
            sol = poly_tri( n_moment,element(iele)%umodal(1:n_moment,0),pg(1:2),pc(1:2),ar,oa(1:5,1:5) )

 
            errori = max(errori,  sol-fun_exact( pg(1:2),tfinal ) )

            error_distri(iele) = sol-fun_exact( pg(1:2),tfinal )
            
            !************************
            etmp =  abs(sol-fun_exact( pg(1:2),tfinal ) )
            error1 = error1 +etmp*w(ig)*2.*ar 
            error2 = error2 + etmp**2*w(ig)*2.*ar 
        enddo
        
        
    enddo
    error1=error1/( (xright-xleft)*(ytop-ybottom) )
    error2=sqrt(error2/( (xright-xleft)*(ytop-ybottom) ))
    !********************************
    ! error distribution

    !********************************

    if(kkkk.eq.1) write(123,103) nx,ny,error1,error2,errori
    write(*,*) error1,error2,errori,'nx',nx

    if(kkkk.gt.1) then
        rr1=log(er1/error1)/log( 2. )
        rr2=log(er2/error2)/log( 2. )
        !rri=log(eri/errori)/log( real(kkkk)/real(kkkk-1) )
        rri=log(eri/errori)/log( 2. )
        write(123,102) nx,ny,error1,rr1,error2, rr2,errori, rri
        write(*,*) nx,ny,rr1,rr2,rri
    endif

    close(111)


111 format(4(1x,f12.4))
102 format(i6,'*',i6,1x,3('&',1x, es12.2e2,1x,'&' 1x,f8.2 ,1x),'&',i8,'&',i8,'\\',1x,'\hline')
103 format(i6,'*',i6,1x,3('&',1x,es12.2E2,1x,'&',1x) ,'&',i8,'&',i8,'\\',1x,'\hline')

    eri = errori
    er1 = error1
    er2 = error2

    !call solution_nodes(0)
    !call output_tri(n_node,n_element)

    end subroutine order_dg
    !*************************************************************************
 