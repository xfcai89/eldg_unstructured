    !*************************************************************************
    !
    !  MODULE: read function for Gmsh
    !
    !  PURPOSE:  for read the Gmsh in the MSH file format version 1 
    !
    !                    By Xiaofeng Cai, University of Delaware, 03/26/2020
    !*************************************************************************
    !  $NOD
    !  number-of-nodes
    !  node-number x-coord y-coord z-coord
    !  ¡­
    !  $ENDNOD
    !  $ELM
    !  number-of-elements
    !  elm-number elm-type reg-phys reg-elem number-of-nodes node-number-list
    !  ¡­
    !  $ENDELM
    !*************************************************************************
    !  List of Routines:
    !     
    !        gmsh_nodes_elements is for getting the nodes and elements from a 
    !        gmsh mesh file in the MSH file format version 1. 
    !*************************************************************************
    module module_read_gmsh


    implicit none
    
    contains
    !*************************************************************************
    ! Input: 
    !         the file name: gmsh_filename.
    ! Output: 
    !         the number of nodes: n_node.
    !         the number of elements: n_element.
    !         nodes of the unstructured mesh: node.
    !         elements of the unstructured mesh: element.
    !*************************************************************************
    subroutine gmsh_nodes_elements(gmsh_filename,n_node,n_element)
    use module_unstructured_mesh_data, only : node,element
    implicit none
    character(len=255),intent(in) :: gmsh_filename
    integer,intent(out) :: n_node,n_element
    
    ! local variables
    integer :: input_stat,input
    character(len=25) text
    real(kind=8) :: useless
    integer :: i,ii

    integer :: n1,n2,n3,n4,n5
    integer :: n6,n7,n8,n_all

    integer :: istart,iele
    
    input = 11
    open(unit=input,file=gmsh_filename,status='old',iostat = input_stat )
    !¡°old¡± indicates that the file is already there
    if ( input_stat /= 0 ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'GMSH_DATA_READ - Fatal error!'
        write ( *, '(a)' ) 'Could not open input file'
    endif

    read ( input,*) text
 

    read ( input,* ) n_node
 
    allocate( node(1:n_node) )

    do i = 1,n_node
        read(input, *) &
            ii,node(i)%coor(1),node(i)%coor(2),useless

    enddo


    read ( input,*) text
 
    read ( input,*) text
 
    read ( input,* ) n_all
 

    istart =0
    iele = 0

    do i = 1,n_all
        read(input,*) n1,n2,n3,n4,n5
        if(n5==1)then
            BACKSPACE (input)
            read(input, *) n1,n2,n3,n4,n5,n6

        elseif(n5==2)then
            BACKSPACE (input)
            read(input, *) n1,n2,n3,n4,n5,n6,n7

        elseif(n5==3)then

            if( istart == 0 )then
                istart = 1
                n_element = n_all - (n1-1)
                allocate( element(1:n_element) )
            endif
            BACKSPACE (input)
            read(input, *) n1,n2,n3,n4,n5,n6,n7,n8

            iele = iele + 1

            element(iele)%node(1) = n6
            element(iele)%node(2) = n7
            element(iele)%node(3) = n8

        endif

    enddo

    close(input)
    
    end subroutine gmsh_nodes_elements 
    !*************************************************************************
    subroutine deallocate_gmsh
    use module_unstructured_mesh_data, only : node,element,side
    implicit none

    deallocate( node )
    deallocate( element )
    
    deallocate( side )
    
    end subroutine deallocate_gmsh
    !*************************************************************************

    end module module_read_gmsh