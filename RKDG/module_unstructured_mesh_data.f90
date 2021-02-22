    !*************************************************************************
    module module_unstructured_mesh_data
    implicit none
    !*************************************************************************
    ! define the information of nodes

    type :: type_node
        sequence
        real(kind=8) :: coor(1:2)
        !*****************above is a basic information
        integer :: ivalue
        real :: sol
        !************************************************
        ! it get in our routine getting sides
        integer :: nside
        integer :: is2(1:10 ),sside(1:10)
        !*************************************************
    end type

    type(type_node), allocatable :: node(:)
    !*************************************************************************
    ! define the information of elements.
    type :: type_side_of_element
        sequence
        integer :: n
        real :: clock ! clock=1. means anticlockwise
        ! clock =-1. means clockwise
    end type
    type :: type_element
        sequence
        integer :: node(1:3)
        integer :: neighbor(1:3)
        !integer :: side(1:3)
        type(type_side_of_element) :: side(1:3)
        real :: bary(1:2)
        !*****************above is a basic information
        real :: area
        real :: gauss(1:9,1:2)
        
        real :: umodal(1:6,0:3) !        
        
        real :: a(1:5,1:5)
        
        real :: orth_d(1:6)
        
        real :: sol(1:3,1:3)
        
        real :: jac_inv(1:2,1:2)
        
        !************************************************
        integer :: icut ! this is will be used in module_upstream
    end type

    type(type_element), allocatable,target :: element(:)
    
    
    type :: type_side
        sequence
        integer :: node(1:2)
        integer :: be(1:2)
        integer :: mark ! mark=1 is boundary; mark = 2 is interior side.

        integer :: ibe(1:2) ! cooresponding side of the be(1:2)
        !*****************above is a basic information
        ! normal vector
        real :: v_normal(1:2)
        real :: gau(1:3,1:2)
        real :: length


    end type

    type(type_side), allocatable :: side(:)  
    type(type_side), allocatable :: side_tmp(:) 
    
    end module module_unstructured_mesh_data