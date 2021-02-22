    !*************************************************************************
    module module_regular_mesh_data 
    implicit none
    !*************************************************************************
    ! define the information of nodes
    ! follow a standard information such as Easymesh,
    ! in which
    !           nodes includes
    !    .   .                      .                       .
    !    .   .                      .                       .
    !    .   .                      .                       .
    !    87  8.175887722469600e-002 9.181926502205810e-001  0
    !    88  1.250000000000000e-001 1.000000000000000e+000  1
    !    89  0.000000000000000e+000 8.750000000000000e-001  1
    !    90  0.000000000000000e+000 1.000000000000000e+000  1
    !-------------------------------------------------------------
    !     n  x                      y                       mark
    !  mark = 1 means it is a boundary node;
    !  mark = 0 means it is an interior node
    type :: type_node
        sequence
        real :: coor(1:2)
        integer :: mark !
        !*****************above is a basic information
        integer :: ivalue
        real :: sol
        !************************************************
        ! for remapping
        integer :: ielement(1:10)
        integer :: ne
        !************************************************
        ! it get in our routine getting sides
        integer :: nside
        integer :: is2(1:10 ),sside(1:10)
        !*************************************************
    end type

    type(type_node), allocatable :: node(:)
    !*************************************************************************
    ! define the information of elements.
    ! Easymesh
    !  .     .    .    .    .    .    .
    !  .     .    .    .    .    .    .
    !  .     .    .    .    .    .    .
    !  141   82   87   84   145  136  143   229  223  225  1.257257364045290e-001 8.742416641198870e-001    -1
    !  142   82   85   88    -1  143  138   230  227  224  1.875000000000000e-001 9.685727724583220e-001    -1
    !  143   82   88   87   146  141  142   233  225  227  1.384652101742670e-001 9.405509479544970e-001    -1
    !  144   89   86   84   140  145   -1   228  231  232  3.132522564571040e-002 8.125000000000000e-001    -1
    !  145   84   87   89   147  144  141   234  231  229  5.939147693375600e-002 8.615550979438770e-001    -1
    !  146   87   88   90    -1  147  143   236  235  233  6.250000000000000e-002 9.807040758660390e-001    -1
    !  147   87   90   89    -1  145  146   237  234  235  1.927030964410030e-002 9.375000000000000e-001    -1
    !----------------------------------------------------------------------------------------------------------
    !     e   i,   j,   k,   ei,  ej,  ek,   si,  sj,  sk   xV,                    yV                       sign
    !
    !          / \
    !         / ei\
    !        /  si \
    !       k-------j
    !      / \      / \
    !     / sj\  e /sk \
    !    / ej  \  / ek  \
    !    -------i--------
    !
    !   ei = -1 means no neighbor element
    !   (xV,yV) 是三角形外接圆圆心？
    !   here we let it be the barycenter of the triangle!
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
        integer :: mark ! mark=1 means inside; mark = 0 means ghost.
        !*****************above is a basic information
        real :: area
        real :: gauss(1:9,1:2)

        real :: umodal(1:6,0:3) !

        real :: orth_d(1:6)
        real :: a(1:5,1:5)

        real :: sol(1:3,1:3)
        
        real :: jac_inv(1:2,1:2)
        !***********************************************
        ! information for periodic boundary conditions
        integer :: i_periodic
        !************************************************

    end type

    type(type_element), allocatable,target :: element(:)
    ! attribute target is used to be the target of pointer of the periodic
    ! ghost elements.

    !*************************************************************************
    ! define the information of sides.
    ! Easymesh
    !  233    88   87  143  146  0
    !  234    89   87  147  145  0
    !  235    90   87  146  147  0
    !  236    88   90  146   -1  1
    !  237    90   89  147   -1  1
    !--------------------------------
    !   s     c    d   ea   eb   mark
    !
    !            d-------
    !           / \ eb  /
    !          / side  /
    !         /  ea \ /
    !         -------c

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

    !*************************************************************************
    ! data structure for forming ghost cells.
    !
    !     ----
    !    |\ie1|
    !    | \  |
    !    |  \ |
    !    |ie2\|
    !     ----
    !      sidelower
    type :: type_side_b
        sequence
        integer :: node(1:2)
        integer :: ie(1:2)
    end type

    type(type_side_b), allocatable :: side_lower(:,:)
    type(type_side_b), allocatable :: side_upper(:,:)
    type(type_side_b), allocatable :: side_left(:,:)
    type(type_side_b), allocatable :: side_right(:,:)

    end module module_regular_mesh_data 