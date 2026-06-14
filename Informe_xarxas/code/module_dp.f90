module dp
    implicit none

    public :: exponential_distribution, append_int,rem_element,create_degree_list,create_pointers, create_neighbors

contains

SUBROUTINE append_int(vec, val, val2)
    !***********************************************************************
    ! Appends val in vec if not already present
    ! Evaluate if it is a 2D or 1D array
    !***********************************************************************
    INTEGER, DIMENSION(..), ALLOCATABLE, INTENT(INOUT)   :: vec
    INTEGER, DIMENSION(:), ALLOCATABLE                 :: vec_aux1
    INTEGER, DIMENSION(:,:), ALLOCATABLE                 :: vec_aux
    INTEGER, INTENT(IN)                                  :: val
    INTEGER, OPTIONAL, INTENT(IN)                        :: val2
    
    select rank(vec)
        rank(1)
            IF (.NOT. ANY(vec .EQ. val)) THEN
                ALLOCATE(vec_aux1(size(vec, 1)+1))
                vec_aux1(:size(vec)) = vec;  vec_aux1(size(vec)+1) = val
                call move_alloc(vec_aux1,vec)
                
            END IF
        rank(2)
            ALLOCATE(vec_aux(size(vec, 1)+1,size(vec, 2)))
            vec_aux(:size(vec, 1), :)=vec; vec_aux(size(vec, 1)+1,:)=(/val,val2/)
            call move_alloc(vec_aux,vec)
    end select
    END SUBROUTINE

    subroutine rem_element(vec, val, swicht_case)

    !***********************************************************************
    ! Remove a value in a vec
    ! Evaluate if it is a 2D or 1D array and can be choosen if it remove a value 
    ! or a index
    ! It also remove multiple values if a 2d array is given as val
    !***********************************************************************
        INTEGER, DIMENSION(..), ALLOCATABLE, INTENT(INOUT)   :: vec
        INTEGER, DIMENSION(:), ALLOCATABLE                   :: vec_aux1
        INTEGER, DIMENSION(:,:), ALLOCATABLE                 :: vec_aux
        INTEGER,DIMENSION(..), INTENT(IN)                    :: val
        logical, dimension(size(vec,1))                      :: mask
        CHARACTER(*), INTENT(IN)                             :: swicht_case
        integer                                              :: i
        
        select rank(vec) ! Case 1D vector
            rank(1)
                select rank(val)
                    rank(0)
                    allocate(vec_aux1(size(vec)-1))
                    select case(swicht_case)
                        case('index')
                            vec_aux1= pack(vec, mask=((/(i, i=1,size(vec))/).ne.val))
                        case('value')
                            vec_aux1= pack(vec, mask=(vec.ne.val))
                    end select
                    call move_alloc(vec_aux1, vec)
                endselect
            
            rank(2) ! Case 2D vector
                select rank(val)
                    rank(0)

                    allocate(vec_aux(size(vec,1)-1,size(vec,2)))
                    select case(swicht_case)
                        case('index')
                            vec_aux(:,1)=pack(vec(:,1), mask=((/(i, i=1,size(vec,1))/).ne.val))
                            vec_aux(:,2)=pack(vec(:,2), mask=((/(i, i=1,size(vec,1))/).ne.val))

                        case('value')
                            vec_aux(:,1)=pack(vec(:,1), mask=(vec(:,1).ne.val .and.vec(:,2).ne.val))
                            vec_aux(:,2)=pack(vec(:,2), mask=(vec(:,1).ne.val .and.vec(:,2).ne.val))
                    end select
                    call move_alloc(vec_aux, vec)

                    rank(1)
                    allocate(vec_aux(size(vec,1)-size(val),size(vec,2)))
                    select case(swicht_case)
                        case('index')
                            mask=.True.; mask(val)=.False.
                            vec_aux(:,1)=pack(vec(:,1), mask=mask)
                            vec_aux(:,2)=pack(vec(:,2), mask=mask)

                        case('value')
                            mask=.True.
                            do i=1, size(val)
                            where (vec(:,1).eq. val(i) .or. vec(:,2).eq. val(i)) mask=.False.
                            enddo
                            vec_aux(:,1)=pack(vec(:,1), mask=mask)
                            vec_aux(:,2)=pack(vec(:,2), mask=mask)
                    end select
                    call move_alloc(vec_aux, vec)

                endselect
        end select

    endsubroutine



    function exponential_distribution(lambda)result(exp_dis)
        !Function to create a exponentian distribution using a uniform distribution.
        !The expression is -1/lambda*ln(1-R) where R is a random variable of a uniform
        !distribution between [0,1)
        !Input:
        !   lambda    real. Parameter in the exponential distribution -> distribution lambda exp(-lambda x)
        !Output:
        !   exp_dis   real. Number distributed in a exponential distribution.
        use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
        implicit none
    
        !Input variable
        real(kind=8), intent(in)                        :: lambda
        !Output 
        real(kind=8)                                    :: exp_dis
        !Internal
        real(kind=8)                                   :: ran
    
        call random_number(ran)
        exp_dis = -1.d0/lambda*dlog((1.d0-mod(ran,1.d0)))
    
    end function exponential_distribution
    
    !THE FOLLOWING SUBROUTINES ARE ALREADY IN THE FILE SUBROUTINE.F90
    !THEY ARE HERE TOO BECAUSE IT IS DIFFICULT TO ACCES TO A FORTRAN
    !FILE WITHOUT A MODULE

    subroutine create_degree_list(N_nodes,edge_list, degree_list)
        ! Creates a degree list where the index correspondes to the label
        ! of the node and the value the number of neighbors of this node.
        !-----------------------------------------------------------------------
        ! INPUT
        !-----------------------------------------------------------------------
        ! N_nodes     Number of nodes of the network
        ! edge_list   Array with shape (N_nodes, 2) where all the edges are stored
        !-----------------------------------------------------------------------
        ! OUTPUT
        !-----------------------------------------------------------------------
        ! degree_list One dimensional array with length N_nodes
        ! 

        implicit none
        integer, intent(in)                 :: N_nodes
        integer, dimension(:,:), intent(in) :: edge_list
        integer, dimension(N_nodes), intent(out) :: degree_list
        integer                             :: i
    
        degree_list=0
        do i=1, N_nodes
            degree_list(i)=count(edge_list(:,:).eq.i)
        enddo   
    
    end subroutine create_degree_list
    
    subroutine create_pointers(N_nodes, degree_list, pointers)
        ! Creates pointers to make possible read the neighbors array. This array has
        ! the first position of the neighbors array where the neighbors of a certain
        ! node begins and the position where the neighbors of this node ends. In fact,
        ! it is where the neighbors ends plus one to make easy use them in Python
        !-----------------------------------------------------------------------
        ! INPUT
        !-----------------------------------------------------------------------
        ! N_nodes     Number of nodes of the network
        ! degree_list One dimensional array with length N_nodes and information of the 
        !             edges of each node
        !-----------------------------------------------------------------------
        ! OUTPUT
        !-----------------------------------------------------------------------
        ! pointers    Two dimensional array with shape (N_nodes, 2)
        ! 
        implicit none
        integer, intent(in)                         :: N_nodes
        integer, dimension(:), intent(in)     :: degree_list
        integer, dimension(N_nodes, 2), intent(out) ::  pointers
        integer                                     :: i
    
        pointers(1,1)=1
        pointers(1,2)=degree_list(1)+1
        do i=2, N_nodes
            pointers(i, 1)=pointers(i-1,2)
            pointers(i,2)=pointers(i, 1)+ degree_list(i)
        enddo
    
    endsubroutine create_pointers
    
    subroutine create_neighbors(nlines,edge_list, pointers, neighbors)
        ! Creates neighbors array where all neighbors for all the nodes are stored
        ! The neighbors of node i will be neighbors(pointers(i,1):ponters(i,2)-1)
        !-----------------------------------------------------------------------
        ! INPUT
        !-----------------------------------------------------------------------
        ! nlines      Number of edges of the network
        ! edge_list   Array with shape (N_nodes, 2) where all the edges are stored
        ! pointers    Two dimensional array with shape (N_nodes, 2)
        !-----------------------------------------------------------------------
        ! OUTPUT
        !-----------------------------------------------------------------------
        ! neighbors   One dimensional array length 2*nlines
        ! 
        implicit none
        integer, intent(in)                       :: nlines
        integer, dimension(:,:), intent(in)  :: edge_list
        integer, dimension(:,:), intent(in) :: pointers
        integer, dimension(2*nlines), intent(out) ::  neighbors
        integer                                   :: i
        integer                                   :: position_1, position_2, n_1, n_2
        neighbors=0
    
        do i=1, nlines
            n_1=edge_list(i,1)
            position_1=pointers(n_1, 1)
            
            do while (.True.)
                if (neighbors(position_1).eq. 0) then
                    neighbors(position_1)=edge_list(i,2)
                    exit
                else
                    position_1=position_1+1
                endif
            enddo
    
            n_2=edge_list(i,2)
            position_2=pointers(n_2, 1)
            do while (.True.)
                if (neighbors(position_2).eq. 0) then
                    neighbors(position_2)=edge_list(i,1)
                    exit
                else
                    position_2=position_2+1
                endif
            enddo
        enddo
    endsubroutine create_neighbors
    

end module