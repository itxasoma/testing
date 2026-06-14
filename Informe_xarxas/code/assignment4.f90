subroutine SW_model(N_nodes, probability, edge_list)
    ! Creates networks using the small world
    !-----------------------------------------------------------------------
    ! INPUT
    !-----------------------------------------------------------------------
    ! N_nodes     Number of nodes of the network
    ! probability Probability of doing the rewiring
    !-----------------------------------------------------------------------
    ! OUTPUT
    !-----------------------------------------------------------------------
    ! edge_list   Edge list of the network
    ! 
    implicit none
    integer, intent(in) :: N_nodes
    real(kind=8), intent(in) :: probability
    integer, dimension(N_nodes):: first_array, second_array, third_array
    integer, dimension(N_nodes*2, 2), intent(out) :: edge_list
    integer, dimension(N_nodes*2) :: test1, test2

    integer                          :: i, node, j
    real(kind=8)                     :: ran, ran_node

    ! Network with <k>=4
    first_array= (/(i, i= 1,N_nodes)/)
    second_array=cshift(first_array, shift=1)
    third_array=cshift(first_array, shift=2)

    edge_list(:N_nodes, 1)=first_array
    edge_list(:N_nodes, 2)=second_array
    edge_list(N_nodes+1:, 1)=first_array
    edge_list(N_nodes+1:, 2)=third_array


    do i=1,2*N_nodes
        call random_number(ran)
        if (ran.lt. probability) then
            call random_number(ran_node)
            node= mod(int(N_nodes*ran_node),N_nodes)+1

            forall (j=1:N_nodes*2) test1(j) = sum(abs(edge_list(j,:) -(/edge_list(i,1),node /)))
            forall (j=1:N_nodes*2) test2(j) = sum(abs(edge_list(j,:) -(/node, edge_list(i,2) /)))
        
            if (edge_list(i, 1).ne. node .and. all(test1.ne.0)& !Check self-loop and multiconnexions
             .and. all(test2.ne.0) ) then
                edge_list(i, 2) = node
            endif
        endif
    enddo



endsubroutine