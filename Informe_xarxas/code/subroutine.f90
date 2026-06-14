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
    integer, intent(in)                      :: N_nodes
    integer, dimension(:,:), intent(in)      :: edge_list
    integer, dimension(N_nodes), intent(out) :: degree_list
    integer                                  :: i

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
    integer, dimension(:,:), intent(in)       :: edge_list
    integer, dimension(:,:), intent(in)       :: pointers
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

subroutine create_knn(max_k, N_nodes, neighbors, degree_list, pointers,k_nn)
    ! Computation of the ANND function normalized by <k²>/<k> using the defintion
    ! of the Ekk' matrix and the probabilities P(k,k')
    !-----------------------------------------------------------------------
    ! INPUT
    !-----------------------------------------------------------------------
    ! max_k       Integer. Maximum value of neighbors in the network
    ! N_nodes     Integer. Number of nodes of the network
    ! degree_list One dimensional array with length N_nodes and information of the 
    !             edges of each node
    ! pointers    Two dimensional array with shape (N_nodes, 2)
    ! neighbors   One dimensional array length 2*nlines
    !-----------------------------------------------------------------------
    ! OUTPUT
    !-----------------------------------------------------------------------
    ! k_nn        One dimensional array length max_k
    ! 
    
    implicit none
    integer, intent(in)                           :: max_k, N_nodes
    integer, dimension(:,:), intent(in)           :: pointers
    integer, dimension(:), intent(in)             :: neighbors, degree_list
    integer, dimension(max_k, max_k)              :: E
    real(kind=8), dimension(max_k,max_k)          :: P_kprima_k
    real(kind=8), dimension(max_k)                :: Pk
    integer, dimension(max_k)                     :: N_k
    real(kind=8), dimension(max_k), intent(out)   :: k_nn
    integer                                       :: nod_neig, i,j, m, k_prima, k
    real(kind=8)                                  :: mean_val_k, mean_val_k_sq, aux

    E=0
    
    do i=1, size(neighbors)
        nod_neig=neighbors(i)
        k_prima=degree_list(nod_neig)
        m=1
        do while (.True.)
            
            if (i.ge. pointers(m, 1) .and. i.lt. pointers(m, 2)) then
                k=degree_list(m)
                
                E(k, k_prima)= E(k, k_prima)+1
                
                exit
            endif
            m=m+1
        enddo
       
    enddo

    do i=1, max_k
        N_k(i)=count(degree_list(:).eq.i)
    enddo
    P_kprima_k=0.d0
    do i=1, max_k
        do j= 1,max_k
            if (E(j,i).ne.0) then
                P_kprima_k(j,i)= dble(E(j,i))/dble(i*N_k(i))
                
            endif
        enddo
    enddo

    do i=1, max_k
    k_nn(i)=0.d0
        do j=1,max_k
            k_nn(i)=k_nn(i)+j *P_kprima_k(j,i)
        enddo
    enddo

    !Normalization
    !<k>
    mean_val_k=dble(sum(E(:,:)))/dble(N_nodes)
    do i=1, max_k
        Pk(i)=dble(count(degree_list(:).eq.i))/dble(size(degree_list))
    enddo
    
    !<k²>
    mean_val_k_sq=0.d0
    do i=1, max_k
        aux=0.d0
        do j=1, max_k
            aux=aux+j*P_kprima_k(j,i)
        enddo
        mean_val_k_sq=mean_val_k_sq+i*Pk(i)*aux
    enddo

    k_nn=k_nn/(mean_val_k_sq/mean_val_k)
endsubroutine create_knn

subroutine knn(max_k, N_nodes, neighbors, degree_list, pointers,k_nn)
    ! Computation of the ANND function normalized by <k²>/<k>. It have been 
    ! checked that it give the same result as the previous function
    !-----------------------------------------------------------------------
    ! INPUT
    !-----------------------------------------------------------------------
    ! max_k       Integer. Maximum value of neighbors in the network
    ! N_nodes     Integer. Number of nodes of the network
    ! degree_list One dimensional array with length N_nodes and information of the 
    !             edges of each node
    ! pointers    Two dimensional array with shape (N_nodes, 2)
    ! neighbors   One dimensional array length 2*nlines
    !-----------------------------------------------------------------------
    ! OUTPUT
    !-----------------------------------------------------------------------
    ! k_nn        One dimensional array length max_k
    ! 
    implicit none
    integer, intent(in)                           :: max_k, N_nodes
    integer, dimension(:,:), intent(in)           :: pointers
    integer, dimension(:), intent(in)             :: neighbors, degree_list
    integer                                       :: i, k,n
    integer, dimension(N_nodes)                   :: nodes
    integer, dimension(:), allocatable            :: nodes_k, neig, k_sum
    real(kind=8)                                  :: sumation, k_mean, k_sq_mean
    real(kind=8), dimension(max_k), intent(out)   :: k_nn

    nodes=(/(i, i=1,N_nodes)/)
    
    do k=1, max_k
        nodes_k=pack(nodes, degree_list.eq.k)

        sumation=0.d0
        if (size(nodes_k).gt.0) then
            do n=1, size(nodes_k)
                allocate(neig(k))
                allocate( k_sum(k))
                neig(:)=neighbors(pointers(nodes_k(n), 1):pointers(nodes_k(n), 2)-1)
                forall (i=1:size(neig)) k_sum(i)=degree_list(neig(i)) 
                sumation=sumation + dble(sum(k_sum))/dble(k)
                deallocate(neig, k_sum)
            enddo
        
    
        k_nn(k)=sumation/dble(size(nodes_k))
        else
            k_nn(k)=0.d0
        endif
    enddo

    !Normalization
    k_mean=dble(sum(degree_list))/dble(size(degree_list))
    k_sq_mean=dble(sum(degree_list**2))/dble(size(degree_list))
    k_nn=k_nn/(k_sq_mean/k_mean)
    
endsubroutine


subroutine global_cluster_coef(N_nodes, pointers, neighbors, global_c)
    ! Computation of the global clustering coefficient computing all the 
    ! triangles and all the connected triads (edges that share a vertex).
    ! This is also known as transitivity
    !-----------------------------------------------------------------------
    ! INPUT
    !-----------------------------------------------------------------------
    ! max_k       Integer. Maximum value of neighbors in the network
    ! N_nodes     Integer. Number of nodes of the network
    ! pointers    Two dimensional array with shape (N_nodes, 2)
    ! neighbors   One dimensional array length 2*nlines
    !-----------------------------------------------------------------------
    ! OUTPUT
    !-----------------------------------------------------------------------
    ! global_c    Global cluster coefficient
    ! 
    implicit none
    integer, intent(in)                 :: N_nodes
    integer, dimension(:,:), intent(in) :: pointers
    integer, dimension(:), intent(in)   :: neighbors
    real(kind=8), intent(out)           :: global_c
    integer                             :: triangles, i,j,k, triad


    triangles=0; triad=0
    
    do i=1, N_nodes-2
        do j=i+1, N_nodes-1
            do k=j+1, N_nodes
                ! Check all the triangles
                if (any(neighbors(pointers(i,1):pointers(i,2)-1).eq. j) .and.&
                any(neighbors(pointers(j,1):pointers(j,2)-1).eq. k) .and. &
                any(neighbors(pointers(k,1):pointers(k,2)-1).eq. i)) then
                    triangles=triangles+1
                else
                    ! Check all the triads
                
                    if (any(neighbors(pointers(i,1):pointers(i,2)-1).eq. j) .and.&
                    any(neighbors(pointers(j,1):pointers(j,2)-1).eq. k)) then
                        triad=triad+1
                
                    else if (any(neighbors(pointers(i,1):pointers(i,2)-1).eq. k) .and.&
                    any(neighbors(pointers(k,1):pointers(k,2)-1).eq. j)) then
                        triad=triad+1

                    else if (any(neighbors(pointers(j,1):pointers(j,2)-1).eq. i) .and.&
                    any(neighbors(pointers(i,1):pointers(i,2)-1).eq. k) ) then
                        triad=triad+1
                        
                     endif

            
                endif
            enddo
        enddo
    enddo

    global_c = 3.d0*dble(triangles)/dble(3*triangles+(triad))
   
endsubroutine

subroutine cluster_coef(max_k,degree_list, neighbors, pointers, c_k)
    ! Computation of the clustering coefficient dependent of the deighbors
    ! computing all the triangles with the nearest neighbors for all the
    ! nodes of degree class k
    !-----------------------------------------------------------------------
    ! INPUT
    !-----------------------------------------------------------------------
    ! max_k       Integer. Maximum value of neighbors in the network
    ! degree_list One dimensional array with length N_nodes and information of the 
    !             edges of each node
    ! pointers    Two dimensional array with shape (N_nodes, 2)
    ! neighbors   One dimensional array length 2*nlines
    !-----------------------------------------------------------------------
    ! OUTPUT
    !-----------------------------------------------------------------------
    ! c_k         One dimensional array length max_k with information of the 
    !             clustering
    ! 
    implicit none
    integer, intent(in)                         :: max_k
    integer, dimension(:), intent(in)           :: degree_list
    integer, dimension(:,:), intent(in)         :: pointers
    integer, dimension(:), intent(in)           :: neighbors
    integer, dimension(:), allocatable          :: mask, nodes
    real(kind=8), dimension(max_k), intent(out) :: c_k
    integer                                     :: k, j, i, l, ni, nj, nl, triangles



    c_k=0.d0
    do k=2, max_k
        mask= merge(1,0, degree_list.eq.k)
       
        if (sum(mask).eq. 0) then
            cycle
        endif
        allocate(nodes(sum(mask)))

        j=1
        do i=1, size(mask)
            if (mask(i).eq. 1) then
                nodes(j)=i
                j=j+1
            endif
        enddo

        triangles=0
        do i=1, size(nodes) 
            ni=nodes(i)
            ! Check the triangles
            do j=1, k-1
                do l=j+1, k
                    
                    nj= neighbors(pointers(ni,1)+j-1)
                    nl= neighbors(pointers(ni,1)+l-1)
                    if (any(neighbors(pointers(nj,1):pointers(nj,2)-1).eq. nl)) then
                        triangles=triangles+1
                    endif

                enddo
            enddo
        enddo
        c_k(k)= 2.d0/dble(k*(k-1)* size(nodes))*triangles
        
        deallocate(nodes)
        

enddo

endsubroutine cluster_coef

subroutine cm_creation(nodes, link_list)
    ! Creates networks using the configuration model
    !-----------------------------------------------------------------------
    ! INPUT
    !-----------------------------------------------------------------------
    ! nodes       One dimendional array of length 2*E where E is the number of
    !             links of the network. It has a speacial construction, each node 
    !             label is repited as meny times as its degree. See the python 
    !             script assignment5.py for more information
    !-----------------------------------------------------------------------
    ! OUTPUT
    !-----------------------------------------------------------------------
    ! link_list   Edge list of the network
    ! 
    implicit none
    integer, dimension(:), intent(in)                  :: nodes
    integer, dimension(:), allocatable                 :: nodes_aux1, nodes_aux2
    integer, dimension((size(nodes)/2),2), intent(out) :: link_list
    integer                                            :: nod_ran1, nod_ran2, link_num,n, nod_aux,&
     nod_pos1, nod_pos2,i, m
    real(kind=8)                                       :: ran1, ran2
    integer, dimension(size(link_list(:,1))) :: test

    link_num=size(nodes)
    link_list=0
    n=1; m=0
    allocate(nodes_aux1(link_num))
    nodes_aux1=nodes


    do while (n.le. size(link_list(:,1)))
        allocate(nodes_aux2(link_num))
        nodes_aux2=nodes_aux1
        call random_number(ran1);call random_number(ran2)
    
        nod_pos1= mod(int(link_num*ran1),link_num)+1
        nod_pos2= mod(int(link_num*ran2),link_num)+1

        nod_ran1=nodes_aux1(nod_pos1)
        nod_ran2=nodes_aux1(nod_pos2)
        deallocate(nodes_aux1)

        if (nod_ran1 .gt. nod_ran2) then
             nod_aux=nod_ran2 
             nod_ran2=nod_ran1
             nod_ran1=nod_aux
        endif
        m=m+1
        forall (i=1:size(link_list(:,1))) test(i) = sum(abs(link_list(i,:) -(/nod_ran1,nod_ran2 /)))
        
        if (nod_ran1.ne. nod_ran2  .and. all(test.ne.0)) then !Check self-loop and multiconnexions
            link_list(n,:)=(/nod_ran1, nod_ran2 /)
            where ((/(i, i=1, link_num)/).eq. nod_pos1 .or.(/(i, i=1, link_num)/).eq. nod_pos2) nodes_aux2=-1
            n=n+1
            link_num=link_num-2
            
            allocate(nodes_aux1(link_num))
            nodes_aux1=pack(nodes_aux2, nodes_aux2.gt. 0)
            deallocate(nodes_aux2)
            
        else
            allocate(nodes_aux1(link_num))
            nodes_aux1=nodes_aux2
            deallocate(nodes_aux2)
            
        endif

        if (m.gt.100000) exit
    enddo


endsubroutine cm_creation