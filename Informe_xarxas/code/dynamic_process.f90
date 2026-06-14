module SIS
    implicit none

    public :: SIS_process, creation_link_act

    contains
subroutine SIS_process(t_fin, lambda, pointers, neighbors, act_nodes, act_links,sim)
    use ::dp, only: exponential_distribution, append_int,rem_element

    implicit none
    integer, dimension(:,:), intent(in)                 :: pointers
    integer, dimension(:), intent(in)                   :: neighbors
    integer, dimension(:), intent(in)                   :: act_nodes
    integer, dimension(:,:), intent(in)                 :: act_links
    integer,allocatable, dimension(:)                   :: act_nodes_aux
    integer,allocatable, dimension(:,:)                 :: act_links_aux
    integer, allocatable, dimension(:)                  :: neig,vec,brk_neigh,pos_links_rem,nod_add
    real(kind=8), intent(in)                            :: t_fin
    integer, intent(in)                                 :: sim
    real(kind=8), parameter                             :: delta=1
    real(kind=8), intent(in)                            :: lambda
    integer                                             :: c, Elinks, n_i, i, val,val2,j,loc, unit_au
    integer, dimension(size(neighbors))                 :: point_act, resta_point_act
    real(kind=8)                                        :: prob_rec, ran, ran_2,t,rho
    logical, dimension(size(neighbors))                 :: mask
    logical, dimension(size(pointers,1))                :: mask_out
    integer, dimension(size(pointers,1))                :: n

    
    point_act=0
    t=0.d0
    allocate(act_links_aux(size(act_links,1), size(act_links,2))); allocate(act_nodes_aux(size(act_nodes)))
   
    !################################################################
    ! Create allocatable array to append and remove terms
    !################################################################

    act_links_aux= act_links; act_nodes_aux=act_nodes

    !################################################################
    ! Create array of pointers for the active links array
    ! It pount out to the positions of the array fo active nodes where
    ! a certain infected-susceptible node is located. It is realed to 
    ! the neighbors array 
    ! If the infected node is 2 and the active links are:
    !  2  4
    !  2  5
    !   :
    ! and the neighbors array is:
    ! ... 11 [4 5] 8...     where [] denote where the pointers of node 2 are
    ! the act_poit array will be:
    ! ... 0 [1 2] 0...  
    !################################################################
    do i=1, size(act_nodes_aux)
        neig=neighbors(pointers(act_nodes_aux(i),1):pointers(act_nodes_aux(i),2)-1)
        
        do j=1, size(act_links_aux,1)
            if (any(act_links_aux(j,:).eq. act_nodes_aux(i))) then
                val=act_links_aux(j,1); val2=act_links_aux(j,2)
                if (val.eq.act_nodes_aux(i)) val=val2
                point_act(findloc(neig,val,dim=1)+pointers(act_nodes_aux(i),1)-1)=j
            endif
            
        enddo
    enddo

  
    call name_ouput_file(int(1000.d0*lambda), unit_au, sim)
    write(unit_au,*) '#time, rho'
    c=0
    !##################################################################
    ! Main loop
    !##################################################################
    do while (t.lt. t_fin)
        Elinks=size(act_links_aux,1) !number of active links
        n_i=size(act_nodes_aux) !number of infected nodes

        !##################################################################
        ! Update of the time
        !##################################################################
    
        t=t+exponential_distribution(n_i+lambda*Elinks)


        !##################################################################
        ! Probability of recovery
        !##################################################################
        prob_rec=dble(n_i)/dble(n_i+lambda*Elinks)

        call random_number(ran)
          !------------------------------------------------------------------------------------
          ! The node infects
        if (ran.ge.prob_rec) then
            
            call random_number(ran_2)
            i= int(ran_2*dble(size(act_links_aux,1)-1))+1
            ! Choose of the active link that will create the infection
            val=act_links_aux(i,1);val2=act_links_aux(i,2); if (any(act_nodes_aux.eq. val)) val=val2

            !##################################################################
            !Broken active links due to that infection
            !##################################################################

            mask=.False.
            where (neighbors.eq.val) mask=.True.
            !print*, mask
            brk_neigh=pack(point_act, mask.and. point_act.ne.0)
            
            !##################################################################
            !Elimination of the old active links
            !##################################################################
           
            call rem_element(act_links_aux, brk_neigh, 'index')
            !print*, mask
            where (mask) point_act=0
            resta_point_act=0
            do i=1, size(brk_neigh)
                where (point_act.gt.brk_neigh(i)) resta_point_act=resta_point_act-1
            enddo
            point_act=point_act+resta_point_act

            
            !------------------------------------------------------------------------------------
            !##################################################################
            ! Append the node as infected
            !##################################################################

            call append_int(act_nodes_aux, val)

            !################################################################
            ! Neighbours of the new infected node
            !################################################################

            
            
            allocate(vec(size(neighbors(pointers(val, 1):pointers(val, 2)-1))))
            vec=neighbors(pointers(val, 1):pointers(val, 2)-1)

            !################################################################
            ! Append new active links
            !################################################################


            do i=1, size(vec)
                if (all(act_nodes_aux.ne.vec(i))) then
                    call append_int(act_links_aux, val, vec(i))
                    point_act(pointers(val, 1)+i-1)=size(act_links_aux,1)
                endif
            enddo
            deallocate(vec)

            

        else 
            !################################################################
            ! Recovery of a node
            !################################################################
            
            call random_number(ran_2)
            ! Choose node
            i= int(ran_2*dble(size(act_nodes_aux)-1))+1
            val=act_nodes_aux(i)

            !################################################################
            ! Remove the node from the infected nodes
            !################################################################

            call rem_element(act_nodes_aux,val, 'value')

            !################################################################
            ! Locate old active links
            !################################################################
            pos_links_rem=pack(point_act(pointers(val, 1):pointers(val, 2)-1), point_act(pointers(val, 1):pointers(val, 2)-1).ne.0)
           
            !################################################################
            ! Choose new active links that will be append
            !################################################################
            nod_add=pack(neighbors(pointers(val, 1):pointers(val, 2)-1), point_act(pointers(val, 1):pointers(val, 2)-1).eq.0)
            

            !################################################################
            !Remove old active links
            !################################################################
            call rem_element(act_links_aux, pos_links_rem, 'index')
            point_act(pointers(val, 1):pointers(val, 2)-1)=0
            resta_point_act=0
            do i=1, size(pos_links_rem)
                where (point_act.gt.pos_links_rem(i)) resta_point_act=resta_point_act-1
            enddo
            point_act=point_act+resta_point_act
            
           
            !################################################################
            ! Append new active links to the list
            !################################################################

            do i=1, size(nod_add)
                val2=nod_add(i)
                call append_int(act_links_aux,val, val2)
                
                loc=findloc(neighbors(pointers(val2, 1):pointers(val2, 2)-1), val, dim=1)
                point_act(pointers(val2, 1)+loc-1)=size(act_links_aux,1)
            enddo

        endif
        
        mask_out=.False.
        
        mask_out(act_nodes_aux)=.True.

        !##################################################################
        !Calculation of the prevalence for all nodes as rho=<n>
        !##################################################################
    
        n=mask_out ; rho=dble(sum(n))/dble(size(n) )
        

        
        c=c+1
        write(unit_au,*) t, rho
        if (rho<1d-14) exit ! The loop stops if rho is near zero to avoid Infinity
        if (c.gt.100000) stop
    enddo

 close(unit_au)

endsubroutine

subroutine creation_link_act(pointers, neighbors, act_nodes, link_act)
    ! Creates the array of active links it has a shape (E_{act}, 2)
    !-----------------------------------------------------------------------
    ! INPUT
    !-----------------------------------------------------------------------
    ! pointers    Two dimensional array with shape (N_nodes, 2)
    ! neighbors   One dimensional array length 2*nlines
    ! act_nodes   Array of infected nodes
    !-----------------------------------------------------------------------
    ! OUTPUT
    !-----------------------------------------------------------------------
    ! act_list    Two dimensional array of active links shape (E_{act}, 2)
    ! 
    use ::dp, only:  append_int
    implicit none
    integer, dimension(:,:), intent(in)                 :: pointers
    integer, dimension(:), intent(in)                   :: neighbors
    integer, dimension(:), intent(in)                   :: act_nodes
    integer, dimension(:), allocatable                  :: nod_add,neig, nodes_inf
    integer                                             :: i,j,node, val2,n
    integer, dimension(:,:),allocatable, intent(out)    :: link_act

        allocate(link_act(1,2))
        link_act=0               ! active edges array
        n=0
        ! Creation of the array needed for the simulation
        ! act_nodes array has the label of the infected nodes
        do j=1, size(act_nodes)
            node=act_nodes(j)
            neig=neighbors(pointers(node, 1):pointers(node, 2)-1)      !neighbors of each the active nodes
            
            allocate(nodes_inf(size(neig)));nodes_inf=0                ! neighbors that can be or not infected
            do i=1, size(neig)
                if (all(act_nodes.ne. neig(i))) nodes_inf(i)= neig(i)
            enddo
            nod_add=pack(nodes_inf, nodes_inf.ne.0)                    ! no infected nodes in an active link
           
            deallocate(nodes_inf)

            do i=1, size(nod_add)
                val2=nod_add(i)
                if (n.ne.0) then; call append_int(link_act,node, val2) ! append to the array of active links
                else; link_act(1,:)=(/node, val2/); endif              ! initial creation of the array
                
                n=n+1
            enddo
        enddo
endsubroutine

subroutine name_ouput_file(lambda, fu, sim,str)
    ! Gives names to the output file depending on the value of lambda parameter
    ! and the repetition of the simulation. Also a possible extra variable str
    !-----------------------------------------------------------------------
    ! INPUT
    !-----------------------------------------------------------------------
    ! lambda      Rate of infection
    ! fu          Unit of the file
    ! sim         Optional number of the repetition of the simulation 
    ! str         Optional, used to produce files with the Runge-Kutta labe
    ! 
    integer(kind=4), intent(in)       :: lambda, fu
    integer(kind=4), optional,intent(in)          :: sim
    character(*), optional,intent(in)          :: str
    character(len=1024)               :: filename
    character(len=10)                 :: file_id, file_id2

    write(file_id, '(i0)') lambda
   
    if (present(str)) then
        filename = 'output/'//'sim_out'//'_'//str//'_lambda' // trim(adjustl(file_id))//'.dat'
    else
        if (present(sim))then
            write(file_id2, '(i0)') sim
        filename = 'output/'//'sim_out_lambda' // trim(adjustl(file_id))//'_sim'//trim(adjustl(file_id2))//'.dat'
        endif
    endif

    filename= trim(filename)

    open(unit=fu, file=filename, action='write', status='replace')
end subroutine name_ouput_file

subroutine func_rk( lambda, rho, neighbors,pointers, rho_out)
    !Function used for the Runge-Kutta integration
    !-----------------------------------------------------------------------
    ! INPUT
    !-----------------------------------------------------------------------
    ! lambda      Rate of infection
    ! rho         Perseverance array with length N_nodes
    ! pointers    Two dimensional array with shape (N_nodes, 2)
    ! neighbors   One dimensional array length 2*nlines
    !-----------------------------------------------------------------------
    ! OUTPUT
    !-----------------------------------------------------------------------
    ! rho_out     Value of d rho/dt
    ! 
    implicit none
    real(kind=8), intent(in)                 :: lambda
    real(kind=8),dimension(:), intent(in)    :: rho
    integer,dimension(:), intent(in)         :: neighbors
    integer,dimension(:,:), intent(in)       :: pointers
    real(kind=8),dimension(:), intent(out)   :: rho_out
    real(kind=8)                             :: suma
    integer                                  :: i,j
    
    do i=1, size(rho)
        suma=0.d0
        do j= pointers(i, 1),  pointers(i, 2)-1
            suma=suma + rho(neighbors(j))
        enddo

        rho_out(i)= -rho(i)+lambda*suma*(1.d0-rho(i))
    enddo
endsubroutine


subroutine runge_kutta(t_fin,dt, lambda, rho, neighbors,pointers )
    ! Runge Kutta integration
    !-----------------------------------------------------------------------
    ! INPUT
    !-----------------------------------------------------------------------
    ! t_fin       Final time of the integraion
    ! dt          Time step of the integration
    ! rho         Perseverance array with length N_nodes
    ! pointers    Two dimensional array with shape (N_nodes, 2)
    ! neighbors   One dimensional array length 2*nlines
    !-----------------------------------------------------------------------
    ! OUTPUT
    !-----------------------------------------------------------------------
    ! rho_out     Value of d rho/dt
    ! 
    implicit none
    real(kind=8), intent(in)                 :: t_fin, dt, lambda
    real(kind=8),dimension(:), intent(in)    :: rho
    integer,dimension(:), intent(in)         :: neighbors
    integer,dimension(:,:), intent(in)       :: pointers
    ! Internal variables
    real(kind=8),dimension(size(rho))        :: u, u_n,aux
    real(kind=8),dimension(size(rho))        :: rho2, rho3, rho4, k1,k2,k3,k4
    real(kind=8)                             :: t
    integer                                  :: c,unit_au

        u_n=rho
        call name_ouput_file(int(1000.d0*lambda), unit_au, str='rk')
        t=0.d0
        c=0
        do while (t.lt.t_fin)
            call func_rk( lambda, u_n, neighbors,pointers, k1)
            rho2= u_n+ dt*k1*0.5d0
            
            call func_rk( lambda, rho2, neighbors,pointers, k2)
            rho3= u_n+ dt*k2*0.5d0
            
            call func_rk( lambda, rho3, neighbors,pointers, k3)
            rho4= u_n+ dt*k3

            call func_rk( lambda, rho4, neighbors,pointers, k4)

            u= u_n+ dt/6.d0 * (k1+2.d0*k2+2.d0*k3+k4)
            
            aux=u; u_n=aux; u=u_n
            t=t+dt
            c=c+1
            !Write in the file
            write(unit_au,*) t, sum(u_n)/dble(size(u_n))
            if (c.gt.100000) stop
        enddo
        close(unit_au)
endsubroutine


endmodule SIS

 program dynamic_process
      use ::dp, only: create_degree_list,create_pointers, create_neighbors
      use ::SIS, only:  SIS_process, creation_link_act,runge_kutta
      implicit none
      character(*), parameter    :: filename='edges_list.txt'
      integer, dimension(:,:), allocatable                 :: pointers
      integer, dimension(:), allocatable                   :: neighbors, degree_list
      integer, dimension(:), allocatable                   :: act_nodes
      integer, dimension(:,:), allocatable                 :: act_links, edge_list
      logical, dimension(:), allocatable                   :: mask_out
      real(kind=8), dimension(:), allocatable              :: rho
      integer                                              :: n_inf,nlines, N_nodes,io, i,j,sim, l
      real(kind=8)                                         :: ran, lambda, t_fin,dt

      ! Reading of the edge list

      open(10,file= filename, iostat=io, status='old')
      if (io/=0) stop 'Cannot open file! '
      nlines = 0
     
      do
          read(10,*,iostat=io)
          if (io/=0) exit
          nlines = nlines + 1
      end do
      close(10)
    
      allocate(edge_list(nlines,2))
   
      open(1, file = filename, action='read', status='old')
     
      do i = 1, nlines
          read(1, *) (edge_list(i,j), j=1,2)
      end do
      close(1)

      t_fin=15.d0
      dt=1d-2
     

    !-----------------------------------------------------------------------------------
    ! Creation of an initial configuration
    ! Number of initially infected nodes
      n_inf=50
      nlines=size(edge_list,1)
      N_nodes=int(maxval(edge_list))
    
      allocate(act_nodes(n_inf))
      i=1
      do while (i.le.n_inf)
          call random_number(ran)

            !Creation of the array of infected nodes avoiding repetitions          
            if (all(act_nodes.ne. int(ran*dble(N_nodes-1))+1))then
                act_nodes(i)= int(ran*dble(N_nodes-1))+1
                i=i+1
            endif
      enddo
     
      allocate( mask_out(N_nodes))
      mask_out=.False.; mask_out(act_nodes)=.True.
      

      rho=merge(1.d0,0.d0,mask_out) ! prevalence for the Runge Kuta
      
      allocate(degree_list(N_nodes), pointers(N_nodes,2), neighbors(2*nlines))

      print*, sum(rho) !Check the number i¡of infected nodes
      call create_degree_list(N_nodes,edge_list, degree_list)

      call create_pointers(N_nodes, degree_list, pointers)

      call create_neighbors(nlines,edge_list, pointers,neighbors)


      call creation_link_act(pointers, neighbors, act_nodes, act_links)   


    ! Loop over the possible valueas of lambda 
      do l=1, 40

      lambda=0.05d0+ dble(l*(4.0d0-0.05d0))/10.d0 !Lineary distributed lambda
      call runge_kutta(t_fin,dt, lambda, rho, neighbors,pointers ) !Runge Kutta integrator
      
          do sim=1, 10             
              call SIS_process(t_fin, lambda, pointers, neighbors, act_nodes, act_links, sim) !Stocastic process
          enddo
            
      enddo
  endprogram