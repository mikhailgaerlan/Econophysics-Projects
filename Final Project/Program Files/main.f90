program main
  use rands
  use convert
  use gen_tools
  use omp_lib
  implicit none
  
  integer :: popsize,indsize,objectives,ngen
  logical :: prints,seeded
  real :: start,finish
  real(kind=8) :: cxpb,distance,unifpb
  real(kind=8) :: mutpb,flipprob,mu,sigma
  real(kind=8) :: indtotal
  real(kind=8),dimension(:,:),allocatable :: pop,fitness,offspring,offfitness,total,totalfitness
  real(kind=8) :: rand
  real(kind=8),dimension(10) :: returns,sds
  real(kind=8),dimension(90) :: correlations
  integer :: i,j,k,n,prntrows
  integer :: c1,c2,cr
  character(len=32) :: numstr,prntformat,genstr
  character(len=32) :: funcname,crossover,mutation
  namelist /parameters/prints,seeded,popsize,ngen,funcname,&
       crossover,cxpb,distance,unifpb,mutation,mutpb,flipprob,mu,sigma,&
       returns,sds,correlations
  call cpu_time(start)
  call system_clock(c1)
  call system_clock(count_rate=cr)
  
  open(20,file='generations.dat',status='replace')
  
  !=============================
  !  Parameter Initialization
  !=============================
  open(10,file='parameters.txt',status='old')
  read(10,nml=parameters)
  
  if (seeded) call set_seed
  
  if (funcname.eq.'ackley') then
     indsize = 2
     objectives = 1
  elseif (funcname.eq.'binhkorn') then
     indsize = 2
     objectives = 2
  elseif (funcname.eq.'kursawe') then
     indsize = 3
     objectives = 2
  elseif (funcname.eq.'econo') then
     indsize=0
     do i=1,10
        if(returns(i).ne.0.0d0) indsize=indsize+1
     enddo
     objectives = 2
  endif

  allocate(pop(popsize,indsize),fitness(popsize,objectives+3))
  allocate(offspring(popsize,indsize),offfitness(popsize,objectives+3))
  allocate(total(2*popsize,indsize),totalfitness(2*popsize,objectives+3))
  prntrows = indsize+objectives+3
  write(numstr,'(i4)') prntrows
  prntformat = '('//trim(adjustl(numstr))//'f10.4)'
  
  !=============================
  !  Population Initialization
  !=============================
  if (prints) write(*,*) "Initializing population..."
  do j=1,indsize
     !$omp parallel shared(pop)
     !$omp do
     do i=1,popsize
        call random_number(rand)
        pop(i,j) = 10.0d0*(rand-0.5d0)
        if (funcname.eq.'econo') pop(i,j) = rand
     enddo
     !$omp end do
     !$omp end parallel
  enddo
  !$omp parallel shared(pop,fitness)
  !$omp do
  do i=1,popsize
     if ((funcname.eq.'binkhorn').or.(funcname.eq.'kursawe')) then
        fitness(i,1) = f1(pop(i,:))
        fitness(i,2) = f2(pop(i,:))
        fitness(i,objectives+3) = 1.0d0
     elseif (funcname.eq.'econo') then
        indtotal = 0.0d0
        do j=1,indsize
           indtotal = indtotal + pop(i,j)
        enddo
        if (indtotal.le.1.0d0) then
           fitness(i,1) = f1(pop(i,:))
           fitness(i,2) = f2(pop(i,:))
           fitness(i,objectives+3) = 1.0d0
        else
           fitness(i,1) = f1(pop(i,:))+1
           fitness(i,2) = f2(pop(i,:))+1
           fitness(i,objectives+3) = 1.0d0
        endif
     endif
  enddo
  !$omp end do
  !$omp end parallel
  if (prints) write(*,*) "Initializing done."
  
  call rank(fitness)
  call distancing(fitness)
  call sort(pop,fitness)
  do i=1,popsize
     write(20,prntformat) pop(i,:),fitness(i,:)
  enddo
  !write(*,*) ''

  !=============================
  !      Begin Algorithm
  !=============================
  do n=1,ngen
     if (prints) write(*,*) "Generation ", n
     offspring = pop
     offfitness = fitness
     
     !=============================
     !    Crossover Operation
     !=============================
     if (prints) write(*,*) "Crossing..."
     call distancing(offfitness)
     !$omp parallel shared(offspring,offfitness)
     !$omp do
     do i=1,popsize
        do j=i,popsize
           call random_number(rand)
           if ((rand<cxpb).and.(abs(offfitness(i,objectives+1)-offfitness(j,objectives+1)).ge.distance)) then
              do k=1,indsize
                 if (crossover.eq.'oneptcross') call oneptcross(pop(i,k),pop(j,k),offspring(i,k),offspring(j,k))
                 if (crossover.eq.'twoptcross') call twoptcross(pop(i,k),pop(j,k),offspring(i,k),offspring(j,k))
                 if (crossover.eq.'unifcross') call unifcross(pop(i,k),pop(j,k),offspring(i,k),offspring(j,k),unifpb)
                 if ((funcname.eq.'binhkorn').or.(funcname.eq.'kursawe')) then
                    if (offspring(i,k)<-5.0d0) offspring(i,k)=-5.0d0
                    if (offspring(i,k)>5.0d0) offspring(i,k)=5.0d0
                    if (offspring(j,k)<-5.0d0) offspring(j,k)=-5.0d0
                    if (offspring(j,k)>5.0d0) offspring(j,k)=5.0d0
                 elseif (funcname.eq.'econo') then
                    if (offspring(i,k)<0.0d0) offspring(i,k)=0.0d0
                    if (offspring(i,k)>1.0d0) offspring(i,k)=1.0d0
                    if (offspring(j,k)<0.0d0) offspring(j,k)=0.0d0
                    if (offspring(j,k)>1.0d0) offspring(j,k)=1.0d0
                 endif
                 enddo
              offfitness(i,4) = 0.0d0
              offfitness(j,4) = 0.0d0
           endif
        enddo
     enddo
     !$omp end do
     !$omp end parallel
     if (prints) write(*,*) "Crossing done."
     
     !=============================
     !           Mutation
     !=============================
     !$omp parallel shared(offspring,offfitness)
     !$omp do
     do i=1,popsize
        call random_number(rand)
        if (rand<mutpb) then
           do k=1,indsize
              if (mutation.eq.'mutflipbit') call mutgauss(offspring(i,k),mu,sigma)
              if (mutation.eq.'mutgauss') call mutflipbit(offspring(i,k),flipprob)
              if ((funcname.eq.'binhkorn').or.(funcname.eq.'kursawe')) then
                 if (offspring(i,k)<-5.0d0) offspring(i,k)=-5.0d0
                 if (offspring(i,k)>5.0d0) offspring(i,k)=5.0d0
              elseif (funcname.eq.'econo') then
                 if (offspring(i,k)<0.0d0) offspring(i,k)=0.0d0
                 if (offspring(i,k)>1.0d0) offspring(i,k)=1.0d0
              endif
           enddo
           offfitness(i,4) = 0.0d0
        endif
     enddo
     !$omp end do
     !$omp end parallel
     
     !=============================
     !     Calculate Fitness
     !=============================
     if (prints) write(*,*) "Evaluating fitnesses..."
     !$omp parallel shared(offspring,offfitness)
     !$omp do
     do i=1,popsize
        if ((funcname.eq.'binkhorn').or.(funcname.eq.'kursawe')) then
           offfitness(i,1) = f1(pop(i,:))
           offfitness(i,2) = f2(pop(i,:))
           offfitness(i,objectives+3) = 1.0d0
        elseif (funcname.eq.'econo') then
           indtotal = 0.0d0
           do j=1,indsize
              indtotal = indtotal + offspring(i,j)
           enddo
           if (indtotal.le.1.0d0) then
              offfitness(i,1) = f1(pop(i,:))
              offfitness(i,2) = f2(pop(i,:))
              offfitness(i,objectives+3) = 1.0d0
           else
              offfitness(i,1) = f1(pop(i,:))+1
              offfitness(i,2) = f2(pop(i,:))+1
              offfitness(i,objectives+3) = 1.0d0
           endif
        endif
     enddo
     !$omp end do
     !$omp end parallel
     if (prints) write(*,*) "Evaluating done."
     
     !=============================
     !          Sorting
     !=============================
     total(:popsize,:) = pop
     totalfitness(:popsize,:) = fitness
     total(popsize+1:,:) = offspring
     totalfitness(popsize+1:,:) = offfitness
     if (prints) write(*,*) "Ranking..."
     call rank(totalfitness)
     if (prints) write(*,*) "Ranking done"
     if (prints) write(*,*) "Sorting..."
     call dominsort(total,totalfitness)
     if (prints) write(*,*) "Sorting done."
     pop = total(:popsize,:)
     fitness = totalfitness(:popsize,:)
     if (prints) write(*,*) ""
     do i=1,popsize
        write(20,prntformat) pop(i,:),fitness(i,:)
     enddo
  enddo
  close(20)
  call cpu_time(finish)
  call system_clock(c2)
  if (prints) write(*,*) "System Time = ", (c2 - c1)/real(cr), " seconds."
  if (prints) write(*,*) "CPU Time = ", finish-start, " seconds."
  
contains

  function f1(x)
    real(kind=8) :: f1
    real(kind=8),intent(in) :: x(:)
    real(kind=8),parameter :: pi=4.0d0*atan(1.0d0)
    integer :: i
    
    f1=0.0d0
    if (funcname.eq.'ackley') then
       f1 = -20.0d0*exp(-0.2d0*sqrt(0.5d0*(x(1)**2.0d0+x(2)**2.0d0)))-&
            exp(0.5d0*(cos(2*pi*x(1))+cos(2*pi*x(2))))+exp(1.0d0)+20
    elseif (funcname.eq.'binhkorn') then
       f1 = 4.0d0*(x(1)**2.0d0)+4.0d0*(x(2)**2.0d0) !Binh-Korn 1
    elseif (funcname.eq.'kursawe') then
       do i=1,2
          f1 = f1-10.0d0*exp(-0.2d0*sqrt(x(i)**2.0d0+x(i+1)**2.0d0)) !Kursawe 1
       enddo
    elseif (funcname.eq.'econo') then
       do i=1,indsize
          f1 = f1-x(i)*returns(i)
       enddo
    endif
    
  end function f1

  function f2(x)
    real(kind=8) :: f2
    real(kind=8),intent(in) :: x(:)
    integer :: i,j,k
    
    f2=0.0d0
    if (funcname.eq.'binhkorn') then
       f2 = (x(1)-5.0d0)**2.0d0+(x(2)-5.0d0)**2.0d0 !Binh-Korn 2
    elseif (funcname.eq.'kursawe') then
       do i=1,3
          f2 = f2+abs(x(i))**(0.8d0)+5.0d0*sin(x(i)**3.0d0) !Kursawe 2
       enddo
    elseif (funcname.eq.'econo') then
       do i=1,indsize
          f2 = f2+(x(i)**(2.0d0))*(sds(i)**(2.0d0))
       enddo
       do i=1,indsize
          do j=1,indsize
             if (j.gt.i) then
                k = (i-1)*(indsize-1)+j-1
             elseif (j.lt.i) then
                k = (i-1)*(indsize-1)+j
             endif
             if (i.ne.j) f2 = f2+x(i)*x(j)*sds(i)*sds(j)*correlations(k)
          enddo
       enddo
       f2 = sqrt(f2)
    endif
  end function f2
  
end program main
