program modeliiimain
use m_vib
implicit none

real(8),parameter :: pi = 3.1415926535d0, twopi = 2d0*pi

character(len=2) :: c_ng,c_nt
character(len=9) :: fmt1,fmt2

complex(8) :: coeff,fact,a1,a2,et
complex(8),dimension(:),allocatable :: pol_tot,x,p,rm,pm,f
complex(8),dimension(:,:),allocatable :: pol,hm

integer :: a,b,i,j,is,it,cnt,p_i,p_j,p_k,ib,nmap,ng,nb,nd,basispc
integer :: np,nmcs,mcs,nmds,seed_dimension,nosc,step1,bath,init,nfile
integer,dimension(:),allocatable :: seed,g

real(8) :: gauss,dt,dt2,kondo,delta,beta,ome_max,lumda_d,eg,eb,ed,mu,e0,e1,sij,vomega
real(8) :: step2,dnmcs,tau1,omega1,time3,lambdacheck
real(8),dimension(:),allocatable :: tau,time,omega,c2,kosc,ome
real(8),dimension(:,:),allocatable :: lambda,lmd,ug,ub,ud,hc
real(8),dimension(:,:),allocatable :: sgg,sgb,sgd,sbg,sbb,sbd,sdg,sdb,sdd
real(8),dimension(:,:),allocatable :: llg,llb,lld,llgb,llbg,llbd,lldb,hs

call iniconc()

nmap = ng + nb + nd

allocate(c2(1:nosc),kosc(1:nosc),ome(1:nosc),x(1:nosc),p(1:nosc),f(1:nosc))
allocate(tau(1:np),omega(1:np),time(1:np),g(1:np))
allocate(rm(1:nmap),pm(1:nmap))

allocate(sgg(1:ng,1:ng),sgb(1:ng,1:nb),sgd(1:ng,1:nd))
allocate(sbg(1:nb,1:ng),sbb(1:nb,1:nb),sbd(1:nb,1:nd))
allocate(sdg(1:nd,1:ng),sdb(1:nd,1:nb),sdd(1:nd,1:nd))
allocate(hm(1:nmap,1:nmap),hc(1:nmap,1:nmap),hs(1:nmap,1:nmap))
allocate(ug(1:nmap,1:nmap),ub(1:nmap,1:nmap),ud(1:nmap,1:nmap))
allocate(lambda(1:nmap,1:nmap),llg(1:nmap,1:nmap),llb(1:nmap,1:nmap),lld(1:nmap,1:nmap))
allocate(llgb(1:nmap,1:nmap),llbg(1:nmap,1:nmap))
allocate(llbd(1:nmap,1:nmap),lldb(1:nmap,1:nmap))

call iniconq_d()

dt  = twopi*dt
dt2 = 0.5d0*dt

do i = 1, np
   tau(i)   = tau1
   omega(i) = omega1
end do

time(1) = 0.3d0
time(2) = 0.3d0
time(3) = (time3 + nfile*step2)

nmds = nmds + nfile*step1

allocate(pol(1:nmds+1,1:2**np))
allocate(pol_tot(1:nmds+1))

pol = cmplx(0d0,0d0)

if (ng > 9) then
   write(c_ng,'(i2)') ng
else
   write(c_ng,'(i1)') ng
end if

if (nmap > 9) then
   write(c_nt,'(i2)') nmap
else
   write(c_nt,'(i1)') nmap
end if

fmt1 = '('//trim(c_ng)//'f10.5)'
fmt2 = '('//trim(c_nt)//'f10.5)'

call get_lambda_eigenvectors(ng,nb,nd,eg,eb,ed,delta,vomega, &
                              sgg,sgb,sgd,sbg,sbb,sbd,sdg,sdb,sdd,lambda,hs)


!check orthonormality of lambdas
do i = 1, nmap
   do j = 1, nmap
      lambdacheck = dot_product(lambda(1:nmap,i),lambda(1:nmap,j))
      if (lambdacheck > 1d-7) then
         if (i /= j) then
            print *, 'warning in orthonormality of these lambda eigenvectors'
            print *, i, j, lambdacheck
            print *, ''
            do is = 1, nmap
               print *, is
               print *, lambda(1:nmap,is)
            end do
         end if
      end if
   end do
end do

!parameters for faster calculation
llg = 0d0
do i = 1, nmap
   do j = 1, nmap
      llg(i,j) = dot_product(lambda(1:ng,i),lambda(1:ng,j))
   end do
end do

!write(*,*) 'products'
!write(*,fmt2) llg

llb = 0d0
do i = 1, nmap
   do j = 1, nmap
      llb(i,j) = dot_product(lambda(ng+1:ng+nb,i),lambda(ng+1:ng+nb,j))
   end do
end do

!write(*,*) 'products'
!write(*,fmt2) llb

lld = 0d0
do i = 1, nmap
   do j = 1, nmap
      lld(i,j) = dot_product(lambda(ng+nb+1:nmap,i),lambda(ng+nb+1:nmap,j))
   end do
end do

!write(*,*) 'products'
!write(*,fmt2) lld

llgb = 0d0
do i = 1, nmap
   do j = 1, nmap
      do a = 1, ng
         do b = ng+1, ng+nb
            llgb(i,j) = llgb(i,j) + lambda(a,i)*lambda(b,j)*sgb(a,b-ng)
         end do
      end do
   end do
end do

!write(*,*) 'products-overlaps gb'
!write(*,fmt2) llgb

llbg = 0d0
do i = 1, nmap
   do j = 1, nmap
      do a = ng+1, ng+nb
         do b = 1, ng
            llbg(i,j) = llbg(i,j) + lambda(a,i)*lambda(b,j)*sbg(a-ng,b)
         end do
      end do
   end do
end do

!write(*,*) 'products-overlaps gb'
!write(*,fmt2) llbg

llbd = 0d0
do i = 1, nmap
   do j = 1, nmap
      do a = ng+1, ng+nb
         do b = ng+nb+1, nmap
            llbd(i,j) = llbd(i,j) + lambda(a,i)*lambda(b,j)*sbd(a-ng,b-ng-nb)
         end do
      end do
   end do
end do

!write(*,*) 'products-overlaps bd'
!write(*,fmt2) llbd

lldb = 0d0
do i = 1, nmap
   do j = 1, nmap
      do a = ng+nb+1, nmap
         do b = ng+1, ng+nb
            lldb(i,j) = lldb(i,j) + lambda(a,i)*lambda(b,j)*sdb(a-ng-nb,b-ng)
         end do
      end do
   end do
end do

!write(*,*) 'products-overlaps db'
!write(*,fmt2) lldb

cnt = 1

g(1) = p_i
g(2) = p_j
g(3) = p_k

MonteCarlo: do mcs = 1, nmcs
   call sampling_class(bath,beta,kosc,c2,x,p)

   call sampling_mapng(init,rm,pm)
   
   call get_coeff(ng,beta,vomega,rm,pm,coeff)
!   coeff = (rm(1)**2 + pm(1)**2 - 0.5d0)
   call get_fact(nmap,llgb,llbg,rm,pm,fact)
   fact = fact*coeff*mu
   
   ib = 1
   pol(ib,cnt) = pol(ib,cnt) + fact

   call get_a(c2,ome,x,a1,a2)
   call get_force(nmap,ng,nb,lld,kosc,x,c2,rm,pm,f)

   MolecularDynamics: do it = 1, nmds
      call get_pulsefield(np,tau,it,dt,time,g,E0,E1,omega,et)
      
      call get_hm2(nmap,mu,et,a1,a2,hs,llgb,llbg,lld,hm)
      !call make_hm_traceless(nmap,hm)
      !write(*,*) 'hm'
      !write(*,fmt2) dble(hm)

      call update_p(dt2,f,p)

      call update_pm(dt2,hm,rm,pm)
      
      call update_x(dt,p,x)
      
      call update_a2(c2,x,a2)

      call get_hm2(nmap,mu,et,a1,a2,hs,llgb,llbg,lld,hm)
      !call make_hm_traceless(nmap,hm)

      call update_rm(dt,hm,pm,rm)

      call update_pm(dt2,hm,rm,pm)

      call get_force(nmap,ng,nb,lld,kosc,x,c2,rm,pm,f)

      call update_p(dt2,f,p)

      ib = it + 1
      call get_fact(nmap,llgb,llbg,rm,pm,fact)
      fact = fact*coeff*mu
      pol(ib,cnt) = pol(ib,cnt) + fact
   end do MolecularDynamics
end do MonteCarlo

dnmcs = dble(nmcs)
open(333,file="polariz.out")
do ib = 1, nmds + 1
!   pol_tot(ib) = -pol(ib,1) + pol(ib,2) + pol(ib,3) - pol(ib,4) + pol(ib,5)
!   pol_tot(ib) = pol_tot(ib) - pol(ib,6) - pol(ib,7) + pol(ib,8)
   pol_tot(ib) = pol(ib,1)/dnmcs
   write(333,*) time(3), ib-1, dble(pol_tot(ib)), aimag(pol_tot(ib))
end do

deallocate(c2)
deallocate(kosc)
deallocate(ome)
deallocate(tau)
deallocate(omega)
deallocate(time)
deallocate(pol)
deallocate(pol_tot)
deallocate(x)
deallocate(p)
deallocate(g)
deallocate(rm)
deallocate(pm)
deallocate(hm)

contains

subroutine iniconc()
implicit none

integer :: i

open(666,file="map.in")

read(666,*)
read(666,*) np,delta,kondo,nosc,ome_max
read(666,*)
read(666,*) nmcs,nmds,seed_dimension,dt,lumda_d
read(666,*)
read(666,*) eg,eb,ed,mu,e0,e1,beta,vomega
read(666,*)
read(666,*) tau1,omega1,time3,step1,step2
read(666,*)
read(666,*) bath,init,nfile
read(666,*)
read(666,*) basispc,ng,nb,nd
read(666,*)
read(666,*) p_i, p_j, p_k
close(666)

call random_seed(size=seed_dimension)
allocate (seed(seed_dimension))
do i = 1, seed_dimension
   seed(i) = 3*2**i - 1
end do
call random_seed(put=seed)
deallocate(seed)
end subroutine iniconc
end program modeliiimain
