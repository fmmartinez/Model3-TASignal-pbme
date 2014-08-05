program modeliiimain
use m_map, only: iniconq_d,get_preh,sampling_class,sampling_mapng,get_coeff,  &
                  get_fact,get_a,get_force_traceless,get_pulsefield,get_hm2,  &
                  make_hm_traceless,update_p,update_x,update_pm,update_rm,    &
                  update_a2
implicit none

real(8),parameter :: pi=3.1415926535d0, twopi = 2d0*pi

character(len=2) :: c_ng,c_nt
character(len=9) :: fmt1,fmt2

complex(8) :: coeff,fact,a1,a2,et
complex(8),dimension(:),allocatable :: pol_tot,x,p,rm,pm,f
complex(8),dimension(:,:),allocatable :: pol,hm

integer :: a,b,i,j,is,it,cnt,p_i,p_j,p_k,ib,nmap,ng,nb,nd,basispc
integer :: np,nmcs,mcs,nmds,seed_dimension,nosc,step1,bath,init,nfile
integer,dimension(:),allocatable :: seed1,g

real(8) :: gauss,dt,dt2,kondo,delta,beta,ome_max,lumda_d,eg,eb,ed,mu,e0,e1,sij,vomega
real(8) :: step2,dnmcs,tau1,omega1,tau2,omega2,time3,lambdacheck
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

call iniconq_d(nosc,lumda_d,ome_max,ome,c2,kosc)

dt  = twopi*dt
dt2 = 0.5d0*dt

do i = 1, np-1
   tau(i)   = tau1
   omega(i) = omega1
end do
tau(np) = tau2
omega(np) = omega2

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

!call get_lambda_eigenvectors(ng,nb,nd,eg,eb,ed,delta,vomega, &
!                              sgg,sgb,sgd,sbg,sbb,sbd,sdg,sdb,sdd,lambda,hs)

call get_preh(ng,nb,nd,eg,eb,ed,delta,vomega,hs)


!declare variables for faster calculations
llg = 0d0
do i = 1, ng
   llg(i,i) = 1d0
end do
llb = 0d0
do i = ng+1, ng+nb
   llb(i,i) = 1d0
end do
lld = 0d0
do i = ng+nb+1,nmap
   lld(i,i) = 1d0
end do

llgb = 0d0
llgb(1:ng,ng+1:ng+nb) = hs(1:ng,ng+1:ng+nb)

llbg = 0d0
llbg(ng+1:ng+nb,1:ng) = hs(ng+1:ng+nb,1:ng)

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
   call get_force_traceless(nmap,ng,nb,lld,kosc,x,c2,rm,pm,f)

   MolecularDynamics: do it = 1, nmds
      call get_pulsefield(np,tau,it,dt,time,g,E0,E1,omega,et)
      
      call get_hm2(nmap,mu,et,a1,a2,hs,llgb,llbg,lld,hm)
      call make_hm_traceless(nmap,hm)
      !write(*,*) 'hm'
      !write(*,fmt2) dble(hm)

      call update_p(dt2,f,p)

      call update_pm(dt2,hm,rm,pm)
      
      call update_x(dt,p,x)
      
      call update_a2(c2,x,a2)

      call get_hm2(nmap,mu,et,a1,a2,hs,llgb,llbg,lld,hm)
      call make_hm_traceless(nmap,hm)

      call update_rm(dt,hm,pm,rm)

      call update_pm(dt2,hm,rm,pm)

      call get_force_traceless(nmap,ng,nb,lld,kosc,x,c2,rm,pm,f)

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
read(666,*) tau1,omega1,tau2,omega2,time3,step1,step2
read(666,*)
read(666,*) bath,init,nfile
read(666,*)
read(666,*) ng,nb,nd
read(666,*)
read(666,*) p_i, p_j, p_k
close(666)

call random_seed(size=seed_dimension)
allocate (seed1(seed_dimension))
do i = 1, seed_dimension
   seed1(i) = 3*2**i - 1
end do
call random_seed(put=seed1)
deallocate(seed1)
end subroutine iniconc

end program modeliiimain
