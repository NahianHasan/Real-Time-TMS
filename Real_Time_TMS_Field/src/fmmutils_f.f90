subroutine computeEprimary(rs,js,robs,Eprimary,ntarget,nsource,precEP)
  implicit none
  integer, parameter :: dp=kind(1.0d0)
  integer, parameter :: spint=kind(8)
  real(kind=dp),parameter :: pi=3.141592653589793_dp,epsthresh=1d-12
  real(kind=dp),parameter ::eps0=8.854187820000000d-12
  real(kind=dp) ::  mu0over4pi=-1.000000000000000d-7
real(kind=dp),dimension(3,nsource),intent(in) :: rs
real(kind=dp),dimension(3,ntarget),intent(in) :: robs
real(kind=dp),dimension(3,nsource),intent(in) :: js
real(kind=dp),dimension(3,ntarget),intent(out) :: Eprimary
integer,intent(in) :: nsource,ntarget
real(kind=dp),intent(in) :: precEP

integer(kind=spint) :: i,j,ier


call lfmm3d_t_c_p_vec(3,precEP,nsource,rs,js,ntarget,robs,Eprimary)
Eprimary=mu0over4pi*Eprimary



end subroutine


subroutine computeHprimary(rs,js,robs,Hprimary,ntarget,nsource,precEP)
  implicit none
  integer, parameter :: dp=kind(1.0d0)
  integer, parameter :: spint=kind(8)
  real(kind=dp),parameter :: pi=3.141592653589793_dp,epsthresh=1d-12
  real(kind=dp),parameter ::eps0=8.854187820000000d-12
  real(kind=dp) ::  mu0over4pi=-1.000000000000000d-7
real(kind=dp),dimension(3,nsource),intent(in) :: rs
real(kind=dp),dimension(3,ntarget),intent(in) :: robs
real(kind=dp),dimension(3,nsource),intent(in) :: js
real(kind=dp),dimension(3,ntarget),intent(out) :: Hprimary
integer,intent(in) :: nsource,ntarget
real(kind=dp),intent(in) :: precEP

integer(kind=spint) :: i,j,ier
real(kind=dp),dimension(:),allocatable :: pottarg,tempcharge
real(kind=dp),dimension(:,:),allocatable :: fldtarg1



allocate(tempcharge(nsource))
allocate(fldtarg1(3,ntarget))
allocate(pottarg(ntarget))


tempcharge=reshape(mu0over4pi*js(1,:),(/nsource/))
call lfmm3d_t_c_g(precEP,nsource,rs,tempcharge,ntarget,robs,pottarg,fldtarg1)


	Hprimary(2,:)=-reshape(fldtarg1(3,:),(/ntarget/))
	Hprimary(3,:)=reshape(fldtarg1(2,:),(/ntarget/))

tempcharge=reshape(mu0over4pi*js(2,:),(/nsource/))
call lfmm3d_t_c_g(precEP,nsource,rs,tempcharge,ntarget,robs,pottarg,fldtarg1)

	Hprimary(1,:)=reshape(fldtarg1(3,:),(/ntarget/))
	Hprimary(3,:)=Hprimary(3,:)-reshape(fldtarg1(1,:),(/ntarget/))

tempcharge=reshape(mu0over4pi*js(3,:),(/nsource/))
call lfmm3d_t_c_g(precEP,nsource,rs,tempcharge,ntarget,robs,pottarg,fldtarg1)

	Hprimary(1,:)=Hprimary(1,:)-reshape(fldtarg1(2,:),(/ntarget/))
	Hprimary(2,:)=Hprimary(2,:)+reshape(fldtarg1(1,:),(/ntarget/))


deallocate(fldtarg1)
deallocate(pottarg)
deallocate(tempcharge)
end subroutine


subroutine computeEphiprimary(rs,rho,robs,Ephiprimary,ntarget,nsource,precEP)
  implicit none
  integer, parameter :: dp=kind(1.0d0)
  integer, parameter :: spint=kind(8)
  real(kind=dp),parameter :: pi=3.141592653589793_dp,epsthresh=1d-12
  real(kind=dp),parameter ::eps0=8.854187820000000d-12
  real(kind=dp) ::  mu0over4pi=-1.000000000000000d-7
real(kind=dp),dimension(3,nsource),intent(in) :: rs
real(kind=dp),dimension(3,ntarget),intent(in) :: robs
real(kind=dp),dimension(nsource),intent(in) :: rho
real(kind=dp),dimension(3,ntarget),intent(out) :: Ephiprimary
integer,intent(in) :: nsource,ntarget
real(kind=dp),intent(in) :: precEP

integer(kind=spint) :: i,j,ier
real(kind=dp),dimension(:),allocatable :: pottarg,tempcharge
real(kind=dp),dimension(:,:),allocatable :: fldtarg1



allocate(pottarg(ntarget))
call lfmm3d_t_c_g(precEP,nsource,rs,rho,ntarget,robs,pottarg,Ephiprimary)
Ephiprimary=-Ephiprimary/(4.0_dp*pi)
deallocate(pottarg)


end subroutine



subroutine computeAprimary(rs,js,robs,Aprimary,curlAprimary,ntarget,nsource,precEP)
  implicit none
  integer, parameter :: dp=kind(1.0d0)
  integer, parameter :: spint=kind(8)
  real(kind=dp),parameter :: pi=3.141592653589793_dp,epsthresh=1d-12
  real(kind=dp),parameter ::eps0=8.854187820000000d-12
  real(kind=dp) ::  mu0over4pi=1.000000000000000d-7
real(kind=dp),dimension(3,nsource),intent(in) :: rs
real(kind=dp),dimension(3,ntarget),intent(in) :: robs
real(kind=dp),dimension(3,nsource),intent(in) :: js
real(kind=dp),dimension(3,ntarget),intent(out) :: Aprimary,curlAprimary
integer,intent(in) :: nsource,ntarget
real(kind=dp),intent(in) :: precEP

integer(kind=spint) :: i,j,ier
real(kind=dp),dimension(:),allocatable :: pottarg,tempcharge
real(kind=dp),dimension(:,:),allocatable :: fldtarg1



allocate(tempcharge(nsource))
allocate(fldtarg1(3,ntarget))
allocate(pottarg(ntarget))



tempcharge=reshape(mu0over4pi*js(1,:),(/nsource/))
call lfmm3d_t_c_g(precEP,nsource,rs,tempcharge,ntarget,robs,pottarg,fldtarg1)
 Aprimary(1,:)=reshape(pottarg(:),(/ntarget/))
	curlAprimary(2,:)=-fldtarg1(3,:)
	curlAprimary(3,:)=fldtarg1(2,:)

tempcharge=reshape(mu0over4pi*js(2,:),(/nsource/))
call lfmm3d_t_c_g(precEP,nsource,rs,tempcharge,ntarget,robs,pottarg,fldtarg1)
 Aprimary(2,:)=reshape(pottarg(:),(/ntarget/))
	curlAprimary(1,:)=fldtarg1(3,:)
	curlAprimary(3,:)=curlAprimary(3,:)-fldtarg1(1,:)

tempcharge=reshape(mu0over4pi*js(3,:),(/nsource/))
call lfmm3d_t_c_g(precEP,nsource,rs,tempcharge,ntarget,robs,pottarg,fldtarg1)
Aprimary(3,:)=reshape(pottarg(:),(/ntarget/))
	curlAprimary(1,:)=curlAprimary(1,:)-fldtarg1(2,:)
	curlAprimary(2,:)=curlAprimary(2,:)+fldtarg1(1,:)


deallocate(fldtarg1)
deallocate(pottarg)
deallocate(tempcharge)



end subroutine
