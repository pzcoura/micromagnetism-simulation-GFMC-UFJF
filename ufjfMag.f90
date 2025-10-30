    Module VARIAVEIS
!    
    real(kind=8), Allocatable,dimension (:) :: uxi, uyi,uzi
    real(kind=8), Allocatable,dimension (:):: uxi0, uyi0,uzi0    
    real(kind=8), Allocatable,dimension (:) :: rvij3,uvrijx,uvrijy,uvrijz   
    real(kind=8), Allocatable,dimension (:) :: campo_efetivo_x,campo_efetivo_y,campo_efetivo_z 
    real(kind=8), Allocatable,dimension (:) :: mdx,mdy,mdz
    real(kind=8), Allocatable,dimension (:) :: mdx0,mdy0,mdz0
    real(kind=8), Allocatable,dimension (:) :: md1x,md1y,md1z
    real(kind=8), Allocatable,dimension (:) :: md2x,md2y,md2z
    real(kind=8), Allocatable,dimension (:) :: md3x,md3y,md3z
    real(kind=8), Allocatable,dimension (:) :: md0x,md0y,md0z
    real(kind=8), Allocatable,dimension (:) :: k1x,k1y,k1z
    real(kind=8), Allocatable,dimension (:) :: k2x,k2y,k2z
    real(kind=8), Allocatable,dimension (:) :: k3x,k3y,k3z
    real(kind=8), Allocatable,dimension (:) :: k4x,k4y,k4z
    real(kind=8), Allocatable,dimension (:) :: jotaty,dmty,anisoty,magty
    real(kind=8), Allocatable,dimension (:) :: vecnormalx,vecnormaly,vecnormalz
    real(kind=8), Allocatable,dimension (:) :: vectangxxd,vectangyxd,vectangzxd
    real(kind=8), Allocatable,dimension (:) :: vectangxxe,vectangyxe,vectangzxe
    real(kind=8), Allocatable,dimension (:) :: vectangxyd,vectangyyd,vectangzyd
    real(kind=8), Allocatable,dimension (:) :: vectangxye,vectangyye,vectangzye
    real(kind=8), Allocatable,dimension (:) :: vectangxzd,vectangyzd,vectangzzd
    real(kind=8), Allocatable,dimension (:) :: vectangxze,vectangyze,vectangzze
    real(kind=8), Allocatable,dimension (:) :: modtgxd,modtgxe,modtgyd,modtgye,modtgzd,modtgze,volcel
    real(kind=8), Allocatable,dimension (:) :: jotacelxd,jotacelxe,jotacelyd,jotacelye,jotacelzd,jotacelze,magsavol
    real(kind=8), Allocatable,dimension (:) :: dmcelxd,dmcelxe,dmcelyd,dmcelye,dmcelzd,dmcelze,anisovol,dist_x
    real(kind=8), Allocatable,dimension (:) :: tangxxd,tangyxd,tangzxd
    real(kind=8), Allocatable,dimension (:) :: tangxxe,tangyxe,tangzxe
    real(kind=8), Allocatable,dimension (:) :: tangxyd,tangyyd,tangzyd
    real(kind=8), Allocatable,dimension (:) :: tangxye,tangyye,tangzye
    real(kind=8), Allocatable,dimension (:) :: tangxzd,tangyzd,tangzzd
    real(kind=8), Allocatable,dimension (:) :: tangxze,tangyze,tangzze    
    real(kind=8), Allocatable,dimension (:) :: dvecnx,dvecny,dvecnz
    real(kind=8), Allocatable,dimension (:) :: pgx,pgy,pgz,grho1,incli
    real(kind=8), Allocatable,dimension (:) :: exchangx,exchangy,exchangz,anisotroz,anisotrox,anisotroy 
    real(kind=8), Allocatable,dimension (:) :: dzymoriax,dzymoriay,dzymoriaz
    real(kind=8), Allocatable,dimension (:) :: dipx,dipy,dipz,dipvol
    integer, Allocatable,dimension (:) :: usa,nvtp
    integer, Allocatable,dimension (:) :: gx,gy,gz,gvx,gvy,gvz,gvv
    integer, Allocatable,dimension (:) :: gvxd1,gvxe1,gvyd1,gvye1,gvzd1,gvze1
    integer, Allocatable,dimension (:) :: ty
    integer, Allocatable,dimension (:) :: sinal,lj,nli,nlf,amp
    real(kind=8) :: raio, exps,lim_raio,lim_raio1
    real(kind=8) :: energiatot,energdip,dzym,troca,anisotro
    real(kind=8) :: camp_bx,camp_by,camp_bz,camp_bxi,camp_byi,camp_bzi,energcamp
    real(kind=8) :: ap,mags,jota,dip,aniso,dm,gama_alfa,gama_alfa1,alfa
    real(kind=8) :: magty_1,jotaty_1,jotaty_2,dmty_1,dmty_2,anisoty_1
    real(kind=8) :: aniso1,aniso_ef,aniso_efy
    real(kind=8) :: magsa3,dt,dt4,pol,je,adiaba,vj,x0,y0
    real(kind=8) :: vj_ap,vj_alfa3,vj_qui_alfa2,vj_qui1,ct_sot,teta_SHE
    real(kind=8) :: raio_imp,raio_imp1,pxic,pyjc
    real(kind=8) :: beta,passo,gama_tip
    real(kind=8) :: uzipe,xcve,ycve,zcve,uzipd,xcvd,xcvm,ycvm,ycvd,zcvm,zcvd,uziped
    real(kind=8) :: xcvei,ycvej,xcvdi,ycvdj,xcvmi,ycvmj
    real(kind=8) :: nlx1,nlx2,nlx3,nlx4,nlx5,nlx6
    real(kind=8) :: pgxmax1,pgxmax2,pgxmax3,pgxmax4,pgxmax5,pgxmax6,pgzmax1,beta1
    integer(kind=8) :: xmax,ymax,zmax,npt,npt1,npt2,flag_cut,cut,dy,ncx,ncy
    integer :: flag_cor,num_it,npequi,ntpassos
    integer :: nj,it,itv,tipo_dm,adams,polar_sot,tipo
    integer :: impu,concav
    integer :: forma_imp,pxi,pxf,pyi,pyf,nano,m_tip,direc,direc_tipo,nl1,nl2,nl12,nl122
    real(kind=8), PARAMETER :: pi = 3.14159265359, umfento=1D+15, mi_0=4*pi*1D-07,gama=1.76D+11
    real(kind=8), PARAMETER :: qe=1.60217662D-19, magneton=9.27400899D-24  
    real(kind=8), PARAMETER :: cte_h=6.62607004D-34
!
!   qe=1.60217662D-19  ! carga do eletron (Unit: Coulomb)
!   magneton=9.27400899D-24  ! magneton de bohr (Unit: J/T - T => Tesla)
!   mi_0=4*pi*1D-7  ! permeabilidade magnética no vácuo (Unit: Newton/Ampere^2)
!   gama=1,76D+11   !  razão giromagnética do elétron (Unit: 1/(T.S) - T => Tesla - S => segundo)
!   6,62607004 × 10-34  ! constante de planck (unit: m^2 kg / s)
!
    End Module
!
                        program skyrmions
!
    use VARIAVEIS
    real :: start_time, end_time,je1
!
    it=0
!
    call ler_dados
! 
    call configuracao_fita
!
    call tabela
!
    if(flag_cut.eq.0) call energia_campo_local
    if(flag_cut.eq.1) call energia_campo_local_corte
!
    call escreve(it)
!
    do nj=1,3
!
    do itv=0,3
!
    call runge_kuta(nj,itv)
!
    if(flag_cut.eq.0) call energia_campo_local
    if(flag_cut.eq.1) call energia_campo_local_corte
!
    end do
!
    end do
!
!!!!! configuração de equilíbrio !!!!!!!!!!!!!!!!!!!
!
!
    do it=1,npequi
!
    if(mod(it,num_it).eq.0) then
!
    write(*,*)it
    call escreve(it)
!    
    end if
!
    if(adams.eq.2) then
!
    do nj=1,3
!
    do itv=0,3
!
    call runge_kuta(nj,itv)
!
    if(flag_cut.eq.0) call energia_campo_local
    if(flag_cut.eq.1) call energia_campo_local_corte
!
    end do
!
    end do
!
    else
!
    call adams_bashforth
!
    if(flag_cut.eq.0) call energia_campo_local
    if(flag_cut.eq.1) call energia_campo_local_corte
!
    end if
!
    end do
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    call cpu_time( start_time )
!
    do it=npequi+1,ntpassos
!    
    if(mod(it,num_it).eq.0) then
! 
    write(*,*)it
    call escreve(it)
!
    end if
!
    if(adams.eq.2) then
!
    do itv=0,3
!
    if(flag_cor.eq.0) call runge_kuta_flag_cor0(nj,itv)
    if(flag_cor.eq.1) call runge_kuta_flag_cor1(nj,itv)
    if(flag_cor.eq.3) call runge_kuta(nj,itv)
!
    if(flag_cut.eq.0) call energia_campo_local
    if(flag_cut.eq.1) call energia_campo_local_corte
!
    end do
!
    else
!
    if(flag_cut.eq.0) call energia_campo_local
    if(flag_cut.eq.1) call energia_campo_local_corte
!
    if(flag_cor.eq.0) call adams_bashforth_flag_cor0
    if(flag_cor.eq.1) call adams_bashforth_flag_cor1
    if(flag_cor.eq.3) call adams_bashforth
!
    end if
!
    end do
!
    call cpu_time( end_time )
    write(*,*)'tempo de execução=',end_time,'segundos'
    write(*,*)'tempo de execução=',end_time/60d0,'minutos'
    write(*,*)'tempo de execução=',end_time/3600d0,'horas'
    end 
!
!
!
    subroutine ler_dados
    Use VARIAVEIS
    integer :: allocateStatus,ip,i,j,k
!
    open(20,file='dados_entrada_artigo')
    read(20,*)xmax,ymax,ncx,ncy,zmax,flag_cut,cut,flag_cor,tipo,x0,y0,raio,exps,lim_raio,lim_raio1
    read(20,*)ap,mags,jota,aniso,dm,tipo_dm,alfa,teta_SHE,polar_sot,dt,ntpassos,npequi,num_it
    read(20,*)pol,je,adiaba
    read(20,*)impu
    read(20,*)magty_1,jotaty_1,jotaty_2,dmty_1,dmty_2,anisoty_1,forma_imp
    read(20,*)pxi,pxf,pyi,pyf,pxic,pyjc,raio_imp,raio_imp1,nano,m_tip,gama_tip
    read(20,*)direc,direc_tipo,nl1,nl2,beta,passo,concav,camp_bxi,camp_byi,camp_bzi,adams
    close(20)
!
    if(nano.eq.1.or.nano.eq.3)dy=0
    if(nano.eq.2)then
    dy=1
    write(*,*)'Since the chosen option was nano=2, one layer was added to the ymax variable &
     due to the periodic boundary condition along the cylinder circumference.'
    end if
!
    npt=(2*xmax+ncx+1)*(2*ymax+ncy+dy+1)*zmax
    npt1=(2*xmax+ncx+3)*(2*ymax+ncy+dy+3)*zmax
    if(flag_cut.eq.0) npt2=npt*(npt-1)/2+1
    if(flag_cut.eq.1.and.zmax.eq.1) npt2=npt*4*(cut**2)
    if(flag_cut.eq.1.and.zmax.gt.1) npt2=zmax*4*npt*(cut**2)
!
!   Lattice parameter (ap - Unit: meters)
!   Saturation magnetization (mags - Unit: Ampere/meter)
!   Stiffness constant (jota - unit: Joule/meter)
!   Anisotropic interaction (aniso - Unit: Joule/m³)
!   Dzyaloshinskii-Moriya interaction (dm - Unit: Joule/m²)
!   Dipolar constant for the local dipolar field (dip - Not an input parameter - Unit: Tesla)
!   Damping constant (alfa - Unit: dimensionless)
!   Time step (dt - Unit: seconds)
!   Spin polarization current (pol - dimensionless)
!   Polarization current intensity (je - Unit: A/m²)
!   Adiabatic constant (adiaba - dimensionless)
!   Effective Spin Hall angle constant (teta_SHE - dimensionless)
!
!!!!!!!!   The constants below refer to those described above   !!!!!
!
    magsa3=mags*ap**3
    jota = 2d0*jota
    dm=dm
    aniso1=aniso
    aniso=2d0*aniso
!!!! efective field !!!!!!!!!!!!!!!!!!!!!!!!!
    If(tipo.eq.2) aniso1=2*aniso1
    aniso_ef=2d0*(aniso1-mi_0*((mags)**2)/2d0)
    aniso_efy=2d0*(aniso1*anisoty_1-mi_0*((magty_1*mags)**2)/2d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!  external field !!!!!!!!!!!!!!!!!!!!!!!!!!!
    camp_bx=camp_bxi
    camp_by=camp_byi
    camp_bz=camp_bzi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dip=mi_0*mags/(4d0*pi)
    gama_alfa=-gama/(1+alfa**2)
!
    dt4=dt/24d0
!
!!!!!!! Constants for the current term - STT !!!!!!!!!!!!!!!!!!
!
    vj=pol*je*magneton/(qe*mags*(1d0+adiaba**2))
    quisi=adiaba
    vj_ap=vj/ap
    vj_alfa3=vj_ap*alfa
    vj_qui_alfa2=(quisi*alfa+1d0)*vj_ap
    vj_qui1=vj_ap*quisi
    gama_alfa1=-1d0/(1d0+alfa**2)
!
!!!!!!! Constants for the current term - SOT  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! journal of Applied Physics 117, 17E505 (2015) - Kyung-Jin Lee !!!!!!!!!!!
!
    ct_sot=-gama*cte_h*teta_SHE*je/((1d0+alfa**2)*(4d0*pi*qe*ap*mags))
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    
    Allocate (uxi(0:npt1), uyi(0:npt1),uzi(0:npt1),STAT = AllocateStatus)
    Allocate (uxi0(0:npt1), uyi0(0:npt1),uzi0(0:npt1))    
    IF (AllocateStatus /= 0) STOP "*** Not enough memory uxi,uyi,uzi***"
    Allocate (rvij3(npt2),STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory rvij rvij3***" 
    Allocate (mdx(npt),mdy(npt),mdz(npt))
    Allocate (mdx0(npt),mdy0(npt),mdz0(npt))
    Allocate (md1x(npt),md1y(npt),md1z(npt))
    Allocate (md2x(npt),md2y(npt),md2z(npt))
    Allocate (md3x(npt),md3y(npt),md3z(npt))
    Allocate (md0x(npt),md0y(npt),md0z(npt))
    Allocate (k1x(npt),k1y(npt),k1z(npt))
    Allocate (k2x(npt),k2y(npt),k2z(npt))
    Allocate (k3x(npt),k3y(npt),k3z(npt))
    Allocate (k4x(npt),k4y(npt),k4z(npt))
    Allocate (campo_efetivo_x(npt),campo_efetivo_y(npt),campo_efetivo_z(npt))
    Allocate (uvrijx(npt2),uvrijy(npt2),uvrijz(npt2))
    Allocate (gvxd1(npt),gvxe1(npt),gvyd1(npt),gvye1(npt),gvzd1(npt),gvze1(npt))
    Allocate ( gx(npt),gy(npt),gz(npt),pgx(npt),pgy(npt),pgz(npt))
    Allocate (gvv(npt2))
    Allocate (usa(npt),nvtp(0:npt))
    Allocate (ty(0:npt1))
    Allocate (sinal(npt),lj(npt),nli(npt),nlf(npt))
    Allocate (jotaty(0:3),dmty(0:3),anisoty(0:3),magty(0:3))
    Allocate (vecnormalx(npt),vecnormaly(npt),vecnormalz(npt))
    Allocate (vectangxxd(0:npt1),vectangyxd(0:npt1),vectangzxd(0:npt1))
    Allocate (vectangxxe(0:npt1),vectangyxe(0:npt1),vectangzxe(0:npt1))
    Allocate (vectangxyd(0:npt1),vectangyyd(0:npt1),vectangzyd(0:npt1))
    Allocate (vectangxye(0:npt1),vectangyye(0:npt1),vectangzye(0:npt1))
    Allocate (vectangxzd(0:npt1),vectangyzd(0:npt1),vectangzzd(0:npt1))
    Allocate (vectangxze(0:npt1),vectangyze(0:npt1),vectangzze(0:npt1))
    Allocate (tangxxd(0:npt1),tangyxd(0:npt1),tangzxd(0:npt1))
    Allocate (tangxxe(0:npt1),tangyxe(0:npt1),tangzxe(0:npt1))
    Allocate (tangxyd(0:npt1),tangyyd(0:npt1),tangzyd(0:npt1))
    Allocate (tangxye(0:npt1),tangyye(0:npt1),tangzye(0:npt1))
    Allocate (tangxzd(0:npt1),tangyzd(0:npt1),tangzzd(0:npt1))
    Allocate (tangxze(0:npt1),tangyze(0:npt1),tangzze(0:npt1))
    Allocate (modtgxd(0:npt1),modtgxe(0:npt1),modtgyd(0:npt1),modtgye(0:npt1),modtgzd(0:npt1),modtgze(0:npt1),volcel(0:npt1))
    Allocate (magsavol(0:npt1),dist_x(0:npt1))
    Allocate (jotacelxd(0:npt1),jotacelxe(0:npt1),jotacelyd(0:npt1),jotacelye(0:npt1),jotacelzd(0:npt1),jotacelze(0:npt1))
    Allocate (dmcelxd(0:npt1),dmcelxe(0:npt1),dmcelyd(0:npt1),dmcelye(0:npt1),dmcelzd(0:npt1),dmcelze(0:npt1),anisovol(0:npt1))
    Allocate (dvecnx(0:npt1),dvecny(0:npt1),dvecnz(0:npt1))
    Allocate (exchangx(npt1),exchangy(npt1),exchangz(npt1),anisotroz(npt1),anisotrox(npt1),anisotroy(npt1))
    Allocate (dzymoriax(npt1),dzymoriay(npt1),dzymoriaz(npt1))
    Allocate (dipx(npt1),dipy(npt1),dipz(npt1),dipvol(npt1))
    Allocate (grho1(npt),incli(npt))
!
    return
    end
!
!
!
    subroutine configuracao_fita
    Use VARIAVEIS
    integer :: i,j,k,ip,nic,is,isx,isy,im,jdif,idif
    real(kind=8) :: r,cosphi,senophi,teta,senoteta,costeta,a,s,l,il
    real(kind=8) :: rho,rho1,raio_circun,ux,uy,uz,cosphi1,senophi1,abertura
    real(kind=8) :: isr,deltas,deltarx,deltary,deltarz,pgz_1
    real(kind=8), dimension (-1000:1000) :: pgxi,pgyi,pgzi,pgxa,pgya,pgza
    open(11,file='config.xyz')
    open(2,file='abertura')
    open(1, file='normal.xyz')
    write(1,*)npt
    write(1,*)'vetor normal'
!
!!!!! We identified the host with index 0 and the impurity with index 1 via the vector ty !!
!
      ty=0                  !!!! All are initialized with index 0
      magty(0)=1d0
      magty(1)=magty_1
!
!!!!! We have 3 types of possibilities:Host-host interaction             0+0=0
!!!!!                                  Host-impurity interaction         0+1=1
!!!!!                                  Impurity-impurity interaction     1+1=2
!
!!!! host-host interaction 0+0=0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      jotaty(0)=jota
      anisoty(0)=aniso
      if(cut.eq.-1.and.flag_cut.eq.1) anisoty(0)=aniso_ef
      vecnormalx=0d0
      vecnormaly=0d0
      vecnormalz=0d0
      dmty(0)=dm
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!! If there are no impurities, all have the host value. !!!!!!!!!!!!!!!!!!!
!
     if(impu.eq.0) then
!
      magty_1=1d0          !!!! ratio mag_impurity/mag_host
      jotaty_1=jota        !!!! ratio exchange-impurity-host/exchange-host
      jotaty(1)=jotaty_1
      dmty_1=dm            !!!! ratio DM_impurity-host/DM-host
      dmty(1)=dmty_1
      anisoty_1=aniso      !!!! ratio anisotropy-impurity/anisotropy-host
      anisoty(1)=anisoty_1
      jotaty_2=jota        !!!! ratio exchange_impurity-impurity/exchange-host
      dmty_2=dm            !!!! ratio DM_impurity-impurity/DM-host
      jotaty(2)=jotaty_2
      dmty(2)=dmty_2
!
    end if
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
     if(impu.eq.1) then
!
!!!! Host-impurity interaction 0+1=1  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      magty_1=magty_1
      jotaty_1=jotaty_1*jota
      jotaty(1)=jotaty_1
      dmty_1=dmty_1*dm
      dmty(1)=dmty_1
      if(magty_1.ne.0) anisoty(1)=anisoty_1*aniso/magty_1
      if(magty_1.eq.0) anisoty(1)=0d0

      if(cut.eq.-1.and.flag_cut.eq.1) then

      if(magty_1.ne.0) anisoty(1)=aniso_efy/magty_1
      if(magty_1.eq.0) anisoty(1)=0d0

      end if

!
!!!! Impurity-impurity interaction 1+1=2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      jotaty_2=jotaty_2*jota
      dmty_2=dmty_2*dm
      jotaty(2)=jotaty_2
      dmty(2)=dmty_2
!
      end if
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
	a=raio
	s=exps
!	
    il=x0
    l=y0
    ip=0
!
    uxi=0d0;uyi=0d0;uzi=0d0
    ux=0d0;uy=0d0;uz=0d0
    grho1=0d0
    raio_circun=beta
    beta=pi*beta/180d0
    nl12=nl2-nl1
    nl122=nl12/2
    pgxmax=0;pgymax=0;pgzmax=0
!
    nic=0
!
	do i=-xmax,xmax+ncx
!
 		do j=-ymax,ymax+ncy+dy
!
            do k=1,zmax
!
            ip=ip+1
            gx(ip)=i
            gy(ip)=j
            gz(ip)=k
            pgx(ip)=i
            pgy(ip)=j
            pgz(ip)=k
            vecnormalx(ip)=0d0
            vecnormaly(ip)=0d0
            vecnormalz(ip)=1d0
            grho1(ip)=0d0
!
    if(nano.eq.1) then
!
!!!!!Defect in the X direction !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if(direc.eq.1) then
!
    if(direc_tipo.eq.4) then
!
!              rho=beta/(nl12-1)
!              raio_circun=0.5d0/(1d0-(cos((0.5d0*rho)))**2)**0.5
!              abertura=beta
             
             rho=2d0*acos((1d0-1/(2*raio_circun)**2)**0.5d0)
             abertura=rho*(nl12)
            
             
             if(ip.eq.1) then
             write(2,*)raio_circun,abertura*180/pi
             flush(2)
             end if
             if(ip.eq.1) write(*,*)'Raio=',raio_circun,'abertura=',abertura*180/pi
             if(abertura*180/pi.gt.180) stop
!            
            if(i.le.nl1) then
            pgx(ip)=i
            pgy(ip)=j
            pgz(ip)=k
            pgxi(i)=i
            pgzi(i)=k
            vecnormalx(ip)=0d0
            vecnormaly(ip)=0d0
            vecnormalz(ip)=1d0
            dvecnx(i)=vecnormalx(ip)
            dvecny(i)=vecnormaly(ip)
            dvecnz(i)=vecnormalz(ip)
            grho1(ip)=0d0
            incli(ip)=0d0
            end if
!            
            if(i.gt.nl1.and.i.le.(nl2)) then
            rho1=(i-nl1)*rho
            pgz(ip)=k+concav*raio_circun*(1d0-cos(rho1))
            pgzi(i)=k+concav*raio_circun*(1d0-cos(rho1))
            pgx(ip)=nl1+raio_circun*sin(rho1)
            pgxi(i)=nl1+raio_circun*sin(rho1)
            pgy(ip)=j
            vecnormaly(ip)=0d0
            vecnormalx(ip)=-concav*sin(rho1)
            vecnormalz(ip)=cos(rho1)
            dvecnx(i)=vecnormalx(ip)
            dvecny(i)=vecnormaly(ip)
            dvecnz(i)=vecnormalz(ip)
             grho1(ip)=-concav*(pi/2d0-asin(vecnormalz(ip)))
             incli(ip)=-concav*atan((pgzi(i)-pgzi(i-1))/(pgxi(i)-pgxi(i-1)))
             end if
!      
            if(i.gt.nl2) then
            pgz(ip)=k+concav*raio_circun*(1d0-cos(rho1))+sin(concav*abertura)*(i-nl2)
            pgx(ip)=pgxi(nl2)+cos(concav*abertura)*(i-nl2)
            pgy(ip)=j
            vecnormaly(ip)=0d0
            vecnormalz(ip)=cos(concav*abertura)
            vecnormalx(ip)=-sin(concav*abertura)
            dvecnx(i)=vecnormalx(ip)
            dvecny(i)=vecnormaly(ip)
            dvecnz(i)=vecnormalz(ip)
            grho1(ip)=-concav*(pi/2d0-asin(vecnormalz(ip)))
            incli(ip)=0d0
            end if
!            
    end if
!
    if(direc_tipo.eq.5) then
!
            if(i.eq.-xmax) then
            idif=i+xmax
            pgxi(idif)=0
!
            pgzi(idif)=-concav*sqrt((pgxi(idif)**2/(10.61**2)+1)*(10.61)**2)+k ! Function to generate a curved surface in the X direction: z=f(x)
!
            deltarx=1
            deltarz=0d0
            gx(ip)=idif
            pgz(ip)=pgzi(idif)
            pgx(ip)=pgxi(idif)
            pgy(ip)=j
            vecnormaly(ip)=0d0
            vecnormalz(ip)=-1d0
            vecnormalx(ip)=0d0
            grho1(ip)=asin(vecnormalx(ip))
            end if
!
            if(i.gt.-xmax) then
            if(i.le.0) idif=i+xmax
            if(i.gt.0) idif=-i
            do is=1,100000
            isr=float(is)/100000
            if(i.le.0) pgxi(idif)=pgxi(idif-1)+isr
            if(i.gt.0) pgxi(idif)=pgxi(idif+1)-isr
!
!!!!!!!!!!! The defect consists of a flat section in the intervals [-xmax, nl1] and [nl2, xmax], and a curved surface
!!!!!!!!!!! in the interval [nl1, nl2]. All measurements are relative to the origin (center of the surface).
!!!!!!!!!!! If nl1 < -xmax and nl2 > xmax, the entire surface in the interval [-xmax, xmax] becomes curved.
!!!!!!!!!!! The surface along the y-direction remains unchanged.
!
            if(idif.gt.0.and.idif.le.nl2) pgzi(idif)=-concav*sqrt((pgxi(idif)**2/(10.61**2)+1)*(10.61)**2)+k  ! Function to generate a curved surface in the X direction: z=f(x)
            if(idif.ge.nl2) pgzi(idif)= pgzi(idif-1)
            if(idif.lt.0.and.idif.gt.nl1) pgzi(idif)=-concav*sqrt((pgxi(idif)**2/(10.61**2)+1)*(10.61)**2)+k  ! Function to generate a curved surface in the X direction: z=f(x)
            if(idif.le.nl1) pgzi(idif)= pgzi(idif+1)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(i.le.0) then
            deltarx=(pgxi(idif)-pgxi(idif-1))
            deltarz=-(pgzi(idif)-pgzi(idif-1))
            end if
            if(i.gt.0) then
            deltarx=-(pgxi(idif)-pgxi(idif+1))
            deltarz=(pgzi(idif)-pgzi(idif+1))
            end if
            r=sqrt(deltarx**2+deltarz**2)
            deltas=abs(r-1d0)
            if(deltas.le.0.0001d0) then
            gx(ip)=idif
            pgz(ip)=pgzi(idif)
            pgx(ip)=pgxi(idif)
            pgy(ip)=j
            vecnormaly(ip)=0d0
            vecnormalz(ip)=-deltarx
            vecnormalx(ip)=-deltarz
            grho1(ip)=asin(deltarz)
            go to 111
            end if
            end do
            end if
!
       end if
! 
!    
 111   end if
!
            if(direc_tipo.eq.6) then
!
            rho= beta*abs(float(i+xmax)) /passo
!
            pgx(ip)=float(i)
!
!!!!!!!!! helix around the x-axis !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            pgz(ip)=float(j)*sin(rho)+(k-1)*cos(rho)
            pgy(ip)=float(j)*cos(rho)-(k-1)*sin(rho)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            vecnormalx(ip)=0d0
            vecnormaly(ip)=sin(rho)
            vecnormalz(ip)=cos(rho)
            grho1(ip)=rho
!
            end if



!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!Defect in the Y direction  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if(direc.eq.2) then
!
    if(direc_tipo.eq.4) then
!
!              rho=beta/(nl12-1)
!              raio_circun=0.5d0/(1d0-(cos((0.5d0*rho)))**2)**0.5
!              abertura=beta
             rho=2d0*acos((1d0-1/(2*raio_circun)**2)**0.5d0)
             abertura=rho*(nl12-1)
             
             if(ip.eq.1) then 
             write(2,*)raio_circun,abertura*180/pi
             flush(2)
             end if
             if(ip.eq.1) write(*,*)raio_circun,abertura*180/pi
             if(abertura*180/pi.gt.180) stop
!            
!            
            if(j.le.nl1) then
            pgx(ip)=i
            pgy(ip)=j
            pgz(ip)=k
            pgyi(j)=j
            pgzi(j)=k
            vecnormalx(ip)=0d0
            vecnormaly(ip)=0d0
            vecnormalz(ip)=1d0
            grho1(ip)=0d0
            end if
!            
            if(j.gt.nl1.and.j.lt.(nl2)) then
            rho1=(j-nl1)*rho
            pgz(ip)=k+concav*raio_circun*(1d0-cos(rho1))
            pgzi(j)=k+concav*raio_circun*(1d0-cos(rho1))
            pgy(ip)=nl1+raio_circun*sin(rho1)
            pgyi(j)=nl1+raio_circun*sin(rho1)
            pgx(ip)=i
            vecnormalx(ip)=0d0
            vecnormaly(ip)=-concav*sin(rho1)
            vecnormalz(ip)=cos(rho1)
            grho1(ip)=-concav*(pi/2d0-asin(vecnormalz(ip)))
            end if
!      
            if(j.ge.nl2) then
            pgz(ip)=k+concav*raio_circun*(1d0-cos(rho1))+sin(concav*abertura)*(j-nl2+1)
            pgy(ip)=pgyi(nl2-1)+cos(concav*abertura)*(j-nl2+1) 
            pgx(ip)=i
            vecnormalx(ip)=0d0
            vecnormalz(ip)=cos(concav*abertura)
            vecnormaly(ip)=-sin(concav*abertura)
            grho1(ip)=-concav*(pi/2d0-asin(vecnormalz(ip)))
            end if
!            
    end if
!
    if(direc_tipo.eq.5) then
!
            if(j.eq.-ymax) then
            jdif=j+ymax
            pgyi(jdif)=0d0
            pgzi(idif)=-concav*sqrt(pgyi(jdif)**2+(10.61)**2)+k   ! Function to generate a curved surface in the Y direction: z=f(y)
!             pgzi(jdif)=k
            gy(ip)=jdif
            pgz(ip)=pgzi(jdif)
            pgy(ip)=pgyi(jdif)
            pgx(ip)=i
            vecnormalx(ip)=0d0
            vecnormalz(ip)=1d0
            vecnormaly(ip)=0
            grho1(ip)=asin(vecnormaly(ip))
            end if
!
            if(j.gt.-ymax) then
            if(j.le.0) jdif=j+ymax
            if(j.gt.0) jdif=-j
            do is=1,100000
            isr=float(is)/100000
            if(j.le.0) pgyi(jdif)=pgyi(jdif-1)+isr
            if(j.gt.0) pgyi(jdif)=pgyi(jdif+1)-isr
!
!!!!!!!!!!! The defect consists of a flat section in the intervals [-ymax, nl1] and [nl2, ymax], and a curved surface
!!!!!!!!!!! in the interval [nl1, nl2]. All measurements are relative to the origin (center of the surface).
!!!!!!!!!!! If nl1 < -ymax and nl2 > ymax, the entire surface in the interval [-ymax, ymax] becomes curved.
!!!!!!!!!!! The surface along the x-direction remains unchanged.
!
            if(jdif.gt.0.and.jdif.le.nl2) pgzi(jdif)=-concav*sqrt(pgyi(jdif)**2+(10.61)**2)+k  ! Function to generate a curved surface in the Y direction: z=f(y)
            if(jdif.ge.nl2) pgzi(jdif)= pgzi(jdif-1)
            if(jdif.lt.0.and.jdif.gt.nl1) pgzi(jdif)=-concav*sqrt(pgyi(jdif)**2+(10.61)**2)+k  ! Function to generate a curved surface in the Y direction: z=f(y)
            if(jdif.le.nl1) pgzi(jdif)= pgzi(jdif+1)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(j.le.0) then
            deltary=(pgyi(jdif)-pgyi(jdif-1))
            deltarz=-(pgzi(jdif)-pgzi(jdif-1))
            end if
            if(j.gt.0) then
            deltary=-(pgyi(jdif)-pgyi(jdif+1))
            deltarz=(pgzi(jdif)-pgzi(jdif+1))
            end if
            r=sqrt(deltary**2+deltarz**2)
            deltas=abs(r-1d0)
            if(deltas.le.0.0001) then
            gy(ip)=jdif
            pgz(ip)=pgzi(jdif)
            pgy(ip)=pgyi(jdif)
            pgx(ip)=i
            vecnormalx(ip)=0d0
            vecnormalz(ip)=deltary
            vecnormaly(ip)=deltarz
            grho1(ip)=asin(vecnormaly(ip))
            go to 222
            end if
            end do
            end if
!
       end if
!
 222   end if
!
       end if
!
!!!!!! Cylinder !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
            if(nano.eq.2.or.nano.eq.3) then
            if(nano.eq.2) then
            rho=(2d0*pi)/(2d0*ymax+dy+1d0)
            raio_circun=(1d0/2d0)/(1d0-(cos(rho/2d0))**2d0)**0.5d0+(k-1)
            if(ip.eq.1) write(*,*)'raio_circum=',raio_circun
            end if
            if(nano.eq.3) then
            rho=2d0*acos((1d0-1d0/(2d0*raio_circun)**2)**0.5d0)
            if(ip.eq.1) write(*,*)rho
            end if
            if(raio_circun.ne.0) then
            rho1=-concav*(j+ymax)*rho
            grho1(ip)=rho1-pi*(1d0+concav)/2d0
            pgz(ip)=raio_circun*cos(rho1)+0*(k-1)
            pgx(ip)=i
            pgy(ip)=raio_circun*sin(rho1)
            r=dsqrt(pgy(ip)**2+pgz(ip)**2)
            vecnormaly(ip)=-concav*pgy(ip)/r
            vecnormalx(ip)=0d0
            vecnormalz(ip)=-concav*pgz(ip)/r
            end if
            end if
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Magnetic moments
!
        r=dsqrt(real((gx(ip)-il))**2d0+real((gy(ip)-l)**2d0))
!
        if(r.ne.0) then
!
        cosphi1=real(gx(ip)-il)/r
        senophi1=real(gy(ip)-l)/r
        cosphi=cosphi1*cos(gama_tip)-m_tip*senophi1*sin(gama_tip)
        senophi=m_tip*senophi1*cos(gama_tip)+cosphi1*sin(gama_tip)
        teta=-((pi)*(r/a)**(2d0/(1d0 - s)))/(1d0 + (r/a)**(2d0/(1d0 - s)))+pi
        senoteta=dsin(teta)
        costeta=dcos(teta)
!
        else
!
        cosphi=0
        senophi=0
        teta=-(pi*(r/a)**(2d0/(1d0 - s)))/(1d0 + (r/a)**(2d0/(1d0 - s)))+pi
        senoteta=dsin(teta)
        costeta=dcos(teta)
!
        end if
!
!!!!! Skyrmion tipo Neel !!!!!!!
!
    if(tipo.eq.1) then
         sinal(ip)=1d0
         ux=cosphi*senoteta
         uy=senophi*senoteta
         uz=dcos(teta)
         if(nano.eq.1) then
         if(direc_tipo.eq.0) then
         uyi(ip)=uy
         uzi(ip)=(uz*cos(grho1(ip))-ux*sin(grho1(ip)))
         uxi(ip)=(uz*sin(grho1(ip))+ux*cos(grho1(ip)))
         end if
         if(direc.eq.1) then
         if(direc_tipo.eq.4.or.direc_tipo.eq.5) then
         uyi(ip)=uy
         uzi(ip)=(uz*cos(grho1(ip))-ux*sin(grho1(ip)))
         uxi(ip)=(uz*sin(grho1(ip))+ux*cos(grho1(ip)))
         end if
         if(direc_tipo.eq.6) then
         uxi(ip)=ux
         uzi(ip)=(uz*cos(grho1(ip))+uy*sin(grho1(ip)))
         uyi(ip)=(uz*sin(grho1(ip))-uy*cos(grho1(ip)))
         end if
         end if
         if(direc.eq.2) then
         if(direc_tipo.eq.4.or.direc_tipo.eq.5) then
         uxi(ip)=ux
         uzi(ip)=(uz*cos(grho1(ip))-uy*sin(grho1(ip)))
         uyi(ip)=(uz*sin(grho1(ip))+uy*cos(grho1(ip)))
         end if
         end if
    end if
!
        if(nano.eq.2.or.nano.eq.3) then
        uxi(ip)=ux
        uzi(ip)=(uz*cos(grho1(ip))-uy*sin(grho1(ip)))
        uyi(ip)=(uz*sin(grho1(ip))+uy*cos(grho1(ip)))
        end if
!         
    end if
!
!
!!!! Skyrmion AFM !!!!!!!!!!!!!!
!
    if(tipo.eq.2) then
        sinal(ip)=(-1)**(abs(i))*(-1)**(abs(j))*(-1)**(k)
        ux=cosphi*senoteta*sinal(ip)
        uy=senophi*senoteta*sinal(ip)
        uz=dcos(teta)*sinal(ip)
        if(nano.eq.1) then
        if(direc_tipo.eq.0) then
         uyi(ip)=uy
         uzi(ip)=(uz*cos(grho1(ip))-ux*sin(grho1(ip)))
         uxi(ip)=(uz*sin(grho1(ip))+ux*cos(grho1(ip)))
         end if
         if(direc.eq.1) then
         if(direc_tipo.eq.4.or.direc_tipo.eq.5) then
         uyi(ip)=uy
         uzi(ip)=(uz*cos(grho1(ip))-ux*sin(grho1(ip)))
         uxi(ip)=(uz*sin(grho1(ip))+ux*cos(grho1(ip)))
         end if
         if(direc_tipo.eq.6) then
         uxi(ip)=ux
         uzi(ip)=(uz*cos(grho1(ip))+uy*sin(grho1(ip)))
         uyi(ip)=(uz*sin(grho1(ip))-uy*cos(grho1(ip)))
         end if
         end if
         if(direc.eq.2) then
         if(direc_tipo.eq.4.or.direc_tipo.eq.5) then
         uxi(ip)=ux
         uzi(ip)=(uz*cos(grho1(ip))-uy*sin(grho1(ip)))
         uyi(ip)=(uz*sin(grho1(ip))+uy*cos(grho1(ip)))
         end if
         end if
         end if
!
        if(nano.eq.2.or.nano.eq.3) then
        uxi(ip)=ux
        uzi(ip)=(uz*cos(grho1(ip))-uy*sin(grho1(ip)))
        uyi(ip)=(uz*sin(grho1(ip))+uy*cos(grho1(ip)))
        end if
    end if
!    
!!!! parede de domínio neell !!!!!!
!
    if(tipo.eq.3) then
    sinal(ip)=1
    teta=-(2d0*atan(exp((float(i)-1.0*(il))/raio)))
	if(direc.eq.1.and.direc_tipo.eq.5) then
  	if(i.gt.0) teta=-(2d0*atan(exp((float(-i)-1d0*(il))/raio)))
  	if(i.le.0) teta=-(2d0*atan(exp((float(i+xmax)-1d0*(il))/raio)))
  	end if
  	if(direc.eq.2.and.direc_tipo.eq.5) then
  	if(j.gt.0) teta=-(2d0*atan(exp((float(-j)-1d0*(l))/raio)))
  	if(j.le.0) teta=-(2d0*atan(exp((float(j+ymax)-1d0*(l))/raio)))
  	end if
    ux=dsin(teta)
	uy=0d0
    uz=+dcos(teta)
    if(nano.eq.1) then
         if(direc_tipo.eq.0) then
         uyi(ip)=uy
         uzi(ip)=(uz*cos(grho1(ip))-ux*sin(grho1(ip)))
         uxi(ip)=(uz*sin(grho1(ip))+ux*cos(grho1(ip)))
         end if
         if(direc.eq.1) then
         if(direc_tipo.eq.4.or.direc_tipo.eq.5) then
         uyi(ip)=uy
         uzi(ip)=(uz*cos(grho1(ip))-ux*sin(grho1(ip)))
         uxi(ip)=(uz*sin(grho1(ip))+ux*cos(grho1(ip)))
         end if
         if(direc_tipo.eq.6) then
         uxi(ip)=ux
         uzi(ip)=(uz*cos(grho1(ip))+uy*sin(grho1(ip)))
         uyi(ip)=(uz*sin(grho1(ip))-uy*cos(grho1(ip)))
         end if
         end if
         if(direc.eq.2) then
         if(direc_tipo.eq.4.or.direc_tipo.eq.5) then
         uxi(ip)=ux
         uzi(ip)=(uz*cos(grho1(ip))-uy*sin(grho1(ip)))
         uyi(ip)=(uz*sin(grho1(ip))+uy*cos(grho1(ip)))
         end if
         end if
         end if
!
        if(nano.eq.2.or.nano.eq.3) then
        uxi(ip)=ux
        uzi(ip)=(uz*cos(grho1(ip))-uy*sin(grho1(ip)))
        uyi(ip)=(uz*sin(grho1(ip))+uy*cos(grho1(ip)))
        end if
    end if
!    
!!!! parede de domínio tipo tail to tail ou head to head !!!!!!
!
    if(tipo.eq.4.or.tipo.eq.-4) then
    sinal(ip)=1
!!!!  trocando o sinal de teta inverte o sentido da magnetização na direção y !!!!!!
    teta=-((tipo/abs(tipo)-1)*pi/2+(tipo/abs(tipo))*(2*atan(exp((float(i)-1.0*(il))/raio))))
    if(direc.eq.1) then
    if(direc_tipo.eq.5) then
  	if(i.gt.0) teta=((tipo/abs(tipo)-1)*pi/2+(tipo/abs(tipo))*(2*atan(exp((float(-i)-1.0*(il))/raio))))
  	if(i.le.0) teta=((tipo/abs(tipo)-1)*pi/2+(tipo/abs(tipo))*(2*atan(exp((float(i+xmax)-1.0*(il))/raio))))
  	end if
  	end if
  	if(direc.eq.2.and.direc_tipo.eq.5) then
  	if(j.gt.0) teta=((tipo/abs(tipo)-1)*pi/2+(tipo/abs(tipo))*(2*atan(exp((float(-j)-1.0*(l))/raio))))
  	if(j.le.0) teta=((tipo/abs(tipo)-1)*pi/2+(tipo/abs(tipo))*(2*atan(exp((float(j+ymax)-1.0*(l))/raio))))
  	end if
	ux=dcos(teta)
	uy=dsin(teta)
	uz=0
    if(nano.eq.1) then
         if(direc_tipo.eq.0) then
         uyi(ip)=uy
         uzi(ip)=(uz*cos(grho1(ip))-ux*sin(grho1(ip)))
         uxi(ip)=(uz*sin(grho1(ip))+ux*cos(grho1(ip)))
         end if
         if(direc.eq.1) then
         if(direc_tipo.eq.4.or.direc_tipo.eq.5) then
         uyi(ip)=uy
         uzi(ip)=(uz*cos(grho1(ip))-ux*sin(grho1(ip)))
         uxi(ip)=(uz*sin(grho1(ip))+ux*cos(grho1(ip)))
         end if
         if(direc_tipo.eq.6) then
         uxi(ip)=ux
         uzi(ip)=(uz*cos(grho1(ip))+uy*sin(grho1(ip)))
         uyi(ip)=(uz*sin(grho1(ip))-uy*cos(grho1(ip)))
         end if
         end if
         if(direc.eq.2) then
         if(direc_tipo.eq.4.or.direc_tipo.eq.5) then
         uxi(ip)=uy
         uzi(ip)=(uz*cos(grho1(ip))-ux*sin(grho1(ip)))
         uyi(ip)=(uz*sin(grho1(ip))+ux*cos(grho1(ip)))
         end if
         end if
         end if
!
        if(nano.eq.2.or.nano.eq.3) then
        uxi(ip)=ux
        uzi(ip)=(uz*cos(grho1(ip))-uy*sin(grho1(ip)))
        uyi(ip)=(uz*sin(grho1(ip))+uy*cos(grho1(ip)))
        end if
    end if
!
    if(abs(tipo).eq.6) then
    sinal(ip)=1
    teta=dsqrt((i-x0)**2+(j-y0)**2+raio**2)
    if(teta.ne.0) then
	uy=(i-x0)/teta
	ux=(tipo*(j-y0)/teta)/abs(tipo)
    uz=raio/teta
    else
    uy=0d0
	ux=0d0
    uz=1d0
    end if
    if(nano.eq.1) then
         if(direc_tipo.eq.0) then
         uyi(ip)=uy
         uzi(ip)=(uz*cos(grho1(ip))-ux*sin(grho1(ip)))
         uxi(ip)=(uz*sin(grho1(ip))+ux*cos(grho1(ip)))
         end if
         if(direc.eq.1) then
         if(direc_tipo.eq.4.or.direc_tipo.eq.5) then
         uyi(ip)=uy
         uzi(ip)=(uz*cos(grho1(ip))-ux*sin(grho1(ip)))
         uxi(ip)=(uz*sin(grho1(ip))+ux*cos(grho1(ip)))
         end if
         end if
         if(direc_tipo.eq.6) then
         uxi(ip)=ux
         uzi(ip)=(uz*cos(grho1(ip))+uy*sin(grho1(ip)))
         uyi(ip)=(uz*sin(grho1(ip))-uy*cos(grho1(ip)))
         end if
         if(direc.eq.2) then
         if(direc_tipo.eq.4.or.direc_tipo.eq.5) then
         uxi(ip)=ux
         uzi(ip)=(uz*cos(grho1(ip))-uy*sin(grho1(ip)))
         uyi(ip)=(uz*sin(grho1(ip))+uy*cos(grho1(ip)))
         end if
         end if
         end if
         end if
!
        if(nano.eq.2.or.nano.eq.3) then
        uxi(ip)=ux
        uzi(ip)=(uz*cos(grho1(ip))-uy*sin(grho1(ip)))
        uyi(ip)=(uz*sin(grho1(ip))+uy*cos(grho1(ip)))
        end if
!
    write(1,100)'He',pgx(ip),pgy(ip),pgz(ip),'atom_vector',vecnormalx(ip),vecnormaly(ip),vecnormalz(ip)
!
    end do
!
       end do
!
        end do
!
    ip=0
	do i=-xmax,xmax+ncx
!
 		do j=-ymax,ymax+ncy+dy
!
            do k=1,zmax
!
            ip=ip+1
!
        if(impu.eq.1) then
        if(forma_imp.eq.0) then
        if(k.eq.1.and.gx(ip).ge.pxi.and.gx(ip).le.pxf.and.gy(ip).ge.pyi.and.gy(ip).le.pyf.and.ty(ip).eq.0) then
        ty(ip)=1
        if(magty(1).eq.0) then
        uxi(ip)=0d0
        uyi(ip)=0d0
        uzi(ip)=0d0
        end if
        end if
        end if
        if(forma_imp.eq.1) then
        rmin=((gx(ip)-pxic)**2+(gy(ip)-pyjc)**2)**0.5
        if(rmin.le.raio_imp.and.rmin.ge.raio_imp1.and.ty(ip).eq.0) then
        ty(ip)=1
        if(magty(1).eq.0) then
        uxi(ip)=0d0
        uyi(ip)=0d0
        uzi(ip)=0d0
        end if
        end if
        end if
        end if
!
    end do
!
       end do
!
        end do
!
 100 format(A2,1X,F18.13,1x,F18.13,1X,F18.13,1X,A11,3(2x,F9.5))
 300 format(A2,1X,I3,1x,I3,1X,I3,1X,A11,3(2x,F9.5))
 200 format(I2,2x,3(2x,F18.13))
    flush(1)
    close(1)
    close(45)
    return
    end 
!
!
!
    subroutine escreve(it1)
    use VARIAVEIS
    integer :: i,j,k,ip,it1,rx0,ry0
    real :: somax,somay,somaz,pgxr,pgyr,pgzr,zero
    real :: rds,uzipeds,uzipesx,uzipesy,rcvex,rcvey,uzipdsx,uzipdsy,rcvdx,rcvdy
    real :: rcvmx,rcvmy
!
    open(11, file='config.xyz')
    if(tipo.eq.2) then
    open(42,file='config_afm_a.xyz')
    open(44,file='config_afm_b.xyz')
    write(42,*) (2*xmax+1+ncx)*(2*ymax+dy+1+ncy)*zmax
	write(42,*)'skyrmion_a'
    write(44,*) (2*xmax+1+ncx)*(2*ymax+dy+1+ncy)*zmax
	write(44,*)'skyrmion_b'
    end if
    open(12,file='energia.dat')
    open(13, file='comp_mag_xyz.dat')
    open(18, file='posicao_energia_dw.dat')
    open(19, file='posicao_energia_plano_dw.dat')
	write(11,*) (2*xmax+1+ncx)*(2*ymax+dy+1+ncy)*zmax+zmax
	write(11,*)'skyrmion'
!
    soma=0d0;somax=0d0;somay=0d0;somaz=0d0
    zero=0d0
    uzipe=0
    xcve=0d0;ycve=0d0;zcve=0d0;zcvd=0d0
    uzipd=0d0
    xcvd=0d0;ycvd=0d0
    xcvei=0d0;ycvej=0d0
    xcvdi=0d0;ycvdj=0d0
    rcvex=0d0;rcvey=0d0
    uzipesx=0d0;uzipesy=0d0
    rcvdx=0d0;rcvdy=0d0
    uzipdsx=0d0;uzipdsy=0d0
    rcvmx=0d0;rcvmy=0d0
!    
    do ip=1,npt
!
        i=gx(ip)
        j=gy(ip)
        k=gz(ip)
        pgxr=pgx(ip)
        pgyr=pgy(ip)
        pgzr=pgz(ip)
        rx0=int(x0)
        ry0=int(y0)
        ivxd1=gvxd1(ip)
        ivxe1=gvxe1(ip)
!
        somax=somax+uxi(ip)
        somay=somay+uyi(ip)
        somaz=somaz+uzi(ip)
        soma=soma+sqrt(uxi(ip)**2+uyi(ip)**2+uzi(ip)**2)
!
    if(sinal(ip).eq.-1.and.ty(ip).eq.0) write(11,100)'He',pgxr,pgyr,pgzr,'atom_vector',&
    uxi(ip),uyi(ip),uzi(ip)
    if(sinal(ip).eq.1.and.ty(ip).eq.0)  write(11,100)'Ni',pgxr,pgyr,pgzr,'atom_vector',&
    uxi(ip),uyi(ip),uzi(ip)
    if(sinal(ip).eq.1.and.ty(ip).eq.1)  write(11,100)'Au',pgxr,pgyr,pgzr,'atom_vector',&
    uxi(ip),uyi(ip),uzi(ip)
    if(sinal(ip).eq.-1.and.ty(ip).eq.1)  write(11,100)'Au',pgxr,pgyr,pgzr,'atom_vector',&
    uxi(ip),uyi(ip),uzi(ip)

!
    if(tipo.eq.1.or.tipo.eq.6) then
!
!     if(j.ge.-ymax+5.and.j.le.ymax+ncy-5) then
!     if(i.ge.-xmax+5.and.i.le.xmax+ncx-5) then
!
    rds=sqrt((float(i)-x0)**2+(float(j)-y0)**2)
    if(rds.le.lim_raio.and.rds.ne.0) then
!
     uziped=exp(-20d0*(1d0+(uxi(ip)*vecnormalx(ip)+uyi(ip)*vecnormaly(ip)+uzi(ip)*vecnormalz(ip))))
     if(tipo.eq.6) uziped=exp(-10d0*(1d0-(uxi(ip)*vecnormalx(ip)+uyi(ip)*vecnormaly(ip)+uzi(ip)*vecnormalz(ip))))
     uzipe=uzipe+uziped
     xcve=xcve+pgxr*uziped
     ycve=ycve+pgyr*uziped
     zcve=zcve+pgzr*uziped
!
    xcvei=xcvei+i*uziped
    ycvej=ycvej+j*uziped
!
    end if

!     end if
!     end if
    end if
!
    if(tipo.eq.3) then
!
    rds=sqrt((float(i)-x0)**2)
    if(rds.le.lim_raio.and.j.eq.int(ry0)) then
     uziped=exp(-20d0*(1d0-abs(uxi(ip)*tangxxd(ip)+uyi(ip)*tangyxd(ip)+uzi(ip)*tangzxd(ip))))
     uzipe=uzipe+uziped
     xcve=xcve+pgxr*uziped
     ycve=ycve+pgyr*uziped
     zcve=zcve+pgzr*uziped
!
     xcvei=xcvei+i*uziped
     ycvej=ycvej+j*uziped
!
    end if
    end if
!
    if(abs(tipo).eq.4) then
!
    rds=sqrt((float(i)-x0)**2)
    if(rds.le.lim_raio.and.j.eq.int(ry0)) then
     uziped=exp(-20d0*(1d0-abs(uxi(ip)*tangxyd(ip)+uyi(ip)*tangyyd(ip)+uzi(ip)*tangzyd(ip))))
     uzipe=uzipe+uziped
     xcve=xcve+pgxr*uziped
     ycve=ycve+pgyr*uziped
     zcve=zcve+pgzr*uziped
!
     xcvei=xcvei+i*uziped
     ycvej=ycvej+j*uziped
!
    end if
    end if
!
    if(tipo.eq.2) then
!    
    if(sinal(ip).eq.-1.and.ty(ip).eq.0) write(42,100)'He',pgxr,pgyr,pgzr,'atom_vector',&
    uxi(ip),uyi(ip),uzi(ip)
    if(sinal(ip).eq.-1.and.ty(ip).eq.1)  write(42,100)'Cu',pgxr,pgyr,pgzr,'atom_vector',&
    uxi(ip),uyi(ip),uzi(ip)
    if(sinal(ip).eq.1.and.ty(ip).eq.0)  write(42,100)'Ni',pgxr,pgyr,pgzr,'atom_vector',&
    zero,zero,zero
    if(sinal(ip).eq.1.and.ty(ip).eq.1)  write(42,100)'Cu',pgxr,pgyr,pgzr,'atom_vector',&
    zero,zero,zero
!
    if(sinal(ip).eq.-1) then
!     if(j.ge.-ymax+5.and.j.le.ymax+ncy-5) then
!     if(i.ge.-xmax+5.and.i.le.xmax+ncx-5) then
!
    rds=sqrt((i-x0)**2+(j-y0)**2)
    if(rds.le.lim_raio) then
     uziped=exp(20d0*(1d0+(uxi(ip)*vecnormalx(ip)+uyi(ip)*vecnormaly(ip)+uzi(ip)*vecnormalz(ip))))
    uzipe=uzipe+uziped
    xcve=xcve+pgxr*uziped
    ycve=ycve+pgyr*uziped
    zcve=zcve+pgzr*uziped
!
    xcvei=xcvei+i*uziped
    ycvej=ycvej+j*uziped
!
    end if
!    
!     end if
!     end if
    end if
!
    if(sinal(ip).eq.1) then
!     if(j.ge.-ymax+5.and.j.le.ymax+ncy-5) then
!     if(i.ge.-xmax+5.and.i.le.xmax+ncx-5) then
!
    rds=sqrt((i-x0)**2+(j-y0)**2)
    if(rds.le.lim_raio) then
    uziped=exp(20d0*(1d0-(uxi(ip)*vecnormalx(ip)+uyi(ip)*vecnormaly(ip)+uzi(ip)*vecnormalz(ip))))
    uzipd=uzipd+uziped
    xcvd=xcvd+pgxr*uziped
    ycvd=ycvd+pgyr*uziped
    zcvd=zcvd+pgzr*uziped
!
    xcvdi=xcvdi+i*uziped
    ycvdj=ycvdj+j*uziped
!
    end if
!
!     end if
!     end if
    end if
!
    if(sinal(ip).eq.-1.and.ty(ip).eq.0) write(44,100)'He',pgxr,pgyr,pgzr,'atom_vector',&
    zero,zero,zero
    if(sinal(ip).eq.-1.and.ty(ip).eq.1)  write(44,100)'Cu',pgxr,pgyr,pgzr,'atom_vector',&
    zero,zero,zero
    if(sinal(ip).eq.1.and.ty(ip).eq.0)  write(44,100)'Ni',pgxr,pgyr,pgzr,'atom_vector',&
    uxi(ip),uyi(ip),uzi(ip)
    if(sinal(ip).eq.1.and.ty(ip).eq.1)  write(44,100)'Cu',pgxr,pgyr,pgzr,'atom_vector',&
    uxi(ip),uyi(ip),uzi(ip)
!
    end if
!    
    end do
!
        if(uzipe.ne.0) then
        if(tipo.eq.1.or.tipo.eq.6.or.tipo.eq.2.or.abs(tipo).eq.4.or.tipo.eq.3) then
        xcve=xcve/uzipe
        ycve=ycve/uzipe
        zcve=zcve/uzipe
!        
        xcvei=xcvei/uzipe
        ycvej=ycvej/uzipe
!
        end if
        end if
        if(uzipe.ne.0) then
        if(tipo.eq.2) then
        xcvd=xcvd/uzipd
        ycvd=ycvd/uzipd
        zcvd=zcvd/uzipe
!
        xcvdi=xcvdi/uzipd
        ycvdj=ycvdj/uzipd
!
        end if
        end if
!
        if(tipo.eq.1.or.tipo.eq.6.or.abs(tipo).eq.4.or.tipo.eq.3) then
        xcvm=xcve
        ycvm=ycve
        zcvm=zcve
!
        xcvmi=xcvei
        ycvmj=ycvej
!       
        end if
!        
        if(tipo.eq.2) then
        xcvm=0.5d0*(xcve+xcvd)
        ycvm=0.5d0*(ycve+ycvd)
        zcvm=0.5d0*(zcve+zcvd)
!
        xcvmi=0.5d0*(xcvei+xcvdi)
        ycvmj=0.5d0*(ycvej+ycvdj)
!
        end if
! 
        x0=xcvmi
        y0=ycvmj
        rx0=int(x0)
        ry0=int(y0)
!
    do ip=1,npt
!
        i=gx(ip)
        j=gy(ip)
        k=gz(ip)
        pgxr=pgx(ip)
        pgyr=pgy(ip)
        pgzr=pgz(ip)
!    
    if(tipo.eq.1.or.tipo.eq.6) then
!     if(j.ge.-ymax+5.and.j.le.ymax+ncy-5) then
!     if(i.ge.-xmax+5.and.i.le.xmax+ncx-5) then
!
    rds=sqrt(float(i-rx0)**2+float(j-ry0)**2)

    if(rds.le.lim_raio1.and.rds.ne.0) then
       uzipeds=exp(-20*abs(uxi(ip)*vecnormalx(ip)+uyi(ip)*vecnormaly(ip)+uzi(ip)*vecnormalz(ip)))
       uzipesx=uzipesx+uzipeds
       rcvex=rcvex+rds*uzipeds
!
       uzipeds=exp(-20*abs(uxi(ip)*vecnormalx(ip)+uyi(ip)*vecnormaly(ip)+uzi(ip)*vecnormalz(ip)))
       uzipesy=uzipesy+uzipeds
       rcvey=rcvey+rds*uzipeds
!
        end if

!        end if
!        end if
!    
    end if
!
    if(tipo.eq.3) then
!
    rds=sqrt((float(i)-rx0)**2)
    if(rds.le.lim_raio1) then
!
      if(j.eq.0) then
      uzipeds=exp(-20*abs(0.648-abs(uxi(ip)*tangxxd(ip)+uyi(ip)*tangyxd(ip)+uzi(ip)*tangzxd(ip))))
      uzipesx=uzipesx+uzipeds
      rcvex=rcvex+rds*uzipeds
      end if
!
    end if
    end if
!
    if(abs(tipo).eq.4) then
!
    rds=sqrt((float(i)-rx0)**2)
    if(rds.le.lim_raio1) then
!
      if(j.eq.0) then
      uzipeds=exp(-20*abs(0.648-abs(uxi(ip)*tangxyd(ip)+uyi(ip)*tangyyd(ip)+uzi(ip)*tangzyd(ip))))
      uzipesx=uzipesx+uzipeds
      rcvex=rcvex+rds*uzipeds
      end if
!
    end if
    end if
!
    if(tipo.eq.2) then
!
    if(sinal(ip).eq.-1) then
!     if(j.ge.-ymax+5.and.j.le.ymax+ncy-5) then
!     if(i.ge.-xmax+5.and.i.le.xmax+ncx-5) then
!
    rds=sqrt(float(i-rx0)**2+float(j-ry0)**2)
    if(rds.le.lim_raio1) then
!
      uzipeds=exp(-20*abs(uxi(ip)*vecnormalx(ip)+uyi(ip)*vecnormaly(ip)+uzi(ip)*vecnormalz(ip)))
      uzipesx=uzipesx+uzipeds
      rcvex=rcvex+rds*uzipeds
!
      uzipeds=exp(-20*abs(uxi(ip)*vecnormalx(ip)+uyi(ip)*vecnormaly(ip)+uzi(ip)*vecnormalz(ip)))
      uzipesy=uzipesy+uzipeds
      rcvey=rcvey+rds*uzipeds
!    
    end if
!    
!     end if
!     end if
!
    end if
!
    if(sinal(ip).eq.1) then
!     if(j.ge.-ymax+5.and.j.le.ymax+ncy-5) then
!     if(i.ge.-xmax+5.and.i.le.xmax+ncx-5) then
!
    rds=sqrt(float(i-rx0)**2+float(j-ry0)**2)
    if(rds.le.lim_raio1) then
!
      uzipeds=exp(-20*abs(uxi(ip)*vecnormalx(ip)+uyi(ip)*vecnormaly(ip)+uzi(ip)*vecnormalz(ip)))
      uzipdsx=uzipdsx+uzipeds
      rcvdx=rcvdx+rds*uzipeds
!
      uzipeds=exp(-20*abs(uxi(ip)*vecnormalx(ip)+uyi(ip)*vecnormaly(ip)+uzi(ip)*vecnormalz(ip)))
      uzipdsy=uzipdsy+uzipeds
      rcvdy=rcvdy+rds*uzipeds
!    
    end if
!
!     end if
!     end if
!
    end if
!
    end if
!
    if(i.eq.rx0.and.j.eq.ry0) write(11,100)'Cu',xcvm,ycvm,zcvm,'atom_vector',0.,0.,0.
!    
    end do
!
        if(tipo.eq.1.or.tipo.eq.6.or.tipo.eq.2) then
!
        rcvex=rcvex/uzipesx
        rcvey=rcvey/uzipesy
!
        end if
!
        if(abs(tipo).eq.4.or.tipo.eq.3) then
!
        rcvex=rcvex/uzipesx
!
        end if
!        
        if(tipo.eq.2) then
!
        rcvdx=rcvdx/uzipdsx
        rcvdy=rcvdy/uzipdsy
!
        end if
!
        if(tipo.eq.1.or.tipo.eq.6.or.abs(tipo).eq.4.or.tipo.eq.3) then
!
        rcvmx=rcvex
        rcvmy=rcvey
!       
        end if
!        
        if(tipo.eq.2) then
!
        rcvmx=0.5*(rcvex+rcvdx)
        rcvmy=0.5*(rcvey+rcvdy)
!
        end if
!
        ! energiatot ⇒ Total magnetic energy
        !
        ! energdip ⇒ Dipolar energy
        !
        ! dzym ⇒ Dzyaloshinskii–Moriya interaction (DMI) energy
        !
        ! troca ⇒ Heisenberg exchange energy
        !
        ! anisotro ⇒ Magnetocrystalline anisotropy energy
        !
        ! xcvm ⇒ X-coordinate of the skyrmion (or domain wall) center
        !
        ! ycvm ⇒ Y-coordinate of the skyrmion (or domain wall) center
        !
        ! xcvmi ⇒ ξ-coordinate of the skyrmion (or domain wall) center along the length of the nanowire
        !
        ! ycvmj ⇒ δ-coordinate of the skyrmion (or domain wall) center perpendicular to the ξ-axis
        !
        ! rcvmx ⇒ Skyrmion radius along the x-direction
        !
        ! rcvmy ⇒ Skyrmion radius along the y-direction
!
        ! soma/npt – normalized total magnetization
        !
        ! somax/npt, somay/npt, somaz/npt – normalized total magnetization components along x, y, and z, respectively
!
        write(12,200)energiatot,energdip,dzym,troca,anisotro
        write(*,200)energiatot,energdip,dzym,troca,anisotro
        write(*,*)xcvmi,ycvmj,rcvmx,rcvmy
        write(13,800)soma/npt,somax/npt,somay/npt,somaz/npt
        write(18,900)energiatot,energdip,dzym,troca,anisotro,xcvm,ycvm,rcvmx,rcvmy
        write(19,900)energiatot,energdip,dzym,troca,anisotro,xcvmi,ycvmj,rcvmx,rcvmy
        flush(11)
        flush(12)
        flush(13)
        flush(18)
        flush(19)
!
 100 format(A2,1X,F9.4,1x,F9.4,1X,F9.4,1X,A11,3(2x,F9.4))
 200 format(5(2x,F19.12))
 300 format(I2,2x,3(2x,F18.13))
 400 format(3(2x,F9.5))
 600 format(6(2x,F9.5))
 700 format(I4,2x,I4,2x,I2,2x,4(2x,F11.5))
 800 format(4(1x,F15.10))
 900 format(5(1x,F19.12),4(1x,F19.12))
!
    return
    end 
!
!
!    
    subroutine tabela
    use VARIAVEIS
    integer :: ip,ip1,i,j,k,iv,jv,kv,ivp,ivvp,nvp,rvijx,rvijy,rvijz
    real(kind=8) :: rdij,ipg,jpg,kpg,vipg,vjpg,vkpg,rvijxx,rvijyy,rvijzz,modtg,corte
    real(kind=8) :: dmxxd,dmyxd,dmzxd,dmxxe,dmyxe,dmzxe,dmxyd,dmyyd,dmzyd,dmxye,dmyye,dmzye
    real(kind=8) :: dmxzd,dmyzd,dmzzd,dmxze,dmyze,dmzze,ap2
    real(kind=8) :: ax,ay,az,bx,by,bz,cx,cy,cz

!
    open(69, file='tangente.xyz')
    write(69,*)npt*6
    write(69,*)'vetor tangente'
!
    gvxd1=0;gvxe1=0;gvyd1=0;gvye1=0;gvzd1=0;gvze1=0
    usa=0
    nvtp=0    
    ivp=0
    volcel=1d0
    dmcelxd=1d0;dmcelxe=1d0;dmcelyd=1d0;dmcelye=1d0;dmcelzd=1d0;dmcelze=1d0
!
    if(flag_cut.eq.0) corte=100000
    if(flag_cut.eq.1) corte=cut
!    
	do ip=1,npt
!
    usa(ip)=1
    i=gx(ip)
    j=gy(ip)
    k=gz(ip)
    ipg=pgx(ip)
    jpg=pgy(ip)
    kpg=pgz(ip)
!
    nvp=0
    ivvp=0
!
    do ip1=1,npt
!
        iv=gx(ip1)
        jv=gy(ip1)
        kv=gz(ip1)
!        
        ivvp=ivvp+1
        ivvp=ip1
		rvijx=i-iv
		rvijy=j-jv
		rvijz=k-kv
!	
        if(usa(ivvp).eq.0) then
!
!!!!!!!! procurando os primeiros vizinhos mais próximos !!!!!!!!
! 
        if(rvijx.eq.-1.and.rvijy.eq.0.and.rvijz.eq.0) then
        gvxd1(ip)=ivvp
        vectangxxd(ip)=pgx(ivvp)-pgx(ip)
        vectangyxd(ip)=pgy(ivvp)-pgy(ip)
        vectangzxd(ip)=pgz(ivvp)-pgz(ip)
        modtg=(vectangxxd(ip)**2+vectangyxd(ip)**2+vectangzxd(ip)**2)**0.5
        modtgxd(ip)=modtg
        vectangxxd(ip)=vectangxxd(ip)/modtg
        vectangyxd(ip)=vectangyxd(ip)/modtg
        vectangzxd(ip)=vectangzxd(ip)/modtg
        gvxe1(ivvp)=ip
        modtgxe(ivvp)=modtg
        vectangxxe(ivvp)=-vectangxxd(ip)
        vectangyxe(ivvp)=-vectangyxd(ip)
        vectangzxe(ivvp)=-vectangzxd(ip)
        end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(rvijx.eq.1.and.rvijy.eq.0.and.rvijz.eq.0)  then
        gvxe1(ip)=ivvp
        vectangxxe(ip)=pgx(ivvp)-pgx(ip)
        vectangyxe(ip)=pgy(ivvp)-pgy(ip)
        vectangzxe(ip)=pgz(ivvp)-pgz(ip)
        modtg=(vectangxxe(ip)**2+vectangyxe(ip)**2+vectangzxe(ip)**2)**0.5
        modtgxe(ip)=modtg
        vectangxxe(ip)=vectangxxe(ip)/modtg
        vectangyxe(ip)=vectangyxe(ip)/modtg
        vectangzxe(ip)=vectangzxe(ip)/modtg
        gvxd1(ivvp)=ip
        modtgxd(ivvp)=modtg
        vectangxxd(ivvp)=-vectangxxe(ip)
        vectangyxd(ivvp)=-vectangyxe(ip)
        vectangzxd(ivvp)=-vectangzxe(ip)
        end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(rvijy.eq.-1.and.rvijx.eq.0.and.rvijz.eq.0) then 
        gvyd1(ip)=ivvp
        vectangxyd(ip)=pgx(ivvp)-pgx(ip)
        vectangyyd(ip)=pgy(ivvp)-pgy(ip)
        vectangzyd(ip)=pgz(ivvp)-pgz(ip)
        modtg=(vectangxyd(ip)**2+vectangyyd(ip)**2+vectangzyd(ip)**2)**0.5
        modtgyd(ip)=modtg
        vectangxyd(ip)=vectangxyd(ip)/modtg
        vectangyyd(ip)=vectangyyd(ip)/modtg
        vectangzyd(ip)=vectangzyd(ip)/modtg
        gvye1(ivvp)=ip
        modtgye(ivvp)=modtg
        vectangxye(ivvp)=-vectangxyd(ip)
        vectangyye(ivvp)=-vectangyyd(ip)
        vectangzye(ivvp)=-vectangzyd(ip)
        end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(rvijy.eq.1.and.rvijx.eq.0.and.rvijz.eq.0)  then
        gvye1(ip)=ivvp
        vectangxye(ip)=pgx(ivvp)-pgx(ip)
        vectangyye(ip)=pgy(ivvp)-pgy(ip)
        vectangzye(ip)=pgz(ivvp)-pgz(ip)
        modtg=(vectangxye(ip)**2+vectangyye(ip)**2+vectangzye(ip)**2)**0.5
        modtgye(ip)=modtg
        vectangxye(ip)=vectangxye(ip)/modtg
        vectangyye(ip)=vectangyye(ip)/modtg
        vectangzye(ip)=vectangzye(ip)/modtg
        gvyd1(ivvp)=ip
        modtgyd(ivvp)=modtg
        vectangxyd(ivvp)=-vectangxye(ip)
        vectangyyd(ivvp)=-vectangyye(ip)
        vectangzyd(ivvp)=-vectangzye(ip)
        end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(nano.eq.2) then
        if(rvijy.eq.2*ymax+dy.and.rvijx.eq.0.and.rvijz.eq.0) then 
        gvyd1(ip)=ivvp
        vectangxyd(ip)=pgx(ivvp)-pgx(ip)
        vectangyyd(ip)=pgy(ivvp)-pgy(ip)
        vectangzyd(ip)=pgz(ivvp)-pgz(ip)
        modtg=(vectangxyd(ip)**2+vectangyyd(ip)**2+vectangzyd(ip)**2)**0.5
        modtgyd(ip)=modtg
        vectangxyd(ip)=vectangxyd(ip)/modtg
        vectangyyd(ip)=vectangyyd(ip)/modtg
        vectangzyd(ip)=vectangzyd(ip)/modtg
        gvye1(ivvp)=ip
        modtgye(ivvp)=modtg
        vectangxye(ivvp)=-vectangxyd(ip)
        vectangyye(ivvp)=-vectangyyd(ip)
        vectangzye(ivvp)=-vectangzyd(ip)
        end if
        if(rvijy.eq.-2*ymax-dy.and.rvijx.eq.0.and.rvijz.eq.0)  then
        gvye1(ip)=ivvp
        vectangxye(ip)=pgx(ivvp)-pgx(ip)
        vectangyye(ip)=pgy(ivvp)-pgy(ip)
        vectangzye(ip)=pgz(ivvp)-pgz(ip)
        modtg=(vectangxye(ip)**2+vectangyye(ip)**2+vectangzye(ip)**2)**0.5
        modtgye(ip)=modtg
        vectangxye(ip)=vectangxye(ip)/modtg
        vectangyye(ip)=vectangyye(ip)/modtg
        vectangzye(ip)=vectangzye(ip)/modtg
        gvyd1(ivvp)=ip
        modtgyd(ivvp)=modtg
        vectangxyd(ivvp)=-vectangxye(ip)
        vectangyyd(ivvp)=-vectangyye(ip)
        vectangzyd(ivvp)=-vectangzye(ip)
        end if
        end if
        
        if(rvijz.eq.-1.and.rvijy.eq.0.and.rvijx.eq.0) then 
        gvzd1(ip)=ivvp
        vectangxzd(ip)=pgx(ivvp)-pgx(ip)
        vectangyzd(ip)=pgy(ivvp)-pgy(ip)
        vectangzzd(ip)=pgz(ivvp)-pgz(ip)
        modtg=(vectangxzd(ip)**2+vectangyzd(ip)**2+vectangzzd(ip)**2)**0.5
        modtgzd(ip)=modtg
        vectangxzd(ip)=vectangxzd(ip)/modtg
        vectangyzd(ip)=vectangyzd(ip)/modtg
        vectangzzd(ip)=vectangzzd(ip)/modtg
        gvze1(ivvp)=ip
        modtgze(ivvp)=modtg
        vectangxze(ivvp)=-vectangxzd(ip)
        vectangyze(ivvp)=-vectangyzd(ip)
        vectangzze(ivvp)=-vectangzzd(ip)
        end if
        if(rvijz.eq.1.and.rvijy.eq.0.and.rvijx.eq.0)  then
        gvze1(ip)=ivvp
        vectangxze(ip)=pgx(ivvp)-pgx(ip)
        vectangyze(ip)=pgy(ivvp)-pgy(ip)
        vectangzze(ip)=pgz(ivvp)-pgz(ip)
        modtg=(vectangxze(ip)**2+vectangyze(ip)**2+vectangzze(ip)**2)**0.5
        modtgze(ip)=modtg
        vectangxze(ip)=vectangxze(ip)/modtg
        vectangyze(ip)=vectangyze(ip)/modtg
        vectangzze(ip)=vectangzze(ip)/modtg
        gvzd1(ivvp)=ip
        modtgzd(ivvp)=modtg
        vectangxzd(ivvp)=-vectangxze(ip)
        vectangyzd(ivvp)=-vectangyze(ip)
        vectangzzd(ivvp)=-vectangzze(ip)
        end if
!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
        vipg=pgx(ip1)
        vjpg=pgy(ip1)
        vkpg=pgz(ip1)
		rvijxx=ipg-vipg
		rvijyy=jpg-vjpg
		rvijzz=kpg-vkpg
		rdij=sqrt(real(rvijxx)**2+real(rvijyy)**2+real(rvijzz)**2)
!		
		if(rdij.ne.0.and.rdij.le.corte) then
!
        ivp=ivp+1
        nvp=nvp+1
        gvv(ivp)=ivvp
!        
        rvij3(ivp)=1d0/rdij**3
!        
        uvrijx(ivp)=real(rvijxx)/rdij
		uvrijy(ivp)=real(rvijyy)/rdij
		uvrijz(ivp)=real(rvijzz)/rdij
!
!
        end if
!
        end if
!
	end do
!
    nvtp(ip)=nvtp(ip-1)+nvp
!
    end do
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!   The program uses a simplified version of the Micromagnetic Finite Element Method (FEM),
!!!!!!   where the volume elements are defined as parallelepipeds whose volume can vary from
!!!!!!   point to point, that is, they can have different edge lengths at each position.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    do ip=1,npt
!
    if(modtgxd(ip).eq.0) modtgxd(ip)=1d0
    if(modtgxe(ip).eq.0) modtgxe(ip)=1d0
    if(modtgyd(ip).eq.0) modtgyd(ip)=1d0
    if(modtgye(ip).eq.0) modtgye(ip)=1d0
    if(modtgzd(ip).eq.0) modtgzd(ip)=1d0
    if(modtgze(ip).eq.0) modtgze(ip)=1d0
    if(vectangxxd(ip).ne.0) then
    ax=(modtgxd(ip))*vectangxxd(ip)
    else
    ax=-(modtgxd(ip))*vectangxxe(ip)
    end if
    if(vectangyxd(ip).ne.0) then
    ay=(modtgxd(ip))*vectangyxd(ip)
    else
    ay=-(modtgxd(ip))*vectangyxe(ip)
    end if
    if(vectangzxd(ip).ne.0) then
    az=(modtgxd(ip))*vectangzxd(ip)
    else
    az=-(modtgxd(ip))*vectangzxe(ip)
    end if
    if(vectangxyd(ip).ne.0) then
    bx=(modtgyd(ip))*vectangxyd(ip)
    else
    bx=-(modtgyd(ip))*vectangxye(ip)
    end if
    if(vectangyyd(ip).ne.0) then
    by=(modtgyd(ip))*vectangyyd(ip)
    else
    by=-(modtgye(ip))*vectangyye(ip)
    end if
    if(vectangzyd(ip).ne.0) then
    bz=(modtgyd(ip))*vectangzyd(ip)
    else
    bz=-(modtgyd(ip))*vectangzye(ip)
    end if
    if(vectangxzd(ip).ne.0) then
    cx=(modtgzd(ip))*vectangxzd(ip)
    else
    cx=-(modtgzd(ip))*vecnormalx(ip)
    end if
    if(vectangyzd(ip).ne.0) then
    cy=(modtgzd(ip))*vectangyzd(ip)
    else
    cy=-(modtgzd(ip))*vecnormaly(ip)
    end if
    if(vectangzzd(ip).ne.0) then
    cz=(modtgzd(ip))*vectangzzd(ip)
    else
    cz=-(modtgzd(ip))*vecnormalz(ip)
    end if

    volcel(ip)=abs(ax*(by*cz-bz*cy)-ay*(bx*cz-bz*cx)+az*(bx*cy-by*cx))
    dist_x(ip)=1d0/(modtgxd(ip)+modtgxe(ip))
    magsavol(ip)=mags*ap**3*volcel(ip)
    dipvol(ip)=dip*volcel(ip)
    anisovol(ip)=ap**3*volcel(ip)/magsavol(ip)
!
    end do
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    do ip=1,npt
!
    write(69,100)'He',pgx(ip),pgy(ip),pgz(ip),'atom_vector',vectangxxd(ip),vectangyxd(ip),vectangzxd(ip)
    write(69,100)'He',pgx(ip),pgy(ip),pgz(ip),'atom_vector',vectangxxe(ip),vectangyxe(ip),vectangzxe(ip)
    write(69,100)'He',pgx(ip),pgy(ip),pgz(ip),'atom_vector',vectangxyd(ip),vectangyyd(ip),vectangzyd(ip)
    write(69,100)'He',pgx(ip),pgy(ip),pgz(ip),'atom_vector',vectangxye(ip),vectangyye(ip),vectangzye(ip)
    write(69,100)'He',pgx(ip),pgy(ip),pgz(ip),'atom_vector',vectangxzd(ip),vectangyzd(ip),vectangzzd(ip)
    write(69,100)'He',pgx(ip),pgy(ip),pgz(ip),'atom_vector',vectangxze(ip),vectangyze(ip),vectangzze(ip)
!
    tangxxd(ip)=vectangxxd(ip)
    tangyxd(ip)=vectangyxd(ip)
    tangzxd(ip)=vectangzxd(ip)
    tangxxe(ip)=vectangxxe(ip)
    tangyxe(ip)=vectangyxe(ip)
    tangzxe(ip)=vectangzxe(ip)
    tangxyd(ip)=vectangxyd(ip)
    tangyyd(ip)=vectangyyd(ip)
    tangzyd(ip)=vectangzyd(ip)
    tangxye(ip)=vectangxye(ip)
    tangyye(ip)=vectangyye(ip)
    tangzye(ip)=vectangzye(ip)
    tangxzd(ip)=vectangxzd(ip)
    tangyzd(ip)=vectangyzd(ip)
    tangzzd(ip)=vectangzzd(ip)
    tangxze(ip)=vectangxze(ip)
    tangyze(ip)=vectangyze(ip)
    tangzze(ip)=vectangzze(ip)
!
    dmxxd=vectangxxd(ip)
    dmyxd=vectangyxd(ip)
    dmzxd=vectangzxd(ip)
    dmxxe=vectangxxe(ip)
    dmyxe=vectangyxe(ip)
    dmzxe=vectangzxe(ip)
    dmxyd=vectangxyd(ip)
    dmyyd=vectangyyd(ip)
    dmzyd=vectangzyd(ip)
    dmxye=vectangxye(ip)
    dmyye=vectangyye(ip)
    dmzye=vectangzye(ip)

    dmxzd=vectangxzd(ip)
    dmyzd=vectangyzd(ip)
    dmzzd=vectangzzd(ip)
    dmxze=vectangxze(ip)
    dmyze=vectangyze(ip)
    dmzze=vectangzze(ip)
!
    ivxd1=gvxd1(ip)
    ivxe1=gvxe1(ip)
    ivyd1=gvyd1(ip)
    ivye1=gvye1(ip)
    ivzd1=gvzd1(ip)
    ivze1=gvze1(ip)
!
    ity=ty(ip)
    isvxd1=ity+ty(ivxd1)
    isvxe1=ity+ty(ivxe1)
    isvyd1=ity+ty(ivyd1)
    isvye1=ity+ty(ivye1)
    isvzd1=ity+ty(ivzd1)
    isvze1=ity+ty(ivze1)
!
    jotacelxd(ip)=ap*jotaty(isvxd1)*0.5d0*(volcel(ip)+volcel(ivxd1))/modtgxd(ip)**2
    jotacelxe(ip)=ap*jotaty(isvxe1)*0.5d0*(volcel(ip)+volcel(ivxe1))/modtgxe(ip)**2
    jotacelyd(ip)=ap*jotaty(isvyd1)*0.5d0*(volcel(ip)+volcel(ivyd1))/modtgyd(ip)**2
    jotacelye(ip)=ap*jotaty(isvye1)*0.5d0*(volcel(ip)+volcel(ivye1))/modtgye(ip)**2
    jotacelzd(ip)=ap*jotaty(isvzd1)*0.5d0*(volcel(ip)+volcel(ivzd1))/modtgzd(ip)**2
    jotacelze(ip)=ap*jotaty(isvze1)*0.5d0*(volcel(ip)+volcel(ivze1))/modtgze(ip)**2
!
    ap2=ap**2
    dmcelxd(ip)=ap2*dmty(isvxd1)*0.5d0*(volcel(ip)+volcel(ivxd1))/modtgxd(ip)
    dmcelxe(ip)=ap2*dmty(isvxe1)*0.5d0*(volcel(ip)+volcel(ivxe1))/modtgxe(ip)
    dmcelyd(ip)=ap2*dmty(isvyd1)*0.5d0*(volcel(ip)+volcel(ivyd1))/modtgyd(ip)
    dmcelye(ip)=ap2*dmty(isvye1)*0.5d0*(volcel(ip)+volcel(ivye1))/modtgye(ip)
    dmcelzd(ip)=ap2*dmty(isvzd1)*0.5d0*(volcel(ip)+volcel(ivzd1))/modtgzd(ip)
    dmcelze(ip)=ap2*dmty(isvze1)*0.5d0*(volcel(ip)+volcel(ivze1))/modtgze(ip)
!
    anisovol(ip)=anisovol(ip)*anisoty(ity)
!    
    vectangxxd(ip)=dmyxd*(vecnormalz(ip)+vecnormalz(ivxd1))/2d0-dmzxd*(vecnormaly(ip)+vecnormaly(ivxd1))/2d0
    vectangyxd(ip)=dmzxd*(vecnormalx(ip)+vecnormalx(ivxd1))/2d0-dmxxd*(vecnormalz(ip)+vecnormalz(ivxd1))/2d0
    vectangzxd(ip)=dmxxd*(vecnormaly(ip)+vecnormaly(ivxd1))/2d0-dmyxd*(vecnormalx(ip)+vecnormalx(ivxd1))/2d0
    
    vectangxxe(ip)=dmyxe*(vecnormalz(ip)+vecnormalz(ivxe1))/2d0-dmzxe*(vecnormaly(ip)+vecnormaly(ivxe1))/2d0
    vectangyxe(ip)=dmzxe*(vecnormalx(ip)+vecnormalx(ivxe1))/2d0-dmxxe*(vecnormalz(ip)+vecnormalz(ivxe1))/2d0
    vectangzxe(ip)=dmxxe*(vecnormaly(ip)+vecnormaly(ivxe1))/2d0-dmyxe*(vecnormalx(ip)+vecnormalx(ivxe1))/2d0
    
    vectangxyd(ip)=dmyyd*(vecnormalz(ip)+vecnormalz(ivyd1))/2d0-dmzyd*(vecnormaly(ip)+vecnormaly(ivyd1))/2d0
    vectangyyd(ip)=dmzyd*(vecnormalx(ip)+vecnormalx(ivyd1))/2d0-dmxyd*(vecnormalz(ip)+vecnormalz(ivyd1))/2d0
    vectangzyd(ip)=dmxyd*(vecnormaly(ip)+vecnormaly(ivyd1))/2d0-dmyyd*(vecnormalx(ip)+vecnormalx(ivyd1))/2d0
    
    vectangxye(ip)=dmyye*(vecnormalz(ip)+vecnormalz(ivye1))/2d0-dmzye*(vecnormaly(ip)+vecnormaly(ivye1))/2d0
    vectangyye(ip)=dmzye*(vecnormalx(ip)+vecnormalx(ivye1))/2d0-dmxye*(vecnormalz(ip)+vecnormalz(ivye1))/2d0
    vectangzye(ip)=dmxye*(vecnormaly(ip)+vecnormaly(ivye1))/2d0-dmyye*(vecnormalx(ip)+vecnormalx(ivye1))/2d0
    
    vectangxzd(ip)=dmyzd*(vecnormalz(ip)+vecnormalz(ivzd1))/2d0-dmzzd*(vecnormaly(ip)+vecnormaly(ivzd1))/2d0
    vectangyzd(ip)=dmzzd*(vecnormalx(ip)+vecnormalx(ivzd1))/2d0-dmxzz*(vecnormalz(ip)+vecnormalz(ivzd1))/2d0
    vectangzzd(ip)=dmxzd*(vecnormaly(ip)+vecnormaly(ivzd1))/2d0-dmyzd*(vecnormalx(ip)+vecnormalx(ivzd1))/2d0
    
    vectangxze(ip)=dmyze*(vecnormalz(ip)+vecnormalz(ivze1))/2d0-dmzze*(vecnormaly(ip)+vecnormaly(ivze1))/2d0
    vectangyze(ip)=dmzze*(vecnormalx(ip)+vecnormalx(ivze1))/2d0-dmxze*(vecnormalz(ip)+vecnormalz(ivze1))/2d0
    vectangzze(ip)=dmxze*(vecnormaly(ip)+vecnormaly(ivze1))/2d0-dmyze*(vecnormalx(ip)+vecnormalx(ivze1))/2d0
!    
    end do
!
 100 format(A2,1X,F18.13,1x,F18.13,1X,F18.13,1X,A11,3(2x,F9.5))
!
    flush(69)
!        
    return
    end
!
!
!    
    subroutine energia_campo_local
    use VARIAVEIS
    implicit none
    real(kind=8) ::rij3,pui,puj,cuxi,cuyi,cuzi,produ,partialx,partialy,partialz
    real(kind=8) ::dmx1,dmx2,dmx3,dmx4,dmy1,dmy2,dmy3,dmy4,dmxy,dmz1,dmz2,dmz3,dmz4,vizty,vizty1
    real(kind=8) ::dmx5,dmx6,dmy5,dmy6,dmz5,dmz6
    integer :: ivp,invp,ip,iv,ivxd1,ivxe1,ivyd1,ivye1,ivzd1,ivze1,ii
    integer :: isvxd1,isvxe1,isvyd1,isvye1,isvzd1,isvze1,ity
!
	energiatot=0d0
	energcamp=0d0
	energdip=0d0
    partialx=0d0
    partialy=0d0
    partialz=0d0
	dzym=0d0
	troca=0d0
	anisotro=0d0
	usa=0
    dipx=0d0
    dipy=0d0
    dipz=0d0
    invp=0
!
	do ip=1,npt
!    
!!!!!! Campo local de troca, DM e de anisotropia  !!!!!!!
!
!       skyrmion do tipo ouriço
!        
    ivxd1=gvxd1(ip)
    ivxe1=gvxe1(ip)
    ivyd1=gvyd1(ip)
    ivye1=gvye1(ip)
    ivzd1=gvzd1(ip)
    ivze1=gvze1(ip)
!    
    ity=ty(ip)
    isvxd1=ity+ty(ivxd1)
    isvxe1=ity+ty(ivxe1)
    isvyd1=ity+ty(ivyd1)
    isvye1=ity+ty(ivye1)
    isvzd1=ity+ty(ivzd1)
    isvze1=ity+ty(ivze1)
!
    exchangx(ip)=(jotacelxd(ip)*uxi(ivxd1)+jotacelxe(ip)*uxi(ivxe1)+&
    jotacelyd(ip)*uxi(ivyd1)+jotacelye(ip)*uxi(ivye1)+jotacelzd(ip)*uxi(ivzd1)+&
    jotacelze(ip)*uxi(ivze1))
!
    exchangy(ip)=(jotacelxd(ip)*uyi(ivxd1)+jotacelxe(ip)*uyi(ivxe1)+&
    jotacelyd(ip)*uyi(ivyd1)+jotacelye(ip)*uyi(ivye1)+jotacelzd(ip)*uyi(ivzd1)+&
    jotacelze(ip)*uyi(ivze1))
!
    exchangz(ip)=(jotacelxd(ip)*uzi(ivxd1)+jotacelxe(ip)*uzi(ivxe1)+&
    jotacelyd(ip)*uzi(ivyd1)+jotacelye(ip)*uzi(ivye1)+jotacelzd(ip)*uzi(ivzd1)+&
    jotacelze(ip)*uzi(ivze1))
!
    if(dm.ne.0) then
!
    SELECT CASE (tipo_dm)
!
    case (1)
!
    dmx1=dmcelxd(ip)*(uyi(ivxd1)*tangzxd(ip)-uzi(ivxd1)*tangyxd(ip))
    dmx2=dmcelxe(ip)*(uyi(ivxe1)*tangzxe(ip)-uzi(ivxe1)*tangyxe(ip))
    dmx3=dmcelyd(ip)*(uyi(ivyd1)*tangzyd(ip)-uzi(ivyd1)*tangyyd(ip))
    dmx4=dmcelye(ip)*(uyi(ivye1)*tangzye(ip)-uzi(ivye1)*tangyye(ip))
    dmx5=dmcelzd(ip)*(uyi(ivzd1)*tangzzd(ip)-uzi(ivzd1)*tangyzd(ip))
    dmx6=dmcelze(ip)*(uyi(ivze1)*tangzze(ip)-uzi(ivze1)*tangyze(ip))
    dzymoriax(ip)=dmx1+dmx2+dmx3+dmx4+dmx5+dmx6
!
    dmy1=dmcelxd(ip)*(uzi(ivxd1)*tangxxd(ip)-uxi(ivxd1)*tangzxd(ip))
    dmy2=dmcelxe(ip)*(uzi(ivxe1)*tangxxe(ip)-uxi(ivxe1)*tangzxe(ip))
    dmy3=dmcelyd(ip)*(uzi(ivyd1)*tangxyd(ip)-uxi(ivyd1)*tangzyd(ip))
    dmy4=dmcelye(ip)*(uzi(ivye1)*tangxye(ip)-uxi(ivye1)*tangzye(ip))
    dmy5=dmcelzd(ip)*(uzi(ivzd1)*tangxzd(ip)-uxi(ivzd1)*tangzzd(ip))
    dmy6=dmcelze(ip)*(uzi(ivze1)*tangxze(ip)-uxi(ivze1)*tangzze(ip))
    dzymoriay(ip)=dmy1+dmy2+dmy3+dmy4+dmy5+dmy6
!
    dmz1=dmcelxd(ip)*(uxi(ivxd1)*tangyxd(ip)-uyi(ivxd1)*tangxxd(ip))
    dmz2=dmcelxe(ip)*(uxi(ivxe1)*tangyxe(ip)-uyi(ivxe1)*tangxxe(ip))
    dmz3=dmcelyd(ip)*(uxi(ivyd1)*tangyyd(ip)-uyi(ivyd1)*tangxyd(ip))
    dmz4=dmcelye(ip)*(uxi(ivye1)*tangyye(ip)-uyi(ivye1)*tangxye(ip))
    dmz5=dmcelzd(ip)*(uxi(ivzd1)*tangyzd(ip)-uyi(ivzd1)*tangxzd(ip))
    dmz6=dmcelze(ip)*(uxi(ivze1)*tangyze(ip)-uyi(ivze1)*tangxze(ip))
    dzymoriaz(ip)=dmz1+dmz2+dmz3+dmz4+dmz5+dmz6
!
    case(2)
!
    dmx1=dmcelxd(ip)*(uyi(ivxd1)*tangzxd(ip)-uzi(ivxd1)*tangyxd(ip))
    dmx2=dmcelxe(ip)*(uyi(ivxe1)*tangzxe(ip)-uzi(ivxe1)*tangyxe(ip))
    dmx3=dmcelyd(ip)*(uyi(ivyd1)*tangzyd(ip)-uzi(ivyd1)*tangyyd(ip))
    dmx4=dmcelye(ip)*(uyi(ivye1)*tangzye(ip)-uzi(ivye1)*tangyye(ip))
    dmx5=dmcelzd(ip)*(uyi(ivzd1)*tangzzd(ip)-uzi(ivzd1)*tangyzd(ip))
    dmx6=dmcelze(ip)*(uyi(ivze1)*tangzze(ip)-uzi(ivze1)*tangyze(ip))
    dzymoriax(ip)=dmx1+dmx2+dmx3+dmx4+dmx5+dmx6
!
    dmy1=dmcelxd(ip)*(-uzi(ivxd1)*tangxxd(ip)-uxi(ivxd1)*tangzxd(ip))
    dmy2=dmcelxe(ip)*(-uzi(ivxe1)*tangxxe(ip)-uxi(ivxe1)*tangzxe(ip))
    dmy3=dmcelyd(ip)*(-uzi(ivyd1)*tangxyd(ip)-uxi(ivyd1)*tangzyd(ip))
    dmy4=dmcelye(ip)*(-uzi(ivye1)*tangxye(ip)-uxi(ivye1)*tangzye(ip))
    dmy5=dmcelzd(ip)*(-uzi(ivzd1)*tangxzd(ip)-uxi(ivzd1)*tangzzd(ip))
    dmy6=dmcelze(ip)*(-uzi(ivze1)*tangxze(ip)-uxi(ivze1)*tangzze(ip))
    dzymoriay(ip)=dmy1+dmy2+dmy3+dmy4+dmy5+dmy6
!
    dmz1=dmcelxd(ip)*(uxi(ivxd1)*tangyxd(ip)+uyi(ivxd1)*tangxxd(ip))
    dmz2=dmcelxe(ip)*(uxi(ivxe1)*tangyxe(ip)+uyi(ivxe1)*tangxxe(ip))
    dmz3=dmcelyd(ip)*(uxi(ivyd1)*tangyyd(ip)+uyi(ivyd1)*tangxyd(ip))
    dmz4=dmcelye(ip)*(uxi(ivye1)*tangyye(ip)+uyi(ivye1)*tangxye(ip))
    dmz5=dmcelzd(ip)*(uxi(ivzd1)*tangyzd(ip)+uyi(ivzd1)*tangxzd(ip))
    dmz6=dmcelze(ip)*(uxi(ivze1)*tangyze(ip)+uyi(ivze1)*tangxze(ip))
    dzymoriaz(ip)=dmz1+dmz2+dmz3+dmz4+dmz5+dmz6
!
    case (0)
!
    dmx1=dmcelxd(ip)*(uyi(ivxd1)*vectangzxd(ip)-uzi(ivxd1)*vectangyxd(ip))
    dmx2=dmcelxe(ip)*(uyi(ivxe1)*vectangzxe(ip)-uzi(ivxe1)*vectangyxe(ip))
    dmx3=dmcelyd(ip)*(uyi(ivyd1)*vectangzyd(ip)-uzi(ivyd1)*vectangyyd(ip))
    dmx4=dmcelye(ip)*(uyi(ivye1)*vectangzye(ip)-uzi(ivye1)*vectangyye(ip))
    dzymoriax(ip)=dmx1+dmx2+dmx3+dmx4
!
    dmy1=dmcelxd(ip)*(uzi(ivxd1)*vectangxxd(ip)-uxi(ivxd1)*vectangzxd(ip))
    dmy2=dmcelxe(ip)*(uzi(ivxe1)*vectangxxe(ip)-uxi(ivxe1)*vectangzxe(ip))
    dmy3=dmcelyd(ip)*(uzi(ivyd1)*vectangxyd(ip)-uxi(ivyd1)*vectangzyd(ip))
    dmy4=dmcelye(ip)*(uzi(ivye1)*vectangxye(ip)-uxi(ivye1)*vectangzye(ip))
    dzymoriay(ip)=dmy1+dmy2+dmy3+dmy4
!
    dmz1=dmcelxd(ip)*(uxi(ivxd1)*vectangyxd(ip)-uyi(ivxd1)*vectangxxd(ip))
    dmz2=dmcelxe(ip)*(uxi(ivxe1)*vectangyxe(ip)-uyi(ivxe1)*vectangxxe(ip))
    dmz3=dmcelyd(ip)*(uxi(ivyd1)*vectangyyd(ip)-uyi(ivyd1)*vectangxyd(ip))
    dmz4=dmcelye(ip)*(uxi(ivye1)*vectangyye(ip)-uyi(ivye1)*vectangxye(ip))
    dzymoriaz(ip)=dmz1+dmz2+dmz3+dmz4
!
    end select
!
    end if
!
    if(aniso.ne.0) then
!
	produ=uxi(ip)*vecnormalx(ip)+uzi(ip)*vecnormalz(ip)+uyi(ip)*vecnormaly(ip)
    anisotroz(ip)=anisovol(ip)*produ*vecnormalz(ip)
	anisotrox(ip)=anisovol(ip)*produ*vecnormalx(ip)
	anisotroy(ip)=anisovol(ip)*produ*vecnormaly(ip)
!
    end if
!
    cuxi=uxi(ip)
    cuyi=uyi(ip)
    cuzi=uzi(ip)
    vizty1=magty(ity)    
!    
    do ivp=ip+1,npt
!    
    invp=invp+1

!!!! Campo dipolar e energia dipolar !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
        rij3=rvij3(invp)
        vizty=magty(ty(ivp))
!
    puj=3.0d0*(uvrijx(invp)*uxi(ivp)+uvrijy(invp)*uyi(ivp)+uvrijz(invp)*uzi(ivp))
    pui=3.0d0*(uvrijx(invp)*cuxi+uvrijy(invp)*cuyi+uvrijz(invp)*cuzi)
!
 	dipx(ip)=dipx(ip)+vizty*dipvol(ivp)*(puj*uvrijx(invp)-uxi(ivp))*rij3
    dipy(ip)=dipy(ip)+vizty*dipvol(ivp)*(puj*uvrijy(invp)-uyi(ivp))*rij3
    dipz(ip)=dipz(ip)+vizty*dipvol(ivp)*(puj*uvrijz(invp)-uzi(ivp))*rij3

!
    dipx(ivp)=dipx(ivp)+vizty1*dipvol(ip)*(pui*uvrijx(invp)-cuxi)*rij3
    dipy(ivp)=dipy(ivp)+vizty1*dipvol(ip)*(pui*uvrijy(invp)-cuyi)*rij3
    dipz(ivp)=dipz(ivp)+vizty1*dipvol(ip)*(pui*uvrijz(invp)-cuzi)*rij3
!
 	end do
!
    exchangx(ip)=exchangx(ip)/magsavol(ip)
    exchangy(ip)=exchangy(ip)/magsavol(ip)
    exchangz(ip)=exchangz(ip)/magsavol(ip)
!
    dzymoriax(ip)=dzymoriax(ip)/magsavol(ip)
    dzymoriay(ip)=dzymoriay(ip)/magsavol(ip)
    dzymoriaz(ip)=dzymoriaz(ip)/magsavol(ip)
!
    if(magty(ity).ne.0) then
    partialx=(dzymoriax(ip)+exchangx(ip))/magty(ity)
    partialy=(dzymoriay(ip)+exchangy(ip))/magty(ity)
    partialz=(dzymoriaz(ip)+exchangz(ip))/magty(ity)
    end if
!!!!!! Campo efetivo local !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    campo_efetivo_x(ip)=(partialx+dipx(ip)+anisotrox(ip)+camp_bx)
    campo_efetivo_y(ip)=(partialy+dipy(ip)+anisotroy(ip)+camp_by)
    campo_efetivo_z(ip)=(partialz+dipz(ip)+anisotroz(ip)+camp_bz)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!! energia DM, troca, anisotropia !!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    dzym=dzym-magsavol(ip)*(uxi(ip)*dzymoriax(ip)+uyi(ip)*dzymoriay(ip)+uzi(ip)*dzymoriaz(ip))
	troca=troca-magsavol(ip)*(uxi(ip)*exchangx(ip)+uyi(ip)*exchangy(ip)+uzi(ip)*exchangz(ip))
	anisotro=anisotro-magty(ity)*magsavol(ip)*anisovol(ip)*produ**2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! calculo da energia dipolar e de Zeeman!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    energdip=energdip-magsavol(ip)*magty(ity)*(uxi(ip)*dipx(ip)+uyi(ip)*dipy(ip)+uzi(ip)*dipz(ip))
    energcamp=energcamp-magty(ity)*magsavol(ip)*(camp_bx*uxi(ip)+camp_by*uyi(ip)+camp_bz*uzi(ip))
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    end do
!
    energdip=energdip*umfento/2d0
    dzym=dzym*umfento/2d0
    troca=troca*umfento/2d0
    anisotro=anisotro*umfento/2d0
    energcamp=energcamp*umfento
    energiatot=(energdip+dzym+troca+anisotro+energcamp)
!    write(*,*)energiatot,dzym,troca,anisotro
! 
	return
	end
!
!
!
    subroutine energia_campo_local_corte
    use VARIAVEIS
    implicit none
    real(kind=8) ::rij3,pui,puj,cuxi,cuyi,cuzi,produ,partialx,partialy,partialz
    real(kind=8) ::dmx1,dmx2,dmx3,dmx4,dmy1,dmy2,dmy3,dmy4,dmz1,dmz2,dmz3,dmz4,vizty,vizty1
    real(kind=8) ::dmx5,dmx6,dmy5,dmy6,dmz5,dmz6
    integer :: ivp,invp,ip,iv,ivxd1,ivxe1,ivyd1,ivye1,ivzd1,ivze1,ii
    integer :: isvxd1,isvxe1,isvyd1,isvye1,isvzd1,isvze1,ity
!
	energiatot=0d0
	energdip=0d0
	energcamp=0d0
    partialx=0d0
    partialy=0d0
    partialz=0d0
	dzym=0d0
	troca=0d0
	anisotro=0d0
	usa=0
    dipx=0d0
    dipy=0d0
    dipz=0d0
    invp=0
!
	do ip=1,npt
!    
!!!!!! Campo local de troca, DM e de anisotropia  !!!!!!!
!
!       skyrmion do tipo ouriço
!        
    ivxd1=gvxd1(ip)
    ivxe1=gvxe1(ip)
    ivyd1=gvyd1(ip)
    ivye1=gvye1(ip)
    ivzd1=gvzd1(ip)
    ivze1=gvze1(ip)
!
    ity=ty(ip)
    isvxd1=ity+ty(ivxd1)
    isvxe1=ity+ty(ivxe1)
    isvyd1=ity+ty(ivyd1)
    isvye1=ity+ty(ivye1)
    isvzd1=ity+ty(ivzd1)
    isvze1=ity+ty(ivze1)
!
    exchangx(ip)=(jotacelxd(ip)*uxi(ivxd1)+jotacelxe(ip)*uxi(ivxe1)+&
    jotacelyd(ip)*uxi(ivyd1)+jotacelye(ip)*uxi(ivye1)+jotacelzd(ip)*uxi(ivzd1)+&
    jotacelze(ip)*uxi(ivze1))
!
    exchangy(ip)=(jotacelxd(ip)*uyi(ivxd1)+jotacelxe(ip)*uyi(ivxe1)+&
    jotacelyd(ip)*uyi(ivyd1)+jotacelye(ip)*uyi(ivye1)+jotacelzd(ip)*uyi(ivzd1)+&
    jotacelze(ip)*uyi(ivze1))
!
    exchangz(ip)=(jotacelxd(ip)*uzi(ivxd1)+jotacelxe(ip)*uzi(ivxe1)+&
    jotacelyd(ip)*uzi(ivyd1)+jotacelye(ip)*uzi(ivye1)+jotacelzd(ip)*uzi(ivzd1)+&
    jotacelze(ip)*uzi(ivze1))
!
    if(dm.ne.0) then
!
    SELECT CASE (tipo_dm)
!
    case (1)
!
    dmx1=dmcelxd(ip)*(uyi(ivxd1)*tangzxd(ip)-uzi(ivxd1)*tangyxd(ip))
    dmx2=dmcelxe(ip)*(uyi(ivxe1)*tangzxe(ip)-uzi(ivxe1)*tangyxe(ip))
    dmx3=dmcelyd(ip)*(uyi(ivyd1)*tangzyd(ip)-uzi(ivyd1)*tangyyd(ip))
    dmx4=dmcelye(ip)*(uyi(ivye1)*tangzye(ip)-uzi(ivye1)*tangyye(ip))
    dmx5=dmcelzd(ip)*(uyi(ivzd1)*tangzzd(ip)-uzi(ivzd1)*tangyzd(ip))
    dmx6=dmcelze(ip)*(uyi(ivze1)*tangzze(ip)-uzi(ivze1)*tangyze(ip))
    dzymoriax(ip)=dmx1+dmx2+dmx3+dmx4+dmx5+dmx6
!
    dmy1=dmcelxd(ip)*(uzi(ivxd1)*tangxxd(ip)-uxi(ivxd1)*tangzxd(ip))
    dmy2=dmcelxe(ip)*(uzi(ivxe1)*tangxxe(ip)-uxi(ivxe1)*tangzxe(ip))
    dmy3=dmcelyd(ip)*(uzi(ivyd1)*tangxyd(ip)-uxi(ivyd1)*tangzyd(ip))
    dmy4=dmcelye(ip)*(uzi(ivye1)*tangxye(ip)-uxi(ivye1)*tangzye(ip))
    dmy5=dmcelzd(ip)*(uzi(ivzd1)*tangxzd(ip)-uxi(ivzd1)*tangzzd(ip))
    dmy6=dmcelze(ip)*(uzi(ivze1)*tangxze(ip)-uxi(ivze1)*tangzze(ip))
    dzymoriay(ip)=dmy1+dmy2+dmy3+dmy4+dmy5+dmy6
!
    dmz1=dmcelxd(ip)*(uxi(ivxd1)*tangyxd(ip)-uyi(ivxd1)*tangxxd(ip))
    dmz2=dmcelxe(ip)*(uxi(ivxe1)*tangyxe(ip)-uyi(ivxe1)*tangxxe(ip))
    dmz3=dmcelyd(ip)*(uxi(ivyd1)*tangyyd(ip)-uyi(ivyd1)*tangxyd(ip))
    dmz4=dmcelye(ip)*(uxi(ivye1)*tangyye(ip)-uyi(ivye1)*tangxye(ip))
    dmz5=dmcelzd(ip)*(uxi(ivzd1)*tangyzd(ip)-uyi(ivzd1)*tangxzd(ip))
    dmz6=dmcelze(ip)*(uxi(ivze1)*tangyze(ip)-uyi(ivze1)*tangxze(ip))
    dzymoriaz(ip)=dmz1+dmz2+dmz3+dmz4+dmz5+dmz6
!
    case(2)
!
    dmx1=dmcelxd(ip)*(uyi(ivxd1)*tangzxd(ip)-uzi(ivxd1)*tangyxd(ip))
    dmx2=dmcelxe(ip)*(uyi(ivxe1)*tangzxe(ip)-uzi(ivxe1)*tangyxe(ip))
    dmx3=dmcelyd(ip)*(uyi(ivyd1)*tangzyd(ip)-uzi(ivyd1)*tangyyd(ip))
    dmx4=dmcelye(ip)*(uyi(ivye1)*tangzye(ip)-uzi(ivye1)*tangyye(ip))
    dmx5=dmcelzd(ip)*(uyi(ivzd1)*tangzzd(ip)-uzi(ivzd1)*tangyzd(ip))
    dmx6=dmcelze(ip)*(uyi(ivze1)*tangzze(ip)-uzi(ivze1)*tangyze(ip))
    dzymoriax(ip)=dmx1+dmx2+dmx3+dmx4+dmx5+dmx6
!
    dmy1=dmcelxd(ip)*(-uzi(ivxd1)*tangxxd(ip)-uxi(ivxd1)*tangzxd(ip))
    dmy2=dmcelxe(ip)*(-uzi(ivxe1)*tangxxe(ip)-uxi(ivxe1)*tangzxe(ip))
    dmy3=dmcelyd(ip)*(-uzi(ivyd1)*tangxyd(ip)-uxi(ivyd1)*tangzyd(ip))
    dmy4=dmcelye(ip)*(-uzi(ivye1)*tangxye(ip)-uxi(ivye1)*tangzye(ip))
    dmy5=dmcelzd(ip)*(-uzi(ivzd1)*tangxzd(ip)-uxi(ivzd1)*tangzzd(ip))
    dmy6=dmcelze(ip)*(-uzi(ivze1)*tangxze(ip)-uxi(ivze1)*tangzze(ip))
    dzymoriay(ip)=dmy1+dmy2+dmy3+dmy4+dmy5+dmy6
!
    dmz1=dmcelxd(ip)*(uxi(ivxd1)*tangyxd(ip)+uyi(ivxd1)*tangxxd(ip))
    dmz2=dmcelxe(ip)*(uxi(ivxe1)*tangyxe(ip)+uyi(ivxe1)*tangxxe(ip))
    dmz3=dmcelyd(ip)*(uxi(ivyd1)*tangyyd(ip)+uyi(ivyd1)*tangxyd(ip))
    dmz4=dmcelye(ip)*(uxi(ivye1)*tangyye(ip)+uyi(ivye1)*tangxye(ip))
    dmz5=dmcelzd(ip)*(uxi(ivzd1)*tangyzd(ip)+uyi(ivzd1)*tangxzd(ip))
    dmz6=dmcelze(ip)*(uxi(ivze1)*tangyze(ip)+uyi(ivze1)*tangxze(ip))
    dzymoriaz(ip)=dmz1+dmz2+dmz3+dmz4+dmz5+dmz6
!
    case (0)
!
    dmx1=dmcelxd(ip)*(uyi(ivxd1)*vectangzxd(ip)-uzi(ivxd1)*vectangyxd(ip))
    dmx2=dmcelxe(ip)*(uyi(ivxe1)*vectangzxe(ip)-uzi(ivxe1)*vectangyxe(ip))
    dmx3=dmcelyd(ip)*(uyi(ivyd1)*vectangzyd(ip)-uzi(ivyd1)*vectangyyd(ip))
    dmx4=dmcelye(ip)*(uyi(ivye1)*vectangzye(ip)-uzi(ivye1)*vectangyye(ip))
    dzymoriax(ip)=dmx1+dmx2+dmx3+dmx4
!
    dmy1=dmcelxd(ip)*(uzi(ivxd1)*vectangxxd(ip)-uxi(ivxd1)*vectangzxd(ip))
    dmy2=dmcelxe(ip)*(uzi(ivxe1)*vectangxxe(ip)-uxi(ivxe1)*vectangzxe(ip))
    dmy3=dmcelyd(ip)*(uzi(ivyd1)*vectangxyd(ip)-uxi(ivyd1)*vectangzyd(ip))
    dmy4=dmcelye(ip)*(uzi(ivye1)*vectangxye(ip)-uxi(ivye1)*vectangzye(ip))
    dzymoriay(ip)=dmy1+dmy2+dmy3+dmy4
!
    dmz1=dmcelxd(ip)*(uxi(ivxd1)*vectangyxd(ip)-uyi(ivxd1)*vectangxxd(ip))
    dmz2=dmcelxe(ip)*(uxi(ivxe1)*vectangyxe(ip)-uyi(ivxe1)*vectangxxe(ip))
    dmz3=dmcelyd(ip)*(uxi(ivyd1)*vectangyyd(ip)-uyi(ivyd1)*vectangxyd(ip))
    dmz4=dmcelye(ip)*(uxi(ivye1)*vectangyye(ip)-uyi(ivye1)*vectangxye(ip))
    dzymoriaz(ip)=dmz1+dmz2+dmz3+dmz4
!
    end select
!
    end if
!
    if(aniso.ne.0) then
!
    produ=uxi(ip)*vecnormalx(ip)+uzi(ip)*vecnormalz(ip)+uyi(ip)*vecnormaly(ip)
    anisotroz(ip)=anisovol(ip)*produ*vecnormalz(ip)
	anisotrox(ip)=anisovol(ip)*produ*vecnormalx(ip)
	anisotroy(ip)=anisovol(ip)*produ*vecnormaly(ip)
!
    end if
!
    cuxi=uxi(ip)
    cuyi=uyi(ip)
    cuzi=uzi(ip)
    vizty1=magty(ity)
!    
    do invp=nvtp(ip-1)+1,nvtp(ip)
!
        ivp=gvv(invp)
        vizty=magty(ty(ivp))

!!!! Campo dipolar e energia dipolar !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       
        rij3=rvij3(invp)
!
    puj=3.0d0*(uvrijx(invp)*uxi(ivp)+uvrijy(invp)*uyi(ivp)+uvrijz(invp)*uzi(ivp))
    pui=3.0d0*(uvrijx(invp)*cuxi+uvrijy(invp)*cuyi+uvrijz(invp)*cuzi)
!
 	dipx(ip)=dipx(ip)+vizty*dipvol(ivp)*(puj*uvrijx(invp)-uxi(ivp))*rij3
    dipy(ip)=dipy(ip)+vizty*dipvol(ivp)*(puj*uvrijy(invp)-uyi(ivp))*rij3
    dipz(ip)=dipz(ip)+vizty*dipvol(ivp)*(puj*uvrijz(invp)-uzi(ivp))*rij3

!        
    dipx(ivp)=dipx(ivp)+vizty1*dipvol(ip)*(pui*uvrijx(invp)-cuxi)*rij3
    dipy(ivp)=dipy(ivp)+vizty1*dipvol(ip)*(pui*uvrijy(invp)-cuyi)*rij3
    dipz(ivp)=dipz(ivp)+vizty1*dipvol(ip)*(pui*uvrijz(invp)-cuzi)*rij3

!
 	end do	
!
    exchangx(ip)=exchangx(ip)/magsavol(ip)
    exchangy(ip)=exchangy(ip)/magsavol(ip)
    exchangz(ip)=exchangz(ip)/magsavol(ip)
!
    dzymoriax(ip)=dzymoriax(ip)/magsavol(ip)
    dzymoriay(ip)=dzymoriay(ip)/magsavol(ip)
    dzymoriaz(ip)=dzymoriaz(ip)/magsavol(ip)
!
    if(magty(ity).ne.0) then
    partialx=(dzymoriax(ip)+exchangx(ip))/magty(ity)
    partialy=(dzymoriay(ip)+exchangy(ip))/magty(ity)
    partialz=(dzymoriaz(ip)+exchangz(ip))/magty(ity)
    end if
!!!!!! Campo efetivo local !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    campo_efetivo_x(ip)=(partialx+dipx(ip)+anisotrox(ip)+camp_bx)
    campo_efetivo_y(ip)=(partialy+dipy(ip)+anisotroy(ip)+camp_by)
    campo_efetivo_z(ip)=(partialz+dipz(ip)+anisotroz(ip)+camp_bz)
!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!! energia DM, troca, anisotropia !!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    dzym=dzym-magsavol(ip)*(uxi(ip)*dzymoriax(ip)+uyi(ip)*dzymoriay(ip)+uzi(ip)*dzymoriaz(ip))
	troca=troca-magsavol(ip)*(uxi(ip)*exchangx(ip)+uyi(ip)*exchangy(ip)+uzi(ip)*exchangz(ip))
	anisotro=anisotro-magty(ity)*magsavol(ip)*anisovol(ip)*produ**2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! calculo da energia dipolar e de Zeeman!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    energdip=energdip-magsavol(ip)*magty(ity)*(uxi(ip)*dipx(ip)+uyi(ip)*dipy(ip)+uzi(ip)*dipz(ip))
    energcamp=energcamp-magty(ity)*magsavol(ip)*(camp_bx*uxi(ip)+camp_by*uyi(ip)+camp_bz*uzi(ip))
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    end do
!
    energdip=energdip*umfento/2d0
    dzym=dzym*umfento/2d0
    troca=troca*umfento/2d0
    anisotro=anisotro*umfento/2d0
    energcamp=energcamp*umfento
    energiatot=(energdip+dzym+troca+anisotro+energcamp)
!    write(*,*)energiatot,dzym,troca,anisotro,energdip
! 
	return
	end
!
!
!
    subroutine runge_kuta(i,ij)
    use VARIAVEIS
    integer :: ip,i,ij
    real(kind=8) :: mbx,mby,mbz,soma
!
    if(ij.eq.0) then
    uxi0=uxi
    uyi0=uyi
    uzi0=uzi
    end if
!
    do ip=1,npt
!
    if(ij.eq.0) then

!!!!!!! Equação LLG dm/dt para ser integrada !!!!!!!!!!!!!!!!!!!!!

    mbx=uyi(ip)*campo_efetivo_z(ip)-uzi(ip)*campo_efetivo_y(ip)
    mby=uzi(ip)*campo_efetivo_x(ip)-uxi(ip)*campo_efetivo_z(ip)
    mbz=uxi(ip)*campo_efetivo_y(ip)-uyi(ip)*campo_efetivo_x(ip)

    k1x(ip)=gama_alfa*(mbx+alfa*(uyi(ip)*mbz-uzi(ip)*mby))
    k1y(ip)=gama_alfa*(mby+alfa*(uzi(ip)*mbx-uxi(ip)*mbz))
    k1z(ip)=gama_alfa*(mbz+alfa*(uxi(ip)*mby-uyi(ip)*mbx))
!
    if(i.eq.1) then
    md3x(ip)=k1x(ip)
    md3y(ip)=k1y(ip)
    md3z(ip)=k1z(ip)
    end if
    if(i.eq.2) then
    md2x(ip)=k1x(ip)
    md2y(ip)=k1y(ip)
    md2z(ip)=k1z(ip)
    end if
    if(i.eq.3) then
    md1x(ip)=k1x(ip)
    md1y(ip)=k1y(ip)
    md1z(ip)=k1z(ip)
    end if

    uxi(ip)=uxi0(ip)+k1x(ip)*dt/2d0
    uyi(ip)=uyi0(ip)+k1y(ip)*dt/2d0
    uzi(ip)=uzi0(ip)+k1z(ip)*dt/2d0
!
    end if
!
    if(ij.eq.1) then
!
    mbx=uyi(ip)*campo_efetivo_z(ip)-uzi(ip)*campo_efetivo_y(ip)
    mby=uzi(ip)*campo_efetivo_x(ip)-uxi(ip)*campo_efetivo_z(ip)
    mbz=uxi(ip)*campo_efetivo_y(ip)-uyi(ip)*campo_efetivo_x(ip)

    k2x(ip)=gama_alfa*(mbx+alfa*(uyi(ip)*mbz-uzi(ip)*mby))
    k2y(ip)=gama_alfa*(mby+alfa*(uzi(ip)*mbx-uxi(ip)*mbz))
    k2z(ip)=gama_alfa*(mbz+alfa*(uxi(ip)*mby-uyi(ip)*mbx))

    uxi(ip)=uxi0(ip)+k2x(ip)*dt/2d0
    uyi(ip)=uyi0(ip)+k2y(ip)*dt/2d0
    uzi(ip)=uzi0(ip)+k2z(ip)*dt/2d0
!
    end if
!
    if(ij.eq.2) then
!
    mbx=uyi(ip)*campo_efetivo_z(ip)-uzi(ip)*campo_efetivo_y(ip)
    mby=uzi(ip)*campo_efetivo_x(ip)-uxi(ip)*campo_efetivo_z(ip)
    mbz=uxi(ip)*campo_efetivo_y(ip)-uyi(ip)*campo_efetivo_x(ip)

    k3x(ip)=gama_alfa*(mbx+alfa*(uyi(ip)*mbz-uzi(ip)*mby))
    k3y(ip)=gama_alfa*(mby+alfa*(uzi(ip)*mbx-uxi(ip)*mbz))
    k3z(ip)=gama_alfa*(mbz+alfa*(uxi(ip)*mby-uyi(ip)*mbx))

    uxi(ip)=uxi0(ip)+k3x(ip)*dt
    uyi(ip)=uyi0(ip)+k3y(ip)*dt
    uzi(ip)=uzi0(ip)+k3z(ip)*dt
!
    end if
!
    if(ij.eq.3) then
!
    mbx=uyi(ip)*campo_efetivo_z(ip)-uzi(ip)*campo_efetivo_y(ip)
    mby=uzi(ip)*campo_efetivo_x(ip)-uxi(ip)*campo_efetivo_z(ip)
    mbz=uxi(ip)*campo_efetivo_y(ip)-uyi(ip)*campo_efetivo_x(ip)

    k4x(ip)=gama_alfa*(mbx+alfa*(uyi(ip)*mbz-uzi(ip)*mby))
    k4y(ip)=gama_alfa*(mby+alfa*(uzi(ip)*mbx-uxi(ip)*mbz))
    k4z(ip)=gama_alfa*(mbz+alfa*(uxi(ip)*mby-uyi(ip)*mbx))

    uxi(ip)=uxi0(ip)+(k1x(ip)+2d0*k2x(ip)+2d0*k3x(ip)+k4x(ip))*dt/6d0
    uyi(ip)=uyi0(ip)+(k1y(ip)+2d0*k2y(ip)+2d0*k3y(ip)+k4y(ip))*dt/6d0
    uzi(ip)=uzi0(ip)+(k1z(ip)+2d0*k2z(ip)+2d0*k3z(ip)+k4z(ip))*dt/6d0
!
    end if
!
    end do
!    
    return
    end
!
    subroutine runge_kuta_flag_cor0(i,ij)
    use VARIAVEIS
    integer :: ip,i,ij
    real(kind=8) :: mbx,mby,mbz,soma,p1vx,p1vy,p1vz,p2vx,p2vy,p2vz,p3vx,p3vy,p3vz
    real(kind=8) :: termo_correntex,termo_correntey,termo_correntez,vizty
!
    if(ij.eq.0) then
    uxi0=uxi
    uyi0=uyi
    uzi0=uzi
    end if
!
    do ip=1,npt
!
    vizty=magty(ty(ip))
    if(vizty.eq.0) vizty=1
!
    if(ij.eq.0) then

!!!!!!! Equação LLG dm/dt para ser integrada !!!!!!!!!!!!!!!!!!!!!

    mbx=uyi(ip)*campo_efetivo_z(ip)-uzi(ip)*campo_efetivo_y(ip)
    mby=uzi(ip)*campo_efetivo_x(ip)-uxi(ip)*campo_efetivo_z(ip)
    mbz=uxi(ip)*campo_efetivo_y(ip)-uyi(ip)*campo_efetivo_x(ip)
!
    ivxd1=gvxd1(ip)
    ivxe1=gvxe1(ip)
!
    p1vx=dist_x(ip)*(uyi(ip)*(uzi(ivxd1)-uzi(ivxe1))-uzi(ip)*(uyi(ivxd1)-uyi(ivxe1)))
    p1vy=dist_x(ip)*(uzi(ip)*(uxi(ivxd1)-uxi(ivxe1))-uxi(ip)*(uzi(ivxd1)-uzi(ivxe1)))
    p1vz=dist_x(ip)*(uxi(ip)*(uyi(ivxd1)-uyi(ivxe1))-uyi(ip)*(uxi(ivxd1)-uxi(ivxe1)))
!
    p2vx=uyi(ip)*p1vz-uzi(ip)*p1vy
    p2vy=uzi(ip)*p1vx-uxi(ip)*p1vz
    p2vz=uxi(ip)*p1vy-uyi(ip)*p1vx
!
    p3vx=uyi(ip)*p2vz-uzi(ip)*p2vy
    p3vy=uzi(ip)*p2vx-uxi(ip)*p2vz
    p3vz=uxi(ip)*p2vy-uyi(ip)*p2vx
!
    termo_correntex=gama_alfa1*(vj_qui1*p1vx+vj_qui_alfa2*p2vx+vj_alfa3*p3vx)
    termo_correntey=gama_alfa1*(vj_qui1*p1vy+vj_qui_alfa2*p2vy+vj_alfa3*p3vy)
    termo_correntez=gama_alfa1*(vj_qui1*p1vz+vj_qui_alfa2*p2vz+vj_alfa3*p3vz)
!
    k1x(ip)=gama_alfa*(mbx+alfa*(uyi(ip)*mbz-uzi(ip)*mby))+termo_correntex/vizty
    k1y(ip)=gama_alfa*(mby+alfa*(uzi(ip)*mbx-uxi(ip)*mbz))+termo_correntey/vizty
    k1z(ip)=gama_alfa*(mbz+alfa*(uxi(ip)*mby-uyi(ip)*mbx))+termo_correntez/vizty
!
    uxi(ip)=uxi0(ip)+k1x(ip)*dt/2d0
    uyi(ip)=uyi0(ip)+k1y(ip)*dt/2d0
    uzi(ip)=uzi0(ip)+k1z(ip)*dt/2d0
!
    end if
!
    if(ij.eq.1) then
!
    mbx=uyi(ip)*campo_efetivo_z(ip)-uzi(ip)*campo_efetivo_y(ip)
    mby=uzi(ip)*campo_efetivo_x(ip)-uxi(ip)*campo_efetivo_z(ip)
    mbz=uxi(ip)*campo_efetivo_y(ip)-uyi(ip)*campo_efetivo_x(ip)
!
    ivxd1=gvxd1(ip)
    ivxe1=gvxe1(ip)
!
    p1vx=dist_x(ip)*(uyi(ip)*(uzi(ivxd1)-uzi(ivxe1))-uzi(ip)*(uyi(ivxd1)-uyi(ivxe1)))
    p1vy=dist_x(ip)*(uzi(ip)*(uxi(ivxd1)-uxi(ivxe1))-uxi(ip)*(uzi(ivxd1)-uzi(ivxe1)))
    p1vz=dist_x(ip)*(uxi(ip)*(uyi(ivxd1)-uyi(ivxe1))-uyi(ip)*(uxi(ivxd1)-uxi(ivxe1)))
!
    p2vx=uyi(ip)*p1vz-uzi(ip)*p1vy
    p2vy=uzi(ip)*p1vx-uxi(ip)*p1vz
    p2vz=uxi(ip)*p1vy-uyi(ip)*p1vx
!
    p3vx=uyi(ip)*p2vz-uzi(ip)*p2vy
    p3vy=uzi(ip)*p2vx-uxi(ip)*p2vz
    p3vz=uxi(ip)*p2vy-uyi(ip)*p2vx
!
    termo_correntex=gama_alfa1*(vj_qui1*p1vx+vj_qui_alfa2*p2vx+vj_alfa3*p3vx)
    termo_correntey=gama_alfa1*(vj_qui1*p1vy+vj_qui_alfa2*p2vy+vj_alfa3*p3vy)
    termo_correntez=gama_alfa1*(vj_qui1*p1vz+vj_qui_alfa2*p2vz+vj_alfa3*p3vz)
!
    k2x(ip)=gama_alfa*(mbx+alfa*(uyi(ip)*mbz-uzi(ip)*mby))+termo_correntex/vizty
    k2y(ip)=gama_alfa*(mby+alfa*(uzi(ip)*mbx-uxi(ip)*mbz))+termo_correntey/vizty
    k2z(ip)=gama_alfa*(mbz+alfa*(uxi(ip)*mby-uyi(ip)*mbx))+termo_correntez/vizty
!
    uxi(ip)=uxi0(ip)+k2x(ip)*dt/2d0
    uyi(ip)=uyi0(ip)+k2y(ip)*dt/2d0
    uzi(ip)=uzi0(ip)+k2z(ip)*dt/2d0
!
    end if
!
    if(ij.eq.2) then
!
    mbx=uyi(ip)*campo_efetivo_z(ip)-uzi(ip)*campo_efetivo_y(ip)
    mby=uzi(ip)*campo_efetivo_x(ip)-uxi(ip)*campo_efetivo_z(ip)
    mbz=uxi(ip)*campo_efetivo_y(ip)-uyi(ip)*campo_efetivo_x(ip)
!
    ivxd1=gvxd1(ip)
    ivxe1=gvxe1(ip)
!
    p1vx=dist_x(ip)*(uyi(ip)*(uzi(ivxd1)-uzi(ivxe1))-uzi(ip)*(uyi(ivxd1)-uyi(ivxe1)))
    p1vy=dist_x(ip)*(uzi(ip)*(uxi(ivxd1)-uxi(ivxe1))-uxi(ip)*(uzi(ivxd1)-uzi(ivxe1)))
    p1vz=dist_x(ip)*(uxi(ip)*(uyi(ivxd1)-uyi(ivxe1))-uyi(ip)*(uxi(ivxd1)-uxi(ivxe1)))
!
    p2vx=uyi(ip)*p1vz-uzi(ip)*p1vy
    p2vy=uzi(ip)*p1vx-uxi(ip)*p1vz
    p2vz=uxi(ip)*p1vy-uyi(ip)*p1vx
!
    p3vx=uyi(ip)*p2vz-uzi(ip)*p2vy
    p3vy=uzi(ip)*p2vx-uxi(ip)*p2vz
    p3vz=uxi(ip)*p2vy-uyi(ip)*p2vx
!
    termo_correntex=gama_alfa1*(vj_qui1*p1vx+vj_qui_alfa2*p2vx+vj_alfa3*p3vx)
    termo_correntey=gama_alfa1*(vj_qui1*p1vy+vj_qui_alfa2*p2vy+vj_alfa3*p3vy)
    termo_correntez=gama_alfa1*(vj_qui1*p1vz+vj_qui_alfa2*p2vz+vj_alfa3*p3vz)
!
    k3x(ip)=gama_alfa*(mbx+alfa*(uyi(ip)*mbz-uzi(ip)*mby))+termo_correntex/vizty
    k3y(ip)=gama_alfa*(mby+alfa*(uzi(ip)*mbx-uxi(ip)*mbz))+termo_correntey/vizty
    k3z(ip)=gama_alfa*(mbz+alfa*(uxi(ip)*mby-uyi(ip)*mbx))+termo_correntez/vizty
!
    uxi(ip)=uxi0(ip)+k3x(ip)*dt
    uyi(ip)=uyi0(ip)+k3y(ip)*dt
    uzi(ip)=uzi0(ip)+k3z(ip)*dt
!
    end if
!
    if(ij.eq.3) then
!
    mbx=uyi(ip)*campo_efetivo_z(ip)-uzi(ip)*campo_efetivo_y(ip)
    mby=uzi(ip)*campo_efetivo_x(ip)-uxi(ip)*campo_efetivo_z(ip)
    mbz=uxi(ip)*campo_efetivo_y(ip)-uyi(ip)*campo_efetivo_x(ip)
!
    ivxd1=gvxd1(ip)
    ivxe1=gvxe1(ip)
!
    p1vx=dist_x(ip)*(uyi(ip)*(uzi(ivxd1)-uzi(ivxe1))-uzi(ip)*(uyi(ivxd1)-uyi(ivxe1)))
    p1vy=dist_x(ip)*(uzi(ip)*(uxi(ivxd1)-uxi(ivxe1))-uxi(ip)*(uzi(ivxd1)-uzi(ivxe1)))
    p1vz=dist_x(ip)*(uxi(ip)*(uyi(ivxd1)-uyi(ivxe1))-uyi(ip)*(uxi(ivxd1)-uxi(ivxe1)))
!
    p2vx=uyi(ip)*p1vz-uzi(ip)*p1vy
    p2vy=uzi(ip)*p1vx-uxi(ip)*p1vz
    p2vz=uxi(ip)*p1vy-uyi(ip)*p1vx
!
    p3vx=uyi(ip)*p2vz-uzi(ip)*p2vy
    p3vy=uzi(ip)*p2vx-uxi(ip)*p2vz
    p3vz=uxi(ip)*p2vy-uyi(ip)*p2vx
!
    termo_correntex=gama_alfa1*(vj_qui1*p1vx+vj_qui_alfa2*p2vx+vj_alfa3*p3vx)
    termo_correntey=gama_alfa1*(vj_qui1*p1vy+vj_qui_alfa2*p2vy+vj_alfa3*p3vy)
    termo_correntez=gama_alfa1*(vj_qui1*p1vz+vj_qui_alfa2*p2vz+vj_alfa3*p3vz)
!
    k4x(ip)=gama_alfa*(mbx+alfa*(uyi(ip)*mbz-uzi(ip)*mby))+termo_correntex/vizty
    k4y(ip)=gama_alfa*(mby+alfa*(uzi(ip)*mbx-uxi(ip)*mbz))+termo_correntey/vizty
    k4z(ip)=gama_alfa*(mbz+alfa*(uxi(ip)*mby-uyi(ip)*mbx))+termo_correntez/vizty
!
    uxi(ip)=uxi0(ip)+(k1x(ip)+2d0*k2x(ip)+2d0*k3x(ip)+k4x(ip))*dt/6d0
    uyi(ip)=uyi0(ip)+(k1y(ip)+2d0*k2y(ip)+2d0*k3y(ip)+k4y(ip))*dt/6d0
    uzi(ip)=uzi0(ip)+(k1z(ip)+2d0*k2z(ip)+2d0*k3z(ip)+k4z(ip))*dt/6d0
!
    end if
!
    end do
!
    return
    end
!
    subroutine runge_kuta_flag_cor1(i,ij)
    use VARIAVEIS
    integer :: ip,i,ij
    real(kind=8) :: mbx,mby,mbz,soma,p1vx,p1vy,p1vz,p2vx,p2vy,p2vz,p3vx,p3vy,p3vz
    real(kind=8) :: termo_correntex,termo_correntey,termo_correntez,vizty
!
    if(ij.eq.0) then
    uxi0=uxi
    uyi0=uyi
    uzi0=uzi
    end if
!
    do ip=1,npt
!
    vizty=magty(ty(ip))
    if(vizty.eq.0) vizty=1
!
    if(ij.eq.0) then

!!!!!!! Equação LLG dm/dt para ser integrada !!!!!!!!!!!!!!!!!!!!!

    mbx=uyi(ip)*campo_efetivo_z(ip)-uzi(ip)*campo_efetivo_y(ip)
    mby=uzi(ip)*campo_efetivo_x(ip)-uxi(ip)*campo_efetivo_z(ip)
    mbz=uxi(ip)*campo_efetivo_y(ip)-uyi(ip)*campo_efetivo_x(ip)
!
    ivxd1=gvxd1(ip)
    ivxe1=gvxe1(ip)
!
    select case (polar_sot)
!
    case(1)
!
!!! polarização na direção x !!!!!!!
!
      p1vx=0d0
      p1vy=-uzi(ip)
      p1vz=uyi(ip)

    case (2)

!! polarização na direção y !!!!!!!

     p1vx=-uzi(ip)
     p1vy=0d0
     p1vz=uxi(ip)
!
    end select
!
    p2vx=uyi(ip)*p1vz-uzi(ip)*p1vy
    p2vy=uzi(ip)*p1vx-uxi(ip)*p1vz
    p2vz=uxi(ip)*p1vy-uyi(ip)*p1vx
!
    termo_correntex=ct_sot*(p2vx-alfa*p1vx)
    termo_correntey=ct_sot*(p2vy-alfa*p1vy)
    termo_correntez=ct_sot*(p2vz-alfa*p1vz)
!
    k1x(ip)=gama_alfa*(mbx+alfa*(uyi(ip)*mbz-uzi(ip)*mby))+termo_correntex/vizty
    k1y(ip)=gama_alfa*(mby+alfa*(uzi(ip)*mbx-uxi(ip)*mbz))+termo_correntey/vizty
    k1z(ip)=gama_alfa*(mbz+alfa*(uxi(ip)*mby-uyi(ip)*mbx))+termo_correntez/vizty
!
    uxi(ip)=uxi0(ip)+k1x(ip)*dt/2d0
    uyi(ip)=uyi0(ip)+k1y(ip)*dt/2d0
    uzi(ip)=uzi0(ip)+k1z(ip)*dt/2d0
!
    end if
!
    if(ij.eq.1) then
!
    mbx=uyi(ip)*campo_efetivo_z(ip)-uzi(ip)*campo_efetivo_y(ip)
    mby=uzi(ip)*campo_efetivo_x(ip)-uxi(ip)*campo_efetivo_z(ip)
    mbz=uxi(ip)*campo_efetivo_y(ip)-uyi(ip)*campo_efetivo_x(ip)
!
    ivxd1=gvxd1(ip)
    ivxe1=gvxe1(ip)
!
    select case (polar_sot)
!
    case(1)
!
!!! polarização na direção x !!!!!!!
!
      p1vx=0d0
      p1vy=-uzi(ip)
      p1vz=uyi(ip)

    case (2)

!! polarização na direção y !!!!!!!

     p1vx=-uzi(ip)
     p1vy=0d0
     p1vz=uxi(ip)
!
    end select
!
    p2vx=uyi(ip)*p1vz-uzi(ip)*p1vy
    p2vy=uzi(ip)*p1vx-uxi(ip)*p1vz
    p2vz=uxi(ip)*p1vy-uyi(ip)*p1vx
!
    termo_correntex=ct_sot*(p2vx-alfa*p1vx)
    termo_correntey=ct_sot*(p2vy-alfa*p1vy)
    termo_correntez=ct_sot*(p2vz-alfa*p1vz)
!
    k2x(ip)=gama_alfa*(mbx+alfa*(uyi(ip)*mbz-uzi(ip)*mby))+termo_correntex/vizty
    k2y(ip)=gama_alfa*(mby+alfa*(uzi(ip)*mbx-uxi(ip)*mbz))+termo_correntey/vizty
    k2z(ip)=gama_alfa*(mbz+alfa*(uxi(ip)*mby-uyi(ip)*mbx))+termo_correntez/vizty
!
    uxi(ip)=uxi0(ip)+k2x(ip)*dt/2d0
    uyi(ip)=uyi0(ip)+k2y(ip)*dt/2d0
    uzi(ip)=uzi0(ip)+k2z(ip)*dt/2d0
!
    end if
!
    if(ij.eq.2) then
!
    mbx=uyi(ip)*campo_efetivo_z(ip)-uzi(ip)*campo_efetivo_y(ip)
    mby=uzi(ip)*campo_efetivo_x(ip)-uxi(ip)*campo_efetivo_z(ip)
    mbz=uxi(ip)*campo_efetivo_y(ip)-uyi(ip)*campo_efetivo_x(ip)
!
    ivxd1=gvxd1(ip)
    ivxe1=gvxe1(ip)
!
    select case (polar_sot)
!
    case(1)
!
!!! polarização na direção x !!!!!!!
!
      p1vx=0d0
      p1vy=-uzi(ip)
      p1vz=uyi(ip)

    case (2)

!! polarização na direção y !!!!!!!

     p1vx=-uzi(ip)
     p1vy=0d0
     p1vz=uxi(ip)
!
    end select
!
    p2vx=uyi(ip)*p1vz-uzi(ip)*p1vy
    p2vy=uzi(ip)*p1vx-uxi(ip)*p1vz
    p2vz=uxi(ip)*p1vy-uyi(ip)*p1vx
!
    termo_correntex=ct_sot*(p2vx-alfa*p1vx)
    termo_correntey=ct_sot*(p2vy-alfa*p1vy)
    termo_correntez=ct_sot*(p2vz-alfa*p1vz)
!
    k3x(ip)=gama_alfa*(mbx+alfa*(uyi(ip)*mbz-uzi(ip)*mby))+termo_correntex/vizty
    k3y(ip)=gama_alfa*(mby+alfa*(uzi(ip)*mbx-uxi(ip)*mbz))+termo_correntey/vizty
    k3z(ip)=gama_alfa*(mbz+alfa*(uxi(ip)*mby-uyi(ip)*mbx))+termo_correntez/vizty
!
    uxi(ip)=uxi0(ip)+k3x(ip)*dt
    uyi(ip)=uyi0(ip)+k3y(ip)*dt
    uzi(ip)=uzi0(ip)+k3z(ip)*dt
!
    end if
!
    if(ij.eq.3) then
!
    mbx=uyi(ip)*campo_efetivo_z(ip)-uzi(ip)*campo_efetivo_y(ip)
    mby=uzi(ip)*campo_efetivo_x(ip)-uxi(ip)*campo_efetivo_z(ip)
    mbz=uxi(ip)*campo_efetivo_y(ip)-uyi(ip)*campo_efetivo_x(ip)
!
    ivxd1=gvxd1(ip)
    ivxe1=gvxe1(ip)
!
    select case (polar_sot)
!
    case(1)
!
!!! polarização na direção x !!!!!!!
!
      p1vx=0d0
      p1vy=-uzi(ip)
      p1vz=uyi(ip)

    case (2)

!! polarização na direção y !!!!!!!

     p1vx=-uzi(ip)
     p1vy=0d0
     p1vz=uxi(ip)
!
    end select
!
    p2vx=uyi(ip)*p1vz-uzi(ip)*p1vy
    p2vy=uzi(ip)*p1vx-uxi(ip)*p1vz
    p2vz=uxi(ip)*p1vy-uyi(ip)*p1vx
!
    termo_correntex=ct_sot*(p2vx-alfa*p1vx)
    termo_correntey=ct_sot*(p2vy-alfa*p1vy)
    termo_correntez=ct_sot*(p2vz-alfa*p1vz)
!
    k4x(ip)=gama_alfa*(mbx+alfa*(uyi(ip)*mbz-uzi(ip)*mby))+termo_correntex/vizty
    k4y(ip)=gama_alfa*(mby+alfa*(uzi(ip)*mbx-uxi(ip)*mbz))+termo_correntey/vizty
    k4z(ip)=gama_alfa*(mbz+alfa*(uxi(ip)*mby-uyi(ip)*mbx))+termo_correntez/vizty
!
    uxi(ip)=uxi0(ip)+(k1x(ip)+2d0*k2x(ip)+2d0*k3x(ip)+k4x(ip))*dt/6d0
    uyi(ip)=uyi0(ip)+(k1y(ip)+2d0*k2y(ip)+2d0*k3y(ip)+k4y(ip))*dt/6d0
    uzi(ip)=uzi0(ip)+(k1z(ip)+2d0*k2z(ip)+2d0*k3z(ip)+k4z(ip))*dt/6d0
!
    end if
!
    end do
!
    return
    end
!
!
    subroutine adams_bashforth
    use VARIAVEIS
    integer :: ip,ivxd1,ivxe1
    real(kind=8) :: mbx,mby,mbz,soma,p1vx,p1vy,p1vz,p2vx,p2vy,p2vz,p3vx,p3vy,p3vz
!
    uxi0=uxi
    uyi0=uyi
    uzi0=uzi    
!
    if(adams.eq.0.or.adams.eq.1) then
!
    do ip=1,npt
!
!!!!!!! Equação LLg dm/dt para ser integrada !!!!!!!!!!!!!!!!!!!!!
!
    mbx=uyi(ip)*campo_efetivo_z(ip)-uzi(ip)*campo_efetivo_y(ip)
    mby=uzi(ip)*campo_efetivo_x(ip)-uxi(ip)*campo_efetivo_z(ip)
    mbz=uxi(ip)*campo_efetivo_y(ip)-uyi(ip)*campo_efetivo_x(ip)
!
!
    mdx(ip)=gama_alfa*(mbx+alfa*(uyi(ip)*mbz-uzi(ip)*mby))
    mdy(ip)=gama_alfa*(mby+alfa*(uzi(ip)*mbx-uxi(ip)*mbz))
    mdz(ip)=gama_alfa*(mbz+alfa*(uxi(ip)*mby-uyi(ip)*mbx))
!
    uxi(ip)=uxi(ip)+dt4*(55d0*mdx(ip)-59d0*md1x(ip)+37d0*md2x(ip)-9d0*md3x(ip))
    uyi(ip)=uyi(ip)+dt4*(55d0*mdy(ip)-59d0*md1y(ip)+37d0*md2y(ip)-9d0*md3y(ip))
    uzi(ip)=uzi(ip)+dt4*(55d0*mdz(ip)-59d0*md1z(ip)+37d0*md2z(ip)-9d0*md3z(ip))
!
    end do
!
    end if
!
    if(adams.eq.1) then
    
    if(flag_cut.eq.0) call energia_campo_local
    if(flag_cut.eq.1) call energia_campo_local_corte
!
    do ip=1,npt
!
!!!!!!! Equação LLg dm/dt para ser integrada !!!!!!!!!!!!!!!!!!!!!
!
    mbx=uyi(ip)*campo_efetivo_z(ip)-uzi(ip)*campo_efetivo_y(ip)
    mby=uzi(ip)*campo_efetivo_x(ip)-uxi(ip)*campo_efetivo_z(ip)
    mbz=uxi(ip)*campo_efetivo_y(ip)-uyi(ip)*campo_efetivo_x(ip)
!
    md0x(ip)=gama_alfa*(mbx+alfa*(uyi(ip)*mbz-uzi(ip)*mby))
    md0y(ip)=gama_alfa*(mby+alfa*(uzi(ip)*mbx-uxi(ip)*mbz))
    md0z(ip)=gama_alfa*(mbz+alfa*(uxi(ip)*mby-uyi(ip)*mbx))
!
    uxi(ip)=uxi0(ip)+dt4*(9d0*md0x(ip)+19d0*mdx(ip)-5d0*md1x(ip)+md2x(ip))
    uyi(ip)=uyi0(ip)+dt4*(9d0*md0y(ip)+19d0*mdy(ip)-5d0*md1y(ip)+md2y(ip))
    uzi(ip)=uzi0(ip)+dt4*(9d0*md0z(ip)+19d0*mdz(ip)-5d0*md1z(ip)+md2z(ip))
!
    end do
! 
    end if
! 
    md3x=md2x
    md2x=md1x
    md1x=mdx
    md3y=md2y
    md2y=md1y
    md1y=mdy
    md3z=md2z
    md2z=md1z
    md1z=mdz
!    
    return
    end
!    
!
!
    subroutine adams_bashforth_flag_cor0
    use VARIAVEIS
    integer :: ip,ij,ivxd1,ivxe1
    real(kind=8) :: mbx,mby,mbz,soma,p1vx,p1vy,p1vz,p2vx,p2vy,p2vz,p3vx,p3vy,p3vz
    real(kind=8) :: termo_correntex,termo_correntey,termo_correntez,vizty
!
    termo_correntex=0
    termo_correntey=0
    termo_correntez=0
    uxi0=uxi
    uyi0=uyi
    uzi0=uzi    
!
    if(adams.eq.0.or.adams.eq.1) then
!
    do ip=1,npt
!
    vizty=magty(ty(ip))
    if(vizty.eq.0) vizty=1
!!!!!!! Equação LLg dm/dt para ser integrada !!!!!!!!!!!!!!!!!!!!!
!
    mbx=uyi(ip)*campo_efetivo_z(ip)-uzi(ip)*campo_efetivo_y(ip)
    mby=uzi(ip)*campo_efetivo_x(ip)-uxi(ip)*campo_efetivo_z(ip)
    mbz=uxi(ip)*campo_efetivo_y(ip)-uyi(ip)*campo_efetivo_x(ip)
!
    ivxd1=gvxd1(ip)
    ivxe1=gvxe1(ip)
!
    p1vx=dist_x(ip)*(uyi(ip)*(uzi(ivxd1)-uzi(ivxe1))-uzi(ip)*(uyi(ivxd1)-uyi(ivxe1)))
    p1vy=dist_x(ip)*(uzi(ip)*(uxi(ivxd1)-uxi(ivxe1))-uxi(ip)*(uzi(ivxd1)-uzi(ivxe1)))
    p1vz=dist_x(ip)*(uxi(ip)*(uyi(ivxd1)-uyi(ivxe1))-uyi(ip)*(uxi(ivxd1)-uxi(ivxe1)))
!
    p2vx=uyi(ip)*p1vz-uzi(ip)*p1vy
    p2vy=uzi(ip)*p1vx-uxi(ip)*p1vz
    p2vz=uxi(ip)*p1vy-uyi(ip)*p1vx
!
    p3vx=uyi(ip)*p2vz-uzi(ip)*p2vy
    p3vy=uzi(ip)*p2vx-uxi(ip)*p2vz
    p3vz=uxi(ip)*p2vy-uyi(ip)*p2vx
!
    termo_correntex=gama_alfa1*(vj_qui1*p1vx+vj_qui_alfa2*p2vx+vj_alfa3*p3vx)
    termo_correntey=gama_alfa1*(vj_qui1*p1vy+vj_qui_alfa2*p2vy+vj_alfa3*p3vy)
    termo_correntez=gama_alfa1*(vj_qui1*p1vz+vj_qui_alfa2*p2vz+vj_alfa3*p3vz)
!
    mdx(ip)=gama_alfa*(mbx+alfa*(uyi(ip)*mbz-uzi(ip)*mby))+termo_correntex/vizty
    mdy(ip)=gama_alfa*(mby+alfa*(uzi(ip)*mbx-uxi(ip)*mbz))+termo_correntey/vizty
    mdz(ip)=gama_alfa*(mbz+alfa*(uxi(ip)*mby-uyi(ip)*mbx))+termo_correntez/vizty
    
!
    uxi(ip)=uxi(ip)+dt4*(55d0*mdx(ip)-59d0*md1x(ip)+37d0*md2x(ip)-9d0*md3x(ip))
    uyi(ip)=uyi(ip)+dt4*(55d0*mdy(ip)-59d0*md1y(ip)+37d0*md2y(ip)-9d0*md3y(ip))
    uzi(ip)=uzi(ip)+dt4*(55d0*mdz(ip)-59d0*md1z(ip)+37d0*md2z(ip)-9d0*md3z(ip))
!
    end do
!    
   end if
!
    if(adams.eq.1) then
!
    if(flag_cut.eq.0) call energia_campo_local
    if(flag_cut.eq.1) call energia_campo_local_corte
!
    do ip=1,npt
!
    vizty=magty(ty(ip))
    if(vizty.eq.0) vizty=1
!!!!!!! Equação LLg dm/dt para ser integrada !!!!!!!!!!!!!!!!!!!!!
!
    mbx=uyi(ip)*campo_efetivo_z(ip)-uzi(ip)*campo_efetivo_y(ip)
    mby=uzi(ip)*campo_efetivo_x(ip)-uxi(ip)*campo_efetivo_z(ip)
    mbz=uxi(ip)*campo_efetivo_y(ip)-uyi(ip)*campo_efetivo_x(ip)
!
    ivxd1=gvxd1(ip)
    ivxe1=gvxe1(ip)
!    
    p1vx=dist_x(ip)*(uyi(ip)*(uzi(ivxd1)-uzi(ivxe1))-uzi(ip)*(uyi(ivxd1)-uyi(ivxe1)))
    p1vy=dist_x(ip)*(uzi(ip)*(uxi(ivxd1)-uxi(ivxe1))-uxi(ip)*(uzi(ivxd1)-uzi(ivxe1)))
    p1vz=dist_x(ip)*(uxi(ip)*(uyi(ivxd1)-uyi(ivxe1))-uyi(ip)*(uxi(ivxd1)-uxi(ivxe1)))
!
    p2vx=uyi(ip)*p1vz-uzi(ip)*p1vy
    p2vy=uzi(ip)*p1vx-uxi(ip)*p1vz
    p2vz=uxi(ip)*p1vy-uyi(ip)*p1vx
!
    p3vx=uyi(ip)*p2vz-uzi(ip)*p2vy
    p3vy=uzi(ip)*p2vx-uxi(ip)*p2vz
    p3vz=uxi(ip)*p2vy-uyi(ip)*p2vx
!
    termo_correntex=gama_alfa1*(vj_qui1*p1vx+vj_qui_alfa2*p2vx+vj_alfa3*p3vx)
    termo_correntey=gama_alfa1*(vj_qui1*p1vy+vj_qui_alfa2*p2vy+vj_alfa3*p3vy)
    termo_correntez=gama_alfa1*(vj_qui1*p1vz+vj_qui_alfa2*p2vz+vj_alfa3*p3vz)
!
    md0x(ip)=gama_alfa*(mbx+alfa*(uyi(ip)*mbz-uzi(ip)*mby))+termo_correntex/vizty
    md0y(ip)=gama_alfa*(mby+alfa*(uzi(ip)*mbx-uxi(ip)*mbz))+termo_correntey/vizty
    md0z(ip)=gama_alfa*(mbz+alfa*(uxi(ip)*mby-uyi(ip)*mbx))+termo_correntez/vizty
!
    uxi(ip)=uxi0(ip)+dt4*(9d0*md0x(ip)+19d0*mdx(ip)-5d0*md1x(ip)+md2x(ip))
    uyi(ip)=uyi0(ip)+dt4*(9d0*md0y(ip)+19d0*mdy(ip)-5d0*md1y(ip)+md2y(ip))
    uzi(ip)=uzi0(ip)+dt4*(9d0*md0z(ip)+19d0*mdz(ip)-5d0*md1z(ip)+md2z(ip))
!
    end do
! 
    end if
! 
    md3x=md2x
    md2x=md1x
    md1x=mdx
    md3y=md2y
    md2y=md1y
    md1y=mdy
    md3z=md2z
    md2z=md1z
    md1z=mdz
!    
    return
    end
!
!
!
    subroutine adams_bashforth_flag_cor1
    use VARIAVEIS
    integer :: ip,ij,ivxd1,ivxe1,igx,igy
    real(kind=8) :: mbx,mby,mbz,soma,p1vx,p1vy,p1vz,p2vx,p2vy,p2vz,p3vx,p3vy,p3vz
    real(kind=8) :: termo_correntex,termo_correntey,termo_correntez,vizty
    real(kind=8) :: rxy,pcosx,pseny
!
    termo_correntex=0
    termo_correntey=0
    termo_correntez=0
    uxi0=uxi
    uyi0=uyi
    uzi0=uzi
!
    if(adams.eq.0.or.adams.eq.1) then
!
    do ip=1,npt
!
    vizty=magty(ty(ip))
    if(vizty.eq.0) vizty=1
!    
!!!!!!! Equação LLg dm/dt para ser integrada !!!!!!!!!!!!!!!!!!!!!
!
    mbx=uyi(ip)*campo_efetivo_z(ip)-uzi(ip)*campo_efetivo_y(ip)
    mby=uzi(ip)*campo_efetivo_x(ip)-uxi(ip)*campo_efetivo_z(ip)
    mbz=uxi(ip)*campo_efetivo_y(ip)-uyi(ip)*campo_efetivo_x(ip)
!
    select case (polar_sot)
!
    case(1)
!
!!! polarização na direção x !!!!!!!
!
      p1vx=0d0
      p1vy=-uzi(ip)
      p1vz=uyi(ip)

    case (2)

!! polarização na direção y !!!!!!!

     p1vx=-uzi(ip)
     p1vy=0d0
     p1vz=uxi(ip)
!
    end select
!
    p2vx=uyi(ip)*p1vz-uzi(ip)*p1vy
    p2vy=uzi(ip)*p1vx-uxi(ip)*p1vz
    p2vz=uxi(ip)*p1vy-uyi(ip)*p1vx
!
    termo_correntex=ct_sot*(p2vx-alfa*p1vx)
    termo_correntey=ct_sot*(p2vy-alfa*p1vy)
    termo_correntez=ct_sot*(p2vz-alfa*p1vz)
!
    mdx(ip)=gama_alfa*(mbx+alfa*(uyi(ip)*mbz-uzi(ip)*mby))+termo_correntex/vizty
    mdy(ip)=gama_alfa*(mby+alfa*(uzi(ip)*mbx-uxi(ip)*mbz))+termo_correntey/vizty
    mdz(ip)=gama_alfa*(mbz+alfa*(uxi(ip)*mby-uyi(ip)*mbx))+termo_correntez/vizty
!
    uxi(ip)=uxi(ip)+dt4*(55d0*mdx(ip)-59d0*md1x(ip)+37d0*md2x(ip)-9d0*md3x(ip))
    uyi(ip)=uyi(ip)+dt4*(55d0*mdy(ip)-59d0*md1y(ip)+37d0*md2y(ip)-9d0*md3y(ip))
    uzi(ip)=uzi(ip)+dt4*(55d0*mdz(ip)-59d0*md1z(ip)+37d0*md2z(ip)-9d0*md3z(ip))
!
    end do
!    
    end if
!
    if(adams.eq.1) then
!
    if(flag_cut.eq.0) call energia_campo_local
    if(flag_cut.eq.1) call energia_campo_local_corte
!
    do ip=1,npt
!
    vizty=magty(ty(ip))
    
!!!!!!! Equação LLg dm/dt para ser integrada !!!!!!!!!!!!!!!!!!!!!
!
    mbx=uyi(ip)*campo_efetivo_z(ip)-uzi(ip)*campo_efetivo_y(ip)
    mby=uzi(ip)*campo_efetivo_x(ip)-uxi(ip)*campo_efetivo_z(ip)
    mbz=uxi(ip)*campo_efetivo_y(ip)-uyi(ip)*campo_efetivo_x(ip)
!
    select case (polar_sot)

    case(1)
!
!!! polarização na direção x !!!!!!!
!
      p1vx=0d0
      p1vy=-uzi(ip)
      p1vz=uyi(ip)

    case (2)
!
!!! polarização na direção y !!!!!!!

     p1vx=-uzi(ip)
     p1vy=0d0
     p1vz=uxi(ip)
!
    end select
!
    p2vx=uyi(ip)*p1vz-uzi(ip)*p1vy
    p2vy=uzi(ip)*p1vx-uxi(ip)*p1vz
    p2vz=uxi(ip)*p1vy-uyi(ip)*p1vx
!
    termo_correntex=ct_sot*(p2vx-alfa*p1vx)
    termo_correntey=ct_sot*(p2vy-alfa*p1vy)
    termo_correntez=ct_sot*(p2vz-alfa*p1vz)

!
    md0x(ip)=gama_alfa*(mbx+alfa*(uyi(ip)*mbz-uzi(ip)*mby))+termo_correntex/vizty
    md0y(ip)=gama_alfa*(mby+alfa*(uzi(ip)*mbx-uxi(ip)*mbz))+termo_correntey/vizty
    md0z(ip)=gama_alfa*(mbz+alfa*(uxi(ip)*mby-uyi(ip)*mbx))+termo_correntez/vizty
!
    uxi(ip)=uxi0(ip)+dt4*(9d0*md0x(ip)+19d0*mdx(ip)-5d0*md1x(ip)+md2x(ip))
    uyi(ip)=uyi0(ip)+dt4*(9d0*md0y(ip)+19d0*mdy(ip)-5d0*md1y(ip)+md2y(ip))
    uzi(ip)=uzi0(ip)+dt4*(9d0*md0z(ip)+19d0*mdz(ip)-5d0*md1z(ip)+md2z(ip))
!
    end do
! 
    end if
! 
    md3x=md2x
    md2x=md1x
    md1x=mdx
    md3y=md2y
    md2y=md1y
    md1y=mdy
    md3z=md2z
    md2z=md1z
    md1z=mdz
!    
    return
    end































