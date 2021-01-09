program dvr_pcet

!dvr calculation to calculate the mixed electron-proton adiabats as a function of the solvent coordinate for the new2 golden rule comparison system-bath model (linear solvent coupling, coulombic e-p coupling)
!output will be the lowest four mixed electron-proton eigenstates as long s as well as electron coulombic parameters

implicit none

real*8, parameter :: pi=dacos(-1.d0),hbar=1.d0

integer i,j,k,p,nqe,ns,nx,ii,jj,int1,int2

real*8 :: rd,ra,rcut,rdo,rdi,rao,rai,qd,qa,rcut_ep,dis
real*8 :: ad,bd,cd,aa,ba,ca
real*8 :: me,mx,wx,vo,yx,ms,ws,mues,mups,muep
real*8 :: qemin,qemax,smin,smax,xmin,xmax
real*8 :: dqe,ds,dx
real*8 :: vp,vs,ves,vps,vep,conste,constp
real*8, allocatable :: qe(:),s(:),x(:)
real*8, allocatable :: coulpot(:),dvrh(:,:)
character*10 :: str

logical :: safe=.false.

!variables used in dsyev
integer :: ifail=0, lwork
real*8, allocatable :: egval(:), work(:)

external dsyev

write(*,*) 'reading parameters'

!read in parameters
open(unit=1,file='system')
read(1,*) str,me
read(1,*) str,rcut
read(1,*) str,mx
read(1,*) str,wx
read(1,*) str,vo
read(1,*) str,yx
read(1,*) str,ms
read(1,*) str,ws
read(1,*) str,mues
read(1,*) str,mups
read(1,*) str,muep
read(1,*) str,rcut_ep
close(1)

!read in control
open(unit=2,file='control')
read(2,*)
read(2,*)
read(2,*) str,qemin,nqe,qemax
read(2,*) str,smin,ns,smax
read(2,*) str,xmin,nx,xmax
read(2,*) str,qd
read(2,*) str,ra
close(2)

qa=qd
rd=-ra

!grid spacings
dqe=(qemax-qemin)/dble(nqe)
ds=(smax-smin)/dble(ns)
dx=(xmax-xmin)/dble(nx)

!create arrays of the position of ele, solv, and htr for the grids
allocate(qe(nqe))
allocate(x(nx))
allocate(s(ns))
do i=1,nqe
  qe(i)=qemin+(i-1)*dqe
enddo
do i=1,nx
  x(i)=xmin+(i-1)*dx
enddo
do i=1,ns
  s(i)=smin+(i-1)*ds
enddo


!allocate the components for dysev, the diagonalization
allocate(egval(nqe*nx))
lwork=3*nqe*nx-1
allocate(work(lwork))

!allocate all other arrays
allocate(coulpot(nqe))
allocate(dvrh(nqe*nx,nqe*nx))

write(*,*) 'qd= ',qd,'rd= ',rd

write(*,*) 'calculating coulombic potential'
!calculate capped parabolas for coulombic potential
!note that rd should be less than ra
rdo=rd-rcut
rdi=rd+rcut
rai=ra-rcut
rao=ra+rcut

aa=(qa*(rai-ra)/dabs(ra-rai)**3+qd*(rai-rd)/dabs(rai-rd)**3-((rao**2-rai**2)*(qa*(rai-ra)/dabs(ra-rai)**3+qd*(rai-rd)/dabs(rai-rd)**3)+2.d0*rai*(-qa/dabs(ra-rai)+qa/dabs(ra-rao)+qd*(-1.d0/dabs(rai-rd)+1.d0/dabs(rao-rd))))/(rai-rao)**2)/(2.d0*rai)

ba=((rao**2-rai**2)*(qa*(rai-ra)/dabs(ra-rai)**3+qd*(rai-rd)/dabs(rai-rd)**3)+2.d0*rai*(-qa/dabs(ra-rai)+qa/dabs(ra-rao)+qd*(-1.d0/dabs(rai-rd)+1.d0/dabs(rao-rd))))/(rai-rao)**2

ca=(-qa*(ra-rai)*rai*(rai-rao)*rao*dabs(ra-rao)*abs(rai-rd)**3*dabs(rao-rd)+qa*(2.d0*rai-rao)*rao*dabs(ra-rai)**2*dabs(ra-rao)*dabs(rai-rd)**3*dabs(rao-rd)+dabs(ra-rai)**3*(-qa*rai**2*dabs(rai-rd)**3*dabs(rao-rd)+qd*dabs(ra-rao)*(-rai**2*dabs(rai-rd)**3+rai*(rai-rao)*rao*(rai-rd)*dabs(rao-rd)+(2.d0*rai-rao)*rao*dabs(rai-rd)**2*dabs(rao-rd))))/((rai-rao)**2*dabs(ra-rai)**3*dabs(ra-rao)*dabs(rai-rd)**3*dabs(rao-rd))

ad=((qa*(rdi-ra))/dabs(ra-rdi)**3-((-(rdi**2-rdo**2))*((qa*(rdi-ra))/dabs(ra-rdi)**3+(qd*(rdi-rd))/dabs(rd-rdi)**3)-2*rdi*(qa/dabs(ra-rdi)-qa/dabs(ra-rdo)+qd/dabs(rd-rdi)-qd/dabs(rd-rdo)))/(rdi-rdo)**2+(qd*(rdi-rd))/dabs(rd-rdi)**3)/(2*rdi)

bd=((rdi**2-rdo**2)*(qa*(rdi-ra)/dabs(ra-rdi)**3+qd*(rdi-rd)/dabs(rd-rdi)**3)+2.d0*rdi*(qa/dabs(ra-rdi)+qd/dabs(rd-rdi)-qa/dabs(ra-rdo)-qd/dabs(rd-rdo)))/(-(rdi-rdo)**2)

cd=(-qa*(ra-rdi)*rdi*(rdi-rdo)*rdo*dabs(rd-rdi)**3*dabs(ra-rdo)*dabs(rd-rdo)+qa*(2*rdi-rdo)*rdo*dabs(ra-rdi)**2*dabs(rd-rdi)**3*dabs(ra-rdo)*dabs(rd-rdo)+dabs(ra-rdi)**3*(-qd*(rd-rdi)*rdi*(rdi-rdo)*rdo*dabs(ra-rdo)*dabs(rd-rdo)+qd*(2*rdi-rdo)*rdo*dabs(rd-rdi)**2*dabs(ra-rdo)*dabs(rd-rdo)-rdi**2*dabs(rd-rdi)**3*(qd*dabs(ra-rdo)+qa*dabs(rd-rdo))))/((rdi-rdo)**2*dabs(ra-rdi)**3*dabs(rd-rdi)**3*dabs(ra-rdo)*dabs(rd-rdo))


!calculate the coulombic potential for each electron position
do i=1,nqe

  if((qe(i).gt.rdo).and.(qe(i).lt.rdi)) then
    coulpot(i)=ad*qe(i)**2+bd*qe(i)+cd
  elseif((qe(i).gt.rai).and.(qe(i).lt.rao)) then
    coulpot(i)=aa*qe(i)**2+ba*qe(i)+ca
  else
    coulpot(i)=-qd/dabs(qe(i)-rd)-qa/dabs(qe(i)-ra)
  endif

enddo


!initialize the dvr matrix
dvrh=0.d0

open(unit=110,file='adiabats.dat')

!run through the solvent positions
do i=1,ns

    write(*,*) 'making the dvr matrix for s= ',s(i)

!   construct the dvr matrix
    do k=1,nqe
      do ii=1,nx
        do p=1,nqe
          do jj=1,nx

            conste=(hbar**2)*(-1)**(k-p)/(2.d0*me*dqe**2)
            constp=(hbar**2)*(-1)**(ii-jj)/(2.d0*mx*dx**2)

!           index of the dvr matrix for given electron/proton positions
            int1=nx*(k-1)+ii
            int2=nx*(p-1)+jj

            if(k.eq.p.and.ii.eq.jj) then

!             proton potential
              vp=-mx*wx**2/2.d0*x(ii)**2+mx**2*wx**4/(16.d0*vo)*x(ii)**4-yx*x(ii)**3

!             solvent potential
              vs=0.5d0*ms*ws**2*s(i)**2

!             electron-solvent potential
              ves=-mues*s(i)*qe(k)

!             proton-solvent potential
              vps=-mups*s(i)*x(ii)

!             electron-proton potential
              dis=abs(qe(k)-x(ii))

              if(dis.gt.rcut_ep)then
                vep=-muep/dis
              else
                vep=-muep/rcut_ep
              endif

!             diagonal element of the dvr matrix
              dvrh(int1,int2)=coulpot(k)+vp+vs+ves+vps+vep+conste*(pi**2)/3.d0+constp*(pi**2)/3.d0
          
            elseif(k.ne.p.and.ii.eq.jj) then

!             kinetic contribution from the electron, so off-diagonal in electron, but still diagonal in proton
              dvrh(int1,int2)=conste*2.d0/(k-p)**2

            elseif(k.eq.p.and.ii.ne.jj) then

!             kinetic contribution from the proton, so off-diagonal in proton, but still diagonal in electron
              dvrh(int1,int2)=constp*2.d0/(ii-jj)**2

!           other term, where off diagonal in both electron and proton, is equal to 0

            endif

          enddo
        enddo
      enddo
    enddo

    write(*,*) 'diagonalizing the dvr matrix for s= ',s(i)

!   diagonalize the dvr matrix to get the adiabats using the nag subroutine dsyev
    call dsyev ('v' , 'l' , nqe*nx , dvrh , nqe*nx , egval , work , lwork , ifail )

    write(*,*) 'diagonalization complete, ifail=',ifail

!   write out energy of lowest 4 adiabats
    write(110,"(12e15.6)") s(i),egval(1),egval(2),egval(3),egval(4)

enddo

!write parameters used for the coulombic potential
open(unit=100,file='pos_rep_param.dat')
  write(100,*) 'qd= ',qd
  write(100,*) 'qa= ',qa
  write(100,*) 'rd= ',rd
  write(100,*) 'ra= ',ra
  write(100,*) 'ad= ',ad
  write(100,*) 'bd= ',bd
  write(100,*) 'cd= ',cd
  write(100,*) 'aa= ',aa
  write(100,*) 'ba= ',ba
  write(100,*) 'ca= ',ca


end program dvr_pcet
