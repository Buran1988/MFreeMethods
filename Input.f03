SUBROUTINE Input(x,numd,nx,numnode,ndivx,ndivy,ndivxq,ndivyq,&
    nconn2,nquado,pAlf,Dmat,ALFs,numcell,numq,noCell,ncn,xc,&
    npEBCnum,npEBC,pEBC,npNBCnum,npNBC,pNBC)
    !------------------------------------------------------------------
    ! Input data from outside
    ! Output�all variables are output
    !------------------------------------------------------------------
    implicit real*8 (a-h,o-z)
    common/para/xlength,ylength,p,young,anu,aimo
    COMMON/rpim/ALFC,DC,Q,nRBF
    common /basis/mbasis
    CHARACTER*40 NAM
    dimension npEBC(3,100),pEBC(2,100),npNBC(3,100),pNBC(2,100)
    dimension x(nx,numd),Dmat(3,3),noCell(3,ncn),xc(nx,numd)

    read(4,10)nam
    read(4,*) xlength,ylength,young,anu,p
    read(4,10)nam
    read(4,*)numnode,nconn2
    read(4,10)nam
    read(4,*) ndivx,ndivy
    read(4,10)nam
    read(4,*)numq,numcell
    read(4,10)nam
    read(4,*)ndivxq,ndivyq
    read(4,10)nam
    read(4,*)nquado,pAlf
    read(4,10)nam
    read(4,*)ALFs
    numgauss=nquado*nquado
    read(4,10)nam
    do i=1,numnode
        read(4,*)j,x(1,i),x(2,i)
    enddo
    read(4,10)nam
    do i=1,numq
        read(4,*)j,xc(1,i),xc(2,i)
    enddo
    read(4,10)nam
    do j=1,numcell
        read(4,*)i,noCell(1,j),noCell(2,j),noCell(3,j)
    enddo
    read(4,10)nam
    read(4,*)npEBCnum
    read(4,10)nam
    do i=1,npEBCnum
        read(4,*)npEBC(1,i),npEBC(2,i),npEBC(3,i),pEBC(1,i),pEBC(2,i)
    enddo
    read(4,10)nam
    read(4,*)npNBCnum
    read(4,10)nam
    do i=1,npNBCnum
        read(4,*)npNBC(1,i),npNBC(2,i),npNBC(3,i),pNBC(1,i),pNBC(2,i)
    enddo
    read(4,10)nam
    READ(4,*)nRBF, alfc,dc, q
    read(4,10)nam
    READ(4,*)mbasis
    ! ************* Compute material matrix D[] for the plane stress
!    you=young/(1.-anu*anu)
!    aimo=(1./12.)*ylength**3
!    Dmat(1,1)=you
!    Dmat(1,2)=anu*you
!    Dmat(1,3)=0.
!    Dmat(2,1)=anu*you
!    Dmat(2,2)=you
!    Dmat(2,3)=0.
!    Dmat(3,1)=0.
!    Dmat(3,2)=0.
!    Dmat(3,3)=0.5*(1.-anu)*you
!------------------------------------------------
! Introducing composite material ROB
!([[2.29650146e+04, 6.08142164e+03, 5.68434189e-14],
!       [6.08142164e+03, 1.64983850e+04, 5.68434189e-14],
!       [5.68434189e-14, 5.68434189e-14, 6.49028044e+03]])
    Dmat(1,1)=2.29650146e+04
    Dmat(1,2)=6.08142164e+03
    Dmat(1,3)=0.
    Dmat(2,1)=6.08142164e+03
    Dmat(2,2)=1.64983850e+04
    Dmat(2,3)=0.
    Dmat(3,1)=0.
    Dmat(3,2)=0.
    Dmat(3,3)=6.49028044e+03

10  format(a40)
    RETURN
END
