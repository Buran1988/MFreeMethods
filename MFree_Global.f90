program MFree_Global

    !----------------------------------------------------------------------------
    ! main program--2D FORTRAN 90 CODE-MFree global weak-form methods
    ! Using square support domain and square background cells
    ! input file -- input.dat
    ! output file -- result.dat
    ! include file -- parameter.h, variable.h
    !----------------------------------------------------------------------------
    implicit real*8 (a-h,o-z)
include 'parameters.h'
include 'variables.h'

    open(ninput,file='Input175_55.dat')
    open(noutput,file='result.dat',status='unknown')
    ! ************* Input data
!    call input(x,numd,nx,numnode,ndivx,ndivy,ndivxq,ndivyq,&
!        nconn2,nquado,pAlf,Dmat,ALFs,numcell,numq,noCell,ncn,xc,&
!        npEBCnum,npEBC,pEBC,npNBCnum,npNBC,pNBC)
!
    numgauss=nquado*nquado !total number of Gauss points in a cell
    ! ************* Determine sizes of influence domains -- uniform nodal spacing
    xspace=xlength/ndivx
    yspace=ylength/ndivy
    do i=1,numnode
        ds(1,i)=alfs*xspace
        ds(2,i)=alfs*yspace
    enddo

    ! ************* Coefficients of Gauss points,Weights and Jacobian for a cell
    !call GaussCoefficient(nquado,gauss)

    do ik=1,ng
        do jk=1,numgauss
            gs(ik,jk)=0
        enddo
    enddo

    do ik=1,2*numd
        force(ik)=0.
        do jk=1,2*numd
            ak(ik,jk)=0.
        enddo
    enddo

end program MFree_Global
