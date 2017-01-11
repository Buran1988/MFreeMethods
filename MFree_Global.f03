program MFree_Global

    !----------------------------------------------------------------------------
    ! main program--2D FORTRAN 90 CODE-MFree global weak-form methods
    ! Using square support domain and square background cells
    ! input file -- input.dat
    ! output file -- result.dat
    ! include file -- parameters.h, variables.h
    !----------------------------------------------------------------------------
    implicit real*8 (a-h,o-z)
include 'parameters.h'
include 'variables.h'

    open(ninput,file='Input175_55.dat')
    open(noutput,file='result.dat',status='unknown')

    ! ************* Input data
    call input(x,numd,nx,numnode,ndivx,ndivy,ndivxq,ndivyq,&
        nconn2,nquado,pAlf,Dmat,ALFs,numcell,numq,noCell,ncn,xc,&
        npEBCnum,npEBC,pEBC,npNBCnum,npNBC,pNBC)

    numgauss=nquado*nquado !total number of Gauss points in a cell
    ! ************* Determine sizes of influence domains -- uniform nodal spacing
    xspace=xlength/ndivx
    yspace=ylength/ndivy
    do i=1,numnode
        ds(1,i)=alfs*xspace
        ds(2,i)=alfs*yspace
    enddo
    ! ************* Coefficients of Gauss points,Weights and Jacobian for a cell
    call GaussCoefficient(nquado,gauss)

    do ik=1,ng
        do jk=1,numgauss
            gs(ik,jk)=0
        enddo
    enddo


    force=0.0
    ak=0.0

    ! ************* Loop for background cells
    do 10 ibk=1,numcell
        write(*,*)'Cell No.=',ibk
        ! ************* Set Gauss points for this cell
        call CellGaussPoints(ibk,numcell,nquado,numq,numgauss, &
            xc,noCell,gauss,gs)
        ! ************* Loop over Gauss points to assemble discrete equations
        do 20 ie=1,numgauss
            gpos(1)=gs(1,ie) ! Gauss point x
            gpos(2)=gs(2,ie) ! Gauss point y
            weight=gs(3,ie) ! weight coefficent
            ajac=gs(4,ie) ! Jacobian
            ! ************* Determine the support domain of Gauss point
            ndex=0
            call SupportDomain(numnode,nx,gpos,x,ds,ndex,nv)
            do ik=1,3*ndex
                do jk=1,10
                    ph(jk,ik)=0.
                enddo
            enddo
            ! ************* Construct RPIM shape functions for a Gauss point
            call RPIM_ShapeFunc_2D(gpos,x,nv,ph,nx,numnode,ndex,&
                alfc,dc,q,nRBF, mbasis)
            do ik=1,2*ndex
                ne(ik)=0
            enddo
            do ine=1,ndex
                n1=2*ine-1
                n2=2*ine
                ne(n1)=2*nv(ine)-1
                ne(n2)=2*nv(ine)
            enddo
            mbdb=4*ndex*ndex
            do kbdb=1,mbdb
                GSPk(kbdb)=0.
            enddo
            ! ************* Compute the stiffness matrix for a Gauss point
            call PointStiffnessMatrix(ndex,weight,ajac,ph,Dmat,GSPk)
            nb=2*ndex
            do ikk=1,nb
                do jkk=1,nb
                    m1=ne(ikk)
                    m2=ne(jkk)
                    nbdb=(jkk-1)*nb+ikk
                    ak(m1,m2)=ak(m1,m2)+GSPk(nbdb)
                enddo
            enddo
20      continue ! end of loop for Gauss points
        ! ************* Implement natural BC
        in=0
        jn=0
        nn=noCell(3,ibk)
        if(xc(1,nn)==xlength) in=nn

        nn= noCell(4,ibk)
        if(xc(1,nn).eq.xlength) jn=nn
        if((in.ne.0).and.(jn.ne.0)) then
            call naturalBC_distributed(numnode,numq,in,jn, &
                alfs,x,xc,ds,gauss,nquado,force)
        endif
10  continue ! end of loop for cells
    ! ************* Boundary conditions: essential BC
    write(*,*)' Boundary conditions....'
    nak=2*numd
    call EssentialBC(numnode,pAlf,alfs,x,ds,ak,force,npEBCnum,npEBC,pEBC)
    ! ************* Boundary conditions: concentrated natural BC
    call NaturalBC_concentrated(x,nx,numnode,force,ds,alfs,npNBCnum,npNBC,pNBC)
    nak=2*numd
    b=1.d-10
    ! ************* Solve equation to get the solutions
    write(*,*)' Solving....'
    call SolverBand(ak,force,2*numnode,2*numd)
    nnn=2*numd
    do ik=1,nx
        do jk=1,numnode
            u2(ik,jk)=0.
        enddo
    enddo
    do ik=1,numnode
        jk=2*ik-1
        u2(1,ik)=force(jk)
        u2(2,ik)=force(jk+1)
    enddo
    ! ************* Get the final displacement
    call GetDisplacement(x,ds,u2,disp,alfs,nx,numnode)
    ! ************* Get stress
    call GetStress(x,noCell,ds,Dmat,u2,alfs,nx,numnode,numgauss,&
        xc,gauss,nquado,ng,numq,numcell, ENORM,Stressnode)
    STOP


    write(*,*) 'Success!'


end program MFree_Global
