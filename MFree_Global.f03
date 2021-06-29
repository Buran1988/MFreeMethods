PROGRAM MFree_Global

    !----------------------------------------------------------------------------
    ! main program--2D FORTRAN 90 CODE-MFree global weak-form methods
    ! Using square support domain and square background cells
    ! input file -- input.dat
    ! output file -- result.dat
    ! include file -- parameters.h, variables.h
    !----------------------------------------------------------------------------
    IMPLICIT REAL*8 (A-H,O-Z)
    CHARACTER input_file*140
include 'parameters.h'


include 'variables.h'

    ! ¬вожу переменную дл€ передачи имени файла из параметров запуска
    !    OPEN(ninput,FILE='Input175_55.dat')
    call get_command_argument(1, input_file)
    write(*,*) 'Input file aquired: ', input_file
    OPEN(ninput,FILE=input_file)

    OPEN(noutput,FILE='result.dat',STATUS='unknown')

    ! ************* Input data
    CALL input(x,numd,nx,numnode,ndivx,ndivy,ndivxq,ndivyq,&
        nconn2,nquado,pAlf,Dmat,ALFs,numcell,numq,noCell,ncn,xc,&
        npEBCnum,npEBC,pEBC,npNBCnum,npNBC,pNBC)

    numgauss=nquado*nquado !total number of Gauss points in a cell
    ! ************* Determine sizes of influence domains -- uniform nodal spacing
    xspace=xlength/ndivx
    yspace=ylength/ndivy
    DO i=1,numnode
        ds(1,i)=alfs*xspace
        ds(2,i)=alfs*yspace
    ENDDO
    ! ************* Coefficients of Gauss points,Weights and Jacobian for a cell
    CALL GaussCoefficient(nquado,gauss)


    gs=0.0
    force=0.0
    ak=0.0

    ! ************* Loop for background cells
    DO 10 ibk=1,numcell
        WRITE(*,*)'Cell No.=',ibk
        ! ************* Set Gauss points for this cell
        CALL CellGaussPoints(ibk,numcell,nquado,numq,numgauss, &
            xc,noCell,gauss,gs)
        ! ************* Loop over Gauss points to assemble discrete equations
        DO 20 ie=1,numgauss
            gpos(1)=gs(1,ie) ! Gauss point x
            gpos(2)=gs(2,ie) ! Gauss point y
            weight=gs(3,ie) ! weight coefficent
            ajac=gs(4,ie) ! Jacobian
            ! ************* Determine the support domain of Gauss point
            ndex=0
            CALL SupportDomain(numnode,nx,gpos,x,ds,ndex,nv)
            DO ik=1,3*ndex
                DO jk=1,10
                    ph(jk,ik)=0.
                ENDDO
            ENDDO
            ! ************* Construct RPIM shape functions for a Gauss point
            CALL RPIM_ShapeFunc_2D(gpos,x,nv,ph,nx,numnode,ndex,&
                alfc,dc,q,nRBF, mbasis)
            DO ik=1,2*ndex
                ne(ik)=0
            ENDDO
            DO ine=1,ndex
                n1=2*ine-1
                n2=2*ine
                ne(n1)=2*nv(ine)-1
                ne(n2)=2*nv(ine)
            ENDDO
            mbdb=4*ndex*ndex
            DO kbdb=1,mbdb
                GSPk(kbdb)=0.
            ENDDO
            ! ************* Compute the stiffness matrix for a Gauss point
            CALL PointStiffnessMatrix(ndex,weight,ajac,ph,Dmat,GSPk)
            nb=2*ndex
            DO ikk=1,nb
                DO jkk=1,nb
                    m1=ne(ikk)
                    m2=ne(jkk)
                    nbdb=(jkk-1)*nb+ikk
                    ak(m1,m2)=ak(m1,m2)+GSPk(nbdb)
                ENDDO
            ENDDO
20      CONTINUE ! end of loop for Gauss points
        ! ************* Implement natural BC
        in=0
        jn=0
        nn=noCell(3,ibk)
        IF(xc(1,nn)==xlength) in=nn

        nn= noCell(4,ibk)
        IF(xc(1,nn).EQ.xlength) jn=nn
        IF((in.NE.0).AND.(jn.NE.0)) THEN
            CALL naturalBC_distributed(numnode,numq,in,jn, &
                alfs,x,xc,ds,gauss,nquado,force)
        ENDIF
10  CONTINUE ! end of loop for cells
    ! ************* Boundary conditions: essential BC
    WRITE(*,*)' Boundary conditions....'
    nak=2*numd
    CALL EssentialBC(numnode,pAlf,alfs,x,ds,ak,force,npEBCnum,npEBC,pEBC)
    ! ************* Boundary conditions: concentrated natural BC
    CALL NaturalBC_concentrated(x,nx,numnode,force,ds,alfs,npNBCnum,npNBC,pNBC)
    nak=2*numd
    b=1.d-10
    ! ************* Solve equation to get the solutions
    WRITE(*,*)' Solving....'
    CALL SolverBand(ak,force,2*numnode,2*numd)
    nnn=2*numd
    DO ik=1,nx
        DO jk=1,numnode
            u2(ik,jk)=0.
        ENDDO
    ENDDO
    DO ik=1,numnode
        jk=2*ik-1
        u2(1,ik)=force(jk)
        u2(2,ik)=force(jk+1)
    ENDDO
    ! ************* Get the final displacement
    CALL GetDisplacement(x,ds,u2,disp,alfs,nx,numnode)
    ! ************* Get stress
    CALL GetStress(x,noCell,ds,Dmat,u2,alfs,nx,numnode,numgauss,&
        xc,gauss,nquado,ng,numq,numcell, ENORM,Stressnode)



    WRITE(*,*) 'Success!'

    STOP
END PROGRAM MFree_Global
