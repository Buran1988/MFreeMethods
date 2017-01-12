SUBROUTINE SolverBand(ak,fp,neq,nmat)
    !------------------------------------------------------------------
    ! Sloving linear equations; it calls BandSolver & GaussSolver
    ! input—ak,fp,neq,nmat
    ! output--fp
    !------------------------------------------------------------------
    IMPLICIT REAL*8 (A-H,O-Z)
    DIMENSION ak(nmat,nEq),fp(nmat)
    REAL(8), ALLOCATABLE :: tp(:,:)
    REAL(8), ALLOCATABLE :: stfp(:,:)
    ALLOCATE (tp(1:neq,1:nmat))
    ALLOCATE (stfp(1:neq,1:neq))
    ep=1.e-10


            stfp=0.
            tp=0.

    DO i=1,nEq
        DO j=1,nEq
            stfp(i,j)=ak(i,j)
        ENDDO
    ENDDO
    ni=nEq
    Lp=0 ! half band width
    DO 20 i=1,ni
        DO j=ni,i,-1
            IF(stfp(i,j).NE.0.) THEN ! stfp[,] stifness matrix
                IF(abs(j-i).GT.Lp) Lp=abs(j-i)
                GO TO 21
            ENDIF
        ENDDO
21  CONTINUE
    DO j=1,i
        IF(stfp(i,j).NE.0.) THEN
            IF(abs(j-i).GT.Lp) Lp=abs(j-i)
            GO TO 20
        ENDIF
    ENDDO
20 CONTINUE
   ilp=2*lp+1 ! band width
   nm=nEq
   IF(ilp.LT.nEq) THEN
       CALL BandSolver(stfp,fp,tp,nm,lp,ilp,nmat) ! solver for band matrix
   ELSE
       CALL GaussSolver(nEq,nmat,ak,fp,ep,kkkk) ! standard solver
   ENDIF
   DEALLOCATE (tp)
   DEALLOCATE (stfp)
   END
