SUBROUTINE GaussSolver(n,mk,a,b,ep,kwji)
    !------------------------------------------------------------------
    ! Stnadard Gauss elimination slover for linear equations that are
    ! not suitably solved by BandSolver.
    !------------------------------------------------------------------
    IMPLICIT REAL*8 (A-H,O-Z)
    DIMENSION a(mk,mk),b(mk)
    INTEGER, ALLOCATABLE :: m(:)
    ALLOCATE (m(2*mk))
    ep=1.0E-10
    kwji=0
    DO i=1,n
        m(i)=i
    ENDDO
    DO 20 k=1,n
        p=0.0
        DO 30 i=k,n
            DO 30  j=k,n
                IF(abs(a(i,j)).LE.abs(p)) GOTO 30
                p=a(i,j)
                io=i
                jo=j
30      continue

            IF(abs(p)-ep) 200,200,300
200         kwji=1
            RETURN
300         IF(jo.EQ.k) GOTO 400
            DO i=1,n
                t=a(i,jo)
                a(i,jo)=a(i,k)
                a(i,k)=t
            ENDDO
            j=m(k)
            m(k)=m(jo)
            m(jo)=j
400         IF(io.EQ.k) GOTO 500
            DO j=k,n
                t=a(io,j)
                a(io,j)=a(k,j)
                a(k,j)=t
            ENDDO
            t=b(io)
            b(io)=b(k)
            b(k)=t
500         p=1.0/p
            in=n-1
            IF(k.EQ.n) GOTO 600
            DO j=k,in
                a(k,j+1)=a(k,j+1)*p
            ENDDO
600         b(k)=b(k)*p
            IF(k.EQ.n) GOTO 20
            DO i=k,in
                DO j=k,in
                    a(i+1,j+1)=a(i+1,j+1)-a(i+1,k)*a(k,j+1)
                ENDDO
                b(i+1)=b(I+1)-a(i+1,k)*b(k)
            ENDDO
20      CONTINUE
        DO i1=2,n
            i=n+1-i1
            DO j=i,in
                b(i)=b(i)-a(i,j+1)*b(j+1)
            ENDDO
        ENDDO
        DO k=1,n
            i=m(k)
            a(1,i)=b(k)
        ENDDO
        DO k=1,n
            b(k)=a(1,k)
        ENDDO
        kwji=0
        DEALLOCATE (m)
        RETURN
    END
