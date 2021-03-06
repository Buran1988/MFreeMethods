SUBROUTINE GaussCoefficient(k,v)
    !----------------------------------------------------------------------------
    ! This subroutine returns a matrix with Gauss points and their weights
    ! input--k: k -- number of Gauss points;
    ! output--v(2,k): weight matrix of k Gauss points
    !---------------------------------------------------------------------------
    implicit real*8 (a-h,o-z)
    dimension v(2,k)
    SELECT CASE (k)
        Case (2)
            v(1,1)=-.57735
            v(1,2)=-v(1,1)
            v(2,1)=1.00000
            v(2,2)=v(2,1)
        Case (3)
            v(1,1)=-.77459
            v(1,2)=-.00000
            v(1,3)=-v(1,1)
            v(2,1)=.55555
            v(2,2)=.88888
            v(2,3)=v(2,1)
        Case (4)
            v(1,1)=-.86113
            v(1,2)=-.33998
            v(1,3)=-v(1,2)
            v(1,4)=-v(1,1)
            v(2,1)=.34785
            v(2,2)=.65214
            v(2,3)=v(2,2)
            v(2,4)=v(2,1)
        Case (6)
            v(1,1)=-.93246
            v(1,2)=-.66120
            v(1,3)=-.23861
            v(1,4)=-v(1,3)
            v(1,5)=-v(1,2)
            v(1,6)=-v(1,1)
            v(2,1)=.17132
            v(2,2)=.36076
            v(2,3)=.46791
            v(2,4)=v(2,3)
            v(2,5)=v(2,2)
            v(2,6)=v(2,1)
        Case (8)
            v(1,1)=-.96028
            v(1,2)=-.79666
            v(1,3)=-.52553
            v(1,4)=-.18343
            v(1,5)=-v(1,4)
            v(1,6)=-v(1,3)
            v(1,7)=-v(1,2)
            v(1,8)=-v(1,1)
            v(2,1)=.10122
            v(2,2)=.22238
            v(2,3)=.31370
            v(2,4)=.36268
            v(2,5)=v(2,4)
            v(2,6)=v(2,3)
            v(2,7)=v(2,2)
            v(2,8)=v(2,1)
    end select
    RETURN
END
