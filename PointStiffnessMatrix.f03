SUBROUTINE PointStiffnessMatrix(ndex,weight,ajac,ph,Dmat,GSPk)
    !----------------------------------------------------------------------------
    ! This subroutine to calculate sparse stiff matrix
    ! input--ndex: the number of nodes in the support domain;
    ! weight: weight of Gauss quadrature;
    ! ajac: Jacobian;
    ! dphix: first dirivetive of x of shape function;
    ! dphiy: first dirivetive of y of shape function;
    ! Dmat(3,3): the matrix of strain-stress;
    ! output--GSPk(2ndex,2ndex): sub-stiffness matrix of the Gauss point
    !---------------------------------------------------------------------------
    IMPLICIT REAL*8 (A-H,O-Z)
    DIMENSION ph(10,ndex),Dmat(3,3),GSPk(2*ndex,2*ndex)
    DIMENSION bmat(3,2*ndex),dphix(ndex),dphiy(ndex)
    nb=2*ndex
    DO i=1,ndex
        dphix(i)=ph(2,i)
        dphiy(i)=ph(3,i)
    ENDDO

    Bmat=0.

    DO in=1,ndex
        j=2*in-1
        k=2*in
        Bmat(1,j)=dphix(in)
        Bmat(1,k)=0.
        Bmat(2,j)=0.
        Bmat(2,k)=dphiy(in)
        Bmat(3,j)=dphiy(in)
        Bmat(3,k)=dphix(in)
    ENDDO

    GSPk=0.

    DO ii=1,nb
        DO jj=1,nb
            DO kk=1,3
                DO mm=1,3
                    GSPk(ii,jj)=GSPk(ii,jj)+weight*ajac*Bmat(kk,ii)* &
                        Dmat(kk,mm)*Bmat(mm,jj)
                ENDDO
            ENDDO
        ENDDO
    ENDDO
    RETURN
END
