SUBROUTINE EssentialBC(numnode,pAlf,x,ds,ak,af,npEBCnum,npEBC,pEBC)
    !----------------------------------------------------------------------------
    ! This subroutine to enforce point essential bc's using penalty method;
    ! input--numnode: total number of field nodes;
    ! pAlf: penalty Fac; npEBCnum: number of e. b.c points
    ! alfs: coefficient of support domain
    ! x(nx,numnode): coordinates of all field nodes;
    ! input and output-- ak[]: stiffness matrix;
    ! af{}:force vector.
    !---------------------------------------------------------------------------
    IMPLICIT REAL*8 (A-H,O-Z)
include 'parameters.h'
    COMMON/rpim/ALFC,DC,Q,nRBF
    COMMON/basis/mbasis
    DIMENSION npEBC(3,100),pEBC(2,100)
    DIMENSION x(nx,numnode),ds(2,numnode),ak(2*numd,2*(numnode)),af(2*numnode)
    DIMENSION nv(numnode),ph(10,numnode), x2(2)
    real:: maxak=0.
    DO iebc=1,2*numnode
        IF(abs(ak(iebc,iebc)).GT.maxak) maxak=abs(ak(iebc,iebc))
    ENDDO
    DO 10 iebc=1,npEBCnum
        ie=npEBC(1,iebc)
        x2(1)=x(1,ie)
        x2(2)=x(2,ie)
        ndex=0
        ! call support(x2,x,ds,nv(1),numnode,nx,ndex)
        CALL SupportDomain(numnode,nx,x2,x,ds,ndex,nv)
        DO ik=1,ndex
            DO jk=1,10
                ph(jk,ik)=0.
            ENDDO
        ENDDO
        CALL RPIM_ShapeFunc_2D(x2,x,nv,ph,nx,numnode,ndex,&
            alfc,dc,q,nRBF, mbasis)
        DO iee=1,ndex
            ine=nv(iee)
            DO ii=1,ndex
                jne=nv(ii)
                IF(npEBC(2,iebc).EQ.1) THEN
                    ak((ine*2-1),(jne*2-1))=ak((ine*2-1),(jne*2-1))-pAlf*maxak* &
                        ph(1,iee)*ph(1,ii)
                ENDIF
                IF(npEBC(3,iebc).EQ.1) THEN
                    ak((ine*2),(jne*2))=ak((ine*2),(jne*2))-pAlf*maxak* &
                        ph(1,iee)*ph(1,ii)
                ENDIF
            ENDDO
            IF(npEBC(2,iebc).EQ.1) THEN
                uu=pEBC(1,iebc)
                af(ine*2-1)=af(ine*2-1)-pAlf*uu*maxak*ph(1,iee)
            ENDIF
            IF(npEBC(3,iebc).EQ.1) THEN
                uu=pEBC(2,iebc)
                af(ine*2)=af(ine*2)-pAlf*uu*maxak*ph(1,iee)
            ENDIF
        ENDDO
10  CONTINUE
    RETURN
END
