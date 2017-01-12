SUBROUTINE NaturalBC_concentrated(x,nx,numnode,af,ds,alfs,npNBCnum,npNBC,pNBC)
    IMPLICIT REAL*8 (A-H,O-Z)
    DIMENSION npNBC(3,100),pNBC(2,100)
    COMMON/rpim/ALFC,DC,Q,nRBF
    COMMON/basis/mbasis
    DIMENSION af(2*numnode),x(nx,numnode),ds(nx,numnode)
    DIMENSION ph(10,numnode),gpos(2),nv(numnode)

    DO 10 iebc=1,npNBCnum
        ie=npNBC(1,iebc)
        gpos(1)=x(1,ie)
        gpos(2)=x(2,ie)
        ndex=0
        CALL SupportDomain(numnode,nx,gpos,x,ds,ndex,nv)

        ph=0.

        CALL RPIM_ShapeFunc_2D(gpos,x,nv,ph,nx,numnode,ndex,&
            alfc,dc,q,nRBF, mbasis)
        DO iee=1,ndex
            ie=nv(iee)
            uu=pNBC(1,iebc)
            af(ie*2-1)=af(ie*2-1)+ph(1,iee)*uu
            uu=pNBC(2,iebc)
            af(ie*2)=af(ie*2)+ph(1,iee)*uu
        ENDDO
10  CONTINUE
    RETURN
END
