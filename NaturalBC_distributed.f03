SUBROUTINE naturalBC_distributed(numnode,numq,in,jn,alfs,x,xc,ds, &
    gauss,nquado,force)
    !----------------------------------------------------------------------------
    ! This subroutine to enforce point natural bc's;
    ! input—numnode, numq, in,jn,alfs,x,xc,ds,gauss, nquado.
    ! input and output-- force{}:force vector.
    !---------------------------------------------------------------------------
    IMPLICIT REAL*8 (A-H,O-Z)
include 'parameters.h'
    COMMON/para/xlength,ylength,p,young,anu,aimo
    COMMON/rpim/ALFC,DC,Q,nRBF
    COMMON /basis/mbasis
    DIMENSION xc(nx,numq),gauss(2,nquado),force(2*numnode),x(nx,numnode)
    DIMENSION ph(10,numnode),gpos(2),nv(numnode),ds(nx,numnode)
    ax=0.5*(xc(1,in)-xc(1,jn))
    ay=0.5*(xc(2,in)-xc(2,jn))
    bx=0.5*(xc(1,in)+xc(1,jn))
    by=0.5*(xc(2,in)+xc(2,jn))

    DO il=1,nquado
        gpos(1)=ax*gauss(1,il)+bx
        gpos(2)=ay*gauss(1,il)+by
        weight=gauss(2,il)
        ajac=0.5*sqrt((xc(1,in)-xc(1,jn))**2+(xc(2,in)-xc(2,jn))**2)
        aimo=(1./12.)*ylength**3
        ty=(-1000./(2.*aimo))*(ylength*ylength/4.-gpos(2)*gpos(2))
        CALL SupportDomain(numnode,nx,gpos,x,ds,ndex,nv)

        ph=0.

        CALL RPIM_ShapeFunc_2D(gpos,x,nv,ph,nx,numnode,ndex,&
            alfc,dc,q,nRBF, mbasis)
        DO ie=1,ndex
            nn=nv(ie)
            force(2*nn)=force(2*nn)+weight*ajac*ph(1,ie)*ty
        ENDDO
    ENDDO
END
