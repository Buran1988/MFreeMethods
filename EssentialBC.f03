SUBROUTINE EssentialBC(numnode,pAlf,alfs,x,ds,ak,af,npEBCnum,npEBC,pEBC)
    !----------------------------------------------------------------------------
    ! This subroutine to enforce point essential bc's using penalty method;
    ! input--numnode: total number of field nodes;
    ! pAlf: penalty Fac; npEBCnum: number of e. b.c points
    ! alfs: coefficient of support domain
    ! x(nx,numnode): coordinates of all field nodes;
    ! input and output-- ak[]: stiffness matrix;
    ! af{}:force vector.
    !---------------------------------------------------------------------------
    implicit real*8 (a-h,o-z)
include 'parameters.h'
    COMMON/rpim/ALFC,DC,Q,nRBF
    common/basis/mbasis
    dimension npEBC(3,100),pEBC(2,100)
    dimension x(nx,numnode),ds(2,numnode),ak(2*numd,2*(numnode)),af(2*numnode)
    dimension nv(numnode),ph(10,numnode), x2(2)
    maxak=0.
    do iebc=1,2*numnode
        if(abs(ak(iebc,iebc)).gt.maxak) maxak=abs(ak(iebc,iebc))
    enddo
    do 10 iebc=1,npEBCnum
        ie=npEBC(1,iebc)
        x2(1)=x(1,ie)
        x2(2)=x(2,ie)
        ndex=0
        ! call support(x2,x,ds,nv(1),numnode,nx,ndex)
        call SupportDomain(numnode,nx,x2,x,ds,ndex,nv)
        do ik=1,ndex
            do jk=1,10
                ph(jk,ik)=0.
            enddo
        enddo
        call RPIM_ShapeFunc_2D(x2,x,nv,ph,nx,numnode,ndex,&
            alfc,dc,q,nRBF, mbasis)
        do iee=1,ndex
            ine=nv(iee)
            do ii=1,ndex
                jne=nv(ii)
                if(npEBC(2,iebc).eq.1) then
                    ak((ine*2-1),(jne*2-1))=ak((ine*2-1),(jne*2-1))-pAlf*maxak* &
                        ph(1,iee)*ph(1,ii)
                endif
                if(npEBC(3,iebc).eq.1) then
                    ak((ine*2),(jne*2))=ak((ine*2),(jne*2))-pAlf*maxak* &
                        ph(1,iee)*ph(1,ii)
                endif
            enddo
            if(npEBC(2,iebc).eq.1) then
                uu=pEBC(1,iebc)
                af(ine*2-1)=af(ine*2-1)-pAlf*uu*maxak*ph(1,iee)
            endif
            if(npEBC(3,iebc).eq.1) then
                uu=pEBC(2,iebc)
                af(ine*2)=af(ine*2)-pAlf*uu*maxak*ph(1,iee)
            endif
        enddo
10  continue
    RETURN
END
