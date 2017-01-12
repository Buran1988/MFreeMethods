SUBROUTINE CellGaussPoints(ibk,numcell,k,numq,numgauss,xc,noCell,gauss,gs)
    !----------------------------------------------------------------------------
    ! This subroutine to set up Gauss points,Jacobian and weights for a cell
    ! input--ibk: the No. of the consider cell;
    ! numq: number of points for background cells;
    ! numcell: number of background cells;
    ! numgauss: number of Gauss points in a cell;
    ! k: number of Gauss points used, numgauss=k*k for 2D cell;
    ! xc(nx,numq): coordinates of points for background cells;
    ! noCell(ng,numcell): No. of points to form this cell;
    ! gauss(2,k): coefficients of Gauss points;
    ! nx,ng: parameters are defined in file parameter.h.
    ! output--gs(ng,numgauss): coordinate of the Gauss points, weight and Jacobian
    !---------------------------------------------------------------------------
    implicit real*8 (a-h,o-z)
include 'parameters.h'
    dimension xc(nx,numq),noCell(ng,numcell),gauss(nx,k),gs(ng,numgauss)
    dimension psiJ(ng),etaJ(ng),xe(ng),ye(ng),aN(ng),aNJpsi(ng),aNJeta(ng)
    index=0
    psiJ(1)=-1.
    psiJ(2)=1.
    psiJ(3)=1.
    psiJ(4)=-1.
    etaJ(1)=-1.
    etaJ(2)=-1.
    etaJ(3)=1.
    etaJ(4)=1.
    l=k
    ie=ibk
    do j=1,ng
        je=noCell(j,ie)
        xe(j)=xc(1,je)
        ye(j)=xc(2,je)
    enddo
    do 10 i=1,l
        do 10 j=1,l
            index=index+1
            eta=gauss(1,i)
            psi=gauss(1,j)
            do ik=1,ng
                aN(ik)=.25*(1.+psi*psiJ(ik))*(1.+eta*etaJ(ik))
                aNJpsi(ik)=.25*psiJ(ik)*(1.+eta*etaJ(ik))
                aNJeta(ik)=.25*etaJ(ik)*(1.+psi*psiJ(ik))
            enddo
            xpsi=0.
            ypsi=0.
            xeta=0.
            yeta=0.
            do jk=1,ng
                xpsi=xpsi+aNJpsi(jk)*xe(jk)
                ypsi=ypsi+aNJpsi(jk)*ye(jk)
                xeta=xeta+aNJeta(jk)*xe(jk)
                yeta=yeta+aNJeta(jk)*ye(jk)
            enddo
            ajcob=xpsi*yeta-xeta*ypsi
            xq=0.
            yq=0.
            do kk=1,ng
                xq=xq+aN(kk)*xe(kk)
                yq=yq+aN(kk)*ye(kk)
            enddo
            gs(1,index)=xq
            gs(2,index)=yq
            gs(3,index)=gauss(2,i)*gauss(2,j)
            gs(4,index)=ajcob
10      continue
        RETURN
    END
