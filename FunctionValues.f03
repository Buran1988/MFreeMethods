
module FunctionValues
    IMPLICIT NONE
contains

    function GetGaussPointsData(bg_cell_id, v, xc, noCell ) result(gs)

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


        INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
        integer, intent(in)  :: bg_cell_id, noCell(:,:) ! input

        REAL(iwp), intent(in) :: xc(:,:)
        REAL(iwp):: v(:,:)

        INTEGER::ngauss_points, ig, grids_in_cell, j,  je


        REAL(iwp), dimension(4,3 ) :: gs !gs: x_point, y_point, W, Jacobian
        REAL(iwp) :: xe(ubound(noCell, 1)), ye(ubound(noCell, 1))


        real(iwp):: L1, L2, L3, weight, xq=0.0, yq=0.0, detJ = 0.0


        ngauss_points = ubound(v, 2)
        grids_in_cell = ubound(noCell, 1)


        do j=1,grids_in_cell
            je=noCell(j,bg_cell_id)
            xe(j)=xc(1,je)
            ye(j)=xc(2,je)
        enddo


        gauss_points: DO ig=1,ngauss_points
            L1 =  1 - v(1,ig) - v(2,ig)
            L2 = v(1,ig)
            L3 = v(2,ig)
            xq = L1*xe(1) + L2*xe(2) + L3*xe(3)
            yq = L1*ye(1) + L2*ye(2) + L3*ye(3)

            weight = v(3,ig)

            !            detJ=(xB-xA)(yC-yA)-(xC-xA)(yB-yA)

            detJ = (xe(2) - xe(1))*(ye(3) - ye(1)) - (xe(3) - xe(1) )*(ye(2) - ye(1) )

            gs(1,ig) = xq
            gs(2,ig) = yq
            gs(3,ig) = weight
            gs(4,ig) = detJ
        END DO gauss_points


    end function GetGaussPointsData

end module
