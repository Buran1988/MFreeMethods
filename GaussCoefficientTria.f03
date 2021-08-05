SUBROUTINE GaussTriaCoefficient(k,v)
       !----------------------------------------------------------------------------
       ! This subroutine returns a matrix with Gauss points and their weights
       ! input--k: k -- number of Gauss points;
       ! output--v(3,k): weight matrix of k Gauss points (mu, lambda, wk)
       !---------------------------------------------------------------------------

    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp)::zero=0.0_iwp,half=0.5_iwp,one_to_six = 1.0/6_iwp
    integer,INTENT(in)::k
    REAL(iwp),INTENT(out)::v(3,k)

    select case (k)
        case (3)
            v(1,1)=     half
            v(2,1)=     half


            v(1,2)=    half
            v(2,2)=     zero

            v(1,3)=    zero
            v(2,3)=    half


            v(3,:)=     one_to_six

        case (9)

            !     mu
            v(1,1)=     0.10635080
            v(1,2)=     0.47182460
            v(1,3)=     0.83729830
            v(1,4)=     0.08452624
            v(1,5)=     0.37500000
            v(1,6)=     0.66547380
            v(1,7)=     0.06270166
            v(1,8)=     0.27817540
            v(1,9)=     0.49364920


            !     lambda
            v(2,1)=     0.106350800
            v(2,2)=     0.084526240
            v(2,3)=     0.062701660
            v(2,4)=     0.471824600
            v(2,5)=     0.375000000
            v(2,6)=     0.278175400
            v(2,7)=     0.837298300
            v(2,8)=     0.665473800
            v(2,9)=     0.493649200

            !     wk
            v(3,1)=     0.068464390
            v(3,2)=     0.085635710
            v(3,3)=     0.038580250
            v(3,4)=     0.085635710
            v(3,5)=     0.098765430
            v(3,6)=     0.037821090
            v(3,7)=     0.038580250
            v(3,8)=     0.037821090
            v(3,9)=     0.008696116

        case default

    end select



    RETURN
END
