	!*****************************************************************************
    !                                       ____  _____  
    !           /\                         |___ \|  __ \ 
    !          /  \   _ __  _   _ _ __ __ _  __) | |  | |
    !         / /\ \ | '_ \| | | | '__/ _` ||__ <| |  | |
    !        / ____ \| | | | |_| | | | (_| |___) | |__| |
    !       /_/    \_\_| |_|\__,_|_|  \__,_|____/|_____/ 
    !
    !
	!	Anura3D - Numerical modelling and simulation of large deformations 
    !   and soil–water–structure interaction using the material point method (MPM)
    !
    !	Copyright (C) 2023  Members of the Anura3D MPM Research Community 
    !   (See Contributors file "Contributors.txt")
    !
    !	This program is free software: you can redistribute it and/or modify
    !	it under the terms of the GNU Lesser General Public License as published by
    !	the Free Software Foundation, either version 3 of the License, or
    !	(at your option) any later version.
    !
    !	This program is distributed in the hope that it will be useful,
    !	but WITHOUT ANY WARRANTY; without even the implied warranty of
    !	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    !	GNU Lesser General Public License for more details.
    !
    !	You should have received a copy of the GNU Lesser General Public License
    !	along with this program.  If not, see <https://www.gnu.org/licenses/>.
	!
    !*****************************************************************************


! Module GetStrain
!**********************************************************************
!    Should this be a module 
!
!     $Revision: 8842 $
!     $Date: 2020-07-30 07:58:40 -0400 (Thu, 30 Jul 2020) $
!
!**********************************************************************

  
  
    
    subroutine Get_Strain(IEl, IPoint, ICon, B, Disp, NDofEx, Eps)
    !*************************************************************************************   
    !    subroutine:     Get_Strain
    ! 
    !    DESCRIPTION:        
    !>   Calculates the strains in an integration point
    !
    !>   @param[in] IEl : Element number
    !    @param[in] IMP : Material point number
    !    @param[in] Icon : Element conectivities
    !    @param[in] B : B matrix
    !    @param[in] disp : Displacements of global dof
    !    @param[in] NDofEx : Reduced dof
    !
    !>   @param[out] Eps : Strain
    !
    !*************************************************************************************
      use ModCounters
      use ModMeshInfo
      use ModMPMData
      use ModReadCalculationData
      
      implicit none
      
      integer(INTEGER_TYPE), intent(in) :: IEl, IPoint
      integer(INTEGER_TYPE), dimension(ELEMENTNODES, Counters%Nel), intent(in) :: ICon
      real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES), intent(in) :: B
      real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES):: B_
      real(REAL_TYPE), dimension(Counters%N), intent(in) :: Disp
      integer(INTEGER_TYPE), dimension(Counters%NodTot+1), intent(in) :: NDofEx
      real(REAL_TYPE), dimension(NTENSOR), intent(out) :: Eps
      
      ! local variables
      integer(INTEGER_TYPE) :: J, NN, ParticleIndex, K, I
      real(REAL_TYPE) :: Position,Position_
      real(REAL_TYPE), dimension(NVECTOR) :: U
      real(REAL_TYPE), dimension(ELEMENTNODES) :: ShapeValues,ShapeValues_
      real(REAL_TYPE) :: RJAC(NDIM, NDIM), RJAC1(NDIM, NDIM),DET, DET1, H
      
      Eps = 0.0

      select case(NDIM)
      
        case(2)
          
          if ( ISAXISYMMETRIC ) then
            if ( IsParticleIntegration(IEl) ) then ! MP-integration
              ParticleIndex = GetParticleIndex(IPoint, IEl) ! MP global ID
              Position = GlobPosArray(ParticleIndex, 1) ! index 1 is radial direction
              ShapeValues(:) = ShapeValuesArray(ParticleIndex, :)
            else ! GP-integration
                if( Number_of_GP(IEL) == 4 .AND. ELEMENTTYPE == QUAD4) then
                    Position = GP4GlobalPositionElement(1, IPoint, IEl) ! index 1 is radial direction
                    ShapeValues(:) = GP4ShapeFunction(IPoint, :)
                else
                    Position = GPGlobalPositionElement(1, IPoint, IEl) ! index 1 is radial direction
                    ShapeValues(:) = GPShapeFunction(IPoint, :)
                endif
            end if
          end if 
          
          if (.NOT.IsParticleIntegration(IEl).AND.CalParams%ApplyB_BAR .AND. Number_of_GP(IEL) == 4 .AND. ELEMENTTYPE == QUAD4) then
              RJac = 0.0
              Position_ = GPGlobalPositionElement(1, 1, IEl) ! index 1 is radial direction
              ShapeValues_(:) = GPShapeFunction(1, :)
              do K = 1, ELEMENTNODES
                  NN = ICon(K, IEl)
                  do I = 1, NDIM
                      do J = 1, NDIM
                          RJac(I,J) = RJac(I,J) + GPShapeFunctionDerivative(1,K,I) * NodalCoordinatesUpd(NN,J)
                      end do
                  end do
              end do
              
              ! Invert Jacobi-matrix
              call RJacInv(NDIM, RJac, RJac1, Det, Det1)
              
              do K = 1, ELEMENTNODES
                  do I = 1, NDIM
                      H = 0.0
                      do J = 1, NDIM
                          H = H + RJac1(I,J) * GPShapeFunctionDerivative(1,K,J)
                     end do
                      B_(I,K) = H
                  end do
              end do
              
              do J = 1, ELEMENTNODES 
                  NN = Icon(J, Iel)
                  U(1) = Disp(NDofEx(NN)+1)
                  U(2) = Disp(NDofEx(NN)+2)
                  Eps(1) = Eps(1) + (B(1,J) * 2.0 + B_(1,J) + ShapeValues_(J)/Position_ - ShapeValues(J) / Position) / 3.0 * U(1) + (B_(2,J) - B(2,J))/3.0 * U(2)                 ! Eps_XX = dUx/dX
                  Eps(2) = Eps(2) + (B_(1,J) - B(1,J) + ShapeValues_(J)/Position_ - ShapeValues(J) / Position) / 3.0 * U(1) + (B_(2,J) + 2.0 * B(2,J))/3.0 * U(2)                 ! Eps_YY = dUy/dY
                  if ( ISAXISYMMETRIC ) then
                      Eps(3) = Eps(3) + (B_(1,J) - B(1,J) + ShapeValues_(J)/Position_ + ShapeValues(J) / Position * 2.0) / 3.0 * U(1) + (B_(2,J) - B(2,J))/3.0 * U(2) 
                  else               
                      Eps(3) = 0.0                                     ! E_zz = 0 in 2D
                  end if
                  Eps(4) = Eps(4) + B(2,J) * U(1) + B(1,J) * U(2)  ! Gam_XY = dUx/dY+dUy/dX
              end do
          else
              do J = 1, ELEMENTNODES 
                  NN = Icon(J, Iel)
                  U(1) = Disp(NDofEx(NN)+1)
                  U(2) = Disp(NDofEx(NN)+2)
                  Eps(1) = Eps(1) + B(1,J) * U(1)                  ! Eps_XX = dUx/dX
                  Eps(2) = Eps(2) + B(2,J) * U(2)                  ! Eps_YY = dUy/dY
                  if ( ISAXISYMMETRIC ) then
                      Eps(3) = Eps(3) + U(1) * ShapeValues(J) / Position ! Eps_tt = u_r*N/r
                  else               
                      Eps(3) = 0.0                                     ! E_zz = 0 in 2D
                  end if
                  Eps(4) = Eps(4) + B(2,J) * U(1) + B(1,J) * U(2)  ! Gam_XY = dUx/dY+dUy/dX
              end do  
          endif
              
            
        case(3)
            
          do J = 1, ELEMENTNODES 
            NN = Icon(J, Iel)
            U(1) = Disp(NDofEx(NN)+1)
            U(2) = Disp(NDofEx(NN)+2)
            U(3) = Disp(NDofEx(NN)+3)
            Eps(1) = Eps(1) + B(1,J) * U(1)                  ! Eps_XX = dUx/dX
            Eps(2) = Eps(2) + B(2,J) * U(2)                  ! Eps_YY = dUy/dY
            Eps(3) = Eps(3) + B(3,J) * U(3)                  ! Eps_ZZ = dUz/dZ
            Eps(4) = Eps(4) + B(2,J) * U(1) + B(1,J) * U(2)  ! Gam_XY = dUx/dY+dUy/dX
            Eps(5) = Eps(5) + B(3,J) * U(2) + B(2,J) * U(3)  ! Gam_YZ = dUy/dZ+dUz/dY
            Eps(6) = Eps(6) + B(3,J) * U(1) + B(1,J) * U(3)  ! Gam_ZX = dUz/dX+dUx/dZ
          end do  
          
      end select
          
      end subroutine Get_Strain