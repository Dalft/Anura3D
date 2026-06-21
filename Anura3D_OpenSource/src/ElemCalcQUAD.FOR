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


      module ModElementEvaluationQUAD
      !**********************************************************************
      !
      !    DESCRIPTION:
      !    This module provides specific routines for evaluating quadrilateral elements.
      !
      !    $Revision: 7657 $
      !    $Date: 2018-08-15 09:54:30 +0200 (Wed, 15 Aug 2018) $
      !
      !**********************************************************************

      use ModGeometryMath
      use ModString
      use ModReadCalculationData
      use ModGlobalConstants

      implicit none

      
      contains
    

        !**********************************************************************
        !
        !    SUBROUTINE: DetermineSideDataQUAD 
        !
        !    DESCRIPTION: 
        !>   Returns the normal vector of the side with SideID and a PointInitialLocalCoordinatesQUAD
        !>   on that plane. The normal vector is normalised and points inward.
        !    
        !>   @note : 2D element
        !
        !>   @param[in] ISide : Side ID of the element
        !
        !>   @param[out] PlaneNormal : Normal to the element side
        !>   @param[out] PlanePoint : Point on the element side
        !
        !**********************************************************************
        subroutine DetermineSideDataQUAD(ISide, PlaneNormal, PlanePoint)
        
        implicit none
        
          integer(INTEGER_TYPE), parameter :: IDim = 2 ! fixed dimension as 2D element
          
          integer(INTEGER_TYPE), intent(in) :: ISide
          real(REAL_TYPE), dimension(IDim), intent(out) :: PlaneNormal
          real(REAL_TYPE), dimension(IDim), intent(out) :: PlanePoint
          
          PlaneNormal = 0.0
          PlanePoint = 0.0
          
          select case(ISide)
            case(1)
              PlaneNormal(2) = 1.0
            case(2)
              PlaneNormal(1) = 1.0
            case(3)
              PlaneNormal(2) = -1.0
            case(4)
              PlaneNormal(1) = -1.0
            case default
              call GiveError("Undefined side number in [subroutine DetermineSideDataQUAD()].")
          end select
        
        end subroutine DetermineSideDataQUAD


        !**********************************************************************
        !
        !    FUNCTION: PointSideDistanceQUAD
        !
        !    DESCRIPTION:
        !>   Returns the minimum distance between side SideID and LocPos.
        !
        !>   @param[in] SideID : ID of the considered side (1..3)
        !>   @param[in] LocPos : Local coordinates of considered point
        !
        !>   @return PointSideDistanceQuadrilateral
        !
        !**********************************************************************
        real(REAL_TYPE) function PointSideDistanceQUAD(SideID, LocPos)
        
        implicit none
        
          integer(INTEGER_TYPE), parameter :: IDim = 2 ! fixed dimension as 2D element 
          
          integer(INTEGER_TYPE), intent(in) :: SideID
          real(REAL_TYPE), dimension(IDim), intent(in) :: LocPos
          
          ! local variables
          real(REAL_TYPE), dimension(IDim) :: LineNormal, LinePoint

          call DetermineSideDataQUAD(SideID, LineNormal, LinePoint)
        
          PointSideDistanceQUAD = LinePointDistance(LineNormal, LinePoint, LocPos)
        
        end function PointSideDistanceQUAD


        !**********************************************************************
        !
        !    SUBROUTINE: InitialLocalMaterialPointCoordinatesQUAD
        !
        !    DESCRIPTION: 
        !>   Determines the local coordinates and integration weight assigned to material point with ID
        !>   IParticle which are returned through WeiGP and PosGP. Currently, all material points are placed at the same local positions.
        !
        !>   @param[in] IParticle : Number of the material point
        !>   @param[in] SolidPointsElement : Number of solid material points per element
        !>   @param[in] LiquidPointsElement : Number of liquid material points per element
        !>   @param[inout] WeiGP : Initial weight assigned to material point IParticle
        !>   @param[inout] PosGP : Initial local position of material point IParticle
        !
        !**********************************************************************
        subroutine InitialLocalMaterialPointCoordinatesQUAD(IParticle, SolidPointsElement, LiquidPointsElement, WeiGP, PosGP)

        implicit none

          integer(INTEGER_TYPE), intent(in) :: IParticle, SolidPointsElement, LiquidPointsElement
          real(REAL_TYPE), intent(inout) :: WeiGP
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP

          ! local variables
          integer(INTEGER_TYPE) :: ID

             ! first the solid material points are determined
          if ( (IParticle <= SolidPointsElement) .and. (SolidPointsElement > 0) ) then
            ID = IParticle
            select case(SolidPointsElement)
            case (1) ! 1 material point per element
              call InitialQUAD_MP1(ID, PosGP, WeiGP)
            case (4) ! 4 material points per element
              call InitialQUAD_MP4(ID, PosGP, WeiGP)
            case (9) ! 9 material points per element
              call InitialQUAD_MP9(ID, PosGP, WeiGP)
            case (16) ! 16 material points per element
              call InitialQUAD_MP16(ID, PosGP, WeiGP)
            case (25) ! 25 material points per element
              call InitialQUAD_MP25(ID, PosGP, WeiGP)
            case (36) ! 36 material points per element
              call InitialQUAD_MP36(ID, PosGP, WeiGP)
            case (49) ! 49 material points per element
              call InitialQUAD_MP49(ID, PosGP, WeiGP)
            case (64) ! 64 material points per element
              call InitialQUAD_MP64(ID, PosGP, WeiGP)
            case (81) ! 81 material points per element
              call InitialQUAD_MP81(ID, PosGP, WeiGP)
            case (100) ! 100 material points per element
              call InitialQUAD_MP100(ID, PosGP, WeiGP)
            case default
              call GiveError("Number of solid material points, " // trim(String(SolidPointsElement)) // &
                              ", is not available for quadrilateral elements! Supported numbers are 1 and 4. Error in [subroutine InitialLocalMaterialPointCoordinatesQUAD()].")
            end select
          end if
          
          if ( (IParticle > SolidPointsElement) .and. (LiquidPointsElement > 0) ) then
            ID = IParticle - SolidPointsElement
            select case(LiquidPointsElement)
            case (1) ! 1 material point per element
              call InitialQUAD_MP1(ID, PosGP, WeiGP)
            case (4) ! 4 material points per element
              call InitialQUAD_MP4(ID, PosGP, WeiGP)
            case (9) ! 9 material points per element
              call InitialQUAD_MP9(ID, PosGP, WeiGP)
            case (16) ! 16 material points per element
              call InitialQUAD_MP16(ID, PosGP, WeiGP)
            case (25) ! 25 material points per element
              call InitialQUAD_MP25(ID, PosGP, WeiGP)
            case (36) ! 36 material points per element
              call InitialQUAD_MP36(ID, PosGP, WeiGP)
            case (49) ! 49 material points per element
              call InitialQUAD_MP49(ID, PosGP, WeiGP)
            case (64) ! 64 material points per element
              call InitialQUAD_MP64(ID, PosGP, WeiGP)
            case (81) ! 81 material points per element
              call InitialQUAD_MP81(ID, PosGP, WeiGP)
            case (100) ! 100 material points per element
              call InitialQUAD_MP100(ID, PosGP, WeiGP)
            case default
              call GiveError("Number of solid material points, " // trim(String(LiquidPointsElement)) // &
                              ", is not available for quadrilateral elements! Supported numbers are 1 and 4. Error in [subroutine InitialLocalMaterialPointCoordinatesQUAD()].")
            end select
          end if

        end subroutine InitialLocalMaterialPointCoordinatesQUAD


        !**********************************************************************
        !
        !    SUBROUTINE: InitialQUAD_MP1
        !
        !    DESCRIPTION:  
        !>   Returns the initial local coordinates and weight for IParticle.
        !>   1 material point is placed in each element, whose initial
        !>   location and weight is identical with those of the Gauss Point.
        !
        !>   @ note : 2D element
        !>   @ note : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF
        !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987
        !    
        !>   @param[in] IParticle : Local number of the considered particle inside an element
        !
        !>   @param[inout] PosGP : Returns the initial local coordinates of the material point with local number IParticle
        !>   @param[inout] WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************
        subroutine InitialQUAD_MP1(IParticle, PosGP, WeiGP)
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          select case (IParticle)
            case (1)
              PosGP(1) = 0.0
              PosGP(2) = 0.0
              WeiGP = 4.0
            case default
              call GiveError("Undefined number of material points in [subroutine InitialQUAD_MP1()].")   
          end select
        
        end subroutine InitialQUAD_MP1


        !**********************************************************************
        !
        !    SUBROUTINE: InitialQUAD_MP4
        !
        !>   Returns the initial local coordinates and weight for IParticle.
        !>   4 material points are placed in each element, whose initial
        !>   locations and weights are identical with those of Gauss Points.
        !
        !>   @ note : 2D element
        !>   @ note : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF, page 3
        !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987
        !    
        !>   @param[in] IParticle : Local number of the considered particle inside an element
        !
        !>   @param[inout] PosGP : Returns the initial local coordinates of the material point with local number IParticle
        !>   @param[inout] WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************
        subroutine InitialQUAD_MP4(IParticle, PosGP, WeiGP)
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          ! local variables
          real(REAL_TYPE) :: a = 0.57735026918962576450914878050196 ! = 1 / sqrt(3)
          real(REAL_TYPE) :: b = 1.0
          
          select case (IParticle)
            case (1)
              PosGP(1) = -a
              PosGP(2) = -a
              WeiGP = b
            case (2)
              PosGP(1) = a
              PosGP(2) = -a
              WeiGP = b
            case (3)
              PosGP(1) = a
              PosGP(2) = a
              WeiGP = b
            case (4)
              PosGP(1) = -a
              PosGP(2) = a
              WeiGP = b
            case default
              call GiveError("Undefined number of material points in [subroutine InitialQUAD_MP4()].")      
          end select

          PosGP = PosGP * ( 1.0 - CalParams%shrinkageMateriaPointPositionFactor )

        end subroutine InitialQUAD_MP4
        
       
        subroutine InitialQUAD_MP9(IParticle, PosGP, WeiGP)
          !**********************************************************************
        !
        !    SUBROUTINE: InitialQUAD_MP9
        !
        !>   Returns the initial local coordinates and weight for IParticle.
        !>   4 material points are placed in each element, whose initial
        !>   locations and weights are identical with those of Gauss Points.
        !
        !>   @ note : 2D element
        !>   @ note : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF, page 3
        !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987
        !    
        !>   @param[in] IParticle : Local number of the considered particle inside an element
        !
        !>   @param[inout] PosGP : Returns the initial local coordinates of the material point with local number IParticle
        !>   @param[inout] WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          ! local variables
          real(REAL_TYPE) :: a = 0.7746 ! = 1 / sqrt(3)
          real(REAL_TYPE) :: b = 0.308641975308
          real(REAL_TYPE) :: c = 0.493827160
          real(REAL_TYPE) :: d = 0.790123456

          
          select case (IParticle)
            case (1)
              PosGP(1) = -a
              PosGP(2) = -a
              WeiGP = b
            case (2)
              PosGP(1) = 0
              PosGP(2) = -a
              WeiGP = c
            case (3)
              PosGP(1) = a
              PosGP(2) = -a
              WeiGP = b
            case (4)
              PosGP(1) = -a
              PosGP(2) = 0
              WeiGP = c
            case (5)
              PosGP(1) = 0
              PosGP(2) = 0
              WeiGP = d
            case (6)
              PosGP(1) = a
              PosGP(2) = 0
              WeiGP = c
            case (7)
              PosGP(1) = -a
              PosGP(2) = a
              WeiGP = b  
            case (8)
              PosGP(1) = 0
              PosGP(2) = a
              WeiGP = c
            case (9)
              PosGP(1) = a
              PosGP(2) = a
              WeiGP = b    
            case default
              call GiveError("Undefined number of material points in [subroutine InitialQUAD_MP9()].")      
          end select

          PosGP = PosGP * ( 1.0 - CalParams%shrinkageMateriaPointPositionFactor )

        end subroutine InitialQUAD_MP9
        
       subroutine InitialQUAD_MP16(IParticle, PosGP, WeiGP)
    !**********************************************************************
    !
    !    SUBROUTINE: InitialQUAD_MP16
    !
    !>   Returns the initial local coordinates and weight for IParticle.
    !>   16 material points are placed in each element, corresponding to
    !>   4×4 Gauss quadrature points with direct weight values.
    !
    !>   @ note : 2D quadrilateral element
    !>   @ note : Weights are pre-calculated for 4×4 Gauss quadrature
    !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis.
    !    
    !>   @param[in] IParticle : Local number of the particle inside element (1~16)
    !
    !>   @param[inout] PosGP : Local coordinates (ξ, η) of the material point
    !>   @param[inout] WeiGP : Pre-calculated weight of the material point
    ! 
    !**********************************************************************
    implicit none

    integer(INTEGER_TYPE), intent(in) :: IParticle
    real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
    real(REAL_TYPE), intent(inout) :: WeiGP
    
    ! Local variables: 4-point Gauss coordinates (ξ/η ∈ [-1,1])
    real(REAL_TYPE) :: xi1 = -0.75 !-0.8611363115940526  ! 1st Gauss coordinate
    real(REAL_TYPE) :: xi2 = -0.25 !-0.3399810435848563  ! 2nd Gauss coordinate
    real(REAL_TYPE) :: xi3 =  0.25 !0.3399810435848563  ! 3rd Gauss coordinate
    real(REAL_TYPE) :: xi4 =  0.75 !0.8611363115940526  ! 4th Gauss coordinate


    select case (IParticle)
        ! First row (η = xi1)
        case (1)
            PosGP(1) = xi1
            PosGP(2) = xi1
            WeiGP = 0.25!0.1209344246254516  ! w1*w1 = 0.3478548451374538^2
        case (2)
            PosGP(1) = xi2
            PosGP(2) = xi1
            WeiGP = 0.25!0.2274108193056141  ! w2*w1 = 0.6521451548625461*0.3478548451374538
        case (3)
            PosGP(1) = xi3
            PosGP(2) = xi1
            WeiGP = 0.25!0.2274108193056141  ! w3*w1 = 0.6521451548625461*0.3478548451374538
        case (4)
            PosGP(1) = xi4
            PosGP(2) = xi1
            WeiGP = 0.25!0.1209344246254516  ! w4*w1 = 0.3478548451374538^2
            
        ! Second row (η = xi2)
        case (5)
            PosGP(1) = xi1
            PosGP(2) = xi2
            WeiGP = 0.25!0.2274108193056141  ! w1*w2 = 0.3478548451374538*0.6521451548625461
        case (6)
            PosGP(1) = xi2
            PosGP(2) = xi2
            WeiGP = 0.25!0.4252933334489412  ! w2*w2 = 0.6521451548625461^2
        case (7)
            PosGP(1) = xi3
            PosGP(2) = xi2
            WeiGP = 0.25!0.4252933334489412  ! w3*w2 = 0.6521451548625461^2
        case (8)
            PosGP(1) = xi4
            PosGP(2) = xi2
            WeiGP = 0.25!0.2274108193056141  ! w4*w2 = 0.3478548451374538*0.6521451548625461
            
        ! Third row (η = xi3)
        case (9)
            PosGP(1) = xi1
            PosGP(2) = xi3
            WeiGP = 0.25!0.2274108193056141  ! w1*w3 = 0.3478548451374538*0.6521451548625461
        case (10)
            PosGP(1) = xi2
            PosGP(2) = xi3
            WeiGP = 0.25!0.4252933334489412  ! w2*w3 = 0.6521451548625461^2
        case (11)
            PosGP(1) = xi3
            PosGP(2) = xi3
            WeiGP = 0.25!0.4252933334489412  ! w3*w3 = 0.6521451548625461^2
        case (12)
            PosGP(1) = xi4
            PosGP(2) = xi3
            WeiGP = 0.25!0.2274108193056141  ! w4*w3 = 0.3478548451374538*0.6521451548625461
            
        ! Fourth row (η = xi4)
        case (13)
            PosGP(1) = xi1
            PosGP(2) = xi4
            WeiGP = 0.25!0.1209344246254516  ! w1*w4 = 0.3478548451374538^2
        case (14)
            PosGP(1) = xi2
            PosGP(2) = xi4
            WeiGP = 0.25!0.2274108193056141  ! w2*w4 = 0.6521451548625461*0.3478548451374538
        case (15)
            PosGP(1) = xi3
            PosGP(2) = xi4
            WeiGP = 0.25!0.2274108193056141  ! w3*w4 = 0.6521451548625461*0.3478548451374538
        case (16)
            PosGP(1) = xi4
            PosGP(2) = xi4
            WeiGP = 0.25!0.1209344246254516  ! w4*w4 = 0.3478548451374538^2
            
        case default
            call GiveError("IParticle must be 1~16 in [subroutine InitialQUAD_MP16()].")      
    end select

    ! Apply position shrinkage factor
    PosGP = PosGP * (1.0 - CalParams%shrinkageMateriaPointPositionFactor)

    end subroutine InitialQUAD_MP16
    
   subroutine InitialQUAD_MP25(IParticle, PosGP, WeiGP)
  !**********************************************************************
  !
  !    SUBROUTINE: InitialQUAD_MP25
  !
  !>   Returns the initial local coordinates and weight for IParticle.
  !>   25 material points are placed in each element, corresponding to
  !>   5×5 Gauss quadrature points with direct weight values.
  !
  !>   @ note : 2D quadrilateral element
  !>   @ note : Weights are pre-calculated for 5×5 Gauss quadrature
  !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis.
  !    
  !>   @param[in] IParticle : Local number of the particle inside element (1~25)
  !
  !>   @param[inout] PosGP : Local coordinates (ξ, η) of the material point
  !>   @param[inout] WeiGP : Pre-calculated weight of the material point
  ! 
  !**********************************************************************
  implicit none

  integer(INTEGER_TYPE), intent(in) :: IParticle
  real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
  real(REAL_TYPE), intent(inout) :: WeiGP
  
  ! Local variables: 5-point Gauss coordinates (ξ/η ∈ [-1,1])
  real(REAL_TYPE) :: xi1 = -0.9061798459386640  ! 1st Gauss coordinate
  real(REAL_TYPE) :: xi2 = -0.5384693101056831  ! 2nd Gauss coordinate
  real(REAL_TYPE) :: xi3 =  0.0000000000000000  ! 3rd Gauss coordinate
  real(REAL_TYPE) :: xi4 =  0.5384693101056831  ! 4th Gauss coordinate
  real(REAL_TYPE) :: xi5 =  0.9061798459386640  ! 5th Gauss coordinate

  ! 5-point Gauss weights
  real(REAL_TYPE) :: w1 = 0.2369268850561891  ! Weight for xi1 and xi5
  real(REAL_TYPE) :: w2 = 0.4786286704993665  ! Weight for xi2 and xi4
  real(REAL_TYPE) :: w3 = 0.5688888888888888  ! Weight for xi3


  select case (IParticle)
      ! First row (η = xi1)
      case (1)
          PosGP(1) = xi1
          PosGP(2) = xi1
          WeiGP = w1 * w1  ! 0.0561298701638075
      case (2)
          PosGP(1) = xi2
          PosGP(2) = xi1
          WeiGP = w2 * w1  ! 0.1134804256851841
      case (3)
          PosGP(1) = xi3
          PosGP(2) = xi1
          WeiGP = w3 * w1  ! 0.1347106603888530
      case (4)
          PosGP(1) = xi4
          PosGP(2) = xi1
          WeiGP = w2 * w1  ! 0.1134804256851841
      case (5)
          PosGP(1) = xi5
          PosGP(2) = xi1
          WeiGP = w1 * w1  ! 0.0561298701638075
          
      ! Second row (η = xi2)
      case (6)
          PosGP(1) = xi1
          PosGP(2) = xi2
          WeiGP = w1 * w2  ! 0.1134804256851841
      case (7)
          PosGP(1) = xi2
          PosGP(2) = xi2
          WeiGP = w2 * w2  ! 0.2290616641358642
      case (8)
          PosGP(1) = xi3
          PosGP(2) = xi2
          WeiGP = w3 * w2  ! 0.2694213207777060
      case (9)
          PosGP(1) = xi4
          PosGP(2) = xi2
          WeiGP = w2 * w2  ! 0.2290616641358642
      case (10)
          PosGP(1) = xi5
          PosGP(2) = xi2
          WeiGP = w1 * w2  ! 0.1134804256851841
          
      ! Third row (η = xi3)
      case (11)
          PosGP(1) = xi1
          PosGP(2) = xi3
          WeiGP = w1 * w3  ! 0.1347106603888530
      case (12)
          PosGP(1) = xi2
          PosGP(2) = xi3
          WeiGP = w2 * w3  ! 0.2694213207777060
      case (13)
          PosGP(1) = xi3
          PosGP(2) = xi3
          WeiGP = w3 * w3  ! 0.3237367056887454
      case (14)
          PosGP(1) = xi4
          PosGP(2) = xi3
          WeiGP = w2 * w3  ! 0.2694213207777060
      case (15)
          PosGP(1) = xi5
          PosGP(2) = xi3
          WeiGP = w1 * w3  ! 0.1347106603888530
          
      ! Fourth row (η = xi4)
      case (16)
          PosGP(1) = xi1
          PosGP(2) = xi4
          WeiGP = w1 * w2  ! 0.1134804256851841
      case (17)
          PosGP(1) = xi2
          PosGP(2) = xi4
          WeiGP = w2 * w2  ! 0.2290616641358642
      case (18)
          PosGP(1) = xi3
          PosGP(2) = xi4
          WeiGP = w3 * w2  ! 0.2694213207777060
      case (19)
          PosGP(1) = xi4
          PosGP(2) = xi4
          WeiGP = w2 * w2  ! 0.2290616641358642
      case (20)
          PosGP(1) = xi5
          PosGP(2) = xi4
          WeiGP = w1 * w2  ! 0.1134804256851841
          
      ! Fifth row (η = xi5)
      case (21)
          PosGP(1) = xi1
          PosGP(2) = xi5
          WeiGP = w1 * w1  ! 0.0561298701638075
      case (22)
          PosGP(1) = xi2
          PosGP(2) = xi5
          WeiGP = w2 * w1  ! 0.1134804256851841
      case (23)
          PosGP(1) = xi3
          PosGP(2) = xi5
          WeiGP = w3 * w1  ! 0.1347106603888530
      case (24)
          PosGP(1) = xi4
          PosGP(2) = xi5
          WeiGP = w2 * w1  ! 0.1134804256851841
      case (25)
          PosGP(1) = xi5
          PosGP(2) = xi5
          WeiGP = w1 * w1  ! 0.0561298701638075
          
      case default
          call GiveError("IParticle must be 1~25 in [subroutine InitialQUAD_MP25()].")      
  end select

  ! Apply position shrinkage factor
  PosGP = PosGP * (1.0 - CalParams%shrinkageMateriaPointPositionFactor)

end subroutine InitialQUAD_MP25
    
    subroutine InitialQUAD_MP36(IParticle, PosGP, WeiGP)
  !**********************************************************************
  !
  !    SUBROUTINE: InitialQUAD_MP36
  !
  !>   Returns the initial local coordinates and weight for IParticle.
  !>   36 material points are placed in each element, corresponding to
  !>   6×6 Gauss quadrature points with direct weight values.
  !
  !>   @ note : 2D quadrilateral element
  !>   @ note : Weights are pre-calculated for 6×6 Gauss quadrature
  !>   @ note : T.J.R. Hughes. The Finite Element Method reference.
  !    
  !>   @param[in] IParticle : Local number of the particle inside element (1~36)
  !
  !>   @param[inout] PosGP : Local coordinates (ξ, η) of the material point
  !>   @param[inout] WeiGP : Pre-calculated weight of the material point
  ! 
  !**********************************************************************
  implicit none

  integer(INTEGER_TYPE), intent(in) :: IParticle
  real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
  real(REAL_TYPE), intent(inout) :: WeiGP
  
  ! Local variables: 6-point Gauss coordinates (ξ/η ∈ [-1,1])
  real(REAL_TYPE) :: xi1 = -0.9324695142031521  ! 1st Gauss coordinate
  real(REAL_TYPE) :: xi2 = -0.6612093864662645  ! 2nd Gauss coordinate
  real(REAL_TYPE) :: xi3 = -0.2386191860831969  ! 3rd Gauss coordinate
  real(REAL_TYPE) :: xi4 =  0.2386191860831969  ! 4th Gauss coordinate
  real(REAL_TYPE) :: xi5 =  0.6612093864662645  ! 5th Gauss coordinate
  real(REAL_TYPE) :: xi6 =  0.9324695142031521  ! 6th Gauss coordinate

  ! 6-point Gauss weights
  real(REAL_TYPE) :: w1 = 0.1713244923791704  ! Weight for xi1 and xi6
  real(REAL_TYPE) :: w2 = 0.3607615730481384  ! Weight for xi2 and xi5
  real(REAL_TYPE) :: w3 = 0.4679139345726910  ! Weight for xi3 and xi4


  select case (IParticle)
      ! First row (η = xi1)
      case (1)
          PosGP(1) = xi1; PosGP(2) = xi1; WeiGP = w1*w1  ! 0.0293580050130844
      case (2)
          PosGP(1) = xi2; PosGP(2) = xi1; WeiGP = w2*w1  ! 0.0618733246021280
      case (3)
          PosGP(1) = xi3; PosGP(2) = xi1; WeiGP = w3*w1  ! 0.0791567635003250
      case (4)
          PosGP(1) = xi4; PosGP(2) = xi1; WeiGP = w3*w1  ! 0.0791567635003250
      case (5)
          PosGP(1) = xi5; PosGP(2) = xi1; WeiGP = w2*w1  ! 0.0618733246021280
      case (6)
          PosGP(1) = xi6; PosGP(2) = xi1; WeiGP = w1*w1  ! 0.0293580050130844
          
      ! Second row (η = xi2)
      case (7)
          PosGP(1) = xi1; PosGP(2) = xi2; WeiGP = w1*w2  ! 0.0618733246021280
      case (8)
          PosGP(1) = xi2; PosGP(2) = xi2; WeiGP = w2*w2  ! 0.1299488289148246
      case (9)
          PosGP(1) = xi3; PosGP(2) = xi2; WeiGP = w3*w2  ! 0.1684834349158386
      case (10)
          PosGP(1) = xi4; PosGP(2) = xi2; WeiGP = w3*w2  ! 0.1684834349158386
      case (11)
          PosGP(1) = xi5; PosGP(2) = xi2; WeiGP = w2*w2  ! 0.1299488289148246
      case (12)
          PosGP(1) = xi6; PosGP(2) = xi2; WeiGP = w1*w2  ! 0.0618733246021280
          
      ! Third row (η = xi3)
      case (13)
          PosGP(1) = xi1; PosGP(2) = xi3; WeiGP = w1*w3  ! 0.0791567635003250
      case (14)
          PosGP(1) = xi2; PosGP(2) = xi3; WeiGP = w2*w3  ! 0.1684834349158386
      case (15)
          PosGP(1) = xi3; PosGP(2) = xi3; WeiGP = w3*w3  ! 0.2189558134556210
      case (16)
          PosGP(1) = xi4; PosGP(2) = xi3; WeiGP = w3*w3  ! 0.2189558134556210
      case (17)
          PosGP(1) = xi5; PosGP(2) = xi3; WeiGP = w2*w3  ! 0.1684834349158386
      case (18)
          PosGP(1) = xi6; PosGP(2) = xi3; WeiGP = w1*w3  ! 0.0791567635003250
          
      ! Fourth row (η = xi4)
      case (19)
          PosGP(1) = xi1; PosGP(2) = xi4; WeiGP = w1*w3  ! 0.0791567635003250
      case (20)
          PosGP(1) = xi2; PosGP(2) = xi4; WeiGP = w2*w3  ! 0.1684834349158386
      case (21)
          PosGP(1) = xi3; PosGP(2) = xi4; WeiGP = w3*w3  ! 0.2189558134556210
      case (22)
          PosGP(1) = xi4; PosGP(2) = xi4; WeiGP = w3*w3  ! 0.2189558134556210
      case (23)
          PosGP(1) = xi5; PosGP(2) = xi4; WeiGP = w2*w3  ! 0.1684834349158386
      case (24)
          PosGP(1) = xi6; PosGP(2) = xi4; WeiGP = w1*w3  ! 0.0791567635003250
          
      ! Fifth row (η = xi5)
      case (25)
          PosGP(1) = xi1; PosGP(2) = xi5; WeiGP = w1*w2  ! 0.0618733246021280
      case (26)
          PosGP(1) = xi2; PosGP(2) = xi5; WeiGP = w2*w2  ! 0.1299488289148246
      case (27)
          PosGP(1) = xi3; PosGP(2) = xi5; WeiGP = w3*w2  ! 0.1684834349158386
      case (28)
          PosGP(1) = xi4; PosGP(2) = xi5; WeiGP = w3*w2  ! 0.1684834349158386
      case (29)
          PosGP(1) = xi5; PosGP(2) = xi5; WeiGP = w2*w2  ! 0.1299488289148246
      case (30)
          PosGP(1) = xi6; PosGP(2) = xi5; WeiGP = w1*w2  ! 0.0618733246021280
          
      ! Sixth row (η = xi6)
      case (31)
          PosGP(1) = xi1; PosGP(2) = xi6; WeiGP = w1*w1  ! 0.0293580050130844
      case (32)
          PosGP(1) = xi2; PosGP(2) = xi6; WeiGP = w2*w1  ! 0.0618733246021280
      case (33)
          PosGP(1) = xi3; PosGP(2) = xi6; WeiGP = w3*w1  ! 0.0791567635003250
      case (34)
          PosGP(1) = xi4; PosGP(2) = xi6; WeiGP = w3*w1  ! 0.0791567635003250
      case (35)
          PosGP(1) = xi5; PosGP(2) = xi6; WeiGP = w2*w1  ! 0.0618733246021280
      case (36)
          PosGP(1) = xi6; PosGP(2) = xi6; WeiGP = w1*w1  ! 0.0293580050130844
          
      case default
          call GiveError("IParticle must be 1~36 in [subroutine InitialQUAD_MP36()].")      
  end select

  ! Apply position shrinkage factor
  PosGP = PosGP * (1.0 - CalParams%shrinkageMateriaPointPositionFactor)

end subroutine InitialQUAD_MP36
    
    subroutine InitialQUAD_MP49(IParticle, PosGP, WeiGP)
  !**********************************************************************
  !
  !    SUBROUTINE: InitialQUAD_MP49
  !
  !>   Returns the initial local coordinates and weight for IParticle.
  !>   49 material points are placed in each element, corresponding to
  !>   7×7 Gauss quadrature points with direct weight values.
  !
  !>   @ note : 2D quadrilateral element
  !>   @ note : Weights are pre-calculated for 7×7 Gauss quadrature
  !>   @ note : T.J.R. Hughes. The Finite Element Method reference.
  !    
  !>   @param[in] IParticle : Local number of the particle inside element (1~49)
  !
  !>   @param[inout] PosGP : Local coordinates (ξ, η) of the material point
  !>   @param[inout] WeiGP : Pre-calculated weight of the material point
  ! 
  !**********************************************************************
  implicit none

  integer(INTEGER_TYPE), intent(in) :: IParticle
  real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
  real(REAL_TYPE), intent(inout) :: WeiGP
  
  ! Local variables: 7-point Gauss coordinates (ξ/η ∈ [-1,1])
  real(REAL_TYPE) :: xi1 = -0.9491079123427585  ! 1st Gauss coordinate
  real(REAL_TYPE) :: xi2 = -0.7415311855993945  ! 2nd Gauss coordinate
  real(REAL_TYPE) :: xi3 = -0.4058451513773972  ! 3rd Gauss coordinate
  real(REAL_TYPE) :: xi4 =  0.0000000000000000  ! 4th Gauss coordinate
  real(REAL_TYPE) :: xi5 =  0.4058451513773972  ! 5th Gauss coordinate
  real(REAL_TYPE) :: xi6 =  0.7415311855993945  ! 6th Gauss coordinate
  real(REAL_TYPE) :: xi7 =  0.9491079123427585  ! 7th Gauss coordinate

  ! 7-point Gauss weights
  real(REAL_TYPE) :: w1 = 0.1294849661688697  ! Weight for xi1 and xi7
  real(REAL_TYPE) :: w2 = 0.2797053914892767  ! Weight for xi2 and xi6
  real(REAL_TYPE) :: w3 = 0.3818300505051189  ! Weight for xi3 and xi5
  real(REAL_TYPE) :: w4 = 0.4179591836734694  ! Weight for xi4


  select case (IParticle)
      ! First row (η = xi1)
      case (1)
          PosGP(1) = xi1; PosGP(2) = xi1; WeiGP = w1*w1  ! 0.0167665641133612
      case (2)
          PosGP(1) = xi2; PosGP(2) = xi1; WeiGP = w2*w1  ! 0.0362282663250433
      case (3)
          PosGP(1) = xi3; PosGP(2) = xi1; WeiGP = w3*w1  ! 0.0504589442746111
      case (4)
          PosGP(1) = xi4; PosGP(2) = xi1; WeiGP = w4*w1  ! 0.0541035169265112
      case (5)
          PosGP(1) = xi5; PosGP(2) = xi1; WeiGP = w3*w1  ! 0.0504589442746111
      case (6)
          PosGP(1) = xi6; PosGP(2) = xi1; WeiGP = w2*w1  ! 0.0362282663250433
      case (7)
          PosGP(1) = xi7; PosGP(2) = xi1; WeiGP = w1*w1  ! 0.0167665641133612
          
      ! Second row (η = xi2)
      case (8)
          PosGP(1) = xi1; PosGP(2) = xi2; WeiGP = w1*w2  ! 0.0362282663250433
      case (9)
          PosGP(1) = xi2; PosGP(2) = xi2; WeiGP = w2*w2  ! 0.0782357442268667
      case (10)
          PosGP(1) = xi3; PosGP(2) = xi2; WeiGP = w3*w2  ! 0.1078331923136667
      case (11)
          PosGP(1) = xi4; PosGP(2) = xi2; WeiGP = w4*w2  ! 0.1161338263664066
      case (12)
          PosGP(1) = xi5; PosGP(2) = xi2; WeiGP = w3*w2  ! 0.1078331923136667
      case (13)
          PosGP(1) = xi6; PosGP(2) = xi2; WeiGP = w2*w2  ! 0.0782357442268667
      case (14)
          PosGP(1) = xi7; PosGP(2) = xi2; WeiGP = w1*w2  ! 0.0362282663250433
          
      ! Third row (η = xi3)
      case (15)
          PosGP(1) = xi1; PosGP(2) = xi3; WeiGP = w1*w3  ! 0.0504589442746111
      case (16)
          PosGP(1) = xi2; PosGP(2) = xi3; WeiGP = w2*w3  ! 0.1078331923136667
      case (17)
          PosGP(1) = xi3; PosGP(2) = xi3; WeiGP = w3*w3  ! 0.1457918006230693
      case (18)
          PosGP(1) = xi4; PosGP(2) = xi3; WeiGP = w4*w3  ! 0.1595439325552063
      case (19)
          PosGP(1) = xi5; PosGP(2) = xi3; WeiGP = w3*w3  ! 0.1457918006230693
      case (20)
          PosGP(1) = xi6; PosGP(2) = xi3; WeiGP = w2*w3  ! 0.1078331923136667
      case (21)
          PosGP(1) = xi7; PosGP(2) = xi3; WeiGP = w1*w3  ! 0.0504589442746111
          
      ! Fourth row (η = xi4)
      case (22)
          PosGP(1) = xi1; PosGP(2) = xi4; WeiGP = w1*w4  ! 0.0541035169265112
      case (23)
          PosGP(1) = xi2; PosGP(2) = xi4; WeiGP = w2*w4  ! 0.1161338263664066
      case (24)
          PosGP(1) = xi3; PosGP(2) = xi4; WeiGP = w3*w4  ! 0.1595439325552063
      case (25)
          PosGP(1) = xi4; PosGP(2) = xi4; WeiGP = w4*w4  ! 0.1747052624673085
      case (26)
          PosGP(1) = xi5; PosGP(2) = xi4; WeiGP = w3*w4  ! 0.1595439325552063
      case (27)
          PosGP(1) = xi6; PosGP(2) = xi4; WeiGP = w2*w4  ! 0.1161338263664066
      case (28)
          PosGP(1) = xi7; PosGP(2) = xi4; WeiGP = w1*w4  ! 0.0541035169265112
          
      ! Fifth row (η = xi5)
      case (29)
          PosGP(1) = xi1; PosGP(2) = xi5; WeiGP = w1*w3  ! 0.0504589442746111
      case (30)
          PosGP(1) = xi2; PosGP(2) = xi5; WeiGP = w2*w3  ! 0.1078331923136667
      case (31)
          PosGP(1) = xi3; PosGP(2) = xi5; WeiGP = w3*w3  ! 0.1457918006230693
      case (32)
          PosGP(1) = xi4; PosGP(2) = xi5; WeiGP = w4*w3  ! 0.1595439325552063
      case (33)
          PosGP(1) = xi5; PosGP(2) = xi5; WeiGP = w3*w3  ! 0.1457918006230693
      case (34)
          PosGP(1) = xi6; PosGP(2) = xi5; WeiGP = w2*w3  ! 0.1078331923136667
      case (35)
          PosGP(1) = xi7; PosGP(2) = xi5; WeiGP = w1*w3  ! 0.0504589442746111
          
      ! Sixth row (η = xi6)
      case (36)
          PosGP(1) = xi1; PosGP(2) = xi6; WeiGP = w1*w2  ! 0.0362282663250433
      case (37)
          PosGP(1) = xi2; PosGP(2) = xi6; WeiGP = w2*w2  ! 0.0782357442268667
      case (38)
          PosGP(1) = xi3; PosGP(2) = xi6; WeiGP = w3*w2  ! 0.1078331923136667
      case (39)
          PosGP(1) = xi4; PosGP(2) = xi6; WeiGP = w4*w2  ! 0.1161338263664066
      case (40)
          PosGP(1) = xi5; PosGP(2) = xi6; WeiGP = w3*w2  ! 0.1078331923136667
      case (41)
          PosGP(1) = xi6; PosGP(2) = xi6; WeiGP = w2*w2  ! 0.0782357442268667
      case (42)
          PosGP(1) = xi7; PosGP(2) = xi6; WeiGP = w1*w2  ! 0.0362282663250433
          
      ! Seventh row (η = xi7)
      case (43)
          PosGP(1) = xi1; PosGP(2) = xi7; WeiGP = w1*w1  ! 0.0167665641133612
      case (44)
          PosGP(1) = xi2; PosGP(2) = xi7; WeiGP = w2*w1  ! 0.0362282663250433
      case (45)
          PosGP(1) = xi3; PosGP(2) = xi7; WeiGP = w3*w1  ! 0.0504589442746111
      case (46)
          PosGP(1) = xi4; PosGP(2) = xi7; WeiGP = w4*w1  ! 0.0541035169265112
      case (47)
          PosGP(1) = xi5; PosGP(2) = xi7; WeiGP = w3*w1  ! 0.0504589442746111
      case (48)
          PosGP(1) = xi6; PosGP(2) = xi7; WeiGP = w2*w1  ! 0.0362282663250433
      case (49)
          PosGP(1) = xi7; PosGP(2) = xi7; WeiGP = w1*w1  ! 0.0167665641133612
          
      case default
          call GiveError("IParticle must be 1~49 in [subroutine InitialQUAD_MP49()].")      
  end select

  ! Apply position shrinkage factor
  PosGP = PosGP * (1.0 - CalParams%shrinkageMateriaPointPositionFactor)

end subroutine InitialQUAD_MP49
    
       subroutine InitialQUAD_MP64(IParticle, PosGP, WeiGP)
  !**********************************************************************
  !
  !    SUBROUTINE: InitialQUAD_MP64
  !
  !>   Returns the initial local coordinates and weight for IParticle.
  !>   64 material points are placed in each element, corresponding to
  !>   8×8 Gauss quadrature points with direct weight values.
  !
  !>   @ note : 2D quadrilateral element
  !>   @ note : Weights are pre-calculated for 8×8 Gauss quadrature
  !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis.
  !    
  !>   @param[in] IParticle : Local number of the particle inside element (1~64)
  !
  !>   @param[inout] PosGP : Local coordinates (ξ, η) of the material point
  !>   @param[inout] WeiGP : Pre-calculated weight of the material point
  ! 
  !**********************************************************************
  implicit none

  integer(INTEGER_TYPE), intent(in) :: IParticle
  real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
  real(REAL_TYPE), intent(inout) :: WeiGP
  
  ! Local variables: 8-point Gauss coordinates (ξ/η ∈ [-1,1])
  real(REAL_TYPE) :: xi(8)
  real(REAL_TYPE) :: w(8)
  integer :: i, j, temp
  
  ! 8-point Gauss quadrature coordinates and weights
  !xi(1) = -0.9602898564975363  
  !xi(2) = -0.7966664774136267  
  !xi(3) = -0.5255324099163290  
  !xi(4) = -0.1834346424956498  
  !xi(5) =  0.1834346424956498  
  !xi(6) =  0.5255324099163290  
  !xi(7) =  0.7966664774136267  
  !xi(8) =  0.9602898564975363  
  
  ! 8-point Gauss quadrature weights
  !w(1) = 0.1012285362903763   
  !w(2) = 0.2223810344533745   
  !w(3) = 0.3137066458778873   
  !w(4) = 0.3626837833783619   
  !w(5) = 0.3626837833783619   
  !w(6) = 0.3137066458778873   
  !w(7) = 0.2223810344533745   
  !w(8) = 0.1012285362903763
  xi(1) = -0.875  
  xi(2) = -0.625  
  xi(3) = -0.375
  xi(4) = -0.125  
  xi(5) =  0.125
  xi(6) =  0.375  
  xi(7) =  0.625  
  xi(8) =  0.875
  ! 8-point uniform weights (所有权重统一，保证64个点权重相等)
  w(1) = 0.25   
  w(2) = 0.25   
  w(3) = 0.25   
  w(4) = 0.25  
  w(5) = 0.25   
  w(6) = 0.25   
  w(7) = 0.25   
  w(8) = 0.25  

  ! Calculate i and j indices (1-8) from IParticle (1-64)
 select case (IParticle)
 case(1:64)
     j = (IParticle - 1) / 8 + 1  ! η direction index (row)
  
  ! Calculate i without using mod function
     temp = IParticle - (j - 1) * 8
     i = temp  ! ξ direction index (column)
  
  ! Assign coordinates based on indices
     PosGP(1) = xi(i)  ! ξ coordinate
     PosGP(2) = xi(j)  ! η coordinate
  
  ! Calculate weight as product of individual weights
     WeiGP = w(i) * w(j)
  
  ! Apply position shrinkage factor
     PosGP = PosGP * (1.0 - CalParams%shrinkageMateriaPointPositionFactor)
  endselect

 end subroutine InitialQUAD_MP64
    
    subroutine InitialQUAD_MP81(IParticle, PosGP, WeiGP)
  !**********************************************************************
  !
  !    SUBROUTINE: InitialQUAD_MP81
  !
  !>   Returns the initial local coordinates and weight for IParticle.
  !>   81 material points are placed in each element, corresponding to
  !>   9×9 Gauss quadrature points with direct weight values.
  !
  !>   @ note : 2D quadrilateral element
  !>   @ note : Weights are pre-calculated for 9×9 Gauss quadrature
  !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis.
  !    
  !>   @param[in] IParticle : Local number of the particle inside element (1~81)
  !
  !>   @param[inout] PosGP : Local coordinates (ξ, η) of the material point
  !>   @param[inout] WeiGP : Pre-calculated weight of the material point
  ! 
  !**********************************************************************
  implicit none

  integer(INTEGER_TYPE), intent(in) :: IParticle
  real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
  real(REAL_TYPE), intent(inout) :: WeiGP
  
  ! Local variables: 9-point Gauss coordinates (ξ/η ∈ [-1,1]) and weights
  real(REAL_TYPE) :: xi(9)  ! Gauss coordinates (symmetric distribution)
  real(REAL_TYPE) :: w(9)   ! Corresponding Gauss weights
  integer :: i, j, temp     ! Indices for ξ(column)、η(row) and temp calculation
  
  ! 9-point Gauss quadrature coordinates (symmetric about 0)
  xi(1) = -0.9602898564975363   ! 1st symmetric coordinate
  xi(2) = -0.7966664774136267   ! 2nd symmetric coordinate
  xi(3) = -0.5255324099163290   ! 3rd symmetric coordinate
  xi(4) = -0.1834346424956498   ! 4th symmetric coordinate
  xi(5) =  0.0000000000000000   ! Central coordinate (symmetric center)
  xi(6) =  0.1834346424956498   ! 4th symmetric coordinate (mirror of xi(4))
  xi(7) =  0.5255324099163290   ! 3rd symmetric coordinate (mirror of xi(3))
  xi(8) =  0.7966664774136267   ! 2nd symmetric coordinate (mirror of xi(2))
  xi(9) =  0.9602898564975363   ! 1st symmetric coordinate (mirror of xi(1))
  
  ! 9-point Gauss quadrature weights (matching coordinates, symmetric)
  w(1) =  0.0812743883863977    ! Weight for xi(1) and xi(9)
  w(2) =  0.1806481606948574    ! Weight for xi(2) and xi(8)
  w(3) =  0.2606106964029354    ! Weight for xi(3) and xi(7)
  w(4) =  0.3123475953466389    ! Weight for xi(4) and xi(6)
  w(5) =  0.3302393550012598    ! Weight for central xi(5)
  w(6) =  0.3123475953466389    ! Weight for xi(6) (same as xi(4))
  w(7) =  0.2606106964029354    ! Weight for xi(7) (same as xi(3))
  w(8) =  0.1806481606948574    ! Weight for xi(8) (same as xi(2))
  w(9) =  0.0812743883863977    ! Weight for xi(9) (same as xi(1))

  ! Calculate i(ξ direction, 1-9) and j(η direction, 1-9) from IParticle(1-81)
  select case (IParticle)
  case(1:81)
      j = (IParticle - 1) / 9 + 1   ! η direction index (row: 1-9, 9 points per row)
  
      ! Calculate ξ direction index without mod function (avoid negative results)
      temp = IParticle - (j - 1) * 9
      i = temp                      ! ξ direction index (column: 1-9)
  
      ! Assign local coordinates (ξ from i, η from j)
      PosGP(1) = xi(i)  ! ξ coordinate (horizontal direction)
      PosGP(2) = xi(j)  ! η coordinate (vertical direction)
  
      ! Calculate 2D Gauss weight (product of 1D weights for ξ and η)
      WeiGP = w(i) * w(j)
  
      ! Apply material point position shrinkage factor (consistent with MP100 version)
      PosGP = PosGP * (1.0 - CalParams%shrinkageMateriaPointPositionFactor)
  end select

end subroutine InitialQUAD_MP81

subroutine InitialQUAD_MP100(IParticle, PosGP, WeiGP)
  !**********************************************************************
  !
  !    SUBROUTINE: InitialQUAD_MP100
  !
  !>   Returns the initial local coordinates and weight for IParticle.
  !>   100 material points are placed in each element, corresponding to
  !>   10×10 Gauss quadrature points with direct weight values.
  !
  !>   @ note : 2D quadrilateral element
  !>   @ note : Weights are pre-calculated for 10×10 Gauss quadrature
  !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis.
  !    
  !>   @param[in] IParticle : Local number of the particle inside element (1~100)
  !
  !>   @param[inout] PosGP : Local coordinates (ξ, η) of the material point
  !>   @param[inout] WeiGP : Pre-calculated weight of the material point
  ! 
  !**********************************************************************
  implicit none

  integer(INTEGER_TYPE), intent(in) :: IParticle
  real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
  real(REAL_TYPE), intent(inout) :: WeiGP
  
  ! Local variables: 10-point Gauss coordinates (ξ/η ∈ [-1,1])
  real(REAL_TYPE) :: xi(10)
  real(REAL_TYPE) :: w(10)
  integer :: i, j, temp
  
  ! 10-point Gauss quadrature coordinates and weights
  xi(1) = -0.9739065285171717   
  xi(2) = -0.8650633666889845   
  xi(3) = -0.6794095682990244   
  xi(4) = -0.4333953941292472   
  xi(5) = -0.1488743389816312   
  xi(6) =  0.1488743389816312   
  xi(7) =  0.4333953941292472   
  xi(8) =  0.6794095682990244   
  xi(9) =  0.8650633666889845   
  xi(10) = 0.9739065285171717   
  
  ! 10-point Gauss quadrature weights
  w(1) = 0.0666713443086881     
  w(2) = 0.1494513491505806     
  w(3) = 0.2190863625159820     
  w(4) = 0.2692667193099963     
  w(5) = 0.2955242247147529     
  w(6) = 0.2955242247147529     
  w(7) = 0.2692667193099963     
  w(8) = 0.2190863625159820     
  w(9) = 0.1494513491505806     
  w(10) = 0.0666713443086881    

  ! Calculate i and j indices (1-10) from IParticle (1-100)
  select case (IParticle)
  case(1:100)
      j = (IParticle - 1) / 10 + 1  ! η direction index (row)
  
      ! Calculate i without using mod function
      temp = IParticle - (j - 1) * 10
      i = temp  ! ξ direction index (column)
  
      ! Assign coordinates based on indices
      PosGP(1) = xi(i)  ! ξ coordinate
      PosGP(2) = xi(j)  ! η coordinate
  
      ! Calculate weight as product of individual weights
      WeiGP = w(i) * w(j)
  
      ! Apply position shrinkage factor
      PosGP = PosGP * (1.0 - CalParams%shrinkageMateriaPointPositionFactor)
  end select

end subroutine InitialQUAD_MP100


        
        !**********************************************************************
        !
        !    SUBROUTINE: GaussQUAD_Q1
        !
        !    DESCRIPTION:
        !>   Returns the local coordinates and weight for IGaussPoint.
        !>   1 Gauss Point is located in the quadrilaters element.
        !
        !>   @ note : 2D element
        !>   @ note : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF
        !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987
        !
        !>   @param[in] IGaussPoint : Local number of the considered Gauss Point inside the element
        !
        !>   @param[inout] PosGP : Returns the xi, eta local coordinates of the Gauss Point with local number IGaussPoint
        !>   @param[inout] WeiGP : Returns the weight of IGaussPoint
        ! 
        !**********************************************************************
        subroutine GaussQUAD_Q1(IGaussPoint, PosGP, WeiGP)
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IGaussPoint
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          select case (IGaussPoint)
            case (1)
              PosGP(1) = 0.0
              PosGP(2) = 0.0
              WeiGP = 4.0
            case default
              call GiveError("Undefined number of Gauss points in [subroutine GaussQUAD_Q1()].")    
          end select

        end subroutine GaussQUAD_Q1


        !**********************************************************************
        !
        !    SUBROUTINE: GaussQUAD_Q4
        !
        !    DESCRIPTION:
        !>   Returns the local coordinates and weight for IGaussPoint.
        !>   4 Gauss Points are located in the quadrilateral element.
        !
        !>   @ note : 2D element
        !>   @ note : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF
        !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987
        !
        !>   @param[in] IGaussPoint : Local number of the considered Gauss Point inside the element
        !
        !>   @param[inout] PosGP : Returns the xi, eta local coordinates of the Gauss Point with local number IGaussPoint
        !>   @param[inout] WeiGP : Return the weight of IGaussPoint
        ! 
        !**********************************************************************
        subroutine GaussQUAD_Q4(IGaussPoint, PosGP, WeiGP)
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IGaussPoint
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          ! local variables
          real(REAL_TYPE) :: a = 0.57735026918962576450914878050196 ! = 1 / sqrt(3)
          real(REAL_TYPE) :: b = 1.0
          
          select case (IGaussPoint)
            case (1)
              PosGP(1) = -a
              PosGP(2) = -a
              WeiGP = b
            case (2)
              PosGP(1) = a
              PosGP(2) = -a
              WeiGP = b
            case (3)
              PosGP(1) = a
              PosGP(2) = a
              WeiGP = b
            case (4)
              PosGP(1) = -a
              PosGP(2) = a
              WeiGP = b
            case default
              call GiveError("Undefined number of Gauss points in [subroutine GaussQUAD_Q4()].")  
          end select
         
        end subroutine GaussQUAD_Q4
        
        !**********************************************************************
        !
        !    SUBROUTINE: ShapeLocPosQUAD4
        !
        !    DESCRIPTION:
        !>   To calculate the values of shape functions and their
        !>   derivatives at LocPos for a 4-noded 2D quadrilateral element.
        !
        !>   @note : 2D element
        !>   @note : https://ses.library.usyd.edu.au/bitstream/2123/709/8/adt-NU20060210.15574814appendixD.pdf
        !>   @note : R. K. Livesley, Finite Elements: An Introduction for Engineers, CUP Archive 1983
        !
        !>   @param[in] LocPos : Local coordinates of a point inside an element
        !
        !>   @param[out] HS(i) : Value of shape function i at LocPos
        !>   @param[out] dHS(i,j) : Value of derivative of shape function i at LocPos with respect to direction j
        !
        !                         ^ Eta
        !                 4       |
        !                +---------------+ 3
        !                |        |      |
        !                |        |      |
        !                |        |      |
        !                |        -------|---> Xi
        !                |               |
        !                |               |
        !                |1              | 2
        !                +---------------+-
        !
        !**********************************************************************
        
        subroutine InitialiseShapeFunctionsQUAD4(HS, dHS, Wt)
        real(REAL_TYPE), dimension(:, :), intent(inout) :: HS      ! Shape functions (Nint,4)
        real(REAL_TYPE), dimension(:, :, :), intent(inout) :: dHS  ! Derivatives (Nint,4,2)
        real(REAL_TYPE), dimension(:), intent(inout) :: Wt         ! Integration weights (Nint)

        ! --- 中心点积分参数 ---
        real(REAL_TYPE), parameter :: xi = 0.0, eta = 0.0  ! 自然坐标中心点
        real(REAL_TYPE), parameter :: w = 4.0              ! 总权重=4.0 (单元面积)

        ! --- 计算形函数值 ---
        HS(1,1) = (1 - xi) * (1 - eta) / 4.0  ! N1 = 0.25
        HS(1,2) = (1 + xi) * (1 - eta) / 4.0  ! N2 = 0.25
        HS(1,3) = (1 + xi) * (1 + eta) / 4.0  ! N3 = 0.25
        HS(1,4) = (1 - xi) * (1 + eta) / 4.0  ! N4 = 0.25

        ! --- 计算形函数导数（∂N/∂ξ）---
        dHS(1,1,1) = -(1 - eta) / 4.0  ! ∂N1/∂ξ = -0.25
        dHS(1,2,1) =  (1 - eta) / 4.0  ! ∂N2/∂ξ = 0.25
        dHS(1,3,1) =  (1 + eta) / 4.0  ! ∂N3/∂ξ = 0.25
        dHS(1,4,1) = -(1 + eta) / 4.0  ! ∂N4/∂ξ = -0.25

        ! --- 计算形函数导数（∂N/∂η）---
        dHS(1,1,2) = -(1 - xi) / 4.0   ! ∂N1/∂η = -0.25
        dHS(1,2,2) = -(1 + xi) / 4.0   ! ∂N2/∂η = -0.25
        dHS(1,3,2) =  (1 + xi) / 4.0   ! ∂N3/∂η = 0.25
        dHS(1,4,2) =  (1 - xi) / 4.0   ! ∂N4/∂η = 0.25

        ! --- 设置积分权重 ---
        Wt(1) = w  ! 权重=4.0
        
        end subroutine InitialiseShapeFunctionsQUAD4
        
        subroutine InitialiseShapeFunctionsQUAD4_4(HS, dHS, Wt)
        real(REAL_TYPE), dimension(:, :), intent(inout) :: HS      ! 形函数 (Nint,4)，Nint=4
        real(REAL_TYPE), dimension(:, :, :), intent(inout) :: dHS  ! 导数 (Nint,4,2)
        real(REAL_TYPE), dimension(:), intent(inout) :: Wt         ! 积分权重 (Nint)
        real(REAL_TYPE), parameter :: xi(4) = [ -1.0/sqrt(3.0),  1.0/sqrt(3.0),  1.0/sqrt(3.0), -1.0/sqrt(3.0) ]
        real(REAL_TYPE), parameter :: eta(4) = [ -1.0/sqrt(3.0), -1.0/sqrt(3.0),  1.0/sqrt(3.0),  1.0/sqrt(3.0) ]
        real(REAL_TYPE), parameter :: w = 1.0  ! 单个积分点权重
        integer :: i  ! 积分点索引

        do i = 1, 4
        ! 形函数值
            HS(i,1) = (1.0 - xi(i)) * (1.0 - eta(i)) / 4.0  ! N1
            HS(i,2) = (1.0 + xi(i)) * (1.0 - eta(i)) / 4.0  ! N2
            HS(i,3) = (1.0 + xi(i)) * (1.0 + eta(i)) / 4.0  ! N3
            HS(i,4) = (1.0 - xi(i)) * (1.0 + eta(i)) / 4.0  ! N4
            
            ! 形函数对ξ的导数
            dHS(i,1,1) = -(1.0 - eta(i)) / 4.0  ! ∂N1/∂ξ
            dHS(i,2,1) =  (1.0 - eta(i)) / 4.0  ! ∂N2/∂ξ
            dHS(i,3,1) =  (1.0 + eta(i)) / 4.0  ! ∂N3/∂ξ
            dHS(i,4,1) = -(1.0 + eta(i)) / 4.0  ! ∂N4/∂ξ
            
            ! 形函数对η的导数
            dHS(i,1,2) = -(1.0 - xi(i)) / 4.0   ! ∂N1/∂η
            dHS(i,2,2) = -(1.0 + xi(i)) / 4.0   ! ∂N2/∂η
            dHS(i,3,2) =  (1.0 + xi(i)) / 4.0   ! ∂N3/∂η
            dHS(i,4,2) =  (1.0 - xi(i)) / 4.0   ! ∂N4/∂η
            
            ! 积分权重
            Wt(i) = w
        
        end do
        
        end subroutine InitialiseShapeFunctionsQUAD4_4
        
        


        !**********************************************************************
        !
        !    SUBROUTINE: ShapeLocPosQUAD4
        !
        !    DESCRIPTION:
        !>   To calculate the values of shape functions and their
        !>   derivatives at LocPos for a 4-noded 2D quadrilateral element.
        !
        !>   @note : 2D element
        !>   @note : https://ses.library.usyd.edu.au/bitstream/2123/709/8/adt-NU20060210.15574814appendixD.pdf
        !>   @note : R. K. Livesley, Finite Elements: An Introduction for Engineers, CUP Archive 1983
        !
        !>   @param[in] LocPos : Local coordinates of a point inside an element
        !
        !>   @param[out] HS(i) : Value of shape function i at LocPos
        !>   @param[out] dHS(i,j) : Value of derivative of shape function i at LocPos with respect to direction j
        !
        !                         ^ Eta
        !                 4       |
        !                +---------------+ 3
        !                |        |      |
        !                |        |      |
        !                |        |      |
        !                |        -------|---> Xi
        !                |               |
        !                |               |
        !                |1              | 2
        !                +---------------+-
        !
        !**********************************************************************
        subroutine ShapeLocPosQUAD4(LocPos, HS, dHS)

        implicit none
      
          real(REAL_TYPE), dimension(:), intent(in) :: LocPos
          real(REAL_TYPE), dimension(:), intent(out) :: HS
          real(REAL_TYPE), dimension(:, :), intent(out) :: dHS

          ! local variables
          real(REAL_TYPE) :: Xi, Eta

          Xi = LocPos(1)
          Eta = LocPos(2)

          ! HS(i)
          HS(1) = (1.0 - Xi) * (1.0 - Eta) / 4.0
          HS(2) = (1.0 + Xi) * (1.0 - Eta) / 4.0
          HS(3) = (1.0 + Xi) * (1.0 + Eta) / 4.0
          HS(4) = (1.0 - Xi) * (1.0 + Eta) / 4.0

          ! dHS(i,1) = dHS / dXi
          dHS(1,1) =  - (1.0 - Eta) / 4.0
          dHS(2,1) =    (1.0 - Eta) / 4.0
          dHS(3,1) =    (1.0 + Eta) / 4.0
          dHS(4,1) =  - (1.0 + Eta) / 4.0

          ! dHS(i,2) = dHS / dEta
          dHS(1,2) =  - (1.0 - Xi) / 4.0
          dHS(2,2) =  - (1.0 + Xi) / 4.0
          dHS(3,2) =    (1.0 + Xi) / 4.0
          dHS(4,2) =    (1.0 - Xi) / 4.0

        end subroutine ShapeLocPosQUAD4


        !**********************************************************************
        !
        !    SUBROUTINE: ShapeLocPosQUAD8
        !
        !    DESCRIPTION:
        !>   To calculate the values of shape functions and their
        !>   derivatives at LocPos for a 8-noded 2D quadrilateral element.
        !
        !>   @note : 2D element
        !>   @note : https://ses.library.usyd.edu.au/bitstream/2123/709/8/adt-NU20060210.15574814appendixD.pdf
        !>   @note : R. K. Livesley, Finite Elements: An Introduction for Engineers, CUP Archive 1983, p.81
        !
        !>   @param[in] LocPos : Local coordinates of a point inside an element
        !
        !>   @param[out] HS(i) : Value of shape function i at LocPos
        !>   @param[out] dHS(i,j) : Value of derivative of shape function i at LocPos with respect to direction j
        !
        !                         ^ Eta
        !                 4       | 7
        !                +--------+-------+ 3
        !                |        |       |
        !                |        |       |
        !                | 8      |       | 6
        !                +        --------|---> Xi
        !                |                |
        !                |                |
        !                |1       5       | 2
        !                +--------+-------+
        !
        !**********************************************************************
        subroutine ShapeLocPosQUAD8(LocPos, HS, dHS)

        implicit none
      
          real(REAL_TYPE), dimension(:), intent(in) :: LocPos
          real(REAL_TYPE), dimension(:), intent(out) :: HS
          real(REAL_TYPE), dimension(:, :), intent(out) :: dHS

          ! local variables
          real(REAL_TYPE) :: Xi, Eta

          Xi = LocPos(1)
          Eta = LocPos(2)

          ! HS(i)
          ! 1..4 corner nodes, 5..8 middle nodes
          HS(1) = - (1.0 - Xi) * (1.0 - Eta) * (1.0 + Xi + Eta) / 4.0   
          HS(2) = - (1.0 + Xi) * (1.0 - Eta) * (1.0 - Xi + Eta) / 4.0 
          HS(3) = - (1.0 + Xi) * (1.0 + Eta) * (1.0 - Xi - Eta) / 4.0
          HS(4) = - (1.0 - Xi) * (1.0 + Eta) * (1.0 + Xi - Eta) / 4.0
          HS(5) =   (1.0 - Xi) * (1.0 + Xi) * (1.0 - Eta) / 2.0
          HS(6) =   (1.0 + Xi) * (1.0 + Eta) * (1.0 - Eta) / 2.0
          HS(7) =   (1.0 - Xi) * (1.0 + Xi) * (1.0 + Eta) / 2.0
          HS(8) =   (1.0 - Xi) * (1.0 + Eta) * (1.0 - Eta) / 2.0

          ! dHS(i,1) = dHS / dXi
          ! 1..4 corner nodess, 5..8 middle nodes
          dHS(1, 1) = - (-1.0) * (1.0 - Eta) * (1.0 + Xi + Eta) / 4.0  - (1.0 - Xi) * (1.0 - Eta) * (1.0) / 4.0
          dHS(2, 1) = - ( 1.0) * (1.0 - Eta) * (1.0 - Xi + Eta) / 4.0  - (1.0 + Xi) * (1.0 - Eta) * (-1.0) / 4.0
          dHS(3, 1) = - ( 1.0) * (1.0 + Eta) * (1.0 - Xi - Eta) / 4.0  - (1.0 + Xi) * (1.0 + Eta) * (-1.0) / 4.0
          dHS(4, 1) = - (-1.0) * (1.0 + Eta) * (1.0 + Xi - Eta) / 4.0  - (1.0 - Xi) * (1.0 + Eta) * (1.0) / 4.0
          dHS(5, 1) = (-1.0) * (1.0 + Xi) * (1.0 - Eta) / 2.0 + (1.0 - Xi) * (1.0) * (1.0 - Eta) / 2.0
          dHS(6, 1) = (1.0) * (1.0 + Eta) * (1.0 - Eta) / 2.0 
          dHS(7, 1) = (-1.0) * (1.0 + Xi) * (1.0 + Eta) / 2.0 + (1.0 - Xi) * (1.0) * (1.0 + Eta) / 2.0
          dHS(8, 1) = (-1.0) * (1.0 + Eta) * (1.0 - Eta) / 2.0 

          ! dHS(i,2) = dHS / dEta
          ! 1..4 corner nodes, 5..8 middle nodes
          dHS(1, 2) = - (1.0 - Xi) * (-1.0) * (1.0 + Xi + Eta) / 4.0 - (1.0 - Xi) * (1.0 - Eta) * (1.0) / 4.0 
          dHS(2, 2) = - (1.0 + Xi) * (-1.0) * (1.0 - Xi + Eta) / 4.0 - (1.0 + Xi) * (1.0 - Eta) * (1.0) / 4.0 
          dHS(3, 2) = - (1.0 + Xi) * (1.0) * (1.0 - Xi - Eta) / 4.0 - (1.0 + Xi) * (1.0 + Eta) * (-1.0) / 4.0
          dHS(4, 2) = - (1.0 - Xi) * (1.0) * (1.0 + Xi - Eta) / 4.0 - (1.0 - Xi) * (1.0 + Eta) * (-1.0) / 4.0
          dHS(5, 2) = (1.0 - Xi) * (1.0 + Xi) * (-1.0) / 2.0
          dHS(6, 2) = (1.0 + Xi) * (1.0) * (1.0 - Eta) / 2.0 + (1.0 + Xi) * (1.0 + Eta) * (- 1.0) / 2.0
          dHS(7, 2) = (1.0 - Xi) * (1.0 + Xi) * (1.0) / 2.0
          dHS(8, 2) = (1.0 - Xi) * (1.0) * (1.0 - Eta) / 2.0 + (1.0 - Xi) * (1.0 + Eta) * (-1.0) / 2.0

        end subroutine ShapeLocPosQUAD8


        !**********************************************************************
        !
        !    SUBROUTINE: GetLocalCoordinatesQUAD4
        !
        !    DESCRIPTION:
        !>   Determination of local coordinates LocPos from global coordinates GlobPos, 
        !>   assuming that the point lies inside the quadrilateral element.
        !>   OutsideElement returns .false. if the local position is in the element.
        !>   CrossedSide returns the number of the side that has been crossed.
        !
        !>   @param[in] GlobPos : Global coordinates of a point inside the element
        !>   @param[in] MInv : Element matrix
        !>   @param[in] MIX1 : Element vector of first element node
        !
        !>   @param[out] LocPos : Local coordinates of the considered point inside IElement
        !>   @param[out] OutsideElement : True, if the local coordinate lie outside IElement
        !>   @param[out] CrossedSide : ID of crossed side, if local coordinate lie outside IElement
        !
        !**********************************************************************
        subroutine GetLocalCoordinatesQUAD4(GlobPos, LocPos, OutsideElement, MInv, MIX1, CrossedSide)

        implicit none

          integer(INTEGER_TYPE), parameter :: IDim = 2 ! fixed dimension as 2D element
          
          real(REAL_TYPE), dimension(IDim), intent(in) :: GlobPos
          real(REAL_TYPE), dimension(IDim), intent(out):: LocPos
          logical, intent(out) :: OutsideElement
          real(REAL_TYPE), dimension(IDim, IDim), intent(in) :: MInv
          real(REAL_TYPE), dimension(IDim), intent(in) :: MIX1
          integer(INTEGER_TYPE), intent(out) :: CrossedSide
          
          ! local variables
          integer(INTEGER_TYPE) :: J
          
          OutsideElement = .true.

          LocPos = 0.0
          do J = 1, 2
            LocPos(1) = LocPos(1) + MInv(1, J) *  GlobPos(J)
            LocPos(2) = LocPos(2) + MInv(2, J) *  GlobPos(J)
          end do
         
          LocPos(1) = LocPos(1) - MIX1(1)
          LocPos(2) = LocPos(2) - MIX1(2)

          CrossedSide = -1
          if (LocPos(1) > 1.0) then
            CrossedSide = 2
          elseif (LocPos(1) < -1.0) then
            CrossedSide = 4
          elseif (LocPos(2) > 1.0) then
            CrossedSide = 3
          elseif (LocPos(2) < -1.0) then
            CrossedSide = 1
          else
            OutsideElement = .false.
          end if
        
        end subroutine GetLocalCoordinatesQUAD4


        !**********************************************************************
        !
        !    FUNCTION: IsInsideElementLocPosQUAD
        !  
        !    DESCRIPTION:
        !>   Returns .true. if LocPos (local coordinates) lies inside the 
        !>   area of the quadrilateral element.
        !
        !>   @param[in] LocPos : Local coordinates of the considered point inside the quadrilateral element
        !
        !>   @return IsInsideElementLocPosQUAD : True, if the point lies inside the quadrilateral element
        !
        !**********************************************************************
        logical function IsInsideElementLocPosQUAD(LocPos)
        
        implicit none
        
          real(REAL_TYPE), dimension(:), intent(in) :: LocPos
        
          if ( (LocPos(1) > 1.0) .or. (LocPos(1) < -1.0) .or. (LocPos(2) > 1.0) .or. (LocPos(2) < -1.0) ) then
            IsInsideElementLocPosQUAD = .false.
          else
            IsInsideElementLocPosQUAD = .true.
          end if
        
        end function IsInsideElementLocPosQUAD


        !**********************************************************************
        !
        !    FUNCTION: IsInsideElementGlobPosQUAD
        !
        !    DESCRIPTION:
        !>   Returns .true. if GlobPos (global coordinates) lies inside the 
        !>   area of the quadrilateral element.
        !
        !>   @param[in] GlobPos : Global coordinates of the considered point inside an element
        !>   @param[in] ElementID : ID of the considered element
        !>   @param[in] NodTot : Total number of nodes
        !>   @param[in] IElTyp : Number of nodes per element
        !>   @param[in] NEl : Total number of elements
        !>   @param[in] NodeCoord : Nodal coordinates
        !>   @param[in] ICon : Element connectivities
        !
        !>   @return IsInsideElementGlobPosQUAD : True, if the point lies inside the element
        !
        !**********************************************************************
        logical function IsInsideElementGlobPosQUAD(GlobPos, ElementID, NodTot, IElTyp, NEl, NodeCoord, ICon)
        
        implicit none
        
          integer(INTEGER_TYPE), parameter :: IDim = 2 ! fixed dimension as 2D element
          
          real(REAL_TYPE), dimension(IDim), intent(in) :: GlobPos
          integer(INTEGER_TYPE), intent(in) :: ElementID
          integer(INTEGER_TYPE), intent(in) :: NodTot, IElTyp, NEl
          real(REAL_TYPE), dimension(NodTot, IDim), intent(in) :: NodeCoord
          integer(INTEGER_TYPE), dimension(IElTyp, NEl), intent(in) :: ICon
          
          ! local variables
          real(REAL_TYPE), dimension(IDim) :: A, B, C, D ! vertices of quadrilateral
          
          A = NodeCoord(ICon(1, ElementID), 1:2)
          B = NodeCoord(ICon(2, ElementID), 1:2)
          C = NodeCoord(ICon(3, ElementID), 1:2)
          D = NodeCoord(ICon(4, ElementID), 1:2)

          IsInsideElementGlobPosQUAD = CheckInsideQuadrilateral(A, B, C, D, GlobPos) 

        end function IsInsideElementGlobPosQUAD
        
        subroutine DetermineAdjacentParticlesQUAD_MP1(NElementParticles, ParticleStatus)
        !**********************************************************************
        !
        !    Function:  Determines which particles of an element lie next to side ISide
        !               (triangular element with initially 1 material point).
        !               Note: ParticleStatus is not initialised to .false. in order
        !                     to allow for a more flexible usage!
        !
        ! I  NElementParticles : Initial number of particles per element
        ! O  ParticleStatus : Set to .true.
        !
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: NElementParticles
          logical, dimension(NElementParticles), intent(inout) :: ParticleStatus

          ParticleStatus(1) = .true.

        end subroutine DetermineAdjacentParticlesQUAD_MP1
        
    subroutine DetermineAdjacentParticlesQUAD_MP4(ISide, NElementParticles, ParticleStatus)
    !**********************************************************************
    !
    !    Function:  Determines which particles of a QUADRILATERAL element lie next to side ISide
    !               (quadrilateral element with 4 material points, one at each node).
    !               Note: ParticleStatus is not initialised to .false. in order
    !                     to allow for a more flexible usage!
    !
    ! I  ISide : Local number of the element side (1 to 4)
    ! I  NElementParticles : Number of particles per element (must be 4)
    ! O  ParticleStatus : Set to .true. for particles adjacent to ISide
    !
    !**********************************************************************
    implicit none

    integer(INTEGER_TYPE), intent(in) :: ISide, NElementParticles
    logical, dimension(NElementParticles), intent(inout) :: ParticleStatus
    ! Determine particles adjacent to the specified side (ISide)
    select case (ISide)
        case (1)  ! Side 1: Between node 1 and node 2 (particles 1 and 2)
            ParticleStatus(1) = .true.
            ParticleStatus(2) = .true.
        case (2)  ! Side 2: Between node 2 and node 3 (particles 2 and 3)
            ParticleStatus(2) = .true.
            ParticleStatus(3) = .true.
        case (3)  ! Side 3: Between node 3 and node 4 (particles 3 and 4)
            ParticleStatus(3) = .true.
            ParticleStatus(4) = .true.
        case (4)  ! Side 4: Between node 4 and node 1 (particles 4 and 1)
            ParticleStatus(4) = .true.
            ParticleStatus(1) = .true.
        end select
    
    end subroutine DetermineAdjacentParticlesQUAD_MP4
    
    subroutine DetermineAdjacentParticlesQUAD_MP9(ISide, NElementParticles, ParticleStatus)
    !**********************************************************************
    !
    !    Function:  Determines which particles of a QUADRILATERAL element lie next to side ISide
    !               (quadrilateral element with 9 material points: 4 corners, 4 edge centers, 1 center).
    !               Note: ParticleStatus is not initialised to .false. by default!
    !
    ! I  ISide : Local number of the element side (1 to 4)
    ! I  NElementParticles : Number of particles per element (must be 9)
    ! O  ParticleStatus : Set to .true. for particles adjacent to ISide
    !
    !     Particle numbering convention (example):
    !       4-----7-----3
    !       |           |
    !       8           6
    !       |           |
    !       1-----5-----2
    !
    !**********************************************************************
    implicit none

    integer(INTEGER_TYPE), intent(in) :: ISide, NElementParticles
    logical, dimension(NElementParticles), intent(inout) :: ParticleStatus

    ! Determine particles adjacent to the specified side (ISide)
    select case (ISide)
        case (1)  ! Side 1: Between node 1 and node 2 (particles 1, 5, 2)
            ParticleStatus(1) = .true.  ! Corner 1
            ParticleStatus(5) = .true.  ! Edge center
            ParticleStatus(2) = .true.  ! Corner 2
        case (2)  ! Side 2: Between node 2 and node 3 (particles 2, 6, 3)
            ParticleStatus(2) = .true.
            ParticleStatus(6) = .true.
            ParticleStatus(3) = .true.
        case (3)  ! Side 3: Between node 3 and node 4 (particles 3, 7, 4)
            ParticleStatus(3) = .true.
            ParticleStatus(7) = .true.
            ParticleStatus(4) = .true.
        case (4)  ! Side 4: Between node 4 and node 1 (particles 4, 8, 1)
            ParticleStatus(4) = .true.
            ParticleStatus(8) = .true.
            ParticleStatus(1) = .true.
    end select

    end subroutine DetermineAdjacentParticlesQUAD_MP9
    
    subroutine DetermineAdjacentParticlesQUAD_MP16(ISide, NElementParticles, ParticleStatus)
    !**********************************************************************
    !
    !    Function:  Determines which particles of a QUADRILATERAL element lie next to side ISide
    !               (quadrilateral element with 16 material points in 4×4 Gauss layout).
    !               Note: ParticleStatus is not initialised to .false. in order
    !                     to allow for a more flexible usage!
    !
    ! I  ISide : Local number of the element side (1 to 4)
    ! I  NElementParticles : Number of particles per element (must be 16)
    ! O  ParticleStatus : Set to .true. for particles adjacent to ISide
    !
    !    Side Definitions:
    !    - Side 1 (Bottom): η = -1 (Particles 1, 2, 3, 4)
    !    - Side 2 (Right):  ξ = +1 (Particles 4, 8, 12, 16)
    !    - Side 3 (Top):    η = +1 (Particles 13, 14, 15, 16)
    !    - Side 4 (Left):   ξ = -1 (Particles 1, 5, 9, 13)
    !
    !**********************************************************************
    implicit none

    integer(INTEGER_TYPE), intent(in) :: ISide, NElementParticles
    logical, dimension(NElementParticles), intent(inout) :: ParticleStatus

    ! Determine particles adjacent to the specified side (ISide)
    select case (ISide)
        case (1)  ! Side 1: Bottom edge (Particles 1, 2, 3, 4)
            ParticleStatus(1) = .true.
            ParticleStatus(2) = .true.
            ParticleStatus(3) = .true.
            ParticleStatus(4) = .true.
        case (2)  ! Side 2: Right edge (Particles 4, 8, 12, 16)
            ParticleStatus(4) = .true.
            ParticleStatus(8) = .true.
            ParticleStatus(12) = .true.
            ParticleStatus(16) = .true.
        case (3)  ! Side 3: Top edge (Particles 13, 14, 15, 16)
            ParticleStatus(13) = .true.
            ParticleStatus(14) = .true.
            ParticleStatus(15) = .true.
            ParticleStatus(16) = .true.
        case (4)  ! Side 4: Left edge (Particles 1, 5, 9, 13)
            ParticleStatus(1) = .true.
            ParticleStatus(5) = .true.
            ParticleStatus(9) = .true.
            ParticleStatus(13) = .true.
        case default
            call GiveError("ISide must be 1~4 in [subroutine DetermineAdjacentParticlesQUAD_MP16()].")
    end select

    end subroutine DetermineAdjacentParticlesQUAD_MP16
    
    subroutine DetermineAdjacentParticlesQUAD_MP25(ISide, NElementParticles, ParticleStatus)
    !**********************************************************************
    !
    !    Function:  Determines which particles of a QUADRILATERAL element lie next to side ISide
    !               (quadrilateral element with 25 material points in 5×5 Gauss layout).
    !               Note: ParticleStatus is not initialised to .false. in order
    !                     to allow for a more flexible usage!
    !
    ! I  ISide : Local number of the element side (1 to 4)
    ! I  NElementParticles : Number of particles per element (must be 25)
    ! O  ParticleStatus : Set to .true. for particles adjacent to ISide
    !
    !    Side Definitions:
    !    - Side 1 (Bottom): η = -1 (Particles 1, 2, 3, 4, 5)
    !    - Side 2 (Right):  ξ = +1 (Particles 5, 10, 15, 20, 25)
    !    - Side 3 (Top):    η = +1 (Particles 21, 22, 23, 24, 25)
    !    - Side 4 (Left):   ξ = -1 (Particles 1, 6, 11, 16, 21)
    !
    !**********************************************************************
    implicit none

    integer(INTEGER_TYPE), intent(in) :: ISide, NElementParticles
    logical, dimension(NElementParticles), intent(inout) :: ParticleStatus

    ! Determine particles adjacent to the specified side (ISide)
    select case (ISide)
        case (1)  ! Side 1: Bottom edge (Particles 1-5, 1st row of 5×5)
            ParticleStatus(1) = .true.
            ParticleStatus(2) = .true.
            ParticleStatus(3) = .true.
            ParticleStatus(4) = .true.
            ParticleStatus(5) = .true.
        case (2)  ! Side 2: Right edge (Particles 5,10,15,20,25, last column of 5×5)
            ParticleStatus(5) = .true.
            ParticleStatus(10) = .true.
            ParticleStatus(15) = .true.
            ParticleStatus(20) = .true.
            ParticleStatus(25) = .true.
        case (3)  ! Side 3: Top edge (Particles 21-25, last row of 5×5)
            ParticleStatus(21) = .true.
            ParticleStatus(22) = .true.
            ParticleStatus(23) = .true.
            ParticleStatus(24) = .true.
            ParticleStatus(25) = .true.
        case (4)  ! Side 4: Left edge (Particles 1,6,11,16,21, first column of 5×5)
            ParticleStatus(1) = .true.
            ParticleStatus(6) = .true.
            ParticleStatus(11) = .true.
            ParticleStatus(16) = .true.
            ParticleStatus(21) = .true.
        case default
            call GiveError("ISide must be 1~4 in [subroutine DetermineAdjacentParticlesQUAD_MP25()].")
    end select

    end subroutine DetermineAdjacentParticlesQUAD_MP25
    
    subroutine DetermineAdjacentParticlesQUAD_MP36(ISide, NElementParticles, ParticleStatus)
    !**********************************************************************
    !
    !    Function:  Determines which particles of a QUADRILATERAL element lie next to side ISide
    !               (quadrilateral element with 36 material points in 6×6 Gauss layout).
    !               Note: ParticleStatus is not initialised to .false. in order
    !                     to allow for a more flexible usage!
    !
    ! I  ISide : Local number of the element side (1 to 4)
    ! I  NElementParticles : Number of particles per element (must be 36)
    ! O  ParticleStatus : Set to .true. for particles adjacent to ISide
    !
    !    Side Definitions:
    !    - Side 1 (Bottom): η = -1 (Particles 1, 2, 3, 4, 5, 6)
    !    - Side 2 (Right):  ξ = +1 (Particles 6, 12, 18, 24, 30, 36)
    !    - Side 3 (Top):    η = +1 (Particles 31, 32, 33, 34, 35, 36)
    !    - Side 4 (Left):   ξ = -1 (Particles 1, 7, 13, 19, 25, 31)
    !
    !**********************************************************************
    implicit none

    integer(INTEGER_TYPE), intent(in) :: ISide, NElementParticles
    logical, dimension(NElementParticles), intent(inout) :: ParticleStatus

    ! Determine particles adjacent to the specified side (ISide)
    select case (ISide)
        case (1)  ! Side 1: Bottom edge (Particles 1-6, 1st row of 6×6)
            ParticleStatus(1) = .true.
            ParticleStatus(2) = .true.
            ParticleStatus(3) = .true.
            ParticleStatus(4) = .true.
            ParticleStatus(5) = .true.
            ParticleStatus(6) = .true.
        case (2)  ! Side 2: Right edge (Particles 6,12,18,24,30,36, last column of 6×6)
            ParticleStatus(6) = .true.
            ParticleStatus(12) = .true.
            ParticleStatus(18) = .true.
            ParticleStatus(24) = .true.
            ParticleStatus(30) = .true.
            ParticleStatus(36) = .true.
        case (3)  ! Side 3: Top edge (Particles 31-36, last row of 6×6)
            ParticleStatus(31) = .true.
            ParticleStatus(32) = .true.
            ParticleStatus(33) = .true.
            ParticleStatus(34) = .true.
            ParticleStatus(35) = .true.
            ParticleStatus(36) = .true.
        case (4)  ! Side 4: Left edge (Particles 1,7,13,19,25,31, first column of 6×6)
            ParticleStatus(1) = .true.
            ParticleStatus(7) = .true.
            ParticleStatus(13) = .true.
            ParticleStatus(19) = .true.
            ParticleStatus(25) = .true.
            ParticleStatus(31) = .true.
        case default
            call GiveError("ISide must be 1~4 in [subroutine DetermineAdjacentParticlesQUAD_MP36()].")
    end select

    end subroutine DetermineAdjacentParticlesQUAD_MP36
    
    subroutine DetermineAdjacentParticlesQUAD_MP49(ISide, NElementParticles, ParticleStatus)
    !**********************************************************************
    !
    !    Function:  Determines which particles of a QUADRILATERAL element lie next to side ISide
    !               (quadrilateral element with 49 material points in 7×7 Gauss layout).
    !               Note: ParticleStatus is not initialised to .false. in order
    !                     to allow for a more flexible usage!
    !
    ! I  ISide : Local number of the element side (1 to 4)
    ! I  NElementParticles : Number of particles per element (must be 49)
    ! O  ParticleStatus : Set to .true. for particles adjacent to ISide
    !
    !    Side Definitions:
    !    - Side 1 (Bottom): η = -1 (Particles 1, 2, 3, 4, 5, 6, 7)
    !    - Side 2 (Right):  ξ = +1 (Particles 7, 14, 21, 28, 35, 42, 49)
    !    - Side 3 (Top):    η = +1 (Particles 43, 44, 45, 46, 47, 48, 49)
    !    - Side 4 (Left):   ξ = -1 (Particles 1, 8, 15, 22, 29, 36, 43)
    !
    !**********************************************************************
    implicit none

    integer(INTEGER_TYPE), intent(in) :: ISide, NElementParticles
    logical, dimension(NElementParticles), intent(inout) :: ParticleStatus

    ! Determine particles adjacent to the specified side (ISide)
    select case (ISide)
        case (1)  ! Side 1: Bottom edge (Particles 1-7, 1st row of 7×7)
            ParticleStatus(1) = .true.
            ParticleStatus(2) = .true.
            ParticleStatus(3) = .true.
            ParticleStatus(4) = .true.
            ParticleStatus(5) = .true.
            ParticleStatus(6) = .true.
            ParticleStatus(7) = .true.
        case (2)  ! Side 2: Right edge (Particles 7,14,21,28,35,42,49, last column of 7×7)
            ParticleStatus(7) = .true.
            ParticleStatus(14) = .true.
            ParticleStatus(21) = .true.
            ParticleStatus(28) = .true.
            ParticleStatus(35) = .true.
            ParticleStatus(42) = .true.
            ParticleStatus(49) = .true.
        case (3)  ! Side 3: Top edge (Particles 43-49, last row of 7×7)
            ParticleStatus(43) = .true.
            ParticleStatus(44) = .true.
            ParticleStatus(45) = .true.
            ParticleStatus(46) = .true.
            ParticleStatus(47) = .true.
            ParticleStatus(48) = .true.
            ParticleStatus(49) = .true.
        case (4)  ! Side 4: Left edge (Particles 1,8,15,22,29,36,43, first column of 7×7)
            ParticleStatus(1) = .true.
            ParticleStatus(8) = .true.
            ParticleStatus(15) = .true.
            ParticleStatus(22) = .true.
            ParticleStatus(29) = .true.
            ParticleStatus(36) = .true.
            ParticleStatus(43) = .true.
        case default
            call GiveError("ISide must be 1~4 in [subroutine DetermineAdjacentParticlesQUAD_MP49()].")
    end select

end subroutine DetermineAdjacentParticlesQUAD_MP49
    
    subroutine DetermineAdjacentParticlesQUAD_MP64(ISide, NElementParticles, ParticleStatus)
    !**********************************************************************
    !
    !    Function:  Determines which particles of a QUADRILATERAL element lie next to side ISide
    !               (quadrilateral element with 64 material points in 8×8 Gauss layout).
    !               Note: ParticleStatus is not initialised to .false. in order
    !                     to allow for a more flexible usage!
    !
    ! I  ISide : Local number of the element side (1 to 4)
    ! I  NElementParticles : Number of particles per element (must be 64)
    ! O  ParticleStatus : Set to .true. for particles adjacent to ISide
    !
    !    Side Definitions:
    !    - Side 1 (Bottom): η = -1 (Particles 1, 2, 3, 4, 5, 6, 7, 8)
    !    - Side 2 (Right):  ξ = +1 (Particles 8, 16, 24, 32, 40, 48, 56, 64)
    !    - Side 3 (Top):    η = +1 (Particles 57, 58, 59, 60, 61, 62, 63, 64)
    !    - Side 4 (Left):   ξ = -1 (Particles 1, 9, 17, 25, 33, 41, 49, 57)
    !
    !**********************************************************************
    implicit none

    integer(INTEGER_TYPE), intent(in) :: ISide, NElementParticles
    logical, dimension(NElementParticles), intent(inout) :: ParticleStatus

    ! Determine particles adjacent to the specified side (ISide)
    select case (ISide)
        case (1)  ! Side 1: Bottom edge (Particles 1, 2, 3, 4, 5, 6, 7, 8)
            ParticleStatus(1) = .true.
            ParticleStatus(2) = .true.
            ParticleStatus(3) = .true.
            ParticleStatus(4) = .true.
            ParticleStatus(5) = .true.
            ParticleStatus(6) = .true.
            ParticleStatus(7) = .true.
            ParticleStatus(8) = .true.
        case (2)  ! Side 2: Right edge (Particles 8, 16, 24, 32, 40, 48, 56, 64)
            ParticleStatus(8) = .true.
            ParticleStatus(16) = .true.
            ParticleStatus(24) = .true.
            ParticleStatus(32) = .true.
            ParticleStatus(40) = .true.
            ParticleStatus(48) = .true.
            ParticleStatus(56) = .true.
            ParticleStatus(64) = .true.
        case (3)  ! Side 3: Top edge (Particles 57, 58, 59, 60, 61, 62, 63, 64)
            ParticleStatus(57) = .true.
            ParticleStatus(58) = .true.
            ParticleStatus(59) = .true.
            ParticleStatus(60) = .true.
            ParticleStatus(61) = .true.
            ParticleStatus(62) = .true.
            ParticleStatus(63) = .true.
            ParticleStatus(64) = .true.
        case (4)  ! Side 4: Left edge (Particles 1, 9, 17, 25, 33, 41, 49, 57)
            ParticleStatus(1) = .true.
            ParticleStatus(9) = .true.
            ParticleStatus(17) = .true.
            ParticleStatus(25) = .true.
            ParticleStatus(33) = .true.
            ParticleStatus(41) = .true.
            ParticleStatus(49) = .true.
            ParticleStatus(57) = .true.
        case default
            call GiveError("ISide must be 1~4 in [subroutine DetermineAdjacentParticlesQUAD_MP64()].")
    end select

    end subroutine DetermineAdjacentParticlesQUAD_MP64
    
    subroutine DetermineAdjacentParticlesQUAD_MP81(ISide, NElementParticles, ParticleStatus)
    !**********************************************************************
    !
    !    Function:  Determines which particles of a QUADRILATERAL element lie next to side ISide
    !               (quadrilateral element with 81 material points in 9×9 Gauss layout).
    !               Note: ParticleStatus is not initialised to .false. in order
    !                     to allow for a more flexible usage!
    !
    ! I  ISide : Local number of the element side (1 to 4)
    ! I  NElementParticles : Number of particles per element (must be 81)
    ! O  ParticleStatus : Set to .true. for particles adjacent to ISide
    !
    !    Side Definitions:
    !    - Side 1 (Bottom): η = -1 (Particles 1, 2, 3, 4, 5, 6, 7, 8, 9)
    !    - Side 2 (Right):  ξ = +1 (Particles 9, 18, 27, 36, 45, 54, 63, 72, 81)
    !    - Side 3 (Top):    η = +1 (Particles 73, 74, 75, 76, 77, 78, 79, 80, 81)
    !    - Side 4 (Left):   ξ = -1 (Particles 1, 10, 19, 28, 37, 46, 55, 64, 73)
    !
    !**********************************************************************
    implicit none

    integer(INTEGER_TYPE), intent(in) :: ISide, NElementParticles
    logical, dimension(NElementParticles), intent(inout) :: ParticleStatus

    ! Determine particles adjacent to the specified side (ISide)
    select case (ISide)
        case (1)  ! Side 1: Bottom edge (Particles 1-9, 1st row of 9×9)
            ParticleStatus(1) = .true.
            ParticleStatus(2) = .true.
            ParticleStatus(3) = .true.
            ParticleStatus(4) = .true.
            ParticleStatus(5) = .true.
            ParticleStatus(6) = .true.
            ParticleStatus(7) = .true.
            ParticleStatus(8) = .true.
            ParticleStatus(9) = .true.
        case (2)  ! Side 2: Right edge (Particles 9,18,27,36,45,54,63,72,81, last column of 9×9)
            ParticleStatus(9) = .true.
            ParticleStatus(18) = .true.
            ParticleStatus(27) = .true.
            ParticleStatus(36) = .true.
            ParticleStatus(45) = .true.
            ParticleStatus(54) = .true.
            ParticleStatus(63) = .true.
            ParticleStatus(72) = .true.
            ParticleStatus(81) = .true.
        case (3)  ! Side 3: Top edge (Particles 73-81, last row of 9×9)
            ParticleStatus(73) = .true.
            ParticleStatus(74) = .true.
            ParticleStatus(75) = .true.
            ParticleStatus(76) = .true.
            ParticleStatus(77) = .true.
            ParticleStatus(78) = .true.
            ParticleStatus(79) = .true.
            ParticleStatus(80) = .true.
            ParticleStatus(81) = .true.
        case (4)  ! Side 4: Left edge (Particles 1,10,19,28,37,46,55,64,73, first column of 9×9)
            ParticleStatus(1) = .true.
            ParticleStatus(10) = .true.
            ParticleStatus(19) = .true.
            ParticleStatus(28) = .true.
            ParticleStatus(37) = .true.
            ParticleStatus(46) = .true.
            ParticleStatus(55) = .true.
            ParticleStatus(64) = .true.
            ParticleStatus(73) = .true.
        case default
            call GiveError("ISide must be 1~4 in [subroutine DetermineAdjacentParticlesQUAD_MP81()].")
    end select

end subroutine DetermineAdjacentParticlesQUAD_MP81

subroutine DetermineAdjacentParticlesQUAD_MP100(ISide, NElementParticles, ParticleStatus)
    !**********************************************************************
    !
    !    Function:  Determines which particles of a QUADRILATERAL element lie next to side ISide
    !               (quadrilateral element with 100 material points in 10×10 Gauss layout).
    !               Note: ParticleStatus is not initialised to .false. in order
    !                     to allow for a more flexible usage!
    !
    ! I  ISide : Local number of the element side (1 to 4)
    ! I  NElementParticles : Number of particles per element (must be 100)
    ! O  ParticleStatus : Set to .true. for particles adjacent to ISide
    !
    !    Side Definitions:
    !    - Side 1 (Bottom): η = -1 (Particles 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
    !    - Side 2 (Right):  ξ = +1 (Particles 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
    !    - Side 3 (Top):    η = +1 (Particles 91, 92, 93, 94, 95, 96, 97, 98, 99, 100)
    !    - Side 4 (Left):   ξ = -1 (Particles 1, 11, 21, 31, 41, 51, 61, 71, 81, 91)
    !
    !**********************************************************************
    implicit none

    integer(INTEGER_TYPE), intent(in) :: ISide, NElementParticles
    logical, dimension(NElementParticles), intent(inout) :: ParticleStatus

    ! Determine particles adjacent to the specified side (ISide)
    select case (ISide)
        case (1)  ! Side 1: Bottom edge (Particles 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
            ParticleStatus(1) = .true.
            ParticleStatus(2) = .true.
            ParticleStatus(3) = .true.
            ParticleStatus(4) = .true.
            ParticleStatus(5) = .true.
            ParticleStatus(6) = .true.
            ParticleStatus(7) = .true.
            ParticleStatus(8) = .true.
            ParticleStatus(9) = .true.
            ParticleStatus(10) = .true.
        case (2)  ! Side 2: Right edge (Particles 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
            ParticleStatus(10) = .true.
            ParticleStatus(20) = .true.
            ParticleStatus(30) = .true.
            ParticleStatus(40) = .true.
            ParticleStatus(50) = .true.
            ParticleStatus(60) = .true.
            ParticleStatus(70) = .true.
            ParticleStatus(80) = .true.
            ParticleStatus(90) = .true.
            ParticleStatus(100) = .true.
        case (3)  ! Side 3: Top edge (Particles 91, 92, 93, 94, 95, 96, 97, 98, 99, 100)
            ParticleStatus(91) = .true.
            ParticleStatus(92) = .true.
            ParticleStatus(93) = .true.
            ParticleStatus(94) = .true.
            ParticleStatus(95) = .true.
            ParticleStatus(96) = .true.
            ParticleStatus(97) = .true.
            ParticleStatus(98) = .true.
            ParticleStatus(99) = .true.
            ParticleStatus(100) = .true.
        case (4)  ! Side 4: Left edge (Particles 1, 11, 21, 31, 41, 51, 61, 71, 81, 91)
            ParticleStatus(1) = .true.
            ParticleStatus(11) = .true.
            ParticleStatus(21) = .true.
            ParticleStatus(31) = .true.
            ParticleStatus(41) = .true.
            ParticleStatus(51) = .true.
            ParticleStatus(61) = .true.
            ParticleStatus(71) = .true.
            ParticleStatus(81) = .true.
            ParticleStatus(91) = .true.
        case default
            call GiveError("ISide must be 1~4 in [subroutine DetermineAdjacentParticlesQUAD_MP100()].")
    end select

end subroutine DetermineAdjacentParticlesQUAD_MP100


        !**********************************************************************
        !
        !    SUBROUTINE: DetermineCheckEdgeNodesQUAD8
        !
        !    DESCRIPTION:
        !>   Determines which edge node of each side to check in order
        !>   to detect an adjacent element for each side.
        !
        !>   @param[inout] CheckEdgeNodes : Array containing for each side of the element, the local number of the edge node
        !
        !**********************************************************************
        subroutine DetermineCheckEdgeNodesQUAD8(CheckEdgeNodes)
                
        implicit none
        
          integer(INTEGER_TYPE), dimension(:, :), intent(inout) :: CheckEdgeNodes ! size(nCheckNodes, NumberOfElementSides)

          CheckEdgeNodes(1,1) = 5 ! Edge node to check for side 1 spanned by nodes 1-2
          CheckEdgeNodes(1,2) = 6 ! Edge node to check for side 2 spanned by nodes 2-3
          CheckEdgeNodes(1,3) = 7 ! Edge node to check for side 3 spanned by nodes 3-4
          CheckEdgeNodes(1,4) = 8 ! Edge node to check for side 3 spanned by nodes 4-1

        end subroutine DetermineCheckEdgeNodesQUAD8
        
      
        !**********************************************************************
        !
        !    FUNCTION: DetermineSideNodesQUAD8
        !
        !    DESCRIPTION:
        !>   Returns the element connectivity of LocalNodeID of side SideID.
        !>     Side 1 is spanned by nodes 1, 2, 5
        !>     Side 2 is spanned by nodes 2, 3, 6
        !>     Side 3 is spanned by nodes 3, 4, 7
        !>     Side 4 is spanned by nodes 4, 1, 8
        !
        !>   @param[in] SideID : ID of the considered side (1 .. 4)
        !>   @param[in] LocalNodeID : ID of the side node  (1 .. 8)
        !
        !>   @return DetermineSideNodesQUAD8 : Local ID of the considered node (1 .. 8)
        !
        !**********************************************************************
        integer function DetermineSideNodesQUAD8(SideID, LocalNodeID)
        
        implicit none

        integer(INTEGER_TYPE), intent(in) :: SideID, LocalNodeID
        
        ! local variables
        integer(INTEGER_TYPE), dimension(3, 4) :: SideConnectivities

        SideConnectivities = reshape( (/  1, 2, 5, &
                                          2, 3, 6, &
                                          3, 4, 7, &
                                          4, 1, 8/), &
                                       (/ 3, 4 /) )

        DetermineSideNodesQUAD8 = SideConnectivities(LocalNodeID, SideID)

      end function DetermineSideNodesQUAD8

      
      subroutine CheckQUADForGlobPos(GlobPos, ElementID, CentrePoint, NodeCoord, ICon, CrossedSide, IsInside)
      !**********************************************************************
      !
      !    SUBROUTINE: CheckQUADForGlobPos
      !
      !    DESCRIPTION:
      !>   Determines whether GlobPos lies inside the element with ElementID
      !>   (result written to IsInside) and, which side of the quadrilateral is
      !>   crossed by the line between the centrepoint of ElementID and GlobPos
      !>   if GlobPos lies in another element (CrossedSide).
      !>     Side 1 (nodes 1-2 & 5) at xi line
      !>     Side 2 (nodes 2-3 & 6) at eta line
      !>     Side 3 (nodes 3-4 & 7) at xi line
      !>     Side 4 (nodes 4-1 & 8) at eta line
      !
      !>   @param[in] GlobPos : Global coordinates of a point inside the mesh
      !>   @param[in] ElementID : Considered element
      !>   @param[in] CentrePoint : Centrepoint of ElementID
      !>   @param[in] NodTot : Total number of nodes
      !>   @param[in] IElTyp : Number of node connectivities of IElement
      !>   @param[in] NEl : Number of elements
      !>   @param[in] NodeCoord : Global nodal coordinates
      !>   @param[in] ICon : Element connectivities ICon(I, J): global node number of local node I in element J
      !
      !>   @param[out] CrossedSide : Side which contains the intersection point of the above mentioned line
      !>   @param[out] IsInside : True, if GlobPos lies inside ElementID
      !
      !**********************************************************************

      implicit none
      
        integer(INTEGER_TYPE), parameter :: IDim = 2 ! fixed dimension as only 2D element
        
        real(REAL_TYPE), dimension(:), intent(in) :: GlobPos
        integer(INTEGER_TYPE), intent(in) :: ElementID
        real(REAL_TYPE), dimension(:), intent(in) :: CentrePoint
        real(REAL_TYPE), dimension(:, :), intent(in) :: NodeCoord
        integer(INTEGER_TYPE), dimension(:, :), intent(in) :: ICon
        integer(INTEGER_TYPE), intent(out) :: CrossedSide
        logical, intent(out) :: IsInside
        
        ! local variables
        integer(INTEGER_TYPE), dimension(2, 4) :: VerticesID ! size (Number of node on each side, Number of sides of the quadrilateral)
        real(REAL_TYPE), dimension(IDim) :: A, B ! coordinates of nodes of one side of a quadrilateral
        integer(INTEGER_TYPE) :: Success, I
        real(REAL_TYPE) :: distanceP_CentrePoint

        VerticesID = reshape( (/ 1, 2,   & ! Corner nodes to check for side 1 spanned by nodes 1-2
                                 2, 3,   & ! Corner nodes to check for side 2 spanned by nodes 2-3
                                 3, 4,   & ! Corner nodes to check for side 3 spanned by nodes 3-4
                                 4, 1/), & ! Corner nodes to check for side 4 spanned by nodes 4-1
                              (/ 2, 4 /) ) ! 8 elements which should be reshaped to a 2x4 matrix. N_SIDE is always 4 for a quadrilateral.

        do I = 1, 4 ! loop over sides of quadrilateral
            
          A = NodeCoord(ICon(VerticesID(1, I), ElementID), 1:NDIM)
          B = NodeCoord(ICon(VerticesID(2, I), ElementID), 1:NDIM)

          ! compute distance P-CentrePoint
          distanceP_CentrePoint = Distance(GlobPos, CentrePoint, NDIM)
          
          if(distanceP_CentrePoint > SMALL) then
              Success = CheckInsideSubTriangle(A, B, CentrePoint, GlobPos)
          else
              Success = 2 ! point is exactly on top of the center node 
          endif
          if (Success == 1) then
            CrossedSide = I
            EXIT
          else if (Success == 2) then
            IsInside = .true.
            EXIT
          end if

        end do

      end subroutine CheckQUADForGlobPos
                                               
    end module ModElementEvaluationQUAD
    
