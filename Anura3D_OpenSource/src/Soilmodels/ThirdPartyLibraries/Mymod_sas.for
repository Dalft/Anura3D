! Copyright (C)  2007  P.-A. von Wolffersdorff, D. Masin
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
!  USA.

C     Last change:  VW   08 April 2009    17:18 pm
      Subroutine MyMod_sas ( IDTask, iMod, IsUndr,
     *                      iStep, iTer, iEl, Int,
     *                      X, Y, Z,
     *                      Time0, dTime,
     *                      Props, Sig0, Swp0, StVar0,
     *                      dEps, D, BulkW,
     *                      Sig, Swp, StVar, ipl,
     *                      nStat,
     *                      NonSym, iStrsDep, iTimeDep, iTang,
     *                      iAbort )
!
! Purpose: User supplied soil model 
!  Depending on IDTask, 1 : Initialize state variables
!                       2 : calculate stresses,
!                       3 : calculate material stiffness matrix
!                       4 : return number of state variables
!                       5 : inquire matrix properties
!                       6 : calculate elastic material stiffness matrix
!                           return switch for non-symmetric D-matrix
!                           stress/time dependent matrix
! Arguments:
!          I/O  Type
!  IDTask   I   I    : see above
!  iMod     I   I    : model number (1..10)
!  IsUndr   I   I    : =1 for undrained, 0 otherwise
!  iStep    I   I    : Global step number
!  iter     I   I    : Global iteration number
!  iel      I   I    : Global element number
!  Int      I   I    : Global integration point number
!  X        I   R    : X-Position of integration point
!  Y        I   R    : Y-Position of integration point
!  Z        I   R    : Z-Position of integration point
!  Time0    I   R    : Time at start of step
!  dTime    I   R    : Time increment
!  Props    I   R()  : List with model parameters
!  Sig0     I   R()  : Stresses at start of step
!  Swp0     I   R    : Excess pore pressure start of step
!  StVar0   I   R()  : State variable at start of step
!  dEps     I   R()  : Strain increment
!  D       I/O  R(,) : Material stiffness matrix
!  BulkW   I/O  R    : Bulkmodulus for water (undrained only)
!  Sig      O   R()  : Resulting stresses
!  Swp      O   R    : Resulting excess pore pressure
!  StVar    O   R()  : Resulting values state variables
!  ipl      O   I    : Plasticity indicator
!  nStat    O   I    : Number of state variables
!  NonSym   O   I    : Non-Symmetric D-matrix ?
!  iStrsDep O   I    : =1 for stress dependent D-matrix
!  iTimeDep O   I    : =1 for time dependent D-matrix
!  iAbort   O   I    : =1 to force stopping of calculation

      Implicit Double Precision (A-H, O-Z)

      Dimension Props(*), Sig0(*), StVar0(*), dEps(*), D(6,6),
     *          Sig(*),   StVar(*)



!---  Local variables

      Character*8 cmname
      Dimension stress(6), statev(36), ddsdde(6,6), dstran(6),
     *          drot(3,3), gstran(6)

! added dimension
      Dimension ddsddt(6), drplde(6), stran(6), time(2),
     *          predef(2), dpred(1), dfgrd0(3,3), dfgrd1(3,3),
     *          coords(3)

      do i=1,6
            stran(i) = 0.0d0  
      enddo!    
      

      nStatV = SIZE(statev)
      
c     Cases of UMAT call

      If (IDTask .Eq. 1 .OR.
     *    IDTask .Eq. 2 .OR.
     *    IDTask .Eq. 3 .OR.
     *    IDtask .Eq. 6 )    Then

         nprops = 50          ! Dimension of array Props (1..50)
         ntens  = 6           ! Dimension of vectors
         ndi    = ntens - 3   ! Number of mean components
         noel   = iel         ! Global element number
         npt    = Int         ! Global integration point number

         kstep  = iStep       ! Global step number
         kinc   = iTer        ! Global iteration number
         iabort = 0
         ipl    = 0        
         
         dtime  = 1.0d0


c     Initialize the rotation tensor as Kronecker tensor

         do i=1,3
           do j=1,3
              if (i .eq. j) then
                drot(i,j) = 1.d0
              else
                drot(i,j) = 0.d0
              end if 
           enddo
         enddo
         
         if(IsUndr.eq.1) then
           if(Props(17) .gt. 1.d-10) then
       	     Write(1,*)'Cannot use Kw>0 and material type undrained'
       	     Write(1,*)'Kw=',Props(15),' IsUndr=',IsUndr
             iAbort = 123
           end if
         end if 	 

!         if(Props(2) .lt. 1.d-10) Props(2)=1

        If (IDTask .Eq. 1) Then ! initialize state variables

            call StVarIni (StVar0, Props, iel, int)
            
         End If ! IDTask = 1

         call CopyRVec (sig0, stress, ntens)
         call CopyRVec (dEps, dstran, ntens)
         call CopyRVec (StVar0, statev, nStatV)

         call ABQTOPX(stress)
         call ABQTOPX(dstran)
         call ABQTOPX(statev)
         
         call UMAT( stress, statev, ddsdde, sse, spd, scd,
     &              rpl, ddsddt, drplde, drpldt,
     &              stran, dstran, time, dtime,
     &              temp, dtemp, predef, dpred, cmname,
     &              ndi, nshr, ntens, nstatv, props, nprops,
     &              coords, drot, pnewdt, celent, dfgrd0, dfgrd1,
     &              noel, npt, layer, kspt, kstep, kinc)
     
c --------------------- for PLAXIS output

c         if (statev(15) .gt. 1.5d0 ) then
c            ipl = 2       ! tension cut-off point 
c         endif 
c         if (statev(16) .gt. 0.0d0 ) then
c            ipl = 1       ! failure point   
c         endif 

c --------------------- for PLAXIS output                  

         If (IDTask .Eq. 1) Then ! initialize state variables
      
            call CopyRVec( statev, StVar0, nStatV )

            do i = 1, ntens
              gstran(i) = statev(i)
            end do

         End If ! IDTask = 1

         If (IDTask .Eq. 2) Then ! Calculate stresses

            call CopyRVec( stress, Sig, 6 )
            call CopyRVec( statev, StVar, nStatV )
            call ABQTOPX(Sig)
            call ABQTOPX(StVar)

              If (IsUndr.Eq.1) Then
              dEpsV = dEps(1) + dEps(2) + dEps(3)
              dSwp = BulkW * dEpsV
              Swp = Swp0 + dSwp
            Else
              Swp = Swp0
            End If
	  
	        do i = 1, ntens
              gstran(i) = statev(i)
            end do

         End If ! IDTask = 2

         
         If (IDTask .Eq. 3 .Or.
     *       IDTask .Eq. 6     ) Then ! Calculate D-Matrix

             do i = 1, 6
               do j = 1, 6
                  D(i,j) = ddsdde(i,j)
               end do
             end do

            do i = 1, ntens
              gstran(i) = statev(i)
            end do

         End If  ! IDTask = 3, 6

      End If  ! IDTask = 1, 2, 3, 6

      If (IDTask .Eq. 4) Then ! Number of state parameters
              nStat    = nStatV
      End If  ! IDTask = 4

      If (IDTask .Eq. 5) Then ! matrix type
        NonSym   = 1  ! 1 for non-symmetric D-matrix
        iStrsDep = 1  ! 1 for stress dependent D-matrix
        iTang    = 1  ! 1 for tangent D-matrix
        iTimeDep = 0  ! 1 for time dependent D-matrix
      End If  ! IDTask = 5

      Return
      End ! MyMod_HC
      
      
