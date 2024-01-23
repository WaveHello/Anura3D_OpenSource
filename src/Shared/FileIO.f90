!*****************************************************************************
!                                       ____  _____
!           /\                         |___ \|  __ \
!          /  \   _ __  _   _ _ __ __ _  __) | |  | |
!         / /\ \ | '_ \| | | | '__/ _` ||__ <| |  | |
!        / ____ \| | | | |_| | | | (_| |___) | |__| |
!       /_/    \_\_| |_|\__,_|_|  \__,_|____/|_____/
!
!
! Anura3D - Numerical modelling and simulation of large deformations
!   and soil�water�structure interaction using the material point method (MPM)
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


module ModFileIO
   !**********************************************************************
   !
   !     $Revision: 9707 $
   !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr 2022) $
   !
   !**********************************************************************
   use ModGlobalConstants

   implicit none

contains


   !*************************************************************************************
   !    SUBROUTINE: FileOpenAction
   !
   !    DESCRIPTION:
   !>   Opens a file in binary format and carries out a read/write action on it.
   !
   !>   @param[in] FileUnit : file unit identifier
   !>   @param[in] FileName : name of the file that is opened
   !>   @param[in] Action : action to be carried out on the file ('R'=read, 'W'=write, 'RW'=readwrite)
   !
   !*************************************************************************************
   subroutine FileOpenAction(FileUnit, FileName, Action)

      implicit none

      integer(INTEGER_TYPE), intent(in) :: FileUnit
      character(len=*), intent(in) :: FileName, Action

      ! local variables
      integer(INTEGER_TYPE) :: ios
      logical :: IsOpen
      character(len=10) :: ToDo

      ToDo = 'ReadWrite' ! default
      if ( trim(Action) == 'R' ) ToDo = 'Read'
      if ( trim(Action) == 'W' ) ToDo = 'Write'
      if ( trim(Action) == 'RW' ) ToDo = 'ReadWrite'
      inquire(unit = FileUnit, opened = IsOpen)
      if ( IsOpen ) close(FileUnit)
      if ( trim(Action) == 'W' ) call EraseFile(FileName)

    #ifdef __INTEL_COMPILER
          open(FileUnit, FILE = FileName, FORM = 'Binary', BLOCKSIZE = 4096, RECORDTYPE = 'stream', BUFFERCOUNT = 64, BUFFERED = 'no', ACTION = ToDo, IOSTAT = ios)
    #else
          open(FileUnit, FILE = FileName, FORM = "UNFORMATTED", ACCESS = "SEQUENTIAL", STATUS = "REPLACE", IOSTAT = ios)
    #endif

      call Assert( ios == 0, 'Error opening file: ' // trim(FileName) // ' ' // trim(ToDo) )

   end subroutine FileOpenAction

   !*************************************************************************************
   !    SUBROUTINE: FileOpen
   !
   !    DESCRIPTION:
   !>   Opens a file.
   !
   !>   @param[in] FileUnit : file unit identifier
   !>   @param[in] FileName : name of the file that is opened
   !
   !*************************************************************************************
   subroutine FileOpen(FileUnit, FileName)
      implicit none

      integer(INTEGER_TYPE), intent(in) :: FileUnit
      character(len=*), intent(in) :: FileName

      ! local variables
      integer(INTEGER_TYPE) :: ios
      character(len=1023) :: NameIn, NameOut

      NameIn = Trim(FileName) // ' '
      call GetOrMakeFileName(NameIn, NameOut)
      open(FileUnit, FILE = NameOut, IOSTAT = ios)
      call Assert( ios == 0, 'Error opening file: ' // trim(FileName) )

   end subroutine FileOpen

   !*************************************************************************************
   !    SUBROUTINE: FileOpenWriteBinary
   !
   !    DESCRIPTION:
   !>   Opens a file in binary format for writing content.
   !
   !>   @param[in] FileUnit : file unit identifier
   !>   @param[in] FileName : name of the file that is opened
   !
   !*************************************************************************************
   subroutine FileOpenWriteBinary(FileUnit, FileName)

      implicit none

      integer(INTEGER_TYPE), intent(in) :: FileUnit
      character(len=*), intent(in) :: FileName

      ! local variables
      integer(INTEGER_TYPE) :: ios
      character(len=1023) :: NameIn, NameOut

      NameIn = Trim(FileName) // ' '
      call GetOrMakeFileName(NameIn, NameOut)
      open(FileUnit, FILE = NameOut, FORM = 'UNFORMATTED', ACCESS = 'SEQUENTIAL', ACTION = 'WRITE', CONVERT = 'BIG_ENDIAN', RECORDTYPE = 'STREAM', BUFFERED = 'YES', IOSTAT = ios)
      call Assert( ios == 0, 'Error opening file: ' // trim(FileName) )

   end subroutine FileOpenWriteBinary

   !*************************************************************************************
   !    SUBROUTINE: FileOpenAppend
   !
   !    DESCRIPTION:
   !>   Opens a file for appending its content.
   !
   !>   @param[in] FileUnit : file unit identifier
   !>   @param[in] FileName : name of the file that is opened
   !
   !*************************************************************************************
   subroutine FileOpenAppend(FileUnit, FileName)

      implicit none

      integer(INTEGER_TYPE), intent(in) :: FileUnit
      character(len=*), intent(in) :: FileName

      ! local variables
      integer(INTEGER_TYPE) :: ios
      character(len=1023) :: NameIn, NameOut

      NameIn = trim(FileName) // ' '
      call GetOrMakeFileName(NameIn, NameOut)
      open(FileUnit, FILE = NameOut, ACCESS = "APPEND", IOSTAT = ios)
      call Assert( ios == 0, 'Error opening file: ' // trim(FileName) )

   end subroutine FileOpenAppend

   !*************************************************************************************
   !    Function: FExist
   !
   !    DESCRIPTION: Check if File "FName" exists
   !
   !*************************************************************************************
   Logical Function FExist(FName)

      Character FName*(*)
      Logical Tmp
      Tmp=.False.
      goto 2

      Goto 1
      Open( 1,file=FName,Status='OLD',Err=1)
      Tmp=.True.
      Close(1)
1     Continue
      Inquire(File=FName,Exist=Tmp)
      FExist=Tmp
      Return

2     Call UniFExist( fName, Tmp )
      FExist=Tmp
      Return
   end function

   !*************************************************************************************
   !    Subroutine: FindIO
   !
   !    DESCRIPTION: Returns information of an Open File
   !
   !*************************************************************************************
   Subroutine FindIO(ioMin,io)


      logical  :: IsOpen
      ! integer(INTEGER_TYPE), parameter, intent(in) :: ioMax=99
      integer(INTEGER_TYPE), intent(in)    :: ioMin
      integer(INTEGER_TYPE), intent(inout) :: io
      io=ioMin-1
1     io=io+1
      Inquire(io,Opened=IsOpen)

   end subroutine

   !*************************************************************************************
   !    Subroutine: SKIP
   !
   !    DESCRIPTION: Skip reading bytes of block/unit
   !
   !*************************************************************************************
   SUBROUTINE SKIP(IUNIT,NBYTS)

      implicit none
      Character c*1

      integer(INTEGER_TYPE), parameter :: n1 = 1024, nT1 = 8*n1
      integer(INTEGER_TYPE), parameter :: n2 = 128, nT2 = 8*n2
      integer(INTEGER_TYPE), parameter :: n3 = 32, nT3 = 8*n3
      integer(INTEGER_TYPE), parameter :: n4 = 8, nT4 = 8*n4

      real(REAL_TYPE), dimension(n1) :: rBuf1
      real(REAL_TYPE), dimension(n2) :: rBuf2
      real(REAL_TYPE), dimension(n3) :: rBuf3
      real(REAL_TYPE), dimension(n4) :: rBuf4

      integer(INTEGER_TYPE), dimension(:), allocatable :: iBigBuf
      integer(INTEGER_TYPE):: I, IBUF, IERR, IDUM, IUNIT, NDUM, NBYTS

      If (-nByts > 100*nT1) Then
         iBuf = nByts
         Allocate( iBigBuf(iBuf), stat=ierr )
         Read(iUnit) iBigBuf
         nByts = nByts - iBuf
         DeAllocate( iBigBuf )
         Write(2,*)'iBuf : ',iBuf
      End If
      If (nByts > nT1) Then
         nDum = nByts / nT1
         Do i=1,nDum
            Read(iUnit) rBuf1
         End Do
         nByts = nByts - nT1 * nDum
      End If
      If (nByts > nT2) Then
         nDum = nByts / nT2
         Do i=1,nDum
            Read(iUnit) rBuf2
         End Do
         nByts = nByts - nT2 * nDum
      End If
      If (nByts > nT3) Then
         nDum = nByts / nT3
         Do i=1,nDum
            Read(iUnit) rBuf3
         End Do
         nByts = nByts - nT3 * nDum
      End If

      If (nByts > nT4) Then
         nDum = nByts / nT4
         Do i=1,nDum
            Read(iUnit) rBuf4
         End Do
         nByts = nByts - nT4 * nDum
      End If

      NDUM=NBYTS/4
      DO I=1,NDUM
         READ(IUNIT) IDUM
      End Do

      NByts=NByts-4*NDum
      If (NByts > 0) Then
         Do I=1,NByts
            READ(IUNIT) C
         End Do
      End If
   end subroutine

   !*************************************************************************************
   !    Subroutine: EraseFile
   !
   !    DESCRIPTION: Erase fName file
   !
   !*************************************************************************************
   Subroutine EraseFile(fName)

      Character fName*(*)
      
      ! local variables
      logical :: fExist
      integer(INTEGER_TYPE) :: io

      Call UniEraseFile ( fName )
      Return
      If (fExist(fName)) Then
         Call FindIO(10,io)
         Open(io,file=fName)
         Close(io,Status='Delete')
      End If
   end subroutine

   !*************************************************************************************
   !    Function: UpCase
   !
   !    DESCRIPTION: Write name in uppercase
   !
   !*************************************************************************************
   Character*255 Function UpCase(Lower)

      Character Lower*(*),Tmp*255,C*1
      integer(INTEGER_TYPE) :: i, LT, j
      Tmp=Lower
      LT=Len_Trim(Lower)
      Do i=1,LT
         C=Tmp(i:i)
         j=ichar(c)
         If (j >=  97 .And. j <= 122) Then
            Tmp(i:i)=Char(j-32)
         End If
      End Do
      Upcase=Tmp
   end function

   !*************************************************************************************
   !    Subroutine UniFExist
   !
   !    DESCRIPTION: check whether file with possible unicode name exists
   !
   !*************************************************************************************
   Subroutine UniFExist( Name, DoesExist )
      Character*(*) Name
      Logical DoesExist
      Character*1023 NameIn, NameOut
      NameIn = Name
      NameIn = Trim(Name)//' '
      Call Get_FileNameIfExists( NameIn, NameOut )
      DoesExist = Len_trim(NameOut) > 1
   end subroutine

   !*************************************************************************************
   !     Subroutine UniEraseFile
   !
   !    DESCRIPTION:  check whether file with possible unicode name exists,
   !                   when exists, delete it
   !
   !*************************************************************************************
   Subroutine UniEraseFile( Name )

      Character*(*) Name

      ! Local variables
      Logical DoesExist
      Character*1023 NameIn, NameOut
      integer(INTEGER_TYPE) :: io

      NameIn = Name
      NameIn = Trim(Name)//' '

      Call Get_FileNameIfExists( NameIn, NameOut )

      DoesExist = Len_trim(NameOut) > 1
      If (DoesExist) Then
         Call FindIO(10,io)
         Open(io,file=NameOut)
         Close(io,Status='Delete',err=1)
1        Continue
      End If
   end subroutine

   Subroutine Get_FileNameIfExists( Name, ShortName )
      !*************************************************************************************
      !     Subroutine Get_FileNameIfExists
      !
      !    DESCRIPTION:  interface routine to DLL routine to return the 'dos'-name
      !                  for an existing file with utf8-name
      !
      !*************************************************************************************

      Logical Tmp
      Character*(*) Name
      Character*1023 ShortName
      ShortName = ''
      ShortName = Name
      Inquire(File=Name,Exist=Tmp)
      If (.Not.Tmp) ShortName =  ' '

   end subroutine

   !*************************************************************************************
   !     Subroutine GetOrMakeFileName
   !
   !    DESCRIPTION:  Get or make the name of a File
   !
   !*************************************************************************************
   Subroutine GetOrMakeFileName( Name, ShortName )

      Character*(*) Name
      Character*1023 ShortName

      ShortName = Name

   end subroutine

  ! subroutine read_CPS_integer_value(FileUnit, )
  

  ! end subroutine

  ! subroutine read_logical_value()
  ! end subroutine

  ! subroutine read_real_value()
  ! end subroutine

  ! subroutine read_integer_value()
  ! ! Read an integer value from a file and store it in the passed variable
  ! ! The purpose of this subroutine is to read a line of file data it should
  ! ! For the time being it'll be based off the read 
  ! end subroutine
end module ModFileIO
