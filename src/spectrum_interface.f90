! \file simulation.hh
! \brief Declares application level classes to perform a complete simulation from start to finish
! \author Vladimir Florinski

! This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.

subroutine spectrum_init_mpi(comm) bind(C)

   use BATL_mpi, ONLY: init_mpi

   implicit none

   integer, value, intent(in) :: comm

   call init_mpi(comm)

end subroutine spectrum_init_mpi

!=====================================================================================================================================================

subroutine spectrum_get_node(pos, node) bind(C)

   use BATL_lib, ONLY: nDim, CoordMin_D, DomainSize_D
   use BATL_tree, ONLY: find_tree_node
   use iso_c_binding, ONLY: c_int, c_double

   implicit none

   real(c_double), intent(in) :: pos(nDim)
   integer(c_int), intent(out) :: node
   real :: CoordTree_D(nDim)

! Normalize the distance to the domain size and call the tree search routine. The result is available on all CPUs.
   CoordTree_D = (pos - CoordMin_D) / (DomainSize_D)
   
!   write(*,*) CoordTree_D
!   write(*,*) CoordMin_D
!   write(*,*) DomainSize_D
   
   call find_tree_node(CoordTree_D, node)

end subroutine spectrum_get_node

!=====================================================================================================================================================

subroutine spectrum_get_neighbor_node(node, i, j, k, neighbor_node, neighbor_level) bind(C)

   use BATL_lib, ONLY: iTree_IA, iNodeNei_IIIB, DiLevelNei_IIIB, iProc, Block_, Proc_, Unset_
   use iso_c_binding, ONLY: c_int

   implicit none

   integer(c_int), value, intent(in) :: node, i, j, k
   integer(c_int), intent(out) :: neighbor_node, neighbor_level
   integer :: iBlock

! Make sure the input is correct
   if (i < 0 .or. i > 3 .or. j < 0 .or. j > 3 .or. k < 0 .or. k > 3) then
      neighbor_node = -1
      neighbor_level = 0
      return
   end if

! Neighbor information is only available locally
   iBlock = iTree_IA(Block_, node)
   if (iProc == iTree_IA(Proc_, node)) then
      neighbor_node = iNodeNei_IIIB(i, j, k, iBlock)
      neighbor_level = DiLevelNei_IIIB((i + 1) / 2 - 1, (j + 1) / 2 - 1, (k + 1) / 2 - 1, iBlock)
   else
      neighbor_node = -1
      neighbor_level = 0
   end if

end subroutine spectrum_get_neighbor_node

!=====================================================================================================================================================

!subroutine spectrum_get_all_neighbor_nodes(node, neighbor_nodes, neighbor_levels) bind(C)

!   use BATL_lib, ONLY: iProc, Block_, Proc_
!   use BATL_tree, ONLY: iTree_IA, iNodeNei_IIIB, DiLevelNei_IIIB
!   
!   use iso_c_binding, ONLY: c_int, c_ptr, c_null_ptr, c_loc

!   implicit none

!   integer(c_int), value, intent(in) :: node
!   type(c_ptr), intent(out) :: neighbor_nodes, neighbor_levels
!   integer :: iBlock

!! Pass a pointer to the neighbor structure. This is always 4 bytes long, so no conversion in necessary.
!   iBlock = iTree_IA(Block_, node)
!   if (iProc == iTree_IA(Proc_, node)) then
!      neighbor_nodes = c_loc(iNodeNei_IIIB(0, 0, 0, iBlock))
!      neighbor_levels = c_loc(DiLevelNei_IIIB(-1, -1, -1, iBlock))
!   else
!      neighbor_nodes = c_null_ptr
!      neighbor_levels = c_null_ptr
!   end if

!end subroutine spectrum_get_all_neighbor_nodes

!=====================================================================================================================================================

subroutine spectrum_get_all_neighbor_copies(node, neighbor_nodes, neighbor_levels) bind(C)

   use BATL_lib, ONLY: iTree_IA, iNodeNei_IIIB, DiLevelNei_IIIB, iProc, Block_, Proc_
   use iso_c_binding, ONLY: c_int, c_ptr, c_null_ptr, c_loc

   implicit none

   integer(c_int), value, intent(in) :: node
   integer(c_int), intent(out) :: neighbor_nodes(0 : 3, 0 : 3, 0 : 3)
   integer(c_int), intent(out) :: neighbor_levels(-1 : 1, -1 : 1, -1 : 1)
   integer :: iBlock

   iBlock = iTree_IA(Block_, node)
   if (iProc == iTree_IA(Proc_, node)) then
      neighbor_nodes = iNodeNei_IIIB(:, :, :, iBlock)
      neighbor_levels = DiLevelNei_IIIB(:, :, :, iBlock)
   else
      neighbor_nodes = -1
      neighbor_levels = 0
   end if

end subroutine spectrum_get_all_neighbor_copies

!=====================================================================================================================================================

subroutine spectrum_get_block_corners(node, face_min, face_max) bind(C)

   use BATL_lib, ONLY: iTree_IA, nDim, CoordMin_DB, CoordMax_DB, iProc, Block_, Proc_
   use iso_c_binding, ONLY: c_int, c_double

   implicit none

   integer(c_int), value, intent(in) :: node
   real(c_double), intent(out) :: face_min(nDim), face_max(nDim)
   integer :: iBlock
   
   iBlock = iTree_IA(Block_, node)
   if (iProc == iTree_IA(Proc_, node)) then
      face_min = CoordMin_DB(:, iBlock)
      face_max = CoordMax_DB(:, iBlock)
   else
      face_min = 0.0
      face_max = 0.0
   end if

end subroutine spectrum_get_block_corners

!=====================================================================================================================================================

subroutine spectrum_get_block_data(node, variables) bind(C)

   use BATL_lib, ONLY: iTree_IA, iProc, Proc_, Block_, MinI, MaxI, MinJ, MaxJ, MinK, MaxK
   use ModReadAmr, ONLY: State_VGB, nVar
   use iso_c_binding, ONLY: c_int, c_double

   implicit none

   integer(c_int), value, intent(in) :: node
   real(c_double), intent(out) :: variables(nVar, MinI : MaxI, MinJ : MaxJ, MinK : MaxK)
   integer :: iBlock, i, j, k, var

! This performs a copy of the portion of the data array.
   iBlock = iTree_IA(Block_, node)

   if (iProc == iTree_IA(Proc_, node)) then
      variables = State_VGB(:, :, :, :, iBlock)
   else
      variables = 0.0
   end if

end subroutine spectrum_get_block_data

!=====================================================================================================================================================

subroutine spectrum_get_interpolation_stencil(pos, n_zones, stencil_nodes, stencil_zones, stencil_weights) bind(C)

! FIXME: Interpolate_grid_amr does not work properly and neither does interpolate_grid_amr_gc.
!        Data obtained from requesting stencil and interpolating does not match data from BATL.
   use BATL_grid, ONLY: interpolate_grid_amr, interpolate_grid
   use BATL_tree, ONLY: iNode_B
   use BATL_size, ONLY: nDim, MaxDim
   
   use iso_c_binding, ONLY: c_int, c_double

   implicit none

   real(c_double), intent(in) :: pos(nDim)
   integer(c_int), intent(out) :: n_zones, stencil_nodes(2**nDim), stencil_zones(nDim * 2**nDim)
   real(c_double), intent(out) :: stencil_weights(2**nDim)

   integer :: nCell, iDim, iCell, idx
   integer :: iCell_II(0 : nDim, 2**nDim)
   real :: Xyz_D(MaxDim)
   real :: Weight_I(2**nDim)

   Xyz_D(1 : nDim) = pos(1 : nDim)
   ! call interpolate_grid_amr(Xyz_D, nCell, iCell_II, Weight_I)
   call interpolate_grid(Xyz_D, nCell, iCell_II, Weight_I)

   idx = 1
   do iCell = 1, nCell
      stencil_nodes(iCell) = iNode_B(iCell_II(0, iCell))
      do iDim = 1, nDim
         stencil_zones(idx) = iCell_II(iDim, iCell) - 1
         idx = idx + 1
      end do
   end do

   n_zones = nCell
   stencil_weights = Weight_I

end subroutine spectrum_get_interpolation_stencil

!=====================================================================================================================================================

! TODO this is not needed - use wrapamr
subroutine spectrum_read_header(NameFile_I, l, lVerbose) bind(C)

  ! Read header information from .h or .info
  ! This sets number and names of variables, domain size, etc.
  ! l is the maximum length of the filename
  ! If lVerbose is 1, write out verbose information.

  use ModUtilities, ONLY: char_array_to_string
  use ModReadAmr, ONLY:readamr_init

  use iso_c_binding, ONLY: c_char, c_int

  implicit none

  character(c_char), intent(in):: NameFile_I(*)
  integer(c_int), value, intent(in):: l
  integer(c_int), value, intent(in):: lVerbose

  character(len=l):: NameFile
  integer:: i
  !----------------------------------------------------------------------------
  call char_array_to_string(NameFile_I, NameFile)

  ! Cut off extension from the file name (if any)
  i = index(NameFile, ".", BACK=.true.)

  call readamr_init(NameFile(1:i-1), IsVerbose=lVerbose==1_c_int)

end subroutine spectrum_read_header

