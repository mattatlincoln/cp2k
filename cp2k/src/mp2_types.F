!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types needed for MP2 calculations
!> \par History
!>       2011.05 created [Mauro Del Ben]
!> \author MDB
! **************************************************************************************************
MODULE mp2_types
   USE cp_eri_mme_interface,            ONLY: cp_eri_mme_finalize,&
                                              cp_eri_mme_param
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type,&
                                              dbcsr_type
   USE hfx_types,                       ONLY: hfx_release,&
                                              hfx_type,&
                                              pair_list_element_type
   USE input_constants,                 ONLY: &
        do_eri_mme, mp2_method_direct, mp2_method_gpw, mp2_method_laplace, mp2_method_none, &
        mp2_ri_optimize_basis, ri_mp2_laplace, ri_mp2_method_gpw, ri_rpa_method_gpw
   USE kinds,                           ONLY: default_path_length,&
                                              dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_types'

   PUBLIC :: mp2_type, &
             integ_mat_buffer_type, &
             integ_mat_buffer_type_2D, &
             mp2_method_none, &
             mp2_method_laplace, &
             mp2_method_direct, &
             mp2_method_gpw, &
             mp2_ri_optimize_basis, &
             ri_mp2_method_gpw, &
             ri_rpa_method_gpw, &
             ri_mp2_laplace, &
             init_TShPSC_lmax

   PUBLIC :: mp2_env_create, mp2_env_release, mp2_biel_type, &
             pair_list_type_mp2, &
             one_dim_int_array, &
             two_dim_int_array, &
             two_dim_real_array

   INTEGER, SAVE                                         :: init_TShPSC_lmax = -1

! TYPE definitions

   TYPE mp2_biel_type
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: index_table
   END TYPE mp2_biel_type

   TYPE mp2_laplace_type
      INTEGER       :: n_quadrature
      INTEGER       :: integ_group_size
      LOGICAL       :: mo_sos
      REAL(dp)      :: threshold
   END TYPE

   TYPE mp2_direct_type
      LOGICAL  :: big_send
   END TYPE

   TYPE mp2_potential_type
      INTEGER                                  :: potential_type !! 1/r, TShPSC, slab truncation
      REAL(dp)                                 :: truncation_radius !! truncation radius if cutoff potential in use
      CHARACTER(default_path_length)           :: filename
   END TYPE

   TYPE mp2_gpw_type
      REAL(KIND=dp)            :: eps_grid, eps_filter, eps_pgf_orb_S
      INTEGER                  :: print_level
      REAL(KIND=dp)            :: cutoff
      REAL(KIND=dp)            :: relative_cutoff
      LOGICAL                  :: do_mult_2c, do_chol_subgrp
      INTEGER                  :: size_chol_from_input
   END TYPE mp2_gpw_type

   TYPE ri_mp2_type
      INTEGER                  :: block_size
      REAL(dp)                 :: eps_canonical
      LOGICAL                  :: free_hfx_buffer
   END TYPE

   TYPE ri_rpa_type
      INTEGER                  :: rpa_num_quad_points
      INTEGER                  :: rpa_integ_group_size
      INTEGER                  :: mm_style
      TYPE(hfx_type), DIMENSION(:, :), POINTER &
         :: x_data
      LOGICAL                  :: minimax_quad
      LOGICAL                  :: do_ri_g0w0
      LOGICAL                  :: do_admm
      LOGICAL                  :: do_im_time
      LOGICAL                  :: do_RI_overlap
      LOGICAL                  :: do_ri_axk
      TYPE(dbcsr_type), POINTER             :: mo_coeff_o, &
                                               mo_coeff_v
      REAL(KIND=dp)            :: ener_axk
   END TYPE

   TYPE ri_rpa_im_time_util
      INTEGER                  :: color_sub_row
      INTEGER                  :: color_sub_col
      INTEGER                  :: color_sub_P
      INTEGER                  :: n_group_row
      INTEGER                  :: n_group_col
      INTEGER                  :: n_group_P
      INTEGER                  :: cut_RI
      INTEGER                  :: n_group_RI_orig
      INTEGER, ALLOCATABLE, DIMENSION(:)  :: sizes_array_prim_row, starts_array_prim_row, &
                                             ends_array_prim_row, sizes_array_prim_col, &
                                             starts_array_prim_col, ends_array_prim_col, &
                                             sizes_array_prim_fullrow, starts_array_prim_fullrow, &
                                             ends_array_prim_fullrow, sizes_array_prim_fullcol, &
                                             starts_array_prim_fullcol, ends_array_prim_fullcol, &
                                             start_row_data_block, end_row_data_block, &
                                             start_col_data_block, end_col_data_block, &
                                             my_group_L_starts_im_time, my_group_L_ends_im_time, &
                                             my_group_L_sizes_im_time, mepos_P_from_RI_row, &
                                             row_blk_offset_RI
   END TYPE

   TYPE ri_rpa_im_time_2d_util
      INTEGER, ALLOCATABLE, DIMENSION(:, :)  :: offset_combi_block
   END TYPE

   TYPE ri_rpa_im_time_type
      INTEGER                  :: cut_memory
      LOGICAL                  :: memory_info, do_mao, opt_sc_dm_occ, opt_sc_dm_virt, do_lr
      REAL(KIND=dp)            :: eps_filter_im_time, eps_grad_occ, eps_grad_virt, lr_fraction
      INTEGER                  :: group_size_P, group_size_3c, num_points_per_magnitude, &
                                  max_iter_occ, max_iter_virt, group_size_Q
      INTEGER, ALLOCATABLE, DIMENSION(:) :: sizes_array_cm, starts_array_cm, &
                                            ends_array_cm, sizes_array_cm_mao_occ, &
                                            starts_array_cm_mao_occ, ends_array_cm_mao_occ, &
                                            sizes_array_cm_mao_virt, starts_array_cm_mao_virt, &
                                            ends_array_cm_mao_virt
      INTEGER, DIMENSION(:), POINTER     :: nmao_occ, nmao_virt
      LOGICAL                  :: do_gw_im_time
      REAL(KIND=dp)            :: stabilize_exp
   END TYPE

   TYPE ri_g0w0_type
      INTEGER                  :: corr_mos_occ
      INTEGER                  :: corr_mos_virt
      INTEGER                  :: corr_mos_occ_beta
      INTEGER                  :: corr_mos_virt_beta
      REAL(KIND=dp)            :: scaling
      INTEGER                  :: num_poles
      INTEGER                  :: nparam_pade
      INTEGER                  :: analytic_continuation
      REAL(KIND=dp)            :: omega_max_fit
      REAL(KIND=dp)            :: stop_crit
      INTEGER                  :: max_iter_fit
      LOGICAL                  :: check_fit
      INTEGER                  :: crossing_search
      REAL(KIND=dp)            :: fermi_level_offset
      INTEGER                  :: iter_ev_sc
      REAL(KIND=dp)            :: eps_ev_sc_iter
      LOGICAL                  :: hf_like_ev_start, &
                                  print_gw_details, &
                                  do_ri_Sigma_x, &
                                  remove_neg_virt_energies
      INTEGER                  :: ri_metric
      LOGICAL                  :: mix_exchange, do_periodic
      REAL(KIND=dp)            :: frac_exx
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: vec_Sigma_x_minus_vxc_gw
      INTEGER                  :: contour_def_start, contour_def_end
      REAL(KIND=dp)            :: contour_def_offset
      INTEGER, DIMENSION(:), POINTER    :: kp_grid
      INTEGER                  :: num_kp_grids
      REAL(KIND=dp)            :: eps_kpoint
      LOGICAL                  :: do_mo_coeff_gamma, do_average_deg_levels
      REAL(KIND=dp)            :: eps_eigenval
      LOGICAL                  :: do_extra_kpoints, do_aux_bas_gw
      REAL(KIND=dp)            :: frac_aux_mos
      INTEGER                  :: num_omega_points
      LOGICAL                  :: do_bse, normalize_sigma, neglect_norm
      INTEGER, DIMENSION(:), POINTER    :: atoms_for_projection, atom_range_for_projection
      REAL(KIND=dp)            :: eps_charge
      INTEGER                  :: num_z_vectors, max_iter_bse
      REAL(KIND=dp)            :: eps_min_trans
      LOGICAL                  :: do_ic_model, print_ic_values, do_opt_homo_lumo
      REAL(KIND=dp)            :: eps_dist
      REAL(KIND=dp), DIMENSION(:), POINTER               :: ic_corr_list, ic_corr_list_beta, &
                                                            gw_eigenvalues, gw_eigenvalues_beta
   END TYPE

   TYPE ri_basis_opt
      REAL(KIND=dp)            :: DI_rel
      REAL(KIND=dp)            :: DRI
      REAL(KIND=dp)            :: eps_step
      INTEGER                  :: max_num_iter
      INTEGER                  :: basis_quality
      INTEGER, DIMENSION(:), ALLOCATABLE :: RI_nset_per_l
   END TYPE

   TYPE grad_util
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: P_ij, P_ab, P_ij_beta, P_ab_beta
      REAL(KIND=dp), DIMENSION(:, :, :), ALLOCATABLE :: Gamma_P_ia, Gamma_P_ia_beta
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: PQ_half, Gamma_PQ, Gamma_PQ_beta
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER :: G_P_ia, G_P_ia_beta
      TYPE(dbcsr_type), POINTER             :: mo_coeff_o, mo_coeff_o_beta, &
                                               mo_coeff_v, mo_coeff_v_beta
      TYPE(cp_fm_type), POINTER                :: P_mo, P_mo_beta, W_mo, W_mo_beta, &
                                                  L_jb, L_jb_beta
      REAL(KIND=dp) :: cphf_eps_conv
      INTEGER :: cphf_max_num_iter, num_sing_ij, num_sing_ij_beta
   END TYPE

   TYPE mp2_type
      REAL(KIND=dp)            :: energy
      LOGICAL                  :: compute_reference
      INTEGER                  :: method
      TYPE(mp2_laplace_type)   :: ri_laplace
      TYPE(mp2_direct_type)    :: direct_canonical
      TYPE(mp2_potential_type) :: potential_parameter
      TYPE(mp2_gpw_type)       :: mp2_gpw
      TYPE(ri_mp2_type)        :: ri_mp2
      TYPE(ri_rpa_type)        :: ri_rpa
      TYPE(ri_rpa_im_time_util), &
         DIMENSION(:), &
         ALLOCATABLE            :: ri_rpa_im_time_util
      TYPE(ri_rpa_im_time_2d_util), &
         DIMENSION(:, :), &
         ALLOCATABLE            :: ri_rpa_im_time_2d_util
      TYPE(ri_rpa_im_time_type) &
         :: ri_rpa_im_time
      TYPE(ri_g0w0_type)       :: ri_g0w0
      TYPE(ri_basis_opt)       :: ri_opt_param
      TYPE(grad_util)          :: ri_grad
      REAL(dp) :: mp2_memory
      REAL(dp) :: scale_S
      REAL(dp) :: scale_T
      REAL(dp) :: minimal_gap
      INTEGER  :: mp2_num_proc
      INTEGER  :: block_size_row
      INTEGER  :: block_size_col
      INTEGER  :: ri_metric
      LOGICAL  :: calc_PQ_cond_num
      LOGICAL  :: hf_fail
      LOGICAL  :: p_screen
      LOGICAL  :: not_last_hfx
      INTEGER  :: eri_method
      TYPE(cp_eri_mme_param), POINTER  :: eri_mme_param
   END TYPE

   TYPE integ_mat_buffer_type
      REAL(KIND=dp), DIMENSION(:), POINTER  :: msg
      INTEGER, DIMENSION(:), POINTER  :: sizes
      INTEGER, DIMENSION(:, :), POINTER  :: indx
      INTEGER :: proc
      INTEGER :: msg_req
   END TYPE

   TYPE integ_mat_buffer_type_2D
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: msg
      INTEGER :: proc
      INTEGER :: msg_req
   END TYPE

   TYPE pair_list_type_mp2
      TYPE(pair_list_element_type), DIMENSION(:), ALLOCATABLE :: elements
      INTEGER :: n_element
   END TYPE pair_list_type_mp2

   TYPE offset_type
      INTEGER, DIMENSION(:), POINTER :: offset_occ
      INTEGER, DIMENSION(:), POINTER :: LLL
   END TYPE

   TYPE one_dim_int_array
      INTEGER, DIMENSION(:), ALLOCATABLE    :: array
   END TYPE

   TYPE two_dim_int_array
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: array
   END TYPE

   TYPE two_dim_real_array
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: array
   END TYPE

   TYPE offset_row_col_from_LLL
      INTEGER, DIMENSION(:), POINTER :: rows
      INTEGER, DIMENSION(:), POINTER :: cols
      INTEGER, DIMENSION(:), POINTER :: offset
      INTEGER, DIMENSION(:), POINTER :: blocksize
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param mp2_env ...
! **************************************************************************************************
   SUBROUTINE mp2_env_release(mp2_env)
      TYPE(mp2_type), POINTER                            :: mp2_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_env_release', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(mp2_env))

      ! release the HFX section for the EXX calculation
      IF (ASSOCIATED(mp2_env%ri_rpa%x_data)) THEN
         CALL hfx_release(mp2_env%ri_rpa%x_data)
      END IF

      IF (mp2_env%eri_method .EQ. do_eri_mme) CALL cp_eri_mme_finalize(mp2_env%eri_mme_param)
      IF (ASSOCIATED(mp2_env%eri_mme_param)) DEALLOCATE (mp2_env%eri_mme_param)

      DEALLOCATE (mp2_env)

      CALL timestop(handle)

   END SUBROUTINE mp2_env_release

! **************************************************************************************************
!> \brief ...
!> \param mp2_env ...
! **************************************************************************************************
   SUBROUTINE mp2_env_create(mp2_env)
      TYPE(mp2_type), POINTER                            :: mp2_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_env_create', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CPASSERT(.NOT. ASSOCIATED(mp2_env))

      ALLOCATE (mp2_env)

      NULLIFY (mp2_env%ri_rpa%x_data)

      CALL timestop(handle)

   END SUBROUTINE mp2_env_create

END MODULE mp2_types

