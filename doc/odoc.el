(setq petscfem-option-list '(
("A")
("A0")
("A1")
("A2")
("A3")
("A_0")
("A_scale")
("A_van_Driest")
("CN_ctff")
("CO_ctff")
("C_1")
("C_2")
("C_mu")
("C_smag")
("CdN_ctff")
("CdO_ctff")
("Chezy")
("Courant")
("Cp")
("Cr")
("D")
("Dt")
("G_body")
("KN")
("KO")
("KSP_method")
("Krylov_dim")
("LES")
("Nb_ctff")
("Nb_scale")
("Poisson_ratio")
("Pr_t")
("Rgas")
("Sc")
("T")
("T0")
("Tgas")
("Tinfty")
("Young_modulus")
("a_bar")
("activate_debug")
("activate_debug_memory_usage")
("activate_debug_print")
("additional_iprops")
("additional_props")
("additional_tau_pspg")
("advective_jacobians")
("advective_jacobians_type")
("alpha")
("aquifer_phi_dof")
("atol")
("auto_time_step")
("axisymmetric")
("base")
("beta_supg")
("betath")
("block_uploading")
("bottom_length")
("c_distor")
("c_relax")
("c_volume")
("cache_grad_div_u")
("called_from_rosi")
("check_dofmap_id")
("chunk_size")
("coef_file")
("compact_profile_graph_chunk_size")
("compute_fd_adv_jacobian")
("compute_fd_adv_jacobian_eps")
("compute_fd_adv_jacobian_random")
("compute_fd_adv_jacobian_rel_err_threshold")
("conductivity")
("consistent_supg_matrix")
("cyclic_fs")
("cyclic_length")
("debug_compute_prof")
("debug_element_partitioning")
("density")
("dev_comp_mask")
("diffusive_jacobians")
("diffusive_jacobians_type")
("diffusivity")
("displ_factor")
("distor_exp")
("double_layer")
("dry_aquifer_stop")
("dtol")
("dx")
("dx_auto_combine")
("dx_cache_connectivities")
("dx_cache_coords")
("dx_coords_scale_factor")
("dx_coords_scale_factor0")
("dx_do_make_command")
("dx_indices")
("dx_node_coordinates")
("dx_port")
("dx_read_state_from_file")
("dx_split_state")
("dx_state_all_fields")
("dx_steps")
("element_weight")
("enthalpy_jacobians")
("enthalpy_jacobians_type")
("eps_ctff_val")
("eps_min")
("epsilon_fdj")
("epsilon_x")
("expo")
("ext_filename")
("fic_dof")
("filename")
("flux_law_coefficient")
("fractional_step")
("free_surface_damp")
("free_surface_set_level_factor")
("friction_law")
("fs_debug")
("fs_eq_factor")
("fs_relax")
("fs_smoothing_coef")
("function_name")
("g_dir")
("gamma")
("gather_file")
("gather_length")
("gather_pos")
("geometry")
("gmres_orthogonalization")
("gravity")
("h_min")
("half_width")
("hfilm_coeff")
("hfilm_source")
("hm_fac")
("id_cn")
("id_cn1")
("id_fac")
("id_lumped_cn")
("id_lumped_cn1")
("identify_volume_elements")
("iisdmat_print_statistics")
("interface_full_preco_fill")
("interface_full_preco_maxits")
("interface_full_preco_pc")
("interface_full_preco_relax_factor")
("inwt_stop")
("jacobian_factor")
("kap_ctff_val")
("kappa")
("ket_min")
("lagrange_diagonal_factor")
("lagrange_residual_factor")
("lagrange_row_scale_factor")
("lagrange_scale_factor")
("launch_mesh_move")
("layers")
("local_solver")
("local_store")
("local_time_step")
("low_pass_filter")
("lumped_mass")
("lumped_wallke")
("mat_new_nonzero_allocation_err")
("max_partgraph_vertices")
("max_partgraph_vertices_proc")
("maxits")
("measure_performance")
("measure_performance_weight")
("moment_center")
("ndim")
("ndimel")
("ndimelf")
("nel_surf")
("newton_relaxation_factor")
("nfile")
("nfilt")
("ngather")
("nnwt")
("nnwt_liq")
("non_inertial_frame")
("non_inertial_frame_reverse_sign")
("normal_1d")
("npg")
("nrec")
("ns_id_cn")
("ns_id_cn1")
("ns_id_fac")
("ns_id_lumped_cn")
("ns_id_lumped_cn1")
("nsave")
("nsaverot")
("nsc")
("nsome")
("nstep")
("nstep_cpu_stat")
("nu_t")
("omega")
("omega_newton")
("p_thrsh")
("part_include_fic")
("partitioning_method")
("pc_lu_fill")
("peak")
("phieq")
("preco_side")
("preco_type")
("pressure_comp_mask")
("pressure_control_coef")
("print_Schur_matrix")
("print_dofmap_id")
("print_fsm_transition_info")
("print_interface_full_preco_conv")
("print_internal_loop_conv")
("print_linear_system_and_stop")
("print_local_chunk_size")
("print_partitioning_statistics")
("print_residual")
("print_some_file")
("pspg_advection_factor")
("pspg_factor")
("rain")
("reactive_jacobians")
("reactive_jacobians_type")
("real_nodes")
("report_assembly_time")
("report_consumed_time")
("report_consumed_time_stat")
("report_option_access")
("residual_factor")
("restart")
("reuse_mat")
("rho")
("rho_Cp")
("rho_thrsh")
("roughness")
("rtol")
("save_file")
("save_file_pattern")
("save_file_some")
("save_file_some_append")
("shape")
("shocap")
("shock_capturing")
("shock_capturing_factor")
("shock_capturing_threshold")
("sigma")
("sigma_e")
("sigma_k")
("solve_system")
("solver")
("solver_mom")
("source term")
("source_term_type")
("start_comp_time")
("start_time")
("state_ref")
("stdout_file")
("steady")
("stop_on_neg_val")
("t0")
("t_0")
("t_scale")
("tau_fac")
("temporal_stability_factor")
("time_step_stop")
("tol_linear")
("tol_mass")
("tol_newton")
("tol_steady")
("turb_prod_coef")
("turbulence_coef")
("u0")
("update_jacobian_iters")
("update_jacobian_start_iters")
("update_jacobian_start_steps")
("update_jacobian_steps")
("use_compact_profile")
("use_exterior_normal")
("use_iisd")
("use_interface_full_preco")
("use_interface_full_preco_nlay")
("use_log_vars")
("vel_min")
("viscosity")
("volume_elemset")
("volume_gather_pos")
("volume_ref")
("volume_relax_coef")
("von_Karman_cnst")
("weak_form")
("wet_aquifer_width_min")
("y_wall")
("y_wall_plus")
))
