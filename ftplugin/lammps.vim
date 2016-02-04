if exists('*s:ShowLMPDoc') && g:lmpdoc_perform_mappings
    call s:PerformMappings()
    finish
endif

if !exists('g:lmpdoc_perform_mappings')
    let g:lmpdoc_perform_mappings = 1
endif

if !exists('g:lmpdoc_highlight')
    let g:lmpdoc_highlight = 1
endif

if !exists('g:lmpdoc_path')
    let g:lmpdoc_path = "/opt/lammps/doc/"
endif

if !exists('g:lmpdoc_cmd')
    let g:lmpdoc_cmd = 'tail -n +13 '. g:lmpdoc_path . ''
endif

if !exists('g:lmpdoc_open_cmd')
    let g:lmpdoc_open_cmd = 'split'
endif

setlocal switchbuf=useopen
highlight lmpdoc cterm=reverse gui=reverse
"{{{ 1
let b:lmpdoc_map = {
      \ "angle_style charmm": "angle_charmm",
      \ "angle_style class2": "angle_class2",
      \ "angle_coeff": "angle_coeff",
      \ "angle_style cosine": "angle_cosine",
      \ "angle_style cosine/delta": "angle_cosine_delta",
      \ "angle_style cosine/periodic": "angle_cosine_periodic",
      \ "angle_style cosine/shift": "angle_cosine_shift",
      \ "angle_style cosine/shift_exp": "angle_cosine_shift_exp",
      \ "angle_style cosine/squared": "angle_cosine_squared",
      \ "angle_style dipole": "angle_dipole",
      \ "angle_style fourier": "angle_fourier",
      \ "angle_style fourier_simple": "angle_fourier_simple",
      \ "angle_style harmonic": "angle_harmonic",
      \ "angle_style hybrid": "angle_hybrid",
      \ "angle_style none": "angle_none",
      \ "angle_style quartic": "angle_quartic",
      \ "angle_style sdk": "angle_sdk",
      \ "angle_style style": "angle_style",
      \ "angle_style table": "angle_table",
      \ "bond_style class2": "bond_class2",
      \ "bond_coeff": "bond_coeff",
      \ "bond_style fene": "bond_fene",
      \ "bond_style fene/expand": "bond_fene_expand",
      \ "bond_style harmonic": "bond_harmonic",
      \ "bond_style harmonic/shift": "bond_harmonic_shift",
      \ "bond_style harmonic/shift/cut": "bond_harmonic_shift_cut",
      \ "bond_style hybrid": "bond_hybrid",
      \ "bond_style morse": "bond_morse",
      \ "bond_style none": "bond_none",
      \ "bond_style nonlinear": "bond_nonlinear",
      \ "bond_style quartic": "bond_quartic",
      \ "bond_style table": "bond_table",
      \ "compute ackland/atom": "compute_ackland_atom",
      \ "compute angle/local": "compute_angle_local",
      \ "compute angmom/chunk": "compute_angmom_chunk",
      \ "compute basal/atom": "compute_basal_atom",
      \ "compute body/local": "compute_body_local",
      \ "compute bond/local": "compute_bond_local",
      \ "compute centro/atom": "compute_centro_atom",
      \ "compute chunk/atom": "compute_chunk_atom",
      \ "compute cluster/atom": "compute_cluster_atom",
      \ "compute cna/atom": "compute_cna_atom",
      \ "compute com": "compute_com",
      \ "compute com/chunk": "compute_com_chunk",
      \ "compute contact/atom": "compute_contact_atom",
      \ "compute coord/atom": "compute_coord_atom",
      \ "compute damage/atom": "compute_damage_atom",
      \ "compute dihedral/local": "compute_dihedral_local",
      \ "compute dilatation/atom": "compute_dilatation_atom",
      \ "compute displace/atom": "compute_displace_atom",
      \ "compute erotate/asphere": "compute_erotate_asphere",
      \ "compute erotate/rigid": "compute_erotate_rigid",
      \ "compute erotate/sphere": "compute_erotate_sphere",
      \ "compute erotate/sphere/atom": "compute_erotate_sphere_atom",
      \ "compute event/displace": "compute_event_displace",
      \ "compute fep": "compute_fep",
      \ "compute group/group": "compute_group_group",
      \ "compute gyration": "compute_gyration",
      \ "compute gyration/chunk": "compute_gyration_chunk",
      \ "compute heat/flux": "compute_heat_flux",
      \ "compute hexorder/atom": "compute_hexorder_atom",
      \ "compute improper/local": "compute_improper_local",
      \ "compute inertia/chunk": "compute_inertia_chunk",
      \ "compute ke": "compute_ke",
      \ "compute ke/atom": "compute_ke_atom",
      \ "compute ke/atom/eff": "compute_ke_atom_eff",
      \ "compute ke/eff": "compute_ke_eff",
      \ "compute ke/rigid": "compute_ke_rigid",
      \ "compute meso/e/atom": "compute_meso_e_atom",
      \ "compute meso/rho/atom": "compute_meso_rho_atom",
      \ "compute meso/t/atom": "compute_meso_t_atom",
      \ "compute modify": "compute_modify",
      \ "compute msd": "compute_msd",
      \ "compute msd/chunk": "compute_msd_chunk",
      \ "compute msd/nongauss": "compute_msd_nongauss",
      \ "compute omega/chunk": "compute_omega_chunk",
      \ "compute orientorder/atom": "compute_orientorder_atom",
      \ "compute pair": "compute_pair",
      \ "compute pair/local": "compute_pair_local",
      \ "compute pe": "compute_pe",
      \ "compute pe/atom": "compute_pe_atom",
      \ "compute plasticity/atom": "compute_plasticity_atom",
      \ "compute pressure": "compute_pressure",
      \ "compute property/atom": "compute_property_atom",
      \ "compute property/chunk": "compute_property_chunk",
      \ "compute property/local": "compute_property_local",
      \ "compute rdf": "compute_rdf",
      \ "compute reduce": "compute_reduce",
      \ "compute saed": "compute_saed",
      \ "compute slice": "compute_slice",
      \ "compute smd/contact/radius": "compute_smd_contact_radius",
      \ "compute smd/damage": "compute_smd_damage",
      \ "compute smd/hourglass/error": "compute_smd_hourglass_error",
      \ "compute smd/internal/energy": "compute_smd_internal_energy",
      \ "compute smd/plastic/strain": "compute_smd_plastic_strain",
      \ "compute smd/plastic/strain/rate": "compute_smd_plastic_strain_rate",
      \ "compute smd/rho": "compute_smd_rho",
      \ "compute smd/tlsph/defgrad": "compute_smd_tlsph_defgrad",
      \ "compute smd/tlsph/dt": "compute_smd_tlsph_dt",
      \ "compute smd/tlsph/num/neighs": "compute_smd_tlsph_num_neighs",
      \ "compute smd/tlsph/shape": "compute_smd_tlsph_shape",
      \ "compute smd/tlsph/strain": "compute_smd_tlsph_strain",
      \ "compute smd/tlsph/strain/rate": "compute_smd_tlsph_strain_rate",
      \ "compute smd/tlsph/stress": "compute_smd_tlsph_stress",
      \ "compute smd/triangle/mesh/vertices": "compute_smd_triangle_mesh_vertices",
      \ "compute smd/ulsph/num/neighs": "compute_smd_ulsph_num_neighs",
      \ "compute smd/ulsph/strain": "compute_smd_ulsph_strain",
      \ "compute smd/ulsph/strain/rate": "compute_smd_ulsph_strain_rate",
      \ "compute smd/ulsph/stress": "compute_smd_ulsph_stress",
      \ "compute smd/vol": "compute_smd_vol",
      \ "compute sna/atom": "compute_sna_atom",
      \ "compute stress/atom": "compute_stress_atom",
      \ "compute tally": "compute_tally",
      \ "compute temp": "compute_temp",
      \ "compute temp/asphere": "compute_temp_asphere",
      \ "compute temp/chunk": "compute_temp_chunk",
      \ "compute temp/com": "compute_temp_com",
      \ "compute temp/cs": "compute_temp_cs",
      \ "compute temp/deform": "compute_temp_deform",
      \ "compute temp/deform/eff": "compute_temp_deform_eff",
      \ "compute temp/drude": "compute_temp_drude",
      \ "compute temp/eff": "compute_temp_eff",
      \ "compute temp/partial": "compute_temp_partial",
      \ "compute temp/profile": "compute_temp_profile",
      \ "compute temp/ramp": "compute_temp_ramp",
      \ "compute temp/region": "compute_temp_region",
      \ "compute temp/region/eff": "compute_temp_region_eff",
      \ "compute temp/rotate": "compute_temp_rotate",
      \ "compute temp/sphere": "compute_temp_sphere",
      \ "compute ti": "compute_ti",
      \ "compute torque/chunk": "compute_torque_chunk",
      \ "compute vacf": "compute_vacf",
      \ "compute vcm/chunk": "compute_vcm_chunk",
      \ "compute voronoi/atom": "compute_voronoi_atom",
      \ "compute xrd": "compute_xrd",
      \ "dihedral_style": "dihedral_style",
      \ "dihedral_style charmm": "dihedral_charmm",
      \ "dihedral_style class2": "dihedral_class2",
      \ "dihedral_coeff": "dihedral_coeff",
      \ "dihedral_style cosine/shift/exp": "dihedral_cosine_shift_exp",
      \ "dihedral_style fourier": "dihedral_fourier",
      \ "dihedral_style harmonic": "dihedral_harmonic",
      \ "dihedral_style helix": "dihedral_helix",
      \ "dihedral_style hybrid": "dihedral_hybrid",
      \ "dihedral_style multi/harmonic": "dihedral_multi_harmonic",
      \ "dihedral_style nharmonic": "dihedral_nharmonic",
      \ "dihedral_style none": "dihedral_none",
      \ "dihedral_style opls": "dihedral_opls",
      \ "dihedral_style quadratic": "dihedral_quadratic",
      \ "dihedral_style table": "dihedral_table",
      \ "fix adapt": "fix_adapt",
      \ "fix adapt/fep": "fix_adapt_fep",
      \ "fix addforce": "fix_addforce",
      \ "fix addtorque": "fix_addtorque",
      \ "fix append/atoms": "fix_append_atoms",
      \ "fix atc": "fix_atc",
      \ "fix atom/swap": "fix_atom_swap",
      \ "fix ave/atom": "fix_ave_atom",
      \ "fix ave/chunk": "fix_ave_chunk",
      \ "fix ave/correlate": "fix_ave_correlate",
      \ "fix ave/correlate/long": "fix_ave_correlate_long",
      \ "fix ave/histo": "fix_ave_histo",
      \ "fix ave/spatial": "fix_ave_spatial",
      \ "fix ave/spatial/sphere": "fix_ave_spatial_sphere",
      \ "fix ave/time": "fix_ave_time",
      \ "fix aveforce": "fix_aveforce",
      \ "fix balance": "fix_balance",
      \ "fix bond/break": "fix_bond_break",
      \ "fix bond/create": "fix_bond_create",
      \ "fix bond/swap": "fix_bond_swap",
      \ "fix box/relax": "fix_box_relax",
      \ "fix colvars": "fix_colvars",
      \ "fix deform": "fix_deform",
      \ "fix deposit": "fix_deposit",
      \ "fix drag": "fix_drag",
      \ "fix drude": "fix_drude",
      \ "fix drude/transform": "fix_drude_transform",
      \ "fix dt/reset": "fix_dt_reset",
      \ "fix efield": "fix_efield",
      \ "fix enforce2d": "fix_enforce2d",
      \ "fix evaporate": "fix_evaporate",
      \ "fix external": "fix_external",
      \ "fix freeze": "fix_freeze",
      \ "fix gcmc": "fix_gcmc",
      \ "fix gld": "fix_gld",
      \ "fix gle": "fix_gle",
      \ "fix gravity": "fix_gravity",
      \ "fix heat": "fix_heat",
      \ "fix imd": "fix_imd",
      \ "fix indent": "fix_indent",
      \ "fix ipi": "fix_ipi",
      \ "fix langevin": "fix_langevin",
      \ "fix langevin/drude": "fix_langevin_drude",
      \ "fix langevin/eff": "fix_langevin_eff",
      \ "fix lb/fluid": "fix_lb_fluid",
      \ "fix lb/momentum": "fix_lb_momentum",
      \ "fix lb/pc": "fix_lb_pc",
      \ "fix lb/rigid/pc/sphere": "fix_lb_rigid_pc_sphere",
      \ "fix lb/viscous": "fix_lb_viscous",
      \ "fix lineforce": "fix_lineforce",
      \ "fix meso": "fix_meso",
      \ "fix meso/stationary": "fix_meso_stationary",
      \ "fix modify": "fix_modify",
      \ "fix momentum": "fix_momentum",
      \ "fix move": "fix_move",
      \ "fix msst": "fix_msst",
      \ "fix neb": "fix_neb",
      \ "fix nh": "fix_nh",
      \ "fix nvt": "fix_nh",
      \ "fix npt": "fix_nh",
      \ "fix nph": "fix_nh",
      \ "fix nh/eff": "fix_nh_eff",
      \ "fix nph/asphere": "fix_nph_asphere",
      \ "fix nph/sphere": "fix_nph_sphere",
      \ "fix nphug": "fix_nphug",
      \ "fix npt/asphere": "fix_npt_asphere",
      \ "fix npt/sphere": "fix_npt_sphere",
      \ "fix nve": "fix_nve",
      \ "fix nve/asphere": "fix_nve_asphere",
      \ "fix nve/asphere/noforce": "fix_nve_asphere_noforce",
      \ "fix nve/body": "fix_nve_body",
      \ "fix nve/eff": "fix_nve_eff",
      \ "fix nve/limit": "fix_nve_limit",
      \ "fix nve/line": "fix_nve_line",
      \ "fix nve/noforce": "fix_nve_noforce",
      \ "fix nve/sphere": "fix_nve_sphere",
      \ "fix nve/tri": "fix_nve_tri",
      \ "fix nvt/asphere": "fix_nvt_asphere",
      \ "fix nvt/sllod": "fix_nvt_sllod",
      \ "fix nvt/sllod/eff": "fix_nvt_sllod_eff",
      \ "fix nvt/sphere": "fix_nvt_sphere",
      \ "fix oneway": "fix_oneway",
      \ "fix orient/fcc": "fix_orient_fcc",
      \ "fix phonon": "fix_phonon",
      \ "fix pimd": "fix_pimd",
      \ "fix planeforce": "fix_planeforce",
      \ "fix poems": "fix_poems",
      \ "fix pour": "fix_pour",
      \ "fix press/berendsen": "fix_press_berendsen",
      \ "fix print": "fix_print",
      \ "fix property/atom": "fix_property_atom",
      \ "fix qbmsst": "fix_qbmsst",
      \ "fix qeq": "fix_qeq",
      \ "fix qeq/comb": "fix_qeq_comb",
      \ "fix qeq/reax": "fix_qeq_reax",
      \ "fix qmmm": "fix_qmmm",
      \ "fix qtb": "fix_qtb",
      \ "fix reax/bonds": "fix_reax_bonds",
      \ "fix reax/c/species": "fix_reaxc_species",
      \ "fix recenter": "fix_recenter",
      \ "fix restrain": "fix_restrain",
      \ "fix rigid": "fix_rigid",
      \ "fix saed/vtk": "fix_saed_vtk",
      \ "fix setforce": "fix_setforce",
      \ "fix shake": "fix_shake",
      \ "fix smd": "fix_smd",
      \ "fix smd/adjust/dt": "fix_smd_adjust_dt",
      \ "fix smd/integrate/tlsph": "fix_smd_integrate_tlsph",
      \ "fix smd/integrate/ulsph": "fix_smd_integrate_ulsph",
      \ "fix smd/move/triangulated/surface": "fix_smd_move_triangulated_surface",
      \ "fix smd/setvel": "fix_smd_setvel",
      \ "fix smd/tlsph/reference/configuration": "fix_smd_tlsph_reference_configuration",
      \ "fix smd/wall/surface": "fix_smd_wall_surface",
      \ "fix spring": "fix_spring",
      \ "fix spring/rg": "fix_spring_rg",
      \ "fix spring/self": "fix_spring_self",
      \ "fix srd": "fix_srd",
      \ "fix store/force": "fix_store_force",
      \ "fix store/state": "fix_store_state",
      \ "fix temp/berendsen": "fix_temp_berendsen",
      \ "fix temp/csvr": "fix_temp_csvr",
      \ "fix temp/rescale": "fix_temp_rescale",
      \ "fix temp/rescale/eff": "fix_temp_rescale_eff",
      \ "fix tfmc": "fix_tfmc",
      \ "fix thermal/conductivity": "fix_thermal_conductivity",
      \ "fix ti/rs": "fix_ti_rs",
      \ "fix ti/spring": "fix_ti_spring",
      \ "fix tmd": "fix_tmd",
      \ "fix ttm": "fix_ttm",
      \ "fix tune/kspace": "fix_tune_kspace",
      \ "fix vector": "fix_vector",
      \ "fix viscosity": "fix_viscosity",
      \ "fix viscous": "fix_viscous",
      \ "fix wall": "fix_wall",
      \ "fix wall/gran": "fix_wall_gran",
      \ "fix wall/piston": "fix_wall_piston",
      \ "fix wall/reflect": "fix_wall_reflect",
      \ "fix wall/region": "fix_wall_region",
      \ "fix wall/srd": "fix_wall_srd",
      \ "pair_style adp": "pair_adp",
      \ "pair_style airebo": "pair_airebo",
      \ "pair_style awpmd": "pair_awpmd",
      \ "pair_style beck": "pair_beck",
      \ "pair_style body": "pair_body",
      \ "pair_style bop": "pair_bop",
      \ "pair_style born": "pair_born",
      \ "pair_style brownian": "pair_brownian",
      \ "pair_style buck": "pair_buck",
      \ "pair_style buck/long": "pair_buck_long",
      \ "pair_style charmm": "pair_charmm",
      \ "pair_style class2": "pair_class2",
      \ "pair_coeff": "pair_coeff",
      \ "pair_style colloid": "pair_colloid",
      \ "pair_style comb": "pair_comb",
      \ "pair_style coul": "pair_coul",
      \ "pair_style coul/diel": "pair_coul_diel",
      \ "pair_style cs": "pair_cs",
      \ "pair_style dipole": "pair_dipole",
      \ "pair_style dpd": "pair_dpd",
      \ "pair_style dsmc": "pair_dsmc",
      \ "pair_style eam": "pair_eam",
      \ "pair_style edip": "pair_edip",
      \ "pair_style eff": "pair_eff",
      \ "pair_style eim": "pair_eim",
      \ "pair_style gauss": "pair_gauss",
      \ "pair_style gayberne": "pair_gayberne",
      \ "pair_style gran": "pair_gran",
      \ "pair_style gromacs": "pair_gromacs",
      \ "pair_style hbond/dreiding": "pair_hbond_dreiding",
      \ "pair_style hybrid": "pair_hybrid",
      \ "pair_style kim": "pair_kim",
      \ "pair_style lcbop": "pair_lcbop",
      \ "pair_style line/lj": "pair_line_lj",
      \ "pair_style list": "pair_list",
      \ "pair_style lj": "pair_lj",
      \ "pair_style lj96": "pair_lj96",
      \ "pair_style lj/cubic": "pair_lj_cubic",
      \ "pair_style lj/expand": "pair_lj_expand",
      \ "pair_style lj/long": "pair_lj_long",
      \ "pair_style lj/sf": "pair_lj_sf",
      \ "pair_style lj/smooth": "pair_lj_smooth",
      \ "pair_style lj/smooth/linear": "pair_lj_smooth_linear",
      \ "pair_style lj/soft": "pair_lj_soft",
      \ "pair_style lubricate": "pair_lubricate",
      \ "pair_style lubricateU": "pair_lubricateU",
      \ "pair_style meam": "pair_meam",
      \ "pair_style meam/spline": "pair_meam_spline",
      \ "pair_style meam/sw/spline": "pair_meam_sw_spline",
      \ "pair_style mgpt": "pair_mgpt",
      \ "pair_style mie": "pair_mie",
      \ "pair_style modify": "pair_modify",
      \ "pair_style morse": "pair_morse",
      \ "pair_style nb3b/harmonic": "pair_nb3b_harmonic",
      \ "pair_style nm": "pair_nm",
      \ "pair_style none": "pair_none",
      \ "pair_style peri": "pair_peri",
      \ "pair_style polymorphic": "pair_polymorphic",
      \ "pair_style quip": "pair_quip",
      \ "pair_style reax": "pair_reax",
      \ "pair_style reax/c": "pair_reax_c",
      \ "pair_style resquared": "pair_resquared",
      \ "pair_style sdk": "pair_sdk",
      \ "pair_style smd/hertz": "pair_smd_hertz",
      \ "pair_style smd/tlsph": "pair_smd_tlsph",
      \ "pair_style smd/triangulated_surface": "pair_smd_triangulated_surface",
      \ "pair_style smd/ulsph": "pair_smd_ulsph",
      \ "pair_style smtbq": "pair_smtbq",
      \ "pair_style snap": "pair_snap",
      \ "pair_style soft": "pair_soft",
      \ "pair_style sph/heatconduction": "pair_sph_heatconduction",
      \ "pair_style sph/idealgas": "pair_sph_idealgas",
      \ "pair_style sph/lj": "pair_sph_lj",
      \ "pair_style sph/rhosum": "pair_sph_rhosum",
      \ "pair_style sph/taitwater": "pair_sph_taitwater",
      \ "pair_style sph/taitwater/morris": "pair_sph_taitwater_morris",
      \ "pair_style srp": "pair_srp",
      \ "pair_style": "pair_style",
      \ "pair_style sw": "pair_sw",
      \ "pair_style table": "pair_table",
      \ "pair_style tersoff": "pair_tersoff",
      \ "pair_style tersoff/mod": "pair_tersoff_mod",
      \ "pair_style tersoff/zbl": "pair_tersoff_zbl",
      \ "pair_style thole": "pair_thole",
      \ "pair_style tri/lj": "pair_tri_lj",
      \ "pair_style vashishta": "pair_vashishta",
      \ "pair_style write": "pair_write",
      \ "pair_style yukawa": "pair_yukawa",
      \ "pair_style yukawa/colloid": "pair_yukawa_colloid",
      \ "pair_style zbl": "pair_zbl",
      \ "improper_style class2": "improper_class2",
      \ "improper_coeff": "improper_coeff",
      \ "improper_style cossq": "improper_cossq",
      \ "improper_style cvff": "improper_cvff",
      \ "improper_style fourier": "improper_fourier",
      \ "improper_style harmonic": "improper_harmonic",
      \ "improper_style hybrid": "improper_hybrid",
      \ "improper_style none": "improper_none",
      \ "improper_style ring": "improper_ring",
      \ "improper_style umbrella": "improper_umbrella",
      \ "improper_style": "improper_style",
      \ "atom_modify": "atom_modify",
      \ "atom_style": "atom_style",
      \ "balance": "balance",
      \ "body": "body",
      \ "bond_style": "bond_style",
      \ "boundary": "boundary",
      \ "box": "box",
      \ "change_box": "change_box",
      \ "clear": "clear",
      \ "comm_modify": "comm_modify",
      \ "comm_style": "comm_style",
      \ "compute": "compute",
      \ "create_atoms": "create_atoms",
      \ "create_bonds": "create_bonds",
      \ "create_box": "create_box",
      \ "delete_atoms": "delete_atoms",
      \ "delete_bonds": "delete_bonds",
      \ "dielectric": "dielectric",
      \ "dimension": "dimension",
      \ "displace_atoms": "displace_atoms",
      \ "dump": "dump",
      \ "dump_h5md": "dump_h5md",
      \ "dump_image": "dump_image",
      \ "dump_modify": "dump_modify",
      \ "dump_molfile": "dump_molfile",
      \ "echo": "echo",
      \ "fix": "fix",
      \ "group": "group",
      \ "group2ndx": "group2ndx",
      \ "if": "if",
      \ "include": "include",
      \ "info": "info",
      \ "jump": "jump",
      \ "kspace_modify": "kspace_modify",
      \ "kspace_style": "kspace_style",
      \ "label": "label",
      \ "lattice": "lattice",
      \ "log": "log",
      \ "mass": "mass",
      \ "min_modify": "min_modify",
      \ "min_style": "min_style",
      \ "minimize": "minimize",
      \ "molecule": "molecule",
      \ "neb": "neb",
      \ "neigh_modify": "neigh_modify",
      \ "neighbor": "neighbor",
      \ "newton": "newton",
      \ "next": "next",
      \ "package": "package",
      \ "partition": "partition",
      \ "prd": "prd",
      \ "print": "print",
      \ "processors": "processors",
      \ "python": "python",
      \ "quit": "quit",
      \ "read_data": "read_data",
      \ "read_dump": "read_dump",
      \ "read_restart": "read_restart",
      \ "region": "region",
      \ "replicate": "replicate",
      \ "rerun": "rerun",
      \ "reset_timestep": "reset_timestep",
      \ "restart": "restart",
      \ "run": "run",
      \ "run_style": "run_style",
      \ "set": "set",
      \ "shell": "shell",
      \ "special_bonds": "special_bonds",
      \ "suffix": "suffix",
      \ "tad": "tad",
      \ "temper": "temper",
      \ "thermo": "thermo",
      \ "thermo_modify": "thermo_modify",
      \ "thermo_style": "thermo_style",
      \ "timer": "timer",
      \ "timestep": "timestep",
      \ "uncompute": "uncompute",
      \ "undump": "undump",
      \ "unfix": "unfix",
      \ "units": "units",
      \ "variable": "variable",
      \ "velocity": "velocity",
      \ "write_data": "write_data",
      \ "write_dump": "write_dump",
      \ "write_restart": "write_restart"
      \ }
" 1}}}

function! s:GetWindowLine(value)
  if a:value < 1
    return float2nr(winheight(0)*a:value)
  else
    return a:value
  endif
endfunction

function! s:ShowLMPDoc(word, line)

  let s:line = getline(".")

  if s:line =~# 'fix\|compute'
    let s:tmp  = split(s:line)
    echom s:tmp[0] . ' ' . s:tmp[3]
    let s:line = b:lmpdoc_map[s:tmp[0] . ' ' . s:tmp[3]]
  else
    let s:line = b:lmpdoc_map[substitute(matchstr(s:line, '\(\S\+\s\+\)\?'.a:word), '\s\+', ' ', 'g')]
  endif

  if g:lmpdoc_open_cmd == 'split'
    if exists('g:lmpdoc_window_lines')
      let l:lmpdoc_wh = s:GetWindowLine(g:lmpdoc_window_lines)
    else
      let l:lmpdoc_wh = 10
    endif
  endif

  if bufloaded("__doc__")
    let l:buf_is_new = 0
    if bufname("%") == "__doc__"
      " The current buffer is __doc__, thus do not
      " recreate nor resize it
      let l:lmpdoc_wh = -1
    else
      " If the __doc__ buffer is open, jump to it
      if exists("g:lmpdoc_use_drop")
        execute "drop" "__doc__"
      else
        execute "sbuffer" bufnr("__doc__")
      endif
      let l:lmpdoc_wh = -1
    endif
  else
    let l:buf_is_new = 1
    execute g:lmpdoc_open_cmd '__doc__'
    if g:lmpdoc_perform_mappings
      call s:PerformMappings()
    endif
  endif

  setlocal modifiable
  setlocal noswapfile
  setlocal buftype=nofile
  setlocal bufhidden=delete
  setlocal syntax=man
  setlocal nolist

  normal ggdG
  " " Remove function/method arguments
  " let s:name2 = substitute(a:name, '(.*', '', 'g' )
  " Remove all colons
  let s:name = substitute(s:line, '\s\+', '_', 'g' )
  let s:cmd = g:lmpdoc_cmd . s:name . ".txt"

  if &verbose
    echomsg "lmpdoc: calling " s:cmd
  endif
  execute  "silent read !" s:cmd
  normal 1G

  if exists('l:lmpdoc_wh') && l:lmpdoc_wh != -1
    execute "resize" l:lmpdoc_wh
  end

  if g:lmpdoc_highlight == 1
    execute 'syntax match lmpdoc' "'" . s:name . "'"
  endif

  let l:line = getline(2)
  if l:line =~ "^no Lammps documentation found for.*$"
    if l:buf_is_new
      execute "bdelete!"
    else
      normal u
      setlocal nomodified
      setlocal nomodifiable
    endif
    redraw
    echohl WarningMsg | echo l:line | echohl None
  else
    setlocal nomodified
    setlocal nomodifiable
  endif
endfunction

" Mappings
function! s:PerformMappings()
  " remap the K (or 'help') key
  nnoremap <silent> <buffer> K :call <SID>ShowLMPDoc(expand('<cWORD>'), line("."))<CR>
endfunction

if g:lmpdoc_perform_mappings
  call s:PerformMappings()
endif
