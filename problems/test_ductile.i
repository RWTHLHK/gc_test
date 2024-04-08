#This input uses PhaseField-Nonconserved Action to add phase field fracture bulk rate kernels
[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 40
    ny = 20
    ymax = 0.5
  []
  [./noncrack]
    type = BoundingBoxNodeSetGenerator
    new_boundary = noncrack
    bottom_left = '0.5 0 0'
    top_right = '1 0 0'
    input = gen
  [../]
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[AuxVariables]
  [./strain_yy]
    family = MONOMIAL
    order = CONSTANT
  [../]
  # [./stress_yy]
  #   family = MONOMIAL
  #   order = CONSTANT
  # [../]
  [./elastic_strain_yy]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./plastic_strain_yy]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./uncracked_stress_yy]
    family = MONOMIAL
    order = CONSTANT
  [../]
[]

[Modules]
  [./TensorMechanics]
    [./Master]
      [./All]
        add_variables = true
        strain = FINITE
        planar_formulation = PLANE_STRAIN
        additional_generate_output = 'stress_yy vonmises_stress'
        strain_base_name = uncracked
        base_name = ductile_fractured
      [../]
    [../]
  [../]
  [./PhaseField]
    [./Nonconserved]
      [./c]
        free_energy = E_el
        kappa = kappa_op
        mobility = L
      [../]
    [../]
  [../]
[]

[Kernels]
  [./solid_x]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_x
    component = 0
    c = c
  [../]
  [./solid_y]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_y
    component = 1
    c = c
  [../]
  [./off_disp]
    type = AllenCahnElasticEnergyOffDiag
    variable = c
    displacements = 'disp_x disp_y'
    mob_name = L
  [../]
[]

[AuxKernels]
  [./strain_yy]
    type = RankTwoAux
    variable = strain_yy
    rank_two_tensor = uncracked_mechanical_strain
    index_i = 1
    index_j = 1
    execute_on = TIMESTEP_END
  [../]
  # [./stress_yy]
  #   type = RankTwoAux
  #   variable = stress_yy
  #   rank_two_tensor = ductile_fractured_stress
  #   index_i = 1
  #   index_j = 1
  #   execute_on = TIMESTEP_END
  # [../]
  [./elastic_strain_yy]
    type = RankTwoAux
    variable = elastic_strain_yy
    rank_two_tensor = uncracked_elastic_strain
    index_i = 1
    index_j = 1
    execute_on = TIMESTEP_END
  [../]
  [./plastic_strain_yy]
    type = RankTwoAux
    variable = plastic_strain_yy
    rank_two_tensor = uncracked_plastic_strain
    index_i = 1
    index_j = 1
    execute_on = TIMESTEP_END
  [../]
  [./uncracked_stress_yy]
    type = RankTwoAux
    variable = uncracked_stress_yy
    rank_two_tensor = uncracked_stress
    index_i = 1
    index_j = 1
    execute_on = TIMESTEP_END
  [../]
[]

[BCs]
  [./ydisp]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = 't'
  [../]
  [./yfix]
    type = DirichletBC
    variable = disp_y
    boundary = noncrack
    value = 0
  [../]
  [./xfix]
    type = DirichletBC
    variable = disp_x
    boundary = right
    value = 0
  [../]
[]

[Functions]
  [./hf]
    type = PiecewiseLinear
    x = '0    0.001 0.003 0.023'
    y = '0.85 1.0   1.25  1.5'
  [../]
[]

[Materials]
  [./pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'gc_prop l visco'
    prop_values = '1e-2 0.05 5e-3'
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    C_ijkl = '120.0 80.0'
    fill_method = symmetric_isotropic
    base_name = uncracked
  [../]
  [./isotropic_plasticity]
    type = IsotropicPlasticityStressUpdate
    yield_stress = 0.85
    hardening_function = hf
    base_name = uncracked
  [../]
  [./radial_return_stress]
    type = ComputeMultipleInelasticStress
    tangent_operator = elastic
    inelastic_models = 'isotropic_plasticity'
    base_name = uncracked
  [../]
  # [./cracked_stress]
  #   type = ComputeDuctileFracturedStress
  #   base_name = ductile_fractured
  #   c = c
  #   gc_factor = 0.5
  #   c_cri = 0.5
  #   F_name = E_el
  #   use_current_history_variable = true
  #   uncracked_base_name = uncracked
  #   finite_strain_model = true
  # [../]
  [./cracked_stress]
    type = ComputeCrackedStress
    c = c
    kdamage = 1e-5
    F_name = E_el
    use_current_history_variable = true
    uncracked_base_name = uncracked
    finite_strain_model = true
    base_name = ductile_fractured
  [../]
[]

[Postprocessors]
  [./av_stress_yy]
    type = ElementAverageValue
    variable = ductile_fractured_stress_yy
  [../]
  [./av_strain_yy]
    type = SideAverageValue
    variable = disp_y
    boundary = top
  [../]
  [./av_uncracked_stress_yy]
    type = ElementAverageValue
    variable = uncracked_stress_yy
  [../]
  [./max_c]
    type = ElementExtremeValue
    variable = c
  [../]
  # [./G0_pos]
  #   type = ElementAverageMaterialProperty
  #   mat_prop = phi_pos
  # [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient

  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_factor_mat_solving_package'
  petsc_options_value = 'lu superlu_dist'

  nl_rel_tol = 1e-8
  l_tol = 1e-4
  l_max_its = 100
  nl_max_its = 10
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.0002
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 8
  [../]
  # dt = 2.0e-3
  num_steps = 5
  # end_time = 0.005
[]

[Outputs]
  exodus = true
  csv = true
[]
