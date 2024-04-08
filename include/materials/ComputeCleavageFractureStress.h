#pragma once

#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "RotationTensor.h"
#include "DerivativeMaterialInterface.h"

/**
 * Computes energy and modifies the stress for phase field fracture. Can be used with any
 * constitutive model or elastic symmetry.
 */
class ComputeCleavageFractureStress : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  ComputeCleavageFractureStress(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  virtual void initQpStatefulProperties();
  /// Base name of the stress after being modified to include cracks
  const std::string _base_name;

  /// Base name of the uncracked stress and strain
  const std::string _uncracked_base_name;

  /// Mechanical_strain if finite_strain_model = false, otherwise elastic_strain
  const MaterialProperty<RankTwoTensor> & _plastic_strain;

  /// Uncracked stress calculated by another material
  const MaterialProperty<RankTwoTensor> & _uncracked_stress;

  /// Uncracked Jacobian_mult calculated by another material
  const MaterialProperty<RankFourTensor> & _uncracked_Jacobian_mult;
  /// young's modulus
  const Real _E;
  /// micro cracks length
  const Real _a;
  /// stress triaxiality
  const MaterialProperty<Real> & _neta;
  /// lode angle
  const MaterialProperty<Real> & _lode_angle;

  /// Variable defining the phase field damage parameter
  const VariableValue & _c;

  /// initial fracture toughness
  const Real _gc_ini;
  /// final fracture toughness
  const Real _gc_inf;
  /// fracture toughenss degradation factor
  const Real _alpha;
  /// Characteristic length, controls damage zone thickness
  const MaterialProperty<Real> & _l;

  // Viscosity, defining how quickly the crack propagates
  const MaterialProperty<Real> & _visco;

  /// Stress being computed by this kernel
  MaterialProperty<RankTwoTensor> & _stress;

  MaterialProperty<Real> & _F;
  MaterialProperty<Real> & _dFdc;
  MaterialProperty<Real> & _d2Fdc2;
  MaterialProperty<RankTwoTensor> & _d2Fdcdstrain;
  MaterialProperty<RankTwoTensor> & _dstress_dc;

  /// derivative of stress w.r.t. strain (_dstress_dstrain)
  MaterialProperty<RankFourTensor> & _Jacobian_mult;

  /// Property where the value for kappa will be defined
  MaterialProperty<Real> & _kappa;

  /// Property where the value for L will be defined
  MaterialProperty<Real> & _L;
};
