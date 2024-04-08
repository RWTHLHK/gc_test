#pragma once

#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "RotationTensor.h"
#include "DerivativeMaterialInterface.h"
#include "ComputeCrackedStress.h"
/**
 * Computes energy and modifies the stress for phase field fracture. Can be used with any
 * constitutive model or elastic symmetry.
 */
class ComputeDuctileFracturedStress : public ComputeCrackedStress
{
public:
  static InputParameters validParams();

  ComputeDuctileFracturedStress(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  virtual void initQpStatefulProperties();
  // strain
  const MaterialProperty<RankTwoTensor> & _uncracked_stress_old;
  MaterialProperty<Real> & _G0_pos;
  const MaterialProperty<Real> & _G0_pos_old;
  const Real _c_cri;
  const Real _gc_factor;
  const MaterialProperty<RankTwoTensor> & _plastic_strain;
  MaterialProperty<RankTwoTensor> & _stress0pos;
  const MaterialProperty<RankTwoTensor> & _stress0pos_old;
  const MaterialProperty<RankTwoTensor> & _strain_increment;
};
