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
class ComputeBrittleFracturedStress : public ComputeCrackedStress
{
public:
  static InputParameters validParams();

  ComputeBrittleFracturedStress(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  virtual void initQpStatefulProperties();
  const VariableValue &_c_duc;
};
