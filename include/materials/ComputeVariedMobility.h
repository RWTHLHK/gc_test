#pragma once
#include "Material.h"
#include "RankTwoTensor.h"

class ComputeVariedMobility : public Material
{
public:
  static InputParameters validParams();

  ComputeVariedMobility(const InputParameters & parameters);

protected:
virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;
  // member property holds base name
  const std::string _base_name;
  // init mobility
  const Real _L_ini;
  // varied mobility
  MaterialProperty<Real> &_L;
  // strain
  const MaterialProperty<RankTwoTensor> &_strain;
};
