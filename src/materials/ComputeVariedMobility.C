#include "ComputeVariedMobility.h"
#include "RankTwoTensor.h"
#include "RankTwoScalarTools.h"
registerMooseObject("gc_testApp", ComputeVariedMobility);

InputParameters
ComputeVariedMobility::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Computes varied mobility");
  params.addRequiredCoupledVar("c", "Order parameter for damage");

  params.addParam<std::string>("base_name", "The base name used to save the cracked stress");
  params.addRequiredParam<Real>("init_mobility", "value of init mobility");
  return params;
}

void
ComputeVariedMobility::initQpStatefulProperties()
{
  _L[_qp] = 0.0;
}

ComputeVariedMobility::ComputeVariedMobility(const InputParameters & parameters)
  : Material(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _L_ini(getParam<Real>("init_mobility")),
    _L(declareProperty<Real>("ductile_mobility")),
    _strain(getMaterialProperty<RankTwoTensor>((_base_name + "plastic_strain")))
{
}

void
ComputeVariedMobility::computeQpProperties()
{
  // compute effective plastic strain
  Real strain_eq = std::sqrt(2.0 / 3.0 * _strain[_qp].doubleContraction(_strain[_qp]));
  if (strain_eq > 0.0)
  {
    _L[_qp] = _L_ini;
  }
}
