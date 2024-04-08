

#include "ComputeBrittleFracturedStress.h"

registerMooseObject("gc_testApp", ComputeBrittleFracturedStress);

InputParameters
ComputeBrittleFracturedStress::validParams()
{
  InputParameters params = ComputeCrackedStress::validParams();
  params.addClassDescription("Computes energy and modifies the stress for phase field dual fracture");
  params.addRequiredCoupledVar("c_duc", "Order parameter for ductile damage");
  return params;
}

ComputeBrittleFracturedStress::ComputeBrittleFracturedStress(const InputParameters & parameters)
  : ComputeCrackedStress(parameters),
    _c_duc(coupledValue("c_duc"))
{
}

void
ComputeBrittleFracturedStress::initQpStatefulProperties()
{
  ComputeCrackedStress::initQpStatefulProperties();
}

void
ComputeBrittleFracturedStress::computeQpProperties()
{
  const Real c = _c[_qp];
  const Real c_duc = _c_duc[_qp];
  // Zero out values when c > 1
  Real cfactor = 1.0;
  if (c > 1.0)
    cfactor = 0.0;

  // Create the positive and negative projection tensors
  RankFourTensor I4sym(RankFourTensor::initIdentitySymmetricFour);
  std::vector<Real> eigval;
  RankTwoTensor eigvec;
  RankFourTensor Ppos = _uncracked_stress[_qp].positiveProjectionEigenDecomposition(eigval, eigvec);
  RankFourTensor Pneg = I4sym - Ppos;

  // Project the positive and negative stresses
  RankTwoTensor stress0pos = Ppos * _uncracked_stress[_qp];
  RankTwoTensor stress0neg = Pneg * _uncracked_stress[_qp];

  // Compute the positive and negative elastic energies
  Real G0_pos = (stress0pos).doubleContraction(_strain[_qp]) / 2.0;
  Real G0_neg = (stress0neg).doubleContraction(_strain[_qp]) / 2.0;

  // Update the history variable
  if (G0_pos > _hist_old[_qp])
    _hist[_qp] = G0_pos;
  else
    _hist[_qp] = _hist_old[_qp];

  Real hist_variable = _hist_old[_qp];
  if (_use_current_hist)
    hist_variable = _hist[_qp];

  // Compute degredation function and derivatives
  Real h = cfactor * (1.0 - c) * (1.0 - c) * (1.0 - _kdamage) + _kdamage;
  Real dhdc = -2.0 * cfactor * (1.0 - c) * (1.0 - _kdamage);
  Real d2hdc2 = 2.0 * cfactor * (1.0 - _kdamage);

  // Compute stress and its derivatives
  _stress[_qp] = (Ppos * h + Pneg) * _uncracked_stress[_qp];
  _dstress_dc[_qp] = stress0pos * dhdc;
  _Jacobian_mult[_qp] = (Ppos * h + Pneg) * _uncracked_Jacobian_mult[_qp];

  // Compute energy and its derivatives
  _F[_qp] = hist_variable * h/ (1 - c_duc) / (1 - c_duc)- G0_neg + _gc_prop[_qp] * c * c / (2 * _l[_qp]);
  _dFdc[_qp] = hist_variable * dhdc / (1 - c_duc) / (1 - c_duc)+ _gc_prop[_qp] * c / _l[_qp];
  _d2Fdc2[_qp] = hist_variable * d2hdc2 / (1 - c_duc)/ (1 - c_duc) + _gc_prop[_qp] / _l[_qp];

  // 2nd derivative wrt c and strain = 0.0 if we used the previous step's history varible
  if (_use_current_hist)
    _d2Fdcdstrain[_qp] = stress0pos * dhdc;

  // Assign L and kappa
  _kappa[_qp] = _gc_prop[_qp] * _l[_qp];
  _L[_qp] = 1.0 / (_gc_prop[_qp] * _visco[_qp]);
}
