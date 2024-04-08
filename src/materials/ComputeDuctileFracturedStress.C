

#include "ComputeDuctileFracturedStress.h"

registerMooseObject("gc_testApp", ComputeDuctileFracturedStress);

InputParameters
ComputeDuctileFracturedStress::validParams()
{
  InputParameters params = ComputeCrackedStress::validParams();
  params.addClassDescription(
      "Computes energy and modifies the stress for phase field dual fracture(ductile part)");
  params.addRequiredParam<Real>("c_cri",
                                "the maximum critical value that ductile damage can reach");
  params.addRequiredParam<Real>("gc_factor",
                                "ratio of gc_prop of ductile damage compared with total gc_prop");
  return params;
}

ComputeDuctileFracturedStress::ComputeDuctileFracturedStress(const InputParameters & parameters)
  : ComputeCrackedStress(parameters),
    _uncracked_stress_old(getMaterialPropertyOld<RankTwoTensor>(_uncracked_base_name + "stress")),
    _G0_pos(declareProperty<Real>("phi_pos")),
    _G0_pos_old(getMaterialPropertyOld<Real>("phi_pos")),
    _c_cri(getParam<Real>("c_cri")),
    _gc_factor(getParam<Real>("gc_factor")),
    _plastic_strain(getMaterialProperty<RankTwoTensor>(_uncracked_base_name + "plastic_strain")),
    _stress0pos(declareProperty<RankTwoTensor>(_uncracked_base_name + "stress0pos")),
    _stress0pos_old(getMaterialPropertyOld<RankTwoTensor>(_uncracked_base_name + "stress0pos")),
    _strain_increment(
        getMaterialProperty<RankTwoTensor>((_uncracked_base_name + "strain_increment")))
{
}

void
ComputeDuctileFracturedStress::initQpStatefulProperties()
{
  ComputeCrackedStress::initQpStatefulProperties();
  _G0_pos[_qp] = 0.0;
  _stress0pos[_qp].zero();
}

void
ComputeDuctileFracturedStress::computeQpProperties()
{
  const Real c = _c[_qp];
  // Zero out values when c > 1
  Real cfactor = 1.0;
  if (c > _c_cri)
    cfactor = 0.0;

  // Create the positive and negative projection tensors
  RankFourTensor I4sym(RankFourTensor::initIdentitySymmetricFour);
  std::vector<Real> eigval;
  RankTwoTensor eigvec;
  RankFourTensor Ppos = _uncracked_stress[_qp].positiveProjectionEigenDecomposition(eigval, eigvec);
  RankFourTensor Pneg = I4sym - Ppos;
  // Project the positive and negative stresses
  _stress0pos[_qp] = Ppos * _uncracked_stress[_qp];
  RankTwoTensor stress0neg = Pneg * _uncracked_stress[_qp];
  // Compute the positive strain energy
  _G0_pos[_qp] =
      _G0_pos_old[_qp] +
      MetaPhysicL::raw_value(_stress0pos[_qp])
              .doubleContraction(MetaPhysicL::raw_value(_strain_increment[_qp])) /
          2.0 +
      _stress0pos_old[_qp].doubleContraction(MetaPhysicL::raw_value(_strain_increment[_qp])) / 2.0;

  // Update the history variable
  if (_G0_pos[_qp] > _hist_old[_qp])
    _hist[_qp] = _G0_pos[_qp];
  else
    _hist[_qp] = _hist_old[_qp];

  Real hist_variable = _hist_old[_qp];
  if (_use_current_hist)
    hist_variable = _hist[_qp];

  // Compute degredation function and derivatives
  Real h =
      1.0 / _c_cri / _c_cri * cfactor * (_c_cri - c) * (_c_cri - c) * (1.0 - _kdamage) + _kdamage;
  Real dhdc = -2.0 * cfactor / _c_cri / _c_cri * (_c_cri - c) * (1.0 - _kdamage);
  Real d2hdc2 = 2.0 * cfactor * (1.0 - _kdamage) / _c_cri / _c_cri;

  // Compute stress and its derivatives
  _stress[_qp] = (Ppos * h + Pneg) * _uncracked_stress[_qp];
  _dstress_dc[_qp] = _stress0pos[_qp] * dhdc;
  _Jacobian_mult[_qp] = (Ppos * h + Pneg) * _uncracked_Jacobian_mult[_qp];

  // Compute energy and its derivatives
  _F[_qp] = hist_variable * h + _gc_prop[_qp] * c * c / (2 * _l[_qp]);
  _dFdc[_qp] = hist_variable * dhdc + _gc_prop[_qp] * c / _l[_qp];
  _d2Fdc2[_qp] = hist_variable * d2hdc2 + _gc_prop[_qp] / _l[_qp];

  // 2nd derivative wrt c and strain = 0.0 if we used the previous step's history varible
  if (_use_current_hist)
    _d2Fdcdstrain[_qp] = _stress0pos[_qp] * dhdc;

  // Assign L and kappa
  _kappa[_qp] = _gc_factor * _gc_prop[_qp] * _l[_qp];
  // only after plastic deformation ductile damage can occur
  // the gc_prop of ductile damage is lower than brittle one due to decohesion of inclusion-matrix
  // interface
  Real strain_eq =
      std::sqrt(2.0 / 3.0 * _plastic_strain[_qp].doubleContraction(_plastic_strain[_qp]));
  if (strain_eq > 0.0)
  {
    _L[_qp] = 1.0 / (_gc_factor * _gc_prop[_qp] * _visco[_qp]);
  }
  else
  {
    _L[_qp] = 0.0;
  }
}
