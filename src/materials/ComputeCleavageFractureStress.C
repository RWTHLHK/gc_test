#include "ComputeCleavageFractureStress.h"
#include "RankTwoScalarTools.h"

registerMooseObject("gc_testApp", ComputeCleavageFractureStress);

InputParameters
ComputeCleavageFractureStress::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Computes energy and modifies the stress for phase field cleavage fracture");
  params.addRequiredCoupledVar("c", "Order parameter for damage");
  params.addParam<MaterialPropertyName>(
      "F_name", "Griffth", "Name of material property storing the local fracture driving energy");
  params.addParam<MaterialPropertyName>(
      "kappa_name",
      "kappa_op",
      "Name of material property being created to store the interfacial parameter kappa");
  params.addParam<MaterialPropertyName>(
      "mobility_name", "L", "Name of material property being created to store the mobility L");
  params.addParam<std::string>("base_name", "The base name used to save the cracked stress");
  params.addRequiredParam<std::string>("uncracked_base_name",
                                       "The base name used to calculate the original stress");
  params.addRequiredParam<Real>("critical_crack_length", "critical length of micro cracks");
  params.addRequiredParam<Real>("youngs_modulus",
                                "youngs_modulus for local fracture driving energy");
  params.addParam<MaterialPropertyName>(
      "stress_triaxiality", "neta", "stress traxiality for local fracture driving energy");
  params.addParam<MaterialPropertyName>(
      "lode_angle", "theta", "lode angle for local fracture driving energy");
  params.addRequiredParam<Real>("gc_ini", "initial fracture toughness");
  params.addRequiredParam<Real>("gc_inf", "final fracture toughness");
  params.addRequiredParam<Real>("fracture_toughness_degradation_factor",
                                "degradation factor due to plastic strain");
  return params;
}

ComputeCleavageFractureStress::ComputeCleavageFractureStress(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _uncracked_base_name(getParam<std::string>("uncracked_base_name") + "_"),
    _plastic_strain(
        getMaterialPropertyByName<RankTwoTensor>(_uncracked_base_name + "plastic_strain")),
    _uncracked_stress(getMaterialPropertyByName<RankTwoTensor>(_uncracked_base_name + "stress")),
    _uncracked_Jacobian_mult(
        getMaterialPropertyByName<RankFourTensor>(_uncracked_base_name + "Jacobian_mult")),
    _E(getParam<Real>("youngs_modulus")),
    _a(getParam<Real>("critical_crack_length")),
    _neta(getMaterialPropertyByName<Real>(getParam<MaterialPropertyName>("stress_triaxiality"))),
    _lode_angle(getMaterialPropertyByName<Real>(getParam<MaterialPropertyName>("lode_angle"))),
    _c(coupledValue("c")),
    _gc_ini(getParam<Real>("gc_ini")),
    _gc_inf(getParam<Real>("gc_inf")),
    _alpha(getParam<Real>("fracture_toughness_degradation_factor")),
    _l(getMaterialProperty<Real>("l")),
    _visco(getMaterialProperty<Real>("visco")),
    _stress(declareProperty<RankTwoTensor>(_base_name + "stress")),
    _F(declareProperty<Real>(getParam<MaterialPropertyName>("F_name"))),
    _dFdc(declarePropertyDerivative<Real>(getParam<MaterialPropertyName>("F_name"),
                                          coupledName("c", 0))),
    _d2Fdc2(declarePropertyDerivative<Real>(
        getParam<MaterialPropertyName>("F_name"), coupledName("c", 0), coupledName("c", 0))),
    _d2Fdcdstrain(declareProperty<RankTwoTensor>("d2Fdcdstrain")),
    _dstress_dc(declarePropertyDerivative<RankTwoTensor>("stress", coupledName("c", 0))),
    _Jacobian_mult(declareProperty<RankFourTensor>(_base_name + "Jacobian_mult")),
    _kappa(declareProperty<Real>(getParam<MaterialPropertyName>("kappa_name"))),
    _L(declareProperty<Real>(getParam<MaterialPropertyName>("mobility_name")))
{
}

void
ComputeCleavageFractureStress::initQpStatefulProperties()
{
  _stress[_qp].zero();
}

void
ComputeCleavageFractureStress::computeQpProperties()
{
  const Real c = _c[_qp];

  // Zero out values when c > 1
  Real cfactor = 1.0;
  if (c > 1.0)
    cfactor = 0.0;

  // compute equivalent stress
  Real mises_stress = RankTwoScalarTools::vonMisesStress(_uncracked_stress[_qp]);
  // compute maximum principal stress
  Real sigma1 = (_neta[_qp] + 2.0 / 3.0 * std::cos(_lode_angle[_qp])) * mises_stress;
  if (sigma1 < 0)
  {
    sigma1 = 0;
  }
  // compute effective strain energy(Griffth)
  Real phi = sigma1 * sigma1 * M_1_PI * _a / _E;
  // Create the positive and negative projection tensors
  // RankFourTensor I4sym(RankFourTensor::initIdentitySymmetricFour);
  // std::vector<Real> eigval;
  // RankTwoTensor eigvec;
  // RankFourTensor Ppos = _uncracked_stress[_qp].positiveProjectionEigenDecomposition(eigval,
  // eigvec); RankFourTensor Pneg = I4sym - Ppos;

  // // Project the positive and negative stresses
  // RankTwoTensor stress0pos = Ppos * _uncracked_stress[_qp];
  // RankTwoTensor stress0neg = Pneg * _uncracked_stress[_qp];

  // // Compute the positive and negative elastic energies
  // Real G0_pos = (stress0pos).doubleContraction(_strain[_qp]) / 2.0;
  // Real G0_neg = (stress0neg).doubleContraction(_strain[_qp]) / 2.0;

  // Compute degredation function and derivatives
  Real h = cfactor * (1.0 - c) * (1.0 - c);
  Real dhdc = -2.0 * cfactor * (1.0 - c);
  Real d2hdc2 = 2.0 * cfactor;

  // Compute stress and its derivatives
  _stress[_qp] = (1.0 - c) * _uncracked_stress[_qp];
  _dstress_dc[_qp] = -1.0 * _uncracked_stress[_qp];
  _Jacobian_mult[_qp] = h * _uncracked_Jacobian_mult[_qp];
  // compute equivalent plastic strain
  Real strain_eq = RankTwoScalarTools::effectiveStrain(_plastic_strain[_qp]);
  // compute degraded gc_prop
  Real gc_prop = _gc_ini + (_gc_inf - _gc_ini) * std::exp(-_alpha * strain_eq);
  // Compute energy and its derivatives
  _F[_qp] = phi * h + gc_prop * c * c / (2 * _l[_qp]);
  _dFdc[_qp] = phi * dhdc + gc_prop * c / _l[_qp];
  _d2Fdc2[_qp] = phi * d2hdc2 + gc_prop / _l[_qp];
  _d2Fdcdstrain[_qp] = phi * dhdc;
  // Assign L and kappa
  _kappa[_qp] = gc_prop * _l[_qp];
  _L[_qp] = 1.0 / (gc_prop * _visco[_qp]);
}
