#ifndef __PLONK_EXAMPLE_HPP__
#define __PLONK_EXAMPLE_HPP__

namespace libsnark
{

template<typename ppT> class plonk_example
{
public:
  using Field = libff::Fr<ppT>;
  using BaseField = libff::Fq<ppT>;
  template<typename FieldT> using polynomial = std::vector<FieldT>;

  // Hashes of transcript (Fiat-Shamir heuristic)
  Field beta;
  Field gamma;
  Field alpha;
  Field zeta;
  Field nu; // v
  Field u;   
  // Prover Round 3
  std::vector<std::vector<BaseField>> t_poly_at_secret_g1;
  // Prover Round 4
  Field a_zeta;
  Field b_zeta;
  Field c_zeta;
  Field S_0_zeta;
  Field S_1_zeta;
  Field t_zeta;
  Field z_poly_xomega_zeta;
  // Prover Round 5
  Field r_zeta;
  polynomial<Field> W_zeta;
  polynomial<Field> W_zeta_omega;
  // point on the curve as a pair of X,Y coordinates (values in the
  // base field)
  std::vector<BaseField> W_zeta_at_secret;
  std::vector<BaseField> W_zeta_omega_at_secret;
  // Verifier precomputation
  std::vector<std::vector<BaseField>> Q_polys_at_secret_g1;
  std::vector<std::vector<BaseField>> S_polys_at_secret_g1;
  // Verifier Step 5: vanishing polynomial evaluation at zeta
  Field zh_zeta;
  // Verifier Step 6: compute Lagrange polynomial evaluation L1(zeta)
  Field L_0_zeta;
  // Verifier Step7: evaluate public input polynomial at zeta
  Field PI_zeta;
  // Verifier Step 8: compute quotient polynomial evaluation r'(zeta) = r(zeta) - r0, where r0 is a constant term
  Field r_prime_zeta;
  // Verifier Step 9
  std::vector<BaseField> D1;
  // Verifier Step 10: compute full batched polynomial commitment
  std::vector<BaseField> F1;
  // Verifier Step 11: compute group-encoded batch evaluation [E]_1
  std::vector<BaseField> E1;

  void initialize();
};

} // namespace libsnark
  
#include <libsnark/zk_proof_systems/plonk/tests/plonk_example.cpp>

#endif // __PLONK_EXAMPLE_HPP__
