#ifndef __PLONK_EXAMPLE_HPP__
#define __PLONK_EXAMPLE_HPP__

/*

Example Plonk circuit (TODO: give reference)

P(x) = x**3 + x + 5 = 3

circuit has 6 gates + 2 dummy gates (to make power of 2)

gates / constraints

        L   R    O
1 mul   x * x  = v1
2 mul  v1 * x  = v2
3 add  v2 + x  = v3
4 con5  1 * 5  = 5
5 pin   1 * 35 = 35
6 add  v3 + 5  = 35
7 dum   /    /    /
8 dum   /    /    /

wire polynomials

w_L = [ x, v1, v2,  1,  1, v3,  /,  /] = [a1, a2, a3, a4, a5, a6, a7, a8] = a
w_R = [ x,  x,  x,  5, 35,  5,  /,  /] = [b1, b2, b3, b4, b5, b6, b7, b8] = b
w_O = [v1, v2, v3,  5, 35, 35,  /,  /] = [c1, c2, c3, c4, c5, c6, c7, c8] = c

wires = [a1, a2, a3, a4, a5, a6, a7, a8, b1, b2, b3, b4, b5, b6, b7, b8, c1, c2, c3, c4, c5, c6, c7, c8]
index = [ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
perm  = [ 9, 17, 18,  5,  4, 19,  7,  8, 10, 11,  1, 14, 21, 20, 15, 16,  2,  3,  6, 12, 22, 13, 23, 24]

witness

x = 3 => v1 = 9, v2 = 27, v3 = 30

w_L = a = [ 3,  9, 27,  1,  1, 30,  0,  0]
w_R = b = [ 3,  3,  3,  5, 35,  5,  0,  0]
w_O = c = [ 9, 27, 30,  5, 35, 35,  0,  0]

W = w_L + w_R + w_O

encoding of plonk gates

add  = [1, 1,-1, 0, 0]
mul  = [0, 0,-1, 1, 0]
con5 = [0, 1, 0, 0, 5]
pi   = [0, 1, 0, 0, 0]
dum  = [0, 0, 0, 0, 0]

gates matrix

(q_L * a) + (q_R * b) + (q_O * c) + (q_M * a * b) + (q_C) = 0

gates                  q_L q_R q_O q_M q_C
1 mul   x * x  = v1 : [  0,  0, -1,  1,  0]
2 mul  v1 * x  = v2 : [  0,  0, -1,  1,  0]
3 add  v2 + x  = v3 : [  1,  1, -1,  0,  0]
4 con5  1 * 5  = 5  : [  0,  1,  0,  0, -5]
5 pin   1 * 35 = 35 : [  0,  1,  0,  0,  0]
6 add  v2 + 5  = 35 : [  1,  1, -1,  0,  0]
7 dum   /    /    / : [  0,  0,  0,  0,  0]
8 dum   /    /    / : [  0,  0,  0,  0,  0]

q_L = [ 0,  0,  1,  0,  0,  1,  0,  0]
q_R = [ 0,  0,  1,  1,  1,  1,  0,  0]
q_O = [-1, -1, -1,  0,  0,  0,  0,  0]
q_M = [ 1,  1,  0,  0,  0, -1,  0,  0]
q_C = [ 0,  0,  0, -1,  0,  0,  0,  0]

*/
  
namespace libsnark
{

template<typename ppT> class plonk_example
{
public:
  using Field = libff::Fr<ppT>;
  using BaseField = libff::Fq<ppT>;
  template<typename FieldT> using polynomial = std::vector<FieldT>;

  // Circuit data
  
  // number of gates / constraints. we have 6 gates for the example
  // circuit + 2 dummy gates to make it a power of 2 (for the fft)
  size_t num_gates;
  
  // number of q-polynomials
  size_t num_qpolys;
  
  // hard-coded gates matrix for the example circuit
  // P(x) = x**3 + x + 5 = 3
  // Each column is a q-vector
  std::vector<std::vector<Field>>gates_matrix;
  
  // Transposed gates matrix: each row is a q-vector WARN: rows 2
  // q_O and 3 q_M are swapped ti match the Plonk_Py test vectors
  // implementation (reason unclear)
  std::vector<std::vector<Field>>gates_matrix_transpose;
  
  // witness values
  // w_L = a = [ 3,  9, 27,  1,  1, 30,  0,  0]
  // w_R = b = [ 3,  3,  3,  5, 35,  5,  0,  0]
  // w_O = c = [ 9, 27, 30,  5, 35, 35,  0,  0]
  // W = w_L + w_R + w_O
  std::vector<Field> witness;

  // wire permutation (TODO: add function plonk_compute_permutation())
  std::vector<size_t> wire_permutation;

  // public input (PI)
  Field public_input;
  // index of the row of the PI in the non-transposed gates_matrix 
  size_t public_input_index;

  // n-th root of unity omega in Fq (n=8 is the number of constraints
  // in the example). omega is a generator of the multiplicative
  // subgroup H.  Example (2**32)-th primitive root of unity in the
  // base field Fq of bls12-381 i.e. such that omega_base**(2**32) =
  // 1. The bls12-381 prime q is such that any power of 2 divides
  // (q-1). In particular 2**32|(q-1) and so the 2**32-th root of
  // unity exists.
  Field omega_base;
  
  // Constants k1,k2 to generate domains on which to evaluate the witness
  // polynomials. k can be random, but we fix it for debug to match
  // against the test vector values
  Field k1;
  // Similarly, k2 can be random, but we fix it to match the test
  // vectors
  Field k2;

  // H_gen contains the generators of H, k1 H and K2 H in one place
  // ie. omega, omega_k1 and omega_k2
  std::vector<Field> H_gen;
  // H_gen permuted according to the wire permutation
  std::vector<Field> H_gen_permute;  
  // random hidden element secret (toxic waste). we fix it to a
  // constant in order to match against the test vectors
  Field secret; 
  // powers of secret times G1: 1*G1, secret^1*G1, secret^2*G1, ...
  std::vector<std::vector<BaseField>> secret_powers_g1;  
  // powers of secret times G2: 1*G2, secret^1*G2
  std::vector<std::vector<BaseField>> secret_powers_g2;  
  // blinding scalars b1, b2, ..., b9. random but fixed to match
  // the python test vectors
  std::vector<Field> prover_blind_scalars;
  // Hashes of transcript (Fiat-Shamir heuristic)
  Field beta;
  Field gamma;
  Field alpha;
  Field zeta;
  Field nu; // v in the paper
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
  // Verifier Step 12: batch validate all evaluations via pairing
  std::vector<BaseField> pairing_first_lhs;
  std::vector<BaseField> pairing_first_rhs;
  //  std::vector<BaseField> pairing_second_lhs;
  //  std::vector<BaseField> pairing_second_rhs;

  plonk_example();
};

} // namespace libsnark
  
#include <libsnark/zk_proof_systems/plonk/tests/plonk_example.cpp>

#endif // __PLONK_EXAMPLE_HPP__
