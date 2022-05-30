/** @file
*****************************************************************************

Declaration of interfaces for ppzkSNARK proof system Plonk.

This includes:
- class for common preprocessed input

The implementation instantiates the protocol of PlonK \[GWC19],

References:

\[GWC19]: 
"Plonk: Permutations over lagrange-bases for oecumenical noninteractive arguments of knowledge",
Ariel Gabizon, Zachary J. Williamson, and Oana Ciobotaru,
Cryptology ePrint Archive, Report 2019/953, 2019,
<https://eprint.iacr.org/2019/953>

*****************************************************************************
* @author     This file is part of libsnark, developed by SCIPR Lab
*             and contributors (see AUTHORS).
* @copyright  MIT license (see LICENSE file)
*****************************************************************************/

#ifndef PLONK_PPZKSNARK_HPP_
#define PLONK_PPZKSNARK_HPP_

#include <libff/algebra/curves/public_params.hpp>
#include <memory>
#include <libsnark/zk_proof_systems/plonk/tests/example.hpp>

// enable debug checks. in particular enable comparisons to test
// vector values.
#define DEBUG

// number of generators for H, Hk1, Hk2
const size_t NUM_HGEN = 3;

namespace libsnark
{

  enum W_polys_id{a, b, c};
  enum Q_polys_id{L, R, M, O, C};
  enum t_polys_id{lo, mid, hi};
  enum omega_id{base, base_k1, base_k2};

template<typename FieldT> using polynomial = std::vector<FieldT>;

/************************ Common Preprocessed Input ***********************/

/**
 * Plonk common preprocessed input
 */
template<typename ppT> class common_preprocessed_input
{
public:
  using Field = libff::Fr<ppT>;
  // number of gates / constraints
  size_t num_gates;
  // number of selector polynomials (q-polynomials) (= 5 in the
  // vanilla Plonk proposal [GWC19])
  size_t num_qpolys;
  // Lagrange basis
  std::vector<polynomial<Field>> L_basis;
  // Public input polynomial
  polynomial<Field> PI_poly;
  // Circuit selector polynomials (Q-polynomials)
  std::vector<polynomial<Field>> Q_polys;
  // Permutation polynomials S_sigma_1, S_sigma_2, S_sigma_2 (see
  // [GWC19], Sect. 8.1)
  std::vector<polynomial<Field>> S_polys;
  // omega[0] are the n roots of unity, omega[1] are omega[0]*k1,
  // omega[2] are omega[0]*k2
  std::vector<std::vector<Field>> omega_roots;
  // H_gen contains the generators of H, k1 H and k2 H in one place
  // ie. omega, omega_k1 and omega_k2
  std::vector<Field> H_gen;
  // H_gen permuted according to the wire permutation
  std::vector<Field> H_gen_permute;
  // constants for H, k1 H, k2 H
  libff::Fr<ppT> k1;
  libff::Fr<ppT> k2;

  common_preprocessed_input(){};

  // for debug
  void setup_from_example(plonk_example<ppT> example);
  
};  

/***************************** Main algorithms *******************************/


  
} // namespace libsnark

#include "libsnark/zk_proof_systems/plonk/plonk.tcc"

#endif // PLONK_PPZKSNARK_HPP_
