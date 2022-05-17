/** @file
*****************************************************************************

Declaration of interfaces for ppzkSNARK proof system Plonk.

This includes:
- class for proving key
- class for verification key
- class for key pair (proving key & verification key)
- class for proof
- generator algorithm / setup
- prover algorithm
- verifier algorithm 

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

namespace libsnark
{

  enum W_polys_id{a, b, c};
  enum Q_polys_id{L, R, M, O, C};
  enum t_polys_id{lo, mid, hi};
  enum omega_id{base, base_k1, base_k2};

template<typename FieldT> using polynomial = std::vector<FieldT>;

/******************************** Proving key ********************************/
/**
 * A proving key for Plonk
 */
template<typename ppT> class plonk_proving_key
{
public:
  /// Array of powers of secret \alpha, encoded in G1:
  /// [1]_1, [\alpha]_1, [\alpha^2]_1, ..., [\alpha^{n+2}]_1
  std::vector<libff::G1<ppT>> secret_powers_g1;

  plonk_proving_key(){};
  plonk_proving_key(
		    std::vector<libff::G1<ppT>>&& secret_powers_g1)
    :secret_powers_g1(std::move(secret_powers_g1)){};
};

/******************************* Verification key ****************************/

/**
 * A verification key for Plonk
 */
template<typename ppT> class plonk_verification_key
{
public:
  /// Array of powers of secret \alpha, encoded in G2:
  /// [1]_2, [\alpha]_2
  std::vector<libff::G2<ppT>> secret_powers_g2;
  
  plonk_verification_key();
  plonk_verification_key(
		    std::vector<libff::G2<ppT>>&& secret_powers_g2)
    :secret_powers_g2(std::move(secret_powers_g2)){};
};

/********************************** Key pair *********************************/

/**
 * A key pair for Plonk, which consists of a proving key and a
 * verification key.
 */
template<typename ppT> class plonk_keypair
{
public:
    plonk_proving_key<ppT> pk;
    plonk_verification_key<ppT> vk;

    plonk_keypair() = default;
    plonk_keypair(const plonk_keypair<ppT> &other) =
        default;
    plonk_keypair(
        plonk_proving_key<ppT> &&pk,
        plonk_verification_key<ppT> &&vk)
        : pk(std::move(pk)), vk(std::move(vk))
    {
    }

    plonk_keypair(plonk_keypair<ppT> &&other) = default;
};
  
/*********************************** Proof ***********************************/

/**
 * A proof for the Plonk GG-ppzkSNARK.
 */
template<typename ppT> class plonk_proof
{
public:
  // proof elements

  plonk_proof() {};
};
  
/***************************** Main algorithms *******************************/

  // Generate the proving and verification key
  // \see r1cs_gg_ppzksnark_generator_from_secrets
  template<typename ppT>
  plonk_keypair<ppT> plonk_generator_from_secrets(
						  const libff::Fr<ppT> secret,
						  size_t num_gates
						  );

  
} // namespace libsnark

#include "libsnark/zk_proof_systems/plonk/plonk.tcc"

#endif // PLONK_PPZKSNARK_HPP_
