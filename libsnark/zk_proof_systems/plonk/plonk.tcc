/** @file
*****************************************************************************

Implementation of interfaces for a ppzkSNARK for Plonk.

See plonk.hpp .

*****************************************************************************
* @author     This file is part of libsnark, developed by SCIPR Lab
*             and contributors (see AUTHORS).
* @copyright  MIT license (see LICENSE file)
*****************************************************************************/

#ifndef PLONK_PPZKSNARK_TCC_
#define PLONK_PPZKSNARK_TCC_

#include "libsnark/zk_proof_systems/plonk/plonk.hpp"

namespace libsnark
{

  template<typename ppT>
  srs<ppT>::srs(
	   std::vector<libff::G1<ppT>> &&secret_powers_g1,
	   std::vector<libff::G2<ppT>> &&secret_powers_g2)
    : secret_powers_g1(secret_powers_g1), secret_powers_g2(secret_powers_g2)
  {
  }

  // Generate SRS \see r1cs_gg_ppzksnark_generator_from_secrets, \see
  // kzg10<ppT>::setup_from_secret(
  template<typename ppT>
  srs<ppT> plonk_setup_from_secret(
				   const libff::Fr<ppT> secret,
				   size_t num_gates
				   )
  {
    // compute powers of secret times G1: 1*G1, secret^1*G1, secret^2*G1, ...
    const libff::bigint<libff::Fr<ppT>::num_limbs> secret_bigint = secret.as_bigint();
    const size_t window_size =
      std::max(libff::wnaf_opt_window_size<libff::G1<ppT>>(secret_bigint.num_bits()), 1ul);
    const std::vector<long> naf =
      libff::find_wnaf<libff::Fr<ppT>::num_limbs>(window_size, secret_bigint);    
    std::vector<libff::G1<ppT>> secret_powers_g1;
    secret_powers_g1.reserve(num_gates + 3);
    libff::G1<ppT> secret_i_g1 = libff::G1<ppT>::one();
    secret_powers_g1.push_back(secret_i_g1);
    for (size_t i = 1; i < (num_gates + 3); ++i) {
      // secret^i * G1
      secret_i_g1 = libff::fixed_window_wnaf_exp<libff::G1<ppT>>(
								window_size, secret_i_g1, naf);
      secret_powers_g1.push_back(secret_i_g1);
    }
    //    plonk_proving_key<ppT> pk = plonk_proving_key<ppT>(std::move(secret_powers_g1));

    // compute powers of secret times G2: 1*G2, secret^1*G2
    std::vector<libff::G2<ppT>> secret_powers_g2;
    secret_powers_g2.reserve(2);
    // secret^0 * G2 = G2
    libff::G2<ppT> secret_0_g2 = libff::G2<ppT>::one();
    secret_powers_g2.push_back(secret_0_g2);
    // secret^1 * G2
    libff::G2<ppT> secret_1_g2 = secret * libff::G2<ppT>::one();
    secret_powers_g2.push_back(secret_1_g2);
    //    plonk_verification_key<ppT> vk = plonk_verification_key<ppT>(std::move(secret_powers_g2));
    
    srs<ppT> srs(std::move(secret_powers_g1), std::move(secret_powers_g2));
    return srs;
  }
  
} // namespace libsnark

#endif // PLONK_PPZKSNARK_TCC_
