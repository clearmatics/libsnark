/** @file
*****************************************************************************
Implementation of SRS interfaces for a ppzkSNARK for Plonk.

See srs.hpp .

*****************************************************************************
* @author     This file is part of libsnark, developed by SCIPR Lab
*             and contributors (see AUTHORS).
* @copyright  MIT license (see LICENSE file)
*****************************************************************************/

#ifndef PLONK_PPZKSNARK_SRS_TCC_
#define PLONK_PPZKSNARK_SRS_TCC_

namespace libsnark
{

// Compute a universal srs (usrs). It is composed *only* of encoded
// powers of the secret value in the group generator. Therefore a usrs
// is independent of any particular circuit.
//
// Note: the \ref derive_from_secret method is only for debug
template<typename ppT>
void usrs<ppT>::derive_from_secret(const libff::Fr<ppT> secret)
{
    // compute powers of secret times G1: 1*G1, secret^1*G1, secret^2*G1, ...
    const libff::bigint<libff::Fr<ppT>::num_limbs> secret_bigint =
        secret.as_bigint();
    const size_t window_size = std::max(
        libff::wnaf_opt_window_size<libff::G1<ppT>>(secret_bigint.num_bits()),
        1ul);
    const std::vector<long> naf =
        libff::find_wnaf<libff::Fr<ppT>::num_limbs>(window_size, secret_bigint);
    this->secret_powers_g1.reserve(MAX_DEGREE);
    libff::G1<ppT> secret_i_g1 = libff::G1<ppT>::one();
    this->secret_powers_g1.push_back(secret_i_g1);
    for (size_t i = 1; i < MAX_DEGREE; ++i) {
        // secret^i * G1
        secret_i_g1 = libff::fixed_window_wnaf_exp<libff::G1<ppT>>(
            window_size, secret_i_g1, naf);
        this->secret_powers_g1.push_back(secret_i_g1);
    }

    // compute powers of secret times G2: 1*G2, secret^1*G2
    // Note: in Plonk we *always* have 2 encoded elemnts in G2
    this->secret_powers_g2.reserve(2);
    // secret^0 * G2 = G2
    libff::G2<ppT> secret_0_g2 = libff::G2<ppT>::one();
    this->secret_powers_g2.push_back(secret_0_g2);
    // secret^1 * G2
    libff::G2<ppT> secret_1_g2 = secret * libff::G2<ppT>::one();
    this->secret_powers_g2.push_back(secret_1_g2);
}

// Generate (plain) SRS \see r1cs_gg_ppzksnark_generator_from_secrets,
// \see kzg10<ppT>::setup_from_secret
//
// The (plain) srs is a specialization of the usrs for one particular
// circuit and is derived from the usrs e.g.
//
// usrs = <encoded powers of secret>
// srs = (proving_key, verificataion_key) = derive(usrs, circuit_description)
//
template<typename ppT> void srs<ppT>::derive(const usrs<ppT> usrs)
{
    assert(this->circuit.num_gates <= MAX_DEGREE);
    // secret^i * G1
    this->secret_powers_g1.reserve(this->circuit.num_gates + 3);
    for (size_t i = 0; i < (this->circuit.num_gates + 3); ++i) {
        this->secret_powers_g1.push_back(usrs.secret_powers_g1[i]);
    }
    // secret^i * G2
    this->secret_powers_g2.reserve(2);
    for (size_t i = 0; i < 2; ++i) {
        this->secret_powers_g2.push_back(usrs.secret_powers_g2[i]);
    }
}

} // namespace libsnark

#endif // PLONK_PPZKSNARK_SRS_TCC_
