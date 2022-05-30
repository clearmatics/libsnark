/** @file
*****************************************************************************

Declaration of SRS interfaces for ppzkSNARK proof system Plonk.

This includes:
- class srs
- class for proving key (TODO: not implemented)
- class for verification key (TODO: not implemented)
- class for key pair (proving key & verification key) (TODO: not implemented)

References:

\[GWC19]:
"Plonk: Permutations over lagrange-bases for oecumenical noninteractive
arguments of knowledge", Ariel Gabizon, Zachary J. Williamson, and Oana
Ciobotaru, Cryptology ePrint Archive, Report 2019/953, 2019,
<https://eprint.iacr.org/2019/953>

*****************************************************************************
* @author     This file is part of libsnark, developed by SCIPR Lab
*             and contributors (see AUTHORS).
* @copyright  MIT license (see LICENSE file)
*****************************************************************************/

#ifndef PLONK_PPZKSNARK_SRS_HPP_
#define PLONK_PPZKSNARK_SRS_HPP_

namespace libsnark
{
/********************************** SRS ***********************************/

/// The srs generated by the setup step. See kzg10.hpp .
template<typename ppT> class srs
{
public:
    /// Array of powers of secret \alpha, encoded in G1:
    /// [1]_1, [\alpha]_1, [\alpha^2]_1, ..., [\alpha^{n+2}]_1
    std::vector<libff::G1<ppT>> secret_powers_g1;

    /// Array of powers of secret \alpha, encoded in G2:
    /// [1]_2, [\alpha]_2
    std::vector<libff::G2<ppT>> secret_powers_g2;

    srs(){};

    srs(std::vector<libff::G1<ppT>> &&secret_powers_g1,
        std::vector<libff::G2<ppT>> &&secret_powers_g2);

    // for debug only
    void derive_from_secret(const libff::Fr<ppT> secret, size_t num_gates);
};

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
    plonk_proving_key(std::vector<libff::G1<ppT>> &&secret_powers_g1)
        : secret_powers_g1(std::move(secret_powers_g1)){};
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
    plonk_verification_key(std::vector<libff::G2<ppT>> &&secret_powers_g2)
        : secret_powers_g2(std::move(secret_powers_g2)){};
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
    plonk_keypair(const plonk_keypair<ppT> &other) = default;
    plonk_keypair(plonk_proving_key<ppT> &&pk, plonk_verification_key<ppT> &&vk)
        : pk(std::move(pk)), vk(std::move(vk))
    {
    }

    plonk_keypair(plonk_keypair<ppT> &&other) = default;
};

} // namespace libsnark

#include "libsnark/zk_proof_systems/plonk/srs.tcc"

#endif // PLONK_PPZKSNARK_SRS_HPP_
