/** @file
*****************************************************************************
Declaration of SRS interfaces for ppzkSNARK proof system Plonk.

This includes:
- class usrs
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

// description of the circuit. contains only the number of gates for
// now.
template<typename ppT> struct circuit_t {
    // number of gates in the analyzed circuit, denoted by "n" in
    // [GWC19]
    size_t num_gates;
};

//
// A note on the distinction between an srs and a universal srs.
//
// A universal srs (usrs) is composed *only* of monomials i.e. encoded
// powers of the secret value in the group generator. Therefore a usrs
// is independent of any particular circuit.
//
// The (plain) srs is a specialization of the usrs for one particular
// circuit and is derived from the usrs e.g. as
//
// usrs = <encoded powers of secret>
// srs = (proving_key, verificataion_key) = derive(usrs, circuit_description)
//

// Universal srs (usrs). Contains secret encoded monomials with
// maximum degree MAX_DEGREE and is so independent of any particular
// circuit.
template<typename ppT> class usrs
{
public:
    /// Array of powers of secret \alpha, encoded in G1:
    /// [1]_1, [\alpha]_1, [\alpha^2]_1, ..., [\alpha^{n+2}]_1
    std::vector<libff::G1<ppT>> secret_powers_g1;

    /// Array of powers of secret \alpha, encoded in G2:
    /// [1]_2, [\alpha]_2
    std::vector<libff::G2<ppT>> secret_powers_g2;

    usrs(
        std::vector<libff::G1<ppT>> &&secret_powers_g1,
        std::vector<libff::G2<ppT>> &&secret_powers_g2)
        : secret_powers_g1(secret_powers_g1), secret_powers_g2(secret_powers_g2)
    {
    }
};

template<typename ppT>
usrs<ppT> plonk_usrs_derive_from_secret(const libff::Fr<ppT> secret);

// Plain (i.e. non-universal srs). Contains secret encoded monomials
// with maximum degree equal to the number of gates of the analyzed
// circuit + 2. Dependent on the circuit.
template<typename ppT> class srs
{
public:
    // description of the circuit. contains only the number of gates for
    // now. the highest degree of the encoded power monomials will
    // be num_gates+2.
    circuit_t<ppT> circuit;
    /// Array of powers of secret \alpha, encoded in G1:
    /// [1]_1, [\alpha]_1, [\alpha^2]_1, ..., [\alpha^{n+2}]_1
    std::vector<libff::G1<ppT>> secret_powers_g1;

    /// Array of powers of secret \alpha, encoded in G2:
    /// [1]_2, [\alpha]_2
    std::vector<libff::G2<ppT>> secret_powers_g2;

    srs(std::vector<libff::G1<ppT>> &&secret_powers_g1,
        std::vector<libff::G2<ppT>> &&secret_powers_g2)
        : secret_powers_g1(secret_powers_g1), secret_powers_g2(secret_powers_g2)
    {
    }

    srs(const circuit_t<ppT> &circuit) : circuit(circuit){};

    // derive from the usrs
    void derive(const usrs<ppT> usrs);
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
