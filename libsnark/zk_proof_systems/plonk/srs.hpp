/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_SRS_HPP_
#define LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_SRS_HPP_

/// Declaration of SRS interfaces for ppzkSNARK proof system Plonk. This
/// includes:
///
/// - class usrs
/// - class srs
/// - class for proving key (TODO: not implemented)
/// - class for verification key (TODO: not implemented)
/// - class for key pair (proving key & verification key) (TODO: not
///   implemented)
///
/// Reference:
/// - \[GWC19]:
///   Title: "Plonk: Permutations over lagrange-bases for oecumenical
///   noninteractive arguments of knowledge", Ariel Gabizon, Zachary
///   J. Williamson, and Oana Ciobotaru, Cryptology ePrint Archive,
///   Report 2019/953, 2019, <https://eprint.iacr.org/2019/953>

#include "libsnark/common/data_structures/polynomial.hpp"

namespace libsnark
{

///
/// A note on the distinction between an srs and a universal srs.
///
/// A universal srs (usrs) is composed *only* of monomials i.e. encoded
/// powers of the secret value in the group generator. Therefore a usrs
/// is independent of any particular circuit.
///
/// The (plain) srs is a specialization of the usrs for one particular
/// circuit and is derived from the usrs e.g. as
///
/// usrs = <encoded powers of secret>
/// srs = (proving_key, verificataion_key) = derive(usrs, circuit_description)
///

/// Universal srs (usrs). Contains secret encoded monomials with some
/// maximum degree and is so independent of any particular circuit.
template<typename ppT> class usrs
{
public:
    /// Array of powers of secret \secret, encoded in G1:
    /// [1]_1, [\secret]_1, [\secret^2]_1, ..., [\secret^{n+2}]_1
    std::vector<libff::G1<ppT>> secret_powers_g1;

    /// Array of powers of secret \secret, encoded in G2:
    /// [1]_2, [\secret]_2
    std::vector<libff::G2<ppT>> secret_powers_g2;

    usrs(
        std::vector<libff::G1<ppT>> &&secret_powers_g1,
        std::vector<libff::G2<ppT>> &&secret_powers_g2)
        : secret_powers_g1(secret_powers_g1), secret_powers_g2(secret_powers_g2)
    {
    }
};

template<typename ppT>
usrs<ppT> plonk_usrs_derive_from_secret(
    const libff::Fr<ppT> &secret, const size_t max_degree);

/// Plain (i.e. non-universal srs). Contains secret encoded monomials
/// with maximum degree equal to the number of gates of the analyzed
/// circuit + 2, plus circuit description. Dependent on the circuit.
template<typename ppT> class srs
{
public:
    using Field = libff::Fr<ppT>;

    /// number of gates in the analyzed circuit, denoted by "n" in
    /// [GWC19]
    size_t num_gates;

    /// number of selector polynomials (q-polynomials) (= 5 in the
    /// vanilla Plonk proposal [GWC19])
    size_t num_qpolys;

    /// Lagrange basis
    std::vector<polynomial<Field>> L_basis;

    /// Public input polynomial
    polynomial<Field> PI_poly;

    /// Circuit selector polynomials (Q-polynomials)
    std::vector<polynomial<Field>> Q_polys;

    /// Permutation polynomials S_sigma_1, S_sigma_2, S_sigma_2 (see
    /// [GWC19], Sect. 8.1)
    std::vector<polynomial<Field>> S_polys;

    /// omega[0] are the n roots of unity, omega[1] are omega[0]*k1,
    /// omega[2] are omega[0]*k2
    std::vector<std::vector<Field>> omega_roots;

    /// H_gen contains the generators of H, k1 H and k2 H in one place
    /// ie. omega, omega_k1 and omega_k2
    std::vector<Field> H_gen;

    /// H_gen permuted according to the wire permutation
    std::vector<Field> H_gen_permute;

    /// constants for H, k1 H, k2 H
    libff::Fr<ppT> k1;
    libff::Fr<ppT> k2;

    /// Array of powers of secret \alpha, encoded in G1: [1]_1,
    /// [\alpha]_1, [\alpha^2]_1, ..., [\alpha^{n+2}]_1
    std::vector<libff::G1<ppT>> secret_powers_g1;

    /// Array of powers of secret \alpha, encoded in G2: [1]_2,
    /// [\alpha]_2
    std::vector<libff::G2<ppT>> secret_powers_g2;

    srs(const size_t &num_gates,
        const size_t &num_qpolys,
        const std::vector<polynomial<Field>> &L_basis,
        const polynomial<Field> &PI_poly,
        const std::vector<polynomial<Field>> &Q_polys,
        const std::vector<polynomial<Field>> &S_polys,
        const std::vector<std::vector<Field>> &omega_roots,
        const std::vector<Field> &H_gen,
        const std::vector<Field> &H_gen_permute,
        const libff::Fr<ppT> &k1,
        const libff::Fr<ppT> &k2,
        std::vector<libff::G1<ppT>> &&secret_powers_g1,
        std::vector<libff::G2<ppT>> &&secret_powers_g2)
        : num_gates(num_gates)
        , num_qpolys(num_qpolys)
        , L_basis(L_basis)
        , PI_poly(PI_poly)
        , Q_polys(Q_polys)
        , S_polys(S_polys)
        , omega_roots(omega_roots)
        , H_gen(H_gen)
        , H_gen_permute(H_gen_permute)
        , k1(k1)
        , k2(k2)
        , secret_powers_g1(secret_powers_g1)
        , secret_powers_g2(secret_powers_g2)
    {
    }
};

template<typename ppT>
srs<ppT> plonk_srs_derive_from_usrs(
    const usrs<ppT> &usrs, const circuit_t<ppT> &circuit);

/// A proving key for Plonk
template<typename ppT> class plonk_proving_key
{
public:
    /// Array of powers of secret \alpha, encoded in G1: [1]_1,
    /// [\alpha]_1, [\alpha^2]_1, ..., [\alpha^{n+2}]_1
    std::vector<libff::G1<ppT>> secret_powers_g1;

    plonk_proving_key(){};
    plonk_proving_key(std::vector<libff::G1<ppT>> &&secret_powers_g1)
        : secret_powers_g1(std::move(secret_powers_g1)){};
};

/// A verification key for Plonk
template<typename ppT> class plonk_verification_key
{
public:
    /// Array of powers of secret \alpha, encoded in G2: [1]_2,
    /// [\alpha]_2
    std::vector<libff::G2<ppT>> secret_powers_g2;

    plonk_verification_key();
    plonk_verification_key(std::vector<libff::G2<ppT>> &&secret_powers_g2)
        : secret_powers_g2(std::move(secret_powers_g2)){};
};

/// A key pair for Plonk, which consists of a proving key and a
/// verification key.
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

/// Hashes of transcript after prover rounds 1,2,3,4,5
template<typename ppT> struct transcript_hash_t {

    /// - beta: permutation challenge - hashes of transcript after round 1
    libff::Fr<ppT> beta;
    /// - gamma: permutation challenge - hashes of transcript after round 1
    libff::Fr<ppT> gamma;
    /// - alpha: quotient challenge - hash of transcript after rounds 1,2
    libff::Fr<ppT> alpha;
    /// - zeta: evaluation challenge - hash of transcriptafter rounds 1,2,3
    libff::Fr<ppT> zeta;
    /// - nu: opening challenge - hash of transcript after rounds 1,2,3,4
    /// (denoted by v in [GWC19])
    libff::Fr<ppT> nu;
    /// - u: multipoint evaluation challenge -- hash of transcript after rounds
    /// 1,2,3,4,5
    libff::Fr<ppT> u;

    /// stuct constructor
    transcript_hash_t(
        libff::Fr<ppT> &beta,
        libff::Fr<ppT> &gamma,
        libff::Fr<ppT> &alpha,
        libff::Fr<ppT> &zeta,
        libff::Fr<ppT> &nu,
        libff::Fr<ppT> &u)
        : beta(beta), gamma(gamma), alpha(alpha), zeta(zeta), nu(nu), u(u)
    {
    }
};

} // namespace libsnark

#include "libsnark/zk_proof_systems/plonk/srs.tcc"

#endif // LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_SRS_HPP_