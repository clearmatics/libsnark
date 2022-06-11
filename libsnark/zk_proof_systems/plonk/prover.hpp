/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_PROVER_HPP_
#define LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_PROVER_HPP_

#include <libsnark/zk_proof_systems/plonk/srs.hpp>

/// Declaration of Prover interfaces for ppzkSNARK proof system Plonk.This
/// includes:
///
/// - class for proof
/// - class for prover
///
/// Reference:
/// - \[GWC19]:
///   Title: "Plonk: Permutations over lagrange-bases for oecumenical
///   noninteractive arguments of knowledge", Ariel Gabizon, Zachary
///   J. Williamson, and Oana Ciobotaru, Cryptology ePrint Archive,
///   Report 2019/953, 2019, <https://eprint.iacr.org/2019/953>

namespace libsnark
{

/// A proof for the Plonk GG-ppzkSNARK.
template<typename ppT> class plonk_proof
{
public:
    using Field = libff::Fr<ppT>;

    // Plonk proof Pi
    //
    // Pi ([a]_1, [b]_1, [c]_1, [z]_1,
    //     [t_lo]_1, [t_mi]_1, [t_hi]_1,
    //     \bar{a}, \bar{b}, \bar{c},
    //     \bar{S_sigma1}, \bar{S_sigma2}, \bar{z_w},
    //     [W_zeta]_1, [W_{zeta omega_roots}]_1
    //     r_zeta (*))

    // [a]_1, [b]_1, [c]_1
    std::vector<libff::G1<ppT>> W_polys_blinded_at_secret_g1;
    // [z]_1,
    libff::G1<ppT> z_poly_at_secret_g1;
    // [t_lo]_1, [t_mi]_1, [t_hi]_1
    std::vector<libff::G1<ppT>> t_poly_at_secret_g1;
    // \bar{a}, \bar{b}, \bar{c},
    Field a_zeta;
    Field b_zeta;
    Field c_zeta;
    // \bar{S_sigma1}, \bar{S_sigma2},
    Field S_0_zeta;
    Field S_1_zeta;
    // \bar{z_w}
    Field z_poly_xomega_zeta;
    // [W_zeta]_1, [W_{zeta omega_roots}]_1
    libff::G1<ppT> W_zeta_at_secret;
    libff::G1<ppT> W_zeta_omega_at_secret;
    // r_zeta (*)
    Field r_zeta;

    plonk_proof(
        std::vector<libff::G1<ppT>> &W_polys_blinded_at_secret_g1,
        libff::G1<ppT> &z_poly_at_secret_g1,
        std::vector<libff::G1<ppT>> &t_poly_at_secret_g1,
        Field &a_zeta,
        Field &b_zeta,
        Field &c_zeta,
        Field &S_0_zeta,
        Field &S_1_zeta,
        Field &z_poly_xomega_zeta,
        libff::G1<ppT> &W_zeta_at_secret,
        libff::G1<ppT> &W_zeta_omega_at_secret,
        Field &r_zeta)
        : W_polys_blinded_at_secret_g1(W_polys_blinded_at_secret_g1)
        , z_poly_at_secret_g1(z_poly_at_secret_g1)
        , t_poly_at_secret_g1(t_poly_at_secret_g1)
        , a_zeta(a_zeta)
        , b_zeta(b_zeta)
        , c_zeta(c_zeta)
        , S_0_zeta(S_0_zeta)
        , S_1_zeta(S_1_zeta)
        , z_poly_xomega_zeta(z_poly_xomega_zeta)
        , W_zeta_at_secret(W_zeta_at_secret)
        , W_zeta_omega_at_secret(W_zeta_omega_at_secret)
        , r_zeta(r_zeta)
    {
    }
};

// Prover round 0 output
template<typename ppT> struct round_zero_out_t {
    // - zh_poly: vanishing polynomial Zh (from round 0)
    std::vector<libff::Fr<ppT>> zh_poly;
    // - null_poly: 0 polynomial (from round 0)
    polynomial<libff::Fr<ppT>> null_poly;
    // - neg_one_poly: -1 polynomial (from round 0)
    polynomial<libff::Fr<ppT>> neg_one_poly;
};

// Prover round 1 output
template<typename ppT> struct round_one_out_t {
    // - blind_scalars: blinding scalars b1, b2, ..., b9 (only
    //   b1-b6 used in round 1)
    std::vector<libff::Fr<ppT>> blind_scalars;
    // - W_polys: witness polynomials (Lagrange interpolation of the
    //   witness values)
    std::vector<polynomial<libff::Fr<ppT>>> W_polys;
    // - W_polys_blinded: blinded witness polynomials
    std::vector<std::vector<libff::Fr<ppT>>> W_polys_blinded;
    // - W_polys_blinded_at_secret_g1: the blinded witness polynomials
    //   evaluated at the secret input denoted [a]_1, [b]_1, [c]_1 in
    //   [GWC19]
    std::vector<libff::G1<ppT>> W_polys_blinded_at_secret_g1;
};

// Prover round 2 output
template<typename ppT> struct round_two_out_t {
    // - beta: permutation challenge -- hash of transcript
    libff::Fr<ppT> beta;
    // - gamma: permutation challenge -- hash of transcript
    libff::Fr<ppT> gamma;
    // - z_poly: blinded accumulator poly z(x)
    polynomial<libff::Fr<ppT>> z_poly;
    // - z_poly_at_secret_g1: blinded accumulator poly z(x) evaluated at
    // secret
    libff::G1<ppT> z_poly_at_secret_g1;
};

// Prover round 3 output
template<typename ppT> struct round_three_out_t {
    // - alpha: quotient challenge -- hash of transcript
    libff::Fr<ppT> alpha;
    // - z_poly_xomega: the polynomial z(x*w) i.e. z(x) shifted by w
    std::vector<libff::Fr<ppT>> z_poly_xomega;
    // - t_poly: t(x) divided in three parts t(x) = t_lo(x) + t_mid(x) x^n
    //   + t_hi(x) x^{2n}
    std::vector<polynomial<libff::Fr<ppT>>> t_poly;
    // - t_poly_long: the quotient polynomial t(x) (see Round 3, pp28
    //   [GWC19])
    polynomial<libff::Fr<ppT>> t_poly_long;
    // - t_poly_at_secret_g1: t(x) evaluated at the secret input zeta
    //   i.e. t(zeta)
    std::vector<libff::G1<ppT>> t_poly_at_secret_g1;
};

// Prover round 4 output
template<typename ppT> struct round_four_out_t {
    // - zeta: evaluation challenge -- hash of transcript
    libff::Fr<ppT> zeta;
    // - a_zeta, b_zeta, c_zeta: the blinded witness polynomials a(x),
    //   b(x), c(x) (denoted by W_polys_blinded[] output from Round 1)
    //   evaluated at x=zeta i.e. a(z), b(z), c(z)
    libff::Fr<ppT> a_zeta;
    libff::Fr<ppT> b_zeta;
    libff::Fr<ppT> c_zeta;
    // - S_0_zeta, S_1_zeta: the permutation polynomials S_sigma_1(x),
    //   S_sigma_2(x) from the common preprocessed input (see [GWC19],
    //   Sect. 8.1) evaluated at x=zeta i.e. S_sigma_1(z), S_sigma_2(z)
    libff::Fr<ppT> S_0_zeta;
    libff::Fr<ppT> S_1_zeta;
    // - z_poly_xomega_zeta: the polynomial z(x*w) i.e. z(x) shifted by w
    //   (output from Round 3) evaluated at x=zeta i.e. z(zeta*w)
    libff::Fr<ppT> z_poly_xomega_zeta;
    // - t_zeta: the quotient polynomial t(x) output from Round 3, see
    //   pp28 [GWC19]) evaluated at x=zeta i.e. t(z). IMPORTANT! the
    //   original Plonk proposal [GWC19] does not output this parameter
    //   t_zeta. The Python reference implementation does, so we do the
    //   same in order to match the test vectors. TODO can remove t_zeta
    //   in the future
    libff::Fr<ppT> t_zeta;
};

// Prover round 5 output
template<typename ppT> struct round_five_out_t {
    // - nu: opening challenge -- hash of transcript (denoted by v in
    //   [GWC19])
    libff::Fr<ppT> nu;
    // - r_zeta: linearisation polynomial r(x) evaluated at x=zeta
    //   ie. r(zeta)
    libff::Fr<ppT> r_zeta;
    // - W_zeta_at_secret: commitment to opening proof polynomial
    //   W_zeta(x) at secert input i.e. [W_zeta(secret)]_1
    libff::G1<ppT> W_zeta_at_secret;
    // - W_zeta_omega_at_secret: commitment to opening proof polynomial
    //   W_{zeta omega}(x) at secert input i.e. [W_{zeta omega}(secret)]_1
    libff::G1<ppT> W_zeta_omega_at_secret;
    // - u: multipoint evaluation challenge -- hash of transcript
    libff::Fr<ppT> u;
};

/// Plonk prover. Computes object of class plonk_proof.
template<typename ppT> class plonk_prover
{
    using Field = libff::Fr<ppT>;

public:
    static round_zero_out_t<ppT> round_zero(const srs<ppT> srs);

    static round_one_out_t<ppT> round_one(
        const round_zero_out_t<ppT> &round_zero_out,
        const std::vector<libff::Fr<ppT>> &witness,
        const srs<ppT> &srs);

    static round_two_out_t<ppT> round_two(
        const round_zero_out_t<ppT> &round_zero_out,
        const round_one_out_t<ppT> &round_one_out,
        const std::vector<libff::Fr<ppT>> &witness,
        const srs<ppT> &srs);

    static round_three_out_t<ppT> round_three(
        const round_zero_out_t<ppT> &round_zero_out,
        const round_one_out_t<ppT> &round_one_out,
        const round_two_out_t<ppT> &round_two_out,
        const srs<ppT> &srs);

    static round_four_out_t<ppT> round_four(
        const round_one_out_t<ppT> &round_one_out,
        const round_three_out_t<ppT> &round_three_out,
        const srs<ppT> &srs);

    static round_five_out_t<ppT> round_five(
        const round_zero_out_t<ppT> &round_zero_out,
        const round_one_out_t<ppT> &round_one_out,
        const round_two_out_t<ppT> &round_two_out,
        const round_three_out_t<ppT> &round_three_out,
        const round_four_out_t<ppT> &round_four_out,
        const srs<ppT> &srs);

    static plonk_proof<ppT> compute_proof(const srs<ppT> &srs);
};

} // namespace libsnark

#include "libsnark/zk_proof_systems/plonk/prover.tcc"

#endif // LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_PROVER_HPP_
