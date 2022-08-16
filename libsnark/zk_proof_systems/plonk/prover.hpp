/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_PROVER_HPP_
#define LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_PROVER_HPP_

#include "libsnark/zk_proof_systems/plonk/srs.hpp"

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

    /// Plonk proof Pi
    ///
    /// Pi ([a]_1, [b]_1, [c]_1, [z]_1,
    ///     [t_lo]_1, [t_mi]_1, [t_hi]_1,
    ///     \bar{a}, \bar{b}, \bar{c},
    ///     \bar{S_sigma1}, \bar{S_sigma2}, \bar{z_w},
    ///     [W_zeta]_1, [W_{zeta omega_roots}]_1
    ///     r_zeta (*))

    /// [a]_1, [b]_1, [c]_1
    std::vector<libff::G1<ppT>> W_polys_blinded_at_secret_g1;

    /// [z]_1,
    libff::G1<ppT> z_poly_at_secret_g1;

    /// [t_lo]_1, [t_mi]_1, [t_hi]_1
    std::vector<libff::G1<ppT>> t_poly_at_secret_g1;

    /// \bar{a}, \bar{b}, \bar{c},
    Field a_zeta;
    Field b_zeta;
    Field c_zeta;

    /// \bar{S_sigma1}, \bar{S_sigma2},
    Field S_0_zeta;
    Field S_1_zeta;

    /// \bar{z_w}
    Field z_poly_xomega_zeta;

    /// [W_zeta]_1, [W_{zeta omega_roots}]_1
    libff::G1<ppT> W_zeta_at_secret;
    libff::G1<ppT> W_zeta_omega_at_secret;

    /// r_zeta (*)
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
        Field &r_zeta);
};

/// Prover round 0 output
template<typename ppT> struct round_zero_out_t {

    /// - zh_poly: vanishing polynomial Zh (from round 0)
    std::vector<libff::Fr<ppT>> zh_poly;

    /// - null_poly: 0 polynomial (from round 0)
    polynomial<libff::Fr<ppT>> null_poly;

    /// - neg_one_poly: -1 polynomial (from round 0)
    polynomial<libff::Fr<ppT>> neg_one_poly;

    /// stuct constructor
    round_zero_out_t(
        const std::vector<libff::Fr<ppT>> &&zh_poly,
        const polynomial<libff::Fr<ppT>> &&null_poly,
        const polynomial<libff::Fr<ppT>> &&neg_one_poly);
};

/// Prover round 1 output
template<typename ppT> struct round_one_out_t {

    /// - W_polys: witness polynomials (Lagrange interpolation of the
    ///   witness values)
    std::vector<polynomial<libff::Fr<ppT>>> W_polys;

    /// - W_polys_blinded: blinded witness polynomials
    std::vector<std::vector<libff::Fr<ppT>>> W_polys_blinded;

    /// - W_polys_blinded_at_secret_g1: the blinded witness polynomials
    ///   evaluated at the secret input denoted [a]_1, [b]_1, [c]_1 in
    ///   [GWC19]
    std::vector<libff::G1<ppT>> W_polys_blinded_at_secret_g1;

    /// stuct constructor
    round_one_out_t(
        const std::vector<polynomial<libff::Fr<ppT>>> &&W_polys,
        const std::vector<std::vector<libff::Fr<ppT>>> &&W_polys_blinded,
        const std::vector<libff::G1<ppT>> &&W_polys_blinded_at_secret_g1);
};

/// Prover round 2 output
template<typename ppT> struct round_two_out_t {

    /// - z_poly: blinded accumulator poly z(x)
    polynomial<libff::Fr<ppT>> z_poly;

    /// - z_poly_at_secret_g1: blinded accumulator poly z(x) evaluated at
    /// secret
    libff::G1<ppT> z_poly_at_secret_g1;

    /// stuct constructor
    round_two_out_t(
        polynomial<libff::Fr<ppT>> &&z_poly,
        libff::G1<ppT> &&z_poly_at_secret_g1);
};

/// Prover round 3 output
template<typename ppT> struct round_three_out_t {

    /// - z_poly_xomega: the polynomial z(x*w) i.e. z(x) shifted by w
    std::vector<libff::Fr<ppT>> z_poly_xomega;

    /// - t_poly: t(x) divided in three parts t(x) = t_lo(x) + t_mid(x) x^n
    ///   + t_hi(x) x^{2n}
    std::vector<polynomial<libff::Fr<ppT>>> t_poly;

    /// - t_poly_long: the quotient polynomial t(x) (see Round 3, pp28
    ///   [GWC19])
    polynomial<libff::Fr<ppT>> t_poly_long;

    /// - t_poly_at_secret_g1: t(x) evaluated at the secret input zeta
    ///   i.e. t(zeta)
    std::vector<libff::G1<ppT>> t_poly_at_secret_g1;

    /// stuct constructor
    round_three_out_t(
        std::vector<libff::Fr<ppT>> &&z_poly_xomega,
        std::vector<polynomial<libff::Fr<ppT>>> &&t_poly,
        polynomial<libff::Fr<ppT>> &&t_poly_long,
        std::vector<libff::G1<ppT>> &&t_poly_at_secret_g1);
};

/// Prover round 4 output
template<typename ppT> struct round_four_out_t {

    /// - a_zeta, b_zeta, c_zeta: the blinded witness polynomials a(x),
    ///   b(x), c(x) (denoted by W_polys_blinded[] output from Round 1)
    ///   evaluated at x=zeta i.e. a(z), b(z), c(z)
    libff::Fr<ppT> a_zeta;
    libff::Fr<ppT> b_zeta;
    libff::Fr<ppT> c_zeta;

    /// - S_0_zeta, S_1_zeta: the permutation polynomials S_sigma_1(x),
    ///   S_sigma_2(x) from the common preprocessed input (see [GWC19],
    ///   Sect. 8.1) evaluated at x=zeta i.e. S_sigma_1(z), S_sigma_2(z)
    libff::Fr<ppT> S_0_zeta;
    libff::Fr<ppT> S_1_zeta;

    /// - z_poly_xomega_zeta: the polynomial z(x*w) i.e. z(x) shifted by w
    ///   (output from Round 3) evaluated at x=zeta i.e. z(zeta*w)
    libff::Fr<ppT> z_poly_xomega_zeta;

    /// - t_zeta: the quotient polynomial t(x) output from Round 3, see
    ///   pp28 [GWC19]) evaluated at x=zeta i.e. t(z). IMPORTANT! the
    ///   original Plonk proposal [GWC19] does not output this parameter
    ///   t_zeta. The Python reference implementation does, so we do the
    ///   same in order to match the test vectors. TODO can remove t_zeta
    ///   in the future
    libff::Fr<ppT> t_zeta;
    /// stuct constructor
    round_four_out_t(
        libff::Fr<ppT> &&a_zeta,
        libff::Fr<ppT> &&b_zeta,
        libff::Fr<ppT> &&c_zeta,
        libff::Fr<ppT> &&S_0_zeta,
        libff::Fr<ppT> &&S_1_zeta,
        libff::Fr<ppT> &&z_poly_xomega_zeta,
        libff::Fr<ppT> &&t_zeta);
};

/// Prover round 5 output
template<typename ppT> struct round_five_out_t {

    /// - r_zeta: linearisation polynomial r(x) evaluated at x=zeta
    ///   ie. r(zeta)
    libff::Fr<ppT> r_zeta;

    /// - W_zeta_at_secret: commitment to opening proof polynomial
    ///   W_zeta(x) at secert input i.e. [W_zeta(secret)]_1
    libff::G1<ppT> W_zeta_at_secret;

    /// - W_zeta_omega_at_secret: commitment to opening proof polynomial
    ///   W_{zeta omega}(x) at secert input i.e. [W_{zeta omega}(secret)]_1
    libff::G1<ppT> W_zeta_omega_at_secret;

    /// struct constructor
    round_five_out_t(
        libff::Fr<ppT> &&r_zeta,
        libff::G1<ppT> &&W_zeta_at_secret,
        libff::G1<ppT> &&W_zeta_omega_at_secret);
};

/// Plonk prover. Computes object of class plonk_proof.
template<typename ppT> class plonk_prover
{
    using Field = libff::Fr<ppT>;

public:
    /// Prover Round 0 initialization
    ///
    /// Initialization
    ///
    /// INPUT
    /// \param[in] srs: structured reference string containing also
    /// circuit-specific
    ///   information
    ///
    /// OUTPUT
    /// \param[out] W_polys: Lagrange interpolation of the witness values
    /// \param[out] zh_poly: vanishing polynomial
    /// \param[out] null_poly: 0 polynomial
    /// \param[out] neg_one_poly: -1 polynomial
    static round_zero_out_t<ppT> round_zero(const srs<ppT> &srs);

    /// Prover Round 1
    ///
    /// INPUT
    /// \param[in] zh_poly: vanishing polynomial Zh (from round 0)
    /// \param[in] null_poly: 0 polynomial (from round 0)
    /// \param[in] neg_one_poly: -1 polynomial (from round 0)
    /// \param[in] blind_scalars: blinding scalars b1, b2, ..., b9 (only
    ///            b1-b6 used in round 1)
    /// \param[in] witness: witness values
    /// \param[in] srs: structured reference string containing also
    ///            circuit-specific information
    ///
    /// OUTPUT
    /// \param[out] W_polys: witness polynomials (Lagrange interpolation
    ///             of the witness values)
    /// \param[out] W_polys_blinded: blinded witness polynomials
    /// \param[out] W_polys_blinded_at_secret_g1: the blinded witness
    ///             polynomials evaluated at the secret input denoted
    ///             [a]_1, [b]_1, [c]_1 in [GWC19]
    /// \param[out] transcript_hasher: accumulates the communication
    ///             transcript into a buffer to be hashed after prover
    ///             rounds 1,2,3,4,5 (cf. fiat-shamir heuristic).
    static round_one_out_t<ppT> round_one(
        const round_zero_out_t<ppT> &round_zero_out,
        const std::vector<libff::Fr<ppT>> &blind_scalars,
        const std::vector<libff::Fr<ppT>> &witness,
        const srs<ppT> &srs,
        transcript_hasher<ppT> &hasher);

    /// Prover Round 2
    ///
    /// INPUT
    /// \param[in] zh_poly: vanishing polynomial Zh (from round 0)
    /// \param[in] blind_scalars: blinding scalars b1, b2, ..., b9 (only
    ///            b7,b8,b9 used in round 2)
    /// \param[in] witness: witness values
    /// \param[in] srs: structured reference string containing also
    ///            circuit-specific information
    ///
    /// OUTPUT
    /// \param[out] z_poly: blinded accumulator poly z(x)
    /// \param[out] z_poly_at_secret_g1: blinded accumulator poly z(x)
    ///             evaluated at secret
    /// \param[out] transcript_hasher: accumulates the communication
    ///             transcript into a buffer to be hashed after prover
    ///             rounds 1,2,3,4,5 (cf. fiat-shamir heuristic).
    static round_two_out_t<ppT> round_two(
        const libff::Fr<ppT> &beta,
        const libff::Fr<ppT> &gamma,
        const round_zero_out_t<ppT> &round_zero_out,
        const std::vector<libff::Fr<ppT>> blind_scalars,
        const std::vector<libff::Fr<ppT>> &witness,
        const srs<ppT> &srs,
        transcript_hasher<ppT> &hasher);

    /// Prover Round 3
    ///
    /// INPUT
    /// \param[in] zh_poly: vanishing polynomial Zh (from Round 0)
    /// \param[in] W_polys_blinded: blinded witness polynomials (from
    ///            Round 1)
    /// \param[in] beta, gamma: permutation challenges -- hashes of
    ///            transcript (from round 2)
    /// \param[in] z_poly: blinded accumulator poly z(x) (from Round 2)
    /// \param[in] srs: structured reference string containing also
    ///            circuit-specific information
    ///
    /// OUTPUT
    /// \param[out] t_poly_long: the quotient polynomial t(x) (see Round
    ///             3, pp28 [GWC19])
    /// \param[out] t_poly: t(x) divided in three parts t(x) = t_lo(x) +
    ///             t_mid(x) x^n + t_hi(x) x^{2n}
    /// \param[out] t_poly_at_secret_g1: t(x) evaluated at the secret
    ///             input zeta i.e. t(zeta)
    /// \param[out] z_poly_xomega: the polynomial z(x*w) i.e. z(x) shifted
    ///             by w
    /// \param[out] transcript_hasher: accumulates the communication
    ///             transcript into a buffer to be hashed after prover
    ///             rounds 1,2,3,4,5 (cf. fiat-shamir heuristic).
    static round_three_out_t<ppT> round_three(
        const libff::Fr<ppT> &alpha,
        const libff::Fr<ppT> &beta,
        const libff::Fr<ppT> &gamma,
        const round_zero_out_t<ppT> &round_zero_out,
        const round_one_out_t<ppT> &round_one_out,
        const round_two_out_t<ppT> &round_two_out,
        const srs<ppT> &srs,
        transcript_hasher<ppT> &hasher);

    /// Prover Round 4
    ///
    /// INPUT
    /// \param[in] W_polys_blinded: blinded witness polynomials (from
    ///            Round 1)
    /// \param[in] z_poly_xomega: the polynomial z(x*w) i.e. z(x) shifted
    ///            by w (from Round 3)
    /// \param[in] t_poly_long: the quotient polynomial t(x) (see Round 3,
    ///            pp28 [GWC19]) (from Round 3)
    /// \param[in] srs: structured reference string containing also
    ///            circuit-specific information
    ///
    /// OUTPUT
    /// \param[out] a_zeta, b_zeta, c_zeta: the blinded witness
    ///             polynomials a(x), b(x), c(x) (denoted by
    ///             W_polys_blinded[] output from Round 1) evaluated at
    ///             x=zeta i.e. a(z), b(z), c(z)
    /// \param[out] S_0_zeta, S_1_zeta: the permutation polynomials
    ///             S_sigma_1(x), S_sigma_2(x) from the common
    ///             preprocessed input (see [GWC19], Sect. 8.1) evaluated
    ///             at x=zeta i.e. S_sigma_1(z), S_sigma_2(z)
    /// \param[out] z_poly_xomega_zeta: the polynomial z(x*w) i.e. z(x)
    ///             shifted by w (output from Round 3) evaluated at x=zeta
    ///             i.e. z(zeta*w)
    /// \param[out] t_zeta: the quotient polynomial t(x) output from Round
    ///             3, see pp28 [GWC19]) evaluated at x=zeta
    ///             i.e. t(z). IMPORTANT! the original Plonk proposal
    ///             [GWC19] does not output this parameter t_zeta. The
    ///             Python reference implementation does, so we do the
    ///             same in order to match the test vectors. TODO can
    ///             remove t_zeta in the future
    /// \param[out] transcript_hasher: accumulates the communication
    ///             transcript into a buffer to be hashed after prover
    ///             rounds 1,2,3,4,5 (cf. fiat-shamir heuristic).
    static round_four_out_t<ppT> round_four(
        const libff::Fr<ppT> &zeta,
        const round_one_out_t<ppT> &round_one_out,
        const round_three_out_t<ppT> &round_three_out,
        const srs<ppT> &srs,
        transcript_hasher<ppT> &hasher);

    /// Prover Round 5
    ///
    /// INPUT
    /// \param[in] beta, gamma: permutation challenges -- hashes of
    ///            transcript (from round 2)
    /// \param[in] alpha: quotient challenge -- hash of transcript (from
    ///            round 3)
    /// \param[in] zeta: evaluation challenge -- hash of transcript (from
    ///            round 4)
    /// \param[in] a_zeta, b_zeta, c_zeta: the blinded witness polynomials
    ///            a(x), b(x), c(x) (denoted by W_polys_blinded[] output
    ///            from Round 1) evaluated at x=zeta i.e. a(z), b(z), c(z)
    ///            (from round 4)
    /// \param[in] S_0_zeta, S_1_zeta: the permutation polynomials
    ///            S_sigma_1(x), S_sigma_2(x) from the common preprocessed
    ///            input (see [GWC19], Sect. 8.1) evaluated at x=zeta
    ///            i.e. S_sigma_1(z), S_sigma_2(z) (from round 4)
    /// \param[in] t_zeta: the quotient polynomial t(x) output from Round
    ///            3, see pp28 [GWC19]) evaluated at x=zeta
    ///            i.e. t(z). IMPORTANT! the original Plonk proposal
    ///            [GWC19] does not output this parameter t_zeta. The
    ///            Python reference implementation does, so we do the same
    ///            in order to match the test vectors. TODO can remove
    ///            t_zeta in the future (from round 4)
    /// \param[in] z_poly_xomega_zeta: the polynomial z(x*w) i.e. z(x)
    ///            shifted by w (output from Round 3) evaluated at x=zeta
    ///            i.e. z(zeta*w) (from round 4)
    /// \param[in] W_polys_blinded: blinded witness polynomials (from
    ///            round 1)
    /// \param[in] t_poly: t(x) divided in three parts t(x) = t_lo(x) +
    ///            t_mid(x) x^n + t_hi(x) x^{2n} (from round 3)
    /// \param[in] z_poly: blinded accumulator poly z(x) (from round 2)
    /// \param[in] srs: structured reference string containing also
    ///            circuit-specific information
    ///
    /// OUTPUT
    /// \param[out] r_zeta: linearisation polynomial r(x) evaluated at
    ///             x=zeta ie. r(zeta)
    /// \param[out] W_zeta_at_secret: commitment to opening proof
    ///             polynomial W_zeta(x) at secert input
    ///             i.e. [W_zeta(secret)]_1
    /// \param[out] W_zeta_omega_at_secret: commitment to opening proof
    ///             polynomial W_{zeta omega}(x) at secert input
    ///             i.e. [W_{zeta omega}(secret)]_1
    /// \param[out] transcript_hasher: accumulates the communication
    ///             transcript into a buffer to be hashed after prover
    ///             rounds 1,2,3,4,5 (cf. fiat-shamir heuristic).
    static round_five_out_t<ppT> round_five(
        const libff::Fr<ppT> &alpha,
        const libff::Fr<ppT> &beta,
        const libff::Fr<ppT> &gamma,
        const libff::Fr<ppT> &zeta,
        const libff::Fr<ppT> &nu,
        const round_zero_out_t<ppT> &round_zero_out,
        const round_one_out_t<ppT> &round_one_out,
        const round_two_out_t<ppT> &round_two_out,
        const round_three_out_t<ppT> &round_three_out,
        const round_four_out_t<ppT> &round_four_out,
        const srs<ppT> &srs,
        transcript_hasher<ppT> &hasher);

    /// Prover compute SNARK proof
    ///
    /// Pi ([a]_1, [b]_1, [c]_1, [z]_1,
    ///     [t_lo]_1, [t_mi]_1, [t_hi]_1,
    ///     \bar{a}, \bar{b}, \bar{c},
    ///     \bar{S_sigma1}, \bar{S_sigma2}, \bar{z_w},
    ///     [W_zeta]_1, [W_{zeta omega}]_1
    ///     r_zeta)
    ///
    /// \note in the reference Python implementation, r_zeta (the
    /// evaluation of the linearlization polynomial r(X) at zeta from
    /// Prover round 5) is added to the pi-SNARK proof. In the paper this
    /// is omitted, which seems to make the proof shorter by 1 element at
    /// the epxense of a slightly heavier computation on the verifier's
    /// side. Here we follow the reference implementation to make sure we
    /// match the test values. TODO: once all test vectors are verified,
    /// we may remove r_zeta from the proof to be fully compliant with the
    /// paper.
    ///
    /// Mapping code-to-paper quantities
    ///
    /// \param W_polys_blinded_at_secret_g1[a, b, c]: [a]_1, [b]_1, [c]_1
    ///        (from Round 1)
    /// \param z_poly_at_secret_g1: [z]_1 (from Round 2)
    /// \param t_poly_at_secret_g1[lo, mi, hi]: [t_lo]_1, [t_mi]_1,
    ///        [t_hi]_1 (from Round 3)
    /// \param a_zeta, b_zeta, c_zeta, S_0_zeta, S_1_zeta,
    ///        z_poly_xomega_zeta: \bar{a}, \bar{b}, \bar{c},
    ///        \bar{S_sigma1}, \bar{S_sigma2}, \bar{z_w} (from Round 4)
    /// \param W_zeta_at_secret, W_zeta_omega_at_secret: [W_zeta]_1,
    ///        [W_{zeta omega}]_1 (from Round 5)
    ///
    /// INPUT
    /// \param[in] srs: structured reference string containing also
    ///            circuit-specific information
    /// \param[in] witness: all internal values and public input
    ///            corresponding to the given circuit
    /// \param[in] blind_scalars: random blinding scalars b1, b2, ..., b9
    ///            used in prover rounds 1 and 2 (see Sect. 8.3, roumds
    ///            1,2 [GWC19])
    /// \param[in] transcript_hasher: hashes of the communication
    ///            transcript after prover rounds 1,2,3,4,5.
    ///
    /// OUTPUT
    /// \param[out] proof: SNARK proof Pi (see above)
    static plonk_proof<ppT> compute_proof(
        const srs<ppT> &srs,
        const std::vector<Field> &witness,
        const std::vector<libff::Fr<ppT>> &blind_scalars,
        transcript_hasher<ppT> &hasher);
};

} // namespace libsnark

#include "libsnark/zk_proof_systems/plonk/prover.tcc"

#endif // LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_PROVER_HPP_
