/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_VERIFIER_HPP_
#define LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_VERIFIER_HPP_

/// Declaration of Verifier interfaces for ppzkSNARK proof system Plonk. This
/// includes:
///
/// - class for verifier
///
/// References:
///
/// - \[GWC19]:
///   Title: "Plonk: Permutations over lagrange-bases for oecumenical
///   noninteractive arguments of knowledge", Ariel Gabizon, Zachary
///   J. Williamson, and Oana Ciobotaru, Cryptology ePrint Archive,
///   Report 2019/953, 2019, <https://eprint.iacr.org/2019/953>

namespace libsnark
{

/// SNARK proof
///
/// (\see plonk_prover::compute_proof)
///
/// Pi ([a]_1, [b]_1, [c]_1, [z]_1,
///     [t_lo]_1, [t_mi]_1, [t_hi]_1,
///     \bar{a}, \bar{b}, \bar{c},
///     \bar{S_sigma1}, \bar{S_sigma2}, \bar{z_w},
///     [W_zeta]_1, [W_{zeta omega}]_1
///     r_zeta (*))
///
/// Mapping code-to-paper quantities (format code : paper)
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
/// Verifier preprocessed input
///
/// \param Q_polys_at_secret_g1[0..3]: [q_M]_1, [q_L]_1, [q_R]_1,
///        [q_O]_1
/// \param S_polys_at_secret_g1[0..2]: [S_sigma1]_1, [S_sigma2]_1,
///        [S_sigma3]_1
/// \param srs.secret_powers_g2[1]: [secret]_2 = secret^1 * G2
///
/// Public input polynomial
///
/// srs.PI_poly: w_i, 0\le{i}<l<n

/// Verifier preprocessed input
template<typename ppT> struct verifier_preprocessed_input_t {
    /// circuit selector polynomials Q evaluated at the secret input
    std::vector<libff::G1<ppT>> Q_polys_at_secret_g1;
    /// permutation polynomials S evaluated at the secret input
    std::vector<libff::G1<ppT>> S_polys_at_secret_g1;

    verifier_preprocessed_input_t(
        std::vector<libff::G1<ppT>> &&Q_polys_at_secret_g1,
        std::vector<libff::G1<ppT>> &&S_polys_at_secret_g1);
};

/// Verifier step 4 output
template<typename ppT> struct step_four_out_t {
    /// permutation challenge - hashes of transcript (round 1)
    const libff::Fr<ppT> beta;
    /// permutation challenge - hashes of transcript (round 1)
    const libff::Fr<ppT> gamma;
    /// quotinet challenge - hash of transcript (round 2)
    const libff::Fr<ppT> alpha;
    /// evaluation challenge - hash of transcript (round 3)
    const libff::Fr<ppT> zeta;
    /// opening challenge - hash of transcript (round 4) (denoted by v in
    /// [GWC19])
    const libff::Fr<ppT> nu;
    /// multipoint evaluation challenge - hash of transcript (round 5)
    const libff::Fr<ppT> u;
    step_four_out_t(
        const libff::Fr<ppT> &beta,
        const libff::Fr<ppT> &gamma,
        const libff::Fr<ppT> &alpha,
        const libff::Fr<ppT> &zeta,
        const libff::Fr<ppT> &nu,
        const libff::Fr<ppT> &u);
};

/// Verifier step 5 output
template<typename ppT> struct step_five_out_t {
    /// evaluation of vanishing polynomial Zh at x=zeta i.e. Zh(zeta)
    libff::Fr<ppT> zh_zeta;
    step_five_out_t(libff::Fr<ppT> &&zh_zeta);
};

/// Verifier step 6 output
template<typename ppT> struct step_six_out_t {
    /// Lagrange polynomial evaluation of polynomial L1 at x=zeta
    /// i.e. L1(zeta)
    libff::Fr<ppT> L_0_zeta;
    step_six_out_t(libff::Fr<ppT> &&L_0_zeta);
};

/// Verifier step 7 output
template<typename ppT> struct step_seven_out_t {
    /// Public input polynomial PI evaluated at x=zeta .e. PI(zeta)
    libff::Fr<ppT> PI_zeta;
    step_seven_out_t(libff::Fr<ppT> &&PI_zeta);
};

/// Verifier step 8 output
template<typename ppT> struct step_eight_out_t {
    /// compute quotient polynomial evaluation r'(zeta) = r(zeta) - r0,
    /// where r0 is a constant term Note:
    libff::Fr<ppT> r_prime_zeta;
    step_eight_out_t(libff::Fr<ppT> &&r_prime_zeta);
};

/// Verifier step 9 output
template<typename ppT> struct step_nine_out_t {
    /// first part of batched polynomial commitment [D]_1
    libff::G1<ppT> D1;
    step_nine_out_t(libff::G1<ppT> &&D1);
};

/// Verifier step 10 output
template<typename ppT> struct step_ten_out_t {
    /// full batched polynomial commitment [F]_1 = [D]_1 + v [a]_1 + v^2
    /// [b]_1 + v^3 [c]_1 + v^4 [s_sigma_1]_1 + v^5 [s_sigma_2]_1
    libff::G1<ppT> F1;
    step_ten_out_t(libff::G1<ppT> &&F1);
};

/// Verifier step 11 output
template<typename ppT> struct step_eleven_out_t {
    /// group-encoded batch evaluation [E]_1
    libff::G1<ppT> E1;
    step_eleven_out_t(libff::G1<ppT> &&E1);
};

/// Plonk verifier. Verifies object of class plonk_proof.
template<typename ppT, class transcript_hasher> class plonk_verifier
{
    using Field = libff::Fr<ppT>;

public:
    /// Verifier precomputation
    ///
    /// INPUT
    /// \param[in] srs: structured reference string
    ///
    /// OUTPUT
    /// \param[out] verifier_preprocessed_input: see
    ///             verifier_preprocessed_input_t
    static verifier_preprocessed_input_t<ppT> preprocessed_input(
        const srs<ppT> &srs);

    /// Verifier Step 1: validate that elements belong to group G1
    ///
    /// \attention This validation MUST be done by the caller. Empty
    /// function here for consistency with the description in [GWC19]
    static void step_one(const plonk_proof<ppT> &proof);

    /// Verifier Step 2: validate that elements belong to scalar field Fr
    ///
    /// \attention This validation MUST be done by the caller. Empty
    /// function here for consistency with the description in [GWC19]
    static void step_two(const plonk_proof<ppT> &proof);

    /// Verifier Step 3: validate that the public input belongs to scalar
    /// field Fr
    ///
    /// \attention This validation MUST be done by the caller. Empty
    /// function here for consistency with the description in [GWC19]
    static void step_three(const srs<ppT> &srs);

    /// Verifier Step 4: compute challenges hashed transcript as in prover
    /// description, from the common inputs, public input, and elements of
    /// pi-SNARK. TODO: fixed to the test vectors for now
    ///
    /// INPUT
    /// \param[in] proof: SNARK proof produced by the prover
    /// \param[in] transcript_hasher: hashes of the communication
    ///            transcript after prover rounds 1,2,3,4,5.
    ///
    /// OUTPUT
    /// \param[out] step_four_out: see step_four_out
    static step_four_out_t<ppT> step_four(
        const plonk_proof<ppT> &proof, transcript_hasher &hasher);

    /// Verifier Step 5: compute zero polynomial evaluation
    ///
    /// INPUT
    /// \param[in] step_four_out: see step_four_out_t
    /// \param[in] domain: libfqfft evaluation domain (see
    ///            libfqfft/evaluation_domain/evaluation_domain.hpp)
    ///
    /// OUTPUT
    /// \param[out] step_five_out: see step_five_out_t
    static step_five_out_t<ppT> step_five(
        const step_four_out_t<ppT> &step_four_out,
        std::shared_ptr<libfqfft::evaluation_domain<libff::Fr<ppT>>> domain);

    /// Verifier Step 6: Compute Lagrange polynomial evaluation L1(zeta)
    /// Note: the paper counts the L-polynomials from 1; we count from 0
    ///
    /// INPUT
    /// \param[in] step_four_out: see step_four_out_t
    /// \param[in] srs: structured reference string containing also
    ///            circuit-specific information
    ///
    /// OUTPUT
    /// \param[out] step_six_out: see step_six_out_t
    static step_six_out_t<ppT> step_six(
        const step_four_out_t<ppT> &step_four_out, const srs<ppT> &srs);

    /// Verifier Step 7: compute public input polynomial evaluation
    /// PI(zeta)
    ///
    /// INPUT
    /// \param[in] step_four_out: see step_four_out_t
    /// \param[in] srs: structured reference string containing also
    ///            circuit-specific information
    ///
    /// OUTPUT
    /// \param[out] step_seven_out: see step_seven_out_t
    static step_seven_out_t<ppT> step_seven(
        const step_four_out_t<ppT> &step_four_out,
        const std::vector<Field> &PI_value_list,
        const srs<ppT> &srs);

    /// Verifier Step 8: compute quotient polynomial evaluation r'(zeta) =
    /// r(zeta) - r0, where r0 is a constant term \note follows the Python
    /// reference implementation, which slightly deviates from the paper
    /// due to the presence of the r_zeta term in the proof (not present
    /// in the paper).  In particular, the reference code computes and
    /// uses r'(zeta) in step 8, while the paper uses r0. In addition, the
    /// reference code divides r'(zeta) by the vanishing polynomial at
    /// zeta zh_zeta, while the paper does not do that (see also Step 9).
    ///
    /// INPUT
    /// \param[in] step_four_out: see step_four_out_t
    /// \param[in] step_five_out: see step_five_out_t
    /// \param[in] step_six_out: see step_six_out_t
    /// \param[in] step_seven_out: see step_seven_out_t
    /// \param[in] proof: SNARK proof produced by the prover
    ///
    /// OUTPUT
    /// \param[out] step_eight_out: see step_eight_out_t
    static step_eight_out_t<ppT> step_eight(
        const step_four_out_t<ppT> &step_four_out,
        const step_five_out_t<ppT> &step_five_out,
        const step_six_out_t<ppT> &step_six_out,
        const step_seven_out_t<ppT> &step_seven_out,
        const plonk_proof<ppT> &proof);

    /// Verifier Step 9: compute first part of batched polynomial
    /// commitment [D]_1 Note: the reference implemention differs from the
    /// paper -- it does not add the following term to D1, but to F1 (Step
    /// 10): -Zh(zeta)([t_lo]_1 + zeta^n [t_mid]_1 + zeta^2n
    /// [t_hi]_1). Instead ([t_lo]_1 + zeta^n [t_mid]_1 + zeta^2n
    /// [t_hi]_1) is added to F1 in Step 10 and the multiplication by
    /// Zh(zeta) is accounted for by dividing by Zh(zeta) of r'(zeta) in
    /// Step 8.
    ///
    /// INPUT
    /// \param[in] step_four_out: see step_four_out_t
    /// \param[in] step_six_out: see step_six_out_t
    /// \param[in] proof: SNARK proof produced by the prover
    /// \param[in] preprocessed_input: verifier preprocessed input
    /// \param[in] srs: structured reference string containing also
    ///            circuit-specific information
    ///
    /// OUTPUT
    /// \param[out] step_nine_out: see step_nine_out_t
    static step_nine_out_t<ppT> step_nine(
        const step_four_out_t<ppT> &step_four_out,
        const step_six_out_t<ppT> &step_six_out,
        const plonk_proof<ppT> &proof,
        const verifier_preprocessed_input_t<ppT> &preprocessed_input,
        const srs<ppT> &srs);

    /// Verifier Step 10: compute full batched polynomial commitment [F]_1
    /// = [D]_1 + v [a]_1 + v^2 [b]_1 + v^3 [c]_1 + v^4 [s_sigma_1]_1 +
    /// v^5 [s_sigma_2]_1 Note: to [F]_1 the erefernce code also adds the
    /// term ([t_lo]_1 + zeta^n [t_mid]_1 + zeta^2n [t_hi]_1) which is
    /// addedto [D]_1 in the paper (see commenst to Steps 8,9)
    ///
    /// INPUT
    /// \param[in] step_four_out: see step_four_out_t
    /// \param[in] step_nine_out: see step_nine_out_t
    /// \param[in] proof: SNARK proof produced by the prover
    /// \param[in] preprocessed_input: verifier preprocessed input
    /// \param[in] srs: structured reference string containing also
    ///            circuit-specific information
    ///
    /// OUTPUT
    /// \param[out] step_ten_out: see step_ten_out_t
    static step_ten_out_t<ppT> step_ten(
        const step_four_out_t<ppT> &step_four_out,
        const step_nine_out_t<ppT> &step_nine_out,
        const plonk_proof<ppT> &proof,
        const verifier_preprocessed_input_t<ppT> &preprocessed_input,
        const srs<ppT> &srs);

    /// Verifier Step 11: compute group-encoded batch evaluation [E]_1
    ///
    /// INPUT
    /// \param[in] step_four_out: see step_four_out_t
    /// \param[in] step_eight_out: see step_eight_out_t
    /// \param[in] proof: SNARK proof produced by the prover
    ///
    /// OUTPUT
    /// \param[out] step_eleven_out: see step_eleven_out_t
    static step_eleven_out_t<ppT> step_eleven(
        const step_four_out_t<ppT> &step_four_out,
        const step_eight_out_t<ppT> &step_eight_out,
        const plonk_proof<ppT> &proof);

    /// Verifier Step 12: batch validate all evaluations
    ///
    /// Checks the following equality
    ///
    /// e( [W_zeta]_1 + u [W_{zeta srs.omega_roots}]_1, [x]_2 ) * e( -zeta
    /// [W_zeta ]_1 - u zeta srs.omega_roots [W_{zeta srs.omega_roots}]_1
    /// - [F]_1 + [E]_1, [1]_2 ) = Field(1)
    ///
    /// Denoted as:
    ///
    /// e(first_lhs, second_lhs) * e(first_rhs, second_rhs) = 1
    ///
    /// INPUT
    /// \param[in] step_four_out: see step_four_out_t
    /// \param[in] step_ten_out: see step_ten_out_t
    /// \param[in] step_eleven_out: see step_eleven_out_t
    /// \param[in] proof: SNARK proof produced by the prover
    /// \param[in] srs: structured reference string containing also
    ///            circuit-specific information
    ///
    /// OUTPUT
    /// \param[out] boolean 1/0 = valid/invalid proof
    static bool step_twelve(
        const step_four_out_t<ppT> &step_four_out,
        const step_ten_out_t<ppT> &step_ten_out,
        const step_eleven_out_t<ppT> &step_eleven_out,
        const plonk_proof<ppT> &proof,
        const srs<ppT> &srs);

    /// \attention The first three steps (as given in [GWC19] -- see
    /// below) MUST be executed by the caller:
    ///
    /// Verifier Step 1: validate that elements belong to group G1
    /// Verifier Step 2: validate that elements belong to scalar field Fr
    /// Verifier Step 3: validate that the public input belongs to scalar field
    /// Fr
    ///
    /// Therefore verification starts from Step 4 of [GWC19]
    ///
    /// INPUT
    /// \param[in] proof: SNARK proof produced by the prover
    /// \param[in] srs: structured reference string containing also
    ///            circuit-specific information
    /// \param[in] PI_value_list: list of values corresponding to
    ///            public inputs
    /// \param[in] transcript_hasher: hashes of the communication
    ///            transcript after prover rounds 1,2,3,4,5.
    ///
    /// OUTPUT
    /// \param[out] boolean 1/0 = valid/invalid proof
    bool verify_proof(
        const plonk_proof<ppT> &proof,
        const srs<ppT> &srs,
        const std::vector<Field> &PI_value_list,
        transcript_hasher &hasher);
};

} // namespace libsnark

#include "libsnark/zk_proof_systems/plonk/verifier.tcc"

#endif // LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_VERIFIER_HPP_
