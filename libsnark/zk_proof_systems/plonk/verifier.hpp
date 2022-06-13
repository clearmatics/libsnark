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

// Verifier preprocessed input
template<typename ppT> struct verifier_preprocessed_input_t {
    std::vector<libff::G1<ppT>> Q_polys_at_secret_g1;
    std::vector<libff::G1<ppT>> S_polys_at_secret_g1;
};

// Verifier step 4 output
template<typename ppT> struct step_four_out_t {
    libff::Fr<ppT> beta;
    libff::Fr<ppT> gamma;
    libff::Fr<ppT> alpha;
    libff::Fr<ppT> zeta;
    libff::Fr<ppT> nu;
    libff::Fr<ppT> u;
    step_four_out_t(
        libff::Fr<ppT> &&beta,
        libff::Fr<ppT> &&gamma,
        libff::Fr<ppT> &&alpha,
        libff::Fr<ppT> &&zeta,
        libff::Fr<ppT> &&nu,
        libff::Fr<ppT> &&u)
        : beta(beta), gamma(gamma), alpha(alpha), zeta(zeta), nu(nu), u(u)
    {
    }
};

// Verifier step 5 output
template<typename ppT> struct step_five_out_t {
    // evaluation of vanishing polynomial Zh at x=zeta i.e. Zh(zeta)
    libff::Fr<ppT> zh_zeta;
    step_five_out_t(libff::Fr<ppT> &&zh_zeta) : zh_zeta(zh_zeta) {}
};

// Verifier step 6 output
template<typename ppT> struct step_six_out_t {
    // Lagrange polynomial evaluation of polynomial L1 at x=zeta
    // i.e. L1(zeta)
    libff::Fr<ppT> L_0_zeta;
    step_six_out_t(libff::Fr<ppT> &&L_0_zeta) : L_0_zeta(L_0_zeta) {}
};

// Verifier step 7 output
template<typename ppT> struct step_seven_out_t {
    // Public input polynomial PI evaluated at x=zeta .e. PI(zeta)
    libff::Fr<ppT> PI_zeta;
    step_seven_out_t(libff::Fr<ppT> &&PI_zeta) : PI_zeta(PI_zeta) {}
};

// Verifier step 8 output
template<typename ppT> struct step_eight_out_t {
    // compute quotient polynomial evaluation r'(zeta) = r(zeta) - r0,
    // where r0 is a constant term Note:
    libff::Fr<ppT> r_prime_zeta;
    step_eight_out_t(libff::Fr<ppT> &&r_prime_zeta) : r_prime_zeta(r_prime_zeta)
    {
    }
};

// Verifier step 9 output
template<typename ppT> struct step_nine_out_t {
    // first part of batched polynomial commitment [D]_1
    libff::G1<ppT> D1;
    step_nine_out_t(libff::G1<ppT> &&D1) : D1(D1) {}
};

// Verifier step 10 output
template<typename ppT> struct step_ten_out_t {
    // full batched polynomial commitment [F]_1 = [D]_1 + v [a]_1 + v^2
    // [b]_1 + v^3 [c]_1 + v^4 [s_sigma_1]_1 + v^5 [s_sigma_2]_1
    libff::G1<ppT> F1;
    step_ten_out_t(libff::G1<ppT> &&F1) : F1(F1) {}
};

// Verifier step 11 output
template<typename ppT> struct step_eleven_out_t {
    // group-encoded batch evaluation [E]_1
    libff::G1<ppT> E1;
    step_eleven_out_t(libff::G1<ppT> &&E1) : E1(E1) {}
};

/// Plonk verifier. Verifies object of class plonk_proof.
template<typename ppT> class plonk_verifier
{
    using Field = libff::Fr<ppT>;

public:
    static verifier_preprocessed_input_t<ppT> preprocessed_input(
        const srs<ppT> &srs);

    static void step_one(const plonk_proof<ppT> &proof);
    static void step_two(const plonk_proof<ppT> &proof);
    static void step_three(const srs<ppT> &srs);

    static step_four_out_t<ppT> step_four();
    static step_five_out_t<ppT> step_five(
        const step_four_out_t<ppT> &step_four_out, const srs<ppT> &srs);
    static step_six_out_t<ppT> step_six(
        const step_four_out_t<ppT> &step_four_out, const srs<ppT> &srs);
    static step_seven_out_t<ppT> step_seven(
        const step_four_out_t<ppT> &step_four_out, const srs<ppT> &srs);
    static step_eight_out_t<ppT> step_eight(
        const step_four_out_t<ppT> &step_four_out,
        const step_five_out_t<ppT> &step_five_out,
        const step_six_out_t<ppT> &step_six_out,
        const step_seven_out_t<ppT> &step_seven_out,
        const plonk_proof<ppT> &proof);
    static step_nine_out_t<ppT> step_nine(
        const step_four_out_t<ppT> &step_four_out,
        const step_six_out_t<ppT> &step_six_out,
        const plonk_proof<ppT> &proof,
        const verifier_preprocessed_input_t<ppT> &preprocessed_input,
        const srs<ppT> &srs);
    static step_ten_out_t<ppT> step_ten(
        const step_four_out_t<ppT> &step_four_out,
        const step_nine_out_t<ppT> &step_nine_out,
        const plonk_proof<ppT> &proof,
        const verifier_preprocessed_input_t<ppT> &preprocessed_input,
        const srs<ppT> &srs);
    static step_eleven_out_t<ppT> step_eleven(
        const step_four_out_t<ppT> &step_four_out,
        const step_eight_out_t<ppT> &step_eight_out,
        const plonk_proof<ppT> &proof);
    static bool step_twelve(
        const step_four_out_t<ppT> &step_four_out,
        const step_ten_out_t<ppT> &step_ten_out,
        const step_eleven_out_t<ppT> &step_eleven_out,
        const plonk_proof<ppT> &proof,
        const srs<ppT> &srs);
    bool verify_proof(const plonk_proof<ppT> &proof, const srs<ppT> &srs);
};

} // namespace libsnark

#include "libsnark/zk_proof_systems/plonk/verifier.tcc"

#endif // LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_VERIFIER_HPP_
