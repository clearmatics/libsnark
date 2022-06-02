/** @file
*****************************************************************************

Declaration of Verifier interfaces for ppzkSNARK proof system Plonk.

This includes:
- class for verifier

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

#ifndef PLONK_PPZKSNARK_VERIFIER_HPP_
#define PLONK_PPZKSNARK_VERIFIER_HPP_

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
};

// Verifier step 5 output
template<typename ppT> struct step_five_out_t {
    // evaluation of vanishing polynomial Zh at x=zeta i.e. Zh(zeta)
    libff::Fr<ppT> zh_zeta;
};

// Verifier step 6 output
template<typename ppT> struct step_six_out_t {
    // Lagrange polynomial evaluation of polynomial L1 at x=zeta
    // i.e. L1(zeta)
    libff::Fr<ppT> L_0_zeta;
};

// Verifier step 7 output
template<typename ppT> struct step_seven_out_t {
    // Public input polynomial PI evaluated at x=zeta .e. PI(zeta)
    libff::Fr<ppT> PI_zeta;
};

// Verifier step 8 output
template<typename ppT> struct step_eight_out_t {
    // compute quotient polynomial evaluation r'(zeta) = r(zeta) - r0,
    // where r0 is a constant term Note:
    libff::Fr<ppT> r_prime_zeta;
};

// Verifier step 9 output
template<typename ppT> struct step_nine_out_t {
    // first part of batched polynomial commitment [D]_1
    libff::G1<ppT> D1;
};

// Verifier step 10 output
template<typename ppT> struct step_ten_out_t {
    // full batched polynomial commitment [F]_1 = [D]_1 + v [a]_1 + v^2
    // [b]_1 + v^3 [c]_1 + v^4 [s_sigma_1]_1 + v^5 [s_sigma_2]_1
    libff::G1<ppT> F1;
};

// Verifier step 11 output
template<typename ppT> struct step_eleven_out_t {
    // group-encoded batch evaluation [E]_1
    libff::G1<ppT> E1;
};

/**
 * Plonk verifier. Verifies object of class plonk_proof.
 */
template<typename ppT> class plonk_verifier
{
    using Field = libff::Fr<ppT>;

public:
    plonk_verifier(){};

    static verifier_preprocessed_input_t<ppT> preprocessed_input(
        const srs<ppT> srs, const common_preprocessed_input<ppT> common_input);

    static void step_one(const plonk_proof<ppT> proof);
    static void step_two(const plonk_proof<ppT> proof);
    static void step_three(const common_preprocessed_input<ppT> common_input);

    static step_four_out_t<ppT> step_four();
    static step_five_out_t<ppT> step_five(
        const step_four_out_t<ppT> step_four_out,
        const common_preprocessed_input<ppT> common_input);
    static step_six_out_t<ppT> step_six(
        const step_four_out_t<ppT> step_four_out,
        const common_preprocessed_input<ppT> common_input);
    static step_seven_out_t<ppT> step_seven(
        const step_four_out_t<ppT> step_four_out,
        const common_preprocessed_input<ppT> common_input);
    static step_eight_out_t<ppT> step_eight(
        const step_four_out_t<ppT> step_four_out,
        const step_five_out_t<ppT> step_five_out,
        const step_six_out_t<ppT> step_six_out,
        const step_seven_out_t<ppT> step_seven_out,
        const plonk_proof<ppT> proof);
    static step_nine_out_t<ppT> step_nine(
        const step_four_out_t<ppT> step_four_out,
        const step_six_out_t<ppT> step_six_out,
        const plonk_proof<ppT> proof,
        const verifier_preprocessed_input_t<ppT> preprocessed_input,
        const common_preprocessed_input<ppT> common_input);
    static step_ten_out_t<ppT> step_ten(
        const step_four_out_t<ppT> step_four_out,
        const step_nine_out_t<ppT> step_nine_out,
        const plonk_proof<ppT> proof,
        const verifier_preprocessed_input_t<ppT> preprocessed_input,
        const common_preprocessed_input<ppT> common_input);
    static step_eleven_out_t<ppT> step_eleven(
        const step_four_out_t<ppT> step_four_out,
        const step_eight_out_t<ppT> step_eight_out,
        const plonk_proof<ppT> proof);
    static bool step_twelve(
        const step_four_out_t<ppT> step_four_out,
        const step_ten_out_t<ppT> step_ten_out,
        const step_eleven_out_t<ppT> step_eleven_out,
        const plonk_proof<ppT> proof,
        const srs<ppT> srs,
        const common_preprocessed_input<ppT> common_input);
    bool verify_proof(
        const plonk_proof<ppT> proof,
        const srs<ppT> srs,
        const common_preprocessed_input<ppT> common_input);
};

} // namespace libsnark

#include "libsnark/zk_proof_systems/plonk/verifier.tcc"

#endif // PLONK_PPZKSNARK_VERIFIER_HPP_
