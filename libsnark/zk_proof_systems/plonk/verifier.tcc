/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_VERIFIER_TCC_
#define LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_VERIFIER_TCC_

/// Implementation of Verifier interfaces for a ppzkSNARK for Plonk. See
/// verifier.hpp .

namespace libsnark
{

/// Verifier precomputation
///
/// INPUT
/// \param[in] srs: structured reference string
///
/// OUTPUT
/// \param[out] Q_polys_at_secret_g1: circuit selector polynomials Q evaluated
/// at
///   the secret input
/// \param[out] S_polys_at_secret_g1: permutation polynomials S evaluated at the
///   secret input
template<typename ppT>
verifier_preprocessed_input_t<ppT> plonk_verifier<ppT>::preprocessed_input(
    const srs<ppT> &srs)
{
    verifier_preprocessed_input_t<ppT> preprocessed_input;
    preprocessed_input.Q_polys_at_secret_g1.resize(srs.Q_polys.size());
    plonk_evaluate_polys_at_secret_G1<ppT>(
        srs.secret_powers_g1,
        srs.Q_polys,
        preprocessed_input.Q_polys_at_secret_g1);

    preprocessed_input.S_polys_at_secret_g1.resize(srs.S_polys.size());
    plonk_evaluate_polys_at_secret_G1<ppT>(
        srs.secret_powers_g1,
        srs.S_polys,
        preprocessed_input.S_polys_at_secret_g1);
    return preprocessed_input;
}

/// Verifier Step 1: validate that elements belong to group G1
///
/// \attention This validation MUST be done by the caller. Empty
/// function here for consistency with the description in [GWC19]
template<typename ppT>
void plonk_verifier<ppT>::step_one(const plonk_proof<ppT> &proof)
{
}

/// Verifier Step 2: validate that elements belong to scalar field Fr
///
/// \attention This validation MUST be done by the caller. Empty
/// function here for consistency with the description in [GWC19]
template<typename ppT>
void plonk_verifier<ppT>::step_two(const plonk_proof<ppT> &proof)
{
}

/// Verifier Step 3: validate that the public input belongs to scalar
/// field Fr
///
/// \attention This validation MUST be done by the caller. Empty
/// function here for consistency with the description in [GWC19]
template<typename ppT> void plonk_verifier<ppT>::step_three(const srs<ppT> &srs)
{
}

/// Verifier Step 4: compute challenges hashed transcript as in prover
/// description, from the common inputs, public input, and elements of
/// pi-SNARK. TODO: fixed to the test vectors for now
///
/// OUTPUT
/// \param[out] beta, gamma: permutation challenges - hashes of
///             transcript
/// \param[out] alpha: quotinet challenge - hash of transcript
/// \param[out] zeta: evaluation challenge - hash of transcript
/// \param[out] nu: opening challenge - hash of transcript (denoted by
///             v in [GWC19])
/// \param[out] u: multipoint evaluation challenge - hash of
///             transcript
template<typename ppT> step_four_out_t<ppT> plonk_verifier<ppT>::step_four()
{
    // the example class is defined specifically for the BLS12-381
    // curve, so make sure we are using this curve TODO: remove when
    // the implementation is stable and tested
    try {
        plonk_exception_assert_curve_bls12_381<ppT>();
    } catch (const std::domain_error &e) {
        std::cout << "Error: " << e.what() << "\n";
        exit(EXIT_FAILURE);
    }
    // load challenges from example for debug
    plonk_example example;
    // step 4 output
    libff::Fr<ppT> beta;
    libff::Fr<ppT> gamma;
    libff::Fr<ppT> alpha;
    libff::Fr<ppT> zeta;
    libff::Fr<ppT> nu;
    libff::Fr<ppT> u;

    beta = example.beta;
    gamma = example.gamma;
    alpha = example.alpha;
    zeta = example.zeta;
    nu = example.nu;
    u = example.u;

    step_four_out_t<ppT> step_four_out(
        std::move(beta),
        std::move(gamma),
        std::move(alpha),
        std::move(zeta),
        std::move(nu),
        std::move(u));

    return step_four_out;
}

/// Verifier Step 5: compute zero polynomial evaluation
///
/// INPUT
/// \param[in] zeta: evaluation challenge -- hash of transcript (from
///            step 4)
/// \param[in] srs: structured reference string containing also
///            circuit-specific information
///
/// OUTPUT
/// \param[out] zh_zeta: evaluation of vanishing polynomial Zh at
///             x=zeta i.e. Zh(zeta)
template<typename ppT>
step_five_out_t<ppT> plonk_verifier<ppT>::step_five(
    const step_four_out_t<ppT> &step_four_out, const srs<ppT> &srs)
{
    libff::Fr<ppT> zh_zeta;
    std::shared_ptr<libfqfft::evaluation_domain<Field>> domain =
        libfqfft::get_evaluation_domain<Field>(srs.num_gates);
    zh_zeta = domain->compute_vanishing_polynomial(step_four_out.zeta);
    step_five_out_t<ppT> step_five_out(std::move(zh_zeta));
    return step_five_out;
}

/// Verifier Step 6: Compute Lagrange polynomial evaluation L1(zeta)
/// Note: the paper counts the L-polynomials from 1; we count from 0
///
/// INPUT
/// \param[in] zeta: evaluation challenge -- hash of transcript (from
///            step 4)
/// \param[in] srs: structured reference string containing also
///            circuit-specific information
///
/// OUTPUT
/// \param[out] L_0_zeta: Lagrange polynomial evaluation of polynomial
///             L1 at x=zeta i.e. L1(zeta)
template<typename ppT>
step_six_out_t<ppT> plonk_verifier<ppT>::step_six(
    const step_four_out_t<ppT> &step_four_out, const srs<ppT> &srs)
{
    libff::Fr<ppT> L_0_zeta;
    L_0_zeta = libfqfft::evaluate_polynomial<Field>(
        srs.L_basis[0].size(), srs.L_basis[0], step_four_out.zeta);
    step_six_out_t<ppT> step_six_out(std::move(L_0_zeta));
    return step_six_out;
}

/// Verifier Step 7: compute public input polynomial evaluation
/// PI(zeta)
///
/// INPUT
/// \param[in] zeta: evaluation challenge -- hash of transcript (from
///            step 4)
/// \param[in] srs: structured reference string containing also
///            circuit-specific information
///
/// OUTPUT
/// \param[out] PI_zeta: public input polynomial PI evaluated at
///             x=zeta i.e. PI(zeta)
template<typename ppT>
step_seven_out_t<ppT> plonk_verifier<ppT>::step_seven(
    const step_four_out_t<ppT> &step_four_out, const srs<ppT> &srs)
{
    libff::Fr<ppT> PI_zeta;
    PI_zeta = libfqfft::evaluate_polynomial<Field>(
        srs.PI_poly.size(), srs.PI_poly, step_four_out.zeta);
    step_seven_out_t<ppT> step_seven_out(std::move(PI_zeta));
    return step_seven_out;
}

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
/// \param[in] beta, gamma: permutation challenges -- hashes of
///            transcript (from step 4)
/// \param[in] alpha: quotinet challenge -- hash of transcript (from
///            step 4)
/// \param[in] zeta: evaluation challenge -- hash of transcript (from
///            step 4)
/// \param[in] zh_zeta: evaluation of vanishing polynomial Zh at
///            x=zeta i.e. Zh(zeta) (from step 5)
/// \param[in] L_0_zeta: Lagrange polynomial evaluation of polynomial
///            L1 at x=zeta i.e. L1(zeta) (from step 6)
/// \param[in] PI_zeta: public input polynomial PI evaluated at x=zeta
///            i.e. PI(zeta) (from step 7)
/// \param[in] proof: SNARK proof produced by the prover
///
/// OUTPUT
/// \param[out] r_prime_zeta: quotient polynomial evaluation r'(zeta)
///             = r(zeta) - r0, where r0 is a constant term
template<typename ppT>
step_eight_out_t<ppT> plonk_verifier<ppT>::step_eight(
    const step_four_out_t<ppT> &step_four_out,
    const step_five_out_t<ppT> &step_five_out,
    const step_six_out_t<ppT> &step_six_out,
    const step_seven_out_t<ppT> &step_seven_out,
    const plonk_proof<ppT> &proof)
{
    libff::Fr<ppT> r_prime_zeta;
    Field alpha_power2 = libff::power(step_four_out.alpha, libff::bigint<1>(2));

    // compute polynomial r'(zeta) = r(zeta) - r_0
    std::vector<Field> r_prime_parts(5);
    r_prime_parts[0] = proof.r_zeta + step_seven_out.PI_zeta;
    r_prime_parts[1] =
        (proof.a_zeta + (step_four_out.beta * proof.S_0_zeta) +
         step_four_out.gamma);
    r_prime_parts[2] =
        (proof.b_zeta + (step_four_out.beta * proof.S_1_zeta) +
         step_four_out.gamma);
    r_prime_parts[3] = (proof.c_zeta + step_four_out.gamma) *
                       proof.z_poly_xomega_zeta * step_four_out.alpha;
    r_prime_parts[4] = (step_six_out.L_0_zeta * alpha_power2);
    r_prime_zeta = (r_prime_parts[0] -
                    (r_prime_parts[1] * r_prime_parts[2] * r_prime_parts[3]) -
                    r_prime_parts[4]) *
                   step_five_out.zh_zeta.inverse();

    step_eight_out_t<ppT> step_eight_out(std::move(r_prime_zeta));
    return step_eight_out;
}

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
/// \param[in] beta, gamma: permutation challenges -- hashes of
///            transcript (from step 4)
/// \param[in] alpha: quotinet challenge -- hash of transcript (from
///            step 4)
/// \param[in] zeta: evaluation challenge -- hash of transcript (from
///            step 4)
/// \param[in] nu: opening challenge -- hash of transcript (denoted by
///            v in [GWC19]) (from step 4)
/// \param[in] u: multipoint evaluation challenge -- hash of
///            transcript (from step 4)
/// \param[in] L_0_zeta: Lagrange polynomial evaluation of polynomial
///            L1 at x=zeta i.e. L1(zeta) (from step 6)
/// \param[in] Q_polys_at_secret_g1: circuit selector polynomials Q
///            evaluated at the secret input (from verifier
///            preprocessed input)
/// \param[in] S_polys_at_secret_g1: permutation polynomials S
///            evaluated at the secret input (from verifier
///            preprocessed input)
/// \param[in] proof: SNARK proof produced by the prover
/// \param[in] preprocessed_input: verifier preprocessed input
/// \param[in] srs: structured reference string containing also
///            circuit-specific information
///
/// OUTPUT
/// \param[out] D1: first part of batched polynomial commitment [D]_1
template<typename ppT>
step_nine_out_t<ppT> plonk_verifier<ppT>::step_nine(
    const step_four_out_t<ppT> &step_four_out,
    const step_six_out_t<ppT> &step_six_out,
    const plonk_proof<ppT> &proof,
    const verifier_preprocessed_input_t<ppT> &preprocessed_input,
    const srs<ppT> &srs)
{
    libff::G1<ppT> D1;

    // D1 is computed in 3 parts
    std::vector<libff::G1<ppT>> D1_part(3);

    Field alpha_power2 = libff::power(step_four_out.alpha, libff::bigint<1>(2));

    // compute D1_part[0]:
    // (a_bar b_bar [q_M]_1 + a_bar [q_L]_1 + b_bar [q_R]_1 + c_bar [q_O]_1 +
    // [q_C]_1) nu Note: the paper omits the final multiplication by nu
    std::vector<libff::G1<ppT>> curve_points{
        preprocessed_input.Q_polys_at_secret_g1[M],
        preprocessed_input.Q_polys_at_secret_g1[L],
        preprocessed_input.Q_polys_at_secret_g1[R],
        preprocessed_input.Q_polys_at_secret_g1[O],
        preprocessed_input.Q_polys_at_secret_g1[C]};
    std::vector<libff::Fr<ppT>> scalar_elements{
        proof.a_zeta * proof.b_zeta * step_four_out.nu,
        proof.a_zeta * step_four_out.nu,
        proof.b_zeta * step_four_out.nu,
        proof.c_zeta * step_four_out.nu,
        step_four_out.nu};
    D1_part[0] = plonk_multi_exp_G1<ppT>(curve_points, scalar_elements);

    // compute D1_part[1]:
    // ((a_bar + beta zeta + gamma)(b_bar + beta k1 zeta + gamma)(c_bar + beta
    // k2 zeta + gamma) alpha + L1(zeta) alpha^2 + u) [z]_1
    Field D1_part1_scalar =
        (proof.a_zeta + (step_four_out.beta * step_four_out.zeta) +
         step_four_out.gamma) *
            (proof.b_zeta + (step_four_out.beta * srs.k1 * step_four_out.zeta) +
             step_four_out.gamma) *
            (proof.c_zeta + (step_four_out.beta * srs.k2 * step_four_out.zeta) +
             step_four_out.gamma) *
            step_four_out.alpha * step_four_out.nu +
        step_six_out.L_0_zeta * alpha_power2 * step_four_out.nu +
        step_four_out.u;
    D1_part[1] = D1_part1_scalar * proof.z_poly_at_secret_g1;

    // compute D1_part[2]:
    // (a_bar + beta s_sigma1_bar + gamma)(b_bar + beta s_sigma2_bar +
    // gamma)alpha beta z_preprocessed_input.omega_roots_bar [s_sigma3]_1
    Field D1_part2_scalar =
        ((proof.a_zeta + (step_four_out.beta * proof.S_0_zeta) +
          step_four_out.gamma) *
         (proof.b_zeta + (step_four_out.beta * proof.S_1_zeta) +
          step_four_out.gamma) *
         step_four_out.alpha * step_four_out.beta * proof.z_poly_xomega_zeta *
         step_four_out.nu) *
        Field(-1);
    D1_part[2] = D1_part2_scalar * preprocessed_input.S_polys_at_secret_g1[2];

    // Compute D1 = D1_part[0] + D1_part[1] + D1_part[2]
    D1 = D1_part[0] + D1_part[1] + D1_part[2];

    step_nine_out_t<ppT> step_nine_out(std::move(D1));
    return step_nine_out;
}

/// Verifier Step 10: compute full batched polynomial commitment [F]_1
/// = [D]_1 + v [a]_1 + v^2 [b]_1 + v^3 [c]_1 + v^4 [s_sigma_1]_1 +
/// v^5 [s_sigma_2]_1 Note: to [F]_1 the erefernce code also adds the
/// term ([t_lo]_1 + zeta^n [t_mid]_1 + zeta^2n [t_hi]_1) which is
/// addedto [D]_1 in the paper (see commenst to Steps 8,9)
///
/// INPUT
/// \param[in] zeta: evaluation challenge -- hash of transcript (from
///            step 4)
/// \param[in] nu: opening challenge -- hash of transcript (denoted by
///            v in [GWC19]) (from step 4)
/// \param[in] u: multipoint evaluation challenge -- hash of
///            transcript (from step 4)
/// \param[in] D1: first part of batched polynomial commitment [D]_1
///            (from step 9)
/// \param[in] S_polys_at_secret_g1: permutation polynomials S
///            evaluated at the secret input (from verifier
///            preprocessed input)
/// \param[in] proof: SNARK proof produced by the prover
/// \param[in] srs: structured reference string containing also
///            circuit-specific information
///
/// OUTPUT
/// \param[out] F1: full batched polynomial commitment [F]_1
template<typename ppT>
step_ten_out_t<ppT> plonk_verifier<ppT>::step_ten(
    const step_four_out_t<ppT> &step_four_out,
    const step_nine_out_t<ppT> &step_nine_out,
    const plonk_proof<ppT> &proof,
    const verifier_preprocessed_input_t<ppT> &preprocessed_input,
    const srs<ppT> &srs)
{
    libff::G1<ppT> F1;
    Field zeta_power_n =
        libff::power(step_four_out.zeta, libff::bigint<1>(srs.num_gates + 2));
    Field zeta_power_2n = libff::power(
        step_four_out.zeta, libff::bigint<1>(2 * (srs.num_gates + 2)));
    std::vector<Field> nu_power(7);
    for (size_t i = 0; i < nu_power.size(); ++i) {
        nu_power[i] = libff::power(step_four_out.nu, libff::bigint<1>(i));
    }
    std::vector<libff::G1<ppT>> curve_points{
        proof.t_poly_at_secret_g1[lo],              // nu^0
        proof.t_poly_at_secret_g1[mid],             // nu^0
        proof.t_poly_at_secret_g1[hi],              // nu^0
        step_nine_out.D1,                           // nu^1
        proof.W_polys_blinded_at_secret_g1[a],      // nu^2
        proof.W_polys_blinded_at_secret_g1[b],      // nu^3
        proof.W_polys_blinded_at_secret_g1[c],      // nu^4
        preprocessed_input.S_polys_at_secret_g1[0], // nu^5
        preprocessed_input.S_polys_at_secret_g1[1]  // nu^6
    };
    std::vector<libff::Fr<ppT>> scalar_elements{
        Field(1),
        zeta_power_n,  // zeta^(n+2),
        zeta_power_2n, // zeta^(2*(n+2)),
        Field(1),
        nu_power[2], // nu^2,
        nu_power[3], // nu^3,
        nu_power[4], // nu^4,
        nu_power[5], // nu^5,
        nu_power[6]  // nu^6
    };
    F1 = plonk_multi_exp_G1<ppT>(curve_points, scalar_elements);
    step_ten_out_t<ppT> step_ten_out(std::move(F1));
    return step_ten_out;
}

/// Verifier Step 11: compute group-encoded batch evaluation [E]_1
///
/// INPUT
/// \param[in] nu: opening challenge -- hash of transcript (denoted by
///            v in [GWC19]) (from step 4)
/// \param[in] u: multipoint evaluation challenge -- hash of
///            transcript (from step 4)
/// \param[in] r_prime_zeta: quotient polynomial evaluation r'(zeta) =
///            r(zeta) - r0, where r0 is a constant term (from step 8)
/// \param[in] proof: SNARK proof produced by the prover
///
/// OUTPUT
/// \param[out] E1: group-encoded batch evaluation [E]_1
template<typename ppT>
step_eleven_out_t<ppT> plonk_verifier<ppT>::step_eleven(
    const step_four_out_t<ppT> &step_four_out,
    const step_eight_out_t<ppT> &step_eight_out,
    const plonk_proof<ppT> &proof)
{
    libff::G1<ppT> E1;
    std::vector<Field> nu_power(7);
    for (size_t i = 0; i < nu_power.size(); ++i) {
        nu_power[i] = libff::power(step_four_out.nu, libff::bigint<1>(i));
    }
    std::vector<libff::G1<ppT>> curve_points{libff::G1<ppT>::one()};
    std::vector<libff::Fr<ppT>> scalar_elements{
        step_eight_out.r_prime_zeta + nu_power[1] * proof.r_zeta + // v^1
        nu_power[2] * proof.a_zeta +                               // v^2
        nu_power[3] * proof.b_zeta +                               // v^3
        nu_power[4] * proof.c_zeta +                               // v^4
        nu_power[5] * proof.S_0_zeta +                             // v^5
        nu_power[6] * proof.S_1_zeta +                             // v^6
        step_four_out.u * proof.z_poly_xomega_zeta};

    E1 = plonk_multi_exp_G1<ppT>(curve_points, scalar_elements);

    step_eleven_out_t<ppT> step_eleven_out(std::move(E1));
    return step_eleven_out;
}

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
/// \param[in] zeta: evaluation challenge -- hash of transcript (from
///            step 4)
/// \param[in] u: multipoint evaluation challenge -- hash of
///            transcript (from step 4)
/// \param[in] F1: full batched polynomial commitment [F]_1 (from step
///            10)
/// \param[in] E1: group-encoded batch evaluation [E]_1 (from step 11)
/// \param[in] proof: SNARK proof produced by the prover
/// \param[in] srs: structured reference string
/// \param[in] srs: structured reference string containing also
///            circuit-specific information
///
/// OUTPUT
/// \param[out] boolean 1/0 = valid/invalid proof
template<typename ppT>
bool plonk_verifier<ppT>::step_twelve(
    const step_four_out_t<ppT> &step_four_out,
    const step_ten_out_t<ppT> &step_ten_out,
    const step_eleven_out_t<ppT> &step_eleven_out,
    const plonk_proof<ppT> &proof,
    const srs<ppT> &srs)
{
    std::vector<libff::G1<ppT>> curve_points_lhs{
        proof.W_zeta_at_secret, proof.W_zeta_omega_at_secret};
    std::vector<libff::Fr<ppT>> scalar_elements_lhs{Field(1), step_four_out.u};
    libff::G1<ppT> pairing_first_lhs =
        plonk_multi_exp_G1<ppT>(curve_points_lhs, scalar_elements_lhs);
    libff::G2<ppT> pairing_second_lhs = srs.secret_powers_g2[1];

    std::vector<libff::G1<ppT>> curve_points_rhs{
        proof.W_zeta_at_secret,
        proof.W_zeta_omega_at_secret,
        step_ten_out.F1,
        step_eleven_out.E1};
    std::vector<libff::Fr<ppT>> scalar_elements_rhs{
        // Warning! raise to the power of -1 to check e() * e()^-1 = 1
        Field(-1) * step_four_out.zeta,
        Field(-1) * step_four_out.u * step_four_out.zeta *
            srs.omega_roots[base][1],
        Field(-1) * Field(1),
        Field(-1) * Field(-1)};

    libff::G1<ppT> pairing_first_rhs =
        plonk_multi_exp_G1<ppT>(curve_points_rhs, scalar_elements_rhs);
    libff::G2<ppT> pairing_second_rhs = srs.secret_powers_g2[0];

    // TODO: move to unit test for step_twelve
#ifdef DEBUG_PLONK
    // the example class is defined specifically for the BLS12-381
    // curve, so make sure we are using this curve TODO: remove when
    // the implementation is stable and tested
    try {
        plonk_exception_assert_curve_bls12_381<ppT>();
    } catch (const std::domain_error &e) {
        std::cout << "Error: " << e.what() << "\n";
        exit(EXIT_FAILURE);
    }
    // load test vectors for debug
    plonk_example example;
    printf("[%s:%d] pairing_first_lhs\n", __FILE__, __LINE__);
    pairing_first_lhs.print();
    libff::G1<ppT> pairing_first_lhs_aff(pairing_first_lhs);
    pairing_first_lhs_aff.to_affine_coordinates();
    assert(pairing_first_lhs_aff.X == example.pairing_first_lhs[0]);
    assert(pairing_first_lhs_aff.Y == example.pairing_first_lhs[1]);

    printf("[%s:%d] pairing_first_rhs\n", __FILE__, __LINE__);
    pairing_first_rhs.print();
    libff::G1<ppT> pairing_first_rhs_aff(pairing_first_rhs);
    pairing_first_rhs_aff.to_affine_coordinates();
    assert(pairing_first_rhs_aff.X == example.pairing_first_rhs[0]);
    assert(pairing_first_rhs_aff.Y == example.pairing_first_rhs[1]);
#endif // #ifdef DEBUG_PLONK

    const libff::G1_precomp<ppT> _A = ppT::precompute_G1(pairing_first_lhs);
    const libff::G2_precomp<ppT> _B = ppT::precompute_G2(pairing_second_lhs);
    const libff::G1_precomp<ppT> _C = ppT::precompute_G1(pairing_first_rhs);
    const libff::G2_precomp<ppT> _D = ppT::precompute_G2(pairing_second_rhs);
    const libff::Fqk<ppT> miller_result =
        ppT::double_miller_loop(_A, _B, _C, _D);
    const libff::GT<ppT> result = ppT::final_exponentiation(miller_result);
    bool b_accept = (result == libff::GT<ppT>::one());
    return b_accept;
}

/// \attention The first three steps (as given in [GWC19] -- see
/// below) MUST be executed by the caller:
///
/// - Verifier Step 1: validate that elements belong to group G1
/// - Verifier Step 2: validate that elements belong to scalar field
///   Fr
/// - Verifier Step 3: validate that the public input belongs to
///   scalar field Fr
/// .
/// Therefore verification starts from Step 4 of [GWC19]
///
/// INPUT
/// \param[in] proof: SNARK proof produced by the prover
/// \param[in] srs: structured reference string containing also circuit-specific
///   information
///
/// OUTPUT
/// \param[out] boolean 1/0 = valid/invalid proof
template<typename ppT>
bool plonk_verifier<ppT>::verify_proof(
    const plonk_proof<ppT> &proof, const srs<ppT> &srs)
{
    // compute verifier preprocessed input
    const verifier_preprocessed_input_t<ppT> preprocessed_input =
        plonk_verifier::preprocessed_input(srs);

    // TODO: move to unit test for verify_proof
#ifdef DEBUG_PLONK
    // the example class is defined specifically for the BLS12-381
    // curve, so make sure we are using this curve TODO: remove when
    // the implementation is stable and tested
    try {
        plonk_exception_assert_curve_bls12_381<ppT>();
    } catch (const std::domain_error &e) {
        std::cout << "Error: " << e.what() << "\n";
        exit(EXIT_FAILURE);
    }
    // load test vector values form example for debug
    plonk_example example;
    for (int i = 0; i < (int)srs.Q_polys.size(); ++i) {
        printf("srs.Q_polys_at_secret_G1[%d] \n", i);
        preprocessed_input.Q_polys_at_secret_g1[i].print();
        libff::G1<ppT> Q_poly_at_secret_g1_i(
            preprocessed_input.Q_polys_at_secret_g1[i]);
        Q_poly_at_secret_g1_i.to_affine_coordinates();
        assert(Q_poly_at_secret_g1_i.X == example.Q_polys_at_secret_g1[i][0]);
        assert(Q_poly_at_secret_g1_i.Y == example.Q_polys_at_secret_g1[i][1]);
    }
    for (int i = 0; i < (int)srs.S_polys.size(); ++i) {
        printf("S_polys_at_secret_G1[%d] \n", i);
        preprocessed_input.S_polys_at_secret_g1[i].print();
        libff::G1<ppT> S_poly_at_secret_g1_i(
            preprocessed_input.S_polys_at_secret_g1[i]);
        S_poly_at_secret_g1_i.to_affine_coordinates();
        assert(S_poly_at_secret_g1_i.X == example.S_polys_at_secret_g1[i][0]);
        assert(S_poly_at_secret_g1_i.Y == example.S_polys_at_secret_g1[i][1]);
    }
#endif // #ifdef DEBUG_PLONK

    // Verifier Step 1: validate that elements belong to group G1
    // Verifier Step 2: validate that elements belong to scalar field Fr
    // Verifier Step 3: validate that the public input belongs to scalar field
    // Fr
    //
    // Executed by the caller

    // Verifier Step 4: compute challenges hashed transcript as in
    // prover description, from the common inputs, public input, and
    // elements of pi-SNARK (fixed to the test vectors for now)
    const step_four_out_t<ppT> step_four_out = this->step_four();

    // Verifier Step 5: compute zero polynomial evaluation
    const step_five_out_t<ppT> step_five_out =
        this->step_five(step_four_out, srs);
    // TODO: uni test for step_five
#ifdef DEBUG_PLONK
    printf("[%s:%d] zh_zeta ", __FILE__, __LINE__);
    step_five_out.zh_zeta.print();
    assert(step_five_out.zh_zeta == example.zh_zeta);
#endif // #ifdef DEBUG_PLONK

    // Verifier Step 6: Compute Lagrange polynomial evaluation L1(zeta)
    // Note: the paper counts the L-polynomials from 1; we count from 0
    const step_six_out_t<ppT> step_six_out = this->step_six(step_four_out, srs);
    // TODO: uni test for step_six
#ifdef DEBUG_PLONK
    printf("L_0_zeta ");
    step_six_out.L_0_zeta.print();
    assert(step_six_out.L_0_zeta == example.L_0_zeta);
#endif // #ifdef DEBUG_PLONK

    // Verifier Step 7: compute public input polynomial evaluation PI(zeta)
    const step_seven_out_t<ppT> step_seven_out =
        this->step_seven(step_four_out, srs);
    // TODO: uni test for step_seven
#ifdef DEBUG_PLONK
    printf("PI_zeta ");
    step_seven_out.PI_zeta.print();
    assert(step_seven_out.PI_zeta == example.PI_zeta);
#endif // #ifdef DEBUG_PLONK

    // Verifier Step 8: compute quotient polynomial evaluation
    // r'(zeta) = r(zeta) - r0, where r0 is a constant term
    const step_eight_out_t<ppT> step_eight_out = this->step_eight(
        step_four_out, step_five_out, step_six_out, step_seven_out, proof);
    // TODO: uni test for step_eight
#ifdef DEBUG_PLONK
    assert(step_eight_out.r_prime_zeta == example.r_prime_zeta);
#endif // #ifdef DEBUG_PLONK

    // Verifier Step 9: compute first part of batched polynomial
    // commitment [D]_1
    step_nine_out_t<ppT> step_nine_out = this->step_nine(
        step_four_out, step_six_out, proof, preprocessed_input, srs);
    // TODO: uni test for step_nine
#ifdef DEBUG_PLONK
    step_nine_out.D1.print();
    libff::G1<ppT> D1_aff(step_nine_out.D1);
    D1_aff.to_affine_coordinates();
    assert(D1_aff.X == example.D1[0]);
    assert(D1_aff.Y == example.D1[1]);
#endif // #ifdef DEBUG_PLONK

    // Verifier Step 10: compute full batched polynomial commitment
    // [F]_1
    step_ten_out_t<ppT> step_ten_out = this->step_ten(
        step_four_out, step_nine_out, proof, preprocessed_input, srs);
    // TODO: uni test for step_ten
#ifdef DEBUG_PLONK
    printf("[%s:%d] F1\n", __FILE__, __LINE__);
    step_ten_out.F1.print();
    libff::G1<ppT> F1_aff(step_ten_out.F1);
    F1_aff.to_affine_coordinates();
    assert(F1_aff.X == example.F1[0]);
    assert(F1_aff.Y == example.F1[1]);
#endif // #ifdef DEBUG_PLONK

    // Verifier Step 11: compute group-encoded batch evaluation [E]_1
    const step_eleven_out_t<ppT> step_eleven_out =
        this->step_eleven(step_four_out, step_eight_out, proof);
    // TODO: uni test for step_eleven
#ifdef DEBUG_PLONK
    printf("[%s:%d] E1\n", __FILE__, __LINE__);
    step_eleven_out.E1.print();
    libff::G1<ppT> E1_aff(step_eleven_out.E1);
    E1_aff.to_affine_coordinates();
    assert(E1_aff.X == example.E1[0]);
    assert(E1_aff.Y == example.E1[1]);
#endif // #ifdef DEBUG_PLONK

    // Verifier Step 12: batch validate all evaluations (check
    // pairing)
    bool b_accept = this->step_twelve(
        step_four_out, step_ten_out, step_eleven_out, proof, srs);
    return b_accept;
}

} // namespace libsnark

#endif // LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_VERIFIER_TCC_