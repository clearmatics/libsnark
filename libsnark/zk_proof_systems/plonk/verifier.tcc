/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_VERIFIER_TCC_
#define LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_VERIFIER_TCC_

// Implementation of Verifier interfaces for a ppzkSNARK for Plonk. See
// verifier.hpp .

namespace libsnark
{

template<typename ppT>
verifier_preprocessed_input_t<ppT>::verifier_preprocessed_input_t(
    std::vector<libff::G1<ppT>> &&Q_polys_at_secret_g1,
    std::vector<libff::G1<ppT>> &&S_polys_at_secret_g1)
    : Q_polys_at_secret_g1(Q_polys_at_secret_g1)
    , S_polys_at_secret_g1(S_polys_at_secret_g1)
{
}

template<typename ppT, class transcript_hasher>
verifier_preprocessed_input_t<ppT> plonk_verifier<ppT, transcript_hasher>::
    preprocessed_input(const srs<ppT> &srs)
{
    std::vector<libff::G1<ppT>> Q_polys_at_secret_g1;
    Q_polys_at_secret_g1.resize(srs.Q_polys.size());
    plonk_evaluate_polys_at_secret_G1<ppT>(
        srs.secret_powers_g1, srs.Q_polys, Q_polys_at_secret_g1);

    std::vector<libff::G1<ppT>> S_polys_at_secret_g1;
    S_polys_at_secret_g1.resize(srs.S_polys.size());
    plonk_evaluate_polys_at_secret_G1<ppT>(
        srs.secret_powers_g1, srs.S_polys, S_polys_at_secret_g1);

    verifier_preprocessed_input_t<ppT> preprocessed_input(
        std::move(Q_polys_at_secret_g1), std::move(S_polys_at_secret_g1));
    return preprocessed_input;
}

template<typename ppT>
step_four_out_t<ppT>::step_four_out_t(
    const libff::Fr<ppT> &beta,
    const libff::Fr<ppT> &gamma,
    const libff::Fr<ppT> &alpha,
    const libff::Fr<ppT> &zeta,
    const libff::Fr<ppT> &nu,
    const libff::Fr<ppT> &u)
    : beta(beta), gamma(gamma), alpha(alpha), zeta(zeta), nu(nu), u(u)
{
}

template<typename ppT, class transcript_hasher>
step_four_out_t<ppT> plonk_verifier<ppT, transcript_hasher>::step_four(
    const plonk_proof<ppT> &proof, transcript_hasher &hasher)
{
    // Add outputs from Round 1 to the hash buffer.
    hasher.add_element(proof.W_polys_blinded_at_secret_g1[a]);
    hasher.add_element(proof.W_polys_blinded_at_secret_g1[b]);
    hasher.add_element(proof.W_polys_blinded_at_secret_g1[c]);
    // - beta: permutation challenge - hashes of transcript of round
    // 1 (\attention the original protocl appends a 0, while we don't
    // append anyhting to avoid making a copy of the buffer)
    libff::Fr<ppT> beta = hasher.get_hash();
    // - gamma: permutation challenge - hashes of transcript of round
    // 1 with 1 appended
    hasher.add_element(libff::Fr<ppT>::one());
    libff::Fr<ppT> gamma = hasher.get_hash();

    // Add outputs from Round 2 to the hash buffer.
    hasher.add_element(proof.z_poly_at_secret_g1);
    // - alpha: quotient challenge - hash of transcript of rounds 1,2
    libff::Fr<ppT> alpha = hasher.get_hash();

    // Add outputs from Round 3 to the hash buffer.
    hasher.add_element(proof.t_poly_at_secret_g1[lo]);
    hasher.add_element(proof.t_poly_at_secret_g1[mid]);
    hasher.add_element(proof.t_poly_at_secret_g1[hi]);
    // - zeta: evaluation challenge - hash of transcriptof rounds
    // - 1,2,3
    libff::Fr<ppT> zeta = hasher.get_hash();

    // Add outputs from Round 4 to the hash buffer.
    hasher.add_element(proof.a_zeta);
    hasher.add_element(proof.b_zeta);
    hasher.add_element(proof.c_zeta);
    hasher.add_element(proof.S_0_zeta);
    hasher.add_element(proof.S_1_zeta);
    hasher.add_element(proof.z_poly_xomega_zeta);
    // - nu: opening challenge -- hash of transcript (denoted by v in
    //   [GWC19])
    libff::Fr<ppT> nu = hasher.get_hash();

    // Add outputs from Round 5 to the hash buffer.
    hasher.add_element(proof.r_zeta);
    hasher.add_element(proof.W_zeta_at_secret);
    hasher.add_element(proof.W_zeta_omega_at_secret);
    // u: multipoint evaluation challenge -- hash of transcript from
    // rounds 1,2,3,4,5
    libff::Fr<ppT> u = hasher.get_hash();

    // step 4 output
    step_four_out_t<ppT> step_four_out(beta, gamma, alpha, zeta, nu, u);

    return step_four_out;
}

template<typename ppT>
step_five_out_t<ppT>::step_five_out_t(libff::Fr<ppT> &&zh_zeta)
    : zh_zeta(zh_zeta)
{
}

template<typename ppT, class transcript_hasher>
step_five_out_t<ppT> plonk_verifier<ppT, transcript_hasher>::step_five(
    const step_four_out_t<ppT> &step_four_out,
    std::shared_ptr<libfqfft::evaluation_domain<libff::Fr<ppT>>> domain)
{
    libff::Fr<ppT> zh_zeta =
        domain->compute_vanishing_polynomial(step_four_out.zeta);
    step_five_out_t<ppT> step_five_out(std::move(zh_zeta));
    return step_five_out;
}

/// struct step_six_out_t constructor
template<typename ppT>
step_six_out_t<ppT>::step_six_out_t(libff::Fr<ppT> &&L_0_zeta)
    : L_0_zeta(L_0_zeta)
{
}

template<typename ppT, class transcript_hasher>
step_six_out_t<ppT> plonk_verifier<ppT, transcript_hasher>::step_six(
    const step_four_out_t<ppT> &step_four_out, const srs<ppT> &srs)
{
    libff::Fr<ppT> L_0_zeta = libfqfft::evaluate_polynomial<Field>(
        srs.L_basis_zero.size(), srs.L_basis_zero, step_four_out.zeta);
    step_six_out_t<ppT> step_six_out(std::move(L_0_zeta));
    return step_six_out;
}

template<typename ppT>
step_seven_out_t<ppT>::step_seven_out_t(libff::Fr<ppT> &&PI_zeta)
    : PI_zeta(PI_zeta)
{
}

template<typename ppT, class transcript_hasher>
step_seven_out_t<ppT> plonk_verifier<ppT, transcript_hasher>::step_seven(
    const step_four_out_t<ppT> &step_four_out,
    const std::vector<Field> &PI_value_list,
    const srs<ppT> &srs)
{
    std::shared_ptr<libfqfft::evaluation_domain<Field>> domain =
        libfqfft::get_evaluation_domain<Field>(srs.num_gates);
    // Construct the PI polynomial from the vector of PI values
    // (received as input to the verifier) and the PI wire indices
    // (stored in the srs).
    std::vector<Field> PI_points(srs.num_gates, Field(0));
    for (size_t i = 0; i < PI_value_list.size(); i++) {
        size_t PI_coordinate_x = srs.PI_wire_indices[i] % srs.num_gates;
        PI_points[PI_coordinate_x] = Field(-PI_value_list[i]);
    }

    // Compute PI polynomial.
    polynomial<Field> PI_poly;
    plonk_compute_public_input_polynomial(PI_points, PI_poly, domain);

    libff::Fr<ppT> PI_zeta;
    PI_zeta = libfqfft::evaluate_polynomial<Field>(
        PI_poly.size(), PI_poly, step_four_out.zeta);
    step_seven_out_t<ppT> step_seven_out(std::move(PI_zeta));
    return step_seven_out;
}

template<typename ppT>
step_eight_out_t<ppT>::step_eight_out_t(libff::Fr<ppT> &&r_prime_zeta)
    : r_prime_zeta(r_prime_zeta)
{
}

template<typename ppT, class transcript_hasher>
step_eight_out_t<ppT> plonk_verifier<ppT, transcript_hasher>::step_eight(
    const step_four_out_t<ppT> &step_four_out,
    const step_five_out_t<ppT> &step_five_out,
    const step_six_out_t<ppT> &step_six_out,
    const step_seven_out_t<ppT> &step_seven_out,
    const plonk_proof<ppT> &proof)
{
    libff::Fr<ppT> r_prime_zeta;
    Field alpha_power2 = libff::power(step_four_out.alpha, libff::bigint<1>(2));

    // Compute polynomial r'(zeta) = r(zeta) - r_0 .
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

template<typename ppT>
step_nine_out_t<ppT>::step_nine_out_t(libff::G1<ppT> &&D1) : D1(D1)
{
}

template<typename ppT, class transcript_hasher>
step_nine_out_t<ppT> plonk_verifier<ppT, transcript_hasher>::step_nine(
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

    // Compute D1_part[0]:
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

    // Compute D1_part[1]:
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

    // Compute D1_part[2]:
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

    // Compute D1 = D1_part[0] + D1_part[1] + D1_part[2] .
    D1 = D1_part[0] + D1_part[1] + D1_part[2];

    step_nine_out_t<ppT> step_nine_out(std::move(D1));
    return step_nine_out;
}

template<typename ppT>
step_ten_out_t<ppT>::step_ten_out_t(libff::G1<ppT> &&F1) : F1(F1)
{
}

template<typename ppT, class transcript_hasher>
step_ten_out_t<ppT> plonk_verifier<ppT, transcript_hasher>::step_ten(
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

template<typename ppT>
step_eleven_out_t<ppT>::step_eleven_out_t(libff::G1<ppT> &&E1) : E1(E1)
{
}

template<typename ppT, class transcript_hasher>
step_eleven_out_t<ppT> plonk_verifier<ppT, transcript_hasher>::step_eleven(
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

template<typename ppT, class transcript_hasher>
bool plonk_verifier<ppT, transcript_hasher>::step_twelve(
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

template<typename ppT, class transcript_hasher>
bool plonk_verifier<ppT, transcript_hasher>::verify_proof(
    const plonk_proof<ppT> &proof,
    const srs<ppT> &srs,
    const std::vector<Field> &PI_value_list,
    transcript_hasher &hasher)
{
    std::shared_ptr<libfqfft::evaluation_domain<libff::Fr<ppT>>> domain =
        libfqfft::get_evaluation_domain<libff::Fr<ppT>>(srs.num_gates);

    // Compute verifier preprocessed input.
    const verifier_preprocessed_input_t<ppT> preprocessed_input =
        plonk_verifier::preprocessed_input(srs);

    // Verifier Step 1: validate that elements belong to group G1
    // Verifier Step 2: validate that elements belong to scalar field Fr
    // Verifier Step 3: validate that the public input belongs to scalar field
    // Fr
    //
    // Executed by the caller

    // Verifier Step 4: compute challenges hashed transcript as in
    // prover description, from the common inputs, public input, and
    // elements of pi-SNARK (fixed to the test vectors for now)
    const step_four_out_t<ppT> step_four_out = this->step_four(proof, hasher);

    // Verifier Step 5: compute zero polynomial evaluation
    const step_five_out_t<ppT> step_five_out =
        this->step_five(step_four_out, domain);

    // Verifier Step 6: Compute Lagrange polynomial evaluation L1(zeta)
    // Note: the paper counts the L-polynomials from 1; we count from 0
    const step_six_out_t<ppT> step_six_out = this->step_six(step_four_out, srs);

    // Verifier Step 7: compute public input polynomial evaluation PI(zeta)
    const step_seven_out_t<ppT> step_seven_out =
        this->step_seven(step_four_out, PI_value_list, srs);

    // Verifier Step 8: compute quotient polynomial evaluation
    // r'(zeta) = r(zeta) - r0, where r0 is a constant term
    const step_eight_out_t<ppT> step_eight_out = this->step_eight(
        step_four_out, step_five_out, step_six_out, step_seven_out, proof);

    // Verifier Step 9: compute first part of batched polynomial
    // commitment [D]_1
    step_nine_out_t<ppT> step_nine_out = this->step_nine(
        step_four_out, step_six_out, proof, preprocessed_input, srs);

    // Verifier Step 10: compute full batched polynomial commitment
    // [F]_1
    step_ten_out_t<ppT> step_ten_out = this->step_ten(
        step_four_out, step_nine_out, proof, preprocessed_input, srs);

    // Verifier Step 11: compute group-encoded batch evaluation [E]_1
    const step_eleven_out_t<ppT> step_eleven_out =
        this->step_eleven(step_four_out, step_eight_out, proof);

    // Verifier Step 12: batch validate all evaluations (check
    // pairing)
    bool b_accept = this->step_twelve(
        step_four_out, step_ten_out, step_eleven_out, proof, srs);
    return b_accept;
}

} // namespace libsnark

#endif // LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_VERIFIER_TCC_
