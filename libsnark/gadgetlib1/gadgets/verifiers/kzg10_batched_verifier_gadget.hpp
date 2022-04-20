/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_GADGETLIB1_GADGETS_VERIFIERS_KZG10_BATCHED_VERIFIER_GADGET_HPP_
#define LIBSNARK_GADGETLIB1_GADGETS_VERIFIERS_KZG10_BATCHED_VERIFIER_GADGET_HPP_

#include "libsnark/gadgetlib1/gadgets/verifiers/kzg10_verifier_gadget.hpp"
#include "libsnark/polynomial_commitments/kzg10_batched.hpp"

namespace libsnark
{

// TODO: create an "evaluations variable object"

template<typename ppT> class kzg10_batched_witness_variable
{
public:
    kzg10_witness_variable<ppT> W_1;
    kzg10_witness_variable<ppT> W_2;

    kzg10_batched_witness_variable(
        protoboard<libff::Fr<ppT>> &pb, const std::string &annotationn_prefix);

    void generate_r1cs_witness(
        const typename kzg10_batched_2_point<
            other_curve<ppT>>::evaluation_witness &eval_witness);
};

#if 0
/// Use Horners method to evaluate a polynomial expression using an optimal
/// number of conditions. This gadget is intended to be called with x and
/// coeffs already allocated, and uses references rather than taking a copy.
template<typename ppT, size_t num_entries>
class evaluate_polynomial_horner : gadget<libff::Fr<ppT>>
{
public:
    static_assert(num_entries > 2, "fewer than 3 entries not supported");

    // TODO: switch to shared_pointer here for consistency.

    // Compute \sum_{i=0}^{n-1} coeffs[i] x^i as:
    //
    //   x * ( ... (x * (x * coeffs[n-1] + coeff[n-2]) + coeffs[n-3]) ... )
    //
    // where intermediate sum_terms are:
    //
    //   sum_term[0] = x * coeffs[n-1] + coeffs[n-2]
    //   sum_term[0 < i < n-2] = x * sum_term[i-1] + coeffs[n-(i+2)]
    //   result = x * sum_term[n-3] + coeffs[0]
    //
    const pb_linear_combination<libff::Fr<ppT>> &x;
    const pb_linear_combination_array<libff::Fr<ppT>> &coeffs;
    pb_variable_array<libff::Fr<ppT>> sum_terms;
    pb_variable<libff::Fr<ppT>> result;

    evaluate_polynomial_horner(
        protoboard<libff::Fr<ppT>> &pb,
        const pb_linear_combination<libff::Fr<ppT>> &x,
        const pb_linear_combination_array<libff::Fr<ppT>> &coeffs,
        pb_variable<libff::Fr<ppT>> &result,
        const std::string &annotation_prefix);

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};
#endif

/// Given evals[0...n-1], a value gamma and gamma_powers[0...n-3] where
///
///   gamma_powers[i] = gamma^[i+2]
///
/// compute:
///
///   result = \sum_{i=0}^{n-1} gamma^i * evals[i]
template<typename ppT, size_t n>
class kzg10_batched_compute_gamma_eval_sum : gadget<libff::Fr<ppT>>
{
public:
    static_assert(n > 2, "fewer than 3 entries not supported");

    // TODO: We take a copy of the array (in case the array being passed in is
    // actually a pb_variable_array, which then creates a transient
    // pb_linear_combination_array via the conversion constructor).

    // TODO: if we switch to references to fixed-size arrays, this may not be
    // necesary.

    const pb_linear_combination<libff::Fr<ppT>> gamma;
    const pb_linear_combination_array<libff::Fr<ppT>> evals;

    // gamma_powers[i] = gamma^(i+2)
    pb_variable_array<libff::Fr<ppT>> gamma_powers;

    // sum_terms[0] = gamma * evals[1] + evals[0]
    // sum_terms[1] = gamma_powers[0] * evals[2] + sum_terms[0]
    // ...
    // sum_terms[i=2..n-3] = gamma_powers[i-1] * evals[i+1] + sum_terms[i-1]
    pb_variable_array<libff::Fr<ppT>> sum_terms;

    // result = gamma_powers[n-3] * evals[n-1] + sum_terms[n-3]
    pb_variable<libff::Fr<ppT>> result;

    kzg10_batched_compute_gamma_eval_sum(
        protoboard<libff::Fr<ppT>> &pb,
        const pb_linear_combination<libff::Fr<ppT>> &gamma,
        const pb_variable_array<libff::Fr<ppT>> &gamma_powers,
        const pb_linear_combination_array<libff::Fr<ppT>> &evals,
        pb_variable<libff::Fr<ppT>> &result,
        const std::string &annotation_prefix);

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

/// Given an array of commitments, and array of evaluations, and some field
/// element gamma, compute terms of the form:
///
///   \sum_{i=1}^{t-1} \gamma^{i-1} (cm_i - [s_i]_1)
///
template<typename ppT, size_t num_entries>
class kzg10_batched_compute_commit_minus_eval_sum : gadget<libff::Fr<ppT>>
{
    static_assert(num_entries > 2, "num_entries must be greater that 2");

public:
    using Field = libff::Fr<ppT>;

    pb_linear_combination_array<Field> evals;

    // Array does not contain 1 or gamma, so:
    //
    //   gamma_powers[i] = gamma * gamma_powers[i - 1]
    //                   = gamma^{i + 2}
    //
    // for i = 0 ... num_entries - 3 (len = num_entries - 2)
    pb_linear_combination<Field> gamma;
    pb_variable_array<Field> gamma_powers;

    // gamma_power_times_commit[0] = gamma * commits[1]
    // gamma_power_times_commit[i>0] = gamma^{i+1} * commits[i + 1]
    //
    // for i = 0 ... num_entries - 2 (len = num_entries - 1)
    std::vector<G1_variable_or_identity<ppT>> gamma_power_times_commit;
    std::vector<G1_mul_by_scalar_gadget<ppT>> compute_gamma_power_times_commit;

    // sum_gamma_power_times_commit[0] = cm_0 + gamma_power_times_commit[0]
    // sum_gamma_power_times_commit[i>0] =
    //   sum_gamma_power_times_commit[i] + gamma_power_times_commit[i]
    //
    // for i = 0 ... num_entries - 2 (len = num_entries - 1)
    std::vector<G1_variable<ppT>> sum_gamma_power_times_commit;
    std::vector<G1_add_variable_and_variable_or_identity_gadget<ppT>>
        compute_sum_gamma_power_times_commit;

    // Compute sum_gamma_power_times_commit = \sum_i \gamma^{i-1} s_i
    pb_variable<libff::Fr<ppT>> sum_gamma_power_times_eval;
    kzg10_batched_compute_gamma_eval_sum<ppT, num_entries>
        compute_sum_gamma_power_times_eval;

    // encoded_sum_gamma_power_times_eval = [ sum_gamma_power_times_eval ]_1
    G1_variable_or_identity<ppT> encoded_sum_gamma_power_times_eval;
    G1_mul_by_scalar_gadget<ppT> compute_encoded_sum_gamma_power_times_eval;

    // result =
    //   sum_gamma_power_times_commit - encoded_sum_gamma_power_times_eval;
    G1_add_variable_and_variable_or_identity_gadget<ppT> compute_result;

    kzg10_batched_compute_commit_minus_eval_sum(
        protoboard<libff::Fr<ppT>> &pb,
        const pb_linear_combination<libff::Fr<ppT>> &gamma,
        const std::vector<kzg10_commitment_variable<ppT>> &commits,
        const pb_linear_combination_array<libff::Fr<ppT>> &evals,
        G1_variable<ppT> &result,
        const std::string &annotation_prefix);

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

// Specialization for num_entries == 2
template<typename ppT>
class kzg10_batched_compute_commit_minus_eval_sum<ppT, 2>
    : gadget<libff::Fr<ppT>>
{
public:
    std::vector<G1_variable_or_identity<ppT>> encoded_evals;
    std::vector<G1_mul_by_scalar_gadget<ppT>> compute_encoded_evals;

    // commit_minus_encoded_eval[i] = cm_i - [s_i]_1
    std::vector<G1_variable<ppT>> commit_minus_encoded_eval;
    std::vector<G1_add_variable_and_variable_or_identity_gadget<ppT>>
        compute_commit_minus_encoded_eval;

    // gamma * (cm_1 - [s_1])
    G1_variable_or_identity<ppT> gamma_times_commit_1_minus_encoded_eval_1;
    std::shared_ptr<G1_mul_by_scalar_gadget<ppT>>
        compute_gamma_times_commit_1_minus_encoded_eval_1;

    // result = (cm_0 - [s_0]) + gamma * (cm_1 - [s_1])
    std::shared_ptr<G1_add_variable_and_variable_or_identity_gadget<ppT>>
        compute_result;

    kzg10_batched_compute_commit_minus_eval_sum(
        protoboard<libff::Fr<ppT>> &pb,
        pb_linear_combination<libff::Fr<ppT>> gamma,
        const std::vector<kzg10_commitment_variable<ppT>> commitments,
        const pb_linear_combination_array<libff::Fr<ppT>> &evals,
        G1_variable<ppT> &result,
        const std::string &annotation_prefix);

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

//   F = \sum_{i=1}^{t1} \gamma_1^{i-1} (cm_1[i] - [s_1[i]]_1) +      (G)
//       r \sum_{i=1}^{t2} \gamma_2^{i-1} (cm_2[i] - [s_2[i]]_1)      (H)

/// Gadget version of the native KZG10 batched verifier in
/// libsnark/polynomial_commitments/kzg10_batched.hpp.
///
/// Each polynomials can be evaluated at 1 of 2 points. `polyomials_1`
/// determines the number of polynomials evaluated at the first point `z_1`,
/// and `polynomials_2` determines the number to be evaluated at the second
/// point `z_2`. Hence these also determine the number of commitments and
/// evaluation points.
template<typename ppT, size_t num_polyomials_1, size_t num_polyomials_2>
class kzg10_batched_verifier_gadget : public gadget<libff::Fr<ppT>>
{
public:
    using Field = libff::Fr<ppT>;

    // Matching the native calculations in kzg10_batched.tcc, compute:
    //
    //   F = \sum_{i=1}^{t1} \gamma_1^{i-1} (cm_1[i] - [s_1[i]]_1) +
    //       r \sum_{i=1}^{t2} \gamma_2^{i-1} (cm_2[i] - [s_2[i]]_1)
    //     = G + r * H
    //
    // where:
    //
    //   G = \sum_{i=1}^{t1} \gamma_1^{i-1} (cm_1[i] - [s_1[i]]_1)
    //   H = \sum_{i=1}^{t2} \gamma_2^{i-1} (cm_2[i] - [s_2[i]]_1)
    G1_variable<ppT> G;
    kzg10_batched_compute_commit_minus_eval_sum<ppT, num_polyomials_1>
        compute_G;
    G1_variable<ppT> H;
    kzg10_batched_compute_commit_minus_eval_sum<ppT, num_polyomials_2>
        compute_H;
    G1_variable_or_identity<ppT> rH;
    G1_mul_by_scalar_gadget<ppT> compute_rH;
    G1_variable<ppT> F;
    G1_add_variable_and_variable_or_identity_gadget<ppT> compute_F;

    // Expression to check is:
    //   e(W_1 + r * W_2, srs.alpha_g2) *
    //     e(F + z_1 * W_1 + r * z_2 * W_2, -[1]_2)
    //   = 1
    //   = e(A, srs.alpha_g2) * e(B, -[1]_2)
    //
    // where
    //   A = W_1 + r * W_2
    //   B = F + z_1 * W_1 + r * z_2 * W_2

    G1_variable_or_identity<ppT> r_times_W_2;
    G1_mul_by_scalar_gadget<ppT> compute_r_times_W_2;

    G1_variable<ppT> A;
    G1_add_variable_and_variable_or_identity_gadget<ppT> compute_A;

    G1_variable_or_identity<ppT> r_times_z_2_times_W_2;
    G1_variable_or_identity_mul_by_scalar_gadget<ppT>
        compute_r_times_z_2_times_W_2;

    G1_variable_or_identity<ppT> z_1_times_W_1;
    G1_mul_by_scalar_gadget<ppT> compute_z_1_times_W_1;

    G1_variable<ppT> F_plus_z_1_times_W_1;
    G1_add_variable_and_variable_or_identity_gadget<ppT>
        compute_F_plus_z_1_times_W_1;

    G1_variable<ppT> B;
    G1_add_variable_and_variable_or_identity_gadget<ppT> compute_B;

    kzg10_pairing_check_gadget<ppT> pairing_check;

    // TODO: Since polyomials_1 and polyomials_2 are statically defined, we
    // could use arrays here.
    kzg10_batched_verifier_gadget(
        protoboard<libff::Fr<ppT>> &pb,
        pb_linear_combination<libff::Fr<ppT>> z_1,
        pb_linear_combination<libff::Fr<ppT>> z_2,
        const pb_linear_combination_array<libff::Fr<ppT>> &poly_evals_1,
        const pb_linear_combination_array<libff::Fr<ppT>> &poly_evals_2,
        const kzg10_srs_variable<ppT> &srs,
        pb_linear_combination<libff::Fr<ppT>> gamma_1,
        pb_linear_combination<libff::Fr<ppT>> gamma_2,
        const kzg10_batched_witness_variable<ppT> &eval_witness,
        const std::vector<kzg10_commitment_variable<ppT>> &commitments_1,
        const std::vector<kzg10_commitment_variable<ppT>> &commitments_2,
        pb_linear_combination<libff::Fr<ppT>> r,
        pb_variable<libff::Fr<ppT>> result,
        const std::string &annotation_prefix);

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

} // namespace libsnark

#include "libsnark/gadgetlib1/gadgets/verifiers/kzg10_batched_verifier_gadget.tcc"

#endif // LIBSNARK_GADGETLIB1_GADGETS_VERIFIERS_KZG10_BATCHED_VERIFIER_GADGET_HPP_
