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
        protoboard<libff::Fr<ppT>> &pb, const std::string &annotation_prefix);

    void generate_r1cs_witness(
        const typename kzg10_batched_2_point<
            other_curve<ppT>>::evaluation_witness &eval_witness);
};

/// Given a value `gamma` and a vector `points` of `n` group variables, compute:
///
///   result = \sum_{i=0}^{n-1} gamma^i * points[n]
///
/// by computing:
///
///   intermediate[0] = gamma * points[n] + points[n-1]
///   intermediate[i>0] = gamma * intermediate[i-1] + points[n-1-i]
///   result = gamma * intermediate[n-2] + points[0]
template<typename ppT, size_t n>
class kzg10_batched_compute_gamma_powers_times_points : gadget<libff::Fr<ppT>>
{
public:
    // Full calculation is as follows:
    //
    //   intermediate_mul[0] = gamma * points[n-1]
    //   intermediate_sum[0] = intermediate_mul[0] + points[n-2]
    //   intermediate_mul[i=1..n-2] = gamma * intermediate_sum[i-1]
    //   intermediate_sum[i=1..n-3] = intermediate_mul[i] + points[n-2-i]
    //   intermediate_sum[n-2] = result = intermediate_mul[n-2] + points[0]
    //
    // so intermediate_mul.size() = intermediate_sum.size() = n-1
    std::vector<G1_mul_by_scalar_gadget<ppT>> intermediate_mul;
    std::vector<G1_add_variable_and_variable_or_identity_gadget<ppT>>
        intermediate_sum;

    kzg10_batched_compute_gamma_powers_times_points(
        protoboard<libff::Fr<ppT>> &pb,
        const pb_linear_combination<libff::Fr<ppT>> &gamma,
        const std::vector<G1_variable<ppT>> &points,
        G1_variable<ppT> &result,
        const std::string &annotation_prefix);

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

/// Given an array of commitments, and array of evaluations, and some field
/// element gamma, compute terms of the form:
///
///   \sum_{i=1}^{t_1} \gamma^{i-1} (cm_i - [s_i]_1)
///
/// The `num_entries` parameter here is intended to reflect the fact that this
/// must be statically defined, although all internal structures are currently
/// dynamic. This also allows specialization for the case of 2 entries.
template<typename ppT, size_t num_entries>
class kzg10_batched_compute_commit_minus_eval_sum : gadget<libff::Fr<ppT>>
{
    static_assert(num_entries > 2, "num_entries must be greater that 2");

public:
    using Field = libff::Fr<ppT>;

    // These are the negative encoded evaluations:
    //
    //   encoded_evals[i] = G1_mul(evals[i], -G1::one());
    std::vector<G1_variable_or_identity<ppT>> encoded_evals;
    std::vector<G1_mul_by_scalar_gadget<ppT>> compute_encoded_evals;

    // Negative evals are added to commits to compute cm_i - [s_i]_1:
    //
    //   commit_minus_encoded_eval[i] = cm_i - [s_i]_1
    std::vector<G1_variable<ppT>> commit_minus_encoded_eval;
    std::vector<G1_add_variable_and_variable_or_identity_gadget<ppT>>
        compute_commit_minus_encoded_eval;

    // result = sum_{i=0}^{n-1} gamma^i commit_minus_encoded_eval[i]
    kzg10_batched_compute_gamma_powers_times_points<ppT, num_entries>
        compute_gamma_power_times_commit_minus_encoded_eval;

    kzg10_batched_compute_commit_minus_eval_sum(
        protoboard<libff::Fr<ppT>> &pb,
        const pb_linear_combination<libff::Fr<ppT>> &gamma,
        const std::vector<kzg10_commitment_variable<ppT>> &commitments,
        const pb_linear_combination_array<libff::Fr<ppT>> &evals,
        G1_variable<ppT> &result,
        const std::string &annotation_prefix);

    void generate_r1cs_constraints();
    void generate_r1cs_witness();

    const G1_variable<ppT> &result() const;
};

// Specialization for num_entries == 2 (in which we do not need to compute
// further powers of gamma). This simplifies the generic (num_entries >= 3)
// version, since it does not need to account for special cases.
template<typename ppT>
class kzg10_batched_compute_commit_minus_eval_sum<ppT, 2>
    : gadget<libff::Fr<ppT>>
{
public:
    // encoded_evals[i] = evals[i] * -G1::one()
    std::vector<G1_mul_by_scalar_gadget<ppT>> compute_encoded_evals;

    // cm_minus_encoded_eval[i] = commits[i] - encoded_evals[i]
    std::vector<G1_add_variable_and_variable_or_identity_gadget<ppT>>
        compute_cm_minus_eval;

    // gamma_term = gamma * commit_minus_encoded_eval[1]
    // return = gamma_term + (commit_minus_encoded_eval[0]
    std::shared_ptr<G1_mul_by_scalar_gadget<ppT>> compute_gamma_term;
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

    const G1_variable<ppT> &result() const;
};

/// Gadget version of the native KZG10 batched verifier in
/// libsnark/polynomial_commitments/kzg10_batched.hpp.
///
/// Each polynomials can be evaluated at 1 of 2 points. `polyomials_1`
/// determines the number of polynomials evaluated at the first point `z_1`,
/// and `polynomials_2` determines the number to be evaluated at the second
/// point `z_2`. Hence these also determine the number of commitments and
/// evaluation points.
///
/// The number of polynomials inn each group is intentionally a template
/// parameter, reflecting the fact that it must be statically defined for a
/// given circuit.
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
    // could use statically sized containers here. For now, the interfaces and
    // initialization make this a bit inconvenient (requiring default
    // constructors for the contained types).
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
