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

} // namespace libsnark

#include "libsnark/gadgetlib1/gadgets/verifiers/kzg10_batched_verifier_gadget.tcc"

#endif // LIBSNARK_GADGETLIB1_GADGETS_VERIFIERS_KZG10_BATCHED_VERIFIER_GADGET_HPP_
