/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_GADGETLIB1_GADGETS_VERIFIERS_KZG10_BATCHED_VERIFIER_GADGET_TCC_
#define LIBSNARK_GADGETLIB1_GADGETS_VERIFIERS_KZG10_BATCHED_VERIFIER_GADGET_TCC_

#include "libsnark/gadgetlib1/gadgets/verifiers/kzg10_verifier_gadget.hpp"

namespace libsnark
{

namespace internal
{

// Convenience function, since this scheme uses so many vectors of complex
// variables.
template<typename ppT, typename GroupVarT>
std::vector<GroupVarT> allocate_variable_array(
    protoboard<libff::Fr<ppT>> &pb,
    const size_t length,
    const std::string &annotation_prefix)
{
    std::vector<GroupVarT> variables;
    variables.reserve(length);
    for (size_t i = 0; i < length; ++i) {
        variables.emplace_back(pb, FMT(annotation_prefix, "[%zu]", i));
    }

    return variables;
}

} // namespace internal

//
// kzg10_batched_compute_commit_minus_eval_sum
//

// specialization for 2 entries
template<typename ppT>
kzg10_batched_compute_commit_minus_eval_sum<ppT, 2>::
    kzg10_batched_compute_commit_minus_eval_sum(
        protoboard<libff::Fr<ppT>> &pb,
        pb_linear_combination<libff::Fr<ppT>> gamma,
        const std::vector<kzg10_commitment_variable<ppT>> commitments,
        const pb_linear_combination_array<libff::Fr<ppT>> &evals,
        G1_variable<ppT> &result,
        const std::string &annotation_prefix)
    : gadget<libff::Fr<ppT>>(pb, annotation_prefix)
    , encoded_evals(
          internal::allocate_variable_array<ppT, G1_variable_or_identity<ppT>>(
              pb, 2, FMT(annotation_prefix, " encoded_evals")))
    , commit_minus_encoded_eval(
          internal::allocate_variable_array<ppT, G1_variable<ppT>>(
              pb, 2, FMT(annotation_prefix, " commit_minus_encoded_eval")))
    , gamma_times_commit_1_minus_encoded_eval_1(
          pb,
          FMT(annotation_prefix, " gamma_times_commit_1_minus_encoded_eval_1"))
{
    compute_encoded_evals.emplace_back(
        pb,
        evals[0],
        G1_variable<ppT>(
            pb,
            -libff::G1<other_curve<ppT>>::one(),
            FMT(annotation_prefix, " g1_one_0")),
        encoded_evals[0],
        FMT(annotation_prefix, " compute_encoded_evals[0]"));
    compute_encoded_evals.emplace_back(
        pb,
        evals[1],
        G1_variable<ppT>(
            pb,
            -libff::G1<other_curve<ppT>>::one(),
            FMT(annotation_prefix, " g1_one_1")),
        encoded_evals[1],
        FMT(annotation_prefix, " compute_encoded_evals[1]"));

    compute_commit_minus_encoded_eval.emplace_back(
        pb,
        encoded_evals[0],
        commitments[0],
        commit_minus_encoded_eval[0],
        FMT(annotation_prefix, " compute_commit_minus_encoded_eval[0]"));
    compute_commit_minus_encoded_eval.emplace_back(
        pb,
        encoded_evals[1],
        commitments[1],
        commit_minus_encoded_eval[1],
        FMT(annotation_prefix, " compute_commit_minus_encoded_eval[1]"));

    compute_gamma_times_commit_1_minus_encoded_eval_1.reset(
        new G1_mul_by_scalar_gadget<ppT>(
            pb,
            gamma,
            commit_minus_encoded_eval[1],
            gamma_times_commit_1_minus_encoded_eval_1,
            FMT(annotation_prefix,
                " compute_gamma_times_commit_1_minus_encoded_eval_1")));

    compute_result.reset(
        new G1_add_variable_and_variable_or_identity_gadget<ppT>(
            pb,
            gamma_times_commit_1_minus_encoded_eval_1,
            commit_minus_encoded_eval[0],
            result,
            FMT(annotation_prefix, " compute_result")));
}

template<typename ppT>
void kzg10_batched_compute_commit_minus_eval_sum<ppT, 2>::
    generate_r1cs_constraints()
{
    compute_encoded_evals[0].generate_r1cs_constraints();
    compute_encoded_evals[1].generate_r1cs_constraints();
    compute_commit_minus_encoded_eval[0].generate_r1cs_constraints();
    compute_commit_minus_encoded_eval[1].generate_r1cs_constraints();
    compute_gamma_times_commit_1_minus_encoded_eval_1
        ->generate_r1cs_constraints();
    compute_result->generate_r1cs_constraints();
}

template<typename ppT>
void kzg10_batched_compute_commit_minus_eval_sum<ppT, 2>::
    generate_r1cs_witness()
{
    compute_encoded_evals[0].generate_r1cs_witness();
    compute_encoded_evals[1].generate_r1cs_witness();
    compute_commit_minus_encoded_eval[0].generate_r1cs_witness();
    compute_commit_minus_encoded_eval[1].generate_r1cs_witness();
    compute_gamma_times_commit_1_minus_encoded_eval_1->generate_r1cs_witness();
    compute_result->generate_r1cs_witness();
}

// specialization for >2 entries
template<typename ppT, size_t num_entries>
kzg10_batched_compute_commit_minus_eval_sum<ppT, num_entries>::
    kzg10_batched_compute_commit_minus_eval_sum(
        protoboard<libff::Fr<ppT>> &pb,
        const pb_linear_combination<libff::Fr<ppT>> &gamma,
        const std::vector<kzg10_commitment_variable<ppT>> &commitments,
        const pb_linear_combination_array<libff::Fr<ppT>> &evals,
        G1_variable<ppT> &result,
        const std::string &annotation_prefix)
    : gadget<libff::Fr<ppT>>(pb, annotation_prefix)
    , encoded_evals(
          internal::allocate_variable_array<ppT, G1_variable_or_identity<ppT>>(
              pb, num_entries, FMT(annotation_prefix, " encoded_evals")))
    , compute_encoded_evals()
    , commit_minus_encoded_eval(
          internal::allocate_variable_array<ppT, G1_variable<ppT>>(
              pb,
              num_entries,
              FMT(annotation_prefix, " commit_minus_encoded_eval")))
    , compute_commit_minus_encoded_eval()
    , gamma(gamma)
    , gamma_powers()
    , gamma_power_times_commit_minus_encoded_eval(
          internal::allocate_variable_array<ppT, G1_variable_or_identity<ppT>>(
              pb,
              num_entries - 1,
              FMT(annotation_prefix,
                  " gamma_power_times_commit_minus_encoded_eval")))
    , compute_gamma_power_times_commit_minus_encoded_eval()
    , intermediate_sum(internal::allocate_variable_array<ppT, G1_variable<ppT>>(
          pb, num_entries - 2, FMT(annotation_prefix, " intermediate_sum")))
    , compute_intermediate_sum()
{
    // encoded_eval[i] = G1_mul_by_scalar(evals[i], G1::one())
    // len(encoded_eval) = num_entries
    G1_variable<ppT> g1_minus_one(
        pb,
        -libff::G1<other_curve<ppT>>::one(),
        FMT(annotation_prefix, " g1_one"));
    compute_encoded_evals.reserve(num_entries);
    for (size_t i = 0; i < num_entries; ++i) {
        compute_encoded_evals.emplace_back(
            pb,
            evals[i],
            g1_minus_one,
            encoded_evals[i],
            FMT(annotation_prefix, " compute_encoded_evals[%zu]", i));
    }

    // commit_minus_encoded_eval[i] = cm_i - [s_i]_1
    compute_commit_minus_encoded_eval.reserve(num_entries);
    for (size_t i = 0; i < num_entries; ++i) {
        compute_commit_minus_encoded_eval.emplace_back(
            pb,
            encoded_evals[i],
            commitments[i],
            commit_minus_encoded_eval[i],
            FMT(annotation_prefix, " compute_commit_minus_encoded_eval"));
    }

    // gamma_powers[0] = gamma * gamma
    // gamma_powers[i>0] = gamma * gamma_powers[i-1]
    gamma_powers.allocate(
        pb, num_entries - 2, FMT(annotation_prefix, " gamma_powers"));
    this->pb.add_r1cs_constraint(
        r1cs_constraint<Field>(gamma, gamma, gamma_powers[0]),
        FMT(annotation_prefix, " compute_gamma_power[0](gamma^2)"));
    for (size_t i = 1; i < num_entries - 2; ++i) {
        this->pb.add_r1cs_constraint(
            r1cs_constraint<Field>(gamma, gamma_powers[i - 1], gamma_powers[i]),
            FMT(annotation_prefix,
                " compute_gamma_power[%zu](gamma^%zu)",
                i,
                i + 2));
    }

    // gamma_power_times_commit_minus_encoded_eval[0] =
    //   G1_mul(gamma, compute_commit_minus_encoded_eval[1])
    // gamma_power_times_commit_minus_encoded_eval[i>0] =
    //   G1_mul(gamma_powers[i-1], commit_minus_encoded_eval[i+1]
    compute_gamma_power_times_commit_minus_encoded_eval.reserve(
        num_entries - 1);
    compute_gamma_power_times_commit_minus_encoded_eval.emplace_back(
        pb,
        gamma,
        commit_minus_encoded_eval[1],
        gamma_power_times_commit_minus_encoded_eval[0],
        FMT(annotation_prefix,
            "compute_gamma_power_times_commit_minus_encoded_eval[0]"));
    for (size_t i = 1; i < num_entries - 1; ++i) {
        compute_gamma_power_times_commit_minus_encoded_eval.emplace_back(
            pb,
            gamma_powers[i - 1],
            commit_minus_encoded_eval[i + 1],
            gamma_power_times_commit_minus_encoded_eval[i],
            FMT(annotation_prefix,
                "compute_gamma_power_times_commit_minus_encoded_eval[0]"));
    }

    // intermediate_sum[0] = G1_add(
    //   commit_minus_encoded_eval[0],
    //   gamma_power_times_commit_minus_encoded_eval[0])
    // intermediate_sum[0<i<num_entries - 2] = G1_add(
    //   intermediate_sum[i-1],
    //   gamma_power_times_commit_minus_encoded_eval[i])
    // result = G1_add(
    //   intermediate_sum[num_entries - 2],
    //   gamma_power_times_commit_minus_encoded_eval[num_entries - 2])
    compute_intermediate_sum.reserve(num_entries - 1);
    compute_intermediate_sum.emplace_back(
        pb,
        gamma_power_times_commit_minus_encoded_eval[0],
        commit_minus_encoded_eval[0],
        intermediate_sum[0],
        FMT(annotation_prefix, " compute_intermediate_sum[0]"));
    for (size_t i = 1; i < num_entries - 2; ++i) {
        compute_intermediate_sum.emplace_back(
            pb,
            gamma_power_times_commit_minus_encoded_eval[i],
            intermediate_sum[i - 1],
            intermediate_sum[i],
            FMT(annotation_prefix, " compute_intermediate_sum[%zu]", i));
    }

    compute_intermediate_sum.emplace_back(
        pb,
        gamma_power_times_commit_minus_encoded_eval[num_entries - 2],
        intermediate_sum[num_entries - 3],
        result,
        FMT(annotation_prefix,
            " compute_intermediate_sum[%zu]",
            num_entries - 2));
    assert(intermediate_sum.size() == num_entries - 2);
    assert(compute_intermediate_sum.size() == num_entries - 1);
}

template<typename ppT, size_t num_entries>
void kzg10_batched_compute_commit_minus_eval_sum<ppT, num_entries>::
    generate_r1cs_constraints()
{
    for (auto &gadget : compute_encoded_evals) {
        gadget.generate_r1cs_constraints();
    }

    for (auto &gadget : compute_commit_minus_encoded_eval) {
        gadget.generate_r1cs_constraints();
    }

    for (auto &gadget : compute_gamma_power_times_commit_minus_encoded_eval) {
        gadget.generate_r1cs_constraints();
    }

    for (auto &gadget : compute_intermediate_sum) {
        gadget.generate_r1cs_constraints();
    }
}

template<typename ppT, size_t num_entries>
void kzg10_batched_compute_commit_minus_eval_sum<ppT, num_entries>::
    generate_r1cs_witness()
{
    for (auto &gadget : compute_encoded_evals) {
        gadget.generate_r1cs_witness();
    }

    for (auto &gadget : compute_commit_minus_encoded_eval) {
        gadget.generate_r1cs_witness();
    }

    const Field gamma_val = this->pb.lc_val(gamma);
    Field gamma_power_val = gamma_val;

    for (auto &gamma_power : gamma_powers) {
        gamma_power_val = gamma_power_val * gamma_val;
        this->pb.val(gamma_power) = gamma_power_val;
    }

    for (auto &gadget : compute_gamma_power_times_commit_minus_encoded_eval) {
        gadget.generate_r1cs_witness();
    }

    for (auto &gadget : compute_intermediate_sum) {
        gadget.generate_r1cs_witness();
    }
}

} // namespace libsnark

#endif // LIBSNARK_GADGETLIB1_GADGETS_VERIFIERS_KZG10_BATCHED_VERIFIER_GADGET_TCC_
