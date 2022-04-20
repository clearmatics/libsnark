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
// kzg10_batched_witness_variable
//

template<typename ppT>
kzg10_batched_witness_variable<ppT>::kzg10_batched_witness_variable(
    protoboard<libff::Fr<ppT>> &pb, const std::string &annotationn_prefix)
    : W_1(pb, FMT(annotationn_prefix, " W_1"))
    , W_2(pb, FMT(annotationn_prefix, " W_2"))
{
}

template<typename ppT>
void kzg10_batched_witness_variable<ppT>::generate_r1cs_witness(
    const typename kzg10_batched_2_point<other_curve<ppT>>::evaluation_witness
        &eval_witness)
{
    W_1.generate_r1cs_witness(eval_witness.W_1);
    W_2.generate_r1cs_witness(eval_witness.W_2);
}

#if 0
//
// evaluate_polynomial_horner
//

template<typename ppT, size_t num_entries>
evaluate_polynomial_horner<ppT, num_entries>::evaluate_polynomial_horner(
    protoboard<libff::Fr<ppT>> &pb,
    const pb_linear_combination<libff::Fr<ppT>> &x,
    const pb_linear_combination_array<libff::Fr<ppT>> &coeffs,
    pb_variable<libff::Fr<ppT>> &result,
    const std::string &annotation_prefix)
    : gadget<libff::Fr<ppT>>(pb, annotation_prefix)
    , x(x)
    , coeffs(coeffs)
    , result(result)
{
    assert(coeffs.size() == num_entries);

    sum_terms.allocate(
        pb, num_entries - 2, FMT(annotation_prefix, " sum_terms"));
}

template<typename ppT, size_t num_entries>
void evaluate_polynomial_horner<ppT, num_entries>::generate_r1cs_constraints()
{
    // Conditions A * B = C are:
    //
    //   x * coeffs[n-1] = sum_terms[0] - coeffs[n-2]
    //   x * sum_terms[i-1] = sum_terms[i] - coeffs[n-(i+2)]
    //   x * sum_terms[n-3] = result - coeffs[0]
    //
    // for i = 0 ... n - 3

    pb.add_r1cs_constraint(
        {x, coeffs[num_entries - 1], sum_terms[0] - coeffs[num_entries - 2]},
        FMT(annotation_prefix, " compute_sum_term[0]"));

    for (size_t i = 1; i < num_entries - 2; ++i) {
        pb.add_r1cs_constraint(
            {x, sum_terms[i - 1], sum_terms[i] - coeffs[num_entries - (i + 2)]},
            FMT(annotation_prefix, " compute_sum_term[%zu]", i));
    }

    pb.add_r1cs_constraint(
        {x, sum_terms[num_entries - 3], result - coeffs[0]},
        FMT(annotation_prefix, " compute_sum_term[%zu]", num_entries - 2));
}

template<typename ppT, size_t num_entries>
void evaluate_polynomial_horner<ppT, num_entries>::generate_r1cs_witness()
{
    // sum_term[0] = x * coeffs[n-1] + coeffs[n-2]
    // sum_term[0 < i < n-2] = x * sum_term[i-1] + coeffs[n-(i+2)]
    // result = x * sum_term[n-3] + coeffs[0]

    const Field x = pb.val(x);
    pb.val(sum_term[0]) =
        x * pb.val(coeffs[num_entries - 1]) + pb.val(coeffs[num_entries - 2]);
    for (size_t i = 1; i < num_entries - 2; ++i) {
        pb.val(sum_terms[i]) = x * pb.val(sum_terms[i - 1]) +
                               pb.val(coeffs[num_entries - (i + 2)]);
    }
    pb.val(result) = x * pb.val(sum_terms[num_entries - 3]) + pb.val(coeffs[0]);
}
#endif

//
// kzg10_batched_compute_gamma_eval_sum
//

template<typename ppT, size_t n>
kzg10_batched_compute_gamma_eval_sum<ppT, n>::
    kzg10_batched_compute_gamma_eval_sum(
        protoboard<libff::Fr<ppT>> &pb,
        const pb_linear_combination<libff::Fr<ppT>> &gamma,
        const pb_variable_array<libff::Fr<ppT>> &gamma_powers,
        const pb_linear_combination_array<libff::Fr<ppT>> &evals,
        pb_variable<libff::Fr<ppT>> &result,
        const std::string &annotation_prefix)
    : gadget<libff::Fr<ppT>>(pb, annotation_prefix)
    , gamma(gamma)
    , evals(evals)
    , gamma_powers(gamma_powers)
    , sum_terms(pb_variable_array_allocate<libff::Fr<ppT>>(
          pb, n - 2, FMT(annotation_prefix, " sum_terms")))
    , result(result)
{
    assert(evals.size() == n);
    assert(gamma_powers.size() == n - 2);
    assert(sum_terms.size() == n - 2);
}

template<typename ppT, size_t n>
void kzg10_batched_compute_gamma_eval_sum<ppT, n>::generate_r1cs_constraints()
{
    protoboard<libff::Fr<ppT>> &pb = this->pb;
    const std::string &annotation_prefix = this->annotation_prefix;

    // Constraints A * B = C are:
    //
    //   gamma * eval[1] = sum_terms[0] - eval[0]
    //   for i = 0 ... n-4:
    //     gamma_powers[i] * eval[i + 2] = sum_terms[i+1] - sum_terms[i]
    //
    //   gamma_powers[n-3] * eval[n-1] = result - sum_terms[n-3]

    pb.add_r1cs_constraint(
        {gamma, evals[1], sum_terms[0] - evals[0]},
        FMT(annotation_prefix, " compute_sum_terms[0]"));

    for (size_t i = 0; i < n - 3; ++i) {
        pb.add_r1cs_constraint(
            {gamma_powers[i], evals[i + 2], sum_terms[i + 1] - sum_terms[i]},
            FMT(annotation_prefix, " compute_sum_terms[%zu]", i + 1));
    }

    assert(gamma_powers.size() == n - 2);
    assert(evals.size() == n);
    assert(sum_terms.size() == n - 2);

    pb.add_r1cs_constraint(
        {gamma_powers[n - 3], evals[n - 1], result - sum_terms[n - 3]},
        FMT(annotation_prefix, " compute_result"));
}

template<typename ppT, size_t n>
void kzg10_batched_compute_gamma_eval_sum<ppT, n>::generate_r1cs_witness()
{
    protoboard<libff::Fr<ppT>> &pb = this->pb;

    // sum_terms[0] = gamma * evals[1] + evals[0]
    pb.val(sum_terms[0]) =
        pb.lc_val(gamma) * pb.lc_val(evals[1]) + pb.lc_val(evals[0]);

    // sum_terms[i=1..n-3] = gamma_powers[i-1] * evals[i+1] + sum_terms[i-1]
    for (size_t i = 1; i < n - 2; ++i) {
        pb.val(sum_terms[i]) =
            pb.lc_val(gamma_powers[i - 1]) * pb.lc_val(evals[i + 1]) +
            pb.val(sum_terms[i - 1]);
    }

    // result = gamma_powers[n-3] * evals[n-1] + sum_terms[n-3]
    pb.val(result) = pb.lc_val(gamma_powers[n - 3]) * pb.lc_val(evals[n - 1]) +
                     pb.val(sum_terms[n - 3]);
}

//
// kzg10_batched_compute_commit_minus_eval_sum
//

// general case for >2 entries
template<typename ppT, size_t num_entries>
kzg10_batched_compute_commit_minus_eval_sum<ppT, num_entries>::
    kzg10_batched_compute_commit_minus_eval_sum(
        protoboard<libff::Fr<ppT>> &pb,
        const pb_linear_combination<libff::Fr<ppT>> &gamma,
        const std::vector<kzg10_commitment_variable<ppT>> &commits,
        const pb_linear_combination_array<libff::Fr<ppT>> &evals,
        G1_variable<ppT> &result,
        const std::string &annotation_prefix)
    : gadget<libff::Fr<ppT>>(pb, annotation_prefix)
    , evals(evals)
    , gamma(gamma)
    , gamma_powers(pb_variable_array_allocate<libff::Fr<ppT>>(
          pb, num_entries - 2, FMT(annotation_prefix, " gamma_powers")))
    , gamma_power_times_commit(
          internal::allocate_variable_array<ppT, G1_variable_or_identity<ppT>>(
              pb,
              num_entries - 1,
              FMT(annotation_prefix, " gamma_power_times_commit")))
    , compute_gamma_power_times_commit()
    , sum_gamma_power_times_commit(
          internal::allocate_variable_array<ppT, G1_variable<ppT>>(
              pb,
              num_entries - 1,
              FMT(annotation_prefix, " sum_gamma_power_times_commit")))
    , compute_sum_gamma_power_times_commit()
    , sum_gamma_power_times_eval(pb_variable_allocate<libff::Fr<ppT>>(
          pb, FMT(annotation_prefix, " sum_gamma_power_times_eval")))
    , compute_sum_gamma_power_times_eval(
          pb,
          gamma,
          gamma_powers,
          this->evals,
          sum_gamma_power_times_eval,
          FMT(annotation_prefix, " compute_sum_gamma_power_times_eval"))
    , encoded_sum_gamma_power_times_eval(
          pb, FMT(annotation_prefix, " encoded_sum_gamma_power_times_eval"))
    , compute_encoded_sum_gamma_power_times_eval(
          pb,
          sum_gamma_power_times_eval,
          G1_variable<ppT>(
              pb,
              -libff::G1<other_curve<ppT>>::one(),
              FMT(annotation_prefix, " minus_1_g1")),
          encoded_sum_gamma_power_times_eval,
          FMT(annotation_prefix, " compute_encoded_sum_gamma_power_times_eval"))
    , compute_result(
          pb,
          encoded_sum_gamma_power_times_eval,
          sum_gamma_power_times_commit[num_entries - 2],
          result,
          FMT(annotation_prefix, " compute_result"))
{
    // Array does not contain 1 or gamma, so:
    //
    //   gamma_powers[i] = gamma * gamma_powers[i - 1]
    //                   = gamma^{i + 2}
    //
    // for i = 0 ... num_entries - 3 (len = num_entries - 2)
    this->pb.add_r1cs_constraint(
        r1cs_constraint<Field>(this->gamma, this->gamma, gamma_powers[0]),
        FMT(annotation_prefix, " compute_gamma_power[0](gamma^2)"));
    for (size_t i = 1; i < num_entries - 2; ++i) {
        this->pb.add_r1cs_constraint(
            r1cs_constraint<Field>(
                this->gamma, gamma_powers[i - 1], gamma_powers[i]),
            FMT(annotation_prefix,
                " compute_gamma_power[%zu](gamma^%zu)",
                i,
                i + 2));
    }

    // TODO: wrap these two steps in a single
    // kzg10_batched_compute_gamma_powers_times_commit_sum gadget.

    // gamma_power_times_commit[0] = gamma * commits[1]
    // gamma_power_times_commit[i>0]
    //   = gamma^{i+1} * commits[i+1]
    //   = gamma_powers[i-1] * commits[i+1]
    //
    // for i = 0 ... num_entries - 2 (len = num_entries - 1)
    compute_gamma_power_times_commit.reserve(num_entries - 1);
    compute_gamma_power_times_commit.emplace_back(
        pb,
        gamma,
        commits[1],
        gamma_power_times_commit[0],
        FMT(annotation_prefix, "compute_gamma_power_times_commit[0]"));
    for (size_t i = 1; i < num_entries - 1; ++i) {
        compute_gamma_power_times_commit.emplace_back(
            pb,
            gamma_powers[i - 1],
            commits[i + 1],
            gamma_power_times_commit[i],
            FMT(annotation_prefix,
                "compute_gamma_power_times_commit_minus_encoded_eval[0]"));
    }

    // sum_gamma_power_times_commit[0] = commits[0] +
    // gamma_power_times_commit[0] sum_gamma_power_times_commit[i>0] =
    //   sum_gamma_power_times_commit[i-1] + gamma_power_times_commit[i]
    //
    // for i = 0 ... num_entries - 2 (len = num_entries - 1)
    compute_sum_gamma_power_times_commit.reserve(num_entries - 1);
    compute_sum_gamma_power_times_commit.emplace_back(
        pb,
        gamma_power_times_commit[0],
        commits[0],
        sum_gamma_power_times_commit[0],
        FMT(annotation_prefix, " compute_sum_gamma_power_times_commit[0]"));
    for (size_t i = 1; i < num_entries - 1; ++i) {
        compute_sum_gamma_power_times_commit.emplace_back(
            pb,
            gamma_power_times_commit[i],
            sum_gamma_power_times_commit[i - 1],
            sum_gamma_power_times_commit[i],
            FMT(annotation_prefix,
                " compute_sum_gamma_power_times_commit[%zu]",
                i));
    }
}

template<typename ppT, size_t num_entries>
void kzg10_batched_compute_commit_minus_eval_sum<ppT, num_entries>::
    generate_r1cs_constraints()
{
    for (auto &gadget : compute_gamma_power_times_commit) {
        gadget.generate_r1cs_constraints();
    }

    for (auto &gadget : compute_sum_gamma_power_times_commit) {
        gadget.generate_r1cs_constraints();
    }

    compute_sum_gamma_power_times_eval.generate_r1cs_constraints();
    compute_encoded_sum_gamma_power_times_eval.generate_r1cs_constraints();
    compute_result.generate_r1cs_constraints();
}

template<typename ppT, size_t num_entries>
void kzg10_batched_compute_commit_minus_eval_sum<ppT, num_entries>::
    generate_r1cs_witness()
{
    evals.evaluate(this->pb);

    const Field gamma_val = this->pb.lc_val(gamma);
    Field gamma_power_val = gamma_val;

    for (auto &gamma_power : gamma_powers) {
        gamma_power_val = gamma_power_val * gamma_val;
        this->pb.val(gamma_power) = gamma_power_val;
    }

    for (auto &gadget : compute_gamma_power_times_commit) {
        gadget.generate_r1cs_witness();
    }

    for (auto &gadget : compute_sum_gamma_power_times_commit) {
        gadget.generate_r1cs_witness();
    }

    compute_sum_gamma_power_times_eval.generate_r1cs_witness();
    compute_encoded_sum_gamma_power_times_eval.generate_r1cs_witness();
    compute_result.generate_r1cs_witness();
}

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

//
// kzg10_batched_verifier_gadget
//

template<typename ppT, size_t num_polyomials_1, size_t num_polyomials_2>
kzg10_batched_verifier_gadget<ppT, num_polyomials_1, num_polyomials_2>::
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
        const std::string &annotation_prefix)
    : gadget<libff::Fr<ppT>>(pb, annotation_prefix)
    , G(pb, FMT(annotation_prefix, " G"))
    , compute_G(
          pb,
          gamma_1,
          commitments_1,
          poly_evals_1,
          G,
          FMT(annotation_prefix, " compute_G"))
    , H(pb, FMT(annotation_prefix, " H"))
    , compute_H(
          pb,
          gamma_2,
          commitments_2,
          poly_evals_2,
          H,
          FMT(annotation_prefix, " compute_H"))
    , rH(pb, FMT(annotation_prefix, " rH"))
    // NOTE: this ignores H.is_identity
    , compute_rH(pb, r, H, rH, FMT(annotation_prefix, " compute_rH"))
    , F(pb, FMT(annotation_prefix, " F"))
    , compute_F(pb, rH, G, F, FMT(annotation_prefix, " compute_F"))
    //   A = W_1 + r * W_2
    , r_times_W_2(pb, FMT(annotation_prefix, " r_times_W_2"))
    , compute_r_times_W_2(
          pb,
          r,
          eval_witness.W_2,
          r_times_W_2,
          FMT(annotation_prefix, " compute_r_times_W_2"))
    , A(pb, FMT(annotation_prefix, " A"))
    , compute_A(
          pb,
          r_times_W_2,
          eval_witness.W_1,
          A,
          FMT(annotation_prefix, " compute_A"))
    //   B = F + z_1 * W_1 + r * z_2 * W_2
    , r_times_z_2_times_W_2(
          pb, FMT(annotation_prefix, " r_times_z_2_times_W_2"))
    , compute_r_times_z_2_times_W_2(
          pb,
          z_2,
          r_times_W_2,
          r_times_z_2_times_W_2,
          FMT(annotation_prefix, " compute_r_times_z_2_times_W_2"))
    , z_1_times_W_1(pb, FMT(annotation_prefix, " z_1_times_W_1"))
    , compute_z_1_times_W_1(
          pb,
          z_1,
          eval_witness.W_1,
          z_1_times_W_1,
          FMT(annotation_prefix, " compute_z_1_times_W_1"))
    , F_plus_z_1_times_W_1(pb, FMT(annotation_prefix, " F_plus_z_1_times_W_1"))
    , compute_F_plus_z_1_times_W_1(
          pb,
          z_1_times_W_1,
          F,
          F_plus_z_1_times_W_1,
          FMT(annotation_prefix, " compute_F_plus_z_1_times_W_1"))
    , B(pb, FMT(annotation_prefix, " B"))
    , compute_B(
          pb,
          r_times_z_2_times_W_2,
          F_plus_z_1_times_W_1,
          B,
          FMT(annotation_prefix, " compute_B"))
    , pairing_check(
          pb,
          A,
          srs.alpha_g2,
          B,
          libff::G2<other_curve<ppT>>::one(),
          result,
          FMT(annotation_prefix, " pairing_check"))
{
}

template<typename ppT, size_t num_polyomials_1, size_t num_polyomials_2>
void kzg10_batched_verifier_gadget<ppT, num_polyomials_1, num_polyomials_2>::
    generate_r1cs_constraints()
{
    compute_G.generate_r1cs_constraints();
    compute_H.generate_r1cs_constraints();
    compute_rH.generate_r1cs_constraints();
    compute_F.generate_r1cs_constraints();
    compute_r_times_W_2.generate_r1cs_constraints();
    compute_A.generate_r1cs_constraints();
    compute_r_times_z_2_times_W_2.generate_r1cs_constraints();
    compute_z_1_times_W_1.generate_r1cs_constraints();
    compute_F_plus_z_1_times_W_1.generate_r1cs_constraints();
    compute_B.generate_r1cs_constraints();
    pairing_check.generate_r1cs_constraints();
}

template<typename ppT, size_t num_polyomials_1, size_t num_polyomials_2>
void kzg10_batched_verifier_gadget<ppT, num_polyomials_1, num_polyomials_2>::
    generate_r1cs_witness()
{
    compute_G.generate_r1cs_witness();
    compute_H.generate_r1cs_witness();
    compute_rH.generate_r1cs_witness();
    compute_F.generate_r1cs_witness();
    compute_r_times_W_2.generate_r1cs_witness();
    compute_A.generate_r1cs_witness();
    compute_r_times_z_2_times_W_2.generate_r1cs_witness();
    compute_z_1_times_W_1.generate_r1cs_witness();
    compute_F_plus_z_1_times_W_1.generate_r1cs_witness();
    compute_B.generate_r1cs_witness();
    pairing_check.generate_r1cs_witness();

    // // result = (1 - B.is_identity) * pairing_check_result
    // const Field pairing_check_result_val =
    // this->pb.val(pairing_check_result); const Field B_is_identity_val =
    // this->pb.lc_val(B.is_identity); this->pb.val(result) =
    //     (Field::one() - B_is_identity_val) * pairing_check_result_val;
}

} // namespace libsnark

#endif // LIBSNARK_GADGETLIB1_GADGETS_VERIFIERS_KZG10_BATCHED_VERIFIER_GADGET_TCC_
