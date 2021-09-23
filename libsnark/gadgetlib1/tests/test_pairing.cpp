/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include "libsnark/gadgetlib1/gadgets/pairing/mnt/mnt_pairing_params.hpp"
#include "libsnark/gadgetlib1/protoboard.hpp"

#include <gtest/gtest.h>
#include <libff/algebra/curves/mnt/mnt4/mnt4_pp.hpp>
#include <libff/algebra/curves/mnt/mnt6/mnt6_pp.hpp>

namespace libsnark
{

template<typename ppT>
void test_G1_variable_precomp(const std::string &annotation)
{
    protoboard<libff::Fr<ppT>> pb;
    libff::G1<other_curve<ppT>> g_val =
        libff::Fr<other_curve<ppT>>::random_element() *
        libff::G1<other_curve<ppT>>::one();

    G1_variable<ppT> g(pb, "g");
    mnt_G1_precomputation<ppT> precomp;
    mnt_precompute_G1_gadget<ppT> do_precomp(pb, g, precomp, "do_precomp");
    do_precomp.generate_r1cs_constraints();

    g.generate_r1cs_witness(g_val);
    do_precomp.generate_r1cs_witness();
    ASSERT_TRUE(pb.is_satisfied());

    mnt_G1_precomputation<ppT> const_precomp(pb, g_val, "const_precomp");

    libff::affine_ate_G1_precomp<other_curve<ppT>> native_precomp =
        other_curve<ppT>::affine_ate_precompute_G1(g_val);
    ASSERT_EQ(
        precomp.PY_twist_squared->get_element(),
        native_precomp.PY_twist_squared);
    ASSERT_EQ(
        const_precomp.PY_twist_squared->get_element(),
        native_precomp.PY_twist_squared);

    printf(
        "number of constraints for G1 precomp (Fr is %s)  = %zu\n",
        annotation.c_str(),
        pb.num_constraints());
}

template<typename ppT>
void test_G2_variable_precomp(const std::string &annotation)
{
    protoboard<libff::Fr<ppT>> pb;
    libff::G2<other_curve<ppT>> g_val =
        libff::Fr<other_curve<ppT>>::random_element() *
        libff::G2<other_curve<ppT>>::one();

    G2_variable<ppT> g(pb, "g");
    mnt_G2_precomputation<ppT> precomp;
    mnt_precompute_G2_gadget<ppT> do_precomp(pb, g, precomp, "do_precomp");
    do_precomp.generate_r1cs_constraints();

    g.generate_r1cs_witness(g_val);
    do_precomp.generate_r1cs_witness();
    ASSERT_TRUE(pb.is_satisfied());

    libff::affine_ate_G2_precomp<other_curve<ppT>> native_precomp =
        other_curve<ppT>::affine_ate_precompute_G2(g_val);

    // the last precomp is unused, but remains for convenient
    // programming
    ASSERT_EQ(precomp.coeffs.size() - 1, native_precomp.coeffs.size());
    for (size_t i = 0; i < native_precomp.coeffs.size(); ++i) {
        ASSERT_EQ(
            precomp.coeffs[i]->RX->get_element(),
            native_precomp.coeffs[i].old_RX);
        ASSERT_EQ(
            precomp.coeffs[i]->RY->get_element(),
            native_precomp.coeffs[i].old_RY);
        ASSERT_EQ(
            precomp.coeffs[i]->gamma->get_element(),
            native_precomp.coeffs[i].gamma);
        ASSERT_EQ(
            precomp.coeffs[i]->gamma_X->get_element(),
            native_precomp.coeffs[i].gamma_X);
    }

    printf(
        "number of constraints for G2 precomp (Fr is %s)  = %zu\n",
        annotation.c_str(),
        pb.num_constraints());
}

template<typename ppT> void test_mnt_miller_loop(const std::string &annotation)
{
    protoboard<libff::Fr<ppT>> pb;
    libff::G1<other_curve<ppT>> P_val =
        libff::Fr<other_curve<ppT>>::random_element() *
        libff::G1<other_curve<ppT>>::one();
    libff::G2<other_curve<ppT>> Q_val =
        libff::Fr<other_curve<ppT>>::random_element() *
        libff::G2<other_curve<ppT>>::one();

    G1_variable<ppT> P(pb, "P");
    G2_variable<ppT> Q(pb, "Q");

    mnt_G1_precomputation<ppT> prec_P;
    mnt_G2_precomputation<ppT> prec_Q;

    mnt_precompute_G1_gadget<ppT> compute_prec_P(pb, P, prec_P, "prec_P");
    mnt_precompute_G2_gadget<ppT> compute_prec_Q(pb, Q, prec_Q, "prec_Q");

    Fqk_variable<ppT> result(pb, "result");
    mnt_miller_loop_gadget<ppT> miller(pb, prec_P, prec_Q, result, "miller");

    PROFILE_CONSTRAINTS(pb, "precompute P")
    {
        compute_prec_P.generate_r1cs_constraints();
    }
    PROFILE_CONSTRAINTS(pb, "precompute Q")
    {
        compute_prec_Q.generate_r1cs_constraints();
    }
    PROFILE_CONSTRAINTS(pb, "Miller loop")
    {
        miller.generate_r1cs_constraints();
    }
    PRINT_CONSTRAINT_PROFILING();

    P.generate_r1cs_witness(P_val);
    compute_prec_P.generate_r1cs_witness();
    Q.generate_r1cs_witness(Q_val);
    compute_prec_Q.generate_r1cs_witness();
    miller.generate_r1cs_witness();
    assert(pb.is_satisfied());

    libff::affine_ate_G1_precomp<other_curve<ppT>> native_prec_P =
        other_curve<ppT>::affine_ate_precompute_G1(P_val);
    libff::affine_ate_G2_precomp<other_curve<ppT>> native_prec_Q =
        other_curve<ppT>::affine_ate_precompute_G2(Q_val);
    libff::Fqk<other_curve<ppT>> native_result =
        other_curve<ppT>::affine_ate_miller_loop(native_prec_P, native_prec_Q);

    assert(result.get_element() == native_result);
    printf(
        "number of constraints for Miller loop (Fr is %s)  = %zu\n",
        annotation.c_str(),
        pb.num_constraints());
}

template<typename ppT>
void test_mnt_e_over_e_miller_loop(const std::string &annotation)
{
    protoboard<libff::Fr<ppT>> pb;
    libff::G1<other_curve<ppT>> P1_val =
        libff::Fr<other_curve<ppT>>::random_element() *
        libff::G1<other_curve<ppT>>::one();
    libff::G2<other_curve<ppT>> Q1_val =
        libff::Fr<other_curve<ppT>>::random_element() *
        libff::G2<other_curve<ppT>>::one();

    libff::G1<other_curve<ppT>> P2_val =
        libff::Fr<other_curve<ppT>>::random_element() *
        libff::G1<other_curve<ppT>>::one();
    libff::G2<other_curve<ppT>> Q2_val =
        libff::Fr<other_curve<ppT>>::random_element() *
        libff::G2<other_curve<ppT>>::one();

    G1_variable<ppT> P1(pb, "P1");
    G2_variable<ppT> Q1(pb, "Q1");
    G1_variable<ppT> P2(pb, "P2");
    G2_variable<ppT> Q2(pb, "Q2");

    mnt_G1_precomputation<ppT> prec_P1;
    mnt_precompute_G1_gadget<ppT> compute_prec_P1(
        pb, P1, prec_P1, "compute_prec_P1");
    mnt_G1_precomputation<ppT> prec_P2;
    mnt_precompute_G1_gadget<ppT> compute_prec_P2(
        pb, P2, prec_P2, "compute_prec_P2");
    mnt_G2_precomputation<ppT> prec_Q1;
    mnt_precompute_G2_gadget<ppT> compute_prec_Q1(
        pb, Q1, prec_Q1, "compute_prec_Q1");
    mnt_G2_precomputation<ppT> prec_Q2;
    mnt_precompute_G2_gadget<ppT> compute_prec_Q2(
        pb, Q2, prec_Q2, "compute_prec_Q2");

    Fqk_variable<ppT> result(pb, "result");
    mnt_e_over_e_miller_loop_gadget<ppT> miller(
        pb, prec_P1, prec_Q1, prec_P2, prec_Q2, result, "miller");

    PROFILE_CONSTRAINTS(pb, "precompute P")
    {
        compute_prec_P1.generate_r1cs_constraints();
        compute_prec_P2.generate_r1cs_constraints();
    }
    PROFILE_CONSTRAINTS(pb, "precompute Q")
    {
        compute_prec_Q1.generate_r1cs_constraints();
        compute_prec_Q2.generate_r1cs_constraints();
    }
    PROFILE_CONSTRAINTS(pb, "Miller loop")
    {
        miller.generate_r1cs_constraints();
    }
    PRINT_CONSTRAINT_PROFILING();

    P1.generate_r1cs_witness(P1_val);
    compute_prec_P1.generate_r1cs_witness();
    Q1.generate_r1cs_witness(Q1_val);
    compute_prec_Q1.generate_r1cs_witness();
    P2.generate_r1cs_witness(P2_val);
    compute_prec_P2.generate_r1cs_witness();
    Q2.generate_r1cs_witness(Q2_val);
    compute_prec_Q2.generate_r1cs_witness();
    miller.generate_r1cs_witness();
    assert(pb.is_satisfied());

    libff::affine_ate_G1_precomp<other_curve<ppT>> native_prec_P1 =
        other_curve<ppT>::affine_ate_precompute_G1(P1_val);
    libff::affine_ate_G2_precomp<other_curve<ppT>> native_prec_Q1 =
        other_curve<ppT>::affine_ate_precompute_G2(Q1_val);
    libff::affine_ate_G1_precomp<other_curve<ppT>> native_prec_P2 =
        other_curve<ppT>::affine_ate_precompute_G1(P2_val);
    libff::affine_ate_G2_precomp<other_curve<ppT>> native_prec_Q2 =
        other_curve<ppT>::affine_ate_precompute_G2(Q2_val);
    libff::Fqk<other_curve<ppT>> native_result =
        (other_curve<ppT>::affine_ate_miller_loop(
             native_prec_P1, native_prec_Q1) *
         other_curve<ppT>::affine_ate_miller_loop(
             native_prec_P2, native_prec_Q2)
             .inverse());

    assert(result.get_element() == native_result);
    printf(
        "number of constraints for e over e Miller loop (Fr is %s)  = %zu\n",
        annotation.c_str(),
        pb.num_constraints());
}

template<typename ppT>
void test_mnt_e_times_e_over_e_miller_loop(const std::string &annotation)
{
    protoboard<libff::Fr<ppT>> pb;
    libff::G1<other_curve<ppT>> P1_val =
        libff::Fr<other_curve<ppT>>::random_element() *
        libff::G1<other_curve<ppT>>::one();
    libff::G2<other_curve<ppT>> Q1_val =
        libff::Fr<other_curve<ppT>>::random_element() *
        libff::G2<other_curve<ppT>>::one();

    libff::G1<other_curve<ppT>> P2_val =
        libff::Fr<other_curve<ppT>>::random_element() *
        libff::G1<other_curve<ppT>>::one();
    libff::G2<other_curve<ppT>> Q2_val =
        libff::Fr<other_curve<ppT>>::random_element() *
        libff::G2<other_curve<ppT>>::one();

    libff::G1<other_curve<ppT>> P3_val =
        libff::Fr<other_curve<ppT>>::random_element() *
        libff::G1<other_curve<ppT>>::one();
    libff::G2<other_curve<ppT>> Q3_val =
        libff::Fr<other_curve<ppT>>::random_element() *
        libff::G2<other_curve<ppT>>::one();

    G1_variable<ppT> P1(pb, "P1");
    G2_variable<ppT> Q1(pb, "Q1");
    G1_variable<ppT> P2(pb, "P2");
    G2_variable<ppT> Q2(pb, "Q2");
    G1_variable<ppT> P3(pb, "P3");
    G2_variable<ppT> Q3(pb, "Q3");

    mnt_G1_precomputation<ppT> prec_P1;
    mnt_precompute_G1_gadget<ppT> compute_prec_P1(
        pb, P1, prec_P1, "compute_prec_P1");
    mnt_G1_precomputation<ppT> prec_P2;
    mnt_precompute_G1_gadget<ppT> compute_prec_P2(
        pb, P2, prec_P2, "compute_prec_P2");
    mnt_G1_precomputation<ppT> prec_P3;
    mnt_precompute_G1_gadget<ppT> compute_prec_P3(
        pb, P3, prec_P3, "compute_prec_P3");
    mnt_G2_precomputation<ppT> prec_Q1;
    mnt_precompute_G2_gadget<ppT> compute_prec_Q1(
        pb, Q1, prec_Q1, "compute_prec_Q1");
    mnt_G2_precomputation<ppT> prec_Q2;
    mnt_precompute_G2_gadget<ppT> compute_prec_Q2(
        pb, Q2, prec_Q2, "compute_prec_Q2");
    mnt_G2_precomputation<ppT> prec_Q3;
    mnt_precompute_G2_gadget<ppT> compute_prec_Q3(
        pb, Q3, prec_Q3, "compute_prec_Q3");

    Fqk_variable<ppT> result(pb, "result");
    mnt_e_times_e_over_e_miller_loop_gadget<ppT> miller(
        pb,
        prec_P1,
        prec_Q1,
        prec_P2,
        prec_Q2,
        prec_P3,
        prec_Q3,
        result,
        "miller");

    PROFILE_CONSTRAINTS(pb, "precompute P")
    {
        compute_prec_P1.generate_r1cs_constraints();
        compute_prec_P2.generate_r1cs_constraints();
        compute_prec_P3.generate_r1cs_constraints();
    }
    PROFILE_CONSTRAINTS(pb, "precompute Q")
    {
        compute_prec_Q1.generate_r1cs_constraints();
        compute_prec_Q2.generate_r1cs_constraints();
        compute_prec_Q3.generate_r1cs_constraints();
    }
    PROFILE_CONSTRAINTS(pb, "Miller loop")
    {
        miller.generate_r1cs_constraints();
    }
    PRINT_CONSTRAINT_PROFILING();

    P1.generate_r1cs_witness(P1_val);
    compute_prec_P1.generate_r1cs_witness();
    Q1.generate_r1cs_witness(Q1_val);
    compute_prec_Q1.generate_r1cs_witness();
    P2.generate_r1cs_witness(P2_val);
    compute_prec_P2.generate_r1cs_witness();
    Q2.generate_r1cs_witness(Q2_val);
    compute_prec_Q2.generate_r1cs_witness();
    P3.generate_r1cs_witness(P3_val);
    compute_prec_P3.generate_r1cs_witness();
    Q3.generate_r1cs_witness(Q3_val);
    compute_prec_Q3.generate_r1cs_witness();
    miller.generate_r1cs_witness();
    assert(pb.is_satisfied());

    libff::affine_ate_G1_precomp<other_curve<ppT>> native_prec_P1 =
        other_curve<ppT>::affine_ate_precompute_G1(P1_val);
    libff::affine_ate_G2_precomp<other_curve<ppT>> native_prec_Q1 =
        other_curve<ppT>::affine_ate_precompute_G2(Q1_val);
    libff::affine_ate_G1_precomp<other_curve<ppT>> native_prec_P2 =
        other_curve<ppT>::affine_ate_precompute_G1(P2_val);
    libff::affine_ate_G2_precomp<other_curve<ppT>> native_prec_Q2 =
        other_curve<ppT>::affine_ate_precompute_G2(Q2_val);
    libff::affine_ate_G1_precomp<other_curve<ppT>> native_prec_P3 =
        other_curve<ppT>::affine_ate_precompute_G1(P3_val);
    libff::affine_ate_G2_precomp<other_curve<ppT>> native_prec_Q3 =
        other_curve<ppT>::affine_ate_precompute_G2(Q3_val);
    libff::Fqk<other_curve<ppT>> native_result =
        (other_curve<ppT>::affine_ate_miller_loop(
             native_prec_P1, native_prec_Q1) *
         other_curve<ppT>::affine_ate_miller_loop(
             native_prec_P2, native_prec_Q2) *
         other_curve<ppT>::affine_ate_miller_loop(
             native_prec_P3, native_prec_Q3)
             .inverse());

    assert(result.get_element() == native_result);
    printf(
        "number of constraints for e times e over e Miller loop (Fr is %s)  = "
        "%zu\n",
        annotation.c_str(),
        pb.num_constraints());
}

template<typename ppT>
bool test_mnt_e_times_e_times_e_over_e_miller_loop(
    const std::string &annotation)
{
    // This `using` directive is used to avoid
    // any trouble with the use the of libsnark macros
    // in this function.
    using namespace libsnark;

    libsnark::protoboard<libff::Fr<ppT>> pb;
    libff::G1<other_curve<ppT>> P1_val =
        libff::Fr<other_curve<ppT>>::random_element() *
        libff::G1<other_curve<ppT>>::one();
    libff::G2<other_curve<ppT>> Q1_val =
        libff::Fr<other_curve<ppT>>::random_element() *
        libff::G2<other_curve<ppT>>::one();

    libff::G1<other_curve<ppT>> P2_val =
        libff::Fr<other_curve<ppT>>::random_element() *
        libff::G1<other_curve<ppT>>::one();
    libff::G2<other_curve<ppT>> Q2_val =
        libff::Fr<other_curve<ppT>>::random_element() *
        libff::G2<other_curve<ppT>>::one();

    libff::G1<other_curve<ppT>> P3_val =
        libff::Fr<other_curve<ppT>>::random_element() *
        libff::G1<other_curve<ppT>>::one();
    libff::G2<other_curve<ppT>> Q3_val =
        libff::Fr<other_curve<ppT>>::random_element() *
        libff::G2<other_curve<ppT>>::one();

    libff::G1<other_curve<ppT>> P4_val =
        libff::Fr<other_curve<ppT>>::random_element() *
        libff::G1<other_curve<ppT>>::one();
    libff::G2<other_curve<ppT>> Q4_val =
        libff::Fr<other_curve<ppT>>::random_element() *
        libff::G2<other_curve<ppT>>::one();

    libsnark::G1_variable<ppT> P1(pb, "P1");
    libsnark::G2_variable<ppT> Q1(pb, "Q1");
    libsnark::G1_variable<ppT> P2(pb, "P2");
    libsnark::G2_variable<ppT> Q2(pb, "Q2");
    libsnark::G1_variable<ppT> P3(pb, "P3");
    libsnark::G2_variable<ppT> Q3(pb, "Q3");
    libsnark::G1_variable<ppT> P4(pb, "P4");
    libsnark::G2_variable<ppT> Q4(pb, "Q4");

    libsnark::G1_precomputation<ppT> prec_P1;
    libsnark::precompute_G1_gadget<ppT> compute_prec_P1(
        pb, P1, prec_P1, "compute_prec_P1");
    libsnark::G1_precomputation<ppT> prec_P2;
    libsnark::precompute_G1_gadget<ppT> compute_prec_P2(
        pb, P2, prec_P2, "compute_prec_P2");
    libsnark::G1_precomputation<ppT> prec_P3;
    libsnark::precompute_G1_gadget<ppT> compute_prec_P3(
        pb, P3, prec_P3, "compute_prec_P3");
    libsnark::G1_precomputation<ppT> prec_P4;
    libsnark::precompute_G1_gadget<ppT> compute_prec_P4(
        pb, P4, prec_P4, "compute_prec_P4");
    libsnark::G2_precomputation<ppT> prec_Q1;
    libsnark::precompute_G2_gadget<ppT> compute_prec_Q1(
        pb, Q1, prec_Q1, "compute_prec_Q1");
    libsnark::G2_precomputation<ppT> prec_Q2;
    libsnark::precompute_G2_gadget<ppT> compute_prec_Q2(
        pb, Q2, prec_Q2, "compute_prec_Q2");
    libsnark::G2_precomputation<ppT> prec_Q3;
    libsnark::precompute_G2_gadget<ppT> compute_prec_Q3(
        pb, Q3, prec_Q3, "compute_prec_Q3");
    libsnark::G2_precomputation<ppT> prec_Q4;
    libsnark::precompute_G2_gadget<ppT> compute_prec_Q4(
        pb, Q4, prec_Q4, "compute_prec_Q4");

    Fqk_variable<ppT> result(pb, "result");

    mnt_e_times_e_times_e_over_e_miller_loop_gadget<ppT> miller(
        pb,
        prec_P1,
        prec_Q1,
        prec_P2,
        prec_Q2,
        prec_P3,
        prec_Q3,
        prec_P4,
        prec_Q4,
        result,
        "miller");

    PROFILE_CONSTRAINTS(pb, "precompute P")
    {
        compute_prec_P1.generate_r1cs_constraints();
        compute_prec_P2.generate_r1cs_constraints();
        compute_prec_P3.generate_r1cs_constraints();
        compute_prec_P4.generate_r1cs_constraints();
    }
    PROFILE_CONSTRAINTS(pb, "precompute Q")
    {
        compute_prec_Q1.generate_r1cs_constraints();
        compute_prec_Q2.generate_r1cs_constraints();
        compute_prec_Q3.generate_r1cs_constraints();
        compute_prec_Q4.generate_r1cs_constraints();
    }
    PROFILE_CONSTRAINTS(pb, "Miller loop")
    {
        miller.generate_r1cs_constraints();
    }
    PRINT_CONSTRAINT_PROFILING();

    P1.generate_r1cs_witness(P1_val);
    compute_prec_P1.generate_r1cs_witness();
    Q1.generate_r1cs_witness(Q1_val);
    compute_prec_Q1.generate_r1cs_witness();
    P2.generate_r1cs_witness(P2_val);
    compute_prec_P2.generate_r1cs_witness();
    Q2.generate_r1cs_witness(Q2_val);
    compute_prec_Q2.generate_r1cs_witness();
    P3.generate_r1cs_witness(P3_val);
    compute_prec_P3.generate_r1cs_witness();
    Q3.generate_r1cs_witness(Q3_val);
    compute_prec_Q3.generate_r1cs_witness();
    P4.generate_r1cs_witness(P4_val);
    compute_prec_P4.generate_r1cs_witness();
    Q4.generate_r1cs_witness(Q4_val);
    compute_prec_Q4.generate_r1cs_witness();
    miller.generate_r1cs_witness();

    assert(pb.is_satisfied());

    libff::affine_ate_G1_precomp<other_curve<ppT>> native_prec_P1 =
        other_curve<ppT>::affine_ate_precompute_G1(P1_val);
    libff::affine_ate_G2_precomp<other_curve<ppT>> native_prec_Q1 =
        other_curve<ppT>::affine_ate_precompute_G2(Q1_val);
    libff::affine_ate_G1_precomp<other_curve<ppT>> native_prec_P2 =
        other_curve<ppT>::affine_ate_precompute_G1(P2_val);
    libff::affine_ate_G2_precomp<other_curve<ppT>> native_prec_Q2 =
        other_curve<ppT>::affine_ate_precompute_G2(Q2_val);
    libff::affine_ate_G1_precomp<other_curve<ppT>> native_prec_P3 =
        other_curve<ppT>::affine_ate_precompute_G1(P3_val);
    libff::affine_ate_G2_precomp<other_curve<ppT>> native_prec_Q3 =
        other_curve<ppT>::affine_ate_precompute_G2(Q3_val);
    libff::affine_ate_G1_precomp<other_curve<ppT>> native_prec_P4 =
        other_curve<ppT>::affine_ate_precompute_G1(P4_val);
    libff::affine_ate_G2_precomp<other_curve<ppT>> native_prec_Q4 =
        other_curve<ppT>::affine_ate_precompute_G2(Q4_val);
    libff::Fqk<other_curve<ppT>> native_result =
        (other_curve<ppT>::affine_ate_miller_loop(
             native_prec_P1, native_prec_Q1) *
         other_curve<ppT>::affine_ate_miller_loop(
             native_prec_P2, native_prec_Q2) *
         other_curve<ppT>::affine_ate_miller_loop(
             native_prec_P3, native_prec_Q3) *
         other_curve<ppT>::affine_ate_miller_loop(
             native_prec_P4, native_prec_Q4)
             .inverse());

    printf(
        "number of constraints for e times e times e over e Miller loop (Fr is "
        "%s)  = %zu\n",
        annotation.c_str(),
        pb.num_constraints());

    return result.get_element() == native_result;
}

TEST(Pairing, MNT4_G1_Precompute)
{
    test_G1_variable_precomp<libff::mnt4_pp>("mnt4");
}

TEST(Pairing, MNT6_G1_Precompute)
{
    test_G1_variable_precomp<libff::mnt6_pp>("mnt6");
}

TEST(Pairing, MNT4_G2_Precompute)
{
    test_G2_variable_precomp<libff::mnt4_pp>("mnt4");
}

TEST(Pairing, MNT6_G2_Precompute)
{
    test_G2_variable_precomp<libff::mnt6_pp>("mnt6");
}

TEST(Pairing, MNT4_Miller_Loop)
{
    test_mnt_miller_loop<libff::mnt4_pp>("mnt4");
}

TEST(Pairing, MNT6_Miller_Loop)
{
    test_mnt_miller_loop<libff::mnt6_pp>("mnt6");
}

TEST(Pairing, MNT4_e_over_e_Miller_Loop)
{
    test_mnt_e_over_e_miller_loop<libff::mnt4_pp>("mnt4");
}

TEST(Pairing, MNT6_e_over_e_Miller_Loop)
{
    test_mnt_e_over_e_miller_loop<libff::mnt6_pp>("mnt6");
}

TEST(Pairing, MNT4_e_times_e_times_e_over_e_Miller_Loop)
{
    test_mnt_e_times_e_times_e_over_e_miller_loop<libff::mnt4_pp>("mnt4");
}

TEST(Pairing, MNT6_e_times_e_times_e_over_e_Miller_Loop)
{
    test_mnt_e_times_e_times_e_over_e_miller_loop<libff::mnt6_pp>("mnt6");
}

} // namespace libsnark

int main(int argc, char **argv)
{
    libff::mnt4_pp::init_public_params();
    libff::mnt6_pp::init_public_params();
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
