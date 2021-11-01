/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <algorithm>
#include <gtest/gtest.h>
#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>
#include <libff/algebra/curves/bls12_377/bls12_377_pp.hpp>
#include <libsnark/polynomial_commitments/kzg10.hpp>
#include <libsnark/polynomial_commitments/tests/polynomial_commitment_test_utils.hpp>

const size_t MAX_DEGREE = 254;

namespace libsnark
{

template<typename ppT, typename scheme> void test_basic_commitment()
{
    const polynomial<libff::Fr<ppT>> phi = gen_polynomial<ppT>(MAX_DEGREE);
    const polynomial<libff::Fr<ppT>> phi_2 = gen_polynomial<ppT>(MAX_DEGREE);
    const typename scheme::srs srs = scheme::setup(MAX_DEGREE);
    const typename scheme::commitment C = scheme::commit(srs, phi);

    ASSERT_TRUE(scheme::verify_poly(srs, C, phi));
    ASSERT_FALSE(scheme::verify_poly(srs, C, phi_2));
}

template<typename ppT, typename scheme> void test_eval_commitment()
{
    using Field = libff::Fr<ppT>;

    const polynomial<Field> phi = gen_polynomial<ppT>(MAX_DEGREE);
    const polynomial<Field> phi_2 = gen_polynomial<ppT>(MAX_DEGREE);
    const typename scheme::srs srs = scheme::setup(MAX_DEGREE);
    const typename scheme::commitment C = scheme::commit(srs, phi);
    const Field i = Field::random_element();
    const evaluation_and_witness<ppT, scheme> eval_witness =
        scheme::create_witness(srs, phi, i);

    evaluation_and_witness<ppT, scheme> eval_witness_wrong_i =
        scheme::create_witness(srs, phi, i + Field::one());
    eval_witness_wrong_i.i = i;

    const evaluation_and_witness<ppT, scheme> eval_witness_2 =
        scheme::create_witness(srs, phi_2, i);

    ASSERT_TRUE(scheme::verify_eval(srs, C, eval_witness));
    ASSERT_FALSE(scheme::verify_eval(srs, C, eval_witness_wrong_i));
    ASSERT_FALSE(scheme::verify_eval(srs, C, eval_witness_2));
}

template<typename ppT> void test_kzg10_commitment_with_known_secret()
{
    using Field = libff::Fr<ppT>;
    using scheme = kzg10<ppT>;

    // Dummy polynomial
    //   phi(x) = -1 + x + 2x^2 + 3x^3
    //
    //   phi(x) - phi(i) = (x - i) + 2(x^2 - i^2) + 3(x^3 - i^3)
    //   phi(x) - phi(i)/(x - i) =
    //     (x - i)/(x - i) + 2(x^2 - i^2)/(x - i) + 3(x^3 - i^3)/(x - i)
    //   = 1 + 2(x + i) + 3(x^2 + ix + i^2)
    //   = (1 + 2i + 3i^2) + (2 + 3i)x + 3x^2

    // Dummy secret and evaluation point
    //   alpha = 10
    //   i = 2
    //
    //   phi(alpha) = -1 + 10 + 200 + 3000 = 3209
    //   phi(i) = -1 + 2 + 8 + 24 = 33
    //
    //   psi(x) = (phi(x) - phi(i)/(x - i)
    //          = 17 + 8x + 3x^2
    //   psi(alpha) = 17 + 80 + 300 = 397

    const Field alpha = Field(10);
    const Field i = Field(2);
    const polynomial<Field> phi = {-Field("1"), Field(1), Field(2), Field(3)};

    // Check the srs
    const typename scheme::srs srs = scheme::setup_from_secret(16, alpha);
    ASSERT_EQ(libff::G1<ppT>::one(), srs.alpha_powers_g1[0]);
    ASSERT_EQ(Field(10) * libff::G1<ppT>::one(), srs.alpha_powers_g1[1]);
    ASSERT_EQ(Field(100) * libff::G1<ppT>::one(), srs.alpha_powers_g1[2]);
    ASSERT_EQ(Field(1000) * libff::G1<ppT>::one(), srs.alpha_powers_g1[3]);
    ASSERT_EQ(alpha * libff::G2<ppT>::one(), srs.alpha_g2);

    // Check the commitment
    const typename scheme::commitment C = scheme::commit(srs, phi);
    ASSERT_EQ(Field(3209) * libff::G1<ppT>::one(), C);

    // Check the evaluation with witness
    const evaluation_and_witness<ppT, scheme> eval_witness =
        scheme::create_witness(srs, phi, i);

    ASSERT_EQ(i, eval_witness.i);
    ASSERT_EQ(Field(33), eval_witness.phi_i);
    ASSERT_EQ(Field(397) * libff::G1<ppT>::one(), eval_witness.w);
}

template<typename ppT> void test_for_curve()
{
    // Execute all tests for the given curve.

    ppT::init_public_params();

    // KZG10
    test_kzg10_commitment_with_known_secret<ppT>();
    test_basic_commitment<ppT, kzg10<ppT>>();
    test_eval_commitment<ppT, kzg10<ppT>>();
}

TEST(TestPolynomialCommitments, ALT_BN128)
{
    test_for_curve<libff::alt_bn128_pp>();
}

TEST(TestPolynomialCommitments, BLS12_377)
{
    test_for_curve<libff::bls12_377_pp>();
}

} // namespace libsnark
