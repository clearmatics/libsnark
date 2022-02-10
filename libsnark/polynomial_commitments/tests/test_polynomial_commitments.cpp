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
#include <libsnark/polynomial_commitments/kzg10_batched.hpp>
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
    const Field eval = scheme::evaluate_polynomial(phi, i);
    const typename scheme::evaluation_witness eval_witness =
        scheme::create_evaluation_witness(phi, i, eval, srs);

    ASSERT_TRUE(scheme::verify_evaluation(i, eval, srs, eval_witness, C));
    ASSERT_FALSE(scheme::verify_evaluation(
        i + Field::one(), eval, srs, eval_witness, C));
    ASSERT_FALSE(scheme::verify_evaluation(
        i, eval + Field::one(), srs, eval_witness, C));
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
    const Field eval = scheme::evaluate_polynomial(phi, i);
    const typename scheme::evaluation_witness eval_witness =
        scheme::create_evaluation_witness(phi, i, eval, srs);

    ASSERT_EQ(Field(33), eval);
    ASSERT_EQ(Field(397) * libff::G1<ppT>::one(), eval_witness);
}

template<typename ppT> void test_kzg10_batched_2_point()
{
    using Field = libff::Fr<ppT>;
    using scheme = kzg10<ppT>;
    using batch_scheme = kzg10_batched_2_point<ppT>;

    static const size_t MAX_DEGREE_MULTI = 8;
    const Field secret = Field(7);

    // Generate 2 sets of polynomials
    const std::vector<polynomial<Field>> fs{{
        {{1, 2, 3, 4, 5, 6, 7, 8}},
        {{11, 12, 13, 14, 15, 16, 17, 18}},
        {{21, 22, 23, 24, 25, 26, 27, 28}},
        {{31, 32, 33, 34, 35, 36, 37, 38}},
    }};

    const std::vector<polynomial<Field>> gs{{
        {{71, 72, 73, 74, 75, 76, 77, 78}},
        {{81, 82, 83, 84, 85, 86, 87, 88}},
        {{91, 92, 93, 94, 95, 96, 97, 98}},
    }};

    // srs
    const typename scheme::srs srs =
        scheme::setup_from_secret(MAX_DEGREE_MULTI, secret);

    // commitments
    std::vector<typename scheme::commitment> cm_1s;
    cm_1s.reserve(fs.size());
    for (const polynomial<Field> &f : fs) {
        cm_1s.push_back(scheme::commit(srs, f));
    }
    ASSERT_EQ(fs.size(), cm_1s.size());

    std::vector<typename scheme::commitment> cm_2s;
    cm_2s.reserve(gs.size());
    for (const polynomial<Field> &g : gs) {
        cm_2s.push_back(scheme::commit(srs, g));
    }
    ASSERT_EQ(gs.size(), gs.size());

    // Evaluation points
    const Field z_1("123");
    const Field z_2("456");

    // Evaluations
    const typename batch_scheme::evaluations evaluations =
        batch_scheme::evaluate_polynomials(fs, gs, z_1, z_2);
    ASSERT_EQ(fs.size(), evaluations.s_1s.size());
    ASSERT_EQ(gs.size(), evaluations.s_2s.size());

    // Verifier's random challenges
    const Field gamma_1(54321);
    const Field gamma_2(98760);

    // Witness for evaluations
    const typename batch_scheme::evaluation_witness witness =
        batch_scheme::create_evaluation_witness(
            fs, gs, z_1, z_2, evaluations, srs, gamma_1, gamma_2);

    // Check evaluations are correct.
    {
        // Naively evaluate the h_1(X) and h_2(X) polynomials at the secret x.
        // W_1 and W_2 should be these values encoded in G1.

        Field h_1_x = Field::zero();
        for (size_t i = 0; i < fs.size(); ++i) {
            const polynomial<Field> &f_i = fs[i];
            const Field f_x_minus_f_z_1 =
                libfqfft::evaluate_polynomial(f_i.size(), f_i, secret) -
                libfqfft::evaluate_polynomial(f_i.size(), f_i, z_1);
            const Field gamma_power = gamma_1 ^ i;
            h_1_x += gamma_power * f_x_minus_f_z_1 * ((secret - z_1).inverse());
        }
        ASSERT_EQ(h_1_x * libff::G1<ppT>::one(), witness.W_1);

        Field h_2_x = Field::zero();
        for (size_t i = 0; i < gs.size(); ++i) {
            const polynomial<Field> &g_i = gs[i];
            const Field g_x_minus_g_z_2 =
                libfqfft::evaluate_polynomial(g_i.size(), g_i, secret) -
                libfqfft::evaluate_polynomial(g_i.size(), g_i, z_2);
            const Field gamma_power = gamma_2 ^ i;
            h_2_x += gamma_power * g_x_minus_g_z_2 * ((secret - z_2).inverse());
        }
        ASSERT_EQ(h_2_x * libff::G1<ppT>::one(), witness.W_2);
    }

    // Verify the witnesses
    const Field r = Field(23546);
    ASSERT_TRUE(batch_scheme::verify_evaluations(
        z_1,
        z_2,
        evaluations,
        srs,
        gamma_1,
        gamma_2,
        witness,
        cm_1s,
        cm_2s,
        r));
}

template<typename ppT> void test_for_curve()
{
    // Execute all tests for the given curve.

    ppT::init_public_params();

    // KZG10
    test_basic_commitment<ppT, kzg10<ppT>>();
    test_eval_commitment<ppT, kzg10<ppT>>();
    test_kzg10_commitment_with_known_secret<ppT>();
    test_kzg10_batched_2_point<ppT>();
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
