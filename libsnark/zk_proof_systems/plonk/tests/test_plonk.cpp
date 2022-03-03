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
//#include <libff/algebra/curves/bls12_377/bls12_377_pp.hpp>
//#include <libsnark/polynomial_commitments/kzg10.hpp>
//#include <libsnark/polynomial_commitments/tests/polynomial_commitment_test_utils.hpp>

const size_t MAX_DEGREE = 254;

namespace libsnark
{

template<typename ppT> void test_plonk()
{
    // Execute all tests for the given curve.
    ppT::init_public_params();
    printf("[%s:%d] Test OK\n", __FILE__, __LINE__);
}

TEST(TestPolynomialCommitments, ALT_BN128)
{
    test_plonk<libff::alt_bn128_pp>();
}

} // namespace libsnark
