/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_TESTS_ANEMOI_OUTPUTS_HPP_
#define LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_TESTS_ANEMOI_OUTPUTS_HPP_

#include <functional>
#include <libff/algebra/curves/bls12_381/bls12_381_init.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_pp.hpp>
#include <vector>

// Functions returning the expected outputs from Anemoi for each
// tested curve. An implementation of each function should be provided for each
// tested curve.

namespace libsnark
{

// Returns the expected outputs from 1 round of the Anemoi permutation for
// BLS12_381
std::vector<libff::Fr<libff::bls12_381_pp>> anemoi_expected_output_one_round(
    const size_t &NumStateColumns_L);

template<typename ppT>
using expected_round_values_fn_t =
    std::function<std::vector<libff::Fr<ppT>>(const size_t)>;

// Returns the expected outputs from the full Anemoi permutation for
// BLS12_381
std::vector<libff::Fr<libff::bls12_381_pp>> anemoi_expected_output(
    const size_t &NumStateColumns_L);

template<typename ppT>
using expected_values_fn_t =
    std::function<std::vector<libff::Fr<ppT>>(const size_t)>;

} // namespace libsnark

#endif // LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_TESTS_ANEMOI_OUTPUTS_HPP_
