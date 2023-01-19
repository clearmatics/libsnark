/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_PARAMETERS_HPP_
#define LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_PARAMETERS_HPP_

#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>

namespace libsnark
{

/// Instances of this class expose the following Anemoi parameters for every
/// supported curve
///
/// - prime_field      : is it a prime field (True) or a binary field (False)
/// - mult_generator_g : the smallest generator of the multiplicative
///                      subgroup of the scalar field Fr
/// - alpha            : exponent applied in the Flystel E transformation: E(x)
///                      = x^alpha
/// - alpha_inv        : the inverse of alpha modulo r-1 where r is the
///                      modulus of the scalar field Fr. alpha_inv is
///                      the exponent applied in the inverse mapping
///                      of E: E^{-1}(x) = x^{1/alpha}
/// - beta             : multiplicative constant applied in the quadratic
///                      mappings Q_delta = beta x^quad_exponent + delta and
///                      Q_gamma = beta x^quad_exponent + gamma
/// - gamma             : additive constant applied in the quadratic
///                       mapping Q_gamma = beta x^quad_exponent + gamma
/// - delta             : additive constant applied in the quadratic
///                       mapping Q_delta = beta x^quad_exponent + delta
/// - quad_exponent    : quadratic exponent applied in the mappings Q_gamma,
///                      Q_delta. Note that quad_exponent=2 for prime fields and
///                      quad_exponent=3 for binary fields
/// - C_constants_col_"num_cols", D_constants_col_"num_cols" : the C
///                      and D round constants of the Anemoi
///                      permutation for a state with 1,2,3 or 4
///                      columns with the string placeholder
///                      "num_cols" taking the values resp. "one",
///                      "two", "three" or "four". See [BBCPSVW22],
///                      Sect. 5.1 for more details on the C,D
///                      constants.
///
/// The values for the above parameters for each supported curve were generated
/// with the following Sage script scripts/anemoi-hash/parameters.sage .
template<typename ppT> class anemoi_parameters;

} // namespace libsnark

#include "libsnark/gadgetlib1/gadgets/hashes/anemoi/anemoi_parameters_alt_bn128.tcc"
#include "libsnark/gadgetlib1/gadgets/hashes/anemoi/anemoi_parameters_bls12_377.tcc"
#include "libsnark/gadgetlib1/gadgets/hashes/anemoi/anemoi_parameters_bls12_381.tcc"
#include "libsnark/gadgetlib1/gadgets/hashes/anemoi/anemoi_parameters_bn128.tcc"
#include "libsnark/gadgetlib1/gadgets/hashes/anemoi/anemoi_parameters_bw6_761.tcc"
#include "libsnark/gadgetlib1/gadgets/hashes/anemoi/anemoi_parameters_mnt4.tcc"
#include "libsnark/gadgetlib1/gadgets/hashes/anemoi/anemoi_parameters_mnt6.tcc"

#endif // LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_PARAMETERS_HPP_
