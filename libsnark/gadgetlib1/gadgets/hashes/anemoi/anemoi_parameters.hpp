/** @file
 *****************************************************************************

 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
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
///
/// The values for the above parameters for each supported curve were generated
/// with the following Sage script scripts/anemoi-hash/parameters.sage .
template<typename ppT> class anemoi_parameters;

template<> class anemoi_parameters<libff::bls12_381_pp>
{
public:
    using FieldT = libff::Fr<libff::bls12_381_pp>;
    using BignumT = libff::bigint<FieldT::num_limbs>;
    static const bool b_prime_field = false;
    static constexpr size_t multiplicative_generator_g = 7;
    static constexpr size_t alpha = 5;
    static constexpr size_t beta = multiplicative_generator_g;
    static constexpr size_t gamma = 0;
    static constexpr size_t quad_exponent = 2;
    static const BignumT alpha_inv;
    static const BignumT delta;
};

} // namespace libsnark

#include "libsnark/gadgetlib1/gadgets/hashes/anemoi/anemoi_parameters.tcc"

#endif // LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_PARAMETERS_HPP_
