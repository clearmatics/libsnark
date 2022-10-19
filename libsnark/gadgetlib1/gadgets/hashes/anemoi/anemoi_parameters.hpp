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
    using ppT = libff::bls12_381_pp;
    using FieldT = libff::Fr<ppT>;
    static const bool b_prime_field = false;
    static const libff::bigint<FieldT::num_limbs> multiplicative_generator_g;
    static const libff::bigint<FieldT::num_limbs> alpha;
    static const libff::bigint<FieldT::num_limbs> alpha_inv;
    static const libff::bigint<FieldT::num_limbs> beta;
    static const libff::bigint<FieldT::num_limbs> gamma;
    static const libff::bigint<FieldT::num_limbs> delta;
    static const libff::bigint<FieldT::num_limbs> quad_exponent;
};

const libff::bigint<libff::Fr<libff::bls12_381_pp>::num_limbs>
    anemoi_parameters<libff::bls12_381_pp>::multiplicative_generator_g =
        libff::bigint<libff::Fr<libff::bls12_381_pp>::num_limbs>("7");

const libff::bigint<libff::Fr<libff::bls12_381_pp>::num_limbs>
    anemoi_parameters<libff::bls12_381_pp>::alpha =
        libff::bigint<libff::Fr<libff::bls12_381_pp>::num_limbs>("5");

const libff::bigint<libff::Fr<libff::bls12_381_pp>::num_limbs>
    anemoi_parameters<libff::bls12_381_pp>::alpha_inv =
        libff::bigint<libff::Fr<libff::bls12_381_pp>::num_limbs>(
            "209743500700504761917790962032743863350762210002110551290414634799"
            "75432473805");

const libff::bigint<libff::Fr<libff::bls12_381_pp>::num_limbs>
    anemoi_parameters<libff::bls12_381_pp>::beta = multiplicative_generator_g;

const libff::bigint<libff::Fr<libff::bls12_381_pp>::num_limbs>
    anemoi_parameters<libff::bls12_381_pp>::gamma =
        libff::bigint<libff::Fr<libff::bls12_381_pp>::num_limbs>("0");

const libff::bigint<libff::Fr<libff::bls12_381_pp>::num_limbs>
    anemoi_parameters<libff::bls12_381_pp>::delta =
        libff::bigint<libff::Fr<libff::bls12_381_pp>::num_limbs>(
            "14981678621464625851270783002338847382197300714436467949315"
            "331057125308909861");

const libff::bigint<libff::Fr<libff::bls12_381_pp>::num_limbs>
    anemoi_parameters<libff::bls12_381_pp>::quad_exponent =
        libff::bigint<libff::Fr<libff::bls12_381_pp>::num_limbs>("2");

} // namespace libsnark

#endif // LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_PARAMETERS_HPP_
