/** @file
 *****************************************************************************

 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_PARAMETERS_TCC_
#define LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_PARAMETERS_TCC_

namespace libsnark
{

const libff::bigint<libff::Fr<libff::bls12_381_pp>::num_limbs>
    anemoi_parameters<libff::bls12_381_pp>::alpha_inv =
        libff::bigint<libff::Fr<libff::bls12_381_pp>::num_limbs>(
            "209743500700504761917790962032743863350762210002110551290414634799"
            "75432473805");

const libff::bigint<libff::Fr<libff::bls12_381_pp>::num_limbs>
    anemoi_parameters<libff::bls12_381_pp>::delta =
        libff::bigint<libff::Fr<libff::bls12_381_pp>::num_limbs>(
            "14981678621464625851270783002338847382197300714436467949315"
            "331057125308909861");

} // namespace libsnark

#endif // LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_PARAMETERS_TCC_
