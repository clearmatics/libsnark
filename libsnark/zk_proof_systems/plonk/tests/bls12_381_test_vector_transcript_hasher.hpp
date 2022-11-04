/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_TESTS_BLS12_381_TEST_VECTOR_TRANSCRIPT_HASHER_HPP_
#define LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_TESTS_BLS12_381_TEST_VECTOR_TRANSCRIPT_HASHER_HPP_

#include "libsnark/zk_proof_systems/plonk/utils.hpp"

#include <array>

namespace libsnark
{

/// Implementation of the transcript hasher interface (see
/// transcript_hasher.hpp), returning hard-coded test vector values
/// for the plonk proof system on a specific circuit, using the
/// BLS12-381 curve. Also performs a simple check on the transcript to
/// be hashed to verify correct operation of the plonk prover and
/// verifier.
///
/// the hasher works in the Prover as follows:
///
/// round 1
///
///     add_element(buffer, a_eval_exp)
///     add_element(buffer, b_eval_exp)
///     add_element(buffer, c_eval_exp)
///     // buffer = first_output
///
/// round 2
///
///     beta = hash(str(buffer) + "0")
///     gamma = hash(str(buffer) + "1")
///
///     acc_eval = evaluate_in_exponent(CRS, accumulator_poly.to_coeffs())
///
///     add_element(buffer, acc_eval)
///     // buffer = first_output + second_output
///
/// round 3
///
///     alpha = hash(str(buffer))
///
///     add_element(buffer, t_lo_eval_exp)
///     add_element(buffer, t_mid_eval_exp)
///     add_element(buffer, t_hi_eval_exp)
///     // buffer = first_output + second_output + third_output
///
/// round 4
///
///     zeta = hash(str(buffer))
///
///     add_element(buffer, a_zeta)
///     add_element(buffer, b_zeta)
///     add_element(buffer, c_zeta)
///     add_element(buffer, S_0_zeta)
///     add_element(buffer, S_1_zeta)
///     add_element(buffer, accumulator_shift_zeta)
///     add_element(buffer, t_zeta)
///     add_element(buffer, r_zeta)
///     // buffer = first_output + second_output + third_output + fourth_output
///
/// round 5
///
///     nu = hash(str(buffer))
///
///     W_zeta_eval_exp = evaluate_in_exponent(CRS, W_zeta.to_coeffs())
///     W_zeta_omega_eval_exp = evaluate_in_exponent(CRS,
///     W_zeta_omega.to_coeffs())
///
///     add_element(buffer, W_zeta_eval_exp)
///     add_element(buffer, W_zeta_omega_eval_exp)
///
///     // buffer = first_output + second_output + third_output + fourth_output
///     // + fifth_output
///
///     u = hash(str(buffer))
///
class bls12_381_test_vector_transcript_hasher
{
private:
    // buffer accumulating data to be hashed
    std::vector<uint8_t> buffer;
    // array containing the hash values of the communication transcript
    // i.e. the six challenges (in this order): beta, gamma, alpha, zeta, nu, u
    std::array<libff::Fr<libff::bls12_381_pp>, 6> hash_values;
    // vector of valid lengths (\attention specific to BLS12-381)
    const std::array<size_t, 6> length_array;
    // map the index length=0,1...5 to the challenge string=beta,
    // gamma, ...; used to print explicitly the challenge string for debug
    const std::array<std::string, 6> challenge_array;

public:
    bls12_381_test_vector_transcript_hasher();

    // add an Fr element to the transcript buffer for hashing
    void add_element(const libff::Fr<libff::bls12_381_pp> &element);
    // add the coordinates of a G1 curve point to the transcript buffer for
    // hashing
    void add_element(const libff::G1<libff::bls12_381_pp> &element);
    // add the coordinates of a G2 curve point to the transcript buffer for
    // hashing
    void add_element(const libff::G2<libff::bls12_381_pp> &element);

    // dummy implementation of get_hash that directly returns the
    // expected hard-coded hashes for the purposes of unit testing TODO
    // to be replaced by a call to a proper hash function e.g. SHA2,
    // BLAKE, etc.
    libff::Fr<libff::bls12_381_pp> get_hash();

    // clear the buffer (for now only for testing)
    void buffer_clear();

    // get buffer size
    size_t buffer_size();
};

} // namespace libsnark

#endif // LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_TESTS_BLS12_381_TEST_VECTOR_TRANSCRIPT_HASHER_HPP_