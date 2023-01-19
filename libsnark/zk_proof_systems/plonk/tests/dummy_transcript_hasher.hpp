/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_TESTS_DUMMY_TRANSCRIPT_HASHER_HPP_
#define LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_TESTS_DUMMY_TRANSCRIPT_HASHER_HPP_

#include "libsnark/zk_proof_systems/plonk/utils.hpp"

#include <array>

namespace libsnark
{

/// Implementation of a dummy transcript hasher interface (see
/// transcript_hasher.hpp). It returns the number of the elemnts in
/// the hash buffer as an Fr element. Specialized over the curve
/// field. See also class bls12_381_test_vector_transcript_hasher,
/// which is specific to the BLS12_381 curve.
template<typename ppT> class dummy_transcript_hasher
{
private:
    // buffer accumulating data to be hashed
    std::vector<uint8_t> buffer;

public:
    dummy_transcript_hasher();

    // Add an Fr element to the transcript buffer for hashing.
    void add_element(const libff::Fr<ppT> &element);
    // Add the coordinates of a G1 curve point to the transcript buffer for
    // hashing.
    void add_element(const libff::G1<ppT> &element);
    // Add the coordinates of a G2 curve point to the transcript buffer for
    // hashing.
    void add_element(const libff::G2<ppT> &element);

    // Dummy implementation of get_hash that simply returns the number
    // of elements in the buffer as an Fr value for the purposes of
    // unit testing. TODO: to be replaced by a call to a proper hash
    // function e.g. SHA2, BLAKE, etc.
    libff::Fr<ppT> get_hash();

    // clear the buffer (for now only for testing)
    void buffer_clear();

    // get buffer size
    size_t buffer_size();
};

} // namespace libsnark

#include "libsnark/zk_proof_systems/plonk/tests/dummy_transcript_hasher.tcc"

#endif // LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_DUMMY_TRANSCRIPT_HASHER_HPP_
