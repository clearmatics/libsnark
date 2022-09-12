/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_TRANSCRIPT_HASHER_TRANSCRIPT_HASHER_HPP_
#define LIBSNARK_TRANSCRIPT_HASHER_TRANSCRIPT_HASHER_HPP_

// // interface for a common transcript_hasher class used to implement
// // functionality for hashing the communication transcript in ZK proof
// // systems under ./zk_proof_systems
// template<typename ppT> class transcript_hasher
// {
// public:
//     transcript_hasher();
//
//     // add an Fr element to the transcript buffer for hashing
//     void add_element(const libff::Fr<ppT> &element);
//     // add the coordinates of a G1 curve point to the transcript buffer for
//     // hashing
//     void add_element(const libff::G1<ppT> &element);
//     // add the coordinates of a G2 curve point to the transcript buffer for
//     // hashing
//     void add_element(const libff::G2<ppT> &element);
//     // return the hash value of the communication transcript
//     libff::Fr<ppT> get_hash();
// };

#endif // LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_SRS_HPP_
