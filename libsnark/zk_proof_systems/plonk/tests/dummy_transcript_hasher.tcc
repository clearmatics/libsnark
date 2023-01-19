/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_TESTS_DUMMY_TRANSCRIPT_HASHER_CPP_
#define LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_TESTS_DUMMY_TRANSCRIPT_HASHER_CPP_

#include "libsnark/zk_proof_systems/plonk/tests/dummy_transcript_hasher.hpp"

// Implementation of the dummy transcript hasher interface. See
// dummy_transcript_hasher.hpp.
namespace libsnark
{

template<typename ppT> dummy_transcript_hasher<ppT>::dummy_transcript_hasher()
{
}

template<typename ppT> void dummy_transcript_hasher<ppT>::buffer_clear()
{
    this->buffer.clear();
}

template<typename ppT> size_t dummy_transcript_hasher<ppT>::buffer_size()
{
    return this->buffer.size();
}

template<typename ppT>
void dummy_transcript_hasher<ppT>::add_element(const libff::Fr<ppT> &element)
{
    // convert the Fr element into a string
    std::string str;
    {
        std::ostringstream ss;
        libff::field_write<libff::encoding_binary, libff::form_plain>(
            element, ss);
        str = ss.str();
    }
    // copy the string as a sequence of uint8_t elements at the end of
    // the buffer
    std::copy(str.begin(), str.end(), std::back_inserter(this->buffer));
}

template<typename ppT>
void dummy_transcript_hasher<ppT>::add_element(const libff::G1<ppT> &element)
{
    libff::G1<ppT> element_aff(element);
    element_aff.to_affine_coordinates();

    // convert the affine coordinates of the curve point into a string
    std::string str;
    {
        std::ostringstream ss;
        libff::group_write<
            libff::encoding_binary,
            libff::form_plain,
            libff::compression_off>(element_aff, ss);
        str = ss.str();
    }
    // copy the string as a sequence of uint8_t elements at the end of
    // the buffer
    std::copy(str.begin(), str.end(), std::back_inserter(this->buffer));
}

template<typename ppT>
void dummy_transcript_hasher<ppT>::add_element(const libff::G2<ppT> &element)
{
    libff::G2<ppT> element_aff(element);
    element_aff.to_affine_coordinates();

    // convert the affine coordinates of the curve point into a string
    std::string str;
    {
        std::ostringstream ss;
        libff::group_write<
            libff::encoding_binary,
            libff::form_plain,
            libff::compression_off>(element_aff, ss);
        str = ss.str();
    }
    // copy the string as a sequence of uint8_t elements at the end of
    // the buffer
    std::copy(str.begin(), str.end(), std::back_inserter(this->buffer));
}

template<typename ppT> libff::Fr<ppT> dummy_transcript_hasher<ppT>::get_hash()
{
    libff::Fr<ppT> buffer_len = libff::Fr<ppT>(this->buffer.size());
    return buffer_len;
}

} // namespace libsnark

#endif // LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_TESTS_DUMMY_TRANSCRIPT_HASHER_CPP_
