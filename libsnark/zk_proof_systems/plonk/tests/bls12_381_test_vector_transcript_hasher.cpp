/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_TESTS_BLS12_381_TEST_VECTOR_TRANSCRIPT_HASHER_CPP_
#define LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_TESTS_BLS12_381_TEST_VECTOR_TRANSCRIPT_HASHER_CPP_

// Implementation of the transcript hasher interface specific to the BLS12-381
// curve. See bls12_381_test_vector_transcript_hasher.hpp
namespace libsnark
{

bls12_381_test_vector_transcript_hasher::
    bls12_381_test_vector_transcript_hasher()
{
    plonk_example example;

    // initialize to empty vector
    this->buffer.clear();
    // test array containing the expected hash values of the communication
    // transcript i.e. the communication challenges (in this order): beta,
    // gamma, alpha, zeta, nu, u WARNING! specific to curve BLS12-381
    this->hash_values = {
        example.beta,
        example.gamma,
        example.alpha,
        example.zeta,
        example.nu,
        example.u,
    };
}

void bls12_381_test_vector_transcript_hasher::buffer_clear()
{
    this->buffer.clear();
}

size_t bls12_381_test_vector_transcript_hasher::buffer_size()
{
    return this->buffer.size();
}

void bls12_381_test_vector_transcript_hasher::add_element(
    const libff::Fr<libff::bls12_381_pp> &element)
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

void bls12_381_test_vector_transcript_hasher::add_element(
    const libff::G1<libff::bls12_381_pp> &element)
{
    libff::G1<libff::bls12_381_pp> element_aff(element);
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

void bls12_381_test_vector_transcript_hasher::add_element(
    const libff::G2<libff::bls12_381_pp> &element)
{
    libff::G2<libff::bls12_381_pp> element_aff(element);
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

libff::Fr<libff::bls12_381_pp> bls12_381_test_vector_transcript_hasher::
    get_hash()
{
    size_t buffer_len = this->buffer.size();

    // vector of valid lengths (\attention specific to BLS12-381)
    const std::vector<size_t> length{288, 320, 416, 704, 896, 1120};

    // If we are here, then the hasher buffer has invalid length so throw an
    // exception
    bool b_valid_length =
        (0 != std::count(length.begin(), length.end(), buffer_len));
    if (!b_valid_length) {
        throw std::logic_error(
            "Error: invalid length of transcript hasher buffer");
    }

    libff::Fr<libff::bls12_381_pp> challenge = 0;

    // map the index length=0,1...5 to the challenge string=beta,
    // gamma, ...; used to print explicitly the challenge string for debug
    std::map<size_t, std::string> challenge_str;
    challenge_str[0] = "beta";
    challenge_str[1] = "gamma";
    challenge_str[2] = "alpha";
    challenge_str[3] = "zeta";
    challenge_str[4] = "nu";
    challenge_str[5] = "u";

    // find the mathcing index
    size_t i = 0;
    while (buffer_len != length[i]) {
        ++i;
        if (i >= length.size()) {
            throw std::logic_error(
                "Error: invalid index of transcript hasher buffer");
        }
    }

    printf(
        "[%s:%d] buffer_len %d: %s\n",
        __FILE__,
        __LINE__,
        (int)buffer_len,
        challenge_str[i].c_str());
    challenge = this->hash_values[i]; // beta

    return challenge;
}

} // namespace libsnark

#endif // LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_TESTS_BLS12_381_TEST_VECTOR_TRANSCRIPT_HASHER_CPP_