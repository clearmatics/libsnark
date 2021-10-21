/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef __LIBSNARK_POLYNOMIAL_COMMITMENTS_KZG10_TCC__
#define __LIBSNARK_POLYNOMIAL_COMMITMENTS_KZG10_TCC__

#include "libsnark/polynomial_commitments/kzg10.hpp"

#include <libff/algebra/scalar_multiplication/multiexp.hpp>
#include <libff/algebra/scalar_multiplication/wnaf.hpp>
#include <libfqfft/polynomial_arithmetic/basic_operations.hpp>
#include <libfqfft/polynomial_arithmetic/naive_evaluate.hpp>

namespace libsnark
{

template<typename ppT, typename scheme>
evaluation_and_witness<ppT, scheme>::evaluation_and_witness(
    const libff::Fr<ppT> &i,
    const libff::Fr<ppT> &phi_i,
    const typename scheme::witness &w)
    : i(i), phi_i(phi_i), w(w)
{
}

template<typename ppT>
kzg10<ppT>::srs::srs(
    std::vector<libff::G1<ppT>> &&alpha_powers_g1,
    const libff::G2<ppT> &alpha_g2)
    : alpha_powers_g1(alpha_powers_g1), alpha_g2(alpha_g2)
{
}

template<typename ppT>
typename kzg10<ppT>::srs kzg10<ppT>::setup_from_secret(
    size_t max_degree, const Field &alpha)
{
    const libff::bigint<Field::num_limbs> alpha_bigint = alpha.as_bigint();
    const size_t window_size = std::max(
        libff::wnaf_opt_window_size<libff::G1<ppT>>(alpha_bigint.num_bits()),
        1ul);
    const std::vector<long> naf = libff::find_wnaf<Field::num_limbs>(
        window_size, alpha_bigint);

    // TODO: perform in concurrent batches?
    std::vector<libff::G1<ppT>> alpha_powers_g1;
    alpha_powers_g1.reserve(max_degree + 1);
    libff::G1<ppT> alpha_i_g1 = libff::G1<ppT>::one();
    alpha_powers_g1.push_back(alpha_i_g1);
    for (size_t i = 1 ; i < max_degree + 1 ; ++i)
    {
        alpha_i_g1 = libff::fixed_window_wnaf_exp<libff::G1<ppT>>(
            window_size, alpha_i_g1, naf);
        alpha_powers_g1.push_back(alpha_i_g1);
    }

    // assert((libff::G1<ppT> * alpha.pow(max_degree)) == alpha_powers_g1[max_degree]);
    return srs(std::move(alpha_powers_g1), alpha * libff::G2<ppT>::one());
}

template<typename ppT>
typename kzg10<ppT>::srs kzg10<ppT>::setup(size_t max_degree)
{
    const Field alpha = Field::random_element();
    return setup_from_secret(max_degree, alpha);
}

template<typename ppT>
typename kzg10<ppT>::commitment kzg10<ppT>::commit(
    const srs &srs, const polynomial<libff::Fr<ppT>> &phi)
{
    const size_t num_coefficients = phi.size();

    // Caller is responsible for checking this.
    assert(num_coefficients <= srs.alpha_powers_g1.size());

#ifdef MULTICORE
    const size_t chunks = omp_get_max_threads();
#else
    const size_t chunks = 1;
#endif

    // The commitment is the encoding [ \phi(\alpha) ]_1. This is computed
    // using the elements [ \alpha^t ]_1 for t=0, ..., max_degree from the srs.
    return libff::multi_exp<
        libff::G1<ppT>,
        libff::Fr<ppT>,
        libff::multi_exp_method_BDLO12_signed>(
            srs.alpha_powers_g1.begin(),
            srs.alpha_powers_g1.begin() + num_coefficients,
            phi.begin(),
            phi.end(),
            chunks);
}

template<typename ppT>
bool kzg10<ppT>::verify_poly(
    const typename kzg10<ppT>::srs &srs,
    typename kzg10<ppT>::commitment C,
    const polynomial<libff::Fr<ppT>> &phi)
{
    return C == commit(srs, phi);
}

template<typename ppT>
evaluation_and_witness<ppT, kzg10<ppT>> kzg10<ppT>::create_witness(
    const typename kzg10<ppT>::srs &srs,
    const polynomial<libff::Fr<ppT>> &phi,
    const libff::Fr<ppT> i)
{
    const size_t num_coefficients = phi.size();
    assert(num_coefficients > 1);
    assert(num_coefficients <= srs.alpha_powers_g1.size());

    // Evaluate the polynomial
    const Field phi_i = libfqfft::evaluate_polynomial(num_coefficients, phi, i);

    // The witness is:
    //   w_i = [ \psi_i(\alpha) ]_1
    // where:
    //   \psi_i(x) = (\phi(x) - \phi(i)) / (x - i)
    // (see [KZG10] Section 3.2)

    // TODO: Find a more optimal way to compute the special case of dividing by
    // (x - i).
    std::vector<Field> psi;
    std::vector<Field> remainder;

    std::vector<Field> phi_minus_phi_i = phi;
    phi_minus_phi_i[0] -= phi_i;
    libfqfft::_polynomial_division(psi, remainder, phi_minus_phi_i, {-i, 1});
    assert(libfqfft::_is_zero(remainder));
    libff::G1<ppT> w = commit(srs, psi);

    return evaluation_and_witness<ppT, kzg10<ppT>>(i, phi_i, w);
}

template<typename ppT>
bool kzg10<ppT>::verify_eval(
    const typename kzg10<ppT>::srs &srs,
    typename kzg10<ppT>::commitment C,
    const evaluation_and_witness<ppT, kzg10<ppT>> &witness)
{
    // Verify the equality:
    //
    //   \psi(\alpha) (\alpha - i) = \phi(\alpha) - \phi(i)
    //                             = C - phi_i
    // via the pairing equality:
    //
    //   e([\psi(\alpha)]_1, [\alpha - i]_2) = e(C - [phi_i]_1, [1]_2)
    //
    // which (multiplying both sides by e(C - [phi_i]_1, [1]_2)^{-1}) holds iff:
    //
    //   id_T
    //   = e([\psi(\alpha)]_1, [\alpha - i]_2) * e(C - [phi_i]_1, [1]_2)^{-1}
    //   = e([\psi(\alpha)]_1, [\alpha - i]_2) * e([phi_i]_1 - C, [1]_2)
    //   = e([\psi(\alpha)]_1, [\alpha]_2 - [i]_2) * e([phi_i]_1 - C, [1]_2)
    //
    // where id_T is the identity of the group GT (GT::one()). This
    // last line can now be computed in terms of available variables:
    //
    //   e(_A = witness.w, _B = srs.alpha_g2 - i * G2::one()) *
    //   e(_C = witness.phi_i * G1::one() - C, _D = G2::one())

    // See Section 3: https://eprint.iacr.org/2019/953.pdf for further
    // details.

    const libff::G1_precomp<ppT> _A = ppT::precompute_G1(witness.w);
    const libff::G2_precomp<ppT> _B = ppT::precompute_G2(
        srs.alpha_g2 - witness.i * libff::G2<ppT>::one());
    const libff::G1_precomp<ppT> _C = ppT::precompute_G1(
        witness.phi_i * libff::G1<ppT>::one() - C);
    const libff::G2_precomp<ppT> _D = ppT::precompute_G2(libff::G2<ppT>::one());

    const libff::Fqk<ppT> miller_result = ppT::double_miller_loop(_A, _B, _C, _D);
    const libff::GT<ppT> result = ppT::final_exponentiation(miller_result);
    return result == libff::GT<ppT>::one();
}

} // namespace libsnark

#endif // __LIBSNARK_POLYNOMIAL_COMMITMENTS_KZG10_TCC__
