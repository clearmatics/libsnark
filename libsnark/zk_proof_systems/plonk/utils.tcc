/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_UTILS_TCC_
#define LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_UTILS_TCC_

// Implementation of interfaces for a ppzkSNARK for Plonk. See
// utils.hpp .

namespace libsnark
{

template<typename FieldT> void print_vector(const std::vector<FieldT> &v)
{
    for (size_t i = 0; i < v.size(); ++i) {
        printf("[%2d]: ", (int)i);
        v[i].print();
    }
}

template<typename FieldT>
void plonk_compute_lagrange_basis(
    const size_t npoints, std::vector<polynomial<FieldT>> &L)
{
    assert(L.size() != 0);
    assert(L.size() == L[0].size());
    assert(L.size() == npoints);

    std::shared_ptr<libfqfft::evaluation_domain<FieldT>> domain =
        libfqfft::get_evaluation_domain<FieldT>(npoints);
    for (size_t i = 0; i < npoints; ++i) {
        polynomial<FieldT> u(npoints, FieldT(0));
        u[i] = FieldT(1);
        // compute i-th Lagrange basis vector via inverse FFT
        domain->iFFT(u);
        L[i] = u;
    }
}

template<typename FieldT>
void plonk_interpolate_polynomial_from_points(
    const std::vector<FieldT> &f_points,
    polynomial<FieldT> &f_poly,
    std::shared_ptr<libfqfft::evaluation_domain<FieldT>> domain)
{
    f_poly = f_points;
    domain->iFFT(f_poly);
}

template<typename FieldT>
std::vector<polynomial<FieldT>> plonk_compute_selector_polynomials(
    const size_t &num_gates,
    const size_t &num_qpolys,
    const std::vector<std::vector<FieldT>> &gates_matrix_transpose,
    std::shared_ptr<libfqfft::evaluation_domain<FieldT>> domain)
{
    assert(gates_matrix_transpose.size() == num_qpolys);
    assert(gates_matrix_transpose[0].size() == num_gates);

    std::vector<polynomial<FieldT>> Q_polys;
    Q_polys.resize(num_qpolys, polynomial<FieldT>(num_gates));
    for (size_t i = 0; i < num_qpolys; ++i) {
        std::vector<FieldT> q_vec = gates_matrix_transpose[i];
        plonk_interpolate_polynomial_from_points<FieldT>(
            q_vec, Q_polys[i], domain);
    }
    return Q_polys;
};

template<typename FieldT>
void plonk_compute_public_input_polynomial(
    const std::vector<FieldT> &PI_points,
    polynomial<FieldT> &PI_poly,
    std::shared_ptr<libfqfft::evaluation_domain<FieldT>> domain)
{
    plonk_interpolate_polynomial_from_points<FieldT>(
        PI_points, PI_poly, domain);
};

template<typename FieldT>
void plonk_compute_roots_of_unity_omega(
    const size_t num_gates,
    const FieldT k1,
    const FieldT k2,
    std::vector<std::vector<FieldT>> &omega)
{
    // ensure that num_gates is not 0 and is power of 2
    // TODO: check also that it's less than 2^(ppT::s)
    bool b_nonzero = (num_gates > 0);
    bool b_is_power2 = ((num_gates & (num_gates - 1)) == 0);
    if (!(b_nonzero && b_is_power2)) {
        throw std::invalid_argument(
            "Number of gates not power of 2 or is zero.");
    }
    omega.resize(NUM_HSETS, std::vector<FieldT>(num_gates));

    // Get the n-th root of unity omega in Fq (n=8 is the number of
    // constraints in the example). omega is a generator of the
    // multiplicative subgroup H.  Example (2**32)-th primitive root
    // of unity in the base field Fq of bls12-381 i.e. such that
    // omega_base**(2**32) = 1. The bls12-381 prime q is such that any
    // power of 2 divides (q-1). In particular 2**32|(q-1) and so the
    // 2**32-th root of unity exists.
    FieldT omega_base = libff::get_root_of_unity<FieldT>(num_gates);
    FieldT omega_i = 1;
    for (size_t i = 0; i < num_gates; ++i) {
        omega[base][i] = omega_i;
        FieldT omega_k1_i = omega[base][i] * k1;
        FieldT omega_k2_i = omega[base][i] * k2;
        omega[base_k1][i] = omega_k1_i;
        omega[base_k2][i] = omega_k2_i;
        omega_i *= omega_base;
    }
}

template<typename FieldT>
void plonk_compute_cosets_H_k1H_k2H(
    const size_t num_gates,
    const FieldT k1,
    const FieldT k2,
    std::vector<FieldT> &H_prime)
{
    // ensure that num_gates is not 0 and is power of 2
    bool b_nonzero = (num_gates > 0);
    bool b_is_power2 = ((num_gates & (num_gates - 1)) == 0);
    if (!(b_nonzero && b_is_power2)) {
        throw std::invalid_argument(
            "Number of gates not power of 2 or is zero.");
    }

    // omega[0] are the n roots of unity, omega[1] are omega[0]*k1,
    // omega[2] are omega[0]*k2
    std::vector<std::vector<FieldT>> omega;
    plonk_compute_roots_of_unity_omega(num_gates, k1, k2, omega);

    std::copy(omega[base].begin(), omega[base].end(), back_inserter(H_prime));
    std::copy(
        omega[base_k1].begin(), omega[base_k1].end(), back_inserter(H_prime));
    std::copy(
        omega[base_k2].begin(), omega[base_k2].end(), back_inserter(H_prime));
}

template<typename FieldT>
std::vector<FieldT> plonk_permute_subgroup_H(
    const std::vector<FieldT> &H_prime,
    const std::vector<size_t> &wire_permutation,
    const size_t num_gates)
{
    assert(H_prime.size() > 0);
    std::vector<FieldT> H_prime_permute;
    H_prime_permute.resize(NUM_HSETS * num_gates, FieldT(0));
    for (size_t i = 0; i < H_prime.size(); ++i) {
        H_prime_permute[i] = H_prime[wire_permutation[i] - 1];
    }
    return H_prime_permute;
}

template<typename FieldT>
std::vector<polynomial<FieldT>> plonk_compute_permutation_polynomials(
    const std::vector<FieldT> &H_prime_permute,
    const size_t num_gates,
    std::shared_ptr<libfqfft::evaluation_domain<FieldT>> domain)
{
    assert(H_prime_permute.size() == (NUM_HSETS * num_gates));
    std::vector<polynomial<FieldT>> S_polys;
    S_polys.resize(NUM_HSETS, polynomial<FieldT>(num_gates));
    for (size_t i = 0; i < NUM_HSETS; ++i) {
        typename std::vector<FieldT>::const_iterator begin =
            H_prime_permute.begin() + (i * num_gates);
        typename std::vector<FieldT>::const_iterator end =
            H_prime_permute.begin() + (i * num_gates) + (num_gates);
        std::vector<FieldT> S_points(begin, end);
        plonk_interpolate_polynomial_from_points<FieldT>(
            S_points, S_polys[i], domain);
    }
    return S_polys;
}

template<typename ppT>
libff::G1<ppT> plonk_multi_exp_G1(
    const std::vector<libff::G1<ppT>> &curve_points,
    const std::vector<libff::Fr<ppT>> &scalar_elements)
{
    assert(curve_points.size() == scalar_elements.size());
    const size_t chunks = 1;
    libff::G1<ppT> product = libff::multi_exp<
        libff::G1<ppT>,
        libff::Fr<ppT>,
        libff::multi_exp_method_BDLO12_signed>(
        curve_points.begin(),
        curve_points.end(),
        scalar_elements.begin(),
        scalar_elements.end(),
        chunks);
    return product;
}

template<typename ppT>
libff::G1<ppT> plonk_evaluate_poly_at_secret_G1(
    const std::vector<libff::G1<ppT>> &secret_powers_g1,
    const polynomial<libff::Fr<ppT>> &f_poly)
{
    assert(f_poly.size() <= secret_powers_g1.size());
    const size_t chunks = 1;
    const size_t num_coefficients = f_poly.size();
    libff::G1<ppT> f_poly_at_secret_g1 = libff::multi_exp<
        libff::G1<ppT>,
        libff::Fr<ppT>,
        libff::multi_exp_method_BDLO12_signed>(
        secret_powers_g1.begin(),
        secret_powers_g1.begin() + num_coefficients,
        f_poly.begin(),
        f_poly.end(),
        chunks);
    return f_poly_at_secret_g1;
}

template<typename ppT>
void plonk_evaluate_polys_at_secret_G1(
    const std::vector<libff::G1<ppT>> &secret_powers_g1,
    const std::vector<polynomial<libff::Fr<ppT>>> &Q_polys,
    std::vector<libff::G1<ppT>> &Q_polys_at_secret_g1)
{
    assert(secret_powers_g1.size());
    assert(Q_polys.size());
    assert(Q_polys.size() <= secret_powers_g1.size());
    const size_t npolys = Q_polys.size();
    for (size_t i = 0; i < npolys; ++i) {
        Q_polys_at_secret_g1[i] =
            plonk_evaluate_poly_at_secret_G1<ppT>(secret_powers_g1, Q_polys[i]);
    }
}

template<typename FieldT>
FieldT plonk_compute_accumulator_factor(
    const size_t i,
    const size_t num_gates,
    const FieldT beta,
    const FieldT gamma,
    const std::vector<FieldT> &witness,
    const std::vector<FieldT> &H_prime, // H, Hk1, Hk2
    const std::vector<FieldT> &H_prime_permute,
    const std::vector<FieldT> &A)
{
    assert(num_gates);
    assert(i < num_gates);
    assert(witness.size() == (NUM_HSETS * num_gates));
    assert(H_prime.size() == (NUM_HSETS * num_gates));
    assert(H_prime_permute.size() == (NUM_HSETS * num_gates));
    assert(A.size() == num_gates);
    FieldT res = FieldT(1);
    if (i > 0) {
        FieldT nom_1 = witness[i - 1] + (beta * H_prime[i - 1]) + gamma;
        FieldT den_1 = witness[i - 1] + (beta * H_prime_permute[i - 1]) + gamma;

        FieldT nom_2 = witness[num_gates + i - 1] +
                       (beta * H_prime[num_gates + i - 1]) + gamma;
        FieldT den_2 = witness[num_gates + i - 1] +
                       (beta * H_prime_permute[num_gates + i - 1]) + gamma;

        FieldT nom_3 = witness[2 * num_gates + i - 1] +
                       (beta * H_prime[2 * num_gates + i - 1]) + gamma;
        FieldT den_3 = witness[2 * num_gates + i - 1] +
                       (beta * H_prime_permute[2 * num_gates + i - 1]) + gamma;

        FieldT nom = nom_1 * nom_2 * nom_3;
        FieldT den = den_1 * den_2 * den_3;

        res = nom * den.inverse() * A[i - 1];
    }
    return res;
}

template<typename FieldT>
std::vector<FieldT> plonk_compute_accumulator(
    const size_t num_gates,
    const FieldT beta,
    const FieldT gamma,
    const std::vector<FieldT> &witness,
    const std::vector<FieldT> &H_prime, // H, Hk1, Hk2
    const std::vector<FieldT> &H_prime_permute)
{
    assert(num_gates);
    assert(witness.size() == (NUM_HSETS * num_gates));
    assert(H_prime.size() == (NUM_HSETS * num_gates));
    assert(H_prime_permute.size() == (NUM_HSETS * num_gates));
    std::vector<FieldT> A(num_gates, FieldT(0));
    for (size_t i = 0; i < num_gates; ++i) {
        A[i] = plonk_compute_accumulator_factor(
            i, num_gates, beta, gamma, witness, H_prime, H_prime_permute, A);
    }
    return A;
}

template<typename FieldT>
std::vector<std::vector<FieldT>> plonk_gates_matrix_transpose(
    const std::vector<std::vector<FieldT>> &gates_matrix)
{
    const size_t nrows = gates_matrix.size();
    const size_t ncols = gates_matrix[0].size();
    const size_t nrows_transpose = ncols;
    const size_t ncols_transpose = nrows;
    std::vector<std::vector<FieldT>> gates_matrix_transpose(
        nrows_transpose, std::vector<FieldT>(ncols_transpose));
    for (size_t irow = 0; irow < nrows; ++irow) {
        for (size_t icol = 0; icol < ncols; ++icol) {
            assert(gates_matrix[icol].size() == ncols);
            gates_matrix_transpose[icol][irow] = gates_matrix[irow][icol];
        }
    }
    return gates_matrix_transpose;
}

template<typename FieldT>
void plonk_generate_constants_k1_k2(FieldT &k1_result, FieldT &k2_result)
{
    // n = 2^s: maximum order of the H subgroup that is power of 2
    const size_t n = std::pow(2, FieldT::s);
    FieldT k1, k2;
    // choose k1
    do {
        k1 = FieldT::random_element();
    } while ((k1 ^ n) == FieldT::one());
    // choose k2
    FieldT k1_over_k2;
    do {
        k2 = FieldT::random_element();
        k1_over_k2 = k1 * k2.inverse();
    } while (((k2 ^ n) == FieldT::one()) ||
             (((k1_over_k2) ^ n) == FieldT::one()));
    k1_result = k1;
    k2_result = k2;
}

template<typename FieldT>
bool plonk_are_valid_constants_k1_k2(const FieldT &k1, const FieldT &k2)
{
    // n = 2^s: maximum order of the H subgroup that is power of 2
    const size_t n = std::pow(2, FieldT::s);
    bool b_k1 = ((k1 ^ n) != FieldT::one());
    bool b_k2 = ((k2 ^ n) != FieldT::one());
    bool b_k1_k2 = (((k1 * k2.inverse()) ^ n) != FieldT::one());
    return (b_k1 && b_k2 && b_k1_k2);
}

} // namespace libsnark

#endif // LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_UTILS_TCC_
