/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_SRS_TCC_
#define LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_SRS_TCC_

// Implementation of SRS interfaces for a ppzkSNARK for Plonk. See
// srs.hpp .

namespace libsnark
{

template<typename ppT>
usrs<ppT>::usrs(
    std::vector<libff::G1<ppT>> &&secret_powers_g1,
    std::vector<libff::G2<ppT>> &&secret_powers_g2)
    : secret_powers_g1(secret_powers_g1), secret_powers_g2(secret_powers_g2)
{
}

template<typename ppT>
srs<ppT>::srs(
    const size_t &num_gates,
    const size_t &num_qpolys,
    const std::vector<size_t> &PI_wire_indices,
    const std::vector<polynomial<Field>> &Q_polys,
    const std::vector<polynomial<Field>> &S_polys,
    const std::vector<std::vector<Field>> &omega_roots,
    const std::vector<Field> &H_prime,
    const std::vector<Field> &H_prime_permute,
    const libff::Fr<ppT> &k1,
    const libff::Fr<ppT> &k2,
    std::vector<libff::G1<ppT>> &&secret_powers_g1,
    std::vector<libff::G2<ppT>> &&secret_powers_g2,
    const polynomial<Field> &L_basis_zero)
    : num_gates(num_gates)
    , num_qpolys(num_qpolys)
    , PI_wire_indices(PI_wire_indices)
    , Q_polys(Q_polys)
    , S_polys(S_polys)
    , omega_roots(omega_roots)
    , H_prime(H_prime)
    , H_prime_permute(H_prime_permute)
    , k1(k1)
    , k2(k2)
    , secret_powers_g1(secret_powers_g1)
    , secret_powers_g2(secret_powers_g2)
    , L_basis_zero(L_basis_zero)
{
}

template<typename ppT>
plonk_verification_key<ppT>::plonk_verification_key(
    std::vector<libff::G2<ppT>> &&secret_powers_g2)
    : secret_powers_g2(std::move(secret_powers_g2)){};

template<typename ppT>
plonk_keypair<ppT>::plonk_keypair(
    plonk_proving_key<ppT> &&pk, plonk_verification_key<ppT> &&vk)
    : pk(std::move(pk)), vk(std::move(vk))
{
}

template<typename ppT>
usrs<ppT> plonk_usrs_derive_from_secret(
    const libff::Fr<ppT> &secret, const size_t max_degree)
{
    // compute powers of secret times G1: 1*G1, secret^1*G1, secret^2*G1, ...
    const libff::bigint<libff::Fr<ppT>::num_limbs> secret_bigint =
        secret.as_bigint();
    const size_t window_size = std::max(
        libff::wnaf_opt_window_size<libff::G1<ppT>>(secret_bigint.num_bits()),
        1ul);
    const std::vector<long> naf =
        libff::find_wnaf<libff::Fr<ppT>::num_limbs>(window_size, secret_bigint);

    std::vector<libff::G1<ppT>> secret_powers_g1;
    secret_powers_g1.reserve(max_degree);
    libff::G1<ppT> secret_i_g1 = libff::G1<ppT>::one();
    secret_powers_g1.push_back(secret_i_g1);
    for (size_t i = 1; i < max_degree; ++i) {
        // secret^i * G1
        secret_i_g1 = libff::fixed_window_wnaf_exp<libff::G1<ppT>>(
            window_size, secret_i_g1, naf);
        secret_powers_g1.push_back(secret_i_g1);
    }

    // compute powers of secret times G2: 1*G2, secret^1*G2
    // Note: in Plonk we *always* have 2 encoded elemnts in G2
    std::vector<libff::G2<ppT>> secret_powers_g2;
    secret_powers_g2.reserve(2);
    // secret^0 * G2 = G2
    libff::G2<ppT> secret_0_g2 = libff::G2<ppT>::one();
    secret_powers_g2.push_back(secret_0_g2);
    // secret^1 * G2
    libff::G2<ppT> secret_1_g2 = secret * libff::G2<ppT>::one();
    secret_powers_g2.push_back(secret_1_g2);

    return usrs<ppT>(std::move(secret_powers_g1), std::move(secret_powers_g2));
}

template<typename ppT>
srs<ppT> plonk_srs_derive_from_usrs(
    const usrs<ppT> &usrs,
    const std::vector<std::vector<libff::Fr<ppT>>> gates_matrix,
    const std::vector<size_t> wire_permutation,
    const std::vector<size_t> PI_wire_indices)
{
    using Field = libff::Fr<ppT>;
    const std::vector<std::vector<Field>> gates_matrix_transpose =
        plonk_gates_matrix_transpose(gates_matrix);

    // the number of gates is equal to the number of columns in the transposed
    // gates matrix
    size_t num_gates = gates_matrix_transpose[0].size();
    // ensure that num_gates is not 0
    assert(num_gates > 0);
    // ensure num_gates is power of 2
    assert((num_gates & (num_gates - 1)) == 0);

    // the number of Q-polynomials (aka selector polynomials) is equal to the
    // number of rows in the transposed gates matrix
    size_t num_qpolys = gates_matrix_transpose.size();

    // the constraints q_L, q_R, q_O, q_M, q_C and the
    // witness w_L, w_R, w_O are represented as polynomials in the roots of
    // unity e.g. f_{q_L}(omega_i) = q_L[i], 0\le{i}<8
    std::shared_ptr<libfqfft::evaluation_domain<Field>> domain =
        libfqfft::get_evaluation_domain<Field>(num_gates);

    // compute the selector polynomials (q-polynomials) from the
    // transposed gates matrix over the Lagrange basis q_poly = \sum_i
    // q[i] * L[i] where q[i] is a coefficient (a scalar Field
    // element) and L[i] is a polynomial with Field coefficients
    std::vector<polynomial<Field>> Q_polys =
        plonk_compute_selector_polynomials<Field>(
            num_gates, num_qpolys, gates_matrix_transpose, domain);

    // An explanation of the constants k1,k2 from [GWC19], Section 8,
    // page 26 :
    //
    // "We explicitly define the multiplicative subgroup H as
    // containing the n-th roots of unity in F_p , where w (omega) is
    // a primitive n-th root of unity and a generator of H i.e: H =
    // {1, w, ... , w^{n-1}}. We assume that the number of gates in a
    // circuit is no more than n.
    //
    // For the moment k1,k2 are fixed (see below) to the test vector
    // values from the plonk_example class for debug purpouses. Note
    // that these test values are specific to the BLS12-381 curve and
    // hence they satisfy the requirements for BLS12-381. TODO: choose
    // k1,k2 according to the requirements in [GWC19]
    libff::Fr<ppT> k1 =
        Field("706987411474581393682955260879121390206111740035659671471"
              "3673571023200548519");
    libff::Fr<ppT> k2 = libff::power(k1, libff::bigint<1>(2));

    // From [GWC19], Section 8, page 26:
    //
    // "We also include an optimisation suggested by Vitalik Buterin,
    // to define the identity permutations through degree-1
    // polynomials. The identity permutations must map each wire value
    // to a unique element \in F. This can be done by defining
    // S_ID1(X) = X, S_ID2 (X) = k1 X, S_ID3(X) = k2 X [see below for
    // more on S_ID1, S_ID2, S_ID3], where k1, k2 are quadratic
    // non-residues \in F. This effectively maps each wire value to a
    // root of unity in H, with right and output wires having an
    // additional multiplicative factor of k1, k2 applied
    // respectively. By representing the identity permutation via
    // degree-1 polynomials, their evaluations can be directly
    // computed by the verifier. This reduces the size of the proof by
    // 1 F element, as well as reducing the number of
    // Fast-Fourier-Transforms required by the prover."
    //
    // Further in Sect. 8.1 [GWC19]:
    //
    // "S_ID1(X) = X, S_ID2(X) = k 1 X, S ID3 (X) = k 2 X: the
    // identity permutation applied to a, b, c [the wire polynomials,
    // see Round 1, p.27 [GWC19]]. k1, k2 \in F are chosen such that
    // H, k1 H, k2 H are distinct cosets of H in F*, and thus consist
    // of 3n distinct elements. (For example, when w (omega) is a
    // quadratic residue in F, take k1 to be any quadratic
    // non-residue, and k2 to be a quadratic non-residue not contained
    // in k1 H.)"

    // omega[0] are the n roots of unity; omega[1] are omega[0]*k1;
    // omega[2] are omega[0]*k2
    std::vector<std::vector<Field>> omega_roots;
    plonk_compute_roots_of_unity_omega(num_gates, k1, k2, omega_roots);

    // H_prime contains the generators of H, k1 H and k2 H in one
    // place ie. omega_roots, omega_roots_k1 and omega_roots_k2
    std::vector<Field> H_prime;
    plonk_compute_cosets_H_k1H_k2H(num_gates, k1, k2, H_prime);

    // TODO: write unit test for plonk_compute_cosets_H_k1H_k2H
    // assert(H_prime[i] == example.H_prime[i]);

    // permute H_prime according to the wire permutation
    std::vector<Field> H_prime_permute =
        plonk_permute_subgroup_H<Field>(H_prime, wire_permutation, num_gates);

    // TODO: write unit test for plonk_permute_subgroup_H
    // assert(H_prime_permute[i] == example.H_prime_permute[i]);

    // Compute the permutation polynomials S_sigma_1, S_sigma_2,
    // S_sigma_3 (see [GWC19], Sect. 8.1) (our indexing starts from 0)
    std::vector<polynomial<Field>> S_polys =
        plonk_compute_permutation_polynomials<Field>(
            H_prime_permute, num_gates, domain);

    // secret^i * G1
    std::vector<libff::G1<ppT>> secret_powers_g1(
        usrs.secret_powers_g1.begin(),
        usrs.secret_powers_g1.begin() + num_gates + 3);

    for (size_t i = 0; i < (num_gates + 3); ++i) {
        secret_powers_g1.push_back(usrs.secret_powers_g1[i]);
    }
    // secret^i * G2
    std::vector<libff::G2<ppT>> secret_powers_g2;
    secret_powers_g2.reserve(2);
    for (size_t i = 0; i < 2; ++i) {
        secret_powers_g2.push_back(usrs.secret_powers_g2[i]);
    }

    // compute 0-th Lagrange basis vector via inverse FFT
    polynomial<libff::Fr<ppT>> L_basis_zero(num_gates, libff::Fr<ppT>(0));
    L_basis_zero[0] = libff::Fr<ppT>(1);
    domain->iFFT(L_basis_zero);

    srs<ppT> srs(
        std::move(num_gates),
        std::move(num_qpolys),
        std::move(PI_wire_indices),
        std::move(Q_polys),
        std::move(S_polys),
        std::move(omega_roots),
        std::move(H_prime),
        std::move(H_prime_permute),
        std::move(k1),
        std::move(k2),
        std::move(secret_powers_g1),
        std::move(secret_powers_g2),
        std::move(L_basis_zero));

    return srs;
}

} // namespace libsnark

#endif // LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_SRS_TCC_
