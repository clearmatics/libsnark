/** @file
*****************************************************************************
Implementation of Prover interfaces for a ppzkSNARK for Plonk.

See prover.hpp .

*****************************************************************************
* @author     This file is part of libsnark, developed by SCIPR Lab
*             and contributors (see AUTHORS).
* @copyright  MIT license (see LICENSE file)
*****************************************************************************/

#ifndef PLONK_PPZKSNARK_PROVER_TCC_
#define PLONK_PPZKSNARK_PROVER_TCC_

namespace libsnark
{

// Prover Round 0 initialization
//
// Initialization
//
// INPUT
// - srs: structured reference string containing also circuit-specific
//   information
//
// OUTPUT
// - W_polys: Lagrange interpolation of the witness values
// - zh_poly: vanishing polynomial
// - null_poly: 0 polynomial
// - neg_one_poly: -1 polynomial
template<typename ppT>
round_zero_out_t<ppT> plonk_prover<ppT>::round_zero(const srs<ppT> srs)
{
    using Field = libff::Fr<ppT>;

    // initialize hard-coded values from example circuit
#ifdef DEBUG
    plonk_example<ppT> example;
#endif // #ifdef DEBUG

    // output from round 0
    round_zero_out_t<ppT> round_zero_out;

    // vanishing polynomial zh_poly(X) = x^n-1. vanishes on all n roots of
    // unity srs.omega_roots
    round_zero_out.zh_poly.resize(srs.num_gates + 1, Field(0));
    round_zero_out.zh_poly[0] = Field(-1);
    round_zero_out.zh_poly[srs.num_gates] = Field(1);
#ifdef DEBUG
    printf("[%s:%d] Vanishing polynomial\n", __FILE__, __LINE__);
    print_vector(round_zero_out.zh_poly);
    assert(round_zero_out.zh_poly == example.zh_poly);
#endif // #ifdef DEBUG

    round_zero_out.null_poly = {Field(0)};
    round_zero_out.neg_one_poly = {-Field("1")};

    return round_zero_out;
}

// Prover Round 1
//
// INPUT
// - zh_poly: vanishing polynomial Zh (from round 0)
// - null_poly: 0 polynomial (from round 0)
// - neg_one_poly: -1 polynomial (from round 0)
// - witness: witness values
// - srs: structured reference string containing also circuit-specific
//   information
//
// OUTPUT
// - blind_scalars: blinding scalars b1, b2, ..., b9 (only
//   b1-b6 used in round 1)
// - W_polys: witness polynomials (Lagrange interpolation of the
//   witness values)
// - W_polys_blinded: blinded witness polynomials
// - W_polys_blinded_at_secret_g1: the blinded witness polynomials
//   evaluated at the secret input denoted [a]_1, [b]_1, [c]_1 in
//   [GWC19]
template<typename ppT>
round_one_out_t<ppT> plonk_prover<ppT>::round_one(
    const round_zero_out_t<ppT> round_zero_out,
    const std::vector<libff::Fr<ppT>> witness,
    const srs<ppT> srs)
{
    using Field = libff::Fr<ppT>;
    int nwitness = 3;

    // initialize hard-coded values from example circuit
    plonk_example<ppT> example;

    // output from round 1
    round_one_out_t<ppT> round_one_out;

    // compute witness polynomials via Lagrange interpolation
    round_one_out.W_polys.resize(nwitness, polynomial<Field>(srs.num_gates));
    for (int i = 0; i < nwitness; ++i) {
        typename std::vector<Field>::const_iterator begin =
            witness.begin() + (i * srs.num_gates);
        typename std::vector<Field>::const_iterator end =
            witness.begin() + (i * srs.num_gates) + (srs.num_gates);
        std::vector<Field> W_points(begin, end);
        plonk_interpolate_polynomial_from_points<Field>(
            W_points, round_one_out.W_polys[i]);
    }
#ifdef DEBUG
    for (int i = 0; i < nwitness; ++i) {
        printf("[%s:%d] this->W_polys[%d]\n", __FILE__, __LINE__, i);
        print_vector(round_one_out.W_polys[i]);
        assert(round_one_out.W_polys[i] == example.W_polys[i]);
    }
#endif // #ifdef DEBUG

    // hard-coded values for the "random" blinding constants from
    // example circuit
    round_one_out.blind_scalars = example.prover_blind_scalars;
    // represent the blinding scalars b1, b2, ..., b9 as polynomials
#ifdef DEBUG
    printf("[%s:%d] blind_scalars\n", __FILE__, __LINE__);
    print_vector(round_one_out.blind_scalars);
    assert(round_one_out.blind_scalars == example.prover_blind_scalars);
#endif // #ifdef DEBUG
    std::vector<std::vector<Field>> blind_polys{
        {round_one_out.blind_scalars[1],
         round_one_out.blind_scalars[0]}, // b1 + b0 X
        {round_one_out.blind_scalars[3],
         round_one_out.blind_scalars[2]}, // b3 + b2 X
        {round_one_out.blind_scalars[5],
         round_one_out.blind_scalars[4]} // b5 + b4 X
    };

    // compute blinded witness polynomials e.g. a_poly =
    // blind_polys[0] * zh_poly + W_polys[0]
    round_one_out.W_polys_blinded.resize(nwitness);
    for (int i = 0; i < nwitness; ++i) {
        libfqfft::_polynomial_multiplication<Field>(
            round_one_out.W_polys_blinded[i],
            blind_polys[i],
            round_zero_out.zh_poly);
        libfqfft::_polynomial_addition<Field>(
            round_one_out.W_polys_blinded[i],
            round_one_out.W_polys_blinded[i],
            round_one_out.W_polys[i]);
    }
    // evaluate blinded witness polynomials at the secret input
    round_one_out.W_polys_blinded_at_secret_g1.resize(
        round_one_out.W_polys_blinded.size());
    plonk_evaluate_polys_at_secret_G1<ppT>(
        srs.secret_powers_g1,
        round_one_out.W_polys_blinded,
        round_one_out.W_polys_blinded_at_secret_g1);

#ifdef DEBUG
    printf(
        "[%s:%d] poly %d %d %d secret %d\n",
        __FILE__,
        __LINE__,
        (int)round_one_out.W_polys_blinded[0].size(),
        (int)round_one_out.W_polys_blinded[1].size(),
        (int)round_one_out.W_polys_blinded[2].size(),
        (int)srs.secret_powers_g1.size());
#endif // #ifdef DEBUG

    return round_one_out;
}

// Prover Round 2
//
// INPUT
// - blind_scalars: blinding scalars b1, b2, ..., b9 (only
//   b7,b8,b9 used in round 2) (from round 1)
// - zh_poly: vanishing polynomial Zh (from round 0)
// - witness: witness values
// - srs: structured reference string containing also circuit-specific
//   information
//
// OUTPUT
// - beta, gamma: permutation challenges -- hashes of transcript
// - z_poly: blinded accumulator poly z(x)
// - z_poly_at_secret_g1: blinded accumulator poly z(x) evaluated at
// secret
//
template<typename ppT>
round_two_out_t<ppT> plonk_prover<ppT>::round_two(
    const round_zero_out_t<ppT> round_zero_out,
    const round_one_out_t<ppT> round_one_out,
    const std::vector<libff::Fr<ppT>> witness,
    const srs<ppT> srs)
{
    using Field = libff::Fr<ppT>;
    // initialize hard-coded values from example circuit
    plonk_example<ppT> example;

    // output from round 2
    round_two_out_t<ppT> round_two_out;

    // permutation challenges (hashes of transcript); fixed to test
    // values
    round_two_out.beta = example.beta;
    round_two_out.gamma = example.gamma;

    // compute permutation polynomial

    // blinding polynomial: b8 + b7 X + b6 X^2
    std::vector<Field> z1_blind_poly{round_one_out.blind_scalars[8],
                                     round_one_out.blind_scalars[7],
                                     round_one_out.blind_scalars[6]};
    // multiply by the vanishing polynomial: z1 = z1 * this->zh_poly
    libfqfft::_polynomial_multiplication<Field>(
        z1_blind_poly, z1_blind_poly, round_zero_out.zh_poly);
#ifdef DEBUG
    printf("[%s:%d] z1_blind_poly * zh_poly\n", __FILE__, __LINE__);
    print_vector(z1_blind_poly);
#endif // #ifdef DEBUG

    // A[0] = 1; ... A[i] = computed from (i-1)
    std::vector<Field> A_vector(srs.num_gates, Field(0));
    plonk_compute_accumulator(
        srs.num_gates,
        round_two_out.beta,
        round_two_out.gamma,
        witness,
        srs.H_gen,
        srs.H_gen_permute,
        A_vector);
#ifdef DEBUG
    for (int i = 0; i < (int)srs.num_gates; ++i) {
        printf("A[%d] ", i);
        A_vector[i].print();
    }
#endif // #ifdef DEBUG

    polynomial<Field> A_poly(srs.num_gates);
    plonk_interpolate_polynomial_from_points<Field>(A_vector, A_poly);
#ifdef DEBUG
    printf("[%s:%d] A_poly\n", __FILE__, __LINE__);
    print_vector(A_poly);
    assert(A_poly == example.A_poly);
#endif // #ifdef DEBUG

    // add blinding polynomial z_1 to the accumulator polynomial A_poly
    libfqfft::_polynomial_addition<Field>(
        round_two_out.z_poly, z1_blind_poly, A_poly);
    round_two_out.z_poly_at_secret_g1 = plonk_evaluate_poly_at_secret_G1<ppT>(
        srs.secret_powers_g1, round_two_out.z_poly);

    return round_two_out;
}

// Prover Round 3
//
// Warning! The most computationally intensive part of the prover
// routine
//
// INPUT
// - zh_poly: vanishing polynomial Zh (from Round 0)
// - W_polys_blinded: blinded witness polynomials (from Round 1)
// - beta, gamma: permutation challenges -- hashes of transcript (from
//   round 2)
// - z_poly: blinded accumulator poly z(x) (from Round 2)
// - srs: structured reference string containing also circuit-specific
//   information
//
// OUTPUT
// - alpha: quotinet challenge -- hash of transcript
// - t_poly_long: the quotient polynomial t(x) (see Round 3, pp28
//   [GWC19])
// - t_poly: t(x) divided in three parts t(x) = t_lo(x) + t_mid(x) x^n
//   + t_hi(x) x^{2n}
// - t_poly_at_secret_g1: t(x) evaluated at the secret input zeta
//   i.e. t(zeta)
// - z_poly_xomega: the polynomial z(x*w) i.e. z(x) shifted by w
//
template<typename ppT>
round_three_out_t<ppT> plonk_prover<ppT>::round_three(
    const round_zero_out_t<ppT> round_zero_out,
    const round_one_out_t<ppT> round_one_out,
    const round_two_out_t<ppT> round_two_out,
    const srs<ppT> srs)
{
    using Field = libff::Fr<ppT>;
    int num_hgen = NUM_HGEN;

    // initialize hard-coded values from example circuit
    plonk_example<ppT> example;

    // output from round 3
    round_three_out_t<ppT> round_three_out;

    // quotient challenge (hash of transcript); fixed to test value
    round_three_out.alpha = example.alpha;

    // Computing the polynomial z(x*w) i.e. z(x) shifted by w where
    // w=srs.omega_roots is the base root of unity and z is
    // z_poly. we do this by multiplying the coefficients of z by w
    round_three_out.z_poly_xomega.resize(round_two_out.z_poly.size());
    std::fill(
        round_three_out.z_poly_xomega.begin(),
        round_three_out.z_poly_xomega.end(),
        Field(0));
    for (size_t i = 0; i < round_two_out.z_poly.size(); ++i) {
        // omega_roots^i
        Field omega_roots_i =
            libff::power(srs.omega_roots[base][1], libff::bigint<1>(i));
        round_three_out.z_poly_xomega[i] =
            round_two_out.z_poly[i] * omega_roots_i;
    }

#ifdef DEBUG
    printf("[%s:%d] this->z_poly_xomega\n", __FILE__, __LINE__);
    print_vector(round_three_out.z_poly_xomega);
    for (size_t i = 0; i < round_two_out.z_poly.size(); ++i) {
        assert(round_three_out.z_poly_xomega[i] == example.z_poly_xomega[i]);
    }
#endif // #ifdef DEBUG

    // start computation of polynomial t(X) in round 3. we break t
    // into 4 parts which we compute separately. each of the 4 parts
    // is multiplied by 1/zh_poly in the paper
    std::vector<polynomial<Field>> t_part(4);

    // --- Computation of t_part[0]

    // a(x)b(x)q_M(x)
    polynomial<Field> abqM;
    libfqfft::_polynomial_multiplication<Field>(
        abqM,
        round_one_out.W_polys_blinded[a],
        round_one_out.W_polys_blinded[b]);
    libfqfft::_polynomial_multiplication<Field>(abqM, abqM, srs.Q_polys[M]);
    // a(x)q_L(x)
    polynomial<Field> aqL;
    libfqfft::_polynomial_multiplication<Field>(
        aqL, round_one_out.W_polys_blinded[a], srs.Q_polys[L]);
    // b(x)q_R(x)
    polynomial<Field> bqR;
    libfqfft::_polynomial_multiplication<Field>(
        bqR, round_one_out.W_polys_blinded[b], srs.Q_polys[R]);
    // c(x)q_O(x)
    polynomial<Field> cqO;
    libfqfft::_polynomial_multiplication<Field>(
        cqO, round_one_out.W_polys_blinded[c], srs.Q_polys[O]);
    // t_part[0](x) = a(x)b(x)q_M(x) + a(x)q_L(x) + b(x)q_R(x) + c(x)q_O(x) +
    // PI(x) + q_C(x)
    polynomial<Field> poly_null{Field(0)};
    libfqfft::_polynomial_addition<Field>(t_part[0], poly_null, abqM);
    libfqfft::_polynomial_addition<Field>(t_part[0], t_part[0], aqL);
    libfqfft::_polynomial_addition<Field>(t_part[0], t_part[0], bqR);
    libfqfft::_polynomial_addition<Field>(t_part[0], t_part[0], cqO);
    libfqfft::_polynomial_addition<Field>(t_part[0], t_part[0], srs.PI_poly);
    libfqfft::_polynomial_addition<Field>(t_part[0], t_part[0], srs.Q_polys[C]);

    // --- Computation of t_part[1]

    // X*beta as polynomial in X
    std::vector<polynomial<Field>> xbeta_poly{
        {Field(0), round_two_out.beta},          // X*beta
        {Field(0), round_two_out.beta * srs.k1}, // X*beta*k1
        {Field(0), round_two_out.beta * srs.k2}  // X*beta*k2
    };
    // represent gamma as polynomial in X, needed for prover Round 3
    polynomial<Field> gamma_poly{round_two_out.gamma}; // gamma
    // represent alpha as polynomial in X, needed for prover Round 3
    polynomial<Field> alpha_poly{round_three_out.alpha}; // alpha

    // a(x) + beta*x + gamma
    polynomial<Field> a_xbeta_gamma;
    libfqfft::_polynomial_addition<Field>(
        a_xbeta_gamma, round_one_out.W_polys_blinded[a], xbeta_poly[base]);
    libfqfft::_polynomial_addition<Field>(
        a_xbeta_gamma, a_xbeta_gamma, gamma_poly);
    // b(x) + beta_k1*x + gamma
    polynomial<Field> b_xbeta_gamma_k1;
    libfqfft::_polynomial_addition<Field>(
        b_xbeta_gamma_k1,
        round_one_out.W_polys_blinded[b],
        xbeta_poly[base_k1]);
    libfqfft::_polynomial_addition<Field>(
        b_xbeta_gamma_k1, b_xbeta_gamma_k1, gamma_poly);
    // c(x) + beta_k1*x + gamma
    polynomial<Field> c_xbeta_gamma_k2;
    libfqfft::_polynomial_addition<Field>(
        c_xbeta_gamma_k2,
        round_one_out.W_polys_blinded[c],
        xbeta_poly[base_k2]);
    libfqfft::_polynomial_addition<Field>(
        c_xbeta_gamma_k2, c_xbeta_gamma_k2, gamma_poly);
    // t_part[1] = (a(x) + beta*x + gamma)*(b(x) + beta_k1*x +
    // gamma)*(c(x) + beta_k1*x + gamma)*z(x)*alpha
    libfqfft::_polynomial_multiplication<Field>(
        t_part[1], a_xbeta_gamma, b_xbeta_gamma_k1);
    libfqfft::_polynomial_multiplication<Field>(
        t_part[1], t_part[1], c_xbeta_gamma_k2);
    libfqfft::_polynomial_multiplication<Field>(
        t_part[1], t_part[1], round_two_out.z_poly);
    libfqfft::_polynomial_multiplication<Field>(
        t_part[1], t_part[1], alpha_poly);

    // --- Computation of t_part[2]

    // represent beta as polynomial in X, needed for prover Round 3
    polynomial<Field> beta_poly{round_two_out.beta};
    // S*beta as polynomial
    // S_sigma1(x)*beta, S_sigma2(x)*beta, S_sigma3(x)*beta
    std::vector<polynomial<Field>> sbeta_poly(num_hgen);
    for (int i = 0; i < num_hgen; ++i) {
        libfqfft::_polynomial_multiplication<Field>(
            sbeta_poly[i], srs.S_polys[i], beta_poly);
    }
    // a(x) + S_sigma1(x)*beta + gamma
    polynomial<Field> a_sbeta_gamma;
    libfqfft::_polynomial_addition<Field>(
        a_sbeta_gamma, round_one_out.W_polys_blinded[a], sbeta_poly[base]);
    libfqfft::_polynomial_addition<Field>(
        a_sbeta_gamma, a_sbeta_gamma, gamma_poly);
    // b(x) + S_sigma2(x)*beta + gamma
    polynomial<Field> b_sbeta_gamma_k1;
    libfqfft::_polynomial_addition<Field>(
        b_sbeta_gamma_k1,
        round_one_out.W_polys_blinded[b],
        sbeta_poly[base_k1]);
    libfqfft::_polynomial_addition<Field>(
        b_sbeta_gamma_k1, b_sbeta_gamma_k1, gamma_poly);
    // b(x) + S_sigma2(x)*beta + gamma
    polynomial<Field> c_sbeta_gamma_k2;
    libfqfft::_polynomial_addition<Field>(
        c_sbeta_gamma_k2,
        round_one_out.W_polys_blinded[c],
        sbeta_poly[base_k2]);
    libfqfft::_polynomial_addition<Field>(
        c_sbeta_gamma_k2, c_sbeta_gamma_k2, gamma_poly);
    // t_part[2] = (a(x) + S_sigma1(x)*beta + gamma)*(b(x) +
    // S_sigma2(x)*beta + gamma)*(b(x) + S_sigma2(x)*beta +
    // gamma)*z(x*srs.omega_roots)*alpha
    libfqfft::_polynomial_multiplication<Field>(
        t_part[2], a_sbeta_gamma, b_sbeta_gamma_k1);
    libfqfft::_polynomial_multiplication<Field>(
        t_part[2], t_part[2], c_sbeta_gamma_k2);
    libfqfft::_polynomial_multiplication<Field>(
        t_part[2], t_part[2], round_three_out.z_poly_xomega);
    libfqfft::_polynomial_multiplication<Field>(
        t_part[2], t_part[2], alpha_poly);
    // -t_part[2]
    polynomial<Field> neg_one_poly = {-Field("1")};
    libfqfft::_polynomial_multiplication<Field>(
        t_part[2], t_part[2], neg_one_poly);

    // --- Computation of t_part[3]

    // z(x) - 1
    polynomial<Field> z_neg_one;
    libfqfft::_polynomial_addition<Field>(
        z_neg_one, round_two_out.z_poly, neg_one_poly);
    // (z(x)-1) * L_1(x)
    libfqfft::_polynomial_multiplication<Field>(
        t_part[3], z_neg_one, srs.L_basis[0]);
    // (z(x)-1) * L_1(x) * alpha
    libfqfft::_polynomial_multiplication<Field>(
        t_part[3], t_part[3], alpha_poly);
    // (z(x)-1) * L_1(x) * alpha * alpha
    libfqfft::_polynomial_multiplication<Field>(
        t_part[3], t_part[3], alpha_poly);

    // --- computation of t(x)

    // t(x) = (t[0] + t[1] + (-t[2]) + t[3]) / zh(x)
    round_three_out.t_poly_long = {Field(0)};
    libfqfft::_polynomial_addition<Field>(
        round_three_out.t_poly_long, round_three_out.t_poly_long, t_part[0]);
#if 1 // DEBUG
    libfqfft::_polynomial_addition<Field>(
        round_three_out.t_poly_long, round_three_out.t_poly_long, t_part[1]);
    libfqfft::_polynomial_addition<Field>(
        round_three_out.t_poly_long, round_three_out.t_poly_long, t_part[2]);
    libfqfft::_polynomial_addition<Field>(
        round_three_out.t_poly_long, round_three_out.t_poly_long, t_part[3]);
#endif // #if 0 // DEBUG
    //    t(x) = t(x) / zh(x): A/B = (Q, R) st. A = (Q * B) + R.
    polynomial<Field> remainder;
    libfqfft::_polynomial_division(
        round_three_out.t_poly_long,
        remainder,
        round_three_out.t_poly_long,
        round_zero_out.zh_poly);
#ifdef DEBUG
    printf("[%s:%d] round_three_out.t_poly_long\n", __FILE__, __LINE__);
    print_vector(round_three_out.t_poly_long);
    assert(round_three_out.t_poly_long == example.t_poly_long);
    printf("[%s:%d] remainder\n", __FILE__, __LINE__);
    print_vector(remainder);
#endif // #ifdef DEBUG
    assert(libfqfft::_is_zero(remainder));

    // break this->t_poly_long into three parts: lo, mid, hi, each of degree
    // 7. note: (srs.num_gates+3) is the length of the CRS =
    // (srs.num_gates+2) powers of G1 + 1 power of G2
    round_three_out.t_poly.resize(num_hgen);
    for (int i = 0; i < num_hgen; ++i) {
        typename std::vector<Field>::iterator begin =
            round_three_out.t_poly_long.begin() + (i * (srs.num_gates + 2));
        typename std::vector<Field>::iterator end =
            round_three_out.t_poly_long.begin() + (i * (srs.num_gates + 2)) +
            (srs.num_gates + 2);
        std::vector<Field> tmp(begin, end);
        round_three_out.t_poly[i] = tmp;
    }
#ifdef DEBUG
    for (int i = 0; i < num_hgen; ++i) {
        printf("[%s:%d] t_poly[%d]\n", __FILE__, __LINE__, i);
        print_vector(round_three_out.t_poly[i]);
        assert(round_three_out.t_poly[i] == example.t_poly[i]);
    }
#endif // #ifdef DEBUG
    // evaluate each part of t_poly in the secret input
    round_three_out.t_poly_at_secret_g1.resize(num_hgen);
    for (int i = 0; i < num_hgen; ++i) {
        round_three_out.t_poly_at_secret_g1[i] =
            plonk_evaluate_poly_at_secret_G1<ppT>(
                srs.secret_powers_g1, round_three_out.t_poly[i]);
    }
    return round_three_out;
}

// Prover Round 4
//
// INPUT
// - W_polys_blinded: blinded witness polynomials (from Round 1)
// - z_poly_xomega: the polynomial z(x*w) i.e. z(x) shifted by w (from
//   Round 3)
// - t_poly_long: the quotient polynomial t(x) (see Round 3, pp28
//   [GWC19]) (from Round 3)
// - srs: structured reference string containing also circuit-specific
//   information
//
// OUTPUT
// - zeta: evaluation challenge -- hash of transcript
// - a_zeta, b_zeta, c_zeta: the blinded witness polynomials a(x),
//   b(x), c(x) (denoted by W_polys_blinded[] output from Round 1)
//   evaluated at x=zeta i.e. a(z), b(z), c(z)
// - S_0_zeta, S_1_zeta: the permutation polynomials S_sigma_1(x),
//   S_sigma_2(x) from the common preprocessed input (see [GWC19],
//   Sect. 8.1) evaluated at x=zeta i.e. S_sigma_1(z), S_sigma_2(z)
// - z_poly_xomega_zeta: the polynomial z(x*w) i.e. z(x) shifted by w
//   (output from Round 3) evaluated at x=zeta i.e. z(zeta*w)
// - t_zeta: the quotient polynomial t(x) output from Round 3, see
//   pp28 [GWC19]) evaluated at x=zeta i.e. t(z). IMPORTANT! the
//   original Plonk proposal [GWC19] does not output this parameter
//   t_zeta. The Python reference implementation does, so we do the
//   same in order to match the test vectors. TODO can remove t_zeta
//   in the future
//
template<typename ppT>
round_four_out_t<ppT> plonk_prover<ppT>::round_four(
    const round_one_out_t<ppT> round_one_out,
    const round_three_out_t<ppT> round_three_out,
    const srs<ppT> srs)
{
    using Field = libff::Fr<ppT>;
    // initialize hard-coded values from example circuit
    plonk_example<ppT> example;

    // output from round 4
    round_four_out_t<ppT> round_four_out;

    // evaluation challenge (hash of transcript); fixed to test value
    round_four_out.zeta = example.zeta;

#ifdef DEBUG
    printf("[%s:%d] zeta\n", __FILE__, __LINE__);
    round_four_out.zeta.print();
#endif // #ifdef DEBUG
    round_four_out.a_zeta = libfqfft::evaluate_polynomial<Field>(
        srs.num_gates + 2,
        round_one_out.W_polys_blinded[a],
        round_four_out.zeta);
    round_four_out.b_zeta = libfqfft::evaluate_polynomial<Field>(
        srs.num_gates + 2,
        round_one_out.W_polys_blinded[b],
        round_four_out.zeta);
    round_four_out.c_zeta = libfqfft::evaluate_polynomial<Field>(
        srs.num_gates + 2,
        round_one_out.W_polys_blinded[c],
        round_four_out.zeta);
    round_four_out.S_0_zeta = libfqfft::evaluate_polynomial<Field>(
        srs.num_gates, srs.S_polys[0], round_four_out.zeta);
    round_four_out.S_1_zeta = libfqfft::evaluate_polynomial<Field>(
        srs.num_gates, srs.S_polys[1], round_four_out.zeta);
    round_four_out.t_zeta = libfqfft::evaluate_polynomial<Field>(
        round_three_out.t_poly_long.size(),
        round_three_out.t_poly_long,
        round_four_out.zeta);
    round_four_out.z_poly_xomega_zeta = libfqfft::evaluate_polynomial<Field>(
        round_three_out.z_poly_xomega.size(),
        round_three_out.z_poly_xomega,
        round_four_out.zeta);
    return round_four_out;
}

// Prover Round 5
//
// INPUT
// - beta, gamma: permutation challenges -- hashes of transcript (from
//   round 2)
// - alpha: quotinet challenge -- hash of transcript (from round 3)
// - zeta: evaluation challenge -- hash of transcript (from round 4)
// - a_zeta, b_zeta, c_zeta: the blinded witness polynomials a(x),
//   b(x), c(x) (denoted by W_polys_blinded[] output from Round 1)
//   evaluated at x=zeta i.e. a(z), b(z), c(z) (from round 4)
// - S_0_zeta, S_1_zeta: the permutation polynomials S_sigma_1(x),
//   S_sigma_2(x) from the common preprocessed input (see [GWC19],
//   Sect. 8.1) evaluated at x=zeta i.e. S_sigma_1(z), S_sigma_2(z)
//   (from round 4)
// - t_zeta: the quotient polynomial t(x) output from Round 3, see
//   pp28 [GWC19]) evaluated at x=zeta i.e. t(z). IMPORTANT! the
//   original Plonk proposal [GWC19] does not output this parameter
//   t_zeta. The Python reference implementation does, so we do the
//   same in order to match the test vectors. TODO can remove t_zeta
//   in the future (from round 4)
// - z_poly_xomega_zeta: the polynomial z(x*w) i.e. z(x) shifted by w
//   (output from Round 3) evaluated at x=zeta i.e. z(zeta*w) (from
//   round 4)
// - W_polys_blinded: blinded witness polynomials (from round 1)
// - t_poly: t(x) divided in three parts t(x) = t_lo(x) + t_mid(x) x^n
//   + t_hi(x) x^{2n} (from round 3)
// - z_poly: blinded accumulator poly z(x) (from round 2)
// - srs: structured reference string containing also circuit-specific
//   information
//
// OUTPUT
// - nu: opening challenge -- hash of transcript (denoted by v in
//   [GWC19])
// - u: multipoint evaluation challenge -- hash of transcript
// - r_zeta: linearisation polynomial r(x) evaluated at x=zeta
//   ie. r(zeta)
// - W_zeta_at_secret: commitment to opening proof polynomial
//   W_zeta(x) at secert input i.e. [W_zeta(secret)]_1
// - W_zeta_omega_at_secret: commitment to opening proof polynomial
//   W_{zeta omega}(x) at secert input i.e. [W_{zeta omega}(secret)]_1
//
template<typename ppT>
round_five_out_t<ppT> plonk_prover<ppT>::round_five(
    const round_zero_out_t<ppT> round_zero_out,
    const round_one_out_t<ppT> round_one_out,
    const round_two_out_t<ppT> round_two_out,
    const round_three_out_t<ppT> round_three_out,
    const round_four_out_t<ppT> round_four_out,
    const srs<ppT> srs)
{
    using Field = libff::Fr<ppT>;
    polynomial<Field> remainder;

    // initialize hard-coded values from example circuit
    plonk_example<ppT> example;

    // output from round 5
    round_five_out_t<ppT> round_five_out;

    // openinig challenge (hash of transcript); fixed to test value
    // (note: denoted v in [GWZ19])
    round_five_out.nu = example.nu;

    // compute linerisation polynomial r in five parts
    std::vector<polynomial<Field>> r_part(5);

    // --- Computation of r_part[0]

    // represent values as constant term polynomials in orderto use
    // the functions in the libfqfft library on polynomials
    polynomial<Field> a_zeta_poly{round_four_out.a_zeta};
    polynomial<Field> b_zeta_poly{round_four_out.b_zeta};
    polynomial<Field> c_zeta_poly{round_four_out.c_zeta};
    // a(z)b(z)q_M(x)
    polynomial<Field> abqM_zeta;
    libfqfft::_polynomial_multiplication<Field>(
        abqM_zeta, srs.Q_polys[M], a_zeta_poly);
    libfqfft::_polynomial_multiplication<Field>(
        abqM_zeta, abqM_zeta, b_zeta_poly);
    // a(z)q_L(x)
    polynomial<Field> aqL_zeta;
    libfqfft::_polynomial_multiplication<Field>(
        aqL_zeta, srs.Q_polys[L], a_zeta_poly);
    // b(z)q_R(x)
    polynomial<Field> bqR_zeta;
    libfqfft::_polynomial_multiplication<Field>(
        bqR_zeta, srs.Q_polys[R], b_zeta_poly);
    // c(z)q_O(x)
    polynomial<Field> cqO_zeta;
    libfqfft::_polynomial_multiplication<Field>(
        cqO_zeta, srs.Q_polys[O], c_zeta_poly);
    // a(z)b(z)q_M(x) + a(z)q_L(x) + b(z)q_R(x) + c(z)q_O(x) + q_C(x)
    libfqfft::_polynomial_addition<Field>(
        r_part[0], round_zero_out.null_poly, abqM_zeta);
    libfqfft::_polynomial_addition<Field>(r_part[0], r_part[0], aqL_zeta);
    libfqfft::_polynomial_addition<Field>(r_part[0], r_part[0], bqR_zeta);
    libfqfft::_polynomial_addition<Field>(r_part[0], r_part[0], cqO_zeta);
    libfqfft::_polynomial_addition<Field>(r_part[0], r_part[0], srs.Q_polys[C]);

    // --- Computation of r_part[1]

    polynomial<Field> r1_const_poly{
        (round_four_out.a_zeta + (round_two_out.beta * round_four_out.zeta) +
         round_two_out.gamma) *
        (round_four_out.b_zeta +
         (round_two_out.beta * srs.k1 * round_four_out.zeta) +
         round_two_out.gamma) *
        (round_four_out.c_zeta +
         (round_two_out.beta * srs.k2 * round_four_out.zeta) +
         round_two_out.gamma) *
        round_three_out.alpha};
    libfqfft::_polynomial_multiplication<Field>(
        r_part[1], r1_const_poly, round_two_out.z_poly);

    // --- Computation of r_part[2]

    polynomial<Field> r2_const_poly{
        (round_four_out.a_zeta +
         (round_two_out.beta * round_four_out.S_0_zeta) + round_two_out.gamma) *
        (round_four_out.b_zeta +
         (round_two_out.beta * round_four_out.S_1_zeta) + round_two_out.gamma) *
        (round_three_out.alpha * round_two_out.beta *
         round_four_out.z_poly_xomega_zeta)};
    libfqfft::_polynomial_multiplication<Field>(
        r_part[2], r2_const_poly, srs.S_polys[2]);
    // -r_part[2]
    libfqfft::_polynomial_multiplication<Field>(
        r_part[2], r_part[2], round_zero_out.neg_one_poly);

    // --- Computation of r_part[3]

    //     r3 = accumulator_poly_ext3 * eval_poly(L_1, [zeta])[0] * alpha ** 2
    polynomial<Field> L_0_zeta_poly{libfqfft::evaluate_polynomial<Field>(
        srs.L_basis[0].size(), srs.L_basis[0], round_four_out.zeta)};
    polynomial<Field> alpha_power2_poly{
        libff::power(round_three_out.alpha, libff::bigint<1>(2))};
    libfqfft::_polynomial_multiplication<Field>(
        r_part[3], round_two_out.z_poly, L_0_zeta_poly);
    libfqfft::_polynomial_multiplication<Field>(
        r_part[3], r_part[3], alpha_power2_poly);

    // --- Computation of r_poly = (r0+r1-r2+r3)

    //
    // Note: here the reference Python implementation differs from the
    // paper where:
    //
    // r(x) = r(x) - zh(zeta) (t_lo(x) + zeta^n t_mid(x) + zeta^2n t_hi(x))
    //
    // In the reference implementation, the missing term is added in
    // the computation of the W_zeta(x) polynomial
    //
    // linearisation polynomial r(x)
    polynomial<Field> r_poly;
    libfqfft::_polynomial_addition<Field>(
        r_poly, round_zero_out.null_poly, r_part[0]);
    libfqfft::_polynomial_addition<Field>(r_poly, r_poly, r_part[1]);
    libfqfft::_polynomial_addition<Field>(r_poly, r_poly, r_part[2]);
    libfqfft::_polynomial_addition<Field>(r_poly, r_poly, r_part[3]);

#ifdef DEBUG
    //    printf("abqM_zeta\n");
    //    print_vector(abqM_zeta);
    printf("[%s:%d] r_part[0]\n", __FILE__, __LINE__);
    print_vector(r_part[0]);
    printf("[%s:%d] r_part[1]\n", __FILE__, __LINE__);
    print_vector(r_part[1]);
    printf("[%s:%d] r_part[2]\n", __FILE__, __LINE__);
    print_vector(r_part[2]);
    printf("[%s:%d] r_part[3]\n", __FILE__, __LINE__);
    print_vector(r_part[3]);
    printf("[%s:%d] r_poly\n", __FILE__, __LINE__);
    print_vector(r_poly);
    assert(r_poly == example.r_poly);
#endif // #ifdef DEBUG

    // Evaluate the r-polynomial at zeta. Note: in the reference
    // implementation, r_zeta is added to the pi-SNARK proof. In the
    // paper this is omitted, which makes the proof shorter at the
    // epxense of a slightly heavier computation on the verifier's
    // side
    round_five_out.r_zeta = libfqfft::evaluate_polynomial<Field>(
        r_poly.size(), r_poly, round_four_out.zeta);
#ifdef DEBUG
    printf("r_zeta ");
    round_five_out.r_zeta.print();
    assert(round_five_out.r_zeta == example.r_zeta);
#endif // #ifdef DEBUG

    // W_zeta polynomial is of degree 6 in the random element nu and
    // hence has 7 terms
    std::vector<polynomial<Field>> W_zeta_part(7);

    // --- compute W_zeta_part[0]

    // t_lo(x)
    polynomial<Field> t_lo{round_three_out.t_poly[lo]};
    // t_mid(x) * zeta^(n+2)
    polynomial<Field> t_mid_zeta_n;
    polynomial<Field> zeta_powern_poly{
        libff::power(round_four_out.zeta, libff::bigint<1>(srs.num_gates + 2))};
    libfqfft::_polynomial_multiplication<Field>(
        t_mid_zeta_n, round_three_out.t_poly[mid], zeta_powern_poly);
    // t_hi(x) * zeta^(2(n+1))
    polynomial<Field> t_hi_zeta_2n;
    polynomial<Field> zeta_power2n_poly{libff::power(
        round_four_out.zeta, libff::bigint<1>(2 * (srs.num_gates + 2)))};
    libfqfft::_polynomial_multiplication<Field>(
        t_hi_zeta_2n, round_three_out.t_poly[hi], zeta_power2n_poly);
    // -t_zeta as constant term polynomial
    polynomial<Field> t_zeta_poly{-round_four_out.t_zeta};
    // t_lo(x) + (t_mid(x) * zeta^n) + (t_hi(x) * zeta^2n) + t_zeta_poly
    libfqfft::_polynomial_addition<Field>(
        W_zeta_part[0], round_zero_out.null_poly, t_lo);
    libfqfft::_polynomial_addition<Field>(
        W_zeta_part[0], W_zeta_part[0], t_mid_zeta_n);
    libfqfft::_polynomial_addition<Field>(
        W_zeta_part[0], W_zeta_part[0], t_hi_zeta_2n);
    libfqfft::_polynomial_addition<Field>(
        W_zeta_part[0], W_zeta_part[0], t_zeta_poly);

    // --- compute W_zeta_part[1]

    // -r_zeta as constant term polynomial
    polynomial<Field> r_zeta_poly{-round_five_out.r_zeta};
    // r(x) - r_zeta
    polynomial<Field> r_sub_rzeta;
    libfqfft::_polynomial_addition<Field>(r_sub_rzeta, r_poly, r_zeta_poly);
    // (r(x) - r_zeta) * nu
    polynomial<Field> nu_poly{round_five_out.nu};
    libfqfft::_polynomial_multiplication<Field>(
        W_zeta_part[1], r_sub_rzeta, nu_poly);

    // --- compute W_zeta_part[2]

    // -a_zeta as constant term polynomial
    polynomial<Field> a_zeta_poly_neg;
    libfqfft::_polynomial_multiplication<Field>(
        a_zeta_poly_neg, a_zeta_poly, round_zero_out.neg_one_poly);
    // a(x) - a_zeta
    polynomial<Field> a_sub_azeta;
    libfqfft::_polynomial_addition<Field>(
        a_sub_azeta, round_one_out.W_polys_blinded[a], a_zeta_poly_neg);
    // (a(x) - a_zeta) * nu^2
    Field nu2 = libff::power(round_five_out.nu, libff::bigint<1>(2));
    polynomial<Field> nu2_poly{nu2};
    libfqfft::_polynomial_multiplication<Field>(
        W_zeta_part[2], a_sub_azeta, nu2_poly);

    // -b_zeta as constant term polynomial
    polynomial<Field> b_zeta_poly_neg;
    libfqfft::_polynomial_multiplication<Field>(
        b_zeta_poly_neg, b_zeta_poly, round_zero_out.neg_one_poly);
    // (b(x) - b_zeta)
    polynomial<Field> b_sub_bzeta;
    libfqfft::_polynomial_addition<Field>(
        b_sub_bzeta, round_one_out.W_polys_blinded[b], b_zeta_poly_neg);
    // (b(x) - b_zeta) * nu^3
    Field nu3 = libff::power(round_five_out.nu, libff::bigint<1>(3));
    polynomial<Field> nu3_poly{nu3};
    libfqfft::_polynomial_multiplication<Field>(
        W_zeta_part[3], b_sub_bzeta, nu3_poly);

    // -c_zeta as constant term polynomial
    polynomial<Field> c_zeta_poly_neg;
    libfqfft::_polynomial_multiplication<Field>(
        c_zeta_poly_neg, c_zeta_poly, round_zero_out.neg_one_poly);
    // (c(x) - c_zeta)
    polynomial<Field> c_sub_czeta;
    libfqfft::_polynomial_addition<Field>(
        c_sub_czeta, round_one_out.W_polys_blinded[c], c_zeta_poly_neg);
    // (c(x) - c_zeta) * nu^4
    Field nu4 = libff::power(round_five_out.nu, libff::bigint<1>(4));
    polynomial<Field> nu4_poly{nu4};
    libfqfft::_polynomial_multiplication<Field>(
        W_zeta_part[4], c_sub_czeta, nu4_poly);

    // -S_0_zeta as constant term polynomial
    polynomial<Field> S_0_zeta_poly_neg{-round_four_out.S_0_zeta};
    //    libfqfft::_polynomial_multiplication<Field>(S_0_zeta_poly_neg,
    //    S_0_zeta_poly, round_zero_out.neg_one_poly);
    // (S0(x) - S_0_zeta)
    polynomial<Field> S0_sub_szeta;
    libfqfft::_polynomial_addition<Field>(
        S0_sub_szeta, srs.S_polys[0], S_0_zeta_poly_neg);
    // (S0(x) - S_0_zeta) * nu^5
    Field nu5 = libff::power(round_five_out.nu, libff::bigint<1>(5));
    polynomial<Field> nu5_poly{nu5};
    libfqfft::_polynomial_multiplication<Field>(
        W_zeta_part[5], S0_sub_szeta, nu5_poly);

    // -S_1_zeta as constant term polynomial
    polynomial<Field> S_1_zeta_poly_neg{-round_four_out.S_1_zeta};
    //    libfqfft::_polynomial_multiplication<Field>(S_1_zeta_poly_neg,
    //    S_1_zeta_poly, neg_one_poly);
    // (S1(x) - S_1_zeta)
    polynomial<Field> S1_sub_szeta;
    libfqfft::_polynomial_addition<Field>(
        S1_sub_szeta, srs.S_polys[1], S_1_zeta_poly_neg);
    // (S1(x) - S_1_zeta) * nu^6
    Field nu6 = libff::power(round_five_out.nu, libff::bigint<1>(6));
    polynomial<Field> nu6_poly{nu6};
    libfqfft::_polynomial_multiplication<Field>(
        W_zeta_part[6], S1_sub_szeta, nu6_poly);

    // compute full zeta polynomial W_zeta = \sum W_zeta_part[i]
    int nzeta = 7;
    polynomial<Field> W_zeta(round_zero_out.null_poly);
    for (int i = 0; i < nzeta; ++i) {
        libfqfft::_polynomial_addition<Field>(W_zeta, W_zeta, W_zeta_part[i]);
    }

    // compute 1/(X-zeta) * W_zeta
    polynomial<Field> x_sub_zeta_poly{-round_four_out.zeta, Field(1)};
    libfqfft::_polynomial_division(W_zeta, remainder, W_zeta, x_sub_zeta_poly);
#ifdef DEBUG
    printf("W_zeta\n");
    print_vector(W_zeta);
#endif // #ifdef DEBUG
    assert(libfqfft::_is_zero(remainder));

    //    polynomial<Field> z_poly;

    // Compute opening proof:
    // W_zeta_omega = z(X) - z(zeta*srs.omega_roots) / X -
    // (zeta*srs.omega_roots)
    polynomial<Field> W_zeta_omega{round_zero_out.null_poly};

    // -z(zeta*srs.omega_roots)
    polynomial<Field> z_poly_xomega_zeta_neg{
        -round_four_out.z_poly_xomega_zeta};
    // z(X) - z(zeta*srs.omega_roots)
    libfqfft::_polynomial_addition<Field>(
        W_zeta_omega, round_two_out.z_poly, z_poly_xomega_zeta_neg);
    // -zeta*srs.omega_roots; srs.omega_roots[base][1] =
    // srs.omega_roots_base
    polynomial<Field> x_sub_zeta_omega_roots{
        -(round_four_out.zeta * srs.omega_roots[base][1]), Field(1)};

    // z(X) - z(zeta*srs.omega_roots) / X -
    // (zeta*srs.omega_roots)
    libfqfft::_polynomial_division(
        W_zeta_omega, remainder, W_zeta_omega, x_sub_zeta_omega_roots);
    assert(libfqfft::_is_zero(remainder));

#ifdef DEBUG
    printf("W_zeta_part[0]\n");
    print_vector(W_zeta_part[0]);
    printf("W_zeta_part[1]\n");
    print_vector(W_zeta_part[1]);
    printf("W_zeta_part[2]\n");
    print_vector(W_zeta_part[2]);
    printf("W_zeta_part[3]\n");
    print_vector(W_zeta_part[3]);
    printf("W_zeta_part[4]\n");
    print_vector(W_zeta_part[4]);
    printf("W_zeta_part[5]\n");
    print_vector(W_zeta_part[5]);
    printf("W_zeta_part[6]\n");
    print_vector(W_zeta_part[6]);
    printf("W_zeta\n");
    print_vector(W_zeta);
    printf("W_zeta_omega\n");
    print_vector(W_zeta_omega);
#endif // #ifdef DEBUG

    assert(W_zeta == example.W_zeta);
    assert(W_zeta_omega == example.W_zeta_omega);

    // Evaluate polynomials W_zeta and W_zeta_omega at the seceret
    // input
    round_five_out.W_zeta_at_secret =
        plonk_evaluate_poly_at_secret_G1<ppT>(srs.secret_powers_g1, W_zeta);
    round_five_out.W_zeta_omega_at_secret =
        plonk_evaluate_poly_at_secret_G1<ppT>(
            srs.secret_powers_g1, W_zeta_omega);

    // Hashes of transcript (Fiat-Shamir heuristic) -- fixed to match
    // the test vectors
    round_five_out.u = example.u;

    return round_five_out;
}

// Prover compute SNARK proof
//
// Pi ([a]_1, [b]_1, [c]_1, [z]_1,
//     [t_lo]_1, [t_mi]_1, [t_hi]_1,
//     \bar{a}, \bar{b}, \bar{c},
//     \bar{S_sigma1}, \bar{S_sigma2}, \bar{z_w},
//     [W_zeta]_1, [W_{zeta omega}]_1
//     r_zeta (*))
//
// (*) Note: in the reference Python implementation, r_zeta (the
// evaluation of the linearlization polynomial r(X) at zeta from
// Prover round 5) is added to the pi-SNARK proof. In the paper
// this is omitted, which seems to make the proof shorter by 1
// element at the epxense of a slightly heavier computation on the
// verifier's side. Here we follow the reference implementation to
// make sure we match the test values. TODO: once all test vectors
// are verified, we may remove r_zeta from the proof to be fully
// compliant with the paper.
//
// Mapping code-to-paper quantities
//
// - W_polys_blinded_at_secret_g1[a, b, c]: [a]_1, [b]_1, [c]_1 (from Round 1)
// - z_poly_at_secret_g1: [z]_1 (from Round 2)
// - t_poly_at_secret_g1[lo, mi, hi]: [t_lo]_1, [t_mi]_1, [t_hi]_1 (from Round
// 3)
// - a_zeta, b_zeta, c_zeta, S_0_zeta, S_1_zeta, z_poly_xomega_zeta: \bar{a},
// \bar{b}, \bar{c}, \bar{S_sigma1}, \bar{S_sigma2}, \bar{z_w} (from Round 4)
// - W_zeta_at_secret, W_zeta_omega_at_secret: [W_zeta]_1, [W_{zeta omega}]_1
// (from Round 5)
//
// INPUT
// - srs: structured reference string containing also circuit-specific
//   information
//
// OUTPUT
// - proof: SNARK proof Pi (see above)
//
template<typename ppT>
plonk_proof<ppT> plonk_prover<ppT>::compute_proof(const srs<ppT> srs)
{
    using Field = libff::Fr<ppT>;

    // initialize hard-coded values from example circuit
    plonk_example<ppT> example;
    std::vector<Field> witness = example.witness;
#ifdef DEBUG
    assert(srs.k1 == example.k1);
    assert(srs.k2 == example.k2);
#endif // #ifdef DEBUG

    // Prover Round 0 (initialization)
#if 1 // prover round 0
    printf("[%s:%d] Prover Round 0...\n", __FILE__, __LINE__);
    round_zero_out_t<ppT> round_zero_out = plonk_prover::round_zero(srs);
#ifdef DEBUG
    printf("[%s:%d] Vanishing polynomial\n", __FILE__, __LINE__);
    print_vector(round_zero_out.zh_poly);
    assert(round_zero_out.zh_poly == example.zh_poly);
#endif // #ifdef DEBUG
#endif // #if 1 // prover round 0

    // Prover Round 1
#if 1 // prover round 1
    printf("[%s:%d] Prover Round 1...\n", __FILE__, __LINE__);
    round_one_out_t<ppT> round_one_out =
        plonk_prover::round_one(round_zero_out, witness, srs);
    // Prover Round 1 output check against test vectors
#ifdef DEBUG
    for (int i = 0; i < (int)NUM_HGEN; ++i) {
        printf("[%s:%d] W_polys_blinded[%d]\n", __FILE__, __LINE__, i);
        print_vector(round_one_out.W_polys_blinded[i]);
        assert(round_one_out.W_polys_blinded[i] == example.W_polys_blinded[i]);
    }
    printf("[%s:%d] Output from Round 1\n", __FILE__, __LINE__);
    for (int i = 0; i < (int)NUM_HGEN; ++i) {
        printf("W_polys_at_secret_g1[%d]\n", i);
        round_one_out.W_polys_blinded_at_secret_g1[i].print();
        libff::G1<ppT> W_polys_blinded_at_secret_g1_i(
            round_one_out.W_polys_blinded_at_secret_g1[i]);
        W_polys_blinded_at_secret_g1_i.to_affine_coordinates();
        assert(
            W_polys_blinded_at_secret_g1_i.X ==
            example.W_polys_blinded_at_secret_g1[i][0]);
        assert(
            W_polys_blinded_at_secret_g1_i.Y ==
            example.W_polys_blinded_at_secret_g1[i][1]);
    }
#endif // #ifdef DEBUG
#endif // #if 1 // prover round 1

    printf("[%s:%d] Prover Round 2...\n", __FILE__, __LINE__);
#if 1 // prover round 2
    round_two_out_t<ppT> round_two_out =
        plonk_prover::round_two(round_zero_out, round_one_out, witness, srs);
    // Prover Round 2 output check against test vectors
#ifdef DEBUG
    printf("[%s:%d] z_poly\n", __FILE__, __LINE__);
    print_vector(round_two_out.z_poly);
    assert(round_two_out.z_poly == example.z_poly);
    printf("[%s:%d] Output from Round 2\n", __FILE__, __LINE__);
    printf("[%s:%d] z_poly_at_secret_g1\n", __FILE__, __LINE__);
    round_two_out.z_poly_at_secret_g1.print();
    libff::G1<ppT> z_poly_at_secret_g1_aff(round_two_out.z_poly_at_secret_g1);
    z_poly_at_secret_g1_aff.to_affine_coordinates();
    assert(z_poly_at_secret_g1_aff.X == example.z_poly_at_secret_g1[0]);
    assert(z_poly_at_secret_g1_aff.Y == example.z_poly_at_secret_g1[1]);
#endif // #ifdef DEBUG
#endif // #if 1 // prover round 2

    printf("[%s:%d] Prover Round 3...\n", __FILE__, __LINE__);
#if 1 // prover round 3
    round_three_out_t<ppT> round_three_out = plonk_prover::round_three(
        round_zero_out, round_one_out, round_two_out, srs);
    // Prover Round 3 output check against test vectors
#ifdef DEBUG
    printf("[%s:%d] Output from Round 3\n", __FILE__, __LINE__);
    for (int i = 0; i < (int)NUM_HGEN; ++i) {
        printf("[%s:%d] t_poly_at_secret_g1[%d]\n", __FILE__, __LINE__, i);
        round_three_out.t_poly_at_secret_g1[i].print();
        libff::G1<ppT> t_poly_at_secret_g1_i(
            round_three_out.t_poly_at_secret_g1[i]);
        t_poly_at_secret_g1_i.to_affine_coordinates();
        assert(t_poly_at_secret_g1_i.X == example.t_poly_at_secret_g1[i][0]);
        assert(t_poly_at_secret_g1_i.Y == example.t_poly_at_secret_g1[i][1]);
    }
#endif // #ifdef DEBUG
#endif // #if 1 // prover round 3

    printf("[%s:%d] Prover Round 4...\n", __FILE__, __LINE__);
#if 1 // prover round 4
    round_four_out_t<ppT> round_four_out =
        plonk_prover::round_four(round_one_out, round_three_out, srs);
    // Prover Round 4 output check against test vectors
#ifdef DEBUG
    printf("[%s:%d] Output from Round 4\n", __FILE__, __LINE__);
    printf("a_zeta ");
    round_four_out.a_zeta.print();
    assert(round_four_out.a_zeta == example.a_zeta);
    printf("b_zeta ");
    round_four_out.b_zeta.print();
    assert(round_four_out.b_zeta == example.b_zeta);
    printf("c_zeta ");
    round_four_out.c_zeta.print();
    assert(round_four_out.c_zeta == example.c_zeta);
    printf("S_0_zeta ");
    round_four_out.S_0_zeta.print();
    assert(round_four_out.S_0_zeta == example.S_0_zeta);
    printf("S_1_zeta ");
    round_four_out.S_1_zeta.print();
    assert(round_four_out.S_1_zeta == example.S_1_zeta);
    printf("t_zeta ");
    round_four_out.t_zeta.print();
    assert(round_four_out.t_zeta == example.t_zeta);
    printf("z_poly_xomega_zeta ");
    round_four_out.z_poly_xomega_zeta.print();
    assert(round_four_out.z_poly_xomega_zeta == example.z_poly_xomega_zeta);
#endif // #ifdef DEBUG
#endif // #if 1 // prover round 4

    printf("[%s:%d] Prover Round 5...\n", __FILE__, __LINE__);
#if 1 // prover round 5
    round_five_out_t<ppT> round_five_out = plonk_prover::round_five(
        round_zero_out,
        round_one_out,
        round_two_out,
        round_three_out,
        round_four_out,
        srs);
#endif // #if 1 // prover round 5
#ifdef DEBUG
    printf("[%s:%d] Outputs from Prover round 5\n", __FILE__, __LINE__);
    printf("r_zeta ");
    round_five_out.r_zeta.print();
    assert(round_five_out.r_zeta == example.r_zeta);
    printf("[%s:%d] W_zeta_at_secret \n", __FILE__, __LINE__);
    round_five_out.W_zeta_at_secret.print();
    libff::G1<ppT> W_zeta_at_secret_aff(round_five_out.W_zeta_at_secret);
    W_zeta_at_secret_aff.to_affine_coordinates();
    assert(W_zeta_at_secret_aff.X == example.W_zeta_at_secret[0]);
    assert(W_zeta_at_secret_aff.Y == example.W_zeta_at_secret[1]);
    printf("[%s:%d] W_zeta_omega_at_secret \n", __FILE__, __LINE__);
    round_five_out.W_zeta_omega_at_secret.print();
    libff::G1<ppT> W_zeta_omega_at_secret_aff(
        round_five_out.W_zeta_omega_at_secret);
    W_zeta_omega_at_secret_aff.to_affine_coordinates();
    assert(W_zeta_omega_at_secret_aff.X == example.W_zeta_omega_at_secret[0]);
    assert(W_zeta_omega_at_secret_aff.Y == example.W_zeta_omega_at_secret[1]);
#endif // #ifdef DEBUG

    // construct proof
    plonk_proof<ppT> proof(
        round_one_out.W_polys_blinded_at_secret_g1,
        round_two_out.z_poly_at_secret_g1,
        round_three_out.t_poly_at_secret_g1,
        round_four_out.a_zeta,
        round_four_out.b_zeta,
        round_four_out.c_zeta,
        round_four_out.S_0_zeta,
        round_four_out.S_1_zeta,
        round_four_out.z_poly_xomega_zeta,
        round_five_out.W_zeta_at_secret,
        round_five_out.W_zeta_omega_at_secret,
        round_five_out.r_zeta);

    // return proof
    return proof;
}

} // namespace libsnark

#endif // PLONK_PPZKSNARK_PROVER_TCC_
