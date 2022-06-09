/** @file
*****************************************************************************
Implementation of Common Preprocessed Input interfaces for a ppzkSNARK
for Plonk.

See common_input.hpp .

*****************************************************************************
* @author     This file is part of libsnark, developed by SCIPR Lab
*             and contributors (see AUTHORS).
* @copyright  MIT license (see LICENSE file)
*****************************************************************************/

#ifndef PLONK_PPZKSNARK_COMMON_INPUT_TCC_
#define PLONK_PPZKSNARK_COMMON_INPUT_TCC_

namespace libsnark
{
template<typename ppT>
void common_preprocessed_input<ppT>::setup_from_example(
    plonk_example<ppT> example)
{
    using Field = libff::Fr<ppT>;

    // public input (PI)
    Field PI_value = example.public_input;
    // index of the row of the PI in the non-transposed gates_matrix
    int PI_index = example.public_input_index;
    // Transposed gates matrix: each row is a q-vector
    std::vector<std::vector<Field>> gates_matrix_transpose =
        example.gates_matrix_transpose;
    // wire permutation
    std::vector<size_t> wire_permutation = example.wire_permutation;
    // Generate domains on which to evaluate the witness
    // polynomials. k1,k2 can be random, but we fix them for debug to
    // match against the test vector values
    Field k1 = example.k1;
    Field k2 = example.k2;
    this->k1 = example.k1;
    this->k2 = example.k2;
#ifdef DEBUG
    printf("[%s:%d] k1 ", __FILE__, __LINE__);
    k1.print();
    printf("[%s:%d] k2 ", __FILE__, __LINE__);
    k2.print();
#endif // #ifdef DEBUG

    this->num_gates = example.num_gates;
    // ensure that num_gates is not 0
    assert(this->num_gates);
#ifdef DEBUG
    // ensure num_gates is power of 2
    bool b_is_power2 = ((this->num_gates & (this->num_gates - 1)) == 0);
    assert(b_is_power2);
#endif // #ifdef DEBUG

    this->num_qpolys = example.num_qpolys;

    // We represent the constraints q_L, q_R, q_O, q_M, q_C and the
    // witness w_L, w_R, w_O as polynomials in the roots of unity
    // e.g. f_{q_L}(omega_i) = q_L[i], 0\le{i}<8
    // compute Lagrange basis
    this->L_basis.resize(this->num_gates, polynomial<Field>(this->num_gates));
    std::shared_ptr<libfqfft::evaluation_domain<Field>> domain =
        libfqfft::get_evaluation_domain<Field>(this->num_gates);
    plonk_compute_lagrange_basis<Field>(this->num_gates, this->L_basis);

    // compute public input (PI) polynomial
    std::vector<Field> PI_points(this->num_gates, Field(0));
    PI_points[PI_index] = Field(-PI_value);
    plonk_compute_public_input_polynomial(
        PI_points, this->PI_poly, this->L_basis);
#ifdef DEBUG
    printf("[%s:%d] this->PI_poly\n", __FILE__, __LINE__);
    print_vector(this->PI_poly);
#endif // #ifdef DEBUG

    // compute the selector polynomials (q-polynomials) from the
    // transposed gates matrix over the Lagrange basis q_poly = \sum_i
    // q[i] * L[i] where q[i] is a coefficient (a scalar Field
    // element) and L[i] is a polynomial with Field coefficients
    this->Q_polys.resize(this->num_qpolys, polynomial<Field>(this->num_gates));
    plonk_compute_selector_polynomials<Field>(
        gates_matrix_transpose, this->Q_polys, this->L_basis);
#ifdef DEBUG
    for (int i = 0; i < (int)this->num_qpolys; ++i) {
        printf("\n[%s:%d] this->Q_polys[%2d]\n", __FILE__, __LINE__, i);
        print_vector(this->Q_polys[i]);
    }
#endif // #ifdef DEBUG

    // number of generators for H, Hk1, Hk2
    int num_hgen = NUM_HGEN;
    // omega[0] are the n roots of unity, omega[1] are omega[0]*k1,
    // omega[2] are omega[0]*k2
    //    std::vector<std::vector<Field>> omega;
    this->omega_roots.resize(num_hgen, std::vector<Field>(this->num_gates));
    plonk_compute_roots_of_unity_omega(k1, k2, this->omega_roots);
    // H_gen contains the generators of H, k1 H and k2 H in one place
    // ie. this->omega_roots, this->omega_roots_k1 and this->omega_roots_k2
    plonk_roots_of_unity_omega_to_subgroup_H(this->omega_roots, this->H_gen);
#ifdef DEBUG
    printf("[%s:%d] this->H_gen\n", __FILE__, __LINE__);
    print_vector(this->H_gen);
    for (int i = 0; i < (int)this->H_gen.size(); ++i) {
        assert(this->H_gen[i] == example.H_gen[i]);
    }
#endif // #ifdef DEBUG

    // permute this->H_gen according to the wire permutation
    this->H_gen_permute.resize(num_hgen * this->num_gates, Field(0));
    plonk_permute_subgroup_H<Field>(
        this->H_gen_permute, this->H_gen, wire_permutation);
#ifdef DEBUG
    printf("[%s:%d] this->H_gen_permute\n", __FILE__, __LINE__);
    print_vector(this->H_gen_permute);
    for (size_t i = 0; i < this->H_gen_permute.size(); ++i) {
        assert(this->H_gen_permute[i] == example.H_gen_permute[i]);
    }
#endif // #ifdef DEBUG

    // compute the permutation polynomials S_sigma_1, S_sigma_2,
    // S_sigma_3 (see [GWC19], Sect. 8.1) (our indexing starts from 0)
    this->S_polys.resize(num_hgen, polynomial<Field>(this->num_gates));
    plonk_compute_permutation_polynomials<Field>(
        this->S_polys, this->H_gen_permute, this->num_gates);
#ifdef DEBUG
    for (int i = 0; i < num_hgen; ++i) {
        printf("[%s:%d] this->S_polys[%d]\n", __FILE__, __LINE__, i);
        print_vector(this->S_polys[i]);
    }
#endif // #ifdef DEBUG
}
} // namespace libsnark

#endif // PLONK_PPZKSNARK_COMMON_INPUT_TCC_
