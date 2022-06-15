/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_CIRCUIT_TCC_
#define LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_CIRCUIT_TCC_

/// Implementation of Common Preprocessed Input interfaces for a
/// ppzkSNARK for Plonk. See circuit.hpp .

namespace libsnark
{

/// Compute or fill-in ciruit specific data from example.
template<typename ppT>
circuit_t<ppT> plonk_curcuit_description_from_example(
    const plonk_example<ppT> example)
{
    using Field = libff::Fr<ppT>;

    circuit_t<ppT> circuit;

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
    circuit.k1 = example.k1;
    circuit.k2 = example.k2;
#ifdef DEBUG
    printf("[%s:%d] k1 ", __FILE__, __LINE__);
    circuit.k1.print();
    printf("[%s:%d] k2 ", __FILE__, __LINE__);
    circuit.k2.print();
#endif // #ifdef DEBUG

    circuit.num_gates = example.num_gates;
    // TODO: throw exception
#ifdef DEBUG
    // ensure that num_gates is not 0
    assert(circuit.num_gates);
    // ensure num_gates is power of 2
    bool b_is_power2 = ((circuit.num_gates & (circuit.num_gates - 1)) == 0);
    assert(b_is_power2);
#endif // #ifdef DEBUG

    circuit.num_qpolys = example.num_qpolys;

    // We represent the constraints q_L, q_R, q_O, q_M, q_C and the
    // witness w_L, w_R, w_O as polynomials in the roots of unity
    // e.g. f_{q_L}(omega_i) = q_L[i], 0\le{i}<8
    // compute Lagrange basis
    circuit.L_basis.resize(
        circuit.num_gates, polynomial<Field>(circuit.num_gates));
    std::shared_ptr<libfqfft::evaluation_domain<Field>> domain =
        libfqfft::get_evaluation_domain<Field>(circuit.num_gates);
    plonk_compute_lagrange_basis<Field>(circuit.num_gates, circuit.L_basis);

    // compute public input (PI) polynomial
    std::vector<Field> PI_points(circuit.num_gates, Field(0));
    PI_points[PI_index] = Field(-PI_value);
    plonk_compute_public_input_polynomial(
        PI_points, circuit.L_basis, circuit.PI_poly);

    // compute the selector polynomials (q-polynomials) from the
    // transposed gates matrix over the Lagrange basis q_poly = \sum_i
    // q[i] * L[i] where q[i] is a coefficient (a scalar Field
    // element) and L[i] is a polynomial with Field coefficients
    circuit.Q_polys.resize(
        circuit.num_qpolys, polynomial<Field>(circuit.num_gates));
    plonk_compute_selector_polynomials<Field>(
        gates_matrix_transpose, circuit.L_basis, circuit.Q_polys);

    // number of generators for H, Hk1, Hk2
    int num_hgen = NUM_HGEN;
    // omega[0] are the n roots of unity, omega[1] are omega[0]*k1,
    // omega[2] are omega[0]*k2
    //    std::vector<std::vector<Field>> omega;
    circuit.omega_roots.resize(num_hgen, std::vector<Field>(circuit.num_gates));
    plonk_compute_roots_of_unity_omega(
        circuit.k1, circuit.k2, circuit.omega_roots);
    // H_gen contains the generators of H, k1 H and k2 H in one place
    // ie. circuit.omega_roots, circuit.omega_roots_k1 and
    // circuit.omega_roots_k2
    plonk_roots_of_unity_omega_to_subgroup_H(
        circuit.omega_roots, circuit.H_gen);
    // TODO: write unit test for plonk_roots_of_unity_omega_to_subgroup_H
#ifdef DEBUG
    printf("[%s:%d] circuit.H_gen\n", __FILE__, __LINE__);
    print_vector(circuit.H_gen);
    for (int i = 0; i < (int)circuit.H_gen.size(); ++i) {
        assert(circuit.H_gen[i] == example.H_gen[i]);
    }
#endif // #ifdef DEBUG

    // permute circuit.H_gen according to the wire permutation
    circuit.H_gen_permute.resize(num_hgen * circuit.num_gates, Field(0));
    plonk_permute_subgroup_H<Field>(
        circuit.H_gen, wire_permutation, circuit.H_gen_permute);
    // TODO: write unit test for plonk_permute_subgroup_H
#ifdef DEBUG
    printf("[%s:%d] circuit.H_gen_permute\n", __FILE__, __LINE__);
    print_vector(circuit.H_gen_permute);
    for (size_t i = 0; i < circuit.H_gen_permute.size(); ++i) {
        assert(circuit.H_gen_permute[i] == example.H_gen_permute[i]);
    }
#endif // #ifdef DEBUG

    // compute the permutation polynomials S_sigma_1, S_sigma_2,
    // S_sigma_3 (see [GWC19], Sect. 8.1) (our indexing starts from 0)
    circuit.S_polys.resize(num_hgen, polynomial<Field>(circuit.num_gates));
    plonk_compute_permutation_polynomials<Field>(
        circuit.H_gen_permute, circuit.num_gates, circuit.S_polys);

    return circuit;
}

} // namespace libsnark

#endif // LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_CIRCUIT_TCC_
