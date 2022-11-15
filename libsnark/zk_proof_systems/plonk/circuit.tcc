/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_CIRCUIT_TCC_
#define LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_CIRCUIT_TCC_

// Implementation of Common Preprocessed Input interfaces for a
// ppzkSNARK for Plonk. See circuit.hpp .

namespace libsnark
{

// TODO: add here function for describing the target circuit through
// the circuit_t structure

template<typename ppT>
circuit_t<ppT>::circuit_t(
    size_t num_gates,
    size_t num_qpolys,
    std::vector<size_t> &&PI_wire_index,
    std::vector<polynomial<Field>> &&Q_polys,
    std::vector<polynomial<Field>> &&S_polys,
    std::vector<std::vector<Field>> &&omega_roots,
    std::vector<Field> &&H_gen,
    std::vector<Field> &&H_gen_permute,
    libff::Fr<ppT> &&k1,
    libff::Fr<ppT> &&k2)
    : num_gates(num_gates)
    , num_qpolys(num_qpolys)
    , PI_wire_index(PI_wire_index)
    , Q_polys(Q_polys)
    , S_polys(S_polys)
    , omega_roots(omega_roots)
    , H_gen(H_gen)
    , H_gen_permute(H_gen_permute)
    , k1(k1)
    , k2(k2)
{
}

} // namespace libsnark

#endif // LIBSNARK_ZK_PROOF_SYSTEMS_PLONK_CIRCUIT_TCC_
