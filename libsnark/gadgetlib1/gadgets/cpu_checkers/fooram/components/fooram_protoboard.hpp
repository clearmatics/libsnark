/** @file
 *****************************************************************************

 Declaration of interfaces for a protoboard for the FOORAM CPU.

 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef FOORAM_PROTOBOARD_HPP_
#define FOORAM_PROTOBOARD_HPP_

#include <libsnark/gadgetlib1/gadget.hpp>
#include <libsnark/relations/ram_computations/rams/fooram/fooram_aux.hpp>

namespace libsnark
{

template<typename FieldT> class fooram_protoboard : public protoboard<FieldT>
{
public:
    const fooram_architecture_params ap;

    fooram_protoboard(const fooram_architecture_params &ap);
};

template<typename FieldT> class fooram_gadget : public gadget<FieldT>
{
protected:
    fooram_protoboard<FieldT> &pb;

public:
    fooram_gadget(
        fooram_protoboard<FieldT> &pb,
        const std::string &annotation_prefix = "");
};

} // namespace libsnark

#include <libsnark/gadgetlib1/gadgets/cpu_checkers/fooram/components/fooram_protoboard.tcc>

#endif // FOORAM_PROTOBOARD_HPP_
