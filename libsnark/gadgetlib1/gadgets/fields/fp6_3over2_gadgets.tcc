/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_GADGETLIB1_GADGETS_FIELDS_FP6_3OVER2_GADGETS_TCC_
#define LIBSNARK_GADGETLIB1_GADGETS_FIELDS_FP6_3OVER2_GADGETS_TCC_

#include "libsnark/gadgetlib1/gadgets/fields/fp6_3over2_gadgets.hpp"

namespace libsnark
{

// Fp6_3over2_variable methods

template<typename Fp6T>
Fp6_3over2_variable<Fp6T>::Fp6_3over2_variable(
    protoboard<FieldT> &pb, const std::string &annotation_prefix)
    : gadget<FieldT>(pb, annotation_prefix)
    , _c0(pb, FMT(annotation_prefix, " c0"))
    , _c1(pb, FMT(annotation_prefix, " c1"))
    , _c2(pb, FMT(annotation_prefix, " c2"))
{
}

template<typename Fp6T>
Fp6_3over2_variable<Fp6T>::Fp6_3over2_variable(
    protoboard<FieldT> &pb,
    const Fp6_3over2_variable<Fp6T> &el,
    const std::string &annotation_prefix)
    : gadget<FieldT>(pb, annotation_prefix)
    , _c0(el._c0)
    , _c1(el._c1)
    , _c2(el._c2)
{
}

template<typename Fp6T>
Fp6_3over2_variable<Fp6T>::Fp6_3over2_variable(
    protoboard<FieldT> &pb,
    const Fp6T &el,
    const std::string &annotation_prefix)
    : gadget<FieldT>(pb, annotation_prefix)
    , _c0(pb, el.coeffs[0], FMT(annotation_prefix, " c0"))
    , _c1(pb, el.coeffs[1], FMT(annotation_prefix, " c1"))
    , _c2(pb, el.coeffs[2], FMT(annotation_prefix, " c2"))
{
}

template<typename Fp6T>
Fp6_3over2_variable<Fp6T>::Fp6_3over2_variable(
    protoboard<FieldT> &pb,
    const Fp2_variable<Fp2T> &c0,
    const Fp2_variable<Fp2T> &c1,
    const Fp2_variable<Fp2T> &c2,
    const std::string &annotation_prefix)
    : gadget<FieldT>(pb, annotation_prefix), _c0(c0), _c1(c1), _c2(c2)
{
}

template<typename Fp6T>
Fp6_3over2_variable<Fp6T> Fp6_3over2_variable<Fp6T>::operator*(
    const FieldT &scalar) const
{
    return Fp6_3over2_variable<Fp6T>(
        this->pb,
        _c0 * scalar,
        _c1 * scalar,
        _c2 * scalar,
        FMT(this->annotation_prefix, " fp6_var*scalar"));
}

template<typename Fp6T>
Fp6_3over2_variable<Fp6T> Fp6_3over2_variable<Fp6T>::operator*(
    const Fp2T &fp2_const) const
{
    return Fp6_3over2_variable<Fp6T>(
        this->pb,
        _c0 * fp2_const,
        _c1 * fp2_const,
        _c2 * fp2_const,
        FMT(this->annotation_prefix, " fp6_var*fp2_const"));
}

template<typename Fp6T>
Fp6_3over2_variable<Fp6T> Fp6_3over2_variable<Fp6T>::operator*(
    const Fp6T &fp6_const) const
{
    // c0 = a0*b0 + non_residue * (a1*b2 + a2*b1)
    // c1 = a0*b1 + a1*b0 + non_residue * a2*b2
    // c3 = a0*b2 + a1*b1 + a2*b0
    return Fp6_3over2_variable<Fp6T>(
        this->pb,
        _c0 * fp6_const.coeffs[0] +
            (_c1 * fp6_const.coeffs[2] + _c2 * fp6_const.coeffs[1]) *
                Fp6T::non_residue,
        _c0 * fp6_const.coeffs[1] + _c1 * fp6_const.coeffs[0] +
            _c2 * fp6_const.coeffs[2] * Fp6T::non_residue,
        _c0 * fp6_const.coeffs[2] + _c1 * fp6_const.coeffs[1] +
            _c2 * fp6_const.coeffs[0],
        FMT(this->annotation_prefix, " fp6_var*fp6_const"));
}

template<typename Fp6T>
Fp6_3over2_variable<Fp6T> Fp6_3over2_variable<Fp6T>::operator+(
    const Fp6_3over2_variable<Fp6T> &other) const
{
    return Fp6_3over2_variable<Fp6T>(
        this->pb,
        _c0 + other._c0,
        _c1 + other._c1,
        _c2 + other._c2,
        FMT(this->annotation_prefix.c_str(), " +other"));
}

template<typename Fp6T>
Fp6_3over2_variable<Fp6T> Fp6_3over2_variable<Fp6T>::operator-(
    const Fp6_3over2_variable<Fp6T> &other) const
{
    return Fp6_3over2_variable<Fp6T>(
        this->pb,
        _c0 - other._c0,
        _c1 - other._c1,
        _c2 - other._c2,
        FMT(this->annotation_prefix, " -other"));
}

template<typename Fp6T>
Fp6_3over2_variable<Fp6T> Fp6_3over2_variable<Fp6T>::operator-() const
{
    return Fp6_3over2_variable<Fp6T>(
        this->pb,
        -_c0,
        -_c1,
        -_c2,
        FMT(this->annotation_prefix, " fp6_negate"));
}

template<typename Fp6T>
Fp6_3over2_variable<Fp6T> Fp6_3over2_variable<Fp6T>::frobenius_map(
    size_t power) const
{
    return Fp6_3over2_variable<Fp6T>(
        this->pb,
        _c0.frobenius_map(power),
        _c1.frobenius_map(power) * Fp6T::Frobenius_coeffs_c1[power % 6],
        _c2.frobenius_map(power) * Fp6T::Frobenius_coeffs_c2[power % 6],
        FMT(this->annotation_prefix, " fp6_frobenius_map"));
}

template<typename Fp6T> void Fp6_3over2_variable<Fp6T>::evaluate() const
{
    _c0.evaluate();
    _c1.evaluate();
    _c2.evaluate();
}

template<typename Fp6T>
void Fp6_3over2_variable<Fp6T>::generate_r1cs_witness(const Fp6T &el)
{
    _c0.generate_r1cs_witness(el.coeffs[0]);
    _c1.generate_r1cs_witness(el.coeffs[1]);
    _c2.generate_r1cs_witness(el.coeffs[2]);
}

template<typename Fp6T> Fp6T Fp6_3over2_variable<Fp6T>::get_element() const
{
    return Fp6T(_c0.get_element(), _c1.get_element(), _c2.get_element());
}

// Fp6_3over2_mul_gadget methods

template<typename Fp6T>
Fp6_3over2_mul_gadget<Fp6T>::Fp6_3over2_mul_gadget(
    protoboard<FieldT> &pb,
    const Fp6_3over2_variable<Fp6T> &A,
    const Fp6_3over2_variable<Fp6T> &B,
    const Fp6_3over2_variable<Fp6T> &result,
    const std::string &annotation_prefix)
    : gadget<FieldT>(pb, annotation_prefix)
    , _A(A)
    , _B(B)
    , _result(result)
    , _compute_v1(
          pb,
          A._c1,
          B._c1,
          Fp2_variable<Fp2T>(pb, FMT(annotation_prefix, " v1")),
          FMT(annotation_prefix, " _compute_v1"))
    , _compute_v2(
          pb,
          A._c2,
          B._c2,
          Fp2_variable<Fp2T>(pb, FMT(annotation_prefix, " v2")),
          FMT(annotation_prefix, " _compute_v2"))
    , _compute_a1a2_times_b1b2(
          pb,
          A._c1 + A._c2,
          B._c1 + B._c2,
          Fp2_variable<Fp2T>(pb, FMT(annotation_prefix, " (a1+a2)*(b1+b2)")),
          FMT(annotation_prefix, " _compute_a1a2_times_b1b2"))
    // c0 = a0*b0 + non_residue*((a1 + a2)(b1 + b2) - a1*b1 - a2*b2)
    , _compute_v0(
          pb,
          A._c0,
          B._c0,
          _result._c0 - (_compute_a1a2_times_b1b2.result - _compute_v1.result -
                         _compute_v2.result) *
                            Fp6T::non_residue,
          FMT(annotation_prefix, " _compute_v0"))
    // c1 = (a0 + a1)(b0 + b1) - a0*b0 - a1*b1 + non_residue * a2*b2
    , _compute_a0a1_times_b0b1(
          pb,
          A._c0 + A._c1,
          B._c0 + B._c1,
          _result._c1 + _compute_v0.result + _compute_v1.result -
              _compute_v2.result * Fp6T::non_residue,
          FMT(annotation_prefix, " _compute_a0a1_times_b0b1"))
    // c2 = (a0 + a2)(b0 + b2) - a0*b0 - a2*b2 + a1*b1
    , _compute_a0a2_times_b0b2(
          pb,
          A._c0 + A._c2,
          B._c0 + B._c2,
          _result._c2 + _compute_v0.result + _compute_v2.result -
              _compute_v1.result,
          FMT(annotation_prefix, " _compute_a0a2_times_b0b2"))
{
}

template<typename Fp6T>
void Fp6_3over2_mul_gadget<Fp6T>::generate_r1cs_constraints()
{
    _compute_v1.generate_r1cs_constraints();
    _compute_v2.generate_r1cs_constraints();
    _compute_a1a2_times_b1b2.generate_r1cs_constraints();
    _compute_v0.generate_r1cs_constraints();
    _compute_a0a1_times_b0b1.generate_r1cs_constraints();
    _compute_a0a2_times_b0b2.generate_r1cs_constraints();
}

template<typename Fp6T>
void Fp6_3over2_mul_gadget<Fp6T>::generate_r1cs_witness()
{
    const Fp2T a0 = _A._c0.get_element();
    const Fp2T a1 = _A._c1.get_element();
    const Fp2T a2 = _A._c2.get_element();
    const Fp2T b0 = _B._c0.get_element();
    const Fp2T b1 = _B._c1.get_element();
    const Fp2T b2 = _B._c2.get_element();

    // c0 = v1 + non_residue*((a1 + a2)(b1 + b2) - v1 - v2)
    _compute_v1.generate_r1cs_witness();
    const Fp2T v1 = _compute_v1.result.get_element();
    _compute_v2.generate_r1cs_witness();
    const Fp2T v2 = _compute_v2.result.get_element();
    _compute_a1a2_times_b1b2.A.evaluate();
    _compute_a1a2_times_b1b2.B.evaluate();
    _compute_a1a2_times_b1b2.generate_r1cs_witness();
    const Fp2T a1a2_times_b1b2 = _compute_a1a2_times_b1b2.result.get_element();
    const Fp2T v0 = a0 * b0;
    _result._c0.generate_r1cs_witness(
        v0 + Fp6T::mul_by_non_residue(a1a2_times_b1b2 - v1 - v2));
    _compute_v0.generate_r1cs_witness();

    // c1 = (a0 + a1)(b0 + b1) - v1 - v1 + non_residue * v2
    const Fp2T a0a1_times_b0b1 = (a0 + a1) * (b0 + b1);
    _result._c1.generate_r1cs_witness(
        a0a1_times_b0b1 - v0 - v1 + Fp6T::mul_by_non_residue(v2));
    _compute_a0a1_times_b0b1.A.evaluate();
    _compute_a0a1_times_b0b1.B.evaluate();
    _compute_a0a1_times_b0b1.generate_r1cs_witness();

    // c2 = (a0 + a2)(b0 + b2) - v1 - v2 + v1
    const Fp2T a0a2_times_b0b2 = (a0 + a2) * (b0 + b2);
    _result._c2.generate_r1cs_witness(a0a2_times_b0b2 - v0 - v2 + v1);
    _compute_a0a2_times_b0b2.A.evaluate();
    _compute_a0a2_times_b0b2.B.evaluate();
    _compute_a0a2_times_b0b2.generate_r1cs_witness();
}

} // namespace libsnark

#endif // LIBSNARK_GADGETLIB1_GADGETS_FIELDS_FP6_3OVER2_GADGETS_TCC_
