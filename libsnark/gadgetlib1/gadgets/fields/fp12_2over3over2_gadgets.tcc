/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_GADGETLIB1_GADGETS_FIELDS_FP12_2OVER3OVER2_GADGETS_TCC_
#define LIBSNARK_GADGETLIB1_GADGETS_FIELDS_FP12_2OVER3OVER2_GADGETS_TCC_

#include "libsnark/gadgetlib1/gadgets/fields/fp12_2over3over2_gadgets.hpp"

namespace libsnark
{

// Fp12_2over3over2_variable methods

template<typename Fp12T>
Fp12_2over3over2_variable<Fp12T>::Fp12_2over3over2_variable(
    protoboard<FieldT> &pb, const std::string &annotation_prefix)
    : gadget<FieldT>(pb, annotation_prefix)
    , _c0(pb, FMT(annotation_prefix, " c0"))
    , _c1(pb, FMT(annotation_prefix, " c1"))
{
}

template<typename Fp12T>
Fp12_2over3over2_variable<Fp12T>::Fp12_2over3over2_variable(
    protoboard<FieldT> &pb,
    const Fp12T &el,
    const std::string &annotation_prefix)
    : gadget<FieldT>(pb, annotation_prefix)
    , _c0(pb, el.coeffs[0], FMT(annotation_prefix, " c0"))
    , _c1(pb, el.coeffs[1], FMT(annotation_prefix, " c1"))
{
}

template<typename Fp12T>
Fp12_2over3over2_variable<Fp12T>::Fp12_2over3over2_variable(
    protoboard<FieldT> &pb,
    const Fp6_3over2_variable<Fp6T> &c0,
    const Fp6_3over2_variable<Fp6T> &c1,
    const std::string &annotation_prefix)
    : gadget<FieldT>(pb, annotation_prefix), _c0(c0), _c1(c1)
{
}

template<typename Fp12T>
Fp12_2over3over2_variable<Fp12T> Fp12_2over3over2_variable<Fp12T>::operator*(
    const Fp2T &fp2_const) const
{
    return Fp12_2over3over2_variable(
        this->pb,
        _c0 * fp2_const,
        _c1 * fp2_const,
        FMT(this->annotation_prefix, " fp12_var*fp2_const"));
}

template<typename Fp12T>
Fp12_2over3over2_variable<Fp12T> Fp12_2over3over2_variable<Fp12T>::operator*(
    const Fp12T &fp12_const) const
{
    using Fp6_variable = Fp6_3over2_variable<Fp6T>;
    const Fp6_variable a0_b0 = _c0 * fp12_const.coeffs[0];
    const Fp6_variable a1_b1 = _c1 * fp12_const.coeffs[1];
    const Fp6_variable a0a1_b0b1 =
        (_c0 + _c1) * (fp12_const.coeffs[0] + fp12_const.coeffs[1]);

    return Fp12_2over3over2_variable(
        this->pb,
        a0_b0 + (a1_b1 * Fp6T(Fp2T::zero(), Fp2T::one(), Fp2T::zero())),
        a0a1_b0b1 - a0_b0 - a1_b1,
        FMT(this->annotation_prefix, " fp12_var*fp12_const"));
}

template<typename Fp12T>
Fp12_2over3over2_variable<Fp12T> Fp12_2over3over2_variable<
    Fp12T>::frobenius_map(size_t power) const
{
    return Fp12_2over3over2_variable(
        this->pb,
        _c0.frobenius_map(power),
        _c1.frobenius_map(power) * Fp12T::Frobenius_coeffs_c1[power % 12],
        FMT(this->annotation_prefix, " fp12_frobenius_map"));
}

template<typename Fp12T>
Fp12_2over3over2_variable<Fp12T> Fp12_2over3over2_variable<
    Fp12T>::unitary_inverse() const
{
    return Fp12_2over3over2_variable(
        this->pb, _c0, -_c1, FMT(this->annotation_prefix, " fp12_unitary_inv"));
}

template<typename Fp12T> void Fp12_2over3over2_variable<Fp12T>::evaluate() const
{
    _c0.evaluate();
    _c1.evaluate();
}

template<typename Fp12T>
void Fp12_2over3over2_variable<Fp12T>::generate_r1cs_witness(const Fp12T &el)
{
    _c0.generate_r1cs_witness(el.coeffs[0]);
    _c1.generate_r1cs_witness(el.coeffs[1]);
}

template<typename Fp12T>
Fp12T Fp12_2over3over2_variable<Fp12T>::get_element() const
{
    return Fp12T(_c0.get_element(), _c1.get_element());
}

// Multiplication of Fp6 elements by Fp12::non-residue and
// Fp12::non-residue^{-1}.

template<typename Fp12T>
Fp6_3over2_variable<typename Fp12T::my_Fp6> fp6_mul_by_non_residue(
    protoboard<typename Fp12T::my_Fp> &pb,
    const Fp6_3over2_variable<typename Fp12T::my_Fp6> &c,
    const std::string &annotation_prefix)
{
    return Fp6_3over2_variable<typename Fp12T::my_Fp6>(
        pb, c._c2 * Fp12T::non_residue, c._c0, c._c1, annotation_prefix);
}

template<typename Fp12T>
Fp6_3over2_variable<typename Fp12T::my_Fp6> fp6_mul_by_non_residue_inverse(
    protoboard<typename Fp12T::my_Fp> &pb,
    const Fp6_3over2_variable<typename Fp12T::my_Fp6> &c,
    const std::string &annotation_prefix)
{
    return Fp6_3over2_variable<typename Fp12T::my_Fp6>(
        pb,
        c._c1,
        c._c2,
        c._c0 * Fp12T::non_residue.inverse(),
        annotation_prefix);
}

// Fp12_2over3over2_square_gadget methods

template<typename Fp12T>
Fp12_2over3over2_square_gadget<Fp12T>::Fp12_2over3over2_square_gadget(
    protoboard<FieldT> &pb,
    const Fp12_2over3over2_variable<Fp12T> &A,
    const Fp12_2over3over2_variable<Fp12T> &result,
    const std::string &annotation_prefix)
    : gadget<FieldT>(pb, annotation_prefix)
    , _A(A)
    , _result(result)
    , _compute_alpha(
          pb,
          _A._c0,
          _A._c1,
          _result._c1 * (FieldT("2").inverse()),
          FMT(annotation_prefix, " _compute_alpha"))
    , _compute_beta(
          pb,
          _A._c0 + _A._c1,
          _A._c0 + fp6_mul_by_non_residue<Fp12T>(
                       pb, _A._c1, FMT(annotation_prefix, " a1_times_v")),
          _result._c0 +
              fp6_mul_by_non_residue<Fp12T>(
                  pb,
                  _compute_alpha._result,
                  FMT(annotation_prefix, " alpha_times_v")) +
              _compute_alpha._result,
          FMT(annotation_prefix, " _compute_beta"))
{
}

template<typename Fp12T>
const Fp12_2over3over2_variable<Fp12T>
    &Fp12_2over3over2_square_gadget<Fp12T>::result() const
{
    return _result;
}

template<typename Fp12T>
void Fp12_2over3over2_square_gadget<Fp12T>::generate_r1cs_constraints()
{
    _compute_alpha.generate_r1cs_constraints();
    _compute_beta.generate_r1cs_constraints();
}

template<typename Fp12T>
void Fp12_2over3over2_square_gadget<Fp12T>::generate_r1cs_witness()
{
    const Fp6T a0 = _A._c0.get_element();
    const Fp6T a1 = _A._c1.get_element();
    const Fp6T alpha = a0 * a1;
    _result._c1.generate_r1cs_witness(alpha + alpha);
    _compute_alpha.generate_r1cs_witness();

    const Fp6T beta = (a0 + a1) * (a0 + Fp12T::mul_by_non_residue(a1));
    _result._c0.generate_r1cs_witness(
        beta - Fp12T::mul_by_non_residue(alpha) - alpha);
    _compute_beta._A.evaluate();
    _compute_beta._B.evaluate();
    _compute_beta.generate_r1cs_witness();
}

// Fp12_2over3over2_mul_by_024_gadget methods

template<typename Fp12T>
Fp12_2over3over2_mul_by_024_gadget<Fp12T>::Fp12_2over3over2_mul_by_024_gadget(
    protoboard<FieldT> &pb,
    const Fp12_2over3over2_variable<Fp12T> &Z,
    const Fp2_variable<Fp2T> &X_0,
    const Fp2_variable<Fp2T> &X_2,
    const Fp2_variable<Fp2T> &X_4,
    const Fp12_2over3over2_variable<Fp12T> &result,
    const std::string &annotation_prefix)
    : gadget<FieldT>(pb, annotation_prefix)
    , _Z(Z)
    , _X_0(X_0)
    , _X_2(X_2)
    , _X_4(X_4)
    // out_z0 = z0*x0 + non_residue * (z1*x2 + z4*x4)
    // => z0 * x0 = out_z0 - non_residue * (z1*x2 + z4*x4)
    , _compute_z1_x2(
          pb,
          Z._c0._c1,
          X_2,
          Fp2_variable<Fp2T>(pb, FMT(annotation_prefix, " z1_x2")),
          FMT(annotation_prefix, "_compute_z1_x2"))
    , _compute_z4_x4(
          pb,
          Z._c1._c1,
          X_4,
          Fp2_variable<Fp2T>(pb, FMT(annotation_prefix, " z4_x4")),
          FMT(annotation_prefix, " _compute_z4_x4"))
    , _compute_z0_x0(
          pb,
          Z._c0._c0,
          X_0,
          result._c0._c0 + ((_compute_z1_x2.result + _compute_z4_x4.result) *
                            -Fp6T::non_residue),
          FMT(annotation_prefix, " _compute_z0_x0"))
    // out_z1 = z1*x0 + non_residue * (z2*x2 + z5*x4)
    // => z1 * z0 = out_z1 - non_residue * (z2*x2 + z5*x4)
    , _compute_z2_x2(
          pb,
          Z._c0._c2,
          X_2,
          Fp2_variable<Fp2T>(pb, FMT(annotation_prefix, " z2_x2")),
          FMT(annotation_prefix, " _compute_z2_x2"))
    , _compute_z5_x4(
          pb,
          Z._c1._c2,
          X_4,
          Fp2_variable<Fp2T>(pb, FMT(annotation_prefix, " z5_x4")),
          FMT(annotation_prefix, " _compute_z5_x4"))
    , _compute_z1_x0(
          pb,
          Z._c0._c1,
          X_0,
          result._c0._c1 + ((_compute_z2_x2.result + _compute_z5_x4.result) *
                            -Fp6T::non_residue),
          FMT(annotation_prefix, " _compute_z1_x0"))
    // z0*x2 + z2*x0 = (z0 + z2)*(x0 + x2) - z0*x0 - z2*x2
    // out_z2 = z0*x2 + z2*x0 + z3*x4
    //        = (z0 + z2)*(x0 + x2) - z0*x0 - z2*x2 + z3*x4
    // => (z0 + z2)*(x0 + x2) = out_z2 + z0*x0 + z2*x2 - z3*x4
    , _compute_z3_x4(
          pb,
          Z._c1._c0,
          X_4,
          Fp2_variable<Fp2T>(pb, FMT(annotation_prefix, " z3_x4")),
          FMT(annotation_prefix, " _compute_z3_x4"))
    , _compute_z02_x02(
          pb,
          Z._c0._c0 + Z._c0._c2,
          X_0 + X_2,
          result._c0._c2 + _compute_z0_x0.result + _compute_z2_x2.result +
              (_compute_z3_x4.result * -FieldT::one()),
          FMT(annotation_prefix, " _compute_z02_x02"))
    // z2*x4 + z4*x2 = (z2 + z4)*(x2 + x4) - z2*x2 - z4*x4
    // out_z3 = z3*x0 + non_residue * (z2*x4 + z4*x2)
    //        = z3*x0 + non_residue * ((z2 + z4)*(x2 + x4) - z2*x2 - z4*x4)
    // => (z2 + z4)*(x2 + x4) = (out_z3 - z3*x0) / non_residue + z2*x2 + z4_x4
    , _compute_z3_x0(
          pb,
          Z._c1._c0,
          X_0,
          Fp2_variable<Fp2T>(pb, FMT(annotation_prefix, " z3_x0")),
          FMT(annotation_prefix, " _compute_z3_x0"))
    , _compute_z24_x24(
          pb,
          Z._c0._c2 + Z._c1._c1,
          X_2 + X_4,
          _compute_z2_x2.result + _compute_z4_x4.result +
              (result._c1._c0 + _compute_z3_x0.result * -FieldT::one()) *
                  Fp6T::non_residue.inverse(),
          FMT(annotation_prefix, " _compute_z24_x24"))
    // z0*x4 + z4*x0 = (z0 + z4)*(x0 + x4) - z0*x0 - z4*x4
    // out_z4 = z0*x4 + z4*x0 + non_residue * z5*x2
    //        = (z0 + z4)*(x0 + x4) - z0*x0 - z4*x4 + non_residue * z5*x2
    // => (z0 + z4)*(x0 + x4) = out_z4 + z0*x0 + z4*x4 - non_residue * z5*x2
    , _compute_z5_x2(
          pb,
          Z._c1._c2,
          X_2,
          Fp2_variable<Fp2T>(pb, FMT(annotation_prefix, " z5_x2")),
          FMT(annotation_prefix, " _compute_z5_x2"))
    , _compute_z04_x04(
          pb,
          Z._c0._c0 + Z._c1._c1,
          X_0 + X_4,
          result._c1._c1 + _compute_z0_x0.result + _compute_z4_x4.result +
              _compute_z5_x2.result * -Fp6T::non_residue,
          FMT(annotation_prefix, " _compute_z04_x04"))
    // S = z1_x0 + z1_x2 + z3_x0 + z3*x4 + z5_x2 + z5*x4
    // out_z5 = z1*x4 + z3*x2 + z5*x0
    //        = (z1 + z3 + z5)*(x0 + x2 + x4) - S
    // => (z1 + z3 + z5)*(x0 + x2 + x4) = out_z5 + S
    , _S(_compute_z1_x2.result + _compute_z1_x0.result + _compute_z5_x4.result +
         _compute_z3_x4.result + _compute_z3_x0.result + _compute_z5_x2.result)
    , _compute_out_z5_plus_S(
          pb,
          Z._c0._c1 + Z._c1._c0 + Z._c1._c2,
          X_0 + X_2 + X_4,
          result._c1._c2 + _S,
          FMT(annotation_prefix, " _compute_out_z5_plus_S"))
    , _result(result)
{
}

template<typename Fp12T>
const Fp12_2over3over2_variable<Fp12T>
    &Fp12_2over3over2_mul_by_024_gadget<Fp12T>::result() const
{
    return _result;
}

template<typename Fp12T>
void Fp12_2over3over2_mul_by_024_gadget<Fp12T>::generate_r1cs_constraints()
{
    _compute_z1_x2.generate_r1cs_constraints();
    _compute_z4_x4.generate_r1cs_constraints();
    _compute_z0_x0.generate_r1cs_constraints();
    _compute_z2_x2.generate_r1cs_constraints();
    _compute_z5_x4.generate_r1cs_constraints();
    _compute_z1_x0.generate_r1cs_constraints();
    _compute_z3_x4.generate_r1cs_constraints();
    _compute_z02_x02.generate_r1cs_constraints();
    _compute_z3_x0.generate_r1cs_constraints();
    _compute_z24_x24.generate_r1cs_constraints();
    _compute_z5_x2.generate_r1cs_constraints();
    _compute_z04_x04.generate_r1cs_constraints();
    _compute_out_z5_plus_S.generate_r1cs_constraints();
}

template<typename Fp12T>
void Fp12_2over3over2_mul_by_024_gadget<Fp12T>::generate_r1cs_witness()
{
    const Fp2T z0 = _Z._c0._c0.get_element();
    const Fp2T z1 = _Z._c0._c1.get_element();
    const Fp2T z2 = _Z._c0._c2.get_element();
    const Fp2T z3 = _Z._c1._c0.get_element();
    const Fp2T z4 = _Z._c1._c1.get_element();
    const Fp2T z5 = _Z._c1._c2.get_element();

    const Fp2T x0 = _X_0.get_element();
    const Fp2T x2 = _X_2.get_element();
    const Fp2T x4 = _X_4.get_element();

    // out_z0 = z0*x0 + non_residue * (z1*x2 + z4*x4)
    // => z0 * x0 = out_z0 - non_residue * (z1*x2 + z4*x4)
    _compute_z1_x2.generate_r1cs_witness();
    _compute_z4_x4.generate_r1cs_witness();
    const Fp2T z0_x0 = z0 * x0;
    const Fp2T z4_x4 = _compute_z4_x4.result.get_element();
    _result._c0._c0.generate_r1cs_witness(
        z0_x0 +
        Fp6T::non_residue * (_compute_z1_x2.result.get_element() + z4_x4));
    _compute_z0_x0.generate_r1cs_witness();

    // out_z1 = z1*x0 + non_residue * (z2*x2 + z5*x4)
    // => z1 * z0 = out_z1 - non_residue * (z2*x2 + z5*x4)
    _compute_z2_x2.generate_r1cs_witness();
    _compute_z5_x4.generate_r1cs_witness();
    const Fp2T z2_x2 = _compute_z2_x2.result.get_element();
    const Fp2T z1_x0 = z1 * x0;
    _result._c0._c1.generate_r1cs_witness(
        z1_x0 +
        Fp6T::non_residue * (z2_x2 + _compute_z5_x4.result.get_element()));
    _compute_z1_x0.generate_r1cs_witness();

    // z0*x2 + z2*x0 = (z0 + z2)*(x0 + x2) - z0*x0 - z2*x2
    // out_z2 = z0*x2 + z2*x0 + z3*x4
    //        = (z0 + z2)*(x0 + x2) - z0*x0 - z2*x2 + z3*x4
    // => (z0 + z2)*(x0 + x2) = out_z2 + z0*x0 + z2*x2 - z3*x4
    _compute_z3_x4.generate_r1cs_witness();
    const Fp2T z3_x4 = _compute_z3_x4.result.get_element();
    _result._c0._c2.generate_r1cs_witness(
        (z0 + z2) * (x0 + x2) - z0_x0 - z2_x2 + z3_x4);
    _compute_z02_x02.A.evaluate();
    _compute_z02_x02.B.evaluate();
    _compute_z02_x02.generate_r1cs_witness();

    // z2*x4 + z4*x2 = (z2 + z4)*(x2 + x4) - z2*x2 - z4*x4
    // out_z3 = z3*x0 + non_residue * (z2*x4 + z4*x2)
    //        = z3*x0 + non_residue * ((z2 + z4)*(x2 + x4) - z2*x2 - z4*x4)
    // => (z2 + z4)*(x2 + x4) = (out_z3 - z3*x0) / non_residue + z2*x2 + z4_x4
    _compute_z3_x0.generate_r1cs_witness();
    const Fp2T z3_x0 = _compute_z3_x0.result.get_element();
    _result._c1._c0.generate_r1cs_witness(
        z3_x0 + Fp6T::non_residue * ((z2 + z4) * (x2 + x4) - z2_x2 - z4_x4));
    _compute_z24_x24.A.evaluate();
    _compute_z24_x24.B.evaluate();
    _compute_z24_x24.generate_r1cs_witness();

    // z0*x4 + z4*x0 = (z0 + z4)*(x0 + x4) - z0*x0 - z4*x4
    // out_z4 = z0*x4 + z4*x0 + non_residue * z5*x2
    //        = (z0 + z4)*(x0 + x4) - z0*x0 - z4*x4 + non_residue * z5*x2
    // => (z0 + z4)*(x0 + x4) = out_z4 + z0*x0 + z4*x4 - non_residue * z5*x2
    _compute_z5_x2.generate_r1cs_witness();
    const Fp2T z5_x2 = _compute_z5_x2.result.get_element();
    _result._c1._c1.generate_r1cs_witness(
        (z0 + z4) * (x0 + x4) - z0_x0 - z4_x4 + Fp6T::non_residue * z5_x2);
    _compute_z04_x04.A.evaluate();
    _compute_z04_x04.B.evaluate();
    _compute_z04_x04.generate_r1cs_witness();

    // S = z1_x0 - z1_x2 - z3_x0 - z3*x4 - z5_x2 - z5*x4
    // out_z5 = z1*x4 + z3*x2 + z5*x0
    //        = (z1 + z3 + z5)*(x0 + x2 + x4) - S
    // => (z1 + z3 + z5)*(x0 + x2 + x4) = out_z5 + S
    _S.evaluate();
    const Fp2T S = _S.get_element();
    _result._c1._c2.generate_r1cs_witness((z1 + z3 + z5) * (x0 + x2 + x4) - S);
    _compute_out_z5_plus_S.A.evaluate();
    _compute_out_z5_plus_S.B.evaluate();
    _compute_out_z5_plus_S.generate_r1cs_witness();
}

// Fp12_2over3over2_mul_gadget methods

template<typename Fp12T>
Fp12_2over3over2_mul_gadget<Fp12T>::Fp12_2over3over2_mul_gadget(
    protoboard<FieldT> &pb,
    const Fp12_2over3over2_variable<Fp12T> &A,
    const Fp12_2over3over2_variable<Fp12T> &B,
    const Fp12_2over3over2_variable<Fp12T> &result,
    const std::string &annotation_prefix)
    : gadget<FieldT>(pb, annotation_prefix)
    , _A(A)
    , _B(B)
    , _result(result)
    , _compute_v0(
          pb,
          A._c0,
          B._c0,
          Fp6_3over2_variable<Fp6T>(pb, FMT(annotation_prefix, " v0")),
          " _compute_v0")
    // result0 = a0*b0 + non_residue*a1*b1
    //   <=> a1*b1 = (result0 - a0*b0) * non_residue.inverse
    , _compute_v1(
          pb,
          A._c1,
          B._c1,
          fp6_mul_by_non_residue_inverse<Fp12T>(
              pb,
              _result._c0 - _compute_v0._result,
              FMT(annotation_prefix, " v1")),
          FMT(annotation_prefix, " _compute_v1"))
    // result1 = a0*b1 + a1*b0 = (a0 + a1)*(b0 + b1) - a0*b0 - a1*b1
    //   <=> (a0 + a1)(b0 + b1) = result1 + a0*b0 + a1*b1
    , _compute_a0_plus_a1_times_b0_plus_b1(
          pb,
          A._c0 + A._c1,
          B._c0 + B._c1,
          _result._c1 + _compute_v0._result + _compute_v1._result,
          FMT(annotation_prefix, " _compute_a0_plus_a1_times_b0_plus_b1"))
{
}

template<typename Fp12T>
const Fp12_2over3over2_variable<Fp12T>
    &Fp12_2over3over2_mul_gadget<Fp12T>::result() const
{
    return _result;
}

template<typename Fp12T>
void Fp12_2over3over2_mul_gadget<Fp12T>::generate_r1cs_constraints()
{
    _compute_v0.generate_r1cs_constraints();
    _compute_v1.generate_r1cs_constraints();
    _compute_a0_plus_a1_times_b0_plus_b1.generate_r1cs_constraints();
}

template<typename Fp12T>
void Fp12_2over3over2_mul_gadget<Fp12T>::generate_r1cs_witness()
{
    _compute_v0.generate_r1cs_witness();
    const Fp6T a0 = _compute_v0._A.get_element();
    const Fp6T a1 = _compute_v1._A.get_element();
    const Fp6T b0 = _compute_v0._B.get_element();
    const Fp6T b1 = _compute_v1._B.get_element();
    const Fp6T a0b0 = _compute_v0._result.get_element();
    const Fp6T a1b1 = a1 * b1;

    _result._c0.generate_r1cs_witness(a0b0 + Fp12T::mul_by_non_residue(a1b1));
    _compute_v1._result.evaluate();
    _compute_v1.generate_r1cs_witness();

    _result._c1.generate_r1cs_witness((a0 + a1) * (b0 + b1) - a0b0 - a1b1);

    _compute_a0_plus_a1_times_b0_plus_b1._A.evaluate();
    _compute_a0_plus_a1_times_b0_plus_b1._B.evaluate();
    _compute_a0_plus_a1_times_b0_plus_b1.generate_r1cs_witness();
}

// Fp12_2over3over2_inv_gadget methods

template<typename Fp12T>
Fp12_2over3over2_inv_gadget<Fp12T>::Fp12_2over3over2_inv_gadget(
    protoboard<FieldT> &pb,
    const Fp12_2over3over2_variable<Fp12T> &A,
    const Fp12_2over3over2_variable<Fp12T> &result,
    const std::string &annotation_prefix)
    : gadget<FieldT>(pb, annotation_prefix)
    , _A(A)
    , _result(result)
    // _result == A^{-1}
    //   <=> _result * A == Fp12::one()
    , _compute_A_times_result(
          pb,
          _A,
          _result,
          Fp12_2over3over2_variable<Fp12T>(
              pb,
              Fp6_3over2_variable<Fp6T>(pb, Fp6T::one(), " (A*A.inv).c0"),
              Fp6_3over2_variable<Fp6T>(pb, Fp6T::zero(), " (A*A.inv).c1"),
              FMT(annotation_prefix, " A*A.mult_inv")),
          FMT(annotation_prefix, " _compute_A_times_result"))
{
}

template<typename Fp12T>
const Fp12_2over3over2_variable<Fp12T>
    &Fp12_2over3over2_inv_gadget<Fp12T>::result() const
{
    return _result;
}

template<typename Fp12T>
void Fp12_2over3over2_inv_gadget<Fp12T>::generate_r1cs_constraints()
{
    _compute_A_times_result.generate_r1cs_constraints();
}

template<typename Fp12T>
void Fp12_2over3over2_inv_gadget<Fp12T>::generate_r1cs_witness()
{
    _result.generate_r1cs_witness(_A.get_element().inverse());
    _compute_A_times_result.generate_r1cs_witness();
}

// Fp12_2over3over2_cyclotomic_square_gadget methods

template<typename Fp12T>
Fp12_2over3over2_cyclotomic_square_gadget<Fp12T>::
    Fp12_2over3over2_cyclotomic_square_gadget(
        protoboard<FieldT> &pb,
        const Fp12_2over3over2_variable<Fp12T> &A,
        const Fp12_2over3over2_variable<Fp12T> &result,
        const std::string &annotation_prefix)
    : gadget<FieldT>(pb, annotation_prefix)
    , _A(A)
    , _result(result)
    // z0z4 = 6^{-1} * (result4 - 2*z4)
    , _compute_z0z4(
          pb,
          _A._c0._c0,
          _A._c1._c1,
          (_result._c1._c1 - _A._c1._c1 - _A._c1._c1) * FieldT(6).inverse(),
          FMT(annotation_prefix, " _compute_z0z4"))
    // 3*(z0 + z4) * (z0 + non_residue * z4)
    //       = result0 + 3*(1 + non_residue)*z0z4 + 2*z0
    , _check_result_0(
          pb,
          (_A._c0._c0 + _A._c1._c1) * FieldT(3),
          _A._c0._c0 + _A._c1._c1 * Fp6T::non_residue,
          _result._c0._c0 + _A._c0._c0 + _A._c0._c0 +
              _compute_z0z4.result * (Fp2T::one() + Fp6T::non_residue) *
                  FieldT(3),
          FMT(annotation_prefix, " _check_result_0"))
    // z3z2 = 6^{-1} * (result5 - 2*z5)
    , _compute_z3z2(
          pb,
          _A._c1._c0,
          _A._c0._c2,
          (_result._c1._c2 - _A._c1._c2 - _A._c1._c2) * FieldT(6).inverse(),
          FMT(annotation_prefix, " _compute_z3z2"))
    // 3*(z3 + z2)*(z3 + non_residue * z2)
    //       = result1 + 3*(1 + non_residue)*_z3z2 + 2*z1
    , _check_result_1(
          pb,
          (_A._c1._c0 + _A._c0._c2) * FieldT(3),
          _A._c1._c0 + _A._c0._c2 * Fp6T::non_residue,
          _result._c0._c1 + _A._c0._c1 + _A._c0._c1 +
              _compute_z3z2.result * (Fp2T::one() + Fp6T::non_residue) *
                  FieldT(3),
          FMT(annotation_prefix, " _check_result_1"))
    // z1z5 = 6^{-1} * non_residue^{-1} * (result3 - 2*z3)
    , _compute_z1z5(
          pb,
          _A._c0._c1,
          _A._c1._c2,
          (_result._c1._c0 - _A._c1._c0 - _A._c1._c0) *
              Fp6T::non_residue.inverse() * FieldT(6).inverse(),
          FMT(annotation_prefix, " _compute_z1z5"))
    // 3*(z1 + z5)*(z1 + non_residue * z5)
    //       = result2 + 3*(1 + non_residue)*z1z5 + 2*z2
    , _check_result_2(
          pb,
          (_A._c0._c1 + _A._c1._c2) * FieldT(3),
          _A._c0._c1 + _A._c1._c2 * Fp6T::non_residue,
          _result._c0._c2 + _A._c0._c2 + _A._c0._c2 +
              _compute_z1z5.result * (Fp2T::one() + Fp6T::non_residue) *
                  FieldT(3),
          FMT(annotation_prefix, " _check_result_2"))
{
}

template<typename Fp12T>
const Fp12_2over3over2_variable<Fp12T>
    &Fp12_2over3over2_cyclotomic_square_gadget<Fp12T>::result() const
{
    return _result;
}

template<typename Fp12T>
void Fp12_2over3over2_cyclotomic_square_gadget<
    Fp12T>::generate_r1cs_constraints()
{
    _compute_z0z4.generate_r1cs_constraints();
    _check_result_0.generate_r1cs_constraints();
    _compute_z3z2.generate_r1cs_constraints();
    _check_result_1.generate_r1cs_constraints();
    _compute_z1z5.generate_r1cs_constraints();
    _check_result_2.generate_r1cs_constraints();
}

template<typename Fp12T>
void Fp12_2over3over2_cyclotomic_square_gadget<Fp12T>::generate_r1cs_witness()
{
    const Fp2T z0 = _A._c0._c0.get_element();
    const Fp2T z1 = _A._c0._c1.get_element();
    const Fp2T z2 = _A._c0._c2.get_element();
    const Fp2T z3 = _A._c1._c0.get_element();
    const Fp2T z4 = _A._c1._c1.get_element();
    const Fp2T z5 = _A._c1._c2.get_element();

    // result4 = 6 * z0z4 + 2 * z4
    // <=> z0z4 = 6^{-1} * (result4 - 2*z4)
    const Fp2T z0z4 = z0 * z4;
    const Fp2T z0z4_2 = z0z4 + z0z4;
    _result._c1._c1.generate_r1cs_witness(z0z4_2 + z0z4_2 + z0z4_2 + z4 + z4);
    _compute_z0z4.result.evaluate();
    _compute_z0z4.generate_r1cs_witness();

    // result0 = 3*t0_L - 3*t0_R - 2*z0
    //   where
    //     t0_L = (z0 + z4) * (z0 + non_residue * z4)
    //     t0_R = z0z4 * (my_Fp2::one() + my_Fp6::non_residue)
    // <=> 3*(z0 + z4) * (z0 + non_residue * z4)
    //       = result0 + 3*(1 + non_residue)*z0z4 + 2*z0
    const Fp2T t0_L = (z0 + z4) * (z0 + Fp6T::non_residue * z4);
    const Fp2T t0_R = z0z4 * (Fp2T::one() + Fp6T::non_residue);
    _result._c0._c0.generate_r1cs_witness(
        t0_L + t0_L + t0_L - t0_R - t0_R - t0_R - z0 - z0);
    _check_result_0.A.evaluate();
    _check_result_0.B.evaluate();
    _check_result_0.result.evaluate();
    _check_result_0.generate_r1cs_witness();

    // result5 = 6 * z3z2 + 2 * z5
    // <=> z3z2 = 6^{-1} * (result5 - 2*z5)
    const Fp2T z3z2 = z3 * z2;
    const Fp2T z3z2_2 = z3z2 + z3z2;
    _result._c1._c2.generate_r1cs_witness(z3z2_2 + z3z2_2 + z3z2_2 + z5 + z5);
    _compute_z3z2.result.evaluate();
    _compute_z3z2.generate_r1cs_witness();

    // result1 = 3*t2_L - 3*t2_R - 2*z1
    //   where
    //     t2_L = (z3 + z2) * (z3 + non_residue * z2)
    //     t2_R = z3z2 * (1 + non_residue)
    // <=> 3*(z3 + z2)*(z3 + non_residue * z2)
    //       = result1 + 3*(1 + non_residue)*_z3z2 + 2*z1
    const Fp2T t2_L = (z3 + z2) * (z3 + Fp6T::non_residue * z2);
    const Fp2T t2_R = z3z2 * (Fp2T::one() + Fp6T::non_residue);
    _result._c0._c1.generate_r1cs_witness(
        t2_L + t2_L + t2_L - t2_R - t2_R - t2_R - z1 - z1);
    _check_result_1.A.evaluate();
    _check_result_1.B.evaluate();
    _check_result_1.result.evaluate();
    _check_result_1.generate_r1cs_witness();

    // result3 = 6 * non_residue * z1z5 + 2*z3
    // <=> z1z5 = 6^{-1} * non_residue^{-1} * (out3 - 2*z3)
    const Fp2T z1z5 = z1 * z5;
    const Fp2T z1z5_2 = z1z5 + z1z5;
    _result._c1._c0.generate_r1cs_witness(
        (z1z5_2 + z1z5_2 + z1z5_2) * Fp6T::non_residue + z3 + z3);
    _compute_z1z5.result.evaluate();
    _compute_z1z5.generate_r1cs_witness();

    // result2 = 3*t4_L - 3*t4_R - 2*z2
    //   where
    //     t4_L = (z1 + z5) * (z1 + non_residue * z5)
    //     t4_R = z1z5 * (1 + non_residue);
    // <=> 3*(z1 + z5)*(z1 + non_residue * z5)
    //       = result2 + 3*(1 + non_residue)*z1z5 + 2*z2
    const Fp2T t4_L = (z1 + z5) * (z1 + Fp6T::non_residue * z5);
    const Fp2T t4_R = z1z5 * (Fp2T::one() + Fp6T::non_residue);
    _result._c0._c2.generate_r1cs_witness(
        t4_L + t4_L + t4_L - t4_R - t4_R - t4_R - z2 - z2);
    _check_result_2.A.evaluate();
    _check_result_2.B.evaluate();
    _check_result_2.result.evaluate();
    _check_result_2.generate_r1cs_witness();
}

} // namespace libsnark

#endif // LIBSNARK_GADGETLIB1_GADGETS_FIELDS_FP12_2OVER3OVER2_GADGETS_TCC_
