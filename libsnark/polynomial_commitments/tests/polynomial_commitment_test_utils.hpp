/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

namespace libsnark
{

template<typename ppT>
std::vector<libff::Fr<ppT>> gen_polynomial(const size_t degree)
{
    using Field = libff::Fr<ppT>;

    const size_t num_random_values = std::min(4ul, degree);

    std::vector<Field> coefficients;
    coefficients.reserve(degree);

    // First num_random_values are random (ensuring there are no zeroes)
    for (size_t i = 0; i < num_random_values; ++i) {
        coefficients.push_back(Field::random_element());
        if (coefficients[i] == Field::zero()) {
            coefficients[i] = -Field::one();
        }
    }

    // Remaining values coefficients[i] =
    //   coefficients[i-1] * coefficients[i % num_random_values]
    for (size_t i = num_random_values; i < degree; ++i) {
        coefficients.push_back(
            coefficients[i - 1] * coefficients[i % num_random_values]);
    }

    return coefficients;
}

} // namespace libsnark
