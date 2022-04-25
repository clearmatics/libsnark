/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_COMMON_CONSTRAINTS_TRACKER_HPP_
#define LIBSNARK_COMMON_CONSTRAINTS_TRACKER_HPP_

#include <map>
#include <string>
#include <vector>

namespace libsnark
{

/// Simple class to track a set of measurements (initially the number of
/// constraints required by a given gadget) during test or profiling code, and
/// present them together in a sensible format (rather than interspersed with
/// other output to stdout).
///
/// Intended usage is in a test file, create a static global:
///
///   static constraints_tracker constraints_tracker;
///
/// and then later register measurements:
///
///   TEST(TestSuite, SomeTest)
///   {
///     ...
///     constraints_tracker.add_measurement<libff::bls12_377>(
///       "some_gadget", num_constraints);
///     ...
///   }
///
/// if the unit tests exits cleanly, all measurements should be printed to
/// stdout.
///
/// In the future, this may be expanded to collect more data, or write it in
/// other formats to different places.
class constraints_tracker
{
public:
    ~constraints_tracker();

    template<typename ppT>
    void add_measurement(const std::string &name, size_t num_constraints);

protected:
    using measurement = std::pair<std::string, size_t>;
    using measurements_for_curve = std::vector<measurement>;

    void add_measurement_for_curve(
        const std::string &curve_name,
        const std::string &name,
        size_t num_constraints);

    std::map<std::string, measurements_for_curve> _measurements;
};

template<typename ppT>
void constraints_tracker::add_measurement(
    const std::string &name, size_t num_constraints)
{
    add_measurement_for_curve(ppT::name, name, num_constraints);
}

} // namespace libsnark

#endif // LIBSNARK_COMMON_CONSTRAINTS_TRACKER_HPP_
