/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include "libsnark/common/constraints_tracker/constraints_tracker.hpp"

#include <exception>
#include <iostream>

namespace libsnark
{

constraints_tracker::~constraints_tracker()
{
    // For now, just dump all entries to stdout

    std::cout << "====================\n"
              << "     CONSTRAINTS\n"
              << "====================\n";

    for (const auto &curve_it : _measurements) {
        std::cout << "\nCURVE " << curve_it.first << ":\n";
        for (const auto &entry_it : curve_it.second) {
            std::cout << "  " << entry_it.first << ": " << entry_it.second
                      << "\n";
        }
    }

    // Ensure everything is output before the process terminates.
    std::cout.flush();
}

void constraints_tracker::add_measurement_for_curve(
    const std::string &curve_name,
    const std::string &name,
    size_t num_constraints)
{
    std::map<std::string, measurements_for_curve>::iterator it =
        _measurements.find(curve_name);
    if (it == _measurements.end()) {
        _measurements[curve_name] = {{name, num_constraints}};
        return;
    }

#ifndef NDEBUG
    for (const auto &entry_it : it->second) {
        if (entry_it.first == name) {
            throw std::runtime_error(
                "duplicate entry: " + name + " (curve: " + curve_name + ")");
        }
    }
#endif

    it->second.emplace_back(name, num_constraints);
}

} // namespace libsnark
