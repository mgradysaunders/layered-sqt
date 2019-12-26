/* Copyright (c) 2019 M. Grady Saunders
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 
 *   1. Redistributions of source code must retain the above
 *      copyright notice, this list of conditions and the following
 *      disclaimer.
 * 
 *   2. Redistributions in binary form must reproduce the above
 *      copyright notice, this list of conditions and the following
 *      disclaimer in the documentation and/or other materials
 *      provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
/*+-+*/
#include <fstream>
#include <preform/option_parser.hpp>
#include <layered-sqt/common.hpp>
#include <layered-sqt/file_data.hpp>

int main(int argc, char** argv)
{
    pr::option_parser opt_parser("[OPTIONS] filename");

    // Input filename.
    std::string ifs_filename;
    bool show_incident_dirs = false;
    bool show_bsdf_averages = false;
    int which_slice_index = -1;

    // --show-incident-dirs
    opt_parser.on_option(nullptr, "--show-incident-dirs", 0,
    [&](char**) {
        show_incident_dirs = true;
    })
    << "Show incident directions in report.\n"
       "By default, report only shows number of incident directions.\n";

    // --show-bsdf-averages
    opt_parser.on_option(nullptr, "--show-bsdf-averages", 0,
    [&](char**) {
        show_bsdf_averages = true;
    })
    << "Show BSDF averages in report.\n"
       "By default, report only shows number of BSDF averages.\n";

    // -s/--slice
    opt_parser.on_option("-s", "--slice", 1,
    [&](char** argv) {
        try {
            which_slice_index = std::stoi(argv[0]); // This may throw.
            if (which_slice_index < 0) {
                // Trigger catch block manually.
                throw std::exception(); 
            }
        }
        catch (const std::exception&) {
            throw std::runtime_error(
                    "-s/--slice expects 1 non-negative integer");
        }
    })
    << "Specify slice to show by index.\n"
       "By default, report shows all slices.\n";

    // -h/--help
    opt_parser.on_option("-h", "--help", 0,
    [&](char**) {
        std::cout << opt_parser << std::endl;
        std::exit(EXIT_SUCCESS);
    })
    << "Display this help and exit.\n";

    // Positional.
    opt_parser.on_positional(
    [&](char* argv) {
        if (ifs_filename.empty()) {
            ifs_filename = argv;
        }
        else {
            throw std::runtime_error("more than 1 positional argument");
        }
    });

    try {
        // Parse args.
        opt_parser.parse(argc, argv);
    }
    catch (const std::exception& exception) {
        std::cerr << "Unhandled exception in command line arguments!\n";
        std::cerr << "exception.what(): " << exception.what() << "\n";
        std::exit(EXIT_FAILURE);
    }

    // No input filename?
    if (ifs_filename.empty()) {
        std::cout << opt_parser << std::endl;
        std::exit(EXIT_SUCCESS);
    }

    try {

        // Initialize file data.
        ls::FileData file_data;
        {
            // Open file.
            std::ifstream ifs(
                    ifs_filename,
                    std::ios_base::in |
                    std::ios_base::binary);
            if (!ifs.is_open()) {
                throw std::runtime_error("can't open file");
            }

            // Read file.
            file_data.readLss(ifs);
        }

        // Verify slice index.
        if (which_slice_index >= int(file_data.slices.size())) {
            throw std::runtime_error("slice index out of range");
        }

        // Iterate slices.
        for (const auto& slice : file_data.slices) {
            int slice_index = &slice - &file_data.slices[0];
            if (slice_index != which_slice_index && which_slice_index >= 0) {
                continue;
            }

            // Slice.
            std::cout << "slices[" << slice_index << "]\n";
            std::cout << ".outgoing_angle = " << slice.outgoing_angle << "\n";
            std::cout << ".incident_dirs = ";
            if (show_incident_dirs) {
                if (slice.incident_dirs.empty()) { // Shouldn't happen?
                    std::cout << "[]\n";
                }
                else {
                    // Show array.
                    auto incident_dir = slice.incident_dirs.begin();
                    std::cout << "[\n    ";
                    std::cout << *incident_dir++;
                    while (incident_dir < slice.incident_dirs.end()) {
                        std::cout << ",\n    ";
                        std::cout << *incident_dir++;
                    }
                    std::cout << "\n]";
                    std::cout << "\n";
                }
            }
            else {
                // Show array type only.
                std::cout << "Vec3f[" << slice.incident_dirs.size() << "]\n";
            }
            std::cout << ".bsdf_averages = ";
            if (show_bsdf_averages) { 
                if (slice.bsdf_averages.empty()) { // Shouldn't happen?
                    std::cout << "[]\n";
                }
                else {
                    // Show array.
                    auto bsdf_average = slice.bsdf_averages.begin();
                    std::cout << "[\n    ";
                    std::cout << *bsdf_average++;
                    while (bsdf_average < slice.bsdf_averages.end()) {
                        std::cout << ",\n    ";
                        std::cout << *bsdf_average++;
                    }
                    std::cout << "\n]";
                    std::cout << "\n";
                }
            }
            else {
                // Show array type only.
                std::cout << "Float[" << slice.bsdf_averages.size() << "]\n";
            }
            std::cout << ".path_count = " << slice.path_count << "\n";
            std::cout << "\n";
            std::cout.flush();
        }
    }
    catch (const std::exception& exception) {
        std::cerr << "Unhandled exception!\n";
        std::cerr << "exception.what(): " << exception.what() << "\n";
        std::exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}
