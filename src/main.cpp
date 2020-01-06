/* Copyright (c) 2019-20 M. Grady Saunders
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
#include <preform/thread_pool.hpp>
#include <preform/option_parser.hpp>
#include <layered-sqt/common.hpp>
#include <layered-sqt/medium.hpp>
#include <layered-sqt/layer.hpp>
#include <layered-sqt/layered_assembly.hpp>
#include <layered-sqt/progress_bar.hpp>
#include <layered-sqt/file_data.hpp>

int main(int argc, char** argv)
{
    pr::option_parser opt_parser("[OPTIONS] filename");

    // Options.
    int seed = 0;
    int path_count = 10000;
    int wo_count = 12;
    int wi_count = 80;
    int rrss_oversampling = 4;
    int rrss_path_count = 0; // Not an option, determined by path fraction.
    double rrss_path_frac = 0.8;
    bool restart = false;
    int thread_count = 0;
    std::string ifs_filename;
    std::string ofs_filename;

    // -s/--seed
    opt_parser.on_option("-s", "--seed", 1,
    [&](char** argv) {
        try {
            seed = std::stoi(argv[0]);
        }
        catch (const std::exception&) {
            throw 
                std::runtime_error(
                std::string("-s/--seed expects 1 integer")
                    .append("(can't parse ").append(argv[0])
                    .append(")"));
        }
    })
    << "Specify seed. By default, 0.\n";

    // -p/--path-count
    opt_parser.on_option("-p", "--path-count", 1,
    [&](char** argv) {
        try {
            path_count = std::stoi(argv[0]);
            if (!(path_count > 0)) {
                throw std::exception();
            }
        }
        catch (const std::exception&) {
            throw
                std::runtime_error(
                std::string("-p/--path-count expects 1 positive integer ")
                    .append("(can't parse ").append(argv[0])
                    .append(")"));
        }
    })
    << "Specify path count. By default, 10000.\n"
       "This is the number of paths to use to estimate the emergent\n"
       "BRDF/BSDF for each outgoing direction.\n";

    // -wo/--wo-count
    opt_parser.on_option("-wo", "--wo-count", 1,
    [&](char** argv) {
        try {
            wo_count = std::stoi(argv[0]); // This may throw.
            if (!(wo_count >= 4 &&
                  wo_count <= 1024)) {
                // Trigger catch block manually.
                throw std::exception();
            }
        }
        catch (const std::exception&) {
            throw
                std::runtime_error(
                std::string("-wo/--wo-count expects 1 integer in [4, 1024] ")
                    .append("(can't parse ").append(argv[0])
                    .append(")"));
        }
    })
    << "Specify outgoing direction count. By default, 12.\n"
       "This is the number of outgoing directions in the upper hemisphere,\n"
       "uniformly distributed in zenith. As the emergent BRDF/BSDF must be\n"
       "isotropic, the implementation does not sample in azimuth.\n";

    // -wi/--wi-count
    opt_parser.on_option("-wi", "--wi-count", 1,
    [&](char** argv) {
        try {
            wi_count = std::stoi(argv[0]); // This may throw.
            if (!(wi_count >= 16 &&
                  wi_count <= 1024)) {
                // Trigger catch block manually.
                throw std::exception();
            }
        }
        catch (const std::exception&) {
            throw
                std::runtime_error(
                std::string("-wi/--wi-count expects 1 integer in [16, 1024] ")
                    .append("(can't parse ").append(argv[0])
                    .append(")"));
        }
    })
    << "Specify incident direction count. By default, 80.\n"
       "This is the number of incident directions per outgoing direction,\n"
       "sampled via a Redundancy-Reduced Sample Set (RRSS) to match the\n"
       "scattering distribution of the emergent BRDF/BSDF.\n";

    // -rx/--rrss-oversampling
    opt_parser.on_option("-rx", "--rrss-oversampling", 1,
    [&](char** argv) {
        try {
            rrss_oversampling = std::stoi(argv[0]); // This may throw.
            if (!(rrss_oversampling >= 2 &&
                  rrss_oversampling <= 8)) {
                // Trigger catch block manually.
                throw std::exception();
            }
        }
        catch (const std::exception&) {
            throw
                std::runtime_error(
                std::string("-rx/--rrss-oversampling expects 1 "
                            "integer in [2, 8] ")
                    .append("(can't parse ").append(argv[0])
                    .append(")"));
        }
    })
    << "Specify RRSS oversampling multiplier. By default, 4.\n"
       "This is the multiplier on number of incident directions to use\n"
       "when forming the Redundancy-Reduced Sample Set (RRSS) of incident\n"
       "directions for each outgoing direction.\n";

    // -rp/--rrss-path-frac
    opt_parser.on_option("-rp", "--rrss-path-frac", 1,
    [&](char** argv) {
        try {
            rrss_path_frac = std::stod(argv[0]); // This may throw.
            if (!(rrss_path_frac > 0 &&
                  rrss_path_frac <= 1)) {
                // Trigger catch block manually.
                throw std::exception();
            }
        }
        catch (const std::exception&) {
            throw
                std::runtime_error(
                std::string("-rp/--rrss-path-frac expects 1 float in (0, 1] ")
                    .append("(can't parse ").append(argv[0])
                    .append(")"));
        }
    })
    << "Specify RRSS path fraction. By default, 0.8.\n"
       "This is the number of paths to use to estimate the emergent\n"
       "BRDF/BSDF when forming the Redundancy-Reduced Sample Set (RRSS) of\n"
       "incident directions for each outgoing direction.\n";

    // -R/--restart
    opt_parser.on_option("-R", "--restart", 0,
    [&](char**) {
        restart = true;
    })
    << "Ignore LSS file if it exists, effectively restarting simulation.\n";

    // -t/--thread-count
    opt_parser.on_option("-t", "--thread-count", 1,
    [&](char** argv) {
        try {
            thread_count = std::stoi(argv[0]); // This may throw.
            if (!(thread_count >= 1 &&
                  thread_count <= 16)) {
                // Trigger catch block manually.
                throw std::exception();
            }
        }
        catch (const std::exception&) {
            throw
                std::runtime_error(
                std::string("-t/--thread-count expects 1 integer in ")
                    .append("[1, 16] (can't parse ").append(argv[0])
                    .append(")"));
        }
    })
    << "Specify number of threads. By default, all available threads.\n";

    // -o/--output
    opt_parser.on_option("-o", "--output", 1,
    [&](char** argv) {
        ofs_filename = argv[0];
    })
    << "Specify output filename. By default, formed by adding \".raw\"\n"
       "to input filename. If input filename is \"-\" to indicate standard\n"
       "input, then is \"Layered.raw\" by default.\n";

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
    else if (ofs_filename.empty()) {
        if (ifs_filename == "-") {
            ofs_filename = "Layered.raw";
        }
        else {
            ofs_filename = ifs_filename + ".raw";
        }
    }

    // RRSS path count.
    rrss_path_count = int(rrss_path_frac * path_count);
    if (rrss_path_count < 1) {
        rrss_path_count = 1;
    }

    // Hardware concurrency unavailable?
    if (std::thread::hardware_concurrency() == 0) {

        // Default to 4.
        if (thread_count == 0) {
            thread_count = 4;
        }
    }
    else {

        // Default to hardware concurrency.
        if (thread_count == 0) {
            thread_count = int(std::thread::hardware_concurrency());
        }

        // Clamp.
        if (thread_count > int(std::thread::hardware_concurrency())) {
            thread_count = int(std::thread::hardware_concurrency());
        }
    }

    ls::LayeredAssembly layered_assembly;
    try {
        // Use standard input?
        if (ifs_filename == "-") {
            std::cout << "Initializing from std::cin...\n";
            std::cout.flush();

            // Initialize.
            layered_assembly.init(std::cin);

            std::cout << "Done.\n\n";
            std::cout.flush();
        }
        else {
            std::cout << "Initializing from " << ifs_filename << "...\n";
            std::cout.flush();

            // Open filestream.
            std::ifstream ifs(ifs_filename);
            if (!ifs.is_open()) {
                throw std::runtime_error("can't open file");
            }

            // Initialize.
            layered_assembly.init(ifs);

            std::cout << "Done.\n\n";
            std::cout.flush();
        }
    }
    catch (const std::exception& exception) {
        std::cerr << "Unhandled exception!\n";
        std::cerr << "exception.what(): " << exception.what() << "\n";
        std::exit(EXIT_FAILURE);
    }

    // File data.
    ls::FileData file_data;

    std::string lss_filename;
    if (ifs_filename != "-") {
        lss_filename = ifs_filename + ".lss";
    }

    // Read in-progress simulation if available.
    bool has_lss = false;
    if (!lss_filename.empty() && !restart) {
        std::ifstream ifs(
                lss_filename,
                std::ios_base::in |
                std::ios_base::binary);
        if (ifs.good()) {
            has_lss = true;
            std::cout << "Reading in-progress simulation " << lss_filename;
            std::cout << "...\n";
            std::cout.flush();

            try {
                // Read LSQT Slice binary format.
                file_data.readLss(ifs);

                std::cout << "Done.\n\n";
                std::cout.flush();
            }
            catch (const std::exception& exception) {
                std::cerr << "Unhandled exception!\n";
                std::cerr << "exception.what(): " << exception.what() << "\n";
                std::cout << "Error, re-initializing simulation...\n\n";
                std::cout.flush();
                has_lss = false;
            }
        }
    }

    // No simulation read?
    if (!has_lss) {

        // Initialize.
        file_data.basicInit(wo_count, wi_count, seed);

#if 0
        // Check RRSS path count.
        if (path_count < rrss_path_count) {
            path_count = rrss_path_count;
            std::cerr << "Warning: path count less than RRSS path count\n";
        }
#endif
    }
    else {

        // Update count from file data.
        wo_count = file_data.slices.size();

        // Note.
        std::cout << "Note: using in-progress simulation, so the following\n";
        std::cout << "command line args are ignored:\n";
        std::cout << "    -s/--seed\n";
        std::cout << "    -wo/--wo-count\n";
        std::cout << "    -wi/--wi-count\n";
        std::cout << "    -rx/--rrss-oversampling\n";
        std::cout << "    -rp/--rrss-path-frac\n\n";
        std::cout.flush();
    }

    {
        // Thread pool.
        pr::thread_pool pool(thread_count);

        // Progress bar.
        ls::ProgressBar progress_bar(path_count * wo_count);

        // Task.
        auto task = [&](int wo_index) {
            auto& slice = file_data.slices[wo_index];
            try {
                if (has_lss) {
                    // Compute BSDF averages.
                    slice.computeBsdfAverages(
                            layered_assembly,
                            progress_bar,
                            path_count);
                }
                else {
                    // Compute RRSS incident directions first.
                    slice.computeIncidentDirs(
                            layered_assembly,
                            progress_bar,
                            rrss_oversampling,
                            rrss_path_count);

                    // Compute BSDF averages.
                    slice.computeBsdfAverages(
                            layered_assembly,
                            progress_bar,
                            path_count - rrss_path_count);
                }
            }
            catch (const std::exception& exception) {
                std::cerr << 
                    std::string(
                    "Unhandled exception on thread!\n"
                    "exception.what(): ")
                        .append(exception.what()).append("\n");
            }
        };

        // Submit tasks.
        std::vector<std::future<void>> task_futures;
        task_futures.reserve(wo_count);
        for (int wo_index = 0;
                 wo_index < wo_count; wo_index++) {
            task_futures.emplace_back(pool.submit(task, wo_index));
        }

        // Wait.
        for (int wo_index = 0;
                 wo_index < wo_count; wo_index++) {
            task_futures[wo_index].wait();
        }

        // Finish progress bar.
        progress_bar.finish();
    }

    if (!lss_filename.empty()) {
        std::cout << "Writing " << lss_filename << "...\n";
        std::cout.flush();

        // Output filestream.
        std::ofstream ofs(
                lss_filename, 
                std::ios_base::out |
                std::ios_base::binary);
        if (!ofs.is_open()) {
            std::cerr << "Error opening " << lss_filename << "...\n\n";

            // No renaming.
        }
        else {
            try {
                // Output LSQT Slice binary format.
                file_data.writeLss(ofs);

                std::cout << "Done.\n\n";
                std::cout.flush();
            }
            catch (const std::exception& exception) {
                std::cerr << "Unhandled exception!\n";
                std::cerr << "exception.what(): " << exception.what() << "\n";
                std::cout << "Error writing " << lss_filename << "...\n\n";
                std::cout.flush();
            }
        }
    }

    {
        std::cout << "Writing " << ofs_filename << "...\n";
        std::cout.flush();

        // Output filestream.
        std::ofstream ofs(ofs_filename);
        while (!ofs.is_open()) {
            std::cerr << "Error opening " << ofs_filename << "...\n";
            if (ifs_filename == "-") {
                std::exit(EXIT_FAILURE);
            }
            else {
                std::cerr << "Enter alternative filename: ";
                std::cin >> ofs_filename;
                ofs.open(ofs_filename);
            }
        }

        // Output SQT RAW.
        file_data.writeSqtRaw(ofs);

        std::cout << "Done.\n\n";
        std::cout.flush();
    }

    return EXIT_SUCCESS;
}
