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
#include <preform/misc_float.hpp>
#include <preform/thread_pool.hpp>
#include <preform/option_parser.hpp>
#include <preform/bash_format.hpp>
#include <layered-sqt/common.hpp>
#include <layered-sqt/medium.hpp>
#include <layered-sqt/layer.hpp>
#include <layered-sqt/layered_assembly.hpp>
#include <layered-sqt/rrss.hpp>

int main(int argc, char** argv)
{
    pr::option_parser opt_parser("[OPTIONS] filename");

    // Options.
    int seed = 0;
    int path_count = 10000;
    int wo_count = 12;
    int wi_count = 80;
    int rrss_oversampling = 4;
    int rrss_path_count = 2000;
    int thread_count = 0;
    std::string ifs_filename;
    std::string ofs_filename = "Layered.raw";

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

    // --wo-count
    opt_parser.on_option(nullptr, "--wo-count", 1,
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
                std::string("--wo-count expects 1 integer in [4, 1024] ")
                    .append("(can't parse ").append(argv[0])
                    .append(")"));
        }
    })
    << "Specify outgoing direction count. By default, 12.\n"
       "This is the number of outgoing directions in the upper hemisphere,\n"
       "uniformly distributed in zenith. As the emergent BRDF/BSDF must be\n"
       "isotropic, the implementation does not sample in azimuth.\n";

    // --wi-count
    opt_parser.on_option(nullptr, "--wi-count", 1,
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
                std::string("--wi-count expects 1 integer in [16, 1024] ")
                    .append("(can't parse ").append(argv[0])
                    .append(")"));
        }
    })
    << "Specify incident direction count. By default, 80.\n"
       "This is the number of incident directions per outgoing direction,\n"
       "sampled via a Redundancy-Reduced Sample Set (RRSS) to match the\n"
       "scattering distribution of the emergent BRDF/BSDF.\n";

    // --rrss-oversampling
    opt_parser.on_option(nullptr, "--rrss-oversampling", 1,
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
                std::string("--rrss-oversampling expects 1 integer in [2, 8] ")
                    .append("(can't parse ").append(argv[0])
                    .append(")"));
        }
    })
    << "Specify RRSS oversampling multiplier. By default, 4.\n"
       "This is the multiplier on number of incident directions to use\n"
       "when forming the Redundancy-Reduced Sample Set (RRSS) of incident\n"
       "directions for each outgoing direction.\n";

    // --rrss-path-count
    opt_parser.on_option(nullptr, "--rrss-path-count", 1,
    [&](char** argv) {
        try {
            rrss_path_count = std::stoi(argv[0]); // This may throw.
            if (!(rrss_path_count > 0)) {
                // Trigger catch block manually.
                throw std::exception();
            }
        }
        catch (const std::exception&) {
            throw
                std::runtime_error(
                std::string("--rrss-path-count expects 1 positive integer ")
                    .append("(can't parse ").append(argv[0])
                    .append(")"));
        }
    })
    << "Specify RRSS path count. By default, 2000.\n"
       "This is the number of paths to use to estimate the emergent\n"
       "BRDF/BSDF when forming the Redundancy-Reduced Sample Set (RRSS) of\n"
       "incident directions for each outgoing direction.\n";

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
    << "Specify output filename. By default, Layered.raw.\n";

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

    if (path_count < rrss_path_count) {
        path_count = rrss_path_count;
        std::cerr << "Warning: path count less than RRSS path count\n";
    }

    {
        using namespace ls;

        // Raw outgoing angles.
        Float* raw_thetao = new Float[wo_count];

        // Raw outgoing directions.
        Vec3<Float>* raw_wo = new Vec3<Float>[wo_count];

        // Raw incident directions.
        Vec3<Float>* raw_wi = new Vec3<Float>[wo_count * wi_count];

        // Raw BSDFs.
        Float* raw_f = new Float[wo_count * wi_count];

        for (int wo_index = 0; 
                 wo_index < wo_count; wo_index++) {

            // Initialize outgoing angle.
            raw_thetao[wo_index] = 
                    wo_index / Float(wo_count) * 
                    pr::numeric_constants<Float>::M_pi_2();

            // Initialize outgoing direction.
            raw_wo[wo_index] = {
                pr::sin(raw_thetao[wo_index]), 0,
                pr::cos(raw_thetao[wo_index])
            };
        }

        // Print mutex.
        std::mutex print_mutex;
        std::int64_t paths_complete = 0;
        std::int64_t paths_todo = 
                std::int64_t(wo_count) * 
                std::int64_t(path_count);

        // Print initial 0% progress bar.
        std::cout << "Computing...\n";
        std::cout << pr::terminal_progress_bar{0.0};
        std::cout.flush();

        // Task executed on thread.
        auto task = [&](int wo_index) {

            // Permuted-congruential generator.
            Pcg32 pcg(seed, wo_index);

            // Local outgoing direction.
            Vec3<Float>& wo = raw_wo[wo_index];

            // Local incident directions.
            Vec3<Float>* wi = raw_wi + wo_index * wi_count;

            // Local BSDFs.
            Float* f = raw_f + wo_index * wi_count;

            {
                // RRSS incident direction count.
                int rrss_wi_count = 
                    rrss_oversampling * wi_count;

                // RRSS incident directions.
                Vec3<Float>* rrss_wi = new Vec3<Float>[rrss_wi_count];

                // Initialize incident directions.
                for (int rrss_wi_index = 0;
                         rrss_wi_index < rrss_wi_count;
                         rrss_wi_index++) {
                
                    rrss_wi[rrss_wi_index] = 
                    layered_assembly.randomScatterDirection(pcg, wo); 
                }

                // RRSS BSDFs.
                Float* rrss_f = new Float[rrss_wi_count];

                // RRSS BSDF-PDFs.
                Float* rrss_f_pdf = new Float[rrss_wi_count];

                // Average BSDFs/BSDF-PDFs in blocks of 512.
                for (int curr_path_count = 0;
                         curr_path_count < rrss_path_count;
                         curr_path_count += 512) {

                    int next_path_count = curr_path_count + 512;
                    if (next_path_count > rrss_path_count) {
                        next_path_count = rrss_path_count;
                    }

                    // Update average BSDFs/BSDF-PDFs.
                    layered_assembly.computeAverage(
                            curr_path_count,
                            next_path_count,
                            pcg, wo,
                            rrss_wi,
                            rrss_wi_count,
                            rrss_f, rrss_f_pdf);

                    {
                        // Update terminal progress bar.
                        std::unique_lock<std::mutex> lock(print_mutex);
                        paths_complete += 
                            next_path_count - 
                            curr_path_count;
                        std::cout << '\r';
                        std::cout << pr::terminal_progress_bar{
                            double(paths_complete) / 
                            double(paths_todo)};
                        std::cout.flush();
                    }
                }

                // RRSS samples.
                std::vector<Rrss::Sample> rrss_samples;
                rrss_samples.reserve(rrss_wi_count);

                // RRSS BSDF integral.
                Float rrss_f_int = 0;

                for (int rrss_wi_index = 0;
                         rrss_wi_index < rrss_wi_count; 
                         rrss_wi_index++) {

                    // Push sample.
                    rrss_samples.push_back(
                    Rrss::Sample{
                        rrss_wi[rrss_wi_index],
                        rrss_f [rrss_wi_index] // BSDF as initial PDF.
                    });

                    // BSDF integral term.
                    Float f_int = 
                        rrss_f[rrss_wi_index] /
                        rrss_f_pdf[rrss_wi_index];

                    // BSDF integral term okay?
                    if (pr::isfinite(f_int)) {

                        // Update BSDF integral estimate.
                        rrss_f_int =
                        rrss_f_int + (f_int - rrss_f_int) / 
                                (rrss_wi_index + 1);
                    }
                }

                // BSDF integral okay?
                if (rrss_f_int > 0 &&
                    pr::isfinite(rrss_f_int)) {

                    // Normalize.
                    for (Rrss::Sample& rrss_sample : rrss_samples) {
                        rrss_sample.pdf /= rrss_f_int;
                    }
                }
                else {

                    {
                        // Warn.
                        std::unique_lock<std::mutex> lock(print_mutex);
                        std::cerr << 
                            "Warning: RRSS failed to compute valid BSDF "
                            "integral, defaulting to BSDF-PDF\n";
                    }

                    // Default to BSDF-PDF.
                    int rrss_wi_index = 0;
                    for (Rrss::Sample& rrss_sample : rrss_samples) {
                        rrss_sample.pdf = rrss_f_pdf[rrss_wi_index++];
                    }
                }

                // Redundancy-reduced sample set.
                Rrss rrss(rrss_samples);
                if (!rrss.disableUntil(wi_count)) {
                    // Unreachable?
                }

                int wi_index = 0;
                for (int rrss_wi_index = 0;
                         rrss_wi_index < rrss_wi_count; 
                         rrss_wi_index++) {

                    // Is sample enabled?
                    if (rrss.samples()[rrss_wi_index].is_enabled) {

                        // Use incident direction and BSDF.
                        wi[wi_index] = rrss_wi[rrss_wi_index];
                        f [wi_index] = rrss_f [rrss_wi_index];
                        wi_index++;
                    }
                }

                // Delete RRSS BSDF-PDFs.
                delete[] rrss_f_pdf;

                // Delete RRSS BSDFs.
                delete[] rrss_f;

                // Delete RRSS incident directions.
                delete[] rrss_wi;
            }

            // Average BSDFs in blocks of 4096.
            for (int curr_path_count = rrss_path_count;
                     curr_path_count < path_count;
                     curr_path_count += 4096) {

                int next_path_count = curr_path_count + 4096;
                if (next_path_count > path_count) {
                    next_path_count = path_count;
                }

                // Update average BSDFs.
                layered_assembly.computeAverage(
                            curr_path_count,
                            next_path_count,
                            pcg, wo, wi, wi_count, f, nullptr);

                {
                    // Update terminal progress bar.
                    std::unique_lock<std::mutex> lock(print_mutex);
                    paths_complete += 
                        next_path_count - 
                        curr_path_count;
                    std::cout << '\r';
                    std::cout << pr::terminal_progress_bar{
                        double(paths_complete) / 
                        double(paths_todo)};
                    std::cout.flush();
                }
            }

            // SQT expects BSDF without cosine weighting...
            for (int wi_index = 0;
                     wi_index < wi_count; wi_index++) {
                f[wi_index] /= pr::abs(wi[wi_index][2]);
            }
        };

        {
            // Thread pool.
            pr::thread_pool pool(thread_count);

            // Submit tasks.
            std::vector<std::future<void>> task_futures;
            task_futures.reserve(wo_count);
            for (int wo_index = 0;
                     wo_index < wo_count; wo_index++) {
                task_futures.emplace_back(
                        pool.submit(task, wo_index));
            }

            // Wait.
            for (int wo_index = 0;
                     wo_index < wo_count; wo_index++) {
                task_futures[wo_index].wait();
            }   

            // Newlines after terminal progress bar. 
            std::cout << "\n\n";
            std::cout.flush();
        }


        bool is_bsdf = false;
        for (int wi_index = 0; 
                 wi_index < wo_count * wi_count; wi_index++) {
            if (raw_wi[wi_index][2] < 0) {
                is_bsdf = true;
                break;
            }
        }

        if (is_bsdf == false && layered_assembly.isTransmissive()) {
            std::cerr << "Warning: no incident directions sampled in "
                         "lower hemisphere, though layered assembly is "
                         "potentially transmissive\n\n";
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
            ofs << "RAWB" << (is_bsdf ? 'S' : 'H');
            ofs << "10A Layered-SQT\n";
            ofs << "1 0.5\n";
            ofs << wo_count;
            for (int wo_index = 0; 
                     wo_index < wo_count; wo_index++) {
                ofs << ' ';
                ofs << raw_thetao[wo_index];
            }
            ofs << '\n';
            for (int wo_index = 0; 
                     wo_index < wo_count; wo_index++) {
                for (int wi_index = 0; 
                         wi_index < wi_count; wi_index++) {
                    ofs << wo_index << ' ';
                    ofs << raw_wi[wo_index * wi_count + wi_index][0] << ' ';
                    ofs << raw_wi[wo_index * wi_count + wi_index][1] << ' ';
                    ofs << raw_wi[wo_index * wi_count + wi_index][2] << ' ';
                    ofs << raw_f[wo_index * wi_count + wi_index] << '\n';
                }
            }
        }

        // Delete raw BSDFs.
        delete[] raw_f;

        // Delete raw incident directions.
        delete[] raw_wi;

        // Delete raw outgoing directions.
        delete[] raw_wo;

        // Delete raw outgoing angles.
        delete[] raw_thetao;
    }

    return EXIT_SUCCESS;
}
