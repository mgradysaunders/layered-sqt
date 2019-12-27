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
#include <preform/color.hpp>
#include <preform/option_parser.hpp>
#include <layered-sqt/common.hpp>
#include <layered-sqt/file_data.hpp>
#include <layered-sqt/tri.hpp>

extern "C" {
extern int stbi_write_png(const char*, int, int, int, const void*, int);
} // extern "C"

int main(int argc, char** argv)
{
    pr::option_parser opt_parser("[OPTIONS] filename");

    int image_dim = 512;
    float image_mul = 1.0f;
    std::string ifs_filename;
    std::string ofs_filename = "preview.png";

    // --image-dim
    opt_parser.on_option(nullptr, "--image-dim", 1,
    [&](char** argv) {
        try {
            image_dim = std::stoi(argv[0]);
            if (!(image_dim > 0 && 
                  image_dim <= 1024)) {
                throw std::exception();
            }
        }
        catch (const std::exception&) {
            throw
                std::runtime_error(
                std::string("--image-dim expects 1 integer in [1, 1024] ")
                    .append("(can't parse ").append(argv[0])
                    .append(")"));
        }
    })
    << "Specify image dimension. By default, 512.\n";

    // --image-mul
    opt_parser.on_option(nullptr, "--image-mul", 1,
    [&](char** argv) {
        try {
            image_mul = std::stof(argv[0]);
            if (!(image_mul > 0 && 
                  pr::isfinite(image_mul))) {
                throw std::exception();
            }
        }
        catch (const std::exception&) {
            throw
                std::runtime_error(
                std::string("--image-mul expects 1 postive float ")
                    .append("(can't parse ").append(argv[0])
                    .append(")"));
        }
    })
    << "Specify image brightness multiplier. By default, 1.0.\n";

    // -o/--output
    opt_parser.on_option("-o", "--output", 1,
    [&](char** argv) {
        ofs_filename = argv[0];
    })
    << "Specify output filename. By default, preview.png.\n";

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

    if (ifs_filename.empty()) {
        std::cout << opt_parser << std::endl;
        std::exit(EXIT_SUCCESS);
    }

    ls::FileData file_data;
    try {
        // Reading...
        std::cout << "Reading " << ifs_filename << "...\n";
        std::cout.flush();

        // Open filestream.
        std::ifstream ifs(
                ifs_filename,
                std::ios_base::in |
                std::ios_base::binary);
        if (!ifs.is_open()) {
            throw std::runtime_error("can't open file");
        }

        // Read LSQT Slice binary format.
        file_data.readLss(ifs);

        // Done.
        std::cout << "Done.\n\n";
        std::cout.flush();
    }
    catch (const std::exception& exception) {
        std::cerr << "Unhandled exception!\n";
        std::cerr << "exception.what(): " << exception.what() << "\n";
        std::exit(EXIT_FAILURE);
    }

    // Allocate pixels.
    std::uint8_t* image_pixu8 = 
                new std::uint8_t[image_dim * image_dim];

    {
        // Allocate pixels.
        float* image_pix = new float[image_dim * image_dim];

        // Triangulated file data, render sphere example.
        ls::TriFileData tri;
        tri.init(file_data);
        tri.renderSphereExample(
                image_dim, 
                image_pix);

        // Pack as 8-bit unsigned integers.
        for (int i = 0; i < image_dim; i++)
        for (int j = 0; j < image_dim; j++) {
            image_pixu8[i * image_dim + j] = 
                    pr::pack_uint8(
                    pr::srgbenc_hejl_burgess(
                        image_mul * 
                        image_pix[i * image_dim + j]));
        }

        // Deallocate pixels.
        delete[] image_pix;
    }

    {
        // Writing...
        std::cout << "Writing " << ofs_filename << "...\n";
        std::cout.flush();

#if 0
        // Open filestream.
        std::ofstream ofs(ofs_filename);

        // Write header.
        ofs << "P2\n";
        ofs << image_dim << ' ';
        ofs << image_dim << '\n';
        ofs << "255\n";

        // Write pixels.
        for (int j = 0; j < image_dim; j++)
        for (int i = 0; i < image_dim; i++) {
            ofs << int(image_pixu8[i * image_dim + j]) << ' ';
        }
#else
        // TODO Error handling
        stbi_write_png(
                ofs_filename.c_str(), 
                image_dim, 
                image_dim, 
                1,
                image_pixu8, image_dim);
#endif

        // Done.
        std::cout << "Done.\n\n";
        std::cout.flush();
    }

    delete[] image_pixu8;

    return EXIT_SUCCESS;
}
