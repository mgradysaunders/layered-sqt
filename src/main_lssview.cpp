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
#include <cctype>
#include <fstream>
#include <preform/color.hpp>
#include <preform/option_parser.hpp>
#include <layered-sqt/common.hpp>
#include <layered-sqt/file_data.hpp>
#include <layered-sqt/tri.hpp>

// Get filename extension in lowercase.
std::string getFilenameExtension(const std::string& filename)
{
    std::size_t pos = filename.find_last_of('.');
    if (pos == std::string::npos) {
        return std::string();
    }
    else {
        std::string ext = filename.substr(pos);
        std::transform(
                ext.begin(), ext.end(), 
                ext.begin(),
                [](unsigned char c) {
                    return std::tolower(c);
                });
        return ext;
    }
}

extern "C" {

typedef void stbi_write_func(void*, void*, int);

extern 
int stbi_write_png_to_func(
    stbi_write_func*, void*, int, int, int, const void*, int);

extern 
int stbi_write_bmp_to_func(
    stbi_write_func*, void*, int, int, int, const void*);

extern 
int stbi_write_tga_to_func(
    stbi_write_func*, void*, int, int, int, const void*);

extern 
int stbi_write_jpg_to_func(
    stbi_write_func*, void*, int, int, int, const void*, int);

void ofsWriter(void* ofs, void* data, int data_size)
{
    static_cast<std::ofstream*>(ofs)->write(
    static_cast<char*>(data), data_size);
}

} // extern "C"

void stbiWriteImage(
        const std::string& ofs_filename,
        int image_sizex, 
        int image_sizey,
        int image_comp,
        const std::uint8_t* image_data)
{
    std::ofstream ofs(
            ofs_filename,
            std::ios_base::out |
            std::ios_base::binary);
    if (!ofs.is_open() || 
        !ofs.good()) {
        throw 
            std::runtime_error(
            std::string(__PRETTY_FUNCTION__)
                .append(": can't open file"));

    }

    int err = 0;
    std::string ext = getFilenameExtension(ofs_filename);
    if (ext == ".png") {
        err =
        stbi_write_png_to_func(
                ofsWriter, 
                static_cast<void*>(&ofs),
                image_sizex,
                image_sizey,
                image_comp,
                image_data,
                image_sizex * image_comp);
    }
    else if (ext == ".bmp") {
        err =
        stbi_write_bmp_to_func(
                ofsWriter, 
                static_cast<void*>(&ofs),
                image_sizex,
                image_sizey,
                image_comp,
                image_data);
    }
    else if (ext == ".tga") {
        err =
        stbi_write_tga_to_func(
                ofsWriter, 
                static_cast<void*>(&ofs),
                image_sizex,
                image_sizey,
                image_comp,
                image_data);
    }
    else if (ext == ".jpg") {
        err =
        stbi_write_jpg_to_func(
                ofsWriter, 
                static_cast<void*>(&ofs),
                image_sizex,
                image_sizey,
                image_comp,
                image_data, 90);
    }
    else {
        throw 
            std::runtime_error(
            std::string(__PRETTY_FUNCTION__)
                .append(": unknown extension"));
    }

    if (err == 0 || !ofs.good()) {
        throw
            std::runtime_error(
            std::string(__PRETTY_FUNCTION__)
                .append(": write failure"));
    }
}

int main(int argc, char** argv)
{
    pr::option_parser opt_parser("[OPTIONS] filename");

    int image_dim = 512;
    float image_mul = 1.0f;
    std::string ifs_filename;
    std::string ofs_filename;

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
        try {
            // Must have 5 or more characters. 
            ofs_filename = argv[0];
            std::string ext = getFilenameExtension(ofs_filename);
            if (ext != ".png" &&
                ext != ".tga" &&
                ext != ".bmp" && 
                ext != ".jpg") {
                throw std::exception();
            }
        }
        catch (const std::exception&) {
            throw
                std::runtime_error(
                std::string("-o/--output expects string ending with ")
                    .append("\".png\", \".tga\", \".bmp\", or \".jpg\" ")
                    .append("(can't parse ").append(argv[0])
                    .append(")"));
        }
    })
    << "Specify output filename. By default, formed by adding \".png\"\n"
       "to input filename. Note, if this filename is specified, it must end\n"
       "with \".png\", \".bmp\", \".tga\", or \".jpg\".\n";

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
    else if (ofs_filename.empty()) {
        ofs_filename = ifs_filename + ".png";
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

        try {
            stbiWriteImage(
                    ofs_filename, 
                    image_dim, 
                    image_dim, 
                    1,
                    image_pixu8);
        }
        catch (const std::exception& exception) {
            std::cerr << "Unhandled exception!\n";
            std::cerr << "exception.what(): " << exception.what() << "\n";
        }

        // Done.
        std::cout << "Done.\n\n";
        std::cout.flush();
    }

    delete[] image_pixu8;

    return EXIT_SUCCESS;
}
