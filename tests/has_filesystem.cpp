#include <filesystem>

int main()
{
    std::filesystem::path path = std::filesystem::current_path();
    return 0;
}
