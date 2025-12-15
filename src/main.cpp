#include "core/application.hpp"
#include <stdexcept>
#include <iostream>

int main()
{    
    constexpr uint32_t height = 800, width = 1000;

    core::application app{height, width, "Genesis"};

    try 
    {
        app.run();
    }
    catch (std::exception& e)
    {
        std::cerr << "[MAIN] Expection: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}