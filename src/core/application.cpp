#include "application.hpp"

namespace core
{
    application::application(int32_t height, int32_t width, const std::string name)
    {
        window = std::make_unique<core::window_t>(width, height, name);
    }

    application::~application()
    {
    
    }

    void application::run()
    {
        
    }
} // namespace core
