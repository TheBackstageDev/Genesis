#pragma once

#include <SFML/Graphics.hpp>
#include <imgui/imgui.h>

#include <vector>
#include <map>

namespace core
{
    struct videoMetaData
    {
        std::string title;
        std::string description;

        sf::Vector3f box;
        size_t num_atoms;
        size_t num_frames;
    };

    struct frame
    {
        std::vector<sf::Vector3f> positions;
        std::map<size_t, float> temperatures;

        float global_temperature;
        bool isKeyframe = false;
    };

    struct video
    {
        std::vector<frame> frames; // index is defined by order

        videoMetaData metadata;
    };

    class UIHandler
    {
    public:
        void drawVideoControls();

        bool isrewinding() { return rewinding; }
    private:
        bool rewinding = false;
    };
} // namespace core
