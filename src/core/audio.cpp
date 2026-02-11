
#define MINIAUDIO_IMPLEMENTATION
#include "audio.hpp"

namespace core
{
    AudioEngine::AudioEngine()
    {
        if (ma_engine_init(nullptr, &engine) != MA_SUCCESS)
        {
            throw std::runtime_error("Failed to initialize MiniAudio engine");
        }
    }

    AudioEngine::~AudioEngine()
    {
        stopAll();

        for (auto &kv : sounds)
            ma_sound_uninit(kv.second.get());
        
        ma_engine_uninit(&engine);
    }
}; // namespace core