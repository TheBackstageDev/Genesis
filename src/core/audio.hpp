#pragma once

#define MA_ENABLE_MP3
#include <miniaudio.hpp>

#include <string>
#include <iostream>
#include <stdexcept>
#include <memory>
#include <random>

#include <unordered_map>

namespace core
{
    class AudioEngine
    {
    public:
        AudioEngine();
        ~AudioEngine();

        AudioEngine(const AudioEngine &) = delete;
        AudioEngine &operator=(const AudioEngine &) = delete;

        void playSound(const std::string &filePath)
        {
            if (ma_engine_play_sound(&engine, filePath.c_str(), nullptr) != MA_SUCCESS)
            {
                throw std::runtime_error("Failed to play sound: " + filePath);
            }
        }
        void preloadSound(const std::string &name, const std::string &filePath)
        {
            auto sound = std::make_unique<ma_sound>();
            ma_result result = ma_sound_init_from_file(
                &engine,
                filePath.c_str(),
                MA_SOUND_FLAG_DECODE | MA_SOUND_FLAG_ASYNC,
                NULL, NULL,
                sound.get());

            if (result != MA_SUCCESS)
            {
                std::cerr << "[AudioEngine] Failed to preload sound: " << filePath
                          << " (name: " << name << ")"
                          << " | Error code: " << result << "\n";
                throw std::runtime_error("Failed to preload sound: " + filePath);
            }

            if (name.find("Song") != std::string::npos)
            {
                songs.emplace_back(std::move(sound));
                std::cout << "[AudioEngine] Preloaded song: " << filePath << "\n";
            }
            else
            {
                sounds.emplace(name, std::move(sound));
                std::cout << "[AudioEngine] Preloaded sound: " << name
                          << " from " << filePath << "\n";
            }
        }
        void playPreloaded(const std::string &name, float volume = 0.f)
        {
            if (volume == 0.f)
                volume = getGlobalVolume();

            auto it = sounds.find(name);
            if (it == sounds.end())
            {
                std::cerr << " [SOUND ENGINE]: Sound not preloaded! cannot play! \n";
                return;
            }

            ma_sound_set_volume(it->second.get(), volume);
            ma_sound_start(it->second.get());
        }
        void stopSound(const std::string &name)
        {
            auto it = sounds.find(name);
            if (it == sounds.end())
            {
                std::cerr << " [SOUND ENGINE]: Sound not preloaded! cannot stop! \n";
                return;
            }
            ma_sound_stop(it->second.get());
            ma_sound_reset_stop_time(it->second.get());
            ma_sound_seek_to_pcm_frame(it->second.get(), 0);
        }
        
        void stopAll()
        {
            for (auto& sound : sounds)
            {
                ma_sound_stop(sound.second.get());
                ma_sound_reset_stop_time(sound.second.get());
                ma_sound_seek_to_pcm_frame(sound.second.get(), 0);
            }
        }

        void setLooping(const std::string &name, bool loop)
        {
            auto it = sounds.find(name);
            if (it == sounds.end())
            {
                std::cerr << " [SOUND ENGINE]: Sound not preloaded! cannot loop! \n";
                return;
            }
            ma_sound_set_looping(it->second.get(), loop ? MA_TRUE : MA_FALSE);
        }

        void playRandomSong()
        {
            currentSongIndex = rand() % songs.size();

            ma_sound_reset_stop_time(songs[currentSongIndex].get());
            ma_sound_seek_to_pcm_frame(songs[currentSongIndex].get(), 0);
            ma_sound_start(songs[currentSongIndex].get());
        }

        void stopSong()
        {
            if (currentSongIndex != UINT32_MAX)
                ma_sound_stop(songs[currentSongIndex].get());
        }

        bool songFinished()
        {
            return !ma_sound_is_playing(songs[currentSongIndex].get());
        }

        void setGlobalVolume(float volume)
        {
            volume = std::clamp(volume, 0.0f, 1.0f);
            ma_engine_set_volume(&engine, volume);
        }
        
        float getGlobalVolume()
        {
            return ma_engine_get_volume(&engine);
        }

        size_t getSoundCount() { return sounds.size(); }

    private:
        ma_engine engine{};
        std::unordered_map<std::string, std::unique_ptr<ma_sound>> sounds;
        std::vector<std::unique_ptr<ma_sound>> songs;

        uint32_t currentSongIndex = UINT32_MAX;
    };

} // namespace core
