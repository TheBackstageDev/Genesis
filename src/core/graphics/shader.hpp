#pragma once

#define GLM_FORCE_RADIANS
#define GLM_FORCE_INTRINSICS
#define GLM_FORCE_ALIGNED_GENTYPES

#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <glad/glad.h>
#include <filesystem>

#include <vector>

namespace core
{
   class glShader
    {
    public:
        glShader(GLenum type, const std::filesystem::path shaderPath);
        ~glShader();

        glShader(const glShader&) = delete;
        glShader& operator=(const glShader&) = delete;

        glShader(glShader&& other) noexcept;
        glShader& operator=(glShader&& other) noexcept;

        GLuint id() const { return m_id; }
        GLenum type() const { return m_type; }
    private:
        GLuint m_id;
        GLenum m_type;
    };

    class glProgram
    {
    public:

        glProgram() = default;
        glProgram(const std::vector<glShader>& shaders);
        glProgram(const glShader& vertex, const glShader& fragment);
        explicit glProgram(const glShader& computeShader);

        ~glProgram();

        glProgram(const glProgram&) = delete;
        glProgram& operator=(const glProgram&) = delete;

        glProgram(glProgram&& other) noexcept;
        glProgram& operator=(glProgram&& other) noexcept;
        
        GLuint id() const { return m_id; }

        void use() const { glUseProgram(m_id); }

        void setUniform(const std::string& name, const float v) const { glUniform1f(getUniformLocation(name), v); }
        void setUniform(const std::string& name, const int32_t v) const { glUniform1i(getUniformLocation(name), v); }
        void setUniform(const std::string& name, const uint32_t v) const { glUniform1ui(getUniformLocation(name), v); }
        void setUniform(const std::string& name, const glm::vec3& v) const { glUniform3fv(getUniformLocation(name), 1, glm::value_ptr(v)); }
        void setUniform(const std::string& name, const glm::vec4& v) const { glUniform4fv(getUniformLocation(name), 1, glm::value_ptr(v)); }
        void setUniform(const std::string& name, const glm::mat4& m) const { glUniformMatrix4fv(getUniformLocation(name), 1, GL_FALSE, glm::value_ptr(m)); }

        GLint getUniformLocation(const std::string& name) const;
        GLuint getUniformBlockIndex(const std::string& name) const;
        void uniformBlockBinding(const std::string& name, GLuint bindingPoint) const;

    private:
        void createProgram()
        {
            m_id = 0;
            m_id = glCreateProgram();

            if (m_id == 0) 
            {
                throw std::runtime_error("Failed to create program");
            }
        }

        void attach(const glShader& shader) 
        {
            glAttachShader(m_id, shader.id());
        }

        void link();

        GLuint m_id;
    };
} // namespace core
