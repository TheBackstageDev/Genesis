#include "shader.hpp"

#include <fstream>
#include <iostream>

namespace core
{
    std::string loadShaderSource(const std::filesystem::path &path)
    {
        std::ifstream file(path, std::ios::ate | std::ios::binary);
        if (!file.is_open())
        {
            throw std::runtime_error("Cannot open shader: " + path.string());
        }

        auto size = file.tellg();
        std::string content(size, '\0');
        file.seekg(0);
        file.read(content.data(), size);

        if (file.fail())
        {
            throw std::runtime_error("Failed to read shader: " + path.string());
        }

        return content;
    }

    glShader::glShader(GLenum type, const std::filesystem::path shaderPath)
    {
        m_id = glCreateShader(type);
        if (m_id == 0)
        {
            throw std::runtime_error("glCreateShader failed for type " + std::to_string(type));
        }

        std::string shader_content = loadShaderSource(shaderPath);
        const char* sources[] = { shader_content.c_str() };
        glShaderSource(m_id, 1, sources, nullptr);

        glCompileShader(m_id);

        GLint success = 0;
        glGetShaderiv(m_id, GL_COMPILE_STATUS, &success);
        if (success == GL_FALSE)
        {
            GLint logLength = 0;
            glGetShaderiv(m_id, GL_INFO_LOG_LENGTH, &logLength);

            std::vector<char> log(logLength);
            glGetShaderInfoLog(m_id, logLength, nullptr, log.data());

            throw std::runtime_error(
                "Shader compilation failed [" + shaderPath.string() + "]:\n" +
                std::string(log.data())
            );
        }
    }

    glShader::~glShader()
    {
        if (m_id != 0) glDeleteShader(m_id);
    }

    glProgram::glProgram(const std::vector<glShader>& shaders) 
    {
        createProgram();

        for (const glShader& shader : shaders)
            attach(shader);

        link();
    }

    glProgram::glProgram(const glShader& vertex, const glShader& fragment)
    {
        createProgram();

        attach(vertex);
        attach(fragment);
        link();
    }

    glProgram::glProgram(const glShader& computeShader) 
    {
        createProgram();
        attach(computeShader);
        link();
    }

    glProgram::~glProgram()
    {
        glDeleteProgram(m_id);
    }

    glProgram::glProgram(glProgram&& other) noexcept
        : m_id(other.m_id)
    {
        other.m_id = 0;
    }

    glProgram& glProgram::operator=(glProgram&& other) noexcept
    {
        if (this != &other)
        {
            if (m_id != 0)
            {
                glDeleteProgram(m_id);
            }

            m_id = other.m_id;
            other.m_id = 0;
        }
        return *this;
    }

    void glProgram::link()
    {
        glLinkProgram(m_id);

        GLint success = 0;
        glGetProgramiv(m_id, GL_LINK_STATUS, &success);
        if (!success)
        {
            char info[1024]{};
            glGetProgramInfoLog(m_id, 1024, nullptr, info);
            std::string msg = "[Shader link failed]: " + std::string(info);
            std::cerr << msg << "\n";
            throw std::runtime_error(msg);
        }

        glValidateProgram(m_id);
        glGetProgramiv(m_id, GL_VALIDATE_STATUS, &success);
        if (!success)
        {
            char info[1024]{};
            glGetProgramInfoLog(m_id, 1024, nullptr, info);
            std::cerr << "Program validation warning: " << info << "\n";
        }
    }

    GLint glProgram::getUniformLocation(const std::string& name) const
    {
        return glGetUniformLocation(m_id, name.c_str());
    }

    GLuint glProgram::getUniformBlockIndex(const std::string& name) const
    {
        return glGetUniformBlockIndex(m_id, name.c_str());
    }

    void glProgram::uniformBlockBinding(const std::string& name, GLuint bindingPoint) const
    {
        GLint blockIndex = getUniformBlockIndex(name);
        if (blockIndex == GL_INVALID_INDEX) 
        {
            std::cerr << "Warning: Uniform block '" << name 
                    << "' not found in program " << m_id << std::endl;
            return;
        }

        glUniformBlockBinding(m_id, static_cast<GLuint>(blockIndex), bindingPoint);
    }
}; // namespace core