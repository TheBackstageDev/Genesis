#pragma once

#include <glad/glad.h>
#include <cstddef>
#include <stdexcept>
#include <utility>

namespace core
{
    class glBuffer
    {
    public:
        glBuffer() 
            : m_id(0), m_target(GL_SHADER_STORAGE_BUFFER)
        {
            glGenBuffers(1, &m_id);
            if (m_id == 0)
            {
                throw std::runtime_error("Failed to generate buffer");
            }
        }

        glBuffer(GLenum target, const void* data, GLsizeiptr size, GLenum usage = GL_DYNAMIC_DRAW)
            : glBuffer()
        {
            m_target = target;
            upload(data, size, usage);
        }

        ~glBuffer()
        {
            if (m_id != 0) 
            {
                glDeleteBuffers(1, &m_id);
                m_id = 0;
            }
        }

        glBuffer(glBuffer&& other) noexcept
            : m_id(other.m_id), m_target(other.m_target)
        {
            other.m_id = 0;
        }

        glBuffer& operator=(glBuffer&& other) noexcept
        {
            if (this != &other) {
                if (m_id != 0) glDeleteBuffers(1, &m_id);
                m_id     = other.m_id;
                m_target = other.m_target;
                other.m_id = 0;
            }
            return *this;
        }

        glBuffer(const glBuffer&) = delete;
        glBuffer& operator=(const glBuffer&) = delete;

        void bind(GLenum target = 0) const
        {
            glBindBuffer(target ? target : m_target, m_id);
        }

        static void unbind(GLenum target = GL_SHADER_STORAGE_BUFFER)
        {
            glBindBuffer(target, 0);
        }

        void upload(const void* data, GLsizeiptr size, GLenum usage = GL_DYNAMIC_DRAW)
        {
            bind();
            glBufferData(m_target, size, data, usage);
        }

        void update(const void* data, GLsizeiptr size, GLintptr offset = 0)
        {
            bind();
            glBufferSubData(m_target, offset, size, data);
        }

        void* map(GLbitfield access = GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT)
        {
            bind();
            return glMapBufferRange(m_target, 0, getSize(), access);
        }

        void unmap()
        {
            glUnmapBuffer(m_target);
        }

        GLsizeiptr getSize() const
        {
            GLint size = 0;
            bind();
            glGetBufferParameteriv(m_target, GL_BUFFER_SIZE, &size);
            return static_cast<GLsizeiptr>(size);
        }

        GLuint id() const { return m_id; }
        GLenum target() const { return m_target; }

        void bindBase(GLuint binding_point) const
        {
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, binding_point, m_id);
        }
    private:
        GLuint m_id;
        GLenum m_target;
    };
} // namespace core
