#version 330 core

layout(location = 0) in vec3  a_center;
layout(location = 1) in float a_radius;
layout(location = 2) in vec4  a_color;

uniform mat4 u_view;
uniform mat4 u_proj;

out vec3  v_center;
out float v_radius;
out vec4  v_color;
out vec2  v_uv;

void main()
{
    v_radius = a_radius;
    v_color = a_color;

    if (a_center.z > 0.0 || length(a_center) < a_radius) 
    {
        gl_Position = vec4(0.0);
        return;
    }

    v_center = a_center;
    vec2 uv = vec2(((gl_VertexID << 1) & 2) - 1.0, (gl_VertexID & 2) - 1.0) * 1.2f;

    v_uv = uv;

    vec4 pos_view = vec4(a_center, 1.0) + vec4(uv * a_radius, 0.0, 0.0);
    vec4 pos = u_proj * pos_view;

    gl_Position = pos;
}