#version 330 core

layout(location = 0) in vec3  v_center;
layout(location = 1) in float s_radius;
layout(location = 2) in vec4  color;

uniform mat4 u_view;
uniform mat4 u_projection;
uniform vec3 u_campos;

out vec3  v_center_view;
out float v_radius;
out vec4  v_color;
out vec2  v_uv;

void main()
{
    v_radius = s_radius;
    v_color = color;

    vec4 center_view = u_view * vec4(v_center, 1.0);
    v_center_view = center_view.xyz;
    
    vec2 uv = vec2(((gl_VertexID << 1) & 2) - 1.0, (gl_VertexID & 2) - 1.0);

    v_uv = uv;

    gl_Position = u_projection * (center_view + vec4(uv, 0.0, 0.0));
}