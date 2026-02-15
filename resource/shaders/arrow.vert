#version 330 core

layout(location = 0) in vec4  a_centerA;
layout(location = 1) in vec4  a_centerB;
layout(location = 2) in vec4  a_color;
layout(location = 3) in float a_radius;

uniform mat4 u_view;
uniform mat4 u_proj;

out vec3  v_start;      // view-space start of cylinder
out vec3  v_end;        // view-space end of cylinder
out vec4  v_color;
out float v_radius;
out vec3  v_corner;

void main()
{
    vec4 viewA = u_view * a_centerA;
    vec4 viewB = u_view * a_centerB;
    v_start  = viewA.xyz;
    v_end    = viewB.xyz;
    v_color = a_color;
    v_radius    = a_radius;

    vec3 axis     = v_end - v_start;
    float axisLen = length(axis);

    if (axisLen < 1e-7)
    {
        gl_Position = vec4(0.0);
        v_corner = vec3(0.0);
        return;
    }

    vec3 midView = (v_start + v_end) * 0.5;
    vec3 quadCenter = midView * 0.98;

    vec2 uv = vec2(
        ((gl_VertexID << 1) & 2) - 1.0,
        (gl_VertexID & 2) - 1.0
    );

    const float expansion = 25.0;
    vec3 quadCorner = quadCenter + vec3(uv * v_radius * expansion, 0.0);
    v_corner = quadCorner;

    gl_Position = u_proj * vec4(quadCorner, 1.0);
}