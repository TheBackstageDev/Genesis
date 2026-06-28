#version 330 core

layout(location = 0) in vec4  a_centerA;   // world-space center A
layout(location = 1) in vec4  a_centerB;   // world-space center B
layout(location = 2) in vec4  a_colorA;
layout(location = 3) in vec4  a_colorB;
layout(location = 4) in float a_radius;
layout(location = 5) in float a_radiusA;
layout(location = 6) in float a_radiusB;

uniform vec2 u_screen;
uniform mat4 u_view;
uniform mat4 u_proj;

out vec3  v_start;       // world-space sphere A center
out vec3  v_end;         // world-space sphere B center
out vec4  v_colorStart;
out vec4  v_colorEnd;
out float v_radiusA;
out float v_radiusB;
out vec3  v_corner;      // world-space quad corner

void main()
{
    v_start      = a_centerA.xyz;
    v_end        = a_centerB.xyz;
    v_colorStart = a_colorB;
    v_colorEnd   = a_colorA;
    v_radiusA    = a_radiusA;
    v_radiusB    = a_radiusB;

    vec3 axis     = v_end - v_start;
    float axisLen = length(axis);

    if (axisLen < 1e-7)
    {
        gl_Position = vec4(0.0);
        v_corner = vec3(0.0);
        return;
    }

    vec3 midWorld = (v_start + v_end) * 0.5;
    vec3 quadCenter = midWorld;

    vec2 uv = vec2(
        ((gl_VertexID << 1) & 2) - 1.0,
        (gl_VertexID & 2) - 1.0
    );

    const float expansion = 40.0;

    float aspect = u_screen.x / u_screen.y;
    vec2 scaledUV = vec2(uv.x * aspect, uv.y);

    vec3 quadCorner = quadCenter + vec3(scaledUV * a_radius * expansion, 0.0);
    v_corner = quadCorner;

    gl_Position = u_proj * u_view * vec4(quadCorner, 1.0);
}
