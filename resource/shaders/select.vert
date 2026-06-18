#version 330 core
layout(location = 0) in vec3  a_center;
layout(location = 1) in float a_radius;
layout(location = 2) in vec4  a_color;
layout(location = 3) in uint  a_id;

uniform mat4 u_view;
uniform mat4 u_proj;

flat out uvec3 v_id;

void main()
{
    v_id = uvec3(a_id & 0xFFu, (a_id >> 8) & 0xFFu, (a_id >> 16) & 0xFFu);

    vec4 pos = u_proj * u_view * vec4(a_center, 1.0);
    gl_Position = pos;
}