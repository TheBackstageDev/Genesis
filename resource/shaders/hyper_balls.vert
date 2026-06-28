#version 330 core

out vec2 v_uv;

void main() 
{
    vec2 quad[4] = vec2[4](
        vec2(-1.0, -1.0),
        vec2( 1.0, -1.0),
        vec2(-1.0,  1.0),
        vec2( 1.0,  1.0)
    );

    v_uv = quad[gl_VertexID];
    gl_Position = vec4(v_uv, 0.0, 1.0);
}
