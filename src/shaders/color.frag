#version 330 core

in vec2  v_center;
in float v_radius;
in vec4  v_color;
in vec2  v_uv;

uniform vec3 u_lightDir = normalize(vec3(0.4, 0.8, 1.2));
uniform mat4 u_view;
uniform mat4 u_proj;

out vec4 fragColor;

vec2 sphIntersect(vec3 ro, vec3 rd, vec3 ce, float ra)
{
    vec3 oc = ro - ce;
    float b = dot( oc, rd );
    float c = dot( oc, oc ) - ra*ra;
    float h = b*b - c;
    if( h<0.0 ) return vec2(-1.0);
    h = sqrt( h );
    return vec2( -b-h, -b+h );
}

vec3 sphNormal(vec3 pos, vec4 sph)
{
    return normalize(pos-sph.xyz);
}

void main()
{
    float dist = length(v_uv);
    if (dist > 1.0) discard;

    fragColor = v_color;
}