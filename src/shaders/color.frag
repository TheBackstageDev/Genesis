#version 330 core

in vec3  v_center;
in float v_radius;
in vec4  v_color;
in vec2  v_uv;

uniform vec3 u_lightDir = normalize(vec3(0.4, 0.8, 1.2));
uniform vec3 u_campos;
uniform mat4 u_view;
uniform mat4 u_proj;

out vec4 fragColor;

vec2 sphIntersect(in vec3 ro, in vec3 rd, in vec3 ce, in float ra)
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
    vec3 rayOrig = u_campos;
    vec3 rayDir  = normalize(vec3(v_uv, -1.0));

    vec2 tt = sphIntersect(rayOrig, rayDir, v_center, v_radius * 5.f);

    if (length(v_uv) > 1.0f)
    {
        discard;
    }

    fragColor = v_color;
}