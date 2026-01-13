#version 330 core

in vec3  v_start;
in vec3  v_end;
in vec4  v_colorA;
in vec4  v_colorB;
in float v_radius;
in vec3  v_corner;

uniform mat4 u_view;
uniform mat4 u_proj;
uniform vec3 u_lightDir = normalize(vec3(0.4, 0.8, 1.2));

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

float viewPosToDepth(vec3 viewPos)
{
    vec4 clip   = u_proj * vec4(viewPos, 1.0);
    float ndc_z = clip.z / clip.w;
    return 0.5f * ndc_z + 0.5f;
}

void main()
{
    vec3 ro = vec3(0.0);
    vec3 rd = normalize(v_corner);

    if (length(v_corner) > 500.0 * v_radius) discard;

    vec3 center = (v_start + v_end) * 0.5;
    vec2 tt = sphIntersect(ro, rd, center, v_radius);

    float t = -1.0;
    if (tt.x > 0.001) t = tt.x;
    else if (tt.y > 0.001) t = tt.y;

    if (t < 0.001) 
    {
        discard;
    }

    vec3 hit = ro + t * rd;

    vec3 ba = v_end - v_start;
    float len2 = dot(ba, ba);
    float h = dot(hit - v_start, ba) / len2;

    if (h < -0.05 || h > 1.05) 
    {
        discard;
    }

    gl_FragDepth = viewPosToDepth(hit);

    vec3 hitPos_view = ro + t * rd;
    vec3 hit_normal = normalize(hitPos_view - center);

    vec4 final_color = mix(v_colorA, v_colorB, h);

    vec3 normal = normalize(hit - center);
    float NdotL = max(0.0, dot(normal, u_lightDir));

    fragColor = vec4(vec3(1.0, 1.0, 1.0) * (0.2 + 0.8 * NdotL), final_color.w);
}