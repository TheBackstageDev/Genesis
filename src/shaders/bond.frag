#version 330 core

in vec3  v_start;
in vec3  v_end;
in vec4  v_colorStart;
in vec4  v_colorEnd;
in float v_radius;
in vec3  v_corner;

uniform mat4 u_view;
uniform mat4 u_proj;
uniform vec3 u_lightDir = normalize(vec3(0.4, 0.8, 1.2));

out vec4 fragColor;

vec2 cylIntersect(in vec3 ro, in vec3 rd, in vec3 c0, in vec3 c1, in float ra) 
{
    vec3  ba = c1 - c0;
    vec3  ca = ro - c0;
    float baba = dot(ba, ba);
    float bard = dot(ba, rd);
    float baca = dot(ba, ca);
    float card = dot(ca, rd);
    float caca = dot(ca, ca);

    float k2 = baba - bard*bard;
    if (abs(k2) < 1e-6) return vec2(-1.0);

    float k1 = baba * card - baca * bard;
    float k0 = baba * caca - baca * baca - ra*ra * baba;

    float h = k1*k1 - k0*k2;
    if (h < 0.0) return vec2(-1.0);
    h = sqrt(h);

    float t1 = (-k1 - h) / k2;
    float t2 = (-k1 + h) / k2;

    return vec2(t1, t2);
}

float viewPosToDepth(vec3 viewPos) 
{
    vec4 clip   = u_proj * vec4(viewPos, 1.0);
    float ndc_z = clip.z / clip.w;
    return 0.5 * ndc_z + 0.5;
}

void main()
{
    vec3 ro = vec3(0.0);
    vec3 rd = normalize(v_corner);

    vec3 center = (v_start + v_end) * 0.5;
    vec2 tt = cylIntersect(ro, rd, v_start, v_end, v_radius);

    float t = -1.0;
    if (tt.x > 0.001) t = tt.x;
    else if (tt.y > 0.001) t = tt.y;

    if (t < 0) 
    {
        discard;
    }

    vec3 hit = ro + t * rd;

    vec3 ba = v_end - v_start;
    float len2 = dot(ba, ba);
    float h = dot(hit - v_start, ba) / len2;

    if (h < 0.0 || h > 1.0) 
    {
        discard;
    }

    gl_FragDepth = viewPosToDepth(hit);

    vec3 closest = v_start + h * ba;

    vec3 normal = normalize(hit - closest);
    vec4 final_color = mix(v_colorStart, v_colorEnd, clamp(h, 0.0, 1.0));

    vec3 viewDir = normalize(-closest);
    float fresnel = pow(1.0 - abs(dot(normal, viewDir)), 4.0);
    final_color += vec4(0.9, 0.9, 1.0, 1.0) * fresnel * 0.4;

    float NdotL = max(0.0, dot(normal, u_lightDir));
    fragColor = vec4(final_color.xyz * (0.1 + 0.8 * NdotL), final_color.w);
}