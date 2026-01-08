#version 330 core

in vec3  v_start;
in vec3  v_end;
in vec4  v_colorA;
in vec4  v_colorB;
in float v_radius;
in vec3  v_corner;

uniform mat4 u_view;
uniform mat4 u_proj;
uniform vec3 u_cameraPos;
uniform vec3 u_lightDir = normalize(vec3(0.4, 0.8, 1.2));

out vec4 fragColor;

vec2 cylIntersect(vec3 ro, vec3 rd, vec3 pa, vec3 pb, float ra)
{
    vec3  ba = pb - pa;
    vec3  oa = ro - pa;
    float baba = dot(ba, ba);
    float bard = dot(ba, rd);
    float baoa = dot(ba, oa);
    float a    = baba - bard * bard;
    float b    = baoa * bard - dot(oa, rd) * baba;
    float c    = dot(oa, oa) * baba - baoa * baoa - ra * ra * baba;

    if (abs(a) < 1e-6)
    {
        if (dot(oa - (baoa / baba) * ba, oa - (baoa / baba) * ba) > ra * ra)
            return vec2(-1.0);
        return vec2(-1e6, 1e6);
    }

    float h = b * b - a * c;
    if (h < 0.0) return vec2(-1.0);
    h = sqrt(h);
    return vec2((-b - h) / a, (-b + h) / a);
}

float viewPosToDepth(vec3 viewPos)
{
    vec4 clip   = u_proj * vec4(viewPos, 1.0);
    float ndc_z = clip.z / clip.w;
    return 0.5f * ndc_z + 0.5f;
}

void main()
{
    vec3 rayOrig = u_cameraPos;
    vec3 rayDir  = normalize(v_corner - rayOrig);

    vec3 axis    = v_end - v_start;
    float axisLen = length(axis);
    if (axisLen < 1e-5f) discard;

    vec2 tt = cylIntersect(rayOrig, rayDir, v_start, v_end, v_radius);

    float t = (tt.x > 0.0) ? tt.x : tt.y;
    if (t > 0.0001f) {
        fragColor = mix(v_colorA, v_colorB, 0.5);
    } else {
        fragColor = vec4(1.0, 0.2, 0.2, 1.0);
    }
}