#version 330 core

in vec3  v_start;
in vec3  v_end;
in vec4  v_colorStart;
in vec4  v_colorEnd;
in float v_radiusA;
in float v_radiusB;

uniform mat4 u_view;
uniform mat4 u_proj;
uniform vec3 u_campos;
uniform vec2 u_screen;
uniform vec3 u_lightDir = normalize(vec3(0.4, 0.8, 1.2));
uniform float u_time;

out vec4 fragColor;

float sphereSDF(vec3 p, vec3 center, float radius) 
{
    return length(p - center) - radius;
}

float smoothMin(float a, float b, float k) 
{
    float h = clamp(0.5 + 0.5*(b - a)/k, 0.0, 1.0);
    return mix(b, a, h) - k*h*(1.0 - h);
}

float sceneSDF(vec3 p) 
{
    float d0 = sphereSDF(p, v_start, v_radiusA * 0.9);
    float d1 = sphereSDF(p, v_end,   v_radiusB * 0.9);
    return smoothMin(d0, d1, 0.6);
}

vec3 calcNormal(vec3 p) 
{
    float e = 0.001;
    return normalize(vec3(
        sceneSDF(p + vec3(e,0,0)) - sceneSDF(p - vec3(e,0,0)),
        sceneSDF(p + vec3(0,e,0)) - sceneSDF(p - vec3(0,e,0)),
        sceneSDF(p + vec3(0,0,e)) - sceneSDF(p - vec3(0,0,e))
    ));
}

vec4 blendedColor(vec3 p) 
{
    float d0 = sphereSDF(p, v_start, v_radiusA * 0.9);
    float d1 = sphereSDF(p, v_end,   v_radiusB * 0.9);
    float w = smoothstep(-0.1, 0.1, d1 - d0);
    return mix(v_colorStart, v_colorEnd, w);
}

void main() {
    vec2 uv = (gl_FragCoord.xy / u_screen) * 2.0 - 1.0;
    vec4 rayClip = vec4(uv, -1.0, 1.0);
    vec4 rayWorld = inverse(u_proj * u_view) * rayClip;
    rayWorld /= rayWorld.w;

    vec3 ro = u_campos;
    vec3 rd = normalize(rayWorld.xyz - ro);

    float t = 0.0;
    vec3 hit;
    bool found = false;
    for (int i = 0; i < 128; i++) {
        vec3 p = ro + rd * t;
        float dist = sceneSDF(p);
        if (dist < 0.001) 
        {
            hit = p;
            found = true;
            break;
        }
        t += dist;
        if (t > 100.0) break;
    }

    if (!found) discard;

    vec4 clipPos = u_proj * u_view * vec4(hit, 1.0);
    float ndcDepth = clipPos.z / clipPos.w;
    gl_FragDepth = 0.5 * ndcDepth + 0.5;

    vec3 N = calcNormal(hit);
    vec3 V = normalize(u_campos - N);
    vec3 L = normalize(u_campos - u_lightDir);
    vec3 H = normalize(L + V);

    float specPower = 256.0;
    float specStrength = max(dot(N, H), 0.0);
    float specular = pow(specStrength, specPower);

    float diffuse = max(dot(N, normalize(u_lightDir)), 0.0);

    vec4 baseColor = blendedColor(hit);
    vec3 color = baseColor.rgb * diffuse;

    float glow = 0.5 + 0.5 * sin(u_time * 2.0);
    color += baseColor.rgb * glow * 0.1;

    fragColor = vec4(color + specular, baseColor.a);
}
