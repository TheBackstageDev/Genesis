#version 430

struct Sphere {
    vec4 centerRadius; // xyz = center, w = radius
    vec4 color;        // RGBA
};

layout(std430, binding = 0) buffer SphereBuffer {
    Sphere spheres[];      // array of sphere data
};

uniform mat4 u_view;
uniform mat4 u_proj;
uniform vec3 u_campos;
uniform vec2 u_screen;
uniform vec3 u_lightDir = normalize(vec3(0.4, 0.8, 1.2));
uniform float u_smoothFactor;
uniform float u_time;
uniform uint u_count = 0;

out vec4 fragColor;

float smoothMin(float a, float b, float k) 
{
    float h = clamp(0.5 + 0.5*(b - a)/k, 0.0, 1.0);
    return mix(b, a, h) - k*h*(1.0 - h);
}

float sphereSDF(vec3 p, vec3 center, float radius) 
{
    return length(p - center) - radius;
}

float sceneSDF(vec3 p) 
{
    float d = 1e9;
    for (uint i = 0; i < u_count; ++i) 
    {
        vec3 c = spheres[i].centerRadius.xyz;
        float r = spheres[i].centerRadius.w;
        float di = sphereSDF(p, c, r);
        d = smoothMin(d, di, u_smoothFactor);
    }
    return d;
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
    float minDist = 1e9;
    uint idx = 0u;
    for (uint i = 0; i < u_count; ++i) 
    {
        float di = sphereSDF(p, spheres[i].centerRadius.xyz, spheres[i].centerRadius.w);
        if (di < minDist) 
        {
            minDist = di;
            idx = i;
        }
    }
    return spheres[idx].color;
}

float softShadow(vec3 ro, vec3 rd, float mint, float maxt, int steps) 
{
    float res = 1.0;
    float t = mint;
    for (int i = 0; i < steps; ++i) {
        float h = sceneSDF(ro + rd * t);
        if (h < 0.001) return 0.0;

        res = min(res, 10.0 * h / t);
        t += clamp(h, 0.01, 0.5);
        if (t > maxt) break;
    }
    return clamp(res, 0.0, 1.0);
}

float ambientOcclusion(vec3 p, vec3 N) 
{
    float ao = 0.0;
    float sca = 1.0;
    for (int i = 0; i < 5; ++i) 
    {
        float h = 0.02 + 0.1 * float(i); // sample distance
        float d = sceneSDF(p + N * h);
        ao += (h - d) * sca;
        sca *= 0.5;
    }
    return clamp(1.0 - ao, 0.0, 1.0);
}

void main() 
{
    vec2 uv = (gl_FragCoord.xy / u_screen) * 2.0 - 1.0;
    vec4 rayClip = vec4(uv, -1.0, 1.0);
    vec4 rayWorld = inverse(u_proj * u_view) * rayClip;
    rayWorld /= rayWorld.w;

    vec3 ro = u_campos;
    vec3 rd = normalize(rayWorld.xyz - ro);

    float t = 0.0;
    vec3 hit;
    bool found = false;
    for (int i = 0; i < 128; i++) 
    {
        vec3 p = ro + rd * t;
        float dist = sceneSDF(p);
        if (dist < 0.001) 
        {
            hit = p;
            found = true;
            break;
        }
        t += dist;
        if (t > 150.0) break;
    }

    if (!found) discard;

    vec4 clipPos = u_proj * u_view * vec4(hit, 1.0);
    float ndcDepth = clipPos.z / clipPos.w;
    gl_FragDepth = 0.5 * ndcDepth + 0.5;

    vec3 N = calcNormal(hit);
    vec3 V = normalize(u_campos - hit);
    vec3 L = normalize(u_lightDir);
    vec3 H = normalize(L + V);

    float specular = pow(max(dot(N, H), 0.0), 256.0);
    float diffuse = max(dot(N, L), 0.0);

    vec4 baseColor = blendedColor(hit);
    vec3 color = specular + baseColor.rgb * diffuse;

    float fresnel = pow(1.0 - max(dot(N, V), 0.0), 3.0);
    color += baseColor.rgb * fresnel * 0.3;

    float shadow = softShadow(hit, L, 5.0, 20.0, 64);
    color *= shadow;

    fragColor = vec4(color + 0.1, baseColor.a);
}
