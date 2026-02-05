#version 330 core

in vec3  v_start;       // view-space center A
in vec3  v_end;         // view-space center B
in vec4  v_colorStart;
in vec4  v_colorEnd;
in float v_radius;
in float v_radiusA;
in float v_radiusB;
in vec3  v_corner;

uniform mat4 u_view;
uniform mat4 u_proj;
uniform vec3 u_lightDir = normalize(vec3(0.4, 0.8, 1.2));
uniform float hyperLambda = 0.85;

out vec4 fragColor;

vec2 sphIntersect(vec3 ro, vec3 rd, vec3 ce, float ra)
{
    vec3 oc = ro - ce;
    float b = dot(oc, rd);
    float c = dot(oc, oc) - ra*ra;
    float h = b*b - c;
    
    if (h < 0.0) return vec2(-1.0);
    h = sqrt(h);
    return vec2(-b - h, -b + h);
}

float viewPosToDepth(vec3 viewPos)
{
    vec4 clip = u_proj * vec4(viewPos, 1.0);
    float ndc_z = clip.z / clip.w;
    return 0.5 * ndc_z + 0.5;
}

void main()
{
    vec3 ro = vec3(0.0);
    vec3 rd = normalize(v_corner);

    vec3 C1 = v_start;
    vec3 C2 = v_end;
    float r1 = v_radiusA;
    float r2 = v_radiusB;

    float t = -1.0;
    vec3 hitPos;
    vec3 normal;
    vec4 baseColor;

    // Sphere A
    vec2 ttA = sphIntersect(ro, rd, C1, r1);
    float tA = -1.0;
    if (ttA.x > 0.001) tA = ttA.x;
    else if (ttA.y > 0.001) tA = ttA.y;

    if (tA > 0.001)
    {
        t = tA;
        hitPos = ro + t * rd;
        normal = normalize(hitPos - C1);
        baseColor = v_colorStart;
    }
    else
    {
        // Sphere B
        vec2 ttB = sphIntersect(ro, rd, C2, r2);
        float tB = -1.0;
        if (ttB.x > 0.001) tB = ttB.x;
        else if (ttB.y > 0.001) tB = ttB.y;

        if (tB > 0.001)
        {
            t = tB;
            hitPos = ro + t * rd;
            normal = normalize(hitPos - C2);
            baseColor = v_colorEnd;
        }
    }

    if (t < 0.0)
    {
        vec3 ba     = C2 - C1;
        float lenBA = length(ba);
        if (lenBA < 0.001) discard;

        vec3 dir    = ba / lenBA;
        vec3 oc     = ro - C1;
        float proj  = dot(oc, dir);
        vec3 oc_perp = oc - proj * dir;

        float rd_proj = dot(rd, dir);
        vec3 rd_perp  = rd - rd_proj * dir;

        float ll = hyperLambda * hyperLambda;

        float qa = dot(rd_perp, rd_perp) - ll * rd_proj * rd_proj;
        float qb = 2.0 * (dot(oc_perp, rd_perp) - ll * proj * rd_proj);
        float qc = dot(oc_perp, oc_perp) - ll * proj * proj - r1*r1;

        float disc = qb*qb - 4.0*qa*qc;
        if (disc < 0.0) discard;

        float sqrtD = sqrt(disc);
        float t0 = (-qb - sqrtD) / (2.0 * qa);
        float t1 = (-qb + sqrtD) / (2.0 * qa);

        t = -1.0;
        if (t0 > 0.001) t = t0;
        if (t1 > 0.001 && (t < 0.0 || t1 < t)) t = t1;

        if (t < 0.001) discard;

        hitPos = ro + t * rd;

        float h = dot(hitPos - C1, dir);

        float margin = 0.05 * lenBA;

        if (h < -margin || h > lenBA + margin)
        {
            discard;
        }

        vec3 closestPoint = C1 + h * dir;
        float distToAxis = length(hitPos - closestPoint);
        float expectedR = mix(r1, r2, clamp(h / lenBA, 0.0, 1.0));

        if (abs(distToAxis - expectedR) > 0.25 * expectedR)
        {
            discard;
        }

        vec3 posA = hitPos - C1;
        vec3 posB = hitPos - C2;
        normal = normalize(posA / (r1*r1) + posB / (r2*r2));

        float interp = clamp(h / lenBA, 0.0, 1.0);
        baseColor = mix(v_colorStart, v_colorEnd, interp);
    }

    if (t < 0.0) discard;

    gl_FragDepth = viewPosToDepth(hitPos);

    float NdotL = max(0.0, dot(normal, u_lightDir));
    float VdotN = dot(normalize(-hitPos), normal);
    float rim = smoothstep(0.3, 0.0, abs(VdotN));
    rim = pow(rim, 1.6);

    vec3 rimColor = vec3(0.9, 0.95, 1.0);
    vec3 litColor = baseColor.rgb * (0.1 + 0.8 * NdotL);
    litColor += rimColor * rim * 0.55;

    fragColor = vec4(litColor, baseColor.a);
}