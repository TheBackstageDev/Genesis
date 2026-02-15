#version 330 core

in vec3  v_start;    // in view space
in vec3  v_end;      // in view space (arrow tip)
in vec4  v_color;
in float v_radius;
in vec3  v_corner;   // ray direction in view space for this fragment

uniform mat4 u_proj;
uniform vec3 u_lightDir = normalize(vec3(0.4, 0.8, 1.2));

out vec4 fragColor;

vec2 cylIntersect(in vec3 ro, in vec3 rd, in vec3 c0, in vec3 c1, in float ra)
{
    vec3 ba = c1 - c0;
    vec3 ca = ro - c0;
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

vec2 coneIntersect(in vec3 ro, in vec3 rd, in vec3 tip, in vec3 base, in float angle)
{
    vec3 axis = normalize(base - tip);        // axis pointing from tip toward base
    vec3 co = ro - tip;

    float cosA = cos(angle);
    float cosA2 = cosA * cosA;

    float dv = dot(rd, axis);
    float coAxis = dot(co, axis);

    float a = dv*dv - cosA2;
    float b = 2.0*(dv*coAxis - dot(rd, co)*cosA2);
    float c = coAxis*coAxis - dot(co, co)*cosA2;

    float disc = b*b - 4.0*a*c;
    if (disc < 0.0) return vec2(-1.0);

    float s = sqrt(disc);
    float t1 = (-b - s) / (2.0*a);
    float t2 = (-b + s) / (2.0*a);

    return vec2(t1, t2);
}

float planeIntersect(in vec3 ro, in vec3 rd, in vec3 p, in vec3 n)
{
    float denom = dot(rd, n);
    if (abs(denom) < 1e-6) return -1.0;
    return dot(p - ro, n) / denom;
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

    vec3 shaftA = v_start;
    vec3 shaftB = v_end;

    float totalLen = max(1e-6, length(shaftB - shaftA));
    float head_length = v_radius * 10.f;
    float head_angle  = radians(50.0);

    vec3 cone_tip = shaftB; 
    vec3 cone_base = cone_tip - normalize(shaftB - shaftA) * 0.3f;

    vec2 tc = cylIntersect(ro, rd, shaftA, cone_base, v_radius);
    vec2 th = coneIntersect(ro, rd, cone_tip, cone_base, head_angle);

    vec3 axis = normalize(cone_base - shaftA);
    vec3 cap0_n = -axis;
    vec3 cap1_n = axis;

    float t_cap0 = planeIntersect(ro, rd, shaftA, cap0_n);
    float t_cap1 = planeIntersect(ro, rd, cone_base, cap1_n);

    if (t_cap0 > 0.001) {
        vec3 p = ro + rd * t_cap0;
        if (dot(p - shaftA, p - shaftA) > v_radius * v_radius) t_cap0 = -1.0;
    } else t_cap0 = -1.0;

    if (t_cap1 > 0.001) {
        vec3 p = ro + rd * t_cap1;
        if (dot(p - cone_base, p - cone_base) > v_radius * v_radius) t_cap1 = -1.0;
    } else t_cap1 = -1.0;

    float t_coneBase = planeIntersect(ro, rd, cone_base, cap1_n);
    if (t_coneBase > 0.001) {
        vec3 p = ro + rd * t_coneBase;
        float baseRadius = head_length / 3.f;
        if (dot(p - cone_base, p - cone_base) > baseRadius * baseRadius) t_coneBase = -1.0;
    } else t_coneBase = -1.0;

    float t = 1e30;
    int hitType = 0; // 1=cylinder side, 2=shaft cap, 3=cylinder-cone junction cap, 4=cone side, 5=cone base disk

    if (tc.x > 0.001) {
        vec3 hitc = ro + rd * tc.x;
        float hproj = dot(hitc - shaftA, cone_base - shaftA) / dot(cone_base - shaftA, cone_base - shaftA);
        if (hproj >= 0.0 && hproj <= 1.0 && tc.x < t) { t = tc.x; hitType = 1; }
    }
    if (tc.y > 0.001) {
        vec3 hitc = ro + rd * tc.y;
        float hproj = dot(hitc - shaftA, cone_base - shaftA) / dot(cone_base - shaftA, cone_base - shaftA);
        if (hproj >= 0.0 && hproj <= 1.0 && tc.y < t) { t = tc.y; hitType = 1; }
    }

    if (t_cap0 > 0.001 && t_cap0 < t) { t = t_cap0; hitType = 2; }
    if (t_cap1 > 0.001 && t_cap1 < t) { t = t_cap1; hitType = 3; }

    if (th.x > 0.001) {
        vec3 hit = ro + rd * th.x;
        float proj = dot(hit - cone_tip, cone_base - cone_tip) / dot(cone_base - cone_tip, cone_base - cone_tip);
        if (proj >= 0.0 && proj <= 1.0 && th.x < t) { t = th.x; hitType = 4; }
    }
    if (th.y > 0.001) {
        vec3 hit = ro + rd * th.y;
        float proj = dot(hit - cone_tip, cone_base - cone_tip) / dot(cone_base - cone_tip, cone_base - cone_tip);
        if (proj >= 0.0 && proj <= 1.0 && th.y < t) { t = th.y; hitType = 4; }
    }

    if (t_coneBase > 0.001 && t_coneBase < t) { t = t_coneBase; hitType = 5; }

    if (t > 1e29) discard;

    vec3 hitPos = ro + rd * t;
    gl_FragDepth = viewPosToDepth(hitPos);

    vec3 normal;
    if (hitType == 1) 
    {
        vec3 ba = cone_base - shaftA;
        float len2 = dot(ba, ba);
        float h = dot(hitPos - shaftA, ba) / len2;
        vec3 closest = shaftA + clamp(h, 0.0, 1.0) * ba;
        normal = normalize(hitPos - closest);
    } 
    else if (hitType == 2) 
    {
        normal = normalize(cap0_n);
    } 
    else if (hitType == 3) 
    {
        normal = normalize(cap1_n);
    } 
    else if (hitType == 5) 
    {
        normal = normalize(cap1_n);
    } 
    else 
    { 
        vec3 coneAxis = normalize(cone_base - cone_tip);
        float t_axis = dot(hitPos - cone_tip, coneAxis);
        vec3 centerOnAxis = cone_tip + coneAxis * t_axis;
        vec3 radial = hitPos - centerOnAxis;
        float radialLen = length(radial);
        float invTan = 1.0 / tan(head_angle);
        normal = normalize(radial - coneAxis * (radialLen * invTan));
    }

    vec4 final_color = v_color;
    float VdotN = dot(normalize(-hitPos), normal);
    float rim = smoothstep(0.3, 0.0, abs(VdotN));
    rim = pow(rim, 1.6);
    vec3 rimColor = vec3(0.9, 0.95, 1.0);
    final_color.rgb += rimColor * rim * 0.55;

    float NdotL = max(0.0, dot(normal, normalize(u_lightDir)));
    fragColor = vec4(final_color.rgb * (0.1 + 0.8 * NdotL), final_color.a);
}
