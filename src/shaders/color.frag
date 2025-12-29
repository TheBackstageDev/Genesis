#version 330 core

in vec4  v_center;
in float v_radius;
in vec4  v_color;
in vec2  v_uv;

uniform vec3 u_lightDir = normalize(vec3(0.4, 0.8, 1.2));
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

float viewPosToDepth(vec3 viewPos) 
{
    vec4 clip = u_proj * vec4(viewPos, 1.0);
    float ndc_z = clip.z / clip.w;
    return 0.5 * ndc_z + 0.5;
}

void main()
{
    const vec3 rayOrig = vec3(0);

    float tanHalfFovY = 1.0 / u_proj[1][1];
    float aspect       = u_proj[1][1] / u_proj[0][0];
    vec3 rayDir  = normalize(vec3(v_uv.x * tanHalfFovY * aspect, v_uv.y * tanHalfFovY, -1.0));
    vec2 tt = sphIntersect(rayOrig, rayDir, v_center.xyz, v_radius);

    if (tt.x > tt.y || tt.y < 0.0)
    {
        discard;
    }

    float t = (tt.x >= 0.0) ? tt.x : tt.y;

    vec3 hitPos_view = rayOrig + t * rayDir;
    vec3 hit_normal = normalize(hitPos_view - v_center.xyz);

    gl_FragDepth = viewPosToDepth(hitPos_view);

    float NdotL = max(0.0, dot(hit_normal, u_lightDir));
    fragColor = vec4(v_color.xyz * (0.2 + 0.8 * NdotL), 1.0f);
}