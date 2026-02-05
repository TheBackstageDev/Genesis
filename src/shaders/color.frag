#version 330 core

in vec3  v_center;
in float v_radius;
in vec4  v_color;
in vec2  v_uv;

uniform mat4 u_view;
uniform mat4 u_proj;
uniform vec3 u_lightDir = normalize(vec3(0.4, 0.8, 1.2));
uniform vec3 u_campos = vec3(0, 0, 0);

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
    vec3 rayDir = normalize(vec3((v_uv.x * v_radius) + v_center.x, (v_uv.y * v_radius) + v_center.y, v_center.z));

    vec2 tt = sphIntersect(rayOrig, rayDir, v_center, v_radius);

    if (tt.x > tt.y || tt.y < 0.0)
    {
        discard;
    }

    float t = (tt.x >= 0.0) ? tt.x : tt.y;
    
    vec3 hitPos_view = rayOrig + t * rayDir;
    vec3 hit_normal = normalize(hitPos_view - v_center);

    gl_FragDepth = viewPosToDepth(hitPos_view);

    float NdotL = max(0.0, dot(hit_normal, u_lightDir));
    vec3 final_color = v_color.rgb * (0.2 + 0.8 * NdotL);
    
    vec3 viewDir = normalize(-hitPos_view);
    float fresnel = pow(1.0 - abs(dot(hit_normal, viewDir)), 4.0);
    final_color += vec3(0.9, 0.9, 1.0) * fresnel * 0.4;

    vec3 reflectDir = reflect(-u_lightDir, viewDir);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32);
    vec3 specular = 10.f * spec * v_color.rgb;   

    fragColor.xyz += specular;

    fragColor = vec4(final_color, v_color.w);
}