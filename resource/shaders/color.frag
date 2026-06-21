#version 330 core

in vec3  v_center;
in float v_radius;
in vec4  v_color;
in vec2  v_uv;

uniform mat4 u_view;
uniform mat4 u_proj;
uniform vec3 u_lightDir = normalize(vec3(0.4, 0.8, 1.2));
uniform vec3 u_campos = vec3(0.0);

out vec4 fragColor;

vec2 sphIntersect(in vec3 ro, in vec3 rd, in vec3 ce, in float ra)
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
    const vec3 rayOrig = vec3(0.0);
    vec3 rayDir = normalize(vec3((v_uv.x * v_radius) + v_center.x, 
                                 (v_uv.y * v_radius) + v_center.y, 
                                 v_center.z));

    vec2 tt = sphIntersect(rayOrig, rayDir, v_center, v_radius);
    if (tt.x > tt.y || tt.y < 0.0)
    {
        discard;
    }

    float t = (tt.x >= 0.0) ? tt.x : tt.y;
    t = t < 0.0 ? dot(rayDir, v_center - rayOrig) : t;
   
    vec3 hitPos_view = rayOrig + t * rayDir;
    vec3 hit_normal = normalize(hitPos_view - v_center);

    gl_FragDepth = viewPosToDepth(hitPos_view);

    float ao = 0.5 + 0.5 * hit_normal.z;
    ao = pow(ao, 1.8);

    float edgeAO = 1.0 - smoothstep(0.5, 1.0, length(v_uv));
    ao = min(ao, edgeAO);

    float distAO = 1.0 - smoothstep(0.0, v_radius * 3.f, length(hitPos_view - v_center));
    ao *= distAO;

    float NdotL = max(0.0, dot(hit_normal, u_lightDir));
    
    vec3 ambient = v_color.rgb * 0.25;
    vec3 diffuse = v_color.rgb * 0.75 * NdotL;

    vec3 final_color = ambient + diffuse;

    final_color *= (0.7 + 0.3 * ao);

    vec3 viewDir = normalize(-hitPos_view);
    vec3 reflectDir = reflect(-u_lightDir, hit_normal);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 128.0);
    vec3 specular = 0.6 * spec * vec3(1.0);

    final_color += specular;

    fragColor = vec4(final_color, v_color.w);
}