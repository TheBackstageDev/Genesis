#version 330 core
flat in uvec3 v_id;
out uvec3 fragColor;

void main()
{
    fragColor = v_id;
}