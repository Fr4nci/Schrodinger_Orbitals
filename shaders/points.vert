#version 330 core
layout (location = 0) in vec3 aPos;   // posizione 3D del punto campionato (x, y, z)
layout (location = 1) in float aSign; // segno della funzione d'onda in quel punto (+1 / -1)

uniform mat4 view;
uniform mat4 proj;
uniform float pointSize;

out float vSign;

void main() {
    gl_Position = proj * view * vec4(aPos, 1.0);
    gl_PointSize = pointSize;
    vSign = aSign;
}
