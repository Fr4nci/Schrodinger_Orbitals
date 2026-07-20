#version 330 core
in float vSign;
out vec4 FragColor;

void main() {
    // Rende il punto un cerchietto morbido invece di un quadrato
    vec2 coord = gl_PointCoord - vec2(0.5);
    float dist = length(coord);
    if (dist > 0.5) discard;

    // Falloff piu' ripido: nucleo del punto pieno e luminoso, bordo che sfuma
    float falloff = 1.0 - smoothstep(0.0, 0.5, dist);
    falloff = pow(falloff, 0.6);

    // Lobi con fase positiva in blu, lobi con fase negativa in arancio.
    // Colori piu' saturi e "accesi" rispetto a prima.
    vec3 colorPos = vec3(0.15, 0.55, 1.0);
    vec3 colorNeg = vec3(1.0, 0.35, 0.05);
    vec3 color = (vSign > 0.0) ? colorPos : colorNeg;

    // Boost di intensita': il centro del punto e' sovra-illuminato (>1),
    // cosi' col blending additivo le zone dense diventano brillanti invece
    // che restare di un blu/arancio smorto.
    float intensity = 2.2;
    vec3 litColor = color * intensity;

    float alpha = 0.55 * falloff;
    FragColor = vec4(litColor * alpha, alpha);
}
