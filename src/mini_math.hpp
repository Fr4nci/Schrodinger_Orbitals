#pragma once
#include <cmath>

// --- Libreria minimale di algebra lineare (sostituisce GLM) ---
// Matrici 4x4 in column-major order, compatibili con glUniformMatrix4fv.

struct Vec3 {
    float x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}

    Vec3 operator-(const Vec3& o) const { return Vec3(x - o.x, y - o.y, z - o.z); }
    Vec3 operator+(const Vec3& o) const { return Vec3(x + o.x, y + o.y, z + o.z); }
    Vec3 operator*(float s) const { return Vec3(x * s, y * s, z * s); }
};

inline float dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return Vec3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

inline Vec3 normalize(const Vec3& v) {
    float len = std::sqrt(dot(v, v));
    if (len < 1e-8f) return v;
    return Vec3(v.x / len, v.y / len, v.z / len);
}

struct Mat4 {
    float m[16]; // column-major

    static Mat4 identity() {
        Mat4 r{};
        for (int i = 0; i < 16; ++i) r.m[i] = 0.0f;
        r.m[0] = r.m[5] = r.m[10] = r.m[15] = 1.0f;
        return r;
    }
};

// Camera "guarda verso" (right-handed, come lookAt di GLM)
inline Mat4 lookAt(const Vec3& eye, const Vec3& center, const Vec3& up) {
    Vec3 f = normalize(center - eye);
    Vec3 s = normalize(cross(f, up));
    Vec3 u = cross(s, f);

    Mat4 r = Mat4::identity();
    r.m[0] = s.x;  r.m[4] = s.y;  r.m[8]  = s.z;
    r.m[1] = u.x;  r.m[5] = u.y;  r.m[9]  = u.z;
    r.m[2] = -f.x; r.m[6] = -f.y; r.m[10] = -f.z;
    r.m[12] = -dot(s, eye);
    r.m[13] = -dot(u, eye);
    r.m[14] = dot(f, eye);
    return r;
}

// Proiezione prospettica (fovy in radianti)
inline Mat4 perspective(float fovy, float aspect, float zNear, float zFar) {
    Mat4 r = Mat4::identity();
    float tanHalfFovy = std::tan(fovy / 2.0f);

    for (int i = 0; i < 16; ++i) r.m[i] = 0.0f;
    r.m[0] = 1.0f / (aspect * tanHalfFovy);
    r.m[5] = 1.0f / tanHalfFovy;
    r.m[10] = -(zFar + zNear) / (zFar - zNear);
    r.m[11] = -1.0f;
    r.m[14] = -(2.0f * zFar * zNear) / (zFar - zNear);
    return r;
}
