#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <random>

#include "mini_math.hpp"

// ============================================================
//  Integratore Radiale (Numerov)
// ============================================================
std::pair<std::vector<double>, int> solve_radial(double E, int l, double r_max, int steps) {
    double h = r_max / steps;
    std::vector<double> u(steps + 1, 0.0);

    // Inizializzazione analitica per evitare la singolarita' matematica 1/r^2 nell'origine
    int start_idx = 5;
    for (int i = 0; i <= start_idx; ++i) {
        u[i] = std::pow(i * h, l + 1);
    }

    int nodes = 0;

    for (int i = start_idx; i < steps; ++i) {
        double r_curr = i * h;
        double r_prev = (i - 1) * h;
        double r_next = (i + 1) * h;

        auto g = [&](double ri) { return (double)l * (l + 1) / (ri * ri) - 2.0 / ri - E; };

        double g_curr = g(r_curr);
        double g_next = g(r_next);
        double g_prev = g(r_prev);

        double c_curr = 2.0 * (1.0 + (5.0 * h * h / 12.0) * g_curr);
        double c_prev = 1.0 - (h * h / 12.0) * g_prev;
        double c_next = 1.0 - (h * h / 12.0) * g_next;

        u[i + 1] = (c_curr * u[i] - c_prev * u[i - 1]) / c_next;

        if (std::isnan(u[i + 1])) u[i + 1] = 1e10;

        if (u[i] * u[i + 1] < 0) nodes++;

        if (std::abs(u[i + 1]) > 1e10) {
            for (double& val : u) val /= 1e10;
        }
    }
    return {u, nodes};
}

// --- Metodo di Shooting ---
std::vector<float> calcola_radiale(int n, int l, double& r_max) {
    int target_nodes = n - l - 1;
    // Diamo margine sufficiente per far decadere a zero la coda esponenziale.
    r_max = 20.0 + 2.5 * (n * n);

    // Risoluzione matematica altissima, slegata dalla GPU come facevo prima (non funzionava per alti n)
    int math_steps = 150000;

    double E_min = -1.5 / (n * n);
    double E_max = -0.5 / (n * n);
    std::vector<double> u_final;

    for(int iter = 0; iter < 100; ++iter) {
        double E_mid = (E_min + E_max) / 2.0;

        auto res = solve_radial(E_mid, l, r_max, math_steps);

        if (res.second > target_nodes) {
            E_max = E_mid;
        } else if (res.second < target_nodes) {
            E_min = E_mid;
        } else {
            if (target_nodes % 2 == 0) {
                if (res.first.back() > 0) E_min = E_mid;
                else E_max = E_mid;
            } else {
                if (res.first.back() < 0) E_min = E_mid;
                else E_max = E_mid;
            }
        }
        if (iter == 99) u_final = res.first;
    }

    double E_converged = (E_min + E_max) / 2.0;
    std::cout << "[INFO] Energia convergente E = " << E_converged << std::endl;
    std::cout << "[INFO] Valore teorico atteso = " << -1.0 / (n * n) << std::endl;

    r_max = 20.0 + 25.0 * n; // This is required for rendering orbitals less "bigger" and more visible (for larger n this produces a better result)
    double h = r_max / math_steps;
    double norm_factor = 0.0;
    std::vector<double> R_full(math_steps + 1, 0.0);

    // Calcoliamo R(r) e normalizziamo sull'array ad alta definizione
    for (int i = 1; i <= math_steps; ++i) {
        double r = i * h;
        R_full[i] = u_final[i] / r;
        norm_factor += R_full[i] * R_full[i] * r * r * h;
    }
    norm_factor = std::sqrt(norm_factor);
    std::cout << "[INFO] Fattore di normalizzazione = " << norm_factor << std::endl;

    // Sottocampionamento sicuro per la Texture GPU (limite hardware)
    int tex_steps = 8192;
    std::vector<float> R_gpu(tex_steps, 0.0f);

    for(int i = 0; i < tex_steps; ++i) {
        // Estraiamo in modo proporzionale l'indice dall'array gigantesco
        int idx = (i * math_steps) / tex_steps;
        if(idx == 0) idx = 1; // Evita la singolarità in r=0
        R_gpu[i] = (float)(R_full[idx] / norm_factor);
    }

    return R_gpu;
}
// ============================================================
//  MATEMATICA: Polinomio di Legendre associato (precisione double)
//  Stessa formula usata prima nel fragment shader, portata in C++
//  cosi' il campionamento sulla CPU e' coerente con la fisica del progetto.
// ============================================================
double legendre_assoc(int l, int m, double x) {
    double pmm = 1.0;
    if (m > 0) {
        double somx2 = std::sqrt(std::max(0.0, (1.0 - x) * (1.0 + x)));
        double fact = 1.0;
        for (int i = 1; i <= m; i++) {
            pmm *= -fact * somx2;
            fact += 2.0;
        }
    }
    if (l == m) return pmm;

    double pmmp1 = x * (2.0 * m + 1.0) * pmm;
    if (l == m + 1) return pmmp1;

    double pll = 0.0;
    for (int ll = m + 2; ll <= l; ll++) {
        pll = (x * (2.0 * ll - 1.0) * pmmp1 - (double)(ll + m - 1) * pmm) / (double)(ll - m);
        pmm = pmmp1;
        pmmp1 = pll;
    }
    return pll;
}

// Interpola R_nl(r) dalla tabella calcolata da calcola_radiale
double sample_R(double r, const std::vector<float>& R_table, double r_max) {
    if (r >= r_max || r < 0.0) return 0.0;
    double u = (r / r_max) * (double)(R_table.size() - 1);
    int i0 = (int)u;
    int i1 = std::min(i0 + 1, (int)R_table.size() - 1);
    double frac = u - i0;
    return R_table[i0] * (1.0 - frac) + R_table[i1] * frac;
}

// psi(x,y,z) = R_nl(r) * Y_lm(theta, phi)  (forma reale, stessa convenzione del progetto originale)
double psi_cartesian(const Vec3& p, int l, int m, const std::vector<float>& R_table, double r_max) {
    double r = std::sqrt((double)p.x * p.x + (double)p.y * p.y + (double)p.z * p.z);
    if (r < 1e-9 || r >= r_max) return 0.0;

    double theta = std::acos((double)p.z / r);
    double phi = std::atan2((double)p.y, (double)p.x);

    double R_val = sample_R(r, R_table, r_max);
    double P_lm = legendre_assoc(l, std::abs(m), std::cos(theta));
    double Y_lm = P_lm * std::cos((double)m * phi);

    return R_val * Y_lm;
}

// ============================================================
//  CAMPIONAMENTO MONTE CARLO: genera punti distribuiti secondo |psi|^2
//
//  Idea: campioniamo (x,y,z) uniformemente in un cubo che contiene
//  l'orbitale, calcoliamo la densita' di probabilita' |psi(x,y,z)|^2
//  in quel punto, e accettiamo il campione con probabilita'
//  proporzionale a quella densita' (rejection sampling classico).
//  Il numero di punti accettati in una regione e' quindi proporzionale
//  al modulo quadro della funzione d'onda: la "nuvola" di punti
//  rappresenta direttamente la densita' di probabilita'.
// ============================================================
struct PointCloud {
    std::vector<float> positions; // x,y,z per ogni punto
    std::vector<float> signs;     // segno di psi in quel punto (+1/-1), per colorare i lobi
};

PointCloud generate_orbital_points(int l, int m, const std::vector<float>& R_table,
                                    double r_max, int n_points, double bound) {
    PointCloud cloud;
    cloud.positions.reserve(n_points * 3);
    cloud.signs.reserve(n_points);

    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> uni(-bound, bound);
    std::uniform_real_distribution<double> uni01(0.0, 1.0);

    // Stima la densita' massima con una scansione grossolana, per calibrare
    // la soglia di accettazione del rejection sampling.
    double p_max = 1e-12;
    const int scan_samples = 200000;
    for (int i = 0; i < scan_samples; ++i) {
        Vec3 p(uni(rng), uni(rng), uni(rng));
        double psi = psi_cartesian(p, l, m, R_table, r_max);
        double dens = psi * psi;
        if (dens > p_max) p_max = dens;
    }
    p_max *= 1.15; // margine di sicurezza

    long long attempts = 0;
    const long long max_attempts = (long long)n_points * 4000;

    while ((int)cloud.signs.size() < n_points && attempts < max_attempts) {
        attempts++;
        Vec3 p(uni(rng), uni(rng), uni(rng));
        double psi = psi_cartesian(p, l, m, R_table, r_max);
        double dens = psi * psi;

        if (dens <= 0.0) continue;
        if (uni01(rng) * p_max < dens) {
            cloud.positions.push_back(p.x);
            cloud.positions.push_back(p.y);
            cloud.positions.push_back(p.z);
            cloud.signs.push_back(psi >= 0.0 ? 1.0f : -1.0f);
        }
    }

    std::cout << "[INFO] Generati " << cloud.signs.size() << " punti in " << attempts << " tentativi." << std::endl;
    return cloud;
}

// ============================================================
//  OPENGL: Lettura e compilazione shader
// ============================================================
std::string readFile(const char* filePath) {
    std::ifstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "ERRORE FATALE: Non trovo il file " << filePath << std::endl;
        exit(1);
    }
    std::stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

GLuint compileShaders(const char* vertPath, const char* fragPath) {
    std::string vertCode = readFile(vertPath);
    std::string fragCode = readFile(fragPath);
    const char* vShaderCode = vertCode.c_str();
    const char* fShaderCode = fragCode.c_str();

    GLuint vertex, fragment;
    GLint success;
    GLchar infoLog[1024];

    vertex = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex, 1, &vShaderCode, NULL);
    glCompileShader(vertex);
    glGetShaderiv(vertex, GL_COMPILE_STATUS, &success);
    if (!success) { glGetShaderInfoLog(vertex, 1024, NULL, infoLog); std::cerr << "ERRORE VERTEX:\n" << infoLog; }

    fragment = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment, 1, &fShaderCode, NULL);
    glCompileShader(fragment);
    glGetShaderiv(fragment, GL_COMPILE_STATUS, &success);
    if (!success) { glGetShaderInfoLog(fragment, 1024, NULL, infoLog); std::cerr << "ERRORE FRAGMENT:\n" << infoLog; }

    GLuint ID = glCreateProgram();
    glAttachShader(ID, vertex);
    glAttachShader(ID, fragment);
    glLinkProgram(ID);
    glGetProgramiv(ID, GL_LINK_STATUS, &success);
    if (!success) { glGetProgramInfoLog(ID, 1024, NULL, infoLog); std::cerr << "ERRORE LINK:\n" << infoLog; }

    glDeleteShader(vertex);
    glDeleteShader(fragment);
    return ID;
}

// ============================================================
//  CAMERA ORBITALE controllata col mouse (tasto sinistro premuto)
// ============================================================
struct OrbitCamera {
    float yaw = -90.0f * 3.14159265f / 180.0f; // radianti
    float pitch = 20.0f * 3.14159265f / 180.0f;
    float distance = 40.0f;
    Vec3 target = Vec3(0.0f, 0.0f, 0.0f);

    Vec3 eye() const {
        float cp = std::cos(pitch), sp = std::sin(pitch);
        float cy = std::cos(yaw), sy = std::sin(yaw);
        return target + Vec3(distance * cp * cy, distance * sp, distance * cp * sy);
    }
};

OrbitCamera g_camera;
bool g_dragging = false;
double g_lastX = 0.0, g_lastY = 0.0;

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if (action == GLFW_PRESS) {
            g_dragging = true;
            glfwGetCursorPos(window, &g_lastX, &g_lastY);
        } else if (action == GLFW_RELEASE) {
            g_dragging = false;
        }
    }
}

void cursor_position_callback(GLFWwindow* window, double xpos, double ypos) {
    if (!g_dragging) return;

    double dx = xpos - g_lastX;
    double dy = ypos - g_lastY;
    g_lastX = xpos;
    g_lastY = ypos;

    const float sensitivity = 0.005f;
    g_camera.yaw += (float)dx * sensitivity;
    g_camera.pitch += (float)dy * sensitivity;

    // Evita che la camera "ribalti" ai poli
    const float limit = 1.55f; // ~89 gradi
    if (g_camera.pitch > limit) g_camera.pitch = limit;
    if (g_camera.pitch < -limit) g_camera.pitch = -limit;
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
    g_camera.distance -= (float)yoffset * g_camera.distance * 0.1f;
    if (g_camera.distance < 1.0f) g_camera.distance = 1.0f;
}

int main() {
    int n, l, m;
    std::cout << "Inserisci n, l, m: ";
    if (!(std::cin >> n >> l >> m) || n < 1 || l < 0 || l >= n || std::abs(m) > l) {
        std::cerr << "Parametri quantici non validi." << std::endl;
        return 1;
    }

    double r_max;
    std::vector<float> R_r = calcola_radiale(n, l, r_max);

    // Bounding box del campionamento Monte Carlo: la maggior parte della
    // densita' e' concentrata ben dentro r_max, quindi restringiamo un po'.
    double bound = r_max * 0.6;
    int n_points = 150000;

    PointCloud cloud = generate_orbital_points(l, m, R_r, r_max, n_points, bound);

    // --- SETUP OPENGL / FINESTRA ---
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow* window = glfwCreateWindow(900, 900, "Orbitale - point cloud", NULL, NULL);
    if (!window) {
        std::cerr << "Impossibile creare la finestra GLFW." << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetCursorPosCallback(window, cursor_position_callback);
    glfwSetScrollCallback(window, scroll_callback);

    g_camera.distance = (float)r_max * 0.5f;

    GLuint shaderProgram = compileShaders("shaders/points.vert", "shaders/points.frag");

    // Interleaviamo posizione (3 float) + segno (1 float) per ogni punto
    std::vector<float> vertexData;
    vertexData.reserve(cloud.signs.size() * 4);
    for (size_t i = 0; i < cloud.signs.size(); ++i) {
        vertexData.push_back(cloud.positions[i * 3 + 0]);
        vertexData.push_back(cloud.positions[i * 3 + 1]);
        vertexData.push_back(cloud.positions[i * 3 + 2]);
        vertexData.push_back(cloud.signs[i]);
    }

    unsigned int VBO, VAO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertexData.size() * sizeof(float), vertexData.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    GLint viewLoc = glGetUniformLocation(shaderProgram, "view");
    GLint projLoc = glGetUniformLocation(shaderProgram, "proj");
    GLint pointSizeLoc = glGetUniformLocation(shaderProgram, "pointSize");

    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE); // blending additivo: le zone dense "si accendono" di piu'
    glDisable(GL_DEPTH_TEST);          // evita problemi di ordinamento con la trasparenza additiva

    GLsizei n_verts = (GLsizei)cloud.signs.size();

    // --- RENDER LOOP ---
    while (!glfwWindowShouldClose(window)) {
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        glViewport(0, 0, width, height);
        float aspect = (height > 0) ? (float)width / (float)height : 1.0f;

        glClearColor(0.03f, 0.03f, 0.05f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        Vec3 eye = g_camera.eye();
        Mat4 view = lookAt(eye, g_camera.target, Vec3(0.0f, 1.0f, 0.0f));
        Mat4 proj = perspective(45.0f * 3.14159265f / 180.0f, aspect, 0.1f, 1000.0f);

        glUseProgram(shaderProgram);
        glUniformMatrix4fv(viewLoc, 1, GL_FALSE, view.m);
        glUniformMatrix4fv(projLoc, 1, GL_FALSE, proj.m);
        glUniform1f(pointSizeLoc, 2.5f);

        glBindVertexArray(VAO);
        glDrawArrays(GL_POINTS, 0, n_verts);

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}
