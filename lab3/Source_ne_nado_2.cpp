// Сферы + Тетраэдры + Текстуры + SSAA + 3D-камера
// Сферы + Тетраэдры + Зеркальная плоскость
// обратная трассировка, тени, отражения, enabled, управление

#define _CRT_SECURE_NO_WARNINGS
#define NOMINMAX

#pragma comment(lib, "opengl32.lib")
#pragma comment(lib, "glfw3.lib")

#include <thread>
#include <vector>
#include <GLFW/glfw3.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <algorithm>

// =============== Вспомогательные структуры ===============
struct Vec2 {
    float x, y;
    Vec2(float x = 0, float y = 0) : x(x), y(y) {}
};

struct Vec3 {
    float x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(float x, float y, float z) : x(x), y(y), z(z) {}

    Vec3 operator+(const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
    Vec3 operator-(const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
    Vec3 operator-() const { return Vec3(-x, -y, -z); }
    Vec3 operator*(float t) const { return Vec3(x * t, y * t, z * t); }
    Vec3 operator/(float t) const { return t != 0 ? Vec3(x / t, y / t, z / t) : Vec3(0, 0, 0); }

    float dot(const Vec3& v) const { return x * v.x + y * v.y + z * v.z; }
    Vec3 cross(const Vec3& v) const {
        return Vec3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }

    float length() const { return sqrtf(dot(*this)); }
    Vec3 normalize() const {
        float len = length();
        return len > 1e-6f ? (*this) * (1.0f / len) : Vec3(0, 0, 1);
    }

    Vec3 reflect(const Vec3& n) const {
        return *this - n * (2.0f * this->dot(n));
    }
};

struct Color {
    float r, g, b;
    Color() : r(0), g(0), b(0) {}
    Color(float r, float g, float b) {
        this->r = r < 0 ? 0 : (r > 1 ? 1 : r);
        this->g = g < 0 ? 0 : (g > 1 ? 1 : g);
        this->b = b < 0 ? 0 : (b > 1 ? 1 : b);
    }

    Color operator*(float t) const { return Color(r * t, g * t, b * t); }
    Color operator+(const Color& c) const { return Color(r + c.r, g + c.g, b + c.b); }
    Color blend(const Color& c, float t) const {
        t = t < 0 ? 0 : (t > 1 ? 1 : t);
        return Color(r * (1 - t) + c.r * t, g * (1 - t) + c.g * t, b * (1 - t) + c.b * t);
    }

    unsigned char getR() const { return (unsigned char)(r * 255); }
    unsigned char getG() const { return (unsigned char)(g * 255); }
    unsigned char getB() const { return (unsigned char)(b * 255); }
};

// =============== BMP Loader ===============
struct Texture {
    unsigned char* data = nullptr;
    int width = 0, height = 0;

    ~Texture() { if (data) delete[] data; }

    Color sample(float u, float v) const {
        if (!data || width <= 0 || height <= 0) return Color(1, 0, 1);
        u = u - floorf(u);
        v = v - floorf(v);
        int x = (int)(u * width) % width;
        int y = (int)(v * height) % height;
        if (x < 0) x += width;
        if (y < 0) y += height;
        int idx = (y * width + x) * 3;
        return Color(data[idx + 2] / 255.0f, data[idx + 1] / 255.0f, data[idx] / 255.0f);
    }

    static Texture loadBMP(const char* filename) {
        Texture tex;
        FILE* f = fopen(filename, "rb");
        if (!f) return tex;
        unsigned char header[54];
        if (fread(header, 1, 54, f) != 54) { fclose(f); return tex; }
        if (header[0] != 'B' || header[1] != 'M') { fclose(f); return tex; }
        int dataPos = *(int*)&(header[0x0A]);
        int size = *(int*)&(header[0x22]);
        tex.width = *(int*)&(header[0x12]);
        tex.height = *(int*)&(header[0x16]);
        if (size == 0) size = tex.width * tex.height * 3;
        if (dataPos == 0) dataPos = 54;
        tex.data = new unsigned char[size];
        fread(tex.data, 1, size, f);
        fclose(f);
        return tex;
    }
};

// =============== Material ===============
struct Material {
    Color color;
    float ka, kd, ks, reflect;
    int shininess;
    Texture* texture;

    Material() : color(1, 1, 1), ka(0.1f), kd(0.7f), ks(0.3f), reflect(0.0f), shininess(32), texture(nullptr) {}
    Material(Color c, float a, float d, float s, float r, int sh, Texture* tex = nullptr)
        : color(c), ka(a), kd(d), ks(s), reflect(r), shininess(sh), texture(tex) {
    }
};

// =============== Ray & HitInfo ===============
struct Ray {
    Vec3 origin, direction;
    Ray(Vec3 o, Vec3 d) : origin(o), direction(d.normalize()) {}
    Vec3 at(float t) const { return origin + direction * t; }
};

struct HitInfo {
    bool hit = false;
    float t = 1e30f;
    Vec3 point, normal;
    Material mat;
    float u = 0, v = 0;
    HitInfo() {}
};

// =============== Shape ===============
struct Shape {
    bool enabled = true;
    virtual ~Shape() {}
    virtual HitInfo intersect(const Ray& ray) const = 0;
    virtual bool isEnabled() const { return enabled; }
};

// =============== Сфера ===============
struct Sphere : public Shape {
    Vec3 center; float radius; Material mat;
    Sphere(Vec3 c, float r, Material m, bool en = true) : center(c), radius(r), mat(m) { enabled = en; }

    HitInfo intersect(const Ray& ray) const override {
        if (!isEnabled()) return HitInfo();
        HitInfo info;
        Vec3 oc = ray.origin - center;
        float a = ray.direction.dot(ray.direction);
        float b = 2.0f * oc.dot(ray.direction);
        float c = oc.dot(oc) - radius * radius;
        float disc = b * b - 4 * a * c;
        if (disc < 0) return info;

        float sqrt_disc = sqrtf(disc);
        float t1 = (-b - sqrt_disc) / (2 * a);
        float t2 = (-b + sqrt_disc) / (2 * a);

        float t = 1e30f;
        if (t1 > 0.001f) t = t1;
        if (t2 > 0.001f && t2 < t) t = t2;
        if (t > 1e29f) return info;

        info.hit = true;
        info.t = t;
        info.point = ray.at(t);
        info.normal = (info.point - center).normalize();
        info.mat = mat;

        float theta = acosf(-info.normal.y);
        float phi = atan2f(info.normal.z, info.normal.x);
        info.u = (phi + 3.1415926535f) / (2 * 3.1415926535f);
        info.v = theta / 3.1415926535f;

        return info;
    }
};

// =============== Треугольник ===============
struct Triangle : public Shape {
    Vec3 v0, v1, v2, normal;
    Vec2 uv0, uv1, uv2;
    Material mat;

    Triangle(Vec3 a, Vec3 b, Vec3 c, Vec2 a_uv, Vec2 b_uv, Vec2 c_uv, Material m, bool en = true)
        : v0(a), v1(b), v2(c), uv0(a_uv), uv1(b_uv), uv2(c_uv), mat(m) {
        normal = (v1 - v0).cross(v2 - v0).normalize();
        enabled = en;
    }

    HitInfo intersect(const Ray& ray) const override {
        if (!isEnabled()) return HitInfo();
        HitInfo info;
        Vec3 edge1 = v1 - v0, edge2 = v2 - v0;
        Vec3 h = ray.direction.cross(edge2);
        float a = edge1.dot(h);
        if (fabs(a) < 1e-6f) return info;
        float f = 1.0f / a;
        Vec3 s = ray.origin - v0;
        float u = f * s.dot(h);
        if (u < 0 || u > 1) return info;
        Vec3 q = s.cross(edge1);
        float v = f * ray.direction.dot(q);
        if (v < 0 || u + v > 1) return info;
        float t = f * edge2.dot(q);
        if (t < 0.001f) return info;

        info.hit = true;
        info.t = t;
        info.point = ray.at(t);
        info.normal = normal;
        if (info.normal.dot(ray.direction) > 0) info.normal = info.normal * -1.0f;
        info.mat = mat;

        float w = 1.0f - u - v;
        info.u = w * uv0.x + u * uv1.x + v * uv2.x;
        info.v = w * uv0.y + u * uv1.y + v * uv2.y;

        return info;
    }
};

// =============== Тетраэдр ===============
struct Tetrahedron : public Shape {
    std::vector<Triangle> tris;
    Tetrahedron(Vec3 a, Vec3 b, Vec3 c, Vec3 d, Material m, bool en = true) {
        enabled = en;
        tris.emplace_back(a, b, c, Vec2(0, 0), Vec2(1, 0), Vec2(0, 1), m, en);
        tris.emplace_back(a, c, d, Vec2(0, 0), Vec2(0, 1), Vec2(1, 1), m, en);
        tris.emplace_back(a, d, b, Vec2(0, 0), Vec2(1, 1), Vec2(1, 0), m, en);
        tris.emplace_back(b, d, c, Vec2(1, 0), Vec2(1, 1), Vec2(0, 1), m, en);
    }
    HitInfo intersect(const Ray& ray) const override {
        if (!isEnabled()) return HitInfo();
        HitInfo best;
        for (const auto& tri : tris) {
            HitInfo h = tri.intersect(ray);
            if (h.hit && h.t < best.t) best = h;
        }
        return best;
    }
};

// =============== Плоскость ===============
struct Plane : public Shape {
    Vec3 point, normal; Material mat;
    Plane(Vec3 p, Vec3 n, Material m, bool en = true) : point(p), normal(n.normalize()), mat(m) { enabled = en; }

    HitInfo intersect(const Ray& ray) const override {
        if (!isEnabled()) return HitInfo();
        HitInfo info;
        float denom = normal.dot(ray.direction);
        if (fabs(denom) < 1e-6f) return info;
        float t = (point - ray.origin).dot(normal) / denom;
        if (t < 0.001f) return info;

        info.hit = true;
        info.t = t;
        info.point = ray.at(t);
        info.normal = denom > 0 ? normal * -1.0f : normal;
        info.mat = mat;

        info.u = info.point.x * 0.5f;
        info.v = info.point.z * 0.5f;

        return info;
    }
};

// =============== Глобальные переменные ===============
const int RENDER_WIDTH = 1920;
const int RENDER_HEIGHT = 1080;
const int DISPLAY_WIDTH = 1920;
const int DISPLAY_HEIGHT = 1080;
unsigned char* image = nullptr;
GLuint textureID = 0;
Vec3 camera(0, 1, 5);
Vec3 lightPos(2, 5, 2);
std::vector<Shape*> scene;
const int MAX_DEPTH = 4; // чуть глубже для лучшего RTX-эффекта
const float FOV = 60.0f;

// =============== Трассировка ===============
Color traceRay(const Ray& ray, int depth = 0) {
    if (depth > MAX_DEPTH) return Color(0.05f, 0.05f, 0.1f);

    HitInfo closest;
    for (const auto& s : scene) {
        if (!s->isEnabled()) continue;
        HitInfo h = s->intersect(ray);
        if (h.hit && h.t < closest.t) closest = h;
    }

    if (!closest.hit) return Color(0.05f, 0.05f, 0.05f);

    Color baseColor = closest.mat.color;
    if (closest.mat.texture) {
        baseColor = closest.mat.texture->sample(closest.u, closest.v);
    }

    Vec3 L = (lightPos - closest.point).normalize();
    Vec3 V = (camera - closest.point).normalize();
    Vec3 R = L.reflect(closest.normal);

    bool inShadow = false;
    Ray shadowRay(closest.point + closest.normal * 0.001f, L);
    float lightDist = (lightPos - closest.point).length();
    for (const auto& s : scene) {
        if (!s->isEnabled()) continue;
        HitInfo h = s->intersect(shadowRay);
        if (h.hit && h.t < lightDist - 0.001f) { inShadow = true; break; }
    }

    Color ambient = baseColor * closest.mat.ka;
    float diff = std::max(0.0f, closest.normal.dot(L));
    Color diffuse = baseColor * closest.mat.kd * diff;
    float spec = powf(std::max(0.0f, V.dot(R)), (float)closest.mat.shininess);
    Color specular = Color(1, 1, 1) * closest.mat.ks * spec;

    Color result = ambient;
    if (!inShadow) result = result + diffuse + specular;

    if (closest.mat.reflect > 0 && depth < MAX_DEPTH) {
        Vec3 reflDir = ray.direction.reflect(closest.normal);
        Ray reflRay(closest.point + closest.normal * 0.001f, reflDir);
        Color reflCol = traceRay(reflRay, depth + 1);
        result = result.blend(reflCol, closest.mat.reflect);
    }

    return result;
}

// =============== Рендер с SSAA ===============
void renderImage() {
    const int NUM_THREADS = std::thread::hardware_concurrency();
    const int THREADS = (NUM_THREADS > 0) ? NUM_THREADS : 8;

    float tanHalfFov = tanf(FOV * 0.5f * 3.1415926535f / 180.0f);
    float aspect = (float)RENDER_WIDTH / RENDER_HEIGHT;

    std::vector<std::thread> threads;
    threads.reserve(THREADS);

    int rowsPerThread = RENDER_HEIGHT / THREADS;
    for (int t = 0; t < THREADS; t++) {
        int startRow = t * rowsPerThread;
        int endRow = (t == THREADS - 1) ? RENDER_HEIGHT : (t + 1) * rowsPerThread;

        // Захватываем ТОЛЬКО то, что нужно — по значению
        threads.emplace_back([startRow, endRow, tanHalfFov, aspect]() {
            for (int y = startRow; y < endRow; y++) {
                for (int x = 0; x < RENDER_WIDTH; x++) {
                    Color finalColor(0, 0, 0);

                    // SSAA 2x2
                    for (int sy = 0; sy < 2; sy++) {
                        for (int sx = 0; sx < 2; sx++) {
                            float px = (2.0f * (x + (sx + 0.5f) / 2.0f) / RENDER_WIDTH - 1.0f) * aspect * tanHalfFov;
                            float py = (1.0f - 2.0f * (y + (sy + 0.5f) / 2.0f) / RENDER_HEIGHT) * tanHalfFov;
                            Vec3 rayDir(px, py, -1);
                            Ray ray(camera, rayDir);
                            finalColor = finalColor + traceRay(ray);
                        }
                    }
                    finalColor = finalColor * 0.25f;

                    int idx = ((RENDER_HEIGHT - 1 - y) * RENDER_WIDTH + x) * 3;
                    image[idx] = finalColor.getR();
                    image[idx + 1] = finalColor.getG();
                    image[idx + 2] = finalColor.getB();
                }
            }
            });
    }

    for (auto& t : threads) {
        if (t.joinable()) t.join();
    }

    glBindTexture(GL_TEXTURE_2D, textureID);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, RENDER_WIDTH, RENDER_HEIGHT, GL_RGB, GL_UNSIGNED_BYTE, image);
}
// =============== Управление — 3D камера + свет ===============
void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    if (action != GLFW_PRESS) return;
    float step = 0.3f;
    bool changed = false;

    switch (key) {
        // Камера: 3D перемещение
    case GLFW_KEY_W: camera.z -= step; changed = true; break; // вперёд
    case GLFW_KEY_S: camera.z += step; changed = true; break; // назад
    case GLFW_KEY_A: camera.x -= step; changed = true; break; // лево
    case GLFW_KEY_D: camera.x += step; changed = true; break; // право
    case GLFW_KEY_Q: camera.y += step; changed = true; break; // вверх
    case GLFW_KEY_E: camera.y -= step; changed = true; break; // вниз

        // Свет
    case GLFW_KEY_LEFT:  lightPos.x -= step; changed = true; break;
    case GLFW_KEY_RIGHT: lightPos.x += step; changed = true; break;
    case GLFW_KEY_UP:    lightPos.y += step; changed = true; break;
    case GLFW_KEY_DOWN:  lightPos.y -= step; changed = true; break;

    case GLFW_KEY_R: changed = true; break;
    case GLFW_KEY_ESCAPE: glfwSetWindowShouldClose(window, GLFW_TRUE); return;
    }

    if (changed) {
        printf("Camera: (%.2f, %.2f, %.2f) | Light: (%.2f, %.2f, %.2f)\n",
            camera.x, camera.y, camera.z, lightPos.x, lightPos.y, lightPos.z);
        renderImage();
    }
}

// =============== main ===============
int main() {
    if (!glfwInit()) return -1;

    GLFWwindow* window = glfwCreateWindow(DISPLAY_WIDTH, DISPLAY_HEIGHT, "Лабораторная №4 — RTX-стиль", nullptr, nullptr);
    if (!window) { glfwTerminate(); return -1; }

    glfwMakeContextCurrent(window);
    glfwSetKeyCallback(window, keyCallback);

    // Загрузка текстур
    Texture sphereTex = Texture::loadBMP("sphere.bmp");
    Texture floorTex = Texture::loadBMP("floor.bmp");

    printf("Текстуры: sphere=%s, floor=%s\n",
        sphereTex.data ? "OK" : "FAIL",
        floorTex.data ? "OK" : "FAIL");

    // Сцена — Вариант 2
    scene.push_back(new Sphere(
        Vec3(0, 1, 0), 1.0f,
        Material(Color(1, 0, 0), 0.05f, 0.6f, 0.4f, 0.35f, 128, sphereTex.data ? &sphereTex : nullptr),
        true
    ));

    //scene.push_back(new Sphere(
    //    Vec3(0, 1, 0), 1.0f,
    //    Material(Color(1, 0, 0), 0.05f, 0.8f, 0.3f, 0.2f, 128, sphereTex.data ? &sphereTex : nullptr),
    //    true
    //));


    scene.push_back(new Tetrahedron(
        Vec3(2, 0, 0), Vec3(3, 0, 0), Vec3(2, 1, 0), Vec3(2, 0, 1),
        Material(Color(0, 0, 1), 0.05f, 0.6f, 0.4f, 0.4f, 128, nullptr),
        true
    ));

    // ЗЕРКАЛЬНЫЙ ПОЛ — сочный, но с текстурой
    scene.push_back(new Plane(
        Vec3(0, -1, 0), Vec3(0, 1, 0),
        Material(
            Color(0.95f, 0.95f, 0.95f), // fallback
            0.05f,   // ka — фоновое освещение
            0.75f,   // kd — ↑↑↑ ПОВЫШЕН, чтобы текстура "светилась"
            0.25f,   // ks — зеркальный блик
            0.25f,   // reflect — ↓↓↓ СНИЖЕН, чтобы не перекрывать текстуру
            128,     // shininess
            floorTex.data ? &floorTex : nullptr
        ),
        true
    ));

    // Инициализация
    image = new unsigned char[RENDER_WIDTH * RENDER_HEIGHT * 3];
    glGenTextures(1, &textureID);
    glBindTexture(GL_TEXTURE_2D, textureID);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, RENDER_WIDTH, RENDER_HEIGHT, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);

    printf("Рендер запущен (1024x1024 + SSAA)...\n");
    renderImage();

    // Цикл
    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(0, RENDER_WIDTH, 0, RENDER_HEIGHT, -1, 1);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, textureID);
        glBegin(GL_QUADS);
        glTexCoord2f(0, 0); glVertex2i(0, 0);
        glTexCoord2f(1, 0); glVertex2i(RENDER_WIDTH, 0);
        glTexCoord2f(1, 1); glVertex2i(RENDER_WIDTH, RENDER_HEIGHT);
        glTexCoord2f(0, 1); glVertex2i(0, RENDER_HEIGHT);
        glEnd();
        glDisable(GL_TEXTURE_2D);

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // Очистка
    for (auto s : scene) delete s;
    delete[] image;
    glDeleteTextures(1, &textureID);
    glfwTerminate();
    return 0;
}