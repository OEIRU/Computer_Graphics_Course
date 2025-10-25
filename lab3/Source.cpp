#define _CRT_SECURE_NO_WARNINGS
#define GLFW_EXPOSE_NATIVE_WIN32
#define NOMINMAX

#pragma comment(lib, "opengl32.lib")

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <shellapi.h>

#include <thread>
#include <vector>
#include <memory>
#include <string>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <math.h>
#include <chrono>  
#include <algorithm>
#include <locale.h>
#include <GLFW/glfw3.h>
#include <GLFW/glfw3native.h>

// =============== Вспомогательные структуры ===============

// 2D вектор (используется для UV-координат текстур)
struct Vec2 {
    float x, y;
    Vec2(float x = 0, float y = 0) : x(x), y(y) {}
};

// 3D вектор — основа для геометрии и лучей
struct Vec3 {
    float x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(float x, float y, float z) : x(x), y(y), z(z) {}

    // Операторы для удобной арифметики с векторами
    Vec3 operator+(const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
    Vec3 operator-(const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
    Vec3 operator-() const { return Vec3(-x, -y, -z); }
    Vec3 operator*(float t) const { return Vec3(x * t, y * t, z * t); }
    Vec3 operator/(float t) const { return t != 0 ? Vec3(x / t, y / t, z / t) : Vec3(0, 0, 0); }

    // Скалярное произведение (для вычисления углов)
    float dot(const Vec3& v) const { return x * v.x + y * v.y + z * v.z; }
    // Векторное произведение (для нормалей)
    Vec3 cross(const Vec3& v) const {
        return Vec3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }

    // Квадрат длины
    float length2() const { return dot(*this); }
    // Длина вектора
    float length() const { return sqrtf(length2()); }
    // Нормализация (приведение к единичной длине)
    Vec3 normalize() const {
        float len = length();
        return len > 1e-6f ? (*this) * (1.0f / len) : Vec3(0, 0, 1);
    }

    // Отражение вектора относительно нормали
    Vec3 reflect(const Vec3& n) const {
        return *this - n * (2.0f * this->dot(n));
    }

    // Преломление (рефракция) вектора при переходе между средами
    Vec3 refract(const Vec3& n, float eta) const {
        float cosI = -n.dot(*this);
        float sinT2 = eta * eta * (1.0f - cosI * cosI);
        if (sinT2 >= 1.0f) return Vec3(0, 0, 0); // Полное внутреннее отражение
        float cosT = sqrtf(1.0f - sinT2);
        return (*this) * eta + n * (eta * cosI - cosT);
    }
};

// Цвет в формате RGB с компонентами от 0.0 до 1.0
struct Color {
    float r, g, b;
    Color() : r(0), g(0), b(0) {}
    Color(float r, float g, float b) {
        // Ограничиваем значения в диапазоне [0, 1]
        this->r = std::clamp(r, 0.0f, 1.0f);
        this->g = std::clamp(g, 0.0f, 1.0f);
        this->b = std::clamp(b, 0.0f, 1.0f);
    }

    // Арифметика цветов
    Color operator*(float t) const { return Color(r * t, g * t, b * t); }
    Color operator+(const Color& c) const { return Color(r + c.r, g + c.g, b + c.b); }
    // Линейная интерполяция между двумя цветами
    Color blend(const Color& c, float t) const {
        t = std::clamp(t, 0.0f, 1.0f);
        return Color(r * (1 - t) + c.r * t, g * (1 - t) + c.g * t, b * (1 - t) + c.b * t);
    }

    // Преобразование в 8-битные значения для вывода в изображение
    unsigned char getR() const { return (unsigned char)(r * 255); }
    unsigned char getG() const { return (unsigned char)(g * 255); }
    unsigned char getB() const { return (unsigned char)(b * 255); }
};

// =============== Загрузчик BMP-текстур ===============
struct Texture {
    unsigned char* data = nullptr; // Пиксельные данные (BGR, 24 бита)
    int width = 0, height = 0;

    // Деструктор освобождает память
    ~Texture() { delete[] data; }

    // Получение цвета по UV-координатам (с повторением текстуры)
    Color sample(float u, float v) const {
        if (!data || width <= 0 || height <= 0) return Color(1, 0, 1); // Магента = ошибка
        u = u - floorf(u); // Повторение по горизонтали
        v = v - floorf(v); // Повторение по вертикали
        int x = (int)(u * width) % width;
        int y = (int)(v * height) % height;
        if (x < 0) x += width;
        if (y < 0) y += height;
        int idx = (y * width + x) * 3;
        // BMP хранит BGR, поэтому меняем порядок
        return Color(data[idx + 2] / 255.0f, data[idx + 1] / 255.0f, data[idx] / 255.0f);
    }

    // Загрузка BMP-файла (без поддержки сжатия, только 24-bit)
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

// =============== Материал объекта ===============
struct Material {
    Color color;           // Базовый цвет
    float ka, kd, ks;      // Коэффициенты: ambient, diffuse, specular
    float reflect;         // Отражательная способность (0..1)
    float transparency;    // Прозрачность (0..1)
    int shininess;         // Блеск (для specular)
    Texture* texture;      // Текстура (может быть nullptr)

    Material() : color(1, 1, 1), ka(0.1f), kd(0.7f), ks(0.3f), reflect(0.0f), transparency(0.0f), shininess(32), texture(nullptr) {}
    Material(Color c, float a, float d, float s, float r, float t, int sh, Texture* tex = nullptr)
        : color(c), ka(a), kd(d), ks(s), reflect(r), transparency(t), shininess(sh), texture(tex) {
    }
};

// =============== Луч и информация о пересечении ===============
struct Ray {
    Vec3 origin, direction;
    Ray(Vec3 o, Vec3 d) : origin(o), direction(d.normalize()) {}
    Vec3 at(float t) const { return origin + direction * t; } // Точка на луче
};

struct HitInfo {
    bool hit = false;      // Было ли пересечение
    float t = 1e30f;       // Расстояние до пересечения
    Vec3 point, normal;    // Точка и нормаль в точке пересечения
    Material mat;          // Материал объекта
    float u = 0, v = 0;    // UV-координаты для текстур
    HitInfo() {}
};

// =============== Базовый класс для всех фигур ===============
struct Shape {
    bool enabled = true;   // Включен ли объект в рендер
    virtual ~Shape() = default;
    virtual HitInfo intersect(const Ray& ray) const = 0; // Чисто виртуальный метод пересечения
    virtual bool isEnabled() const { return enabled; }
};

// =============== Сфера ===============
struct Sphere : public Shape {
    Vec3 center; float radius; Material mat;
    Sphere(Vec3 c, float r, Material m, bool en = true) : center(c), radius(r), mat(m) { enabled = en; }

    HitInfo intersect(const Ray& ray) const override {
        if (!isEnabled()) return HitInfo();
        Vec3 oc = ray.origin - center;
        float a = 1.0f;
        float b = 2.0f * oc.dot(ray.direction);
        float c = oc.dot(oc) - radius * radius;
        float disc = b * b - 4 * a * c;
        if (disc < 0) return HitInfo(); // Нет пересечения

        float sqrt_disc = sqrtf(disc);
        float t = (-b - sqrt_disc) / (2 * a);
        if (t < 0.001f) { // Проверяем дальнюю точку, если ближняя за камерой
            t = (-b + sqrt_disc) / (2 * a);
            if (t < 0.001f) return HitInfo();
        }

        HitInfo info;
        info.hit = true;
        info.t = t;
        info.point = ray.at(t);
        info.normal = (info.point - center).normalize();
        info.mat = mat;

        // Преобразуем нормаль в UV-координаты (для текстурирования сферы)
        float theta = acosf(std::clamp(-info.normal.y, -1.0f, 1.0f));
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
        Vec3 edge1 = v1 - v0, edge2 = v2 - v0;
        Vec3 h = ray.direction.cross(edge2);
        float a = edge1.dot(h);
        if (fabs(a) < 1e-6f) return HitInfo(); // Луч параллелен плоскости
        float f = 1.0f / a;
        Vec3 s = ray.origin - v0;
        float u = f * s.dot(h);
        if (u < 0 || u > 1) return HitInfo();
        Vec3 q = s.cross(edge1);
        float v = f * ray.direction.dot(q);
        if (v < 0 || u + v > 1) return HitInfo();
        float t = f * edge2.dot(q);
        if (t < 0.001f) return HitInfo();

        HitInfo info;
        info.hit = true;
        info.t = t;
        info.point = ray.at(t);
        info.normal = normal;
        if (info.normal.dot(ray.direction) > 0) info.normal = -info.normal; // Коррекция направления нормали
        info.mat = mat;

        // Интерполяция UV по барицентрическим координатам
        float w = 1.0f - u - v;
        info.u = w * uv0.x + u * uv1.x + v * uv2.x;
        info.v = w * uv0.y + u * uv1.y + v * uv2.y;
        return info;
    }
};

// =============== Тетраэдр (из 4 треугольников) ===============
struct Tetrahedron : public Shape {
    Vec3 v[4];
    Material mat;
    Tetrahedron(Vec3 a, Vec3 b, Vec3 c, Vec3 d, Material m, bool en = true)
        : v{ a, b, c, d }, mat(m) {
        enabled = en;
    }

    HitInfo intersect(const Ray& ray) const override {
        if (!isEnabled()) return HitInfo();
        HitInfo best;
        // Создаём 4 грани тетраэдра
        Triangle faces[4] = {
            Triangle(v[0], v[1], v[2], Vec2(0,0), Vec2(1,0), Vec2(0,1), mat, true),
            Triangle(v[0], v[2], v[3], Vec2(0,0), Vec2(0,1), Vec2(1,1), mat, true),
            Triangle(v[0], v[3], v[1], Vec2(0,0), Vec2(1,1), Vec2(1,0), mat, true),
            Triangle(v[1], v[3], v[2], Vec2(1,0), Vec2(1,1), Vec2(0,1), mat, true)
        };
        for (const auto& tri : faces) {
            HitInfo h = tri.intersect(ray);
            if (h.hit && h.t < best.t) best = h;
        }
        return best;
    }
};

// =============== Бесконечная плоскость ===============
struct Plane : public Shape {
    Vec3 point, normal; Material mat;
    Plane(Vec3 p, Vec3 n, Material m, bool en = true) : point(p), normal(n.normalize()), mat(m) { enabled = en; }

    HitInfo intersect(const Ray& ray) const override {
        if (!isEnabled()) return HitInfo();
        float denom = normal.dot(ray.direction);
        if (fabs(denom) < 1e-6f) return HitInfo(); // Луч параллелен плоскости
        float t = (point - ray.origin).dot(normal) / denom;
        if (t < 0.001f) return HitInfo();

        HitInfo info;
        info.hit = true;
        info.t = t;
        info.point = ray.at(t);
        info.normal = denom > 0 ? -normal : normal; // Коррекция нормали
        info.mat = mat;
        // Простая UV-проекция по координатам X и Z
        info.u = info.point.x * 2.0f;
        info.v = info.point.z * 2.0f;
        return info;
    }
};

// =============== Глобальные переменные ===============
const int RENDER_WIDTH = 1920;
const int RENDER_HEIGHT = 1080;
const int DISPLAY_WIDTH = 1920;
const int DISPLAY_HEIGHT = 1080;
unsigned char* image = nullptr; // Пиксельный буфер для рендера
GLuint textureID = 0;           // Идентификатор текстуры OpenGL

Vec3 camera(0, 1, 5);           // Позиция камеры
Vec3 lightPos(2, 5, 2);         // Позиция точечного источника света
bool useSSAA = false;           // Включено ли суперсэмплинговое сглаживание
std::vector<std::unique_ptr<Shape>> scene; // Список объектов в сцене
const int MAX_DEPTH = 5;        // Максимальная глубина рекурсии (отражения/преломления)
const float FOV = 60.0f;        // Угол обзора камеры

// Параметры рендеринга (предвычислены для ускорения)
struct RenderParams {
    float invWidth, invHeight, aspect, tanHalfFov;
    void update(int w, int h, float fov) {
        invWidth = 1.0f / w;
        invHeight = 1.0f / h;
        aspect = (float)w / h;
        tanHalfFov = tanf(fov * 0.5f * 3.1415926535f / 180.0f);
    }
} g_renderParams;

// =============== Отладочная информация ===============
void logDebugInfo(const std::string& message) {
    printf("[DEBUG] %s\n", message.c_str());
}

void logCameraInfo() {
    printf("[CAMERA] Pos: (%.2f, %.2f, %.2f)\n", camera.x, camera.y, camera.z);
}

void logLightInfo() {
    printf("[LIGHT] Pos: (%.2f, %.2f, %.2f)\n", lightPos.x, lightPos.y, lightPos.z);
}

void logSceneInfo() {
    printf("[SCENE] Objects: %zu\n", scene.size());
    for (size_t i = 0; i < scene.size(); ++i) {
        printf("  Objects %zu: %s\n", i, scene[i]->enabled ? "ON" : "OFF");
    }
}

void logMaterialInfo(const Material& mat, const std::string& name) {
    printf("[Material %s] ka=%.2f, kd=%.2f, ks=%.2f, reflection=%.2f, transparency=%.2f, blind=%d\n",
        name.c_str(), mat.ka, mat.kd, mat.ks, mat.reflect, mat.transparency, mat.shininess);
}

void logRenderInfo() {
    printf("[RENDER] SSAA: %s\n", useSSAA ? "ON" : "OFF");
}

// =============== Физически корректное освещение: Cook-Torrance ===============
Color cookTorrance(const Vec3& N, const Vec3& V, const Vec3& L, const Color& baseColor, float ks, int shininess) {
    Vec3 H = (L + V).normalize(); // Половинный вектор
    float NdotL = std::max(0.0f, N.dot(L));
    float NdotV = std::max(0.0f, N.dot(V));
    float NdotH = std::max(0.0f, N.dot(H));
    float VdotH = std::max(0.0f, V.dot(H));

    if (NdotL == 0 || NdotV == 0) return Color(0, 0, 0);

    // Френель: отражение зависит от угла
    float F0 = 0.04f;
    float F = F0 + (1.0f - F0) * powf(1.0f - VdotH, 5.0f);

    // Шероховатость из shininess
    float alpha = std::clamp(1.0f - (float)shininess / 256.0f, 0.01f, 1.0f);
    float alpha2 = alpha * alpha;

    // Распределение микрограней (GGX)
    float denom = NdotH * NdotH * (alpha2 - 1.0f) + 1.0f;
    float D = alpha2 / (3.1415926535f * denom * denom);

    // Геометрическое затенение
    float k = (alpha + 2.0f) / 8.0f;
    float G1V = NdotV / (NdotV * (1.0f - k) + k);
    float G1L = NdotL / (NdotL * (1.0f - k) + k);
    float G = G1V * G1L;

    // Итоговый коэффициент зеркального отражения
    float denomSpec = 4.0f * NdotV * NdotL;
    if (denomSpec < 1e-6f) return Color(0, 0, 0);

    float specular = F * D * G / denomSpec;
    return Color(1, 1, 1) * ks * specular;
}

// =============== Трассировка луча с рекурсией ===============
Color traceRay(const Ray& ray, int depth = 0) {
    if (depth >= MAX_DEPTH) return Color(0.05f, 0.05f, 0.1f); // Фон

    HitInfo closest;
    for (const auto& s : scene) {
        if (!s->isEnabled()) continue;
        HitInfo h = s->intersect(ray);
        if (h.hit && h.t < closest.t) closest = h;
    }

    if (!closest.hit) return Color(0.05f, 0.05f, 0.05f); // Фон

    // Получаем базовый цвет (с текстурой или без)
    Color baseColor = closest.mat.color;
    if (closest.mat.texture) {
        baseColor = closest.mat.texture->sample(closest.u, closest.v);
    }

    // Направление к источнику света
    Vec3 L_vec = lightPos - closest.point;
    float lightDist2 = L_vec.length2();
    if (lightDist2 < 1e-6f) {
        return baseColor * closest.mat.ka; // Только ambient
    }
    float lightDist = sqrtf(lightDist2);
    Vec3 L = L_vec / lightDist;
    Vec3 V = (camera - closest.point).normalize();
    Vec3 N = closest.normal;

    // Проверка тени: бросаем луч к источнику света
    bool inShadow = false;
    Ray shadowRay(closest.point + N * 0.001f, L);
    for (const auto& s : scene) {
        if (!s->isEnabled()) continue;
        HitInfo h = s->intersect(shadowRay);
        if (h.hit && h.t * h.t < lightDist2 - 0.01f) {
            inShadow = true;
            break;
        }
    }

    // Ambient — всегда есть
    Color ambient = baseColor * closest.mat.ka;
    Color result = ambient;

    if (!inShadow) {
        float diff = std::max(0.0f, N.dot(L));
        Color diffuse = baseColor * closest.mat.kd * diff;
        Color specular = cookTorrance(N, V, L, baseColor, closest.mat.ks, closest.mat.shininess);
        result = result + diffuse + specular;
    }

    // Рекурсивное отражение
    if (closest.mat.reflect > 0 && depth < MAX_DEPTH) {
        Vec3 reflDir = ray.direction.reflect(N);
        Ray reflRay(closest.point + N * 0.001f, reflDir);
        Color reflCol = traceRay(reflRay, depth + 1);
        result = result.blend(reflCol, closest.mat.reflect);
    }

    // Рекурсивное преломление
    if (closest.mat.transparency > 0 && depth < MAX_DEPTH) {
        float eta = 1.0f / 1.5f; // Воздух -> стекло
        Vec3 N_refr = N;
        if (N.dot(ray.direction) > 0) {
            N_refr = -N; // Луч выходит из объекта
            eta = 1.5f;  // Стекло -> воздух
        }
        Vec3 refrDir = ray.direction.refract(N_refr, eta);
        if (refrDir.length2() > 0) {
            Ray refrRay(closest.point - N_refr * 0.001f, refrDir);
            Color refrCol = traceRay(refrRay, depth + 1);
            result = result * (1.0f - closest.mat.transparency) + refrCol * closest.mat.transparency;
        }
    }

    return result;
}

// =============== Многопоточный рендер изображения ===============
void renderImage() {
    logDebugInfo("Start rendering");
    auto startTime = std::chrono::high_resolution_clock::now();

    g_renderParams.update(RENDER_WIDTH, RENDER_HEIGHT, FOV);

    const int NUM_THREADS = std::thread::hardware_concurrency();
    const int THREADS = (NUM_THREADS > 0) ? NUM_THREADS : 8;

    std::vector<std::thread> threads;
    threads.reserve(THREADS);

    int rowsPerThread = RENDER_HEIGHT / THREADS;
    const bool ssaa = useSSAA;

    for (int t = 0; t < THREADS; ++t) {
        int startRow = t * rowsPerThread;
        int endRow = (t == THREADS - 1) ? RENDER_HEIGHT : (t + 1) * rowsPerThread;

        threads.emplace_back([startRow, endRow, ssaa]() {
            const float& invWidth = g_renderParams.invWidth;
            const float& invHeight = g_renderParams.invHeight;
            const float& aspect = g_renderParams.aspect;
            const float& tanHalfFov = g_renderParams.tanHalfFov;

            for (int y = startRow; y < endRow; ++y) {
                for (int x = 0; x < RENDER_WIDTH; ++x) {
                    Color finalColor(0, 0, 0);
                    int samples = ssaa ? 4 : 1;

                    // SSAA: 2x2 субпиксельных сэмпла
                    for (int sy = 0; sy < (ssaa ? 2 : 1); ++sy) {
                        for (int sx = 0; sx < (ssaa ? 2 : 1); ++sx) {
                            float jitterX = ssaa ? (sx + 0.5f) / 2.0f : 0.5f;
                            float jitterY = ssaa ? (sy + 0.5f) / 2.0f : 0.5f;
                            float px = (2.0f * (x + jitterX) * invWidth - 1.0f) * aspect * tanHalfFov;
                            float py = (1.0f - 2.0f * (y + jitterY) * invHeight) * tanHalfFov;
                            Vec3 rayDir(px, py, -1);
                            Ray ray(camera, rayDir);
                            finalColor = finalColor + traceRay(ray);
                        }
                    }

                    finalColor = finalColor * (1.0f / samples);

                    // Запись в буфер (в формате RGB, с переворотом по Y)
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

    // Обновляем текстуру OpenGL
    glBindTexture(GL_TEXTURE_2D, textureID);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, RENDER_WIDTH, RENDER_HEIGHT, GL_RGB, GL_UNSIGNED_BYTE, image);

    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    logDebugInfo("Render completed Za" + std::to_string(duration.count()) + " ms");
}

// =============== Загрузка сцены из файла ===============
void loadSceneFromFile(const char* filename) {
    logDebugInfo("Load scene from: " + std::string(filename));
    scene.clear();
    std::ifstream file(filename);
    if (!file.is_open()) {
        logDebugInfo("ERROR!: File scene not found");
        return;
    }

    std::string line;
    int objectsLoaded = 0;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue; // Комментарии
        std::istringstream iss(line);
        std::string type;
        iss >> type;

        if (type == "sphere") {
            float cx, cy, cz, r, R, G, B, ka, kd, ks, refl, trans, shininess;
            std::string texName;
            iss >> cx >> cy >> cz >> r >> R >> G >> B >> ka >> kd >> ks >> refl >> trans >> shininess >> texName;
            Color col(R, G, B);
            Texture* tex = (texName != "none") ? new Texture(Texture::loadBMP(texName.c_str())) : nullptr;
            scene.push_back(std::make_unique<Sphere>(Vec3(cx, cy, cz), r, Material(col, ka, kd, ks, refl, trans, (int)shininess, tex)));
            objectsLoaded++;
            logDebugInfo("Sphere loaded in (" + std::to_string(cx) + ", " + std::to_string(cy) + ", " + std::to_string(cz) + ")");
        }
        else if (type == "tetra") {
            float coords[12], R, G, B, ka, kd, ks, refl, trans, shininess;
            std::string texName;
            for (int i = 0; i < 12; ++i) iss >> coords[i];
            iss >> R >> G >> B >> ka >> kd >> ks >> refl >> trans >> shininess >> texName;
            Color col(R, G, B);
            Texture* tex = (texName != "none") ? new Texture(Texture::loadBMP(texName.c_str())) : nullptr;
            scene.push_back(std::make_unique<Tetrahedron>(
                Vec3(coords[0], coords[1], coords[2]),
                Vec3(coords[3], coords[4], coords[5]),
                Vec3(coords[6], coords[7], coords[8]),
                Vec3(coords[9], coords[10], coords[11]),
                Material(col, ka, kd, ks, refl, trans, (int)shininess, tex)
            ));
            objectsLoaded++;
            logDebugInfo("Tetra loaded");
        }
        else if (type == "plane") {
            float px, py, pz, nx, ny, nz, R, G, B, ka, kd, ks, refl, trans, shininess;
            std::string texName;
            iss >> px >> py >> pz >> nx >> ny >> nz >> R >> G >> B >> ka >> kd >> ks >> refl >> trans >> shininess >> texName;
            Color col(R, G, B);
            Texture* tex = (texName != "none") ? new Texture(Texture::loadBMP(texName.c_str())) : nullptr;
            scene.push_back(std::make_unique<Plane>(Vec3(px, py, pz), Vec3(nx, ny, nz), Material(col, ka, kd, ks, refl, trans, (int)shininess, tex)));
            objectsLoaded++;
            logDebugInfo("Plane loaded (" + std::to_string(px) + ", " + std::to_string(py) + ", " + std::to_string(pz) + ")");
        }
    }
    logDebugInfo("Scene completed: " + std::to_string(objectsLoaded));
}

// =============== Действия меню ===============
enum MenuActions {
    // Камера
    MENU_CAM_FORWARD, MENU_CAM_BACK, MENU_CAM_LEFT, MENU_CAM_RIGHT, MENU_CAM_UP, MENU_CAM_DOWN, MENU_RESET_CAMERA,

    // Свет
    MENU_LIGHT_UP, MENU_LIGHT_DOWN, MENU_LIGHT_LEFT, MENU_LIGHT_RIGHT, MENU_RESET_LIGHT,

    // Объекты
    MENU_TOGGLE_SPHERE, MENU_TOGGLE_TETRA, MENU_TOGGLE_PLANE,

    // Материалы (сфера)
    MENU_SPHERE_KA_UP, MENU_SPHERE_KA_DOWN,
    MENU_SPHERE_KD_UP, MENU_SPHERE_KD_DOWN,
    MENU_SPHERE_KS_UP, MENU_SPHERE_KS_DOWN,
    MENU_SPHERE_REFLECT_UP, MENU_SPHERE_REFLECT_DOWN,
    MENU_SPHERE_TRANSP_UP, MENU_SPHERE_TRANSP_DOWN,

    // Рендер
    MENU_TOGGLE_SSAA, MENU_RELOAD_SCENE, MENU_RENDER_NOW
};

// Обновление параметров материала сферы
void updateMaterial(Material& mat, int action) {
    const float delta = 0.05f;
    switch (action) {
    case MENU_SPHERE_KA_UP:     mat.ka = std::min(mat.ka + delta, 1.0f); break;
    case MENU_SPHERE_KA_DOWN:   mat.ka = std::max(mat.ka - delta, 0.0f); break;
    case MENU_SPHERE_KD_UP:     mat.kd = std::min(mat.kd + delta, 1.0f); break;
    case MENU_SPHERE_KD_DOWN:   mat.kd = std::max(mat.kd - delta, 0.0f); break;
    case MENU_SPHERE_KS_UP:     mat.ks = std::min(mat.ks + delta, 1.0f); break;
    case MENU_SPHERE_KS_DOWN:   mat.ks = std::max(mat.ks - delta, 0.0f); break;
    case MENU_SPHERE_REFLECT_UP:    mat.reflect = std::min(mat.reflect + delta, 1.0f); break;
    case MENU_SPHERE_REFLECT_DOWN:  mat.reflect = std::max(mat.reflect - delta, 0.0f); break;
    case MENU_SPHERE_TRANSP_UP:     mat.transparency = std::min(mat.transparency + delta, 1.0f); break;
    case MENU_SPHERE_TRANSP_DOWN:   mat.transparency = std::max(mat.transparency - delta, 0.0f); break;
    }
}

// Выполнение действия из меню
void executeMenuAction(int action) {
    bool changed = false;

    switch (action) {
        // Камера
    case MENU_CAM_FORWARD:
        camera.z -= 0.3f;
        changed = true;
        logDebugInfo("Камера: вперёд");
        logCameraInfo();
        break;
    case MENU_CAM_BACK:
        camera.z += 0.3f;
        changed = true;
        logDebugInfo("Камера: назад");
        logCameraInfo();
        break;
    case MENU_CAM_LEFT:
        camera.x -= 0.3f;
        changed = true;
        logDebugInfo("Камера: влево");
        logCameraInfo();
        break;
    case MENU_CAM_RIGHT:
        camera.x += 0.3f;
        changed = true;
        logDebugInfo("Камера: вправо");
        logCameraInfo();
        break;
    case MENU_CAM_UP:
        camera.y += 0.3f;
        changed = true;
        logDebugInfo("Камера: вверх");
        logCameraInfo();
        break;
    case MENU_CAM_DOWN:
        camera.y -= 0.3f;
        changed = true;
        logDebugInfo("Камера: вниз");
        logCameraInfo();
        break;
    case MENU_RESET_CAMERA:
        camera = Vec3(0, 1, 5);
        changed = true;
        logDebugInfo("Камера: сброс");
        logCameraInfo();
        break;

        // Свет
    case MENU_LIGHT_UP:
        lightPos.y += 0.3f;
        changed = true;
        logDebugInfo("Свет: вверх");
        logLightInfo();
        break;
    case MENU_LIGHT_DOWN:
        lightPos.y -= 0.3f;
        changed = true;
        logDebugInfo("Свет: вниз");
        logLightInfo();
        break;
    case MENU_LIGHT_LEFT:
        lightPos.x -= 0.3f;
        changed = true;
        logDebugInfo("Свет: влево");
        logLightInfo();
        break;
    case MENU_LIGHT_RIGHT:
        lightPos.x += 0.3f;
        changed = true;
        logDebugInfo("Свет: вправо");
        logLightInfo();
        break;
    case MENU_RESET_LIGHT:
        lightPos = Vec3(2, 5, 2);
        changed = true;
        logDebugInfo("Свет: сброс");
        logLightInfo();
        break;

        // Объекты
    case MENU_TOGGLE_SPHERE:
        if (scene.size() >= 1) {
            scene[0]->enabled = !scene[0]->enabled;
            changed = true;
            logDebugInfo("Сфера: " + std::string(scene[0]->enabled ? "ON" : "OFF"));
            logSceneInfo();
        }
        break;
    case MENU_TOGGLE_TETRA:
        if (scene.size() >= 2) {
            scene[1]->enabled = !scene[1]->enabled;
            changed = true;
            logDebugInfo("Тетраэдр: " + std::string(scene[1]->enabled ? "ON" : "OFF"));
            logSceneInfo();
        }
        break;
    case MENU_TOGGLE_PLANE:
        if (scene.size() >= 3) {
            scene[2]->enabled = !scene[2]->enabled;
            changed = true;
            logDebugInfo("Плоскость: " + std::string(scene[2]->enabled ? "ON" : "OFF"));
            logSceneInfo();
        }
        break;

        // Материалы (сфера)
    case MENU_SPHERE_KA_UP: case MENU_SPHERE_KA_DOWN:
    case MENU_SPHERE_KD_UP: case MENU_SPHERE_KD_DOWN:
    case MENU_SPHERE_KS_UP: case MENU_SPHERE_KS_DOWN:
    case MENU_SPHERE_REFLECT_UP: case MENU_SPHERE_REFLECT_DOWN:
    case MENU_SPHERE_TRANSP_UP: case MENU_SPHERE_TRANSP_DOWN:
        if (scene.size() >= 1 && scene[0]->isEnabled()) {
            if (Sphere* s = dynamic_cast<Sphere*>(scene[0].get())) {
                updateMaterial(s->mat, action);
                changed = true;
                logDebugInfo("Материал сферы обновлён");
                logMaterialInfo(s->mat, "СФЕРА");
            }
        }
        break;

        // Рендер
    case MENU_TOGGLE_SSAA:
        useSSAA = !useSSAA;
        changed = true;
        logDebugInfo("SSAA: " + std::string(useSSAA ? "ON" : "OFF"));
        logRenderInfo();
        break;
    case MENU_RELOAD_SCENE:
        loadSceneFromFile("scene.txt");
        changed = true;
        break;
    case MENU_RENDER_NOW:
        changed = true;
        logDebugInfo("Принудительный рендер");
        break;
    }

    if (changed) {
        renderImage();
    }
}

// =============== Создание нативного Win32 контекстного меню (на русском!) ===============
HMENU g_hContextMenu = nullptr;

void createNativeMenu() {
    g_hContextMenu = CreatePopupMenu();

    // --- Камера ---
    HMENU hCam = CreatePopupMenu();
    AppendMenuA(hCam, MF_STRING, MENU_CAM_FORWARD, "Вперёд");
    AppendMenuA(hCam, MF_STRING, MENU_CAM_BACK, "Назад");
    AppendMenuA(hCam, MF_STRING, MENU_CAM_LEFT, "Влево");
    AppendMenuA(hCam, MF_STRING, MENU_CAM_RIGHT, "Вправо");
    AppendMenuA(hCam, MF_STRING, MENU_CAM_UP, "Вверх");
    AppendMenuA(hCam, MF_STRING, MENU_CAM_DOWN, "Вниз");
    AppendMenuA(hCam, MF_SEPARATOR, 0, nullptr);
    AppendMenuA(hCam, MF_STRING, MENU_RESET_CAMERA, "Сбросить камеру");
    AppendMenuA(g_hContextMenu, MF_POPUP, (UINT_PTR)hCam, "Камера");

    // --- Свет ---
    HMENU hLight = CreatePopupMenu();
    AppendMenuA(hLight, MF_STRING, MENU_LIGHT_UP, "Свет вверх");
    AppendMenuA(hLight, MF_STRING, MENU_LIGHT_DOWN, "Свет вниз");
    AppendMenuA(hLight, MF_STRING, MENU_LIGHT_LEFT, "Свет влево");
    AppendMenuA(hLight, MF_STRING, MENU_LIGHT_RIGHT, "Свет вправо");
    AppendMenuA(hLight, MF_STRING, MENU_RESET_LIGHT, "Сбросить свет");
    AppendMenuA(g_hContextMenu, MF_POPUP, (UINT_PTR)hLight, "Свет");

    // --- Объекты ---
    HMENU hObjects = CreatePopupMenu();
    AppendMenuA(hObjects, MF_STRING, MENU_TOGGLE_SPHERE, "Сфера: ON/OFF");
    AppendMenuA(hObjects, MF_STRING, MENU_TOGGLE_TETRA, "Тетраэдр: ON/OFF");
    AppendMenuA(hObjects, MF_STRING, MENU_TOGGLE_PLANE, "Пол: ON/OFF");
    AppendMenuA(g_hContextMenu, MF_POPUP, (UINT_PTR)hObjects, "Объекты");

    // --- Материалы (сфера) ---
    HMENU hMat = CreatePopupMenu();
    AppendMenuA(hMat, MF_STRING, MENU_SPHERE_KA_UP, "Фоновое освещение +");
    AppendMenuA(hMat, MF_STRING, MENU_SPHERE_KA_DOWN, "Фоновое освещение -");
    AppendMenuA(hMat, MF_STRING, MENU_SPHERE_KD_UP, "Диффузное +");
    AppendMenuA(hMat, MF_STRING, MENU_SPHERE_KD_DOWN, "Диффузное -");
    AppendMenuA(hMat, MF_STRING, MENU_SPHERE_KS_UP, "Зеркальное +");
    AppendMenuA(hMat, MF_STRING, MENU_SPHERE_KS_DOWN, "Зеркальное -");
    AppendMenuA(hMat, MF_STRING, MENU_SPHERE_REFLECT_UP, "Отражение +");
    AppendMenuA(hMat, MF_STRING, MENU_SPHERE_REFLECT_DOWN, "Отражение -");
    AppendMenuA(hMat, MF_STRING, MENU_SPHERE_TRANSP_UP, "Прозрачность +");
    AppendMenuA(hMat, MF_STRING, MENU_SPHERE_TRANSP_DOWN, "Прозрачность -");
    AppendMenuA(g_hContextMenu, MF_POPUP, (UINT_PTR)hMat, "Материалы (сфера)");

    // --- Рендер ---
    HMENU hRender = CreatePopupMenu();
    AppendMenuA(hRender, MF_STRING, MENU_TOGGLE_SSAA, "Сглаживание (SSAA)");
    AppendMenuA(hRender, MF_STRING, MENU_RELOAD_SCENE, "Перезагрузить сцену");
    AppendMenuA(hRender, MF_STRING, MENU_RENDER_NOW, "Перерисовать сейчас");
    AppendMenuA(g_hContextMenu, MF_POPUP, (UINT_PTR)hRender, "Рендер");

    // --- Выход ---
    AppendMenuA(g_hContextMenu, MF_SEPARATOR, 0, nullptr);
    AppendMenuA(g_hContextMenu, MF_STRING, 0x1000, "Выход"); // Не обрабатываем — Esc
}

// Показ контекстного меню по правому клику
void showNativeContextMenu(int x, int y, GLFWwindow* window) {
    if (!g_hContextMenu) return;

    HWND hwnd = glfwGetWin32Window(window);
    if (!hwnd) return;

    POINT pt = { (LONG)x, (LONG)y };
    ClientToScreen(hwnd, &pt);

    UINT cmd = TrackPopupMenu(
        g_hContextMenu,
        TPM_RETURNCMD | TPM_NONOTIFY,
        pt.x, pt.y,
        0,
        hwnd,
        nullptr
    );

    if (cmd != 0 && cmd != 0x1000) {
        executeMenuAction(cmd);
    }
}

// =============== Обработчики событий GLFW ===============
void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    if (action != GLFW_PRESS) return;

    switch (key) {
    case GLFW_KEY_ESCAPE:
        glfwSetWindowShouldClose(window, GLFW_TRUE);
        logDebugInfo("Нажат ESC — выход");
        break;

    case GLFW_KEY_F:
        useSSAA = !useSSAA;
        logDebugInfo("Нажата F — SSAA: " + std::string(useSSAA ? "ON" : "OFF"));
        logRenderInfo();
        renderImage();
        break;

    case GLFW_KEY_R:
        logDebugInfo("Нажата R — принудительный рендер");
        renderImage();
        break;

        // Управление камерой
    case GLFW_KEY_W: camera.z -= 0.3f; logDebugInfo("W — камера вперёд"); logCameraInfo(); renderImage(); break;
    case GLFW_KEY_S: camera.z += 0.3f; logDebugInfo("S — камера назад"); logCameraInfo(); renderImage(); break;
    case GLFW_KEY_A: camera.x -= 0.3f; logDebugInfo("A — камера влево"); logCameraInfo(); renderImage(); break;
    case GLFW_KEY_D: camera.x += 0.3f; logDebugInfo("D — камера вправо"); logCameraInfo(); renderImage(); break;
    case GLFW_KEY_Q: camera.y += 0.3f; logDebugInfo("Q — камера вверх"); logCameraInfo(); renderImage(); break;
    case GLFW_KEY_E: camera.y -= 0.3f; logDebugInfo("E — камера вниз"); logCameraInfo(); renderImage(); break;

    case GLFW_KEY_C:
        camera = Vec3(0, 1, 5);
        logDebugInfo("C — сброс камеры");
        logCameraInfo();
        renderImage();
        break;

        // Управление светом
    case GLFW_KEY_I: lightPos.y += 0.3f; logDebugInfo("I — свет вверх"); logLightInfo(); renderImage(); break;
    case GLFW_KEY_K: lightPos.y -= 0.3f; logDebugInfo("K — свет вниз"); logLightInfo(); renderImage(); break;
    case GLFW_KEY_J: lightPos.x -= 0.3f; logDebugInfo("J — свет влево"); logLightInfo(); renderImage(); break;
    case GLFW_KEY_L: lightPos.x += 0.3f; logDebugInfo("L — свет вправо"); logLightInfo(); renderImage(); break;
    case GLFW_KEY_O:
        lightPos = Vec3(2, 5, 2);
        logDebugInfo("O — сброс света");
        logLightInfo();
        renderImage();
        break;

        // Переключение объектов
    case GLFW_KEY_1:
        if (scene.size() >= 1) {
            scene[0]->enabled = !scene[0]->enabled;
            logDebugInfo("1 — сфера: " + std::string(scene[0]->enabled ? "ON" : "OFF"));
            logSceneInfo();
            renderImage();
        }
        break;
    case GLFW_KEY_2:
        if (scene.size() >= 2) {
            scene[1]->enabled = !scene[1]->enabled;
            logDebugInfo("2 — тетраэдр: " + std::string(scene[1]->enabled ? "ON" : "OFF"));
            logSceneInfo();
            renderImage();
        }
        break;
    case GLFW_KEY_3:
        if (scene.size() >= 3) {
            scene[2]->enabled = !scene[2]->enabled;
            logDebugInfo("3 — пол: " + std::string(scene[2]->enabled ? "ON" : "OFF"));
            logSceneInfo();
            renderImage();
        }
        break;

    case GLFW_KEY_F5:
        logDebugInfo("F5 — перезагрузка сцены");
        loadSceneFromFile("scene.txt");
        renderImage();
        break;
    }
}

void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
        double xpos, ypos;
        glfwGetCursorPos(window, &xpos, &ypos);
        showNativeContextMenu((int)xpos, (int)ypos, window);
    }
}

// =============== Точка входа ===============
int main() {
    setlocale(LC_ALL, "Russian");
    // Создаём консоль для отладочного вывода
    AllocConsole();
    freopen("CONOUT$", "w", stdout);
    freopen("CONOUT$", "w", stderr);

    printf("=== Code started===\n");
    printf("Move:\n");
    printf("  Camera: W/S/A/D/Q/E — move, C — restart\n");
    printf("  Light: I/K/J/L — move, O — restart\n");
    printf("  Objects: 1/2/3 — ON/OFF sphere/tetra/plane\n");
    printf("  Render: F — SSAA, R — re-render, F5 — re-render scene\n");
    printf("  Menu: Right mouse button\n");
    printf("  Exit: ESC\n\n");

    if (!glfwInit()) return -1;

    GLFWwindow* window = glfwCreateWindow(DISPLAY_WIDTH, DISPLAY_HEIGHT, "LR3", nullptr, nullptr);
    if (!window) { glfwTerminate(); return -1; }

    glfwMakeContextCurrent(window);
    glfwSetKeyCallback(window, keyCallback);
    glfwSetMouseButtonCallback(window, mouseButtonCallback);

    createNativeMenu();
    loadSceneFromFile("scene.txt");

    // Выделение памяти под изображение и создание OpenGL-текстуры
    image = new unsigned char[RENDER_WIDTH * RENDER_HEIGHT * 3]();
    glGenTextures(1, &textureID);
    glBindTexture(GL_TEXTURE_2D, textureID);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, RENDER_WIDTH, RENDER_HEIGHT, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);

    // Вывод начальной информации
    logCameraInfo();
    logLightInfo();
    logSceneInfo();
    logRenderInfo();

    renderImage();

    // Главный цикл отображения
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

    printf("=== Work finished ===\n");

    if (g_hContextMenu) DestroyMenu(g_hContextMenu);
    delete[] image;
    glDeleteTextures(1, &textureID);
    glfwTerminate();
    FreeConsole();
    return 0;
}