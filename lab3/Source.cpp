// compile
// g++ -O3 -std=c++17 -o main Source.cpp -lglfw -lGL -lpthread 

#define _CRT_SECURE_NO_WARNINGS

#include <thread>
#include <vector>
#include <memory>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <cstdlib>

#include <GL/gl.h>
#include <GLFW/glfw3.h>

// ============================================================================
// ВСПОМОГАТЕЛЬНЫЕ СТРУКТУРЫ
// ============================================================================

// 2D-вектор для UV-координат текстур
struct Vec2 {
    float x, y;
    Vec2(float x = 0, float y = 0) : x(x), y(y) {}
};

// 3D-вектор: основа для точек, направлений, нормалей и лучей
struct Vec3 {
    float x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(float x, float y, float z) : x(x), y(y), z(z) {}

    // Арифметические операции
    Vec3 operator+(const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
    Vec3 operator-(const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
    Vec3 operator-() const { return Vec3(-x, -y, -z); }
    Vec3 operator*(float t) const { return Vec3(x * t, y * t, z * t); }
    Vec3 operator/(float t) const { return t != 0 ? Vec3(x / t, y / t, z / t) : Vec3(0, 0, 0); }

    // Скалярное произведение (для углов и проекций)
    float dot(const Vec3& v) const { return x * v.x + y * v.y + z * v.z; }

    // Векторное произведение (для вычисления нормалей)
    Vec3 cross(const Vec3& v) const {
        return Vec3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }

    // Квадрат длины
    float length2() const { return dot(*this); }

    // Длина вектора
    float length() const { return std::sqrt(length2()); }

    // Нормализация: приведение к единичной длине
    Vec3 normalize() const {
        float len = length();
        return len > 1e-6f ? (*this) * (1.0f / len) : Vec3(0, 0, 1);
    }

    // Отражение вектора относительно нормали (для зеркал)
    Vec3 reflect(const Vec3& n) const {
        return *this - n * (2.0f * this->dot(n));
    }

    // Преломление (закон Снеллиуса) при переходе между средами
    Vec3 refract(const Vec3& n, float eta) const {
        float cosI = -n.dot(*this); // косинус угла падения
        float sinT2 = eta * eta * (1.0f - cosI * cosI); // sin² угла преломления
        if (sinT2 >= 1.0f) return Vec3(0, 0, 0); // Полное внутреннее отражение
        float cosT = std::sqrt(1.0f - sinT2);
        return (*this) * eta + n * (eta * cosI - cosT);
    }
};

// Цвет в формате RGB с компонентами от 0.0 до 1.0
struct Color {
    float r, g, b;
    Color() : r(0), g(0), b(0) {}
    Color(float r, float g, float b) {
        // Ограничение значений в [0, 1]
        this->r = std::clamp(r, 0.0f, 1.0f);
        this->g = std::clamp(g, 0.0f, 1.0f);
        this->b = std::clamp(b, 0.0f, 1.0f);
    }

    // Операции над цветами
    Color operator*(float t) const { return Color(r * t, g * t, b * t); }
    Color operator+(const Color& c) const { return Color(r + c.r, g + c.g, b + c.b); }
    Color operator*(const Color& c) const { return Color(r * c.r, g * c.g, b * c.b); }

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

// ============================================================================
// ЗАГРУЗЧИК BMP-ТЕКСТУР
// Поддерживает только 24-битные несжатые BMP
// ============================================================================

struct Texture {
    unsigned char* data = nullptr; // Пиксельные данные (BGR, 24 бита на пиксель)
    int width = 0, height = 0;     // Размеры изображения
    int imageSize = 0;             // Размер данных в байтах
    int dataPos = 0;               // Смещение начала пиксельных данных

    // Деструктор автоматически освобождает память
    ~Texture() { delete[] data; }

    // Получение цвета по UV-координатам с повторением текстуры (tiling)
    Color sample(float u, float v) const {
        if (!data || width <= 0 || height <= 0) return Color(1, 0, 1); // Магента = ошибка

        // Повторение текстуры: дробная часть UV
        u = u - std::floor(u);
        v = v - std::floor(v);

        // Преобразование UV в индексы пикселей
        int x = (int)(u * width) % width;
        int y = (int)(v * height) % height;
        if (x < 0) x += width;
        if (y < 0) y += height;

        // Индекс в массиве данных (BMP хранит BGR)
        int idx = (y * width + x) * 3;
        return Color(
            data[idx + 2] / 255.0f, // R
            data[idx + 1] / 255.0f, // G
            data[idx] / 255.0f      // B
        );
    }

    // Статический метод загрузки BMP-файла
    static Texture loadBMP(const char* filename) {
        Texture tex;
        FILE* f = std::fopen(filename, "rb");
        if (!f) {
            std::printf("ERROR: Cannot open texture file %s\n", filename);
            return tex;
        }

        // Чтение заголовка BMP (54 байта)
        unsigned char header[54];
        if (std::fread(header, 1, 54, f) != 54) {
            std::printf("ERROR: Not a valid BMP file %s\n", filename);
            std::fclose(f);
            return tex;
        }

        // Проверка сигнатуры "BM"
        if (header[0] != 'B' || header[1] != 'M') {
            std::printf("ERROR: Not a valid BMP file %s\n", filename);
            std::fclose(f);
            return tex;
        }

        // Извлечение параметров из заголовка
        tex.dataPos = *(int*)&(header[0x0A]);  // Смещение пиксельных данных
        tex.imageSize = *(int*)&(header[0x22]); // Размер изображения в байтах
        tex.width = *(int*)&(header[0x12]);     // Ширина
        tex.height = *(int*)&(header[0x16]);    // Высота

        // Уточнение параметров, если они нулевые
        if (tex.imageSize == 0) tex.imageSize = tex.width * tex.height * 3;
        if (tex.dataPos == 0) tex.dataPos = 54;

        // Проверка: только 24-битные изображения
        if (*(short*)&(header[0x1C]) != 24) {
            std::printf("ERROR: Only 24-bit BMP supported %s\n", filename);
            std::fclose(f);
            return tex;
        }

        // Выделение памяти и чтение пиксельных данных
        tex.data = new unsigned char[tex.imageSize];
        std::fseek(f, tex.dataPos, SEEK_SET);
        size_t read = std::fread(tex.data, 1, tex.imageSize, f);
        std::fclose(f);

        if (read != (size_t)tex.imageSize) {
            std::printf("ERROR: Could not read full texture data %s\n", filename);
            delete[] tex.data;
            tex.data = nullptr;
        } else {
            std::printf("Texture loaded: %s (%dx%d)\n", filename, tex.width, tex.height);
        }

        return tex;
    }
};

// ============================================================================
// МАТЕРИАЛ ОБЪЕКТА
// Определяет, как объект взаимодействует со светом
// ============================================================================

struct Material {
    Color color;           // Базовый цвет (если нет текстуры)
    float ka, kd, ks;      // Коэффициенты: ambient, diffuse, specular
    float reflect;         // Отражательная способность (0..1)
    float transparency;    // Прозрачность (0..1)
    int shininess;         // Блеск (для specular highlights)
    Texture* texture;      // Указатель на текстуру (может быть nullptr)

    // Конструктор по умолчанию
    Material() : color(1, 1, 1), ka(0.1f), kd(0.7f), ks(0.3f),
                 reflect(0.0f), transparency(0.0f), shininess(32), texture(nullptr) {}

    // Параметризованный конструктор
    Material(Color c, float a, float d, float s, float r, float t, int sh, Texture* tex = nullptr)
        : color(c), ka(a), kd(d), ks(s), reflect(r), transparency(t), shininess(sh), texture(tex) {}
};

// ============================================================================
// ЛУЧ И ИНФОРМАЦИЯ О ПЕРЕСЕЧЕНИИ
// ============================================================================

// Луч: точка начала и направление
struct Ray {
    Vec3 origin, direction;
    Ray(Vec3 o, Vec3 d) : origin(o), direction(d.normalize()) {}
    // Точка на луче на расстоянии t
    Vec3 at(float t) const { return origin + direction * t; }
};

// Результат пересечения луча с объектом
struct HitInfo {
    bool hit = false;      // Было ли пересечение?
    float t = 1e30f;       // Расстояние до точки пересечения
    Vec3 point, normal;    // Точка и нормаль в точке пересечения
    Material mat;          // Материал объекта
    float u = 0, v = 0;    // UV-координаты для текстурирования
    HitInfo() {}
};

// ============================================================================
// БАЗОВЫЙ КЛАСС ДЛЯ ГЕОМЕТРИЧЕСКИХ ФИГУР
// ============================================================================

struct Shape {
    bool enabled = true;   // Включен ли объект в рендер?
    virtual ~Shape() = default;
    // Чисто виртуальный метод пересечения
    virtual HitInfo intersect(const Ray& ray) const = 0;
    virtual bool isEnabled() const { return enabled; }
};

// ============================================================================
// СФЕРА
// ============================================================================

struct Sphere : public Shape {
    Vec3 center; float radius; Material mat;
    Sphere(Vec3 c, float r, Material m, bool en = true) : center(c), radius(r), mat(m) { enabled = en; }

    HitInfo intersect(const Ray& ray) const override {
        if (!isEnabled()) return HitInfo();

        // Вектор от центра сферы к началу луча
        Vec3 oc = ray.origin - center;

        // Коэффициенты квадратного уравнения |O + tD - C|² = R²
        float a = 1.0f;
        float b = 2.0f * oc.dot(ray.direction);
        float c = oc.dot(oc) - radius * radius;
        float disc = b * b - 4 * a * c;

        if (disc < 0) return HitInfo(); // Нет пересечения

        // Вычисление корней
        float sqrt_disc = std::sqrt(disc);
        float t = (-b - sqrt_disc) / (2 * a);

        // Если ближайшая точка за камерой — пробуем дальнюю
        if (t < 0.001f) {
            t = (-b + sqrt_disc) / (2 * a);
            if (t < 0.001f) return HitInfo();
        }

        HitInfo info;
        info.hit = true;
        info.t = t;
        info.point = ray.at(t);
        info.normal = (info.point - center).normalize();
        info.mat = mat;

        // Преобразование нормали в UV (сферическая проекция)
        float theta = std::acos(std::clamp(-info.normal.y, -1.0f, 1.0f)); // Полярный угол
        float phi = std::atan2(info.normal.z, info.normal.x);             // Азимут
        info.u = (phi + 3.1415926535f) / (2 * 3.1415926535f);            // [0,1]
        info.v = theta / 3.1415926535f;                                   // [0,1]

        return info;
    }
};

// ============================================================================
// ТРЕУГОЛЬНИК (АЛГОРИТМ МЁЛЛЕРА–ТРУМБОРА)
// ============================================================================

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

        Vec3 edge1 = v1 - v0;
        Vec3 edge2 = v2 - v0;
        Vec3 h = ray.direction.cross(edge2);
        float a = edge1.dot(h);

        if (std::abs(a) < 1e-6f) return HitInfo(); // Луч параллелен плоскости

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

        // Коррекция нормали: должна смотреть навстречу лучу
        if (info.normal.dot(ray.direction) > 0) info.normal = -info.normal;

        info.mat = mat;

        // Барицентрическая интерполяция UV
        float w = 1.0f - u - v;
        info.u = w * uv0.x + u * uv1.x + v * uv2.x;
        info.v = w * uv0.y + u * uv1.y + v * uv2.y;

        return info;
    }
};

// ============================================================================
// ТЕТРАЭДР (СОСТОИТ ИЗ 4 ТРЕУГОЛЬНИКОВ)
// ============================================================================

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
        // Определение 4 граней тетраэдра
        Triangle faces[4] = {
            Triangle(v[0], v[1], v[2], Vec2(0,0), Vec2(1,0), Vec2(0,1), mat, true),
            Triangle(v[0], v[2], v[3], Vec2(0,0), Vec2(0,1), Vec2(1,1), mat, true),
            Triangle(v[0], v[3], v[1], Vec2(0,0), Vec2(1,1), Vec2(1,0), mat, true),
            Triangle(v[1], v[3], v[2], Vec2(1,0), Vec2(1,1), Vec2(0,1), mat, true)
        };

        // Проверка пересечения с каждой гранью
        for (auto& tri : faces) {
            HitInfo h = tri.intersect(ray);
            if (h.hit && h.t < best.t) {
                // Коррекция нормали
                if (h.normal.dot(ray.direction) > 0) {
                    h.normal = -h.normal;
                }
                best = h;
            }
        }
        return best;
    }
};

// ============================================================================
// БЕСКОНЕЧНАЯ ПЛОСКОСТЬ
// ============================================================================

struct Plane : public Shape {
    Vec3 point, normal; Material mat;
    Plane(Vec3 p, Vec3 n, Material m, bool en = true) : point(p), normal(n.normalize()), mat(m) { enabled = en; }

    HitInfo intersect(const Ray& ray) const override {
        if (!isEnabled()) return HitInfo();

        float denom = normal.dot(ray.direction);
        if (std::abs(denom) < 1e-6f) return HitInfo(); // Луч параллелен плоскости

        float t = (point - ray.origin).dot(normal) / denom;
        if (t < 0.001f) return HitInfo();

        HitInfo info;
        info.hit = true;
        info.t = t;
        info.point = ray.at(t);
        // Нормаль направлена навстречу лучу
        info.normal = denom > 0 ? -normal : normal;
        info.mat = mat;
        // Простая проекция UV по координатам X и Z
        info.u = info.point.x * 2.0f;
        info.v = info.point.z * 2.0f;

        return info;
    }
};

// ============================================================================
// ГЛОБАЛЬНЫЕ ПЕРЕМЕННЫЕ
// ============================================================================

const int RENDER_WIDTH = 1920;   // Ширина рендеримого изображения
const int RENDER_HEIGHT = 1080;  // Высота рендеримого изображения
const int DISPLAY_WIDTH = 1920; // Размер окна OpenGL
const int DISPLAY_HEIGHT = 1080;
unsigned char* image = nullptr; // Буфер пикселей (RGB, 8 бит на канал)
GLuint textureID = 0;           // Идентификатор текстуры OpenGL

Vec3 camera(0, 1, 5);           // Позиция камеры
Vec3 lightPos(2, 5, 2);         // Позиция точечного источника света
bool useSSAA = false;           // Включено ли сглаживание (Super-Sampling Anti-Aliasing)
std::vector<std::unique_ptr<Shape>> scene; // Список объектов в сцене
const int MAX_DEPTH = 5;        // Максимальная глубина рекурсии (отражения/преломления)
const float FOV = 60.0f;        // Угол обзора в градусах

// Параметры рендеринга (предвычислены для ускорения)
struct RenderParams {
    float invWidth, invHeight, aspect, tanHalfFov;
    void update(int w, int h, float fov) {
        invWidth = 1.0f / w;
        invHeight = 1.0f / h;
        aspect = (float)w / h;
        tanHalfFov = std::tan(fov * 0.5f * 3.1415926535f / 180.0f);
    }
} g_renderParams;

// ============================================================================
// ОТЛАДОЧНЫЕ ФУНКЦИИ
// ============================================================================

void logDebugInfo(const std::string& message) {
    std::printf("[DEBUG] %s\n", message.c_str());
}

void logCameraInfo() {
    std::printf("[CAMERA] Pos: (%.2f, %.2f, %.2f)\n", camera.x, camera.y, camera.z);
}

void logLightInfo() {
    std::printf("[LIGHT] Pos: (%.2f, %.2f, %.2f)\n", lightPos.x, lightPos.y, lightPos.z);
}

void logSceneInfo() {
    std::printf("[SCENE] Objects: %zu\n", scene.size());
    for (size_t i = 0; i < scene.size(); ++i) {
        std::printf("  Object %zu: %s\n", i, scene[i]->enabled ? "ON" : "OFF");
    }
}

void logMaterialInfo(const Material& mat, const std::string& name) {
    std::printf("[Material %s] ka=%.2f, kd=%.2f, ks=%.2f, reflection=%.2f, transparency=%.2f, shininess=%d\n",
        name.c_str(), mat.ka, mat.kd, mat.ks, mat.reflect, mat.transparency, mat.shininess);
}

void logRenderInfo() {
    std::printf("[RENDER] SSAA: %s\n", useSSAA ? "ON" : "OFF");
}

// ============================================================================
// МОДЕЛЬ ОСВЕЩЕНИЯ COOK-TORRANCE
// Физически корректный specular
// ============================================================================

Color cookTorrance(const Vec3& N, const Vec3& V, const Vec3& L, const Color& baseColor, float ks, int shininess) {
    Vec3 H = (L + V).normalize(); // Half-vector

    float NdotL = std::max(0.0f, N.dot(L));
    float NdotV = std::max(0.0f, N.dot(V));
    float NdotH = std::max(0.0f, N.dot(H));
    float VdotH = std::max(0.0f, V.dot(H));

    if (NdotL <= 0 || NdotV <= 0) return Color(0, 0, 0);

    // Френель (приближение Шлика)
    float F0 = 0.04f; // Базовый коэффициент отражения для диэлектриков
    float F = F0 + (1.0f - F0) * std::pow(1.0f - VdotH, 5.0f);

    // Распределение микрограней (GGX)
    float alpha = std::max(0.001f, 1.0f - (float)shininess / 100.0f);
    float alpha2 = alpha * alpha;
    float denomD = NdotH * NdotH * (alpha2 - 1.0f) + 1.0f;
    float D = alpha2 / (3.1415926535f * denomD * denomD);

    // Геометрическое затенение (Smith model)
    float k = (alpha + 1.0f) * (alpha + 1.0f) / 8.0f;
    float G1V = NdotV / (NdotV * (1.0f - k) + k);
    float G1L = NdotL / (NdotL * (1.0f - k) + k);
    float G = G1V * G1L;

    // Итоговый specular
    float specular = (F * D * G) / std::max(0.001f, 4.0f * NdotV * NdotL);

    return Color(1, 1, 1) * ks * specular * NdotL;
}

// ============================================================================
// РЕКУРСИВНАЯ ТРАССИРОВКА ЛУЧЕЙ
// ============================================================================

Color traceRay(const Ray& ray, int depth = 0) {
    // Ограничение глубины рекурсии
    if (depth >= MAX_DEPTH) return Color(0.05f, 0.05f, 0.1f); // Тёмно-синий фон

    // Поиск ближайшего пересечения
    HitInfo closest;
    for (const auto& s : scene) {
        if (!s->isEnabled()) continue;
        HitInfo h = s->intersect(ray);
        if (h.hit && h.t < closest.t && h.t > 0.001f) closest = h;
    }

    if (!closest.hit) return Color(0.05f, 0.05f, 0.05f); // Серый фон

    // Основной цвет: либо базовый, либо из текстуры
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

    float lightDist = std::sqrt(lightDist2);
    Vec3 L = L_vec / lightDist; // Нормализованное направление к свету
    Vec3 V = (camera - closest.point).normalize(); // Направление к камере
    Vec3 N = closest.normal; // Нормаль в точке

    // Проверка тени
    bool inShadow = false;
    Ray shadowRay(closest.point + N * 0.01f, L); // Небольшое смещение, чтобы избежать self-intersection
    for (const auto& s : scene) {
        if (!s->isEnabled()) continue;
        HitInfo h = s->intersect(shadowRay);
        if (h.hit && h.t > 0.001f && h.t < lightDist - 0.01f) {
            inShadow = true;
            break;
        }
    }

    // Ambient
    Color ambient = baseColor * closest.mat.ka;
    Color result = ambient;

    if (!inShadow) {
        // Diffuse (Lambert)
        float diff = std::max(0.0f, N.dot(L));
        Color diffuse = baseColor * closest.mat.kd * diff;

        // Specular (Cook-Torrance)
        Color specular = cookTorrance(N, V, L, baseColor, closest.mat.ks, closest.mat.shininess);

        result = result + diffuse + specular;
    }

    // Рекурсивное отражение
    if (closest.mat.reflect > 0.01f && depth < MAX_DEPTH - 1) {
        Vec3 reflDir = ray.direction.reflect(N);
        Ray reflRay(closest.point + N * 0.01f, reflDir);
        Color reflCol = traceRay(reflRay, depth + 1);
        result = result * (1.0f - closest.mat.reflect) + reflCol * closest.mat.reflect;
    }

    // Рекурсивное преломление
    if (closest.mat.transparency > 0.01f && depth < MAX_DEPTH - 1) {
        // Вычисление коэффициента преломления
        float eta = N.dot(ray.direction) < 0 ? (1.0f / 1.5f) : 1.5f; // 1.5 — IOR стекла
        Vec3 refrDir = ray.direction.refract(N, eta);

        if (refrDir.length2() > 0) {
            Ray refrRay(closest.point - N * 0.01f, refrDir);
            Color refrCol = traceRay(refrRay, depth + 1);
            result = result * (1.0f - closest.mat.transparency) + refrCol * closest.mat.transparency;
        }
    }

    return result;
}

// ============================================================================
// РЕНДЕРИНГ ИЗОБРАЖЕНИЯ
// ============================================================================

void renderImage() {
    logDebugInfo("Start rendering");
    auto startTime = std::chrono::high_resolution_clock::now();

    g_renderParams.update(RENDER_WIDTH, RENDER_HEIGHT, FOV);

    const int NUM_THREADS = std::thread::hardware_concurrency();
    const int THREADS = (NUM_THREADS > 0) ? NUM_THREADS : 4;

    std::vector<std::thread> threads;
    threads.reserve(THREADS);

    int rowsPerThread = RENDER_HEIGHT / THREADS;
    const bool ssaa = useSSAA;

    // Создание потоков для рендеринга полос изображения
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

                            // Преобразование пикселя в направление луча
                            float px = (2.0f * (x + jitterX) * invWidth - 1.0f) * aspect * tanHalfFov;
                            float py = (1.0f - 2.0f * (y + jitterY) * invHeight) * tanHalfFov;
                            Vec3 rayDir(px, py, -1);
                            Ray ray(camera, rayDir);

                            finalColor = finalColor + traceRay(ray);
                        }
                    }

                    finalColor = finalColor * (1.0f / samples);

                    // Запись в буфер (с переворотом по Y для OpenGL)
                    int idx = ((RENDER_HEIGHT - 1 - y) * RENDER_WIDTH + x) * 3;
                    image[idx] = finalColor.getR();
                    image[idx + 1] = finalColor.getG();
                    image[idx + 2] = finalColor.getB();
                }
            }
        });
    }

    // Ожидание завершения всех потоков
    for (auto& t : threads) {
        if (t.joinable()) t.join();
    }

    // Обновление текстуры OpenGL
    glBindTexture(GL_TEXTURE_2D, textureID);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, RENDER_WIDTH, RENDER_HEIGHT, GL_RGB, GL_UNSIGNED_BYTE, image);

    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    logDebugInfo("Render completed in " + std::to_string(duration.count()) + " ms");
}

// ============================================================================
// ЗАГРУЗКА СЦЕНЫ ИЗ ФАЙЛА
// Формат: тип x y z ... параметры материала текстура
// ============================================================================

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
        if (line.empty() || line[0] == '#') continue; // Пропуск комментариев

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
            logDebugInfo("Sphere loaded at (" + std::to_string(cx) + ", " + std::to_string(cy) + ", " + std::to_string(cz) + ")");
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
            logDebugInfo("Plane loaded at (" + std::to_string(px) + ", " + std::to_string(py) + ", " + std::to_string(pz) + ")");
        }
    }
    logDebugInfo("Scene loaded: " + std::to_string(objectsLoaded) + " objects");
}

// ============================================================================
// ОБНОВЛЕНИЕ ПАРАМЕТРОВ МАТЕРИАЛА
// ============================================================================

void updateMaterial(Material& mat, int action) {
    const float delta = 0.05f;
    switch (action) {
    case 100: mat.ka = std::min(mat.ka + delta, 1.0f); break; // Ambient +
    case 101: mat.ka = std::max(mat.ka - delta, 0.0f); break; // Ambient -
    case 102: mat.kd = std::min(mat.kd + delta, 1.0f); break; // Diffuse +
    case 103: mat.kd = std::max(mat.kd - delta, 0.0f); break; // Diffuse -
    case 104: mat.ks = std::min(mat.ks + delta, 1.0f); break; // Specular +
    case 105: mat.ks = std::max(mat.ks - delta, 0.0f); break; // Specular -
    case 106: mat.reflect = std::min(mat.reflect + delta, 1.0f); break; // Отражение +
    case 107: mat.reflect = std::max(mat.reflect - delta, 0.0f); break; // Отражение -
    case 108: mat.transparency = std::min(mat.transparency + delta, 1.0f); break; // Прозрачность +
    case 109: mat.transparency = std::max(mat.transparency - delta, 0.0f); break; // Прозрачность -
    }
}

// ============================================================================
// ОБРАБОТЧИК КЛАВИАТУРЫ
// ============================================================================

void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    if (action != GLFW_PRESS) return;

    if (key == GLFW_KEY_ESCAPE) {
        glfwSetWindowShouldClose(window, GLFW_TRUE);
        logDebugInfo("ESC pressed — exiting");
        return;
    }

    bool changed = false;

    auto updateAndRender = [&]() {
        if (changed) renderImage();
    };

    // Управление камерой
    if (key == GLFW_KEY_W) { camera.z -= 0.3f; logDebugInfo("W: camera forward"); changed = true; }
    else if (key == GLFW_KEY_S) { camera.z += 0.3f; logDebugInfo("S: camera back"); changed = true; }
    else if (key == GLFW_KEY_A) { camera.x -= 0.3f; logDebugInfo("A: camera left"); changed = true; }
    else if (key == GLFW_KEY_D) { camera.x += 0.3f; logDebugInfo("D: camera right"); changed = true; }
    else if (key == GLFW_KEY_Q) { camera.y += 0.3f; logDebugInfo("Q: camera up"); changed = true; }
    else if (key == GLFW_KEY_E) { camera.y -= 0.3f; logDebugInfo("E: camera down"); changed = true; }
    else if (key == GLFW_KEY_C) { camera = Vec3(0, 1, 5); logDebugInfo("C: reset camera"); changed = true; }

    // Управление светом
    else if (key == GLFW_KEY_I) { lightPos.y += 0.3f; logDebugInfo("I: light up"); changed = true; }
    else if (key == GLFW_KEY_K) { lightPos.y -= 0.3f; logDebugInfo("K: light down"); changed = true; }
    else if (key == GLFW_KEY_J) { lightPos.x -= 0.3f; logDebugInfo("J: light left"); changed = true; }
    else if (key == GLFW_KEY_L) { lightPos.x += 0.3f; logDebugInfo("L: light right"); changed = true; }
    else if (key == GLFW_KEY_O) { lightPos = Vec3(2, 5, 2); logDebugInfo("O: reset light"); changed = true; }

    // Включение/выключение объектов
    else if (key == GLFW_KEY_1 && scene.size() >= 1) {
        scene[0]->enabled = !scene[0]->enabled;
        logDebugInfo("1: toggle sphere");
        changed = true;
    }
    else if (key == GLFW_KEY_2 && scene.size() >= 2) {
        scene[1]->enabled = !scene[1]->enabled;
        logDebugInfo("2: toggle tetrahedron");
        changed = true;
    }
    else if (key == GLFW_KEY_3 && scene.size() >= 3) {
        scene[2]->enabled = !scene[2]->enabled;
        logDebugInfo("3: toggle plane");
        changed = true;
    }

    // Управление рендером
    else if (key == GLFW_KEY_F) {
        useSSAA = !useSSAA;
        logDebugInfo("F: toggle SSAA");
        changed = true;
    }
    else if (key == GLFW_KEY_R) {
        logDebugInfo("R: re-render");
        changed = true;
    }
    else if (key == GLFW_KEY_F5) {
        loadSceneFromFile("scene.txt");
        changed = true;
    }

    // Пример управления материалом сферы
    else if (key == GLFW_KEY_U && scene.size() >= 1) {
        if (Sphere* s = dynamic_cast<Sphere*>(scene[0].get())) {
            updateMaterial(s->mat, 100); // ka+
            logDebugInfo("U: sphere ka+");
            changed = true;
        }
    }

    updateAndRender();

    // Вывод информации о материале сферы при изменении
    if (changed && scene.size() >= 1) {
        if (Sphere* s = dynamic_cast<Sphere*>(scene[0].get())) {
            logMaterialInfo(s->mat, "SPHERE");
        }
    }
}

// ============================================================================
// ТОЧКА ВХОДА
// ============================================================================

int main() {
    std::printf("=== Raytracer started ===\n");
    std::printf("Controls:\n");
    std::printf("  Camera: W/S/A/D/Q/E — move, C — reset\n");
    std::printf("  Light: I/K/J/L — move, O — reset\n");
    std::printf("  Objects: 1/2/3 — toggle\n");
    std::printf("  Render: F — SSAA, R — re-render, F5 — reload scene\n");
    std::printf("  Exit: ESC\n\n");

    // Инициализация GLFW
    if (!glfwInit()) {
        std::fprintf(stderr, "Failed to initialize GLFW\n");
        return -1;
    }

    // Настройка контекста OpenGL (Legacy, совместим с glBegin/glEnd)
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);

    // Создание окна
    GLFWwindow* window = glfwCreateWindow(DISPLAY_WIDTH, DISPLAY_HEIGHT, "LR3 - Cross-platform Raytracer", nullptr, nullptr);
    if (!window) {
        std::fprintf(stderr, "Failed to create GLFW window\n");
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);
    glfwSetKeyCallback(window, keyCallback);
    glfwSetMouseButtonCallback(window, [](GLFWwindow*, int, int, int){}); // Игнор кликов

    // Загрузка сцены
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

    // Первый рендер
    renderImage();

    // Главный цикл
    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(0, RENDER_WIDTH, 0, RENDER_HEIGHT, -1, 1);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        // Отображение рендеримого изображения как текстуры
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

    std::printf("=== Raytracer finished ===\n");

    // Освобождение ресурсов
    delete[] image;
    glDeleteTextures(1, &textureID);
    glfwTerminate();
    return 0;
}