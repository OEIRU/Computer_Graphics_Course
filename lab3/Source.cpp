// compile
// Linux: g++ -O3 -std=c++17 -o main Source.cpp -lglut -lGL -lGLU -lpthread -lm
// Windows: cl /EHsc /O2 /std:c++17 Source.cpp /link glut32.lib opengl32.lib glu32.lib

#define _CRT_SECURE_NO_WARNINGS
#define NOMINMAX

#ifdef _WIN32
#include <windows.h>
#endif

#include ".\GL\glut.h"

#include <thread>
#include <vector>
#include <memory>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <cstdlib>
#include <mutex>
#include <atomic>
#include <condition_variable>
// ============================================================================
// ВСПОМОГАТЕЛЬНЫЕ СТРУКТУРЫ
// ============================================================================

// Двумерный вектор (для UV-координат текстур)
struct Vec2 {
    float x, y;
    Vec2(float x = 0, float y = 0) : x(x), y(y) {}
};

// Трёхмерный вектор — основа геометрических вычислений
struct Vec3 {
    float x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(float x, float y, float z) : x(x), y(y), z(z) {}

    // Базовые арифметические операции над векторами
    Vec3 operator+(const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
    Vec3 operator-(const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
    Vec3 operator-() const { return Vec3(-x, -y, -z); }
    Vec3 operator*(float t) const { return Vec3(x * t, y * t, z * t); }
    Vec3 operator/(float t) const { return t != 0 ? Vec3(x / t, y / t, z / t) : Vec3(0, 0, 0); }

    // Скалярное произведение — для углов, проекций, освещения
    float dot(const Vec3& v) const { return x * v.x + y * v.y + z * v.z; }

    // Векторное произведение — для нормалей и ориентации
    Vec3 cross(const Vec3& v) const {
        return Vec3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }

    // Длина и нормализация — критичны для корректного направления лучей и нормалей
    float length2() const { return dot(*this); }
    float length() const { return std::sqrt(length2()); }

    Vec3 normalize() const {
        float len = length();         // Избегаем деления на ноль; fallback — вдоль оси Z
        return len > 1e-6f ? (*this) * (1.0f / len) : Vec3(0, 0, 1);
    }

    // Отражение луча относительно нормали — для зеркальных поверхностей
    Vec3 reflect(const Vec3& n) const {
        return *this - n * (2.0f * this->dot(n));
    }

    // Преломление луча (закон Снеллиуса) — для стеклянных/прозрачных объектов
    Vec3 refract(const Vec3& n, float eta) const {
        float cosI = -n.dot(*this);
        float sinT2 = eta * eta * (1.0f - cosI * cosI);
        if (sinT2 >= 1.0f) return Vec3(0, 0, 0);
        float cosT = std::sqrt(1.0f - sinT2);
        return (*this) * eta + n * (eta * cosI - cosT);
    }
};

// Цвет в линейном пространстве (RGB, значения [0..1])
struct Color {
    float r, g, b;
    Color() : r(0), g(0), b(0) {}
    Color(float r, float g, float b) {
        this->r = std::clamp(r, 0.0f, 1.0f);
        this->g = std::clamp(g, 0.0f, 1.0f);
        this->b = std::clamp(b, 0.0f, 1.0f);
    }

    // Операторы для смешивания цветов и затенения; поэлементное умножение (для текстур)
    Color operator*(float t) const { return Color(r * t, g * t, b * t); }
    Color operator+(const Color& c) const { return Color(r + c.r, g + c.g, b + c.b); }
    Color operator*(const Color& c) const { return Color(r * c.r, g * c.g, b * c.b); }

    Color blend(const Color& c, float t) const {
        t = std::clamp(t, 0.0f, 1.0f);
        return Color(r * (1 - t) + c.r * t, g * (1 - t) + c.g * t, b * (1 - t) + c.b * t);
    }

    // Преобразование в 8-битный формат для вывода на экран
    unsigned char getR() const { return (unsigned char)(r * 255); }
    unsigned char getG() const { return (unsigned char)(g * 255); }
    unsigned char getB() const { return (unsigned char)(b * 255); }
};

// Простая текстура, загружаемая из BMP (24-bit, без сжатия)
struct Texture {
    unsigned char* data = nullptr; // сырые данные в формате BGR (как в BMP)
    int width = 0, height = 0;
    int imageSize = 0;             // размер пиксельных данных в байтах
	int dataPos = 0;			  // смещение начала пиксельных данных в файле

	// Деструктор для очистки памяти
    ~Texture() { delete[] data; }

	// Выборка цвета по UV-координатам (повторяющаяся текстура)
    Color sample(float u, float v) const {
        if (!data || width <= 0 || height <= 0) return Color(1, 0, 1);

		// Повторение текстуры
        u = u - std::floor(u);
        v = v - std::floor(v);

        int x = (int)(u * width) % width;
        int y = (int)(v * height) % height;
        if (x < 0) x += width;
        if (y < 0) y += height;

		// Хранение строки в BMP идет снизу вверх
        int idx = (y * width + x) * 3;
        return Color(
            data[idx + 2] / 255.0f,
            data[idx + 1] / 255.0f,
            data[idx] / 255.0f
        );
    }

	// Загрузка текстуры из BMP-файла
    static Texture loadBMP(const char* filename) {
        Texture tex;
        FILE* f = std::fopen(filename, "rb");
        if (!f) {
            std::printf("ERROR: Cannot open texture file %s\n", filename);
            return tex;
        }

		// Чтение заголовка BMP (54 байта - стандарт для 24-bit)
        unsigned char header[54];
        if (std::fread(header, 1, 54, f) != 54) {
            std::printf("ERROR: Not a valid BMP file %s\n", filename);
            std::fclose(f);
            return tex;
        }

		// Проверка сигнатуры 'BM'
        if (header[0] != 'B' || header[1] != 'M') {
            std::printf("ERROR: Not a valid BMP file %s\n", filename);
            std::fclose(f);
            return tex;
        }

		// Извлечение метаданных из заголовка
        tex.dataPos = *(int*)&(header[0x0A]);
        tex.imageSize = *(int*)&(header[0x22]);
        tex.width = *(int*)&(header[0x12]);
        tex.height = *(int*)&(header[0x16]);

        // Расчет размера, если он не задан явно 
        if (tex.imageSize == 0) tex.imageSize = tex.width * tex.height * 3;
        if (tex.dataPos == 0) tex.dataPos = 54;

		// Проверка на 24-битный формат
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
        }
        else {
            std::printf("Texture loaded: %s (%dx%d)\n", filename, tex.width, tex.height);
        }

        return tex;
    }
};
// Материал поверхности с параметрами освещения и текстурой 
struct Material {
	Color color;        // базовый цвет
	float ka, kd, ks;   // коэффициенты амбиентного, диффузного и зеркального отражения
	float reflect;      // коэффициент отражения
	float transparency; // коэффициент прозрачности
	int shininess;      // жёсткость блика
	Texture* texture;   // указатель на текстуру (если есть)

	// Конструктор по умолчанию
    Material() : color(1, 1, 1), ka(0.1f), kd(0.7f), ks(0.3f),
        reflect(0.0f), transparency(0.0f), shininess(32), texture(nullptr) {
    }
	// Полный конструктор
    Material(Color c, float a, float d, float s, float r, float t, int sh, Texture* tex = nullptr)
        : color(c), ka(a), kd(d), ks(s), reflect(r), transparency(t), shininess(sh), texture(tex) {
    }
};

// Луч с началом и направлением
struct Ray {
    Vec3 origin, direction;
    Ray(Vec3 o, Vec3 d) : origin(o), direction(d.normalize()) {}
    Vec3 at(float t) const { return origin + direction * t; }
};

// Информация о пересечении луча с объектом
struct HitInfo {
	bool hit = false;       // было ли пересечение
	float t = 1e30f;        // расстояние до пересечения
	Vec3 point, normal;     // точка пересечения и нормаль в этой точке
	Material mat;           // материал объекта в точке
	float u = 0, v = 0;     // UV-координаты для текстурирования
    HitInfo() {}
};

// Базовый класс для всех геометрических объектов сцены
struct Shape {
    bool enabled = true;
    virtual ~Shape() = default;
    virtual HitInfo intersect(const Ray& ray) const = 0;
    virtual bool isEnabled() const { return enabled; }
};

// Сфера — простой геометрический объект
struct Sphere : public Shape {
    Vec3 center; float radius; Material mat;
    Sphere(Vec3 c, float r, Material m, bool en = true) : center(c), radius(r), mat(m) { enabled = en; }

    HitInfo intersect(const Ray& ray) const override {
        if (!isEnabled()) return HitInfo();

		// Решение квадратного уравнения для пересечения луча и сферы
        Vec3 oc = ray.origin - center;
        float a = 1.0f;
        float b = 2.0f * oc.dot(ray.direction);
        float c = oc.dot(oc) - radius * radius;
        float disc = b * b - 4 * a * c;

		if (disc < 0) return HitInfo(); // нет пересечения

        float sqrt_disc = std::sqrt(disc);
        float t = (-b - sqrt_disc) / (2 * a);

		// Проверка ближайшей точки, 
        // если она слишком близко - проверяем дальше 
        if (t < 0.001f) {
            t = (-b + sqrt_disc) / (2 * a);
            if (t < 0.001f) return HitInfo(); // за камерой или внутри
        }

        HitInfo info;
        info.hit = true;
        info.t = t;
        info.point = ray.at(t);
        info.normal = (info.point - center).normalize();
        info.mat = mat;

		// Вычисление UV-координат для текстурирования сферы
        float theta = std::acos(std::clamp(-info.normal.y, -1.0f, 1.0f));
        float phi = std::atan2(info.normal.z, info.normal.x);
        info.u = (phi + 3.1415926535f) / (2 * 3.1415926535f);
        info.v = theta / 3.1415926535f;

        return info;
    }
};

// Треугольник с поддержкой UV-координат барицентрического типа
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

		// Алгоритм Мёллера–Трамбора для пересечения луча с треугольником
        Vec3 edge1 = v1 - v0;
        Vec3 edge2 = v2 - v0;
        Vec3 h = ray.direction.cross(edge2);
        float a = edge1.dot(h);

        if (std::abs(a) < 1e-6f) return HitInfo(); // Луч параллелен плоскости треугольника 

        float f = 1.0f / a;
        Vec3 s = ray.origin - v0;
        float u = f * s.dot(h);
        if (u < 0 || u > 1) return HitInfo();

        Vec3 q = s.cross(edge1);
        float v = f * ray.direction.dot(q);
        if (v < 0 || u + v > 1) return HitInfo();

        float t = f * edge2.dot(q);
		if (t < 0.001f) return HitInfo(); // пересечение слишком близко или за камерой

        HitInfo info;
        info.hit = true;
        info.t = t;
        info.point = ray.at(t);
        info.normal = normal;

		// направление нормали и направление взгляда 
        if (info.normal.dot(ray.direction) > 0) info.normal = -info.normal;

        info.mat = mat;

		// Барицентрические координаты для UV-текстурирования
        float w = 1.0f - u - v;
        info.u = w * uv0.x + u * uv1.x + v * uv2.x;
        info.v = w * uv0.y + u * uv1.y + v * uv2.y;

        return info;
    }
};

// Тетраэдр, состоящий из четырёх треугольных граней
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
		// Четыре грани тетраэдра
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
                if (h.normal.dot(ray.direction) > 0) {
                    h.normal = -h.normal;
                }
                best = h;
            }
        }
        return best;
    }
};

// Бесконечная плоскость, определяемая точкой и нормалью
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
		info.normal = denom > 0 ? -normal : normal; // Нормаль направлена к лучу
        info.mat = mat;
		// Простейшее текстурирование по XZ-плоскости
        info.u = info.point.x * 2.0f;
        info.v = info.point.z * 2.0f;

        return info;
    }
};

// ============================================================================
// ГЛОБАЛЬНЫЕ ПЕРЕМЕННЫЕ
// ============================================================================

int RENDER_WIDTH = 1280;
int RENDER_HEIGHT = 720;
unsigned char* renderBuffer = nullptr;  // основной буфер рендеринга
unsigned char* displayBuffer = nullptr; // буфер для отображения
GLuint textureID = 0;                   // OpenGL текстура для вывода
int windowWidth = RENDER_WIDTH;         // размеры окна
int windowHeight = RENDER_HEIGHT; 

Vec3 camera(0, 1, 5);    // позиция камеры
Vec3 lightPos(2, 5, 2);  // позиция точечного источника света
bool useSSAA = false;    // флаг суперсэмплинга
std::vector<std::unique_ptr<Shape>> scene; // объекты сцены

const int MAX_DEPTH = 5; // максимальная глубина рекурсии лучей
const float FOV = 60.0f; // угол обзора камеры

// Параметры рендеринга, зависящие от разрешения и FOV
struct RenderParams {
    float invWidth, invHeight, aspect, tanHalfFov;
    void update(int w, int h, float fov) {
        invWidth = 1.0f / w;
        invHeight = 1.0f / h;
        aspect = (float)w / h;
        tanHalfFov = std::tan(fov * 0.5f * 3.1415926535f / 180.0f);
    }
} g_renderParams;

std::mutex renderMutex;
std::atomic<bool> renderingInProgress(false);
std::atomic<bool> renderComplete(false);
std::atomic<int> currentRenderRow(0); // текущая строка рендеринга
std::atomic<bool> stopRendering(false); // флаг остановки рендеринга
std::thread renderWorker; // поток рендеринга

void stopRender();
void startRender(int w, int h);


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
    std::printf("[RENDER] SSAA: %s, Resolution: %dx%d\n", useSSAA ? "ON" : "OFF", RENDER_WIDTH, RENDER_HEIGHT);
}

// ============================================================================
// МОДЕЛЬ ОСВЕЩЕНИЯ COOK-TORRANCE
// ============================================================================

Color cookTorrance(const Vec3& N, const Vec3& V, const Vec3& L, const Color& baseColor, float ks, int shininess) {
    Vec3 H = (L + V).normalize();

    float NdotL = std::max(0.0f, N.dot(L));
    float NdotV = std::max(0.0f, N.dot(V));
    float NdotH = std::max(0.0f, N.dot(H));
    float VdotH = std::max(0.0f, V.dot(H));

    if (NdotL <= 0 || NdotV <= 0) return Color(0, 0, 0);
    
    // Френель: отражение зависит от угла
    float F0 = 0.04f;
    float F = F0 + (1.0f - F0) * std::pow(1.0f - VdotH, 5.0f);

    // Распределение микрограней (GGX/Trowbridge-Reitz)
    float alpha = std::max(0.001f, 1.0f - (float)shininess / 100.0f);
    float alpha2 = alpha * alpha;
    float denomD = NdotH * NdotH * (alpha2 - 1.0f) + 1.0f;
    float D = alpha2 / (3.1415926535f * denomD * denomD);

    // Геометрическое затенение (Smith)
    float k = (alpha + 1.0f) * (alpha + 1.0f) / 8.0f;
    float G1V = NdotV / (NdotV * (1.0f - k) + k);
    float G1L = NdotL / (NdotL * (1.0f - k) + k);
    float G = G1V * G1L;

    // Итоговый вклад
    float specular = (F * D * G) / std::max(0.001f, 4.0f * NdotV * NdotL);

    return Color(1, 1, 1) * ks * specular * NdotL;
}

// ============================================================================
// РЕКУРСИВНАЯ ТРАССИРОВКА ЛУЧЕЙ
// ============================================================================

Color traceRay(const Ray& ray, int depth = 0) {
    if (depth >= MAX_DEPTH) return Color(0.05f, 0.05f, 0.1f); // Затухающий фон

	// Поиск ближайшего пересечения с объектами сцены
    HitInfo closest;
    for (const auto& s : scene) {
        if (!s->isEnabled()) continue;
        HitInfo h = s->intersect(ray);
        if (h.hit && h.t < closest.t && h.t > 0.001f) closest = h;
    }

    if (!closest.hit) return Color(0.05f, 0.05f, 0.05f);

	// Базовый цвет из текстуры или материала
    Color baseColor = closest.mat.color;
    if (closest.mat.texture) {
        baseColor = closest.mat.texture->sample(closest.u, closest.v);
    }

	// Направление к источнику света
    Vec3 L_vec = lightPos - closest.point;
    float lightDist2 = L_vec.length2();
    if (lightDist2 < 1e-6f) {
        return baseColor * closest.mat.ka;
    }

    float lightDist = std::sqrt(lightDist2);
    Vec3 L = L_vec / lightDist;
    Vec3 V = (camera - closest.point).normalize();
    Vec3 N = closest.normal;

	// Проверка на тени - запускаем луч к источнику и проверяем пересечение
    bool inShadow = false;
    Ray shadowRay(closest.point + N * 0.01f, L); // смещение!!! 
    for (const auto& s : scene) {
        if (!s->isEnabled()) continue;
        HitInfo h = s->intersect(shadowRay);
        if (h.hit && h.t > 0.001f && h.t < lightDist - 0.01f) {
            inShadow = true;
            break;
        }
    }

	// Вычисление освещения
    Color ambient = baseColor * closest.mat.ka;
    Color result = ambient;

    if (!inShadow) {
		// Диффузное и зеркальное освещение
        float diff = std::max(0.0f, N.dot(L)); 
        Color diffuse = baseColor * closest.mat.kd * diff;
        Color specular = cookTorrance(N, V, L, baseColor, closest.mat.ks, closest.mat.shininess);
        result = result + diffuse + specular;
    }

	// Рекурсивные отражения 
    if (closest.mat.reflect > 0.01f && depth < MAX_DEPTH - 1) {
        Vec3 reflDir = ray.direction.reflect(N);
        Ray reflRay(closest.point + N * 0.01f, reflDir);
        Color reflCol = traceRay(reflRay, depth + 1);
        result = result * (1.0f - closest.mat.reflect) + reflCol * closest.mat.reflect;
    }

	// преломление и прозрачность
    if (closest.mat.transparency > 0.01f && depth < MAX_DEPTH - 1) {
        float eta = N.dot(ray.direction) < 0 ? (1.0f / 1.5f) : 1.5f;
        Vec3 refrDir = ray.direction.refract(N, eta);

        if (refrDir.length2() > 0) { // преломление возможно
            Ray refrRay(closest.point - N * 0.01f, refrDir);
            Color refrCol = traceRay(refrRay, depth + 1);
            result = result * (1.0f - closest.mat.transparency) + refrCol * closest.mat.transparency;
        }
    }

    return result;
}

// ============================================================================
// ПОТОК РЕНДЕРИНГА
// ============================================================================

void renderThreadFunction(int renderW, int renderH) {
    renderingInProgress = true;
    renderComplete = false;
    stopRendering = false;
    currentRenderRow = 0;

    const int localWidth = renderW;
    const int localHeight = renderH;
    g_renderParams.update(localWidth, localHeight, FOV);

	// Количество потоков для рендеринга
    const int NUM_THREADS = std::thread::hardware_concurrency();
    const int THREADS = (NUM_THREADS > 0) ? NUM_THREADS : 4;
    const bool ssaa = useSSAA;

	// Параметры камеры и рендеринга
    float invWidth = g_renderParams.invWidth;
    float invHeight = g_renderParams.invHeight;
    float aspect = g_renderParams.aspect;
    float tanHalfFov = g_renderParams.tanHalfFov;

	// Временный буфер для рендеринга
    unsigned char* tempBuffer = new unsigned char[localWidth * localHeight * 3]();

    std::vector<std::thread> threads;
    threads.reserve(THREADS);

	// Запуск потоков рендеринга
    for (int t = 0; t < THREADS; ++t) {
        threads.emplace_back([&, t]() {
            while (true) {
                if (stopRendering) break;

                int row;
                {
					// Получение следующей строки для рендеринга
                    std::lock_guard<std::mutex> lock(renderMutex);
                    if (currentRenderRow >= localHeight) break;
                    row = currentRenderRow++;
                }

				// Рендеринг строки по пикселям
                for (int x = 0; x < localWidth; ++x) {
                    Color finalColor(0, 0, 0);
                    int samples = ssaa ? 4 : 1;

					// Суперсэмплинг 2x2 
                    for (int sy = 0; sy < (ssaa ? 2 : 1); ++sy) {
                        for (int sx = 0; sx < (ssaa ? 2 : 1); ++sx) {
                            float jitterX = ssaa ? (sx + 0.5f) / 2.0f : 0.5f;
                            float jitterY = ssaa ? (sy + 0.5f) / 2.0f : 0.5f;

							// Преобразуем экранные координаты в направление луча
                            float px = (2.0f * (x + jitterX) * invWidth - 1.0f) * aspect * tanHalfFov;
                            float py = (1.0f - 2.0f * (row + jitterY) * invHeight) * tanHalfFov;
                            Vec3 rayDir(px, py, -1);
                            Ray ray(camera, rayDir);

                            finalColor = finalColor + traceRay(ray);
                        }
                    }

                    finalColor = finalColor * (1.0f / samples);

					// Запись цвета в буфер (с учётом инверсии по Y)
                    int idx = ((localHeight - 1 - row) * localWidth + x) * 3;
                    if (idx >= 0 && idx + 2 < localWidth * localHeight * 3) {
                        tempBuffer[idx] = finalColor.getR();
                        tempBuffer[idx + 1] = finalColor.getG();
                        tempBuffer[idx + 2] = finalColor.getB();
                    }
                }
            }
            });
    }

	// Ожидание завершения всех потоков
    for (auto& t : threads) {
        if (t.joinable()) t.join();
    }

    {
        // Защищенная передача результата в основной поток
        std::lock_guard<std::mutex> lock(renderMutex);
        if (RENDER_WIDTH == localWidth && RENDER_HEIGHT == localHeight) {
            if (displayBuffer) {
                std::memcpy(displayBuffer, tempBuffer, localWidth * localHeight * 3);
                renderComplete = true;
            }
        }
    }

    delete[] tempBuffer;
    renderingInProgress = false;
}

// Остановка рендеринга 
void stopRender() {
    if (renderWorker.joinable()) {
        stopRendering = true; // request stop
        // join will wait until renderThreadFunction exits
        renderWorker.join();
    }
}

// Запуск рендеринга в отдельном потоке
void startRender(int w, int h) {
    // ensure any previous worker is stopped
    if (renderWorker.joinable()) stopRender();
    stopRendering = false;
    renderWorker = std::thread(renderThreadFunction, w, h);
}

// ============================================================================
// ЗАГРУЗКА СЦЕНЫ
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
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string type;
        iss >> type;

        if (type == "sphere") {
			// Параметры сферы
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
			// Параметры тетраэдра
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
			// Параметры плоскости
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
// ОБНОВЛЕНИЕ МАТЕРИАЛА
// ============================================================================

void updateMaterial(Material& mat, int action) {
    const float delta = 0.05f;
    switch (action) {
    case 100: mat.ka = std::min(mat.ka + delta, 1.0f); break;
    case 101: mat.ka = std::max(mat.ka - delta, 0.0f); break;
    case 102: mat.kd = std::min(mat.kd + delta, 1.0f); break;
    case 103: mat.kd = std::max(mat.kd - delta, 0.0f); break;
    case 104: mat.ks = std::min(mat.ks + delta, 1.0f); break;
    case 105: mat.ks = std::max(mat.ks - delta, 0.0f); break;
    case 106: mat.reflect = std::min(mat.reflect + delta, 1.0f); break;
    case 107: mat.reflect = std::max(mat.reflect - delta, 0.0f); break;
    case 108: mat.transparency = std::min(mat.transparency + delta, 1.0f); break;
    case 109: mat.transparency = std::max(mat.transparency - delta, 0.0f); break;
    }
}

// ============================================================================
// GLUT CALLBACKS
// ============================================================================

// Прототип
void menuCallback(int option);

// Обработка изменения размера окна
void reshape(int width, int height) {
    const int MAX_WIDTH = 1280;
    const int MAX_HEIGHT = 720;

    int newWidth = std::min(width, MAX_WIDTH);
    int newHeight = std::min(height, MAX_HEIGHT);

    // Guard against zero-size windows
    if (newWidth <= 0) newWidth = 1;
    if (newHeight <= 0) newHeight = 1;

    // Ensure any running render is stopped in a controlled way
    stopRender();

    // Final render buffer size follows the (capped) window size
    RENDER_WIDTH = newWidth;
    RENDER_HEIGHT = newHeight;

    // store actual window size for viewport/quad drawing
    windowWidth = width;
    windowHeight = height;

    {
        // Protect buffer reallocation and texture upload against concurrent access from display()
        std::lock_guard<std::mutex> lock(renderMutex);

        delete[] renderBuffer;
        delete[] displayBuffer;
        renderBuffer = new unsigned char[RENDER_WIDTH * RENDER_HEIGHT * 3]();
        displayBuffer = new unsigned char[RENDER_WIDTH * RENDER_HEIGHT * 3]();

        glBindTexture(GL_TEXTURE_2D, textureID);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, RENDER_WIDTH, RENDER_HEIGHT, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);

        // Viewport should cover the actual window size (so textured quad stretches correctly)
        glViewport(0, 0, width, height);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        gluOrtho2D(0, width, 0, height);
        glMatrixMode(GL_MODELVIEW);
    }

    logRenderInfo();

    stopRendering = false;
    startRender(RENDER_WIDTH, RENDER_HEIGHT);
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT);

    // Only update texture when a new rendered image is ready.
    if (renderComplete.load()) {
        std::lock_guard<std::mutex> lock(renderMutex);
        if (renderComplete.load() && displayBuffer) {
            glBindTexture(GL_TEXTURE_2D, textureID);
            glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, RENDER_WIDTH, RENDER_HEIGHT, GL_RGB, GL_UNSIGNED_BYTE, displayBuffer);
            renderComplete = false;
        }
    }

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, textureID);
    glBegin(GL_QUADS);
    glTexCoord2f(0, 0); glVertex2i(0, 0);
    glTexCoord2f(1, 0); glVertex2i(windowWidth, 0);
    glTexCoord2f(1, 1); glVertex2i(windowWidth, windowHeight);
    glTexCoord2f(0, 1); glVertex2i(0, windowHeight);
    glEnd();
    glDisable(GL_TEXTURE_2D);

    glutSwapBuffers();
}

// Обработка нажатий клавиш
void keyboard(unsigned char key, int x, int y) {
    bool changed = false;

    if (key == 27) {
        logDebugInfo("ESC pressed — exiting");
        exit(0);
    }

    if (key == 'w' || key == 'W') { camera.z -= 0.3f; logDebugInfo("W: camera forward"); changed = true; }
    else if (key == 's' || key == 'S') { camera.z += 0.3f; logDebugInfo("S: camera back"); changed = true; }
    else if (key == 'a' || key == 'A') { camera.x -= 0.3f; logDebugInfo("A: camera left"); changed = true; }
    else if (key == 'd' || key == 'D') { camera.x += 0.3f; logDebugInfo("D: camera right"); changed = true; }
    else if (key == 'q' || key == 'Q') { camera.y += 0.3f; logDebugInfo("Q: camera up"); changed = true; }
    else if (key == 'e' || key == 'E') { camera.y -= 0.3f; logDebugInfo("E: camera down"); changed = true; }
    else if (key == 'c' || key == 'C') { camera = Vec3(0, 1, 5); logDebugInfo("C: reset camera"); changed = true; }
    else if (key == 'i' || key == 'I') { lightPos.y += 0.3f; logDebugInfo("I: light up"); changed = true; }
    else if (key == 'k' || key == 'K') { lightPos.y -= 0.3f; logDebugInfo("K: light down"); changed = true; }
    else if (key == 'j' || key == 'J') { lightPos.x -= 0.3f; logDebugInfo("J: light left"); changed = true; }
    else if (key == 'l' || key == 'L') { lightPos.x += 0.3f; logDebugInfo("L: light right"); changed = true; }
    else if (key == 'o' || key == 'O') { lightPos = Vec3(2, 5, 2); logDebugInfo("O: reset light"); changed = true; }
    else if (key == 'f' || key == 'F') { useSSAA = !useSSAA; logDebugInfo("F: toggle SSAA"); changed = true; }
    else if (key == 'r' || key == 'R') { logDebugInfo("R: re-render"); changed = true; }
    else if (key == '1' && scene.size() >= 1) { scene[0]->enabled = !scene[0]->enabled; logDebugInfo("1: toggle sphere"); changed = true; }
    else if (key == '2' && scene.size() >= 2) { scene[1]->enabled = !scene[1]->enabled; logDebugInfo("2: toggle tetrahedron"); changed = true; }
    else if (key == '3' && scene.size() >= 3) { scene[2]->enabled = !scene[2]->enabled; logDebugInfo("3: toggle plane"); changed = true; }

    if (changed) {
        stopRender();
        stopRendering = false;
        startRender(RENDER_WIDTH, RENDER_HEIGHT);
    }
}

void specialKeys(int key, int x, int y) {
    if (key == GLUT_KEY_F5) {
        stopRender();
        loadSceneFromFile("scene.txt");
        stopRendering = false;
        startRender(RENDER_WIDTH, RENDER_HEIGHT);
    }
}

void mouse(int button, int state, int x, int y) {
    if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN) {
        glutCreateMenu(menuCallback);
        glutAddMenuEntry("Toggle SSAA (F)", 1);
        glutAddMenuEntry("Reload Scene (F5)", 2);
        glutAddMenuEntry("Reset Camera (C)", 3);
        glutAddMenuEntry("Reset Light (O)", 5);
        if (scene.size() >= 1) glutAddMenuEntry(scene[0]->enabled ? "Hide Sphere (1)" : "Show Sphere (1)", 10);
        if (scene.size() >= 2) glutAddMenuEntry(scene[1]->enabled ? "Hide Tetrahedron (2)" : "Show Tetrahedron (2)", 11);
        if (scene.size() >= 3) glutAddMenuEntry(scene[2]->enabled ? "Hide Plane (3)" : "Show Plane (3)", 12);
        glutAddMenuEntry("Exit (ESC)", 99);
        glutAttachMenu(GLUT_RIGHT_BUTTON);
    }
}

void menuCallback(int option) {
    bool changed = false;
    switch (option) {
    case 1: useSSAA = !useSSAA; logDebugInfo("Toggle SSAA via menu"); changed = true; break;
    case 2: loadSceneFromFile("scene.txt"); changed = true; break;
    case 3: camera = Vec3(0, 1, 5); logDebugInfo("Reset Camera via menu"); changed = true; break;
    case 5: lightPos = Vec3(2, 5, 2); logDebugInfo("Reset Light via menu"); changed = true; break;
    case 10: if (scene.size() >= 1) { scene[0]->enabled = !scene[0]->enabled; logDebugInfo("Toggle Sphere via menu"); changed = true; } break;
    case 11: if (scene.size() >= 2) { scene[1]->enabled = !scene[1]->enabled; logDebugInfo("Toggle Tetrahedron via menu"); changed = true; } break;
    case 12: if (scene.size() >= 3) { scene[2]->enabled = !scene[2]->enabled; logDebugInfo("Toggle Plane via menu"); changed = true; } break;
    case 99: logDebugInfo("Exit via menu"); exit(0); break;
    }

    if (changed) {
        stopRender();
        stopRendering = false;
        startRender(RENDER_WIDTH, RENDER_HEIGHT);
    }
}

// Функция вызывается, когда система простаивает
void idle() {
    glutPostRedisplay();
	std::this_thread::sleep_for(std::chrono::milliseconds(16)); // примерно 60 FPS
}

// ============================================================================
// MAIN
// ============================================================================

int main(int argc, char** argv) {
    std::printf("=== Raytracer started ===\n");
    std::printf("Controls:\n");
    std::printf("  Camera: W/S/A/D/Q/E — move, C — reset\n");
    std::printf("  Light: I/K/J/L — move, O — reset\n");
    std::printf("  Objects: 1/2/3 — toggle\n");
    std::printf("  Render: F — SSAA, R — re-render, F5 — reload scene\n");
    std::printf("  Menu: Right Mouse Button\n");
    std::printf("  Exit: ESC\n\n");

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(RENDER_WIDTH, RENDER_HEIGHT);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("LR3 - Cross-platform Raytracer");

    glutReshapeFunc(reshape);
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(specialKeys);
    glutMouseFunc(mouse);
    glutIdleFunc(idle);

    loadSceneFromFile("scene.txt");

    renderBuffer = new unsigned char[RENDER_WIDTH * RENDER_HEIGHT * 3]();
    displayBuffer = new unsigned char[RENDER_WIDTH * RENDER_HEIGHT * 3]();

    glGenTextures(1, &textureID);
    glBindTexture(GL_TEXTURE_2D, textureID);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, RENDER_WIDTH, RENDER_HEIGHT, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);

    logCameraInfo();
    logLightInfo();
    logSceneInfo();
    logRenderInfo();

    startRender(RENDER_WIDTH, RENDER_HEIGHT);

    glutMainLoop();

    delete[] renderBuffer;
    delete[] displayBuffer;
    glDeleteTextures(1, &textureID);
    std::printf("=== Raytracer finished ===\n");
    return 0;
}