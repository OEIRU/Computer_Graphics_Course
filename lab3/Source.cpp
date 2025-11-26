// compile
<<<<<<< HEAD
// Linux: g++ -O3 -std=c++17 -o main Source.cpp -lglut -lGL -lGLU -lpthread -lm
// Windows: cl /EHsc /O2 /std:c++17 Source.cpp /link glut32.lib opengl32.lib glu32.lib

#define _CRT_SECURE_NO_WARNINGS
#define NOMINMAX

#ifdef _WIN32
#include <windows.h>
#endif

#include ".\GL\glut.h"
=======
// g++ -O3 -std=c++17 -o main Source.cpp -lglfw -lGL -lpthread 

#define _CRT_SECURE_NO_WARNINGS
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b

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
<<<<<<< HEAD
#include <mutex>
#include <atomic>
#include <condition_variable>
=======

#include <GL/gl.h>
#include <GLFW/glfw3.h>
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b

// ============================================================================
// ВСПОМОГАТЕЛЬНЫЕ СТРУКТУРЫ
// ============================================================================

<<<<<<< HEAD
=======
// 2D-вектор для UV-координат текстур
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
struct Vec2 {
    float x, y;
    Vec2(float x = 0, float y = 0) : x(x), y(y) {}
};

<<<<<<< HEAD
=======
// 3D-вектор: основа для точек, направлений, нормалей и лучей
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
struct Vec3 {
    float x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(float x, float y, float z) : x(x), y(y), z(z) {}

<<<<<<< HEAD
=======
    // Арифметические операции
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
    Vec3 operator+(const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
    Vec3 operator-(const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
    Vec3 operator-() const { return Vec3(-x, -y, -z); }
    Vec3 operator*(float t) const { return Vec3(x * t, y * t, z * t); }
    Vec3 operator/(float t) const { return t != 0 ? Vec3(x / t, y / t, z / t) : Vec3(0, 0, 0); }

<<<<<<< HEAD
    float dot(const Vec3& v) const { return x * v.x + y * v.y + z * v.z; }

=======
    // Скалярное произведение (для углов и проекций)
    float dot(const Vec3& v) const { return x * v.x + y * v.y + z * v.z; }

    // Векторное произведение (для вычисления нормалей)
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
    Vec3 cross(const Vec3& v) const {
        return Vec3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }

    float length2() const { return dot(*this); }
<<<<<<< HEAD
    float length() const { return std::sqrt(length2()); }

=======

    // Длина вектора
    float length() const { return std::sqrt(length2()); }

    // Нормализация: приведение к единичной длине
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
    Vec3 normalize() const {
        float len = length();
        return len > 1e-6f ? (*this) * (1.0f / len) : Vec3(0, 0, 1);
    }

<<<<<<< HEAD
=======
    // Отражение вектора относительно нормали (для зеркал)
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
    Vec3 reflect(const Vec3& n) const {
        return *this - n * (2.0f * this->dot(n));
    }

<<<<<<< HEAD
    Vec3 refract(const Vec3& n, float eta) const {
        float cosI = -n.dot(*this);
        float sinT2 = eta * eta * (1.0f - cosI * cosI);
        if (sinT2 >= 1.0f) return Vec3(0, 0, 0);
=======
    // Преломление (закон Снеллиуса) при переходе между средами
    Vec3 refract(const Vec3& n, float eta) const {
        float cosI = -n.dot(*this); // косинус угла падения
        float sinT2 = eta * eta * (1.0f - cosI * cosI); // sin² угла преломления
        if (sinT2 >= 1.0f) return Vec3(0, 0, 0); // Полное внутреннее отражение
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
        float cosT = std::sqrt(1.0f - sinT2);
        return (*this) * eta + n * (eta * cosI - cosT);
    }
};

struct Color {
    float r, g, b;
    Color() : r(0), g(0), b(0) {}
    Color(float r, float g, float b) {
<<<<<<< HEAD
=======
        // Ограничение значений в [0, 1]
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
        this->r = std::clamp(r, 0.0f, 1.0f);
        this->g = std::clamp(g, 0.0f, 1.0f);
        this->b = std::clamp(b, 0.0f, 1.0f);
    }

<<<<<<< HEAD
=======
    // Операции над цветами
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
    Color operator*(float t) const { return Color(r * t, g * t, b * t); }
    Color operator+(const Color& c) const { return Color(r + c.r, g + c.g, b + c.b); }
    Color operator*(const Color& c) const { return Color(r * c.r, g * c.g, b * c.b); }

<<<<<<< HEAD
=======
    // Линейная интерполяция между двумя цветами
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
    Color blend(const Color& c, float t) const {
        t = std::clamp(t, 0.0f, 1.0f);
        return Color(r * (1 - t) + c.r * t, g * (1 - t) + c.g * t, b * (1 - t) + c.b * t);
    }

    unsigned char getR() const { return (unsigned char)(r * 255); }
    unsigned char getG() const { return (unsigned char)(g * 255); }
    unsigned char getB() const { return (unsigned char)(b * 255); }
};

<<<<<<< HEAD
struct Texture {
    unsigned char* data = nullptr;
    int width = 0, height = 0;
    int imageSize = 0;
    int dataPos = 0;

    ~Texture() { delete[] data; }

    Color sample(float u, float v) const {
        if (!data || width <= 0 || height <= 0) return Color(1, 0, 1);

        u = u - std::floor(u);
        v = v - std::floor(v);

=======
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
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
        int x = (int)(u * width) % width;
        int y = (int)(v * height) % height;
        if (x < 0) x += width;
        if (y < 0) y += height;

<<<<<<< HEAD
        int idx = (y * width + x) * 3;
        return Color(
            data[idx + 2] / 255.0f,
            data[idx + 1] / 255.0f,
            data[idx] / 255.0f
        );
    }

    static Texture* loadBMP(const char* filename) {
        Texture* tex = new Texture();
        FILE* f = std::fopen(filename, "rb");
        if (!f) {
            std::printf("ERROR: Cannot open texture file %s\n", filename);
            delete tex;
            return nullptr;
        }

=======
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
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
        unsigned char header[54];
        if (std::fread(header, 1, 54, f) != 54) {
            std::printf("ERROR: Not a valid BMP file %s\n", filename);
            std::fclose(f);
<<<<<<< HEAD
            delete tex;
            return nullptr;
        }

        if (header[0] != 'B' || header[1] != 'M') {
            std::printf("ERROR: Not a valid BMP file %s\n", filename);
            std::fclose(f);
            delete tex;
            return nullptr;
        }

        tex->dataPos = *(int*)&(header[0x0A]);
        tex->imageSize = *(int*)&(header[0x22]);
        tex->width = *(int*)&(header[0x12]);
        tex->height = *(int*)&(header[0x16]);

        if (tex->imageSize == 0) tex->imageSize = tex->width * tex->height * 3;
        if (tex->dataPos == 0) tex->dataPos = 54;

        if (*(short*)&(header[0x1C]) != 24) {
            std::printf("ERROR: Only 24-bit BMP supported %s\n", filename);
            std::fclose(f);
            delete tex;
            return nullptr;
        }

        tex->data = new unsigned char[tex->imageSize];
        std::fseek(f, tex->dataPos, SEEK_SET);
        size_t read = std::fread(tex->data, 1, tex->imageSize, f);
        std::fclose(f);

        if (read != (size_t)tex->imageSize) {
            std::printf("ERROR: Could not read full texture data %s\n", filename);
            delete[] tex->data;
            tex->data = nullptr;
            delete tex;
            return nullptr;
        }
        else {
            std::printf("Texture loaded: %s (%dx%d)\n", filename, tex->width, tex->height);
=======
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
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
        }

        return tex;
    }
};

<<<<<<< HEAD
struct Material {
    Color color;
    float ka, kd, ks;
    float reflect;
    float transparency;
    int shininess;
    Texture* texture;

    Material() : color(1, 1, 1), ka(0.1f), kd(0.7f), ks(0.3f),
        reflect(0.0f), transparency(0.0f), shininess(32), texture(nullptr) {
    }

=======
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
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
    Material(Color c, float a, float d, float s, float r, float t, int sh, Texture* tex = nullptr)
        : color(c), ka(a), kd(d), ks(s), reflect(r), transparency(t), shininess(sh), texture(tex) {}
};

<<<<<<< HEAD
struct Ray {
    Vec3 origin, direction;
    Ray(Vec3 o, Vec3 d) : origin(o), direction(d.normalize()) {}
=======
// ============================================================================
// ЛУЧ И ИНФОРМАЦИЯ О ПЕРЕСЕЧЕНИИ
// ============================================================================

// Луч: точка начала и направление
struct Ray {
    Vec3 origin, direction;
    Ray(Vec3 o, Vec3 d) : origin(o), direction(d.normalize()) {}
    // Точка на луче на расстоянии t
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
    Vec3 at(float t) const { return origin + direction * t; }
};

// Результат пересечения луча с объектом
struct HitInfo {
<<<<<<< HEAD
    bool hit = false;
    float t = 1e30f;
    Vec3 point, normal;
    Material mat;
    float u = 0, v = 0;
    HitInfo() {}
};

struct Shape {
    bool enabled = true;
    virtual ~Shape() = default;
=======
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
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
    virtual HitInfo intersect(const Ray& ray) const = 0;
    virtual bool isEnabled() const { return enabled; }
};

<<<<<<< HEAD
=======
// ============================================================================
// СФЕРА
// ============================================================================

>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
struct Sphere : public Shape {
    Vec3 center; float radius; Material mat;
    Sphere(Vec3 c, float r, Material m, bool en = true) : center(c), radius(r), mat(m) { enabled = en; }

    HitInfo intersect(const Ray& ray) const override {
        if (!isEnabled()) return HitInfo();

<<<<<<< HEAD
=======
        // Вектор от центра сферы к началу луча
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
        Vec3 oc = ray.origin - center;

        // Коэффициенты квадратного уравнения |O + tD - C|² = R²
        float a = 1.0f;
        float b = 2.0f * oc.dot(ray.direction);
        float c = oc.dot(oc) - radius * radius;
        float disc = b * b - 4 * a * c;
<<<<<<< HEAD

        if (disc < 0) return HitInfo();

        float sqrt_disc = std::sqrt(disc);
        float t = (-b - sqrt_disc) / (2 * a);

=======

        if (disc < 0) return HitInfo(); // Нет пересечения

        // Вычисление корней
        float sqrt_disc = std::sqrt(disc);
        float t = (-b - sqrt_disc) / (2 * a);

        // Если ближайшая точка за камерой — пробуем дальнюю
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
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

<<<<<<< HEAD
        float theta = std::acos(std::clamp(-info.normal.y, -1.0f, 1.0f));
        float phi = std::atan2(info.normal.z, info.normal.x);
        info.u = (phi + 3.1415926535f) / (2 * 3.1415926535f);
        info.v = theta / 3.1415926535f;
=======
        // Преобразование нормали в UV (сферическая проекция)
        float theta = std::acos(std::clamp(-info.normal.y, -1.0f, 1.0f)); // Полярный угол
        float phi = std::atan2(info.normal.z, info.normal.x);             // Азимут
        info.u = (phi + 3.1415926535f) / (2 * 3.1415926535f);            // [0,1]
        info.v = theta / 3.1415926535f;                                   // [0,1]
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b

        return info;
    }
};

<<<<<<< HEAD
=======
// ============================================================================
// ТРЕУГОЛЬНИК (АЛГОРИТМ МЁЛЛЕРА–ТРУМБОРА)
// ============================================================================

>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
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

<<<<<<< HEAD
        if (std::abs(a) < 1e-6f) return HitInfo();
=======
        if (std::abs(a) < 1e-6f) return HitInfo(); // Луч параллелен плоскости
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b

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

<<<<<<< HEAD
=======
        // Коррекция нормали: должна смотреть навстречу лучу
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
        if (info.normal.dot(ray.direction) > 0) info.normal = -info.normal;

        info.mat = mat;

<<<<<<< HEAD
=======
        // Барицентрическая интерполяция UV
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
        float w = 1.0f - u - v;
        info.u = w * uv0.x + u * uv1.x + v * uv2.x;
        info.v = w * uv0.y + u * uv1.y + v * uv2.y;

        return info;
    }
};

<<<<<<< HEAD
=======
// ============================================================================
// ТЕТРАЭДР (СОСТОИТ ИЗ 4 ТРЕУГОЛЬНИКОВ)
// ============================================================================

>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
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
<<<<<<< HEAD
=======
        // Определение 4 граней тетраэдра
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
        Triangle faces[4] = {
            Triangle(v[0], v[1], v[2], Vec2(0,0), Vec2(1,0), Vec2(0,1), mat, true),
            Triangle(v[0], v[2], v[3], Vec2(0,0), Vec2(0,1), Vec2(1,1), mat, true),
            Triangle(v[0], v[3], v[1], Vec2(0,0), Vec2(1,1), Vec2(1,0), mat, true),
            Triangle(v[1], v[3], v[2], Vec2(1,0), Vec2(1,1), Vec2(0,1), mat, true)
        };

<<<<<<< HEAD
        for (auto& tri : faces) {
            HitInfo h = tri.intersect(ray);
            if (h.hit && h.t < best.t) {
=======
        // Проверка пересечения с каждой гранью
        for (auto& tri : faces) {
            HitInfo h = tri.intersect(ray);
            if (h.hit && h.t < best.t) {
                // Коррекция нормали
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
                if (h.normal.dot(ray.direction) > 0) {
                    h.normal = -h.normal;
                }
                best = h;
            }
        }
        return best;
    }
};

<<<<<<< HEAD
=======
// ============================================================================
// БЕСКОНЕЧНАЯ ПЛОСКОСТЬ
// ============================================================================

>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
struct Plane : public Shape {
    Vec3 point, normal; Material mat;
    Plane(Vec3 p, Vec3 n, Material m, bool en = true) : point(p), normal(n.normalize()), mat(m) { enabled = en; }

    HitInfo intersect(const Ray& ray) const override {
        if (!isEnabled()) return HitInfo();

        float denom = normal.dot(ray.direction);
<<<<<<< HEAD
        if (std::abs(denom) < 1e-6f) return HitInfo();
=======
        if (std::abs(denom) < 1e-6f) return HitInfo(); // Луч параллелен плоскости
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b

        float t = (point - ray.origin).dot(normal) / denom;
        if (t < 0.001f) return HitInfo();

        HitInfo info;
        info.hit = true;
        info.t = t;
        info.point = ray.at(t);
<<<<<<< HEAD
        info.normal = denom > 0 ? -normal : normal;
        info.mat = mat;
=======
        // Нормаль направлена навстречу лучу
        info.normal = denom > 0 ? -normal : normal;
        info.mat = mat;
        // Простая проекция UV по координатам X и Z
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
        info.u = info.point.x * 2.0f;
        info.v = info.point.z * 2.0f;

        return info;
    }
};

// ============================================================================
// ГЛОБАЛЬНЫЕ ПЕРЕМЕННЫЕ
// ============================================================================
<<<<<<< HEAD

int RENDER_WIDTH = 800;
int RENDER_HEIGHT = 600;
unsigned char* renderBuffer = nullptr;
unsigned char* displayBuffer = nullptr;
GLuint textureID = 0;

Vec3 camera(0, 1, 5);
Vec3 lightPos(2, 5, 2);
bool useSSAA = false;
std::vector<std::unique_ptr<Shape>> scene;
const int MAX_DEPTH = 5;
const float FOV = 60.0f;
=======

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
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b

struct RenderParams {
    float invWidth, invHeight, aspect, tanHalfFov;
    void update(int w, int h, float fov) {
        invWidth = 1.0f / w;
        invHeight = 1.0f / h;
        aspect = (float)w / h;
        tanHalfFov = std::tan(fov * 0.5f * 3.1415926535f / 180.0f);
    }
} g_renderParams;

<<<<<<< HEAD
std::mutex renderMutex;
std::atomic<bool> renderingInProgress(false);
std::atomic<bool> renderComplete(false);
std::atomic<int> currentRenderRow(0);
std::atomic<bool> stopRendering(false);

=======
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
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
<<<<<<< HEAD
    std::printf("[RENDER] SSAA: %s, Resolution: %dx%d\n", useSSAA ? "ON" : "OFF", RENDER_WIDTH, RENDER_HEIGHT);
=======
    std::printf("[RENDER] SSAA: %s\n", useSSAA ? "ON" : "OFF");
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
}

// ============================================================================
// МОДЕЛЬ ОСВЕЩЕНИЯ COOK-TORRANCE
<<<<<<< HEAD
// ============================================================================

Color cookTorrance(const Vec3& N, const Vec3& V, const Vec3& L, const Color& baseColor, float ks, int shininess) {
    Vec3 H = (L + V).normalize();
=======
// Физически корректный specular
// ============================================================================

Color cookTorrance(const Vec3& N, const Vec3& V, const Vec3& L, const Color& baseColor, float ks, int shininess) {
    Vec3 H = (L + V).normalize(); // Half-vector
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b

    float NdotL = std::max(0.0f, N.dot(L));
    float NdotV = std::max(0.0f, N.dot(V));
    float NdotH = std::max(0.0f, N.dot(H));
    float VdotH = std::max(0.0f, V.dot(H));

    if (NdotL <= 0 || NdotV <= 0) return Color(0, 0, 0);

<<<<<<< HEAD
    float F0 = 0.04f;
    float F = F0 + (1.0f - F0) * std::pow(1.0f - VdotH, 5.0f);

    float alpha = std::max(0.001f, 1.0f - (float)shininess / 100.0f);
    float alpha2 = alpha * alpha;
    float denomD = NdotH * NdotH * (alpha2 - 1.0f) + 1.0f;
    float D = alpha2 / (3.1415926535f * denomD * denomD);

=======
    // Френель (приближение Шлика)
    float F0 = 0.04f; // Базовый коэффициент отражения для диэлектриков
    float F = F0 + (1.0f - F0) * std::pow(1.0f - VdotH, 5.0f);

    // Распределение микрограней (GGX)
    float alpha = std::max(0.001f, 1.0f - (float)shininess / 100.0f);
    float alpha2 = alpha * alpha;
    float denomD = NdotH * NdotH * (alpha2 - 1.0f) + 1.0f;
    float D = alpha2 / (3.1415926535f * denomD * denomD);

    // Геометрическое затенение (Smith model)
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
    float k = (alpha + 1.0f) * (alpha + 1.0f) / 8.0f;
    float G1V = NdotV / (NdotV * (1.0f - k) + k);
    float G1L = NdotL / (NdotL * (1.0f - k) + k);
    float G = G1V * G1L;

<<<<<<< HEAD
=======
    // Итоговый specular
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
    float specular = (F * D * G) / std::max(0.001f, 4.0f * NdotV * NdotL);

    return Color(1, 1, 1) * ks * specular * NdotL;
}

// ============================================================================
// РЕКУРСИВНАЯ ТРАССИРОВКА ЛУЧЕЙ
// ============================================================================
<<<<<<< HEAD

Color traceRay(const Ray& ray, int depth = 0) {
    if (depth >= MAX_DEPTH) return Color(0.05f, 0.05f, 0.1f);
=======
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b

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

<<<<<<< HEAD
    if (!closest.hit) return Color(0.05f, 0.05f, 0.05f);

=======
    if (!closest.hit) return Color(0.05f, 0.05f, 0.05f); // Серый фон

    // Основной цвет: либо базовый, либо из текстуры
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
    Color baseColor = closest.mat.color;
    if (closest.mat.texture) {
        baseColor = closest.mat.texture->sample(closest.u, closest.v);
    }

    Vec3 L_vec = lightPos - closest.point;
    float lightDist2 = L_vec.length2();
    if (lightDist2 < 1e-6f) {
        return baseColor * closest.mat.ka;
    }
<<<<<<< HEAD

    float lightDist = std::sqrt(lightDist2);
    Vec3 L = L_vec / lightDist;
    Vec3 V = (camera - closest.point).normalize();
    Vec3 N = closest.normal;

    bool inShadow = false;
    Ray shadowRay(closest.point + N * 0.01f, L);
=======

    float lightDist = std::sqrt(lightDist2);
    Vec3 L = L_vec / lightDist; // Нормализованное направление к свету
    Vec3 V = (camera - closest.point).normalize(); // Направление к камере
    Vec3 N = closest.normal; // Нормаль в точке

    // Проверка тени
    bool inShadow = false;
    Ray shadowRay(closest.point + N * 0.01f, L); // Небольшое смещение, чтобы избежать self-intersection
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
    for (const auto& s : scene) {
        if (!s->isEnabled()) continue;
        HitInfo h = s->intersect(shadowRay);
        if (h.hit && h.t > 0.001f && h.t < lightDist - 0.01f) {
            inShadow = true;
            break;
        }
    }

<<<<<<< HEAD
=======
    // Ambient
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
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

<<<<<<< HEAD
=======
    // Рекурсивное отражение
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
    if (closest.mat.reflect > 0.01f && depth < MAX_DEPTH - 1) {
        Vec3 reflDir = ray.direction.reflect(N);
        Ray reflRay(closest.point + N * 0.01f, reflDir);
        Color reflCol = traceRay(reflRay, depth + 1);
        result = result * (1.0f - closest.mat.reflect) + reflCol * closest.mat.reflect;
    }

<<<<<<< HEAD
    if (closest.mat.transparency > 0.01f && depth < MAX_DEPTH - 1) {
        float eta = N.dot(ray.direction) < 0 ? (1.0f / 1.5f) : 1.5f;
=======
    // Рекурсивное преломление
    if (closest.mat.transparency > 0.01f && depth < MAX_DEPTH - 1) {
        // Вычисление коэффициента преломления
        float eta = N.dot(ray.direction) < 0 ? (1.0f / 1.5f) : 1.5f; // 1.5 — IOR стекла
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
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
<<<<<<< HEAD
// ПОТОК РЕНДЕРИНГА — БЕЗОПАСНЫЙ
// ============================================================================
=======
// РЕНДЕРИНГ ИЗОБРАЖЕНИЯ
// ============================================================================

void renderImage() {
    logDebugInfo("Start rendering");
    auto startTime = std::chrono::high_resolution_clock::now();
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b

void renderThreadFunction(int renderW, int renderH) {
    renderingInProgress = true;
    renderComplete = false;
    stopRendering = false;
    currentRenderRow = 0;

    const int localWidth = renderW;
    const int localHeight = renderH;
    g_renderParams.update(localWidth, localHeight, FOV);

    const int NUM_THREADS = std::thread::hardware_concurrency();
    const int THREADS = (NUM_THREADS > 0) ? NUM_THREADS : 4;
<<<<<<< HEAD
    const bool ssaa = useSSAA;

    float invWidth = g_renderParams.invWidth;
    float invHeight = g_renderParams.invHeight;
    float aspect = g_renderParams.aspect;
    float tanHalfFov = g_renderParams.tanHalfFov;

    unsigned char* tempBuffer = new unsigned char[localWidth * localHeight * 3]();
=======
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b

    std::vector<std::thread> threads;
    threads.reserve(THREADS);

<<<<<<< HEAD
=======
    int rowsPerThread = RENDER_HEIGHT / THREADS;
    const bool ssaa = useSSAA;

    // Создание потоков для рендеринга полос изображения
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
    for (int t = 0; t < THREADS; ++t) {
        threads.emplace_back([&, t]() {
            while (true) {
                if (stopRendering) break;

                int row;
                {
                    std::lock_guard<std::mutex> lock(renderMutex);
                    if (currentRenderRow >= localHeight) break;
                    row = currentRenderRow++;
                }

                for (int x = 0; x < localWidth; ++x) {
                    Color finalColor(0, 0, 0);
                    int samples = ssaa ? 4 : 1;

                    for (int sy = 0; sy < (ssaa ? 2 : 1); ++sy) {
                        for (int sx = 0; sx < (ssaa ? 2 : 1); ++sx) {
                            float jitterX = ssaa ? (sx + 0.5f) / 2.0f : 0.5f;
                            float jitterY = ssaa ? (sy + 0.5f) / 2.0f : 0.5f;

<<<<<<< HEAD
=======
                            // Преобразование пикселя в направление луча
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
                            float px = (2.0f * (x + jitterX) * invWidth - 1.0f) * aspect * tanHalfFov;
                            float py = (1.0f - 2.0f * (row + jitterY) * invHeight) * tanHalfFov;
                            Vec3 rayDir(px, py, -1);
                            Ray ray(camera, rayDir);

                            finalColor = finalColor + traceRay(ray);
                        }
                    }

                    finalColor = finalColor * (1.0f / samples);

<<<<<<< HEAD
                    int idx = ((localHeight - 1 - row) * localWidth + x) * 3;
                    if (idx >= 0 && idx + 2 < localWidth * localHeight * 3) {
                        tempBuffer[idx] = finalColor.getR();
                        tempBuffer[idx + 1] = finalColor.getG();
                        tempBuffer[idx + 2] = finalColor.getB();
                    }
=======
                    // Запись в буфер (с переворотом по Y для OpenGL)
                    int idx = ((RENDER_HEIGHT - 1 - y) * RENDER_WIDTH + x) * 3;
                    image[idx] = finalColor.getR();
                    image[idx + 1] = finalColor.getG();
                    image[idx + 2] = finalColor.getB();
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
                }
            }
        });
    }

    // Ожидание завершения всех потоков
    for (auto& t : threads) {
        if (t.joinable()) t.join();
    }

<<<<<<< HEAD
    {
        std::lock_guard<std::mutex> lock(renderMutex);
        if (RENDER_WIDTH == localWidth && RENDER_HEIGHT == localHeight) {
            std::memcpy(displayBuffer, tempBuffer, localWidth * localHeight * 3);
            renderComplete = true;
        }
    }

    delete[] tempBuffer;
    renderingInProgress = false;
}

// ============================================================================
// ЗАГРУЗКА СЦЕНЫ
=======
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
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
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
<<<<<<< HEAD
        if (line.empty() || line[0] == '#') continue;
=======
        if (line.empty() || line[0] == '#') continue; // Пропуск комментариев
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b

        std::istringstream iss(line);
        std::string type;
        iss >> type;

        if (type == "sphere") {
            float cx, cy, cz, r, R, G, B, ka, kd, ks, refl, trans, shininess;
            std::string texName;
            iss >> cx >> cy >> cz >> r >> R >> G >> B >> ka >> kd >> ks >> refl >> trans >> shininess >> texName;
            Color col(R, G, B);
            Texture* tex = (texName != "none") ? Texture::loadBMP(texName.c_str()) : nullptr;
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
            Texture* tex = (texName != "none") ? Texture::loadBMP(texName.c_str()) : nullptr;
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
            Texture* tex = (texName != "none") ? Texture::loadBMP(texName.c_str()) : nullptr;
            scene.push_back(std::make_unique<Plane>(Vec3(px, py, pz), Vec3(nx, ny, nz), Material(col, ka, kd, ks, refl, trans, (int)shininess, tex)));
            objectsLoaded++;
            logDebugInfo("Plane loaded at (" + std::to_string(px) + ", " + std::to_string(py) + ", " + std::to_string(pz) + ")");
        }
    }
    logDebugInfo("Scene loaded: " + std::to_string(objectsLoaded) + " objects");
}

// ============================================================================
<<<<<<< HEAD
// ОБНОВЛЕНИЕ МАТЕРИАЛА
=======
// ОБНОВЛЕНИЕ ПАРАМЕТРОВ МАТЕРИАЛА
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
// ============================================================================

void updateMaterial(Material& mat, int action) {
    const float delta = 0.05f;
    switch (action) {
<<<<<<< HEAD
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
=======
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
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
    }
}

// ============================================================================
<<<<<<< HEAD
// GLUT CALLBACKS
// ============================================================================

// Прототип
void menuCallback(int option);

void reshape(int width, int height) {
    const int MAX_WIDTH = 1280;
    const int MAX_HEIGHT = 720;

    int newWidth = std::min(width, MAX_WIDTH);
    int newHeight = std::min(height, MAX_HEIGHT);

    if (renderingInProgress) {
        stopRendering = true;
        while (renderingInProgress) {
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
    }

    RENDER_WIDTH = newWidth;
    RENDER_HEIGHT = newHeight;

    {
        std::lock_guard<std::mutex> lock(renderMutex);
        delete[] renderBuffer;
        delete[] displayBuffer;
        renderBuffer = new unsigned char[RENDER_WIDTH * RENDER_HEIGHT * 3]();
        displayBuffer = new unsigned char[RENDER_WIDTH * RENDER_HEIGHT * 3]();

        glBindTexture(GL_TEXTURE_2D, textureID);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, RENDER_WIDTH, RENDER_HEIGHT, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);
    }

    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, width, 0, height);
    glMatrixMode(GL_MODELVIEW);

    logRenderInfo();

    stopRendering = false;
    std::thread(renderThreadFunction, RENDER_WIDTH, RENDER_HEIGHT).detach();
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT);

    if (renderComplete) {
        {
            std::lock_guard<std::mutex> lock(renderMutex);
            glBindTexture(GL_TEXTURE_2D, textureID);
            glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, RENDER_WIDTH, RENDER_HEIGHT, GL_RGB, GL_UNSIGNED_BYTE, displayBuffer);
            renderComplete = false;
        }
    }

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, textureID);
    glBegin(GL_QUADS);
    glTexCoord2f(0, 0); glVertex2i(0, 0);
    glTexCoord2f(1, 0); glVertex2i(RENDER_WIDTH, 0);
    glTexCoord2f(1, 1); glVertex2i(RENDER_WIDTH, RENDER_HEIGHT);
    glTexCoord2f(0, 1); glVertex2i(0, RENDER_HEIGHT);
    glEnd();
    glDisable(GL_TEXTURE_2D);

    glutSwapBuffers();
}

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
        if (renderingInProgress) {
            stopRendering = true;
            while (renderingInProgress) {
                std::this_thread::sleep_for(std::chrono::milliseconds(10));
            }
        }
        std::thread(renderThreadFunction, RENDER_WIDTH, RENDER_HEIGHT).detach();
    }
}

void specialKeys(int key, int x, int y) {
    if (key == GLUT_KEY_F5) {
        if (renderingInProgress) {
            stopRendering = true;
            while (renderingInProgress) {
                std::this_thread::sleep_for(std::chrono::milliseconds(10));
            }
        }
        loadSceneFromFile("scene.txt");
        std::thread(renderThreadFunction, RENDER_WIDTH, RENDER_HEIGHT).detach();
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
        if (renderingInProgress) {
            stopRendering = true;
            while (renderingInProgress) {
                std::this_thread::sleep_for(std::chrono::milliseconds(10));
            }
        }
        std::thread(renderThreadFunction, RENDER_WIDTH, RENDER_HEIGHT).detach();
    }
}

void idle() {
    glutPostRedisplay();
    std::this_thread::sleep_for(std::chrono::milliseconds(16));
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
    glutCreateWindow("LR3 - Cross-platform Raytracer with GLUT & Multithreading");

    glutReshapeFunc(reshape);
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(specialKeys);
    glutMouseFunc(mouse);
    glutIdleFunc(idle);

=======
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
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
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

<<<<<<< HEAD
    std::thread(renderThreadFunction, RENDER_WIDTH, RENDER_HEIGHT).detach();

    glutMainLoop();

    delete[] renderBuffer;
    delete[] displayBuffer;
    glDeleteTextures(1, &textureID);
    std::printf("=== Raytracer finished ===\n");
=======
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
>>>>>>> 17812082491c034682996f28a9fa28b691d9244b
    return 0;
}