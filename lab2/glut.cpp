// ===================================================================================
// ПРЕПРОЦЕССОРНЫЕ ДИРЕКТИВЫ
// ===================================================================================

#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES

// ===================================================================================
// СТАНДАРТНЫЕ ЗАГОЛОВОЧНЫЕ ФАЙЛЫ
// ===================================================================================

#include <iostream>      
#include <vector>       
#include <string>        
#include <fstream>      
#include <cmath>         
#include <cstdio>        
#include <algorithm>     
#include <unordered_map> 

// ===================================================================================
// ГРАФИЧЕСКИЕ БИБЛИОТЕКИ
// ===================================================================================

#include "C:\Users\Vladi\source\repos\KG2-REBUILD\GL\glut.h"
#include <GL/glu.h>
#include <GL/gl.h>

// Подключение библиотек Windows (для линковки с OpenGL и GLUT)
#pragma comment(lib, "glut32.lib")
#pragma comment(lib, "opengl32.lib")
#pragma comment(lib, "glu32.lib")

// ===================================================================================
// МАТЕМАТИЧЕСКИЕ СТРУКТУРЫ
// ===================================================================================

// Вектор в 3D-пространстве (используется для вершин, нормалей и т.д.)
struct Vector3 {
    double x = 0.0, y = 0.0, z = 0.0;

    // Оператор < нужен для использования Vector3 в std::map
    bool operator<(const Vector3& o) const {
        if (x != o.x) return x < o.x;
        if (y != o.y) return y < o.y;
        return z < o.z;
    }
};

// Хэш-функция для Vector3 — позволяет использовать его как ключ в unordered_map
struct Vector3Hash {
    std::size_t operator()(const Vector3& v) const {
        std::hash<double> hasher;
        size_t h1 = hasher(v.x);
        size_t h2 = hasher(v.y);
        size_t h3 = hasher(v.z);
        // Простой способ комбинирования хэшей
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

// Функция сравнения векторов с учётом погрешности
struct Vector3Equal {
    bool operator()(const Vector3& a, const Vector3& b) const {
        const double epsilon = 1e-9; // Допустимая погрешность
        return std::abs(a.x - b.x) < epsilon &&
            std::abs(a.y - b.y) < epsilon &&
            std::abs(a.z - b.z) < epsilon;
    }
};

// Треугольная грань (face) — состоит из трёх вершин
struct Face {
    Vector3 verts[3]; 
};

// ===================================================================================
// УТИЛИТЫ: МАТЕМАТИЧЕСКИЕ ФУНКЦИИ
// ===================================================================================

// Сложение двух векторов
Vector3 addVec(const Vector3& a, const Vector3& b) {
    return Vector3{ a.x + b.x, a.y + b.y, a.z + b.z };
}

// Нормализация вектора (приведение к длине 1)
Vector3 normalizeVec(const Vector3& a) {
    double len = std::sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
    if (len < 1e-8) return Vector3{ 0, 0, 0 }; // Защита от деления на ноль
    return Vector3{ a.x / len, a.y / len, a.z / len };
}

// Вычисление нормали к треугольнику по трём его вершинам (векторное произведение)
Vector3 calcNormal(const Vector3& v1, const Vector3& v2, const Vector3& v3) {
    // Векторы двух сторон треугольника
    double ux = v2.x - v1.x, uy = v2.y - v1.y, uz = v2.z - v1.z;
    double vx = v3.x - v1.x, vy = v3.y - v1.y, vz = v3.z - v1.z;
    // Векторное произведение → нормаль
    Vector3 n{
        uy * vz - uz * vy,
        uz * vx - ux * vz,
        ux * vy - uy * vx
    };
    return normalizeVec(n); // Возвращаем единичный вектор
}

// ===================================================================================
// ЗАГРУЗЧИК ФАЙЛОВ
// ===================================================================================

// Класс для загрузки числовых данных из текстовых файлов
class FileManager {
public:
    // Загружает все double-значения из файла в вектор out
    // expectedFormat: "even" — чётное количество чисел (для 2D-сечений), "triples" — кратно 3 (для 3D-траекторий)
    static bool loadData(const std::string& path, std::vector<double>& out, const std::string& expectedFormat = "") {
        FILE* f = std::fopen(path.c_str(), "r");
        if (!f) {
            std::cerr << "Файл не найден: " << path << "\n";
            return false;
        }
        double value;
        int count = 0;
        // Читаем все числа из файла
        while (std::fscanf(f, "%lf", &value) == 1) {
            out.push_back(value);
            count++;
        }
        std::fclose(f);

        // Проверка формата данных
        if (expectedFormat == "even" && count % 2 != 0) {
            std::cerr << "Чётное число координат в " << path << "\n";
            return false;
        }
        if (expectedFormat == "triples" && count % 3 != 0) {
            std::cerr << "Координаты по 3 в " << path << "\n";
            return false;
        }

        std::cout << "Загружено " << count << " значений из " << path << "\n";
        return count > 0;
    }
};

// ===================================================================================
// ГЕОМЕТРИЯ: ПОСТРОЕНИЕ 3D-МОДЕЛИ
// ===================================================================================

class Geometry {
public:
    std::vector<Face> meshFaces;           // Список треугольников модели
    std::vector<Vector3> flatNormals;      // Нормали по граням (плоское освещение)
    std::vector<Vector3> smoothNormals;    // Нормали по вершинам (сглаженное освещение)
    Vector3 center{ 0, 0, 0 };             // Центр модели (среднее по всем вершинам)
    Vector3 sectionCenter{ 0, 0, 0 };      // Центр 2D-сечения (профиля)
    double heightRange = 1.0;              // Высота модели по оси Z (для текстурных координат)

    // Вычисляет центр модели
    void computeCenter();
    // Вычисляет нормали (плоские и сглаженные)
    void computeNormals();
    // Строит 3D-модель по сечению, траектории и масштабам
    void buildFromData(const std::vector<double>& section,
        const std::vector<double>& trajectory,
        const std::vector<double>& scaling);

private:
    // Создаёт поперечное сечение вдоль траектории
    std::vector<Vector3> createRing(const std::vector<double>& section,
        double tx, double ty, double tz, double scale);
    // Соединяет сечения треугольниками и добавляет "крышки" на концах
    void triangulateRings(const std::vector<std::vector<Vector3>>& rings);
};

// Вычисление центра модели как среднего по всем вершинам
void Geometry::computeCenter() {
    if (meshFaces.empty()) {
        center = { 0,0,0 };
        return;
    }
    double sx = 0, sy = 0, sz = 0;
    int cnt = 0;
    for (const auto& f : meshFaces)
        for (const auto& v : f.verts) {
            sx += v.x; sy += v.y; sz += v.z; ++cnt;
        }
    center = { sx / cnt, sy / cnt, sz / cnt };
}

// Вычисление нормалей: плоских и сглаженных
void Geometry::computeNormals() {
    flatNormals.clear();
    smoothNormals.clear();
    flatNormals.resize(meshFaces.size());

    // Хэш-таблица: вершина → сумма нормалей соседних граней
    std::unordered_map<Vector3, Vector3, Vector3Hash, Vector3Equal> tempNormals;

    // Шаг 1: вычисляем плоские нормали и накапливаем сглаженные
    for (size_t i = 0; i < meshFaces.size(); ++i) {
        const auto& f = meshFaces[i];
        Vector3 n = calcNormal(f.verts[0], f.verts[1], f.verts[2]);
        flatNormals[i] = n;
        // Добавляем нормаль к каждой вершине треугольника
        for (int j = 0; j < 3; ++j) {
            tempNormals[f.verts[j]] = addVec(tempNormals[f.verts[j]], n);
        }
    }

    // Шаг 2: нормализуем сглаженные нормали
    smoothNormals.resize(meshFaces.size() * 3);
    for (size_t i = 0; i < meshFaces.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            const Vector3& v = meshFaces[i].verts[j];
            auto it = tempNormals.find(v);
            smoothNormals[i * 3 + j] = (it != tempNormals.end()) ? normalizeVec(it->second) : Vector3{ 0,0,0 };
        }
    }
}

// Создание одного поперечного сечения в точке траектории
std::vector<Vector3> Geometry::createRing(const std::vector<double>& section,
    double tx, double ty, double tz, double scale) {
    int profilePoints = section.size() / 2; // Количество точек в профиле
    std::vector<Vector3> ring;
    ring.reserve(profilePoints);
    for (int k = 0; k < profilePoints; ++k) {
        // Смещаем точку относительно центра сечения, масштабируем и перемещаем вдоль траектории
        double dx = section[2 * k] - sectionCenter.x;
        double dy = section[2 * k + 1] - sectionCenter.y;
        ring.emplace_back(Vector3{
            dx * scale + sectionCenter.x + tx,
            dy * scale + sectionCenter.y + ty,
            tz // Z-координата берётся из траектории
            });
    }
    return ring;
}

// создаём боковую поверхность и "крышки" на концах
void Geometry::triangulateRings(const std::vector<std::vector<Vector3>>& rings) {
    meshFaces.clear();

    // Боковая поверхность: соединяем соседние кольца четырёхугольниками → два треугольника
    for (size_t i = 0; i < rings.size() - 1; ++i) {
        const auto& ring0 = rings[i];
        const auto& ring1 = rings[i + 1];
        int ringSize = ring0.size();
        for (int j = 0; j < ringSize; ++j) {
            int j1 = (j + 1) % ringSize;
            // Два треугольника на один четырёхугольник
            meshFaces.push_back(Face{ ring0[j], ring1[j], ring1[j1] });
            meshFaces.push_back(Face{ ring0[j], ring1[j1], ring0[j1] });
        }
    }

    // Крышка на начале (если есть сечение)
    if (!rings.empty() && !rings[0].empty()) {
        const auto& first = rings[0];
        // Центр крышки — среднее по сечению
        Vector3 center{ 0,0,0 };
        for (const auto& v : first) center = addVec(center, v);
        center.x /= first.size(); center.y /= first.size(); center.z /= first.size();
        // Треугольники от центра к краям
        for (size_t i = 0; i < first.size(); ++i) {
            size_t i1 = (i + 1) % first.size();
            meshFaces.push_back(Face{ center, first[i], first[i1] });
        }
    }

    // Крышка на конце
    if (rings.size() > 1 && !rings.back().empty()) {
        const auto& last = rings.back();
        Vector3 center{ 0,0,0 };
        for (const auto& v : last) center = addVec(center, v);
        center.x /= last.size(); center.y /= last.size(); center.z /= last.size();
        for (size_t i = 0; i < last.size(); ++i) {
            size_t i1 = (i + 1) % last.size();
            meshFaces.push_back(Face{ center, last[i], last[i1] });
        }
    }
}

// Главный метод построения геометрии
void Geometry::buildFromData(const std::vector<double>& section,
    const std::vector<double>& trajectory,
    const std::vector<double>& scaling) {

    // Проверка входных данных
    if (section.size() % 2 != 0 || section.size() < 6)
        throw std::runtime_error("Сечение: минимум 3 точки (6 координат)");
    const int profilePoints = section.size() / 2;
    const int trajPoints = trajectory.size() / 3;
    if (trajPoints < 2)
        throw std::runtime_error("Траектория: минимум 2 точки");
    if (trajPoints != (int)scaling.size())
        throw std::runtime_error("Несовпадение траектории и масштабов");

    // Вычисляем центр 2D-сечения
    sectionCenter = { 0,0,0 };
    for (int k = 0; k < profilePoints; ++k) {
        sectionCenter.x += section[2 * k];
        sectionCenter.y += section[2 * k + 1];
    }
    sectionCenter.x /= profilePoints;
    sectionCenter.y /= profilePoints;

    // Создаём кольца вдоль траектории
    std::vector<std::vector<Vector3>> rings;
    for (int i = 0; i < trajPoints; ++i) {
        rings.push_back(createRing(section,
            trajectory[i * 3], trajectory[i * 3 + 1], trajectory[i * 3 + 2],
            scaling[i]));
    }

    // Триангуляция
    triangulateRings(rings);
    computeCenter();

    // Вычисляем диапазон по Z для текстурных координат
    if (!meshFaces.empty()) {
        double minZ = meshFaces[0].verts[0].z;
        double maxZ = minZ;
        for (const auto& f : meshFaces)
            for (const auto& v : f.verts) {
                if (v.z < minZ) minZ = v.z;
                if (v.z > maxZ) maxZ = v.z;
            }
        heightRange = (maxZ - minZ) > 1e-3 ? (maxZ - minZ) : 1.0;
    }

    // Вычисляем нормали после построения сетки
    computeNormals();
}

// ===================================================================================
// КОНФИГУРАЦИЯ И УПРАВЛЕНИЕ КАМЕРОЙ
// ===================================================================================

// Пути к файлам с данными и текстурами
struct AppConfig {
    std::string sectionPath = "files/section.txt";
    std::string trajectoryPath = "files/trajectory.txt";
    std::string scalingPath = "files/scaling.txt";
    std::vector<std::string> texturePaths = { "files/1.bmp", "files/2.bmp" };
};

// Преобразование градусов в радианы
double degToRad(double deg) {
    return deg * M_PI / 180.0;
}

// Контроллер камеры: вращение, зум, ортогональные виды
class CameraController {
private:
    Vector3 angles = { 45, 45, 10 }; // Углы: азимут, высота, расстояние
    double minDistance = 1.0;
    double maxDistance = 100.0;
    bool isOrthoView = false; // Включён ли ортогональный вид
    int orthoMode = 0;        // 0 — спереди, 1 — сбоку, 2 — сверху

public:
    // Вращение камеры (мышка или клавиши)
    void rotate(double deltaX, double deltaY) {
        if (isOrthoView) return; // В орто-режиме вращение запрещено
        angles.x = fmod(angles.x + deltaX, 360.0);
        if (angles.x < 0) angles.x += 360.0;
        angles.y = std::max(-89.0, std::min(89.0, angles.y + deltaY));
    }

    // Зум (колёсико мыши)
    void zoom(double delta) {
        if (isOrthoView) return;
        angles.z = std::max(minDistance, std::min(maxDistance, angles.z + delta));
    }

    // Сброс камеры
    void reset() {
        angles = { 45, 45, 10 };
        isOrthoView = false;
        orthoMode = 0;
    }

    // Включение ортогонального вида
    void setOrthoView(int mode) {
        isOrthoView = true;
        orthoMode = mode;
    }

    // Выход из ортогонального вида
    void exitOrthoView() {
        isOrthoView = false;
    }

    // Настройка матрицы вида (через gluLookAt)
    void setupView() const {
        if (isOrthoView) {
            switch (orthoMode) {
            case 0: gluLookAt(0, 0, 20, 0, 0, 0, 0, 1, 0); break; // Спереди
            case 1: gluLookAt(20, 0, 0, 0, 0, 0, 0, 1, 0); break; // Сбоку
            case 2: gluLookAt(0, 20, 0, 0, 0, 0, 0, 0, -1); break; // Сверху
            }
        }
        else {
            // Расчёт позиции камеры в сферических координатах
            double cameraX = angles.z * cos(degToRad(angles.y)) * cos(degToRad(angles.x));
            double cameraY = angles.z * sin(degToRad(angles.y));
            double cameraZ = angles.z * cos(degToRad(angles.y)) * sin(degToRad(angles.x));
            gluLookAt(cameraX, cameraY, cameraZ, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
        }
    }

    bool getIsOrthoView() const { return isOrthoView; }
};

// ===================================================================================
// СЦЕНА: ГЛОБАЛЬНОЕ СОСТОЯНИЕ ПРИЛОЖЕНИЯ
// ===================================================================================

class Scene {
public:
    Geometry geometry;                     // 3D-модель
    std::vector<unsigned int> textureIDs;  // Идентификаторы загруженных текстур OpenGL
    int textureIndex = -1;                 // Индекс текущей текстуры (-1 = без текстуры)
    bool showWireframe = false;            // Режим каркаса
    bool showNormals = false;              // Отображение нормалей
    bool smoothShading = false;            // Сглаженное освещение
    bool usePerspective = true;            // Перспективная проекция (иначе — ортогональная)
    bool showGrid = true;                  // Сетка на плоскости Z=0
    unsigned short winWidth = 800, winHeight = 600; // Размер окна
    CameraController camera;               // Камера
    AppConfig config;                      // Конфигурация путей

    // Флаги компонентов освещения
    bool showGlobalAmbient = false;
    bool showLightAmbient = false;
    bool showDiffuse = false;
    bool showSpecular = false;
    bool light0Enabled = true;
    bool light1Enabled = false;
    bool light2Enabled = false;

    void reset();                          // Перезагрузка сцены из файлов
    void printSceneInfo() const;           // Вывод информации в консоль
    void updateWindowTitle() const;        // Обновление заголовка окна
};

// Глобальный экземпляр сцены
Scene gScene;

// ===================================================================================
// РАБОТА С ТЕКСТУРАМИ И ОСВЕЩЕНИЕМ
// ===================================================================================

// Генерация текстурных координат (u,v) по вершине
void applyTextureCoordinates(const Vector3& vertex, const Vector3& center, double heightRange) {
    double dx = vertex.x - center.x;
    double dy = vertex.y - center.y;
    double dz = vertex.z - center.z;
    // Угол вокруг оси Z → u (от 0 до 1)
    double angle = std::atan2(dy, dx);
    double u = 0.5 + 0.5 * angle / M_PI;
    // Высота по Z → v (от 0 до 1)
    double v = heightRange > 1e-3 ? dz / heightRange : 0.0;
    glTexCoord2f(static_cast<float>(u), static_cast<float>(v));
}

// Загрузка BMP-текстуры (без использования сторонних библиотек)
void initializeTextureFromBmp(const std::string& filename) {
    unsigned char header[54];
    FILE* file = fopen(filename.c_str(), "rb");
    if (!file) {
        std::cerr << "Текстура не найдена: " << filename << "\n";
        return;
    }
    if (fread(header, 1, 54, file) != 54) {
        fclose(file);
        return;
    }

    // Извлечение параметров из заголовка BMP
    int dataPos = *(int*)&(header[0x0A]);
    int width = *(int*)&(header[0x12]);
    int height = *(int*)&(header[0x16]);
    int rowSize = (width * 3 + 3) & ~3; // Выравнивание строки до 4 байт
    unsigned char* data = new unsigned char[rowSize * height];

    // Чтение пикселей (снизу вверх, как в BMP)
    fseek(file, dataPos, SEEK_SET);
    for (int i = 0; i < height; ++i)
        fread(data + (height - 1 - i) * rowSize, 1, rowSize, file);
    fclose(file);

    // Загрузка в OpenGL
    unsigned int tex;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_BGR_EXT, GL_UNSIGNED_BYTE, data);
    gScene.textureIDs.push_back(tex);
    delete[] data;
    std::cout << "Текстура: " << filename << " (" << width << "x" << height << ")\n";
}

// Настройка освещения в OpenGL
void updateLighting() {
    // Глобальное фоновое освещение
    GLfloat global_amb[] = { 0.0f, 0.0f, 0.0f, 1.0f };
    if (gScene.showGlobalAmbient) {
        global_amb[0] = 0.2f; global_amb[1] = 0.2f; global_amb[2] = 0.25f;
    }
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_amb);

    // Источник света 0 (основной)
    GLfloat amb0[] = { 0.0f, 0.0f, 0.0f, 1.0f };
    GLfloat diff0[] = { 0.0f, 0.0f, 0.0f, 1.0f };
    GLfloat spec0[] = { 0.0f, 0.0f, 0.0f, 1.0f };
    if (gScene.showLightAmbient) { amb0[0] = 0.3f; amb0[1] = 0.3f; amb0[2] = 0.4f; }
    if (gScene.showDiffuse) { diff0[0] = 0.9f; diff0[1] = 0.9f; diff0[2] = 1.0f; }
    if (gScene.showSpecular) { spec0[0] = 1.0f; spec0[1] = 1.0f; spec0[2] = 1.0f; }
    glLightfv(GL_LIGHT0, GL_AMBIENT, amb0);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diff0);
    glLightfv(GL_LIGHT0, GL_SPECULAR, spec0);
    GLfloat pos0[] = { 2.0f, 2.0f, 2.0f, 1.0f }; // Позиция — точечный источник
    glLightfv(GL_LIGHT0, GL_POSITION, pos0);

    // Источник света 1 (дополнительный, зелёный)
    GLfloat pos1[] = { -3.0f, 3.0f, 3.0f, 1.0f };
    glLightfv(GL_LIGHT1, GL_POSITION, pos1);
    GLfloat diff1[] = { 0.0f, 0.7f, 0.0f, 1.0f };
    glLightfv(GL_LIGHT1, GL_DIFFUSE, diff1);

    // Источник света 2 (направленный, синий, "сверху")
    GLfloat pos2[] = { 0.0f, 5.0f, 0.0f, 0.0f }; // w=0 → направленный свет
    glLightfv(GL_LIGHT2, GL_POSITION, pos2);
    GLfloat diff2[] = { 0.0f, 0.0f, 0.8f, 1.0f };
    glLightfv(GL_LIGHT2, GL_DIFFUSE, diff2);

    // Включение/выключение источников
    glEnable(GL_LIGHT0);
    if (gScene.light1Enabled) glEnable(GL_LIGHT1); else glDisable(GL_LIGHT1);
    if (gScene.light2Enabled) glEnable(GL_LIGHT2); else glDisable(GL_LIGHT2);

    // Настройка материала объекта
    GLfloat mat_amb[] = { 0.2f, 0.2f, 0.3f, 1.0f };
    GLfloat mat_diff[] = { 0.8f, 0.8f, 0.9f, 1.0f };
    GLfloat mat_spec[] = { 1.0f, 1.0f, 1.0f, 1.0f };
    GLfloat shininess[] = { 80.0f }; // Блеск (specular highlight)
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_amb);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diff);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_spec);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shininess);
}

// Перестройка проекционной матрицы при изменении размера окна или режима проекции
void rebuildProjectionMatrix() {
    double ratio = (double)gScene.winWidth / (double)gScene.winHeight;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if (gScene.usePerspective && !gScene.camera.getIsOrthoView())
        gluPerspective(60.0, ratio, 0.1, 100.0); // Перспектива
    else
        glOrtho(-ratio * 15, ratio * 15, -15.0, 15.0, -100.0, 100.0); // Орто
    glMatrixMode(GL_MODELVIEW);
}

// ===================================================================================
// УПРАВЛЕНИЕ СЦЕНОЙ (ОБРАБОТЧИКИ КОМАНД)
// ===================================================================================

// Обновление заголовка окна с текущим состоянием
void Scene::updateWindowTitle() const {
    char title[256];
    snprintf(title, sizeof(title),
        "Проекция: %s | Освещение: GA=%d LA=%d D=%d S=%d | Текстура: %s | Нормали: %s",
        (usePerspective && !camera.getIsOrthoView()) ? "Persp" : "Ortho",
        showGlobalAmbient, showLightAmbient, showDiffuse, showSpecular,
        (textureIndex == -1) ? "Нет" : (textureIndex == 0 ? "1" : "2"),
        smoothShading ? "Сглаженные" : "Плоские"
    );
    glutSetWindowTitle(title);
}

// Переключатели состояний (вызываются с клавиатуры или меню)
void switchGlobalAmbient(bool val) { gScene.showGlobalAmbient = val; updateLighting(); gScene.updateWindowTitle(); glutPostRedisplay(); }
void switchLightAmbient(bool val) { gScene.showLightAmbient = val; updateLighting(); gScene.updateWindowTitle(); glutPostRedisplay(); }
void switchDiffuse(bool val) { gScene.showDiffuse = val; updateLighting(); gScene.updateWindowTitle(); glutPostRedisplay(); }
void switchSpecular(bool val) { gScene.showSpecular = val; updateLighting(); gScene.updateWindowTitle(); glutPostRedisplay(); }
void switchGrid(bool val) { gScene.showGrid = val; glutPostRedisplay(); }
void switchCarcass(bool val) { gScene.showWireframe = val; glutPostRedisplay(); }
void switchNormals(bool val) { gScene.showNormals = val; glutPostRedisplay(); }
void switchSmooth(bool val) { gScene.smoothShading = val; gScene.updateWindowTitle(); glutPostRedisplay(); }

// Переключение текстур по кругу
void switchTexture(int val) {
    gScene.textureIndex = val;
    if (gScene.textureIndex >= (int)gScene.textureIDs.size()) gScene.textureIndex = -1;
    gScene.updateWindowTitle();
    glutPostRedisplay();
}

// Переключение проекции (только если не в орто-режиме)
void switchPerspective(bool val) {
    if (gScene.camera.getIsOrthoView()) return;
    gScene.usePerspective = val;
    rebuildProjectionMatrix();
    gScene.updateWindowTitle();
    glutPostRedisplay();
}

// Установка ортогонального вида
void setOrthoView(int mode) {
    gScene.camera.setOrthoView(mode);
    gScene.usePerspective = false;
    rebuildProjectionMatrix();
    gScene.updateWindowTitle();
    glutPostRedisplay();
}

// Выход из ортогонального вида
void exitOrthoView() {
    gScene.camera.exitOrthoView();
    rebuildProjectionMatrix();
    gScene.updateWindowTitle();
    glutPostRedisplay();
}

// Включение/выключение дополнительных источников света
void toggleLight1() { gScene.light1Enabled = !gScene.light1Enabled; updateLighting(); glutPostRedisplay(); }
void toggleLight2() { gScene.light2Enabled = !gScene.light2Enabled; updateLighting(); glutPostRedisplay(); }

// Изменение положения камеры
void changeCamera(double dx, double dy, double dz) {
    gScene.camera.rotate(dx, dy);
    if (dz != 0) gScene.camera.zoom(dz);
    rebuildProjectionMatrix();
    glutPostRedisplay();
}

// ===================================================================================
// ОТРИСОВКА СЦЕНЫ
// ===================================================================================

// Рисование сетки на плоскости Z=0 и осей координат
void drawGrid() {
    if (!gScene.showGrid) return;
    glColor3f(0.3f, 0.3f, 0.3f);
    glBegin(GL_LINES);
    for (int i = -20; i <= 20; ++i) {
        glVertex3f(-20, 0, i); glVertex3f(20, 0, i);
        glVertex3f(i, 0, -20); glVertex3f(i, 0, 20);
    }
    glEnd();
    glLineWidth(2);
    glBegin(GL_LINES);
    glColor3f(1, 0.3, 0.3); glVertex3f(0, 0, 0); glVertex3f(5, 0, 0); // X — красная
    glColor3f(0.3, 1, 0.3); glVertex3f(0, 0, 0); glVertex3f(0, 5, 0); // Y — зелёная
    glColor3f(0.3, 0.3, 1); glVertex3f(0, 0, 0); glVertex3f(0, 0, 5); // Z — синяя
    glEnd();
    glLineWidth(1);
}

// Основная отрисовка 3D-модели
void drawFigure() {
    const auto& faces = gScene.geometry.meshFaces;
    glEnable(GL_LIGHTING);
    if (gScene.textureIndex != -1) {
        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, gScene.textureIDs[gScene.textureIndex]);
    }
    glBegin(GL_TRIANGLES);
    for (size_t i = 0; i < faces.size(); ++i) {
        const auto& f = faces[i];
        for (int j = 0; j < 3; ++j) {
            const Vector3& p = f.verts[j];
            // Выбор нормали: сглаженная или плоская
            Vector3 n = gScene.smoothShading ? gScene.geometry.smoothNormals[i * 3 + j] : gScene.geometry.flatNormals[i];
            applyTextureCoordinates(p, gScene.geometry.center, gScene.geometry.heightRange);
            glNormal3d(n.x, n.y, n.z); // Указание нормали для освещения
            glVertex3d(p.x, p.y, p.z);  // Указание вершины
        }
    }
    glEnd();
    if (gScene.textureIndex != -1) glDisable(GL_TEXTURE_2D);
    glDisable(GL_LIGHTING);
}

// Рисование каркаса (рёбер модели)
void drawCarcass() {
    glColor3ub(255, 0, 0);
    glBegin(GL_LINES);
    for (const auto& f : gScene.geometry.meshFaces)
        for (int j = 0; j < 3; ++j) {
            const Vector3& a = f.verts[j];
            const Vector3& b = f.verts[(j + 1) % 3];
            glVertex3d(a.x, a.y, a.z);
            glVertex3d(b.x, b.y, b.z);
        }
    glEnd();
}

// Рисование нормалей как коротких зелёных линий
void drawNormals() {
    glColor3ub(0, 255, 0);
    glBegin(GL_LINES);
    const auto& faces = gScene.geometry.meshFaces;
    for (size_t i = 0; i < faces.size(); ++i) {
        const auto& f = faces[i];
        for (int j = 0; j < 3; ++j) {
            const Vector3& p = f.verts[j];
            Vector3 n = gScene.smoothShading ? gScene.geometry.smoothNormals[i * 3 + j] : gScene.geometry.flatNormals[i];
            glVertex3d(p.x, p.y, p.z);
            glVertex3d(p.x + n.x * 0.1, p.y + n.y * 0.1, p.z + n.z * 0.1);
        }
    }
    glEnd();
}

// ===================================================================================
// GLUT CALLBACKS
// ===================================================================================

// Главная функция отрисовки
void Display() {
    glClearColor(0.00f, 0.00f, 0.00f, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // Обновляем освещение (позиции источников задаются в мировых координатах)
    updateLighting();

    // Устанавливаем камеру
    gScene.camera.setupView();

    // Рисуем элементы сцены
    drawGrid();
    if (gScene.showWireframe)
        drawCarcass();
    else
        drawFigure();
    if (gScene.showNormals)
        drawNormals();

    glutSwapBuffers(); // Двойная буферизация
}

// Изменение размера окна
void Reshape(int w, int h) {
    gScene.winWidth = w;
    gScene.winHeight = h;
    glViewport(0, 0, w, h);
    rebuildProjectionMatrix();
}

// Обработка нажатий клавиш
void Keyboard(unsigned char key, int x, int y) {
    switch (key) {
    case 27: // ESC
        if (gScene.camera.getIsOrthoView()) exitOrthoView();
        else exit(0);
        break;
    case ' ': // Пробел — сброс
        gScene.reset();
        gScene.camera.reset();
        exitOrthoView();
        break;
    case '1': switchGlobalAmbient(!gScene.showGlobalAmbient); break;
    case '2': switchLightAmbient(!gScene.showLightAmbient); break;
    case '3': switchDiffuse(!gScene.showDiffuse); break;
    case '4': switchSpecular(!gScene.showSpecular); break;
    case '5': toggleLight1(); break;
    case '6': toggleLight2(); break;
    case 'g': switchGrid(!gScene.showGrid); break;
    case 'r': switchCarcass(!gScene.showWireframe); break;
    case 't': switchNormals(!gScene.showNormals); break;
    case 'e': switchSmooth(!gScene.smoothShading); break;
    case 'y': switchTexture(gScene.textureIndex + 1); break;
    case 'q': switchPerspective(!gScene.usePerspective); break;
    case 'f': setOrthoView(0); break; // Спереди
    case 'b': setOrthoView(1); break; // Сбоку
    case 'u': setOrthoView(2); break; // Сверху
    case 'i': gScene.printSceneInfo(); break;
        // Управление камерой (WASD + ZX для зума)
    case 'a': changeCamera(5, 0, 0); break;
    case 'd': changeCamera(-5, 0, 0); break;
    case 'w': changeCamera(0, 5, 0); break;
    case 's': changeCamera(0, -5, 0); break;
    case 'z': changeCamera(0, 0, -0.5); break;
    case 'x': changeCamera(0, 0, 0.5); break;
    }
}

// Обработка колёсика мыши (только в перспективе)
void Mouse(int button, int state, int x, int y) {
    if (state == GLUT_UP) return;
    if (gScene.camera.getIsOrthoView()) return;
    if (button == 3) changeCamera(0, 0, -1); // Колёсико вниз
    else if (button == 4) changeCamera(0, 0, 1); // Колёсико вверх
}

// Обработка меню (правая кнопка мыши)
void Menu(int opt) {
    // ... (переключение опций через меню)
    // Каждый case вызывает соответствующую функцию управления
    switch (opt) {
    case 0: switchGlobalAmbient(!gScene.showGlobalAmbient); break;
    case 1: switchLightAmbient(!gScene.showLightAmbient); break;
    case 2: switchDiffuse(!gScene.showDiffuse); break;
    case 3: switchSpecular(!gScene.showSpecular); break;
    case 4: switchGrid(true); break;
    case 5: switchGrid(false); break;
    case 6: switchCarcass(true); break;
    case 7: switchCarcass(false); break;
    case 8: switchNormals(true); break;
    case 9: switchNormals(false); break;
    case 10: switchSmooth(true); break;
    case 11: switchSmooth(false); break;
    case 12: switchTexture(-1); break;
    case 13: switchTexture(0); break;
    case 14: switchTexture(1); break;
    case 15: switchPerspective(true); break;
    case 16: switchPerspective(false); break;
    case 17: setOrthoView(0); break;
    case 18: setOrthoView(1); break;
    case 19: setOrthoView(2); break;
    case 20: gScene.reset(); gScene.camera.reset(); exitOrthoView(); break;
    case 21: gScene.printSceneInfo(); break;
    case 22: toggleLight1(); break;
    case 23: toggleLight2(); break;
    }
}

// ===================================================================================
// СЦЕНА: РЕАЛИЗАЦИЯ МЕТОДОВ
// ===================================================================================

// Перезагрузка сцены из файлов
void Scene::reset() {
    std::vector<double> section, trajectory, scaling;
    if (!FileManager::loadData(gScene.config.sectionPath, section, "even") ||
        !FileManager::loadData(gScene.config.trajectoryPath, trajectory, "triples") ||
        !FileManager::loadData(gScene.config.scalingPath, scaling)) {
        throw std::runtime_error("Ошибка загрузки файлов");
    }

    geometry.buildFromData(section, trajectory, scaling);
    printSceneInfo();
    updateWindowTitle();
}

// Вывод информации о сцене в консоль
void Scene::printSceneInfo() const {
    std::cout << "=== Сцена ===\n";
    std::cout << "Полигоны: " << geometry.meshFaces.size() << "\n";
    std::cout << "Центр: (" << geometry.center.x << ", " << geometry.center.y << ", " << geometry.center.z << ")\n";
    std::cout << "==============\n";
}

// ===================================================================================
// ИНИЦИАЛИЗАЦИЯ И ЗАПУСК
// ===================================================================================

// Инициализация GLUT и OpenGL
bool initializeApplication(int argc, char* argv[]) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize(gScene.winWidth, gScene.winHeight);
    glutCreateWindow("ЛР3: Тиражирование сечений");

    // Регистрация callback-функций
    glutDisplayFunc(Display);
    glutReshapeFunc(Reshape);
    glutKeyboardFunc(Keyboard);
    glutMouseFunc(Mouse);

    // Включение необходимых возможностей OpenGL
    glEnable(GL_DEPTH_TEST);          // Z-буфер
    glEnable(GL_COLOR_MATERIAL);      // Цвет влияет на материал
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_NORMALIZE);           // Автонормализация нормалей
    glEnable(GL_LIGHTING);            // Освещение

    gScene.updateWindowTitle();

    // Загрузка геометрии
    try {
        gScene.reset();
    }
    catch (const std::exception& e) {
        std::cerr << "Ошибка: " << e.what() << "\n";
        return false;
    }

    // Загрузка текстур
    for (const auto& p : gScene.config.texturePaths)
        initializeTextureFromBmp(p);

    rebuildProjectionMatrix();
    return true;
}

// Точка входа
int main(int argc, char* argv[]) {
    setlocale(LC_ALL, "Russian");
    std::cout << "Управление:\n";
    std::cout << " 1-4: GlobalAmb / LightAmb / Diffuse / Specular\n";
    std::cout << " G: сетка, R: каркас, T: нормали, E: сглаживание\n";
    std::cout << " Y: текстура, Q: проекция\n";
    std::cout << " F/B/U: вид спереди / сбоку / сверху\n";
    std::cout << " WASD/ZX: камера, Колёсико: зум\n";
    std::cout << " Пробел: сброс, ESC: выход\n";

    if (!initializeApplication(argc, argv)) {
        std::cerr << "Ошибка инициализации\n";
        return 1;
    }

    // Создание иерархического меню (правая кнопка мыши)
    int menuLight = glutCreateMenu(Menu);
    glutAddMenuEntry("Global Ambient", 0);
    glutAddMenuEntry("Light Ambient", 1);
    glutAddMenuEntry("Diffuse", 2);
    glutAddMenuEntry("Specular", 3);

    int menuGrid = glutCreateMenu(Menu);
    glutAddMenuEntry("Включить", 4);
    glutAddMenuEntry("Выключить", 5);

    int menuCarcass = glutCreateMenu(Menu);
    glutAddMenuEntry("Включить", 6);
    glutAddMenuEntry("Выключить", 7);

    int menuNormals = glutCreateMenu(Menu);
    glutAddMenuEntry("Включить", 8);
    glutAddMenuEntry("Выключить", 9);

    int menuSmooth = glutCreateMenu(Menu);
    glutAddMenuEntry("Включить", 10);
    glutAddMenuEntry("Выключить", 11);

    int menuTex = glutCreateMenu(Menu);
    glutAddMenuEntry("Нет", 12);
    glutAddMenuEntry("1", 13);
    glutAddMenuEntry("2", 14);

    int menuProj = glutCreateMenu(Menu);
    glutAddMenuEntry("Перспектива", 15);
    glutAddMenuEntry("Орто", 16);

    int menuViews = glutCreateMenu(Menu);
    glutAddMenuEntry("Спереди", 17);
    glutAddMenuEntry("Сбоку", 18);
    glutAddMenuEntry("Сверху", 19);

    int menuLights = glutCreateMenu(Menu);
    glutAddMenuEntry("Источник 1 (вкл/выкл)", 22);
    glutAddMenuEntry("Источник 2 (вкл/выкл)", 23);

    int menuSys = glutCreateMenu(Menu);
    glutAddMenuEntry("Сброс", 20);
    glutAddMenuEntry("Инфо", 21);

    glutCreateMenu(Menu);
    glutAddSubMenu("Освещение", menuLight);
    glutAddSubMenu("Сетка", menuGrid);
    glutAddSubMenu("Каркас", menuCarcass);
    glutAddSubMenu("Нормали", menuNormals);
    glutAddSubMenu("Сглаживание", menuSmooth);
    glutAddSubMenu("Текстура", menuTex);
    glutAddSubMenu("Проекция", menuProj);
    glutAddSubMenu("Виды", menuViews);
    glutAddSubMenu("Источники", menuLights);
    glutAddSubMenu("Система", menuSys);
    glutAttachMenu(GLUT_RIGHT_BUTTON);

    std::cout << "Запущено.\n";
    glutMainLoop(); 
    return 0;
}
