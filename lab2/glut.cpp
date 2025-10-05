#include "stdafx.h"
#define _USE_MATH_DEFINES
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdio>
#include <algorithm>

// OpenGL и GLUT заголовки
// GLUT используется для создания окна и обработки событий (клавиатура, мышь и т.д.)
// GL и GLU — стандартные библиотеки для рендеринга 3D-графики
#include <E:\git_work\Computer_Graphics_Course\lab3\GL\glut.h>
#include <GL/glu.h>
#include <GL/gl.h>

// ===================================================================================
// КОНФИГУРАЦИЯ ПРИЛОЖЕНИЯ
// ===================================================================================
// Структура AppConfig хранит пути к файлам с данными:
// - section.txt: 2D-сечение (профиль) объекта в виде пар (x, y)
// - trajectory.txt: 3D-траектория — точки (x, y, z), по которым "тянется" профиль
// - scaling.txt: коэффициенты масштабирования для каждой точки траектории
// - texturePaths: пути к BMP-текстурам для наложения на модель
struct AppConfig {
    std::string sectionPath = "files/section.txt";
    std::string trajectoryPath = "files/trajectory.txt";
    std::string scalingPath = "files/scaling.txt";
    std::vector<std::string> texturePaths = { "files/1.bmp", "files/2.bmp" };
};

// Преобразование градусов в радианы (для тригонометрических функций OpenGL)
double degToRad(double deg) {
    return deg * M_PI / 180.0;
}

// ===================================================================================
// МАТЕМАТИЧЕСКИЕ СТРУКТУРЫ И ФУНКЦИИ
// ===================================================================================

// Трёхмерный вектор с операторами сравнения и равенства
// Используется для представления точек и нормалей
struct Vector3 {
    double x, y, z;

    // Оператор < необходим для использования Vector3 в std::map
    bool operator<(const Vector3& o) const {
        if (x != o.x) return x < o.x;
        if (y != o.y) return y < o.y;
        return z < o.z;
    }

    // Оператор == с учётом погрешности (epsilon) для сравнения вещественных чисел
    bool operator==(const Vector3& o) const {
        const double epsilon = 1e-9;
        return std::abs(x - o.x) < epsilon &&
            std::abs(y - o.y) < epsilon &&
            std::abs(z - o.z) < epsilon;
    }
};

// Хеш-функция для Vector3 — необходима для использования в unordered_map
struct Vector3Hash {
    std::size_t operator()(const Vector3& v) const {
        // Простой, но рабочий способ комбинирования хешей
        return std::hash<double>()(v.x) ^
            (std::hash<double>()(v.y) << 1) ^
            (std::hash<double>()(v.z) << 2);
    }
};

// Функтор сравнения для unordered_map с учётом погрешности
struct Vector3Equal {
    bool operator()(const Vector3& a, const Vector3& b) const {
        const double epsilon = 1e-9;
        return std::abs(a.x - b.x) < epsilon &&
            std::abs(a.y - b.y) < epsilon &&
            std::abs(a.z - b.z) < epsilon;
    }
};

// Грань (треугольник) модели — состоит из трёх вершин
struct Face {
    Vector3 verts[3];
};

// ===================================================================================
// УПРАВЛЕНИЕ КАМЕРОЙ
// ===================================================================================
// Камера реализована как "орбитальная": вращается вокруг центра сцены (0,0,0)
// Использует сферические координаты: углы поворота и расстояние до центра
class CameraController {
private:
    // angles.x — азимут (вокруг оси Y), angles.y — угол возвышения (от -89° до +89°),
    // angles.z — расстояние от камеры до центра сцены
    Vector3 angles = { 45, 45, 10 };
    double minDistance = 1.0;  // Минимальное расстояние (ограничение приближения)
    double maxDistance = 100.0; // Максимальное расстояние (ограничение отдаления)

public:
    // Вращение камеры: deltaX — горизонтальное, deltaY — вертикальное
    void rotate(double deltaX, double deltaY) {
        angles.x = fmod(angles.x + deltaX, 360.0);
        if (angles.x < 0) angles.x += 360.0;

        // Ограничение угла возвышения, чтобы избежать "переворота" камеры
        angles.y = std::max(-89.0, std::min(89.0, angles.y + deltaY));
    }

    // Изменение расстояния до центра (зум)
    void zoom(double delta) {
        angles.z = std::max(minDistance, std::min(maxDistance, angles.z + delta));
    }

    // Сброс камеры в начальное положение
    void reset() {
        angles = { 45, 45, 10 };
    }

    // Настройка матрицы вида с помощью gluLookAt
    // Вычисляет позицию камеры в декартовых координатах из сферических
    void setupView() const {
        double cameraX = angles.z * cos(degToRad(angles.y)) * cos(degToRad(angles.x));
        double cameraY = angles.z * sin(degToRad(angles.y));
        double cameraZ = angles.z * cos(degToRad(angles.y)) * sin(degToRad(angles.x));

        // gluLookAt(позиция_камеры, точка_наблюдения, вектор_вверх)
        gluLookAt(cameraX, cameraY, cameraZ, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    }

    // Геттер и сеттер для текущих углов камеры
    const Vector3& getAngles() const { return angles; }
    void setAngles(const Vector3& newAngles) { angles = newAngles; }
};

// ===================================================================================
// ГЕОМЕТРИЯ ОБЪЕКТА
// ===================================================================================
// Класс Geometry отвечает за построение 3D-модели из 2D-сечения, траектории и масштабов
class Geometry {
public:
    std::vector<Face> meshFaces; // Список треугольников модели
    // Карта: вершина -> усреднённая нормаль (для сглаженного освещения)
    std::unordered_map<Vector3, Vector3, Vector3Hash, Vector3Equal> vertexNormals;
    Vector3 center; // Центр масс модели

    // Вычисление центра модели (среднее арифметическое всех вершин)
    void computeCenter();

    // Вычисление нормалей:
    // - сначала нормали граней,
    // - затем усреднение по вершинам для smooth shading
    void computeNormals();

    // Основной метод построения модели:
    // section — 2D-профиль (x0,y0,x1,y1,...),
    // trajectory — 3D-точки (x0,y0,z0,x1,y1,z1,...),
    // scaling — коэффициенты масштаба для каждой точки траектории
    void buildFromData(const std::vector<double>& section,
        const std::vector<double>& trajectory,
        const std::vector<double>& scaling);

private:
    // Создаёт одно "кольцо" вершин вокруг точки траектории
    // Сечение центрируется, масштабируется и смещается вдоль траектории
    std::vector<Vector3> createRing(const std::vector<double>& section,
        const std::vector<double>& trajectory,
        const std::vector<double>& scaling,
        int ringIndex);

    // Соединяет соседние кольца треугольниками (триангуляция поверхности)
    void triangulateRings(const std::vector<std::vector<Vector3>>& rings);
};

// ===================================================================================
// МЕНЕДЖЕР ФАЙЛОВ
// ===================================================================================
// Загружает числовые данные из текстового файла в вектор double
class FileManager {
public:
    static bool loadData(const std::string& path, std::vector<double>& out,
        const std::string& expectedFormat = "") {
        FILE* f = std::fopen(path.c_str(), "r");
        if (!f) {
            std::cerr << "Ошибка: не удалось открыть файл " << path << std::endl;
            return false;
        }

        double value;
        int count = 0;
        while (std::fscanf(f, "%lf", &value) == 1) {
            out.push_back(value);
            count++;
        }
        std::fclose(f);

        // Проверка формата данных:
        // - "even": чётное число (для 2D-сечения: x,y,x,y,...)
        // - "triples": кратно трём (для 3D-траектории: x,y,z,x,y,z,...)
        if (expectedFormat == "even" && count % 2 != 0) {
            std::cerr << "Ошибка: " << path << " должен содержать четное число координат" << std::endl;
            return false;
        }
        if (expectedFormat == "triples" && count % 3 != 0) {
            std::cerr << "Ошибка: " << path << " должен содержать координаты в формате x,y,z" << std::endl;
            return false;
        }

        std::cout << "Загружено " << count << " значений из " << path << std::endl;
        return count > 0;
    }
};

// ===================================================================================
// СЦЕНА
// ===================================================================================
// Основной класс сцены: содержит геометрию, текстуры, настройки отображения и камеру
class Scene {
public:
    Geometry geometry;                    // 3D-модель
    std::vector<unsigned int> textureIDs; // Идентификаторы загруженных текстур OpenGL
    int textureIndex = -1;                // Индекс текущей текстуры (-1 = без текстуры)
    bool showWireframe = false;           // Режим каркаса (линии вместо заполненных полигонов)
    bool showNormals = false;             // Отображение нормалей (векторы из вершин)
    bool smoothShading = false;           // Сглаженное освещение (Gouraud shading)
    bool usePerspective = true;           // Перспективная проекция (иначе — ортографическая)
    unsigned short winWidth = 800, winHeight = 600; // Размер окна
    CameraController camera;              // Управление камерой
    AppConfig config;                     // Конфигурация путей к файлам

    // Перезагрузка сцены: загрузка данных и перестроение модели
    void reset();

    // Вывод информации о сцене в консоль
    void printSceneInfo() const;
};

// Глобальный экземпляр сцены — упрощает доступ из callback-функций GLUT
Scene gScene;

// ===================================================================================
// ВСПОМОГАТЕЛЬНЫЕ МАТЕМАТИЧЕСКИЕ ФУНКЦИИ
// ===================================================================================

// Сложение двух векторов
Vector3 addVec(const Vector3& a, const Vector3& b) {
    return Vector3{ a.x + b.x, a.y + b.y, a.z + b.z };
}

// Нормализация вектора (приведение к единичной длине)
Vector3 normalizeVec(const Vector3& a) {
    double len = std::sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
    if (len < 1e-8) return Vector3{ 0, 0, 0 }; // Защита от деления на ноль
    return Vector3{ a.x / len, a.y / len, a.z / len };
}

// Вычисление нормали к треугольнику по трём вершинам (векторное произведение)
Vector3 calcNormal(const Vector3& v1, const Vector3& v2, const Vector3& v3) {
    double ux = v2.x - v1.x, uy = v2.y - v1.y, uz = v2.z - v1.z;
    double vx = v3.x - v1.x, vy = v3.y - v1.y, vz = v3.z - v1.z;
    Vector3 n{
        uy * vz - uz * vy, // X компонента
        uz * vx - ux * vz, // Y компонента
        ux * vy - uy * vx  // Z компонента
    };
    return normalizeVec(n);
}

// Безопасное получение нормали из карты (с возвратом нулевого вектора при отсутствии)
Vector3 getNormal(const std::unordered_map<Vector3, Vector3, Vector3Hash, Vector3Equal>& normals,
    const Vector3& key) {
    auto it = normals.find(key);
    return it != normals.end() ? it->second : Vector3{ 0, 0, 0 };
}

// ===================================================================================
// РЕАЛИЗАЦИЯ МЕТОДОВ КЛАССА Geometry
// ===================================================================================

void Geometry::computeCenter() {
    if (meshFaces.empty()) {
        center = Vector3{ 0, 0, 0 };
        return;
    }

    double sx = 0, sy = 0, sz = 0;
    int cnt = 0;
    for (const auto& f : meshFaces) {
        for (const auto& v : f.verts) {
            sx += v.x; sy += v.y; sz += v.z; ++cnt;
        }
    }
    center = Vector3{ sx / cnt, sy / cnt, sz / cnt };
}

void Geometry::computeNormals() {
    vertexNormals.clear();
    std::unordered_map<Vector3, Vector3, Vector3Hash, Vector3Equal> tempNormals;

    // Для каждой грани вычисляем нормаль и добавляем её ко всем трём вершинам
    for (const auto& face : meshFaces) {
        Vector3 n = calcNormal(face.verts[0], face.verts[1], face.verts[2]);
        for (int i = 0; i < 3; ++i) {
            const Vector3& v = face.verts[i];
            // Если вершина уже есть — суммируем нормали
            tempNormals[v] = addVec(tempNormals[v], n);
        }
    }

    // Нормализуем усреднённые нормали
    for (auto& it : tempNormals) {
        vertexNormals[it.first] = normalizeVec(it.second);
    }
}

// Создание одного кольца вершин вокруг i-й точки траектории
std::vector<Vector3> Geometry::createRing(const std::vector<double>& section,
    const std::vector<double>& trajectory,
    const std::vector<double>& scaling,
    int ringIndex) {
    const int profilePoints = section.size() / 2;
    std::vector<Vector3> ring;
    ring.reserve(profilePoints);

    // Получаем координаты текущей точки траектории и масштаб
    double tx = trajectory[ringIndex * 3 + 0];
    double ty = trajectory[ringIndex * 3 + 1];
    double tz = trajectory[ringIndex * 3 + 2];
    double scale = scaling[ringIndex];

    // Центрируем сечение: вычисляем его центр (чтобы масштабирование было относительно центра)
    double secCx = 0.0, secCy = 0.0;
    for (int k = 0; k < profilePoints; ++k) {
        secCx += section[2 * k];
        secCy += section[2 * k + 1];
    }
    secCx /= profilePoints;
    secCy /= profilePoints;

    // Для каждой точки сечения:
    // 1. Смещаем относительно центра сечения
    // 2. Масштабируем
    // 3. Возвращаем в исходное положение + смещаем вдоль траектории
    for (int k = 0; k < profilePoints; ++k) {
        double dx = section[2 * k] - secCx;
        double dy = section[2 * k + 1] - secCy;
        double x = dx * scale + secCx + tx;
        double y = dy * scale + secCy + ty;
        double z = tz; // Z-координата берётся из траектории
        ring.emplace_back(Vector3{ x, y, z });
    }

    return ring;
}

// Триангуляция: соединение соседних колец четырёхугольниками, разбитыми на два треугольника
void Geometry::triangulateRings(const std::vector<std::vector<Vector3>>& rings) {
    meshFaces.clear();

    for (size_t i = 0; i < rings.size() - 1; ++i) {
        const std::vector<Vector3>& ring0 = rings[i];
        const std::vector<Vector3>& ring1 = rings[i + 1];
        int ringSize = ring0.size();

        for (int j = 0; j < ringSize; ++j) {
            int j1 = (j + 1) % ringSize;
            // Первый треугольник четырёхугольника
            meshFaces.push_back(Face{ ring0[j], ring1[j], ring1[j1] });
            // Второй треугольник
            meshFaces.push_back(Face{ ring0[j], ring1[j1], ring0[j1] });
        }
    }
}

// Основной метод построения модели
void Geometry::buildFromData(const std::vector<double>& section,
    const std::vector<double>& trajectory,
    const std::vector<double>& scaling) {
    // Валидация входных данных
    if (section.size() % 2 != 0 || section.size() < 6) {
        throw std::runtime_error("Неверный формат сечения: требуется минимум 3 точки (6 координат)");
    }

    const int profilePoints = section.size() / 2;
    const int trajectoryPoints = trajectory.size() / 3;

    if (trajectoryPoints < 2) {
        throw std::runtime_error("Траектория должна содержать минимум 2 точки");
    }

    if (trajectoryPoints != (int)scaling.size()) {
        throw std::runtime_error("Количество точек траектории и коэффициентов масштабирования не совпадает");
    }

    // Построение колец для каждой точки траектории
    std::vector<std::vector<Vector3>> rings;
    rings.reserve(trajectoryPoints);

    for (int i = 0; i < trajectoryPoints; ++i) {
        rings.push_back(createRing(section, trajectory, scaling, i));
    }

    // Триангуляция поверхности
    triangulateRings(rings);

    // Вычисление центра и нормалей
    computeCenter();
    computeNormals();
}

// ===================================================================================
// РЕАЛИЗАЦИЯ МЕТОДОВ КЛАССА Scene
// ===================================================================================

void Scene::reset() {
    std::vector<double> section, trajectory, scaling;

    // Загрузка данных из файлов
    if (!FileManager::loadData(config.sectionPath, section, "even") ||
        !FileManager::loadData(config.trajectoryPath, trajectory, "triples") ||
        !FileManager::loadData(config.scalingPath, scaling)) {
        throw std::runtime_error("Ошибка загрузки файлов данных");
    }

    // Построение модели
    geometry.buildFromData(section, trajectory, scaling);
    printSceneInfo();
}

void Scene::printSceneInfo() const {
    std::cout << "=== Информация о сцене ===" << std::endl;
    std::cout << "Количество полигонов: " << geometry.meshFaces.size() << std::endl;
    std::cout << "Количество вершин: " << geometry.meshFaces.size() * 3 << std::endl;
    std::cout << "Центр объекта: (" << geometry.center.x << ", "
        << geometry.center.y << ", " << geometry.center.z << ")" << std::endl;
    std::cout << "Текстур: " << textureIDs.size() << std::endl;
    std::cout << "Размер окна: " << winWidth << "x" << winHeight << std::endl;
    std::cout << "=========================" << std::endl;
}

// ===================================================================================
// РАБОТА С ТЕКСТУРАМИ
// ===================================================================================

// Генерация текстурных координат (u, v) с помощью цилиндрической проекции
// u — угол вокруг оси Z, v — нормализованная высота по Z
void applyTextureCoordinates(const Vector3& vertex, const Vector3& center,
    double heightRange, double angleOffset = 0.0) {
    double dx = vertex.x - center.x;
    double dy = vertex.y - center.y;
    double dz = vertex.z - center.z;

    // Угол в радианах, нормализованный к [0, 1]
    double angle = std::atan2(dy, dx) + angleOffset;
    double u = 0.5 + 0.5 * angle / M_PI;

    // Высота, нормализованная к диапазону heightRange
    double v = (dz - center.z) / heightRange;

    glTexCoord2f(static_cast<float>(u), static_cast<float>(v));
}

// Загрузка BMP-текстуры (без поддержки сжатия, только 24-bit RGB)
void initializeTextureFromBmp(const std::string& filename) {
    unsigned int tex;
    unsigned char bytes[54]; // BMP заголовок — 54 байта
    FILE* file = std::fopen(filename.c_str(), "rb");
    if (!file) {
        std::cerr << "Ошибка: не удалось открыть текстуру " << filename << std::endl;
        return;
    }

    int err = std::fread(bytes, 1, 54, file);
    if (err != 54) {
        std::cerr << "Ошибка чтения BMP заголовка: " << filename << std::endl;
        fclose(file);
        return;
    }

    // Извлечение параметров из заголовка
    int offset;
    memcpy(&offset, bytes + 10, sizeof(int));
    int bmpw = *(int*)(bytes + 18);
    int bmph = *(int*)(bytes + 22);

    // Чтение пикселей (BGR, bottom-up)
    unsigned char* pixels = new unsigned char[3 * bmpw * bmph];
    fseek(file, offset, SEEK_SET);
    err = std::fread(pixels, 3 * bmpw * bmph, 1, file);
    if (err != 1) {
        std::cerr << "Ошибка чтения BMP данных: " << filename << std::endl;
        delete[] pixels;
        fclose(file);
        return;
    }
    std::fclose(file);

    // Создание текстуры OpenGL
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    // Обратите внимание: формат GL_BGR_EXT (т.к. BMP хранит BGR)
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, bmpw, bmph, 0, GL_BGR_EXT, GL_UNSIGNED_BYTE, pixels);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE); // Смешивание с цветом материала
    gScene.textureIDs.push_back(tex);
    delete[] pixels;

    std::cout << "Загружена текстура: " << filename << " (" << bmpw << "x" << bmph << ")" << std::endl;
}

// ===================================================================================
// ОСВЕЩЕНИЕ
// ===================================================================================

// Настройка источника света GL_LIGHT0
void initializeLighting() {
    float ambientLight0[] = { 0.1f, 0.1f, 0.1f, 1.0f };   // Рассеянный свет
    float diffuseLight0[] = { 0.8f, 0.8f, 0.9f, 1.0f };   // Диффузный свет
    float specularLight0[] = { 0.4f, 0.4f, 0.4f, 1.0f };  // Зеркальный свет
    float lightPosition0[] = { 10.0f, 10.0f, 10.0f, 1.0f }; // Позиция источника

    glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight0);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight0);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight0);
    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition0);

    // Настройка материала объекта
    float specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
    float shininess[] = { 50.0f }; // Блеск (0–128)
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shininess);
}

// ===================================================================================
// ПРОЕКЦИЯ И КАМЕРА
// ===================================================================================

// Перестройка матрицы проекции при изменении окна или режима проекции
void rebuildProjectionMatrix() {
    double ratio = (double)gScene.winWidth / (double)gScene.winHeight;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if (gScene.usePerspective)
        gluPerspective(90.0, ratio, 0.1, 100.0); // Угол обзора 90°
    else
        glOrtho(-ratio, ratio, -ratio, ratio, -100.0, 100.0); // Ортографическая проекция

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gScene.camera.setupView(); // Настройка вида камеры
}

// ===================================================================================
// УПРАВЛЕНИЕ СЦЕНОЙ
// ===================================================================================

// Обёртка для изменения параметров камеры и перерисовки
void changeCamera(double x, double y, double z) {
    gScene.camera.rotate(x, y);
    if (z != 0) gScene.camera.zoom(z);

    rebuildProjectionMatrix();
    glutPostRedisplay(); // Запрос на перерисовку
}

// Включение/выключение источника света
void switchLight(GLenum id, bool val) {
    if (val) glEnable(id); else glDisable(id);
    glutPostRedisplay();
}

// Обработка колеса мыши (кнопки 3 и 4 в GLUT)
void Mouse(int button, int state, int x, int y) {
    if (state == 0) { // Нажатие
        if (button == 3) changeCamera(0, 0, -1.0);  // Вверх — приближение
        else if (button == 4) changeCamera(0, 0, 1.0); // Вниз — отдаление
    }
}

// Переключение проекции
void switchPerspective(bool val) {
    gScene.usePerspective = val;
    rebuildProjectionMatrix();
    glutPostRedisplay();
}

// Переключение сглаженного освещения
void switchSmooth(bool val) {
    gScene.smoothShading = val;
    glutPostRedisplay();
}

// Переключение каркасного режима
void switchCarcass(bool val) {
    gScene.showWireframe = val;
    glutPostRedisplay();
}

// Переключение отображения нормалей
void switchShowNormals(bool val) {
    gScene.showNormals = val;
    glutPostRedisplay();
}

// Переключение текстур (циклически: -1 → 0 → 1 → -1 ...)
void switchTexture(int val) {
    gScene.textureIndex = val;
    if (gScene.textureIndex >= (int)gScene.textureIDs.size())
        gScene.textureIndex = -1;
    glutPostRedisplay();
}

// ===================================================================================
// РИСОВАНИЕ
// ===================================================================================

// Рисование сетки координат (серые линии на плоскости Y=0)
void drawCoordinateNetwork() {
    glColor3ub(127, 127, 127);
    glBegin(GL_LINES);
    for (float i = -100; i <= 100; i++) {
        glVertex3f(-100, 0, i);
        glVertex3f(100, 0, i);
        glVertex3f(i, 0, -100);
        glVertex3f(i, 0, 100);
    }
    glVertex3f(0, -100, 0);
    glVertex3f(0, 100, 0);
    glEnd();
}

// Основное рисование модели с освещением и текстурами
void drawFigure() {
    const auto& faces = gScene.geometry.meshFaces;
    const auto& normals = gScene.geometry.vertexNormals;
    bool smooth = gScene.smoothShading;

    glEnable(GL_LIGHTING);
    if (gScene.textureIndex != -1) {
        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, gScene.textureIDs[gScene.textureIndex]);
    }

    glBegin(GL_TRIANGLES);
    for (const auto& face : faces) {
        for (int j = 0; j < 3; j++) {
            const Vector3& p = face.verts[j];
            Vector3 n;
            if (smooth) {
                n = getNormal(normals, p); // Используем усреднённую нормаль
            }
            else {
                n = calcNormal(face.verts[0], face.verts[1], face.verts[2]); // Нормаль грани
            }

            applyTextureCoordinates(p, gScene.geometry.center, 30.0);
            glNormal3d(n.x, n.y, n.z); // Передаём нормаль в OpenGL
            glVertex3d(p.x, p.y, p.z);
        }
    }
    glEnd();

    if (gScene.textureIndex != -1) glDisable(GL_TEXTURE_2D);
    glDisable(GL_LIGHTING);
}

// Рисование каркаса (красные линии)
void drawCarcass() {
    glColor3ub(255, 0, 0);
    const auto& faces = gScene.geometry.meshFaces;
    glBegin(GL_LINES);
    for (const auto& face : faces) {
        for (int j = 0; j < 3; j++) {
            const Vector3& v0 = face.verts[j];
            const Vector3& v1 = face.verts[(j + 1) % 3];
            glVertex3d(v0.x, v0.y, v0.z);
            glVertex3d(v1.x, v1.y, v1.z);
        }
    }
    glEnd();
}

// Рисование нормалей (зелёные векторы длиной 0.1)
void drawNormals() {
    glColor3ub(0, 255, 0);
    const auto& faces = gScene.geometry.meshFaces;
    const auto& normals = gScene.geometry.vertexNormals;
    bool smooth = gScene.smoothShading;

    glBegin(GL_LINES);
    for (const auto& face : faces) {
        for (int j = 0; j < 3; j++) {
            const Vector3& p = face.verts[j];
            Vector3 n;
            if (smooth) {
                n = getNormal(normals, p);
            }
            else {
                n = calcNormal(face.verts[0], face.verts[1], face.verts[2]);
            }

            glVertex3d(p.x, p.y, p.z);
            glVertex3d(p.x + n.x * 0.1, p.y + n.y * 0.1, p.z + n.z * 0.1);
        }
    }
    glEnd();
}

// ===================================================================================
// GLUT CALLBACKS
// ===================================================================================

// Основная функция отрисовки
void Display(void) {
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    drawCoordinateNetwork();
    if (gScene.showWireframe) drawCarcass();
    else drawFigure();

    if (gScene.showNormals) drawNormals();

    glutSwapBuffers(); // Двойная буферизация
}

// Изменение размера окна
void Reshape(GLint newW, GLint newH) {
    gScene.winWidth = newW;
    gScene.winHeight = newH;
    glViewport(0, 0, gScene.winWidth, gScene.winHeight);
    rebuildProjectionMatrix();
}

// Обработка клавиатуры
void Keyboard(unsigned char key, int x, int y) {
    switch (key) {
    case '\033': exit(0); break; // ESC
    case 'a': changeCamera(5, 0, 0); break; // Влево
    case 'd': changeCamera(-5, 0, 0); break; // Вправо
    case 'w': changeCamera(0, 5, 0); break; // Вверх
    case 's': changeCamera(0, -5, 0); break; // Вниз
    case '+': case '=': changeCamera(0, 0, -0.5); break; // Приближение
    case '-': case '_': changeCamera(0, 0, 0.5); break; // Отдаление
    case '1': switchLight(GL_LIGHT0, !glIsEnabled(GL_LIGHT0)); break; // Свет
    case 'q': switchPerspective(!gScene.usePerspective); break; // Проекция
    case 'e': switchSmooth(!gScene.smoothShading); break; // Сглаживание
    case 'r': switchCarcass(!gScene.showWireframe); break; // Каркас
    case 't': switchShowNormals(!gScene.showNormals); break; // Нормали
    case 'y': switchTexture(gScene.textureIndex + 1); break; // Текстура
    case 'i': gScene.printSceneInfo(); break; // Информация
    case ' ': gScene.reset(); break; // Перезагрузка
    }
}

// Обработка меню (правая кнопка мыши)
void Menu(int option) {
    switch (option) {
    case 0: switchPerspective(true); break;
    case 1: switchPerspective(false); break;
    case 2: switchSmooth(true); break;
    case 3: switchSmooth(false); break;
    case 4: switchCarcass(true); break;
    case 5: switchCarcass(false); break;
    case 6: switchShowNormals(true); break;
    case 7: switchShowNormals(false); break;
    case 8: switchTexture(-1); break; // Без текстуры
    case 9: switchTexture(0); break;  // Текстура 1
    case 10: switchTexture(1); break; // Текстура 2
    case 11: switchLight(GL_LIGHT0, true); break;
    case 12: switchLight(GL_LIGHT0, false); break;
    case 13: gScene.reset(); break;
    case 14: gScene.printSceneInfo(); break;
    }
}

// ===================================================================================
// ИНИЦИАЛИЗАЦИЯ ПРИЛОЖЕНИЯ
// ===================================================================================

bool initializeApplication(int argc, char* argv[]) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize(gScene.winWidth, gScene.winHeight);
    glutCreateWindow("Лабораторная работа №3 — Тиражирование сечений");

    // Регистрация callback-функций
    glutDisplayFunc(Display);
    glutReshapeFunc(Reshape);
    glutKeyboardFunc(Keyboard);
    glutMouseFunc(Mouse);

    // Настройка OpenGL
    glEnable(GL_DEPTH_TEST); // Включить Z-буфер
    glEnable(GL_COLOR_MATERIAL); // Цвет материала = цвет вершины
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

    // Включение освещения
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    initializeLighting();

    // Отключение глобального ambient light (чтобы не "засвечивало")
    GLfloat no_ambient[] = { 0.0f, 0.0f, 0.0f, 1.0f };
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, no_ambient);

    // Загрузка данных и построение модели
    try {
        gScene.reset();
    }
    catch (const std::exception& e) {
        std::cerr << "Ошибка инициализации: " << e.what() << std::endl;
        return false;
    }

    // Загрузка текстур
    for (const auto& texPath : gScene.config.texturePaths) {
        initializeTextureFromBmp(texPath);
    }

    rebuildProjectionMatrix();

    return true;
}

// ===================================================================================
// ГЛАВНАЯ ФУНКЦИЯ
// ===================================================================================

int main(int argc, char* argv[]) {
    std::cout << "Лабораторная работа №3 - Тиражирование сечений" << std::endl;
    std::cout << "Управление:" << std::endl;
    std::cout << "  WASD - вращение камеры" << std::endl;
    std::cout << "  +/- - приближение/отдаление" << std::endl;
    std::cout << "  Q - переключение проекции" << std::endl;
    std::cout << "  E - сглаживание нормалей" << std::endl;
    std::cout << "  R - каркасный режим" << std::endl;
    std::cout << "  T - отображение нормалей" << std::endl;
    std::cout << "  Y - переключение текстур" << std::endl;
    std::cout << "  I - информация о сцене" << std::endl;
    std::cout << "  Пробел - перезагрузка сцены" << std::endl;
    std::cout << "  ESC - выход" << std::endl;

    if (!initializeApplication(argc, argv)) {
        std::cerr << "Не удалось инициализировать приложение" << std::endl;
        return 1;
    }

    // Создание контекстного меню (правая кнопка мыши)
    int menuProjection = glutCreateMenu(Menu);
    glutAddMenuEntry("Перспективная", 0);
    glutAddMenuEntry("Ортографическая", 1);

    int menuSmooth = glutCreateMenu(Menu);
    glutAddMenuEntry("Включить", 2);
    glutAddMenuEntry("Выключить", 3);

    int menuCarcass = glutCreateMenu(Menu);
    glutAddMenuEntry("Включить", 4);
    glutAddMenuEntry("Выключить", 5);

    int menuNormal = glutCreateMenu(Menu);
    glutAddMenuEntry("Включить", 6);
    glutAddMenuEntry("Выключить", 7);

    int menuTexture = glutCreateMenu(Menu);
    glutAddMenuEntry("Без текстуры", 8);
    glutAddMenuEntry("Текстура 1", 9);
    glutAddMenuEntry("Текстура 2", 10);

    int menuLight1 = glutCreateMenu(Menu);
    glutAddMenuEntry("Включить", 11);
    glutAddMenuEntry("Выключить", 12);

    int menuSystem = glutCreateMenu(Menu);
    glutAddMenuEntry("Перезагрузить сцену", 13);
    glutAddMenuEntry("Информация о сцене", 14);

    glutCreateMenu(Menu);
    glutAddSubMenu("Проекция", menuProjection);
    glutAddSubMenu("Сглаживание", menuSmooth);
    glutAddSubMenu("Каркас", menuCarcass);
    glutAddSubMenu("Нормали", menuNormal);
    glutAddSubMenu("Текстура", menuTexture);
    glutAddSubMenu("Свет", menuLight1);
    glutAddSubMenu("Система", menuSystem);
    glutAttachMenu(GLUT_RIGHT_BUTTON);

    std::cout << "Приложение успешно запущено" << std::endl;

    glutMainLoop(); // Главный цикл GLUT
    return 0;
}