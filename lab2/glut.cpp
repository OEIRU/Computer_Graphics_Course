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
#include <cstring>
#include <locale.h>

#include <GL/glut.h>
#include <GL/glu.h>
#include <GL/gl.h>

// ===================================================================================
// КОНФИГУРАЦИЯ ПРИЛОЖЕНИЯ
// ===================================================================================
struct AppConfig {
    std::string sectionPath = "files/section.txt";
    std::string trajectoryPath = "files/trajectory.txt";
    std::string scalingPath = "files/scaling.txt";
    std::vector<std::string> texturePaths = { "files/1.bmp", "files/2.bmp" };
};

double degToRad(double deg) {
    return deg * M_PI / 180.0;
}

// ===================================================================================
// МАТЕМАТИЧЕСКИЕ СТРУКТУРЫ И ФУНКЦИИ
// ===================================================================================

struct Vector3 {
    double x, y, z;

    bool operator<(const Vector3& o) const {
        if (x != o.x) return x < o.x;
        if (y != o.y) return y < o.y;
        return z < o.z;
    }

    bool operator==(const Vector3& o) const {
        const double epsilon = 1e-9;
        return std::abs(x - o.x) < epsilon &&
            std::abs(y - o.y) < epsilon &&
            std::abs(z - o.z) < epsilon;
    }
};

struct Vector3Hash {
    std::size_t operator()(const Vector3& v) const {
        return std::hash<double>()(v.x) ^
            (std::hash<double>()(v.y) << 1) ^
            (std::hash<double>()(v.z) << 2);
    }
};

struct Vector3Equal {
    bool operator()(const Vector3& a, const Vector3& b) const {
        const double epsilon = 1e-9;
        return std::abs(a.x - b.x) < epsilon &&
            std::abs(a.y - b.y) < epsilon &&
            std::abs(a.z - b.z) < epsilon;
    }
};

struct Face {
    Vector3 verts[3];
};

// ===================================================================================
// УПРАВЛЕНИЕ КАМЕРОЙ
// ===================================================================================
class CameraController {
private:
    Vector3 angles = { 45, 45, 10 };
    double minDistance = 1.0;
    double maxDistance = 100.0;

public:
    void rotate(double deltaX, double deltaY) {
        angles.x = fmod(angles.x + deltaX, 360.0);
        if (angles.x < 0) angles.x += 360.0;
        angles.y = std::max(-89.0, std::min(89.0, angles.y + deltaY));
    }

    void zoom(double delta) {
        angles.z = std::max(minDistance, std::min(maxDistance, angles.z + delta));
    }

    void reset() {
        angles = { 45, 45, 10 };
    }

    void setupView() const {
        double cameraX = angles.z * cos(degToRad(angles.y)) * cos(degToRad(angles.x));
        double cameraY = angles.z * sin(degToRad(angles.y));
        double cameraZ = angles.z * cos(degToRad(angles.y)) * sin(degToRad(angles.x));
        gluLookAt(cameraX, cameraY, cameraZ, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    }

    const Vector3& getAngles() const { return angles; }
    void setAngles(const Vector3& newAngles) { angles = newAngles; }
};

// ===================================================================================
// ГЕОМЕТРИЯ ОБЪЕКТА
// ===================================================================================
class Geometry {
public:
    std::vector<Face> meshFaces;
    std::unordered_map<Vector3, Vector3, Vector3Hash, Vector3Equal> vertexNormals;
    Vector3 center;

    void computeCenter();
    void computeNormals();
    void buildFromData(const std::vector<double>& section,
        const std::vector<double>& trajectory,
        const std::vector<double>& scaling);

private:
    std::vector<Vector3> createRing(const std::vector<double>& section,
        const std::vector<double>& trajectory,
        const std::vector<double>& scaling,
        int ringIndex);
    void triangulateRings(const std::vector<std::vector<Vector3>>& rings);
};

// ===================================================================================
// МЕНЕДЖЕР ФАЙЛОВ
// ===================================================================================
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
class Scene {
public:
    Geometry geometry;
    std::vector<unsigned int> textureIDs;
    int textureIndex = -1;
    bool showWireframe = false;
    bool showNormals = false;
    bool smoothShading = false;
    bool usePerspective = true;
    unsigned short winWidth = 800, winHeight = 600;
    CameraController camera;
    AppConfig config;
    double heightRange = 30.0; // для текстур

    void reset();
    void printSceneInfo() const;
};

Scene gScene;

// ===================================================================================
// ВСПОМОГАТЕЛЬНЫЕ МАТЕМАТИЧЕСКИЕ ФУНКЦИИ
// ===================================================================================

Vector3 addVec(const Vector3& a, const Vector3& b) {
    return Vector3{ a.x + b.x, a.y + b.y, a.z + b.z };
}

Vector3 normalizeVec(const Vector3& a) {
    double len = std::sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
    if (len < 1e-8) return Vector3{ 0, 0, 0 };
    return Vector3{ a.x / len, a.y / len, a.z / len };
}

Vector3 calcNormal(const Vector3& v1, const Vector3& v2, const Vector3& v3) {
    double ux = v2.x - v1.x, uy = v2.y - v1.y, uz = v2.z - v1.z;
    double vx = v3.x - v1.x, vy = v3.y - v1.y, vz = v3.z - v1.z;
    Vector3 n{
        uy * vz - uz * vy,
        uz * vx - ux * vz,
        ux * vy - uy * vx
    };
    return normalizeVec(n);
}

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

    for (const auto& face : meshFaces) {
        Vector3 n = calcNormal(face.verts[0], face.verts[1], face.verts[2]);
        for (int i = 0; i < 3; ++i) {
            const Vector3& v = face.verts[i];
            tempNormals[v] = addVec(tempNormals[v], n);
        }
    }

    for (auto& it : tempNormals) {
        vertexNormals[it.first] = normalizeVec(it.second);
    }
}

std::vector<Vector3> Geometry::createRing(const std::vector<double>& section,
    const std::vector<double>& trajectory,
    const std::vector<double>& scaling,
    int ringIndex) {
    const int profilePoints = section.size() / 2;
    std::vector<Vector3> ring;
    ring.reserve(profilePoints);

    double tx = trajectory[ringIndex * 3 + 0];
    double ty = trajectory[ringIndex * 3 + 1];
    double tz = trajectory[ringIndex * 3 + 2];
    double scale = scaling[ringIndex];

    double secCx = 0.0, secCy = 0.0;
    for (int k = 0; k < profilePoints; ++k) {
        secCx += section[2 * k];
        secCy += section[2 * k + 1];
    }
    secCx /= profilePoints;
    secCy /= profilePoints;

    for (int k = 0; k < profilePoints; ++k) {
        double dx = section[2 * k] - secCx;
        double dy = section[2 * k + 1] - secCy;
        double x = dx * scale + secCx + tx;
        double y = dy * scale + secCy + ty;
        double z = tz;
        ring.emplace_back(Vector3{ x, y, z });
    }

    return ring;
}

void Geometry::triangulateRings(const std::vector<std::vector<Vector3>>& rings) {
    meshFaces.clear();

    // Боковая поверхность
    for (size_t i = 0; i < rings.size() - 1; ++i) {
        const std::vector<Vector3>& ring0 = rings[i];
        const std::vector<Vector3>& ring1 = rings[i + 1];
        int ringSize = ring0.size();

        for (int j = 0; j < ringSize; ++j) {
            int j1 = (j + 1) % ringSize;
            meshFaces.push_back(Face{ ring0[j], ring1[j], ring1[j1] });
            meshFaces.push_back(Face{ ring0[j], ring1[j1], ring0[j1] });
        }
    }

    // Торец в начале (обратная нормаль)
    if (!rings.empty() && !rings[0].empty()) {
        const auto& first = rings[0];
        Vector3 center{0,0,0};
        for (const auto& v : first) { center = addVec(center, v); }
        center.x /= first.size(); center.y /= first.size(); center.z /= first.size();
        for (size_t i = 0; i < first.size(); ++i) {
            size_t i1 = (i + 1) % first.size();
            meshFaces.push_back(Face{ center, first[i1], first[i] });
        }
    }

    // Торец в конце
    if (rings.size() > 1 && !rings.back().empty()) {
        const auto& last = rings.back();
        Vector3 center{0,0,0};
        for (const auto& v : last) { center = addVec(center, v); }
        center.x /= last.size(); center.y /= last.size(); center.z /= last.size();
        for (size_t i = 0; i < last.size(); ++i) {
            size_t i1 = (i + 1) % last.size();
            meshFaces.push_back(Face{ center, last[i], last[i1] });
        }
    }
}

void Geometry::buildFromData(const std::vector<double>& section,
    const std::vector<double>& trajectory,
    const std::vector<double>& scaling) {
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

    std::vector<std::vector<Vector3>> rings;
    rings.reserve(trajectoryPoints);

    for (int i = 0; i < trajectoryPoints; ++i) {
        rings.push_back(createRing(section, trajectory, scaling, i));
    }

    triangulateRings(rings);
    computeCenter();
    computeNormals();
}

// ===================================================================================
// РЕАЛИЗАЦИЯ МЕТОДОВ КЛАССА Scene
// ===================================================================================

void Scene::reset() {
    std::vector<double> section, trajectory, scaling;

    if (!FileManager::loadData(config.sectionPath, section, "even") ||
        !FileManager::loadData(config.trajectoryPath, trajectory, "triples") ||
        !FileManager::loadData(config.scalingPath, scaling)) {
        throw std::runtime_error("Ошибка загрузки файлов данных");
    }

    geometry.buildFromData(section, trajectory, scaling);

    // Вычисление heightRange для текстур
    if (!geometry.meshFaces.empty()) {
        double minZ = geometry.meshFaces[0].verts[0].z;
        double maxZ = minZ;
        for (const auto& f : geometry.meshFaces) {
            for (const auto& v : f.verts) {
                if (v.z < minZ) minZ = v.z;
                if (v.z > maxZ) maxZ = v.z;
            }
        }
        heightRange = (maxZ - minZ) > 1e-3 ? (maxZ - minZ) : 1.0;
    }

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

void applyTextureCoordinates(const Vector3& vertex, const Vector3& center, double heightRange) {
    double dx = vertex.x - center.x;
    double dy = vertex.y - center.y;
    double dz = vertex.z - center.z;

    double angle = std::atan2(dy, dx);
    double u = 0.5 + 0.5 * angle / M_PI;
    double v = heightRange > 1e-3 ? (dz) / heightRange : 0.0;

    glTexCoord2f(static_cast<float>(u), static_cast<float>(v));
}

void initializeTextureFromBmp(const std::string& filename) {
    unsigned int tex;
    unsigned char bytes[54];
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

    int offset;
    memcpy(&offset, bytes + 10, sizeof(int));
    int bmpw = *(int*)(bytes + 18);
    int bmph = *(int*)(bytes + 22);

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

    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, bmpw, bmph, 0, GL_BGR_EXT, GL_UNSIGNED_BYTE, pixels);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    gScene.textureIDs.push_back(tex);
    delete[] pixels;

    std::cout << "Загружена текстура: " << filename << " (" << bmpw << "x" << bmph << ")" << std::endl;
}

// ===================================================================================
// ОСВЕЩЕНИЕ — УЛУЧШЕННОЕ
// ===================================================================================

void initializeLighting() {
    // GL_LIGHT0 — мягкий белый, спереди
    float pos0[] = { 8.0f, 8.0f, 8.0f, 1.0f };
    float amb0[] = { 0.2f, 0.2f, 0.25f, 1.0f };
    float dif0[] = { 0.6f, 0.6f, 0.7f, 1.0f };
    float spec0[] = { 0.3f, 0.3f, 0.3f, 1.0f };

    glLightfv(GL_LIGHT0, GL_POSITION, pos0);
    glLightfv(GL_LIGHT0, GL_AMBIENT, amb0);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, dif0);
    glLightfv(GL_LIGHT0, GL_SPECULAR, spec0);

    // GL_LIGHT1 — тёплый оранжевый, снизу-справа
    float pos1[] = { 6.0f, -4.0f, 6.0f, 1.0f };
    float amb1[] = { 0.0f, 0.0f, 0.0f, 1.0f };
    float dif1[] = { 0.4f, 0.25f, 0.15f, 1.0f };
    float spec1[] = { 0.1f, 0.05f, 0.02f, 1.0f };

    glLightfv(GL_LIGHT1, GL_POSITION, pos1);
    glLightfv(GL_LIGHT1, GL_AMBIENT, amb1);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, dif1);
    glLightfv(GL_LIGHT1, GL_SPECULAR, spec1);

    // GL_LIGHT2 — холодный синий, сверху-сзади
    float pos2[] = { -5.0f, 6.0f, -8.0f, 1.0f };
    float amb2[] = { 0.0f, 0.0f, 0.0f, 1.0f };
    float dif2[] = { 0.2f, 0.25f, 0.4f, 1.0f };
    float spec2[] = { 0.1f, 0.1f, 0.2f, 1.0f };

    glLightfv(GL_LIGHT2, GL_POSITION, pos2);
    glLightfv(GL_LIGHT2, GL_AMBIENT, amb2);
    glLightfv(GL_LIGHT2, GL_DIFFUSE, dif2);
    glLightfv(GL_LIGHT2, GL_SPECULAR, spec2);

    // Материал
    float specular[] = { 0.8f, 0.8f, 0.9f, 1.0f };
    float shininess[] = { 32.0f };
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shininess);

    GLfloat no_ambient[] = { 0.0f, 0.0f, 0.0f, 1.0f };
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, no_ambient);
}

// ===================================================================================
// ПРОЕКЦИЯ И КАМЕРА
// ===================================================================================

void rebuildProjectionMatrix() {
    double ratio = (double)gScene.winWidth / (double)gScene.winHeight;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if (gScene.usePerspective)
        gluPerspective(60.0, ratio, 0.1, 100.0); // уменьшен угол до 60° для естественности
    else
        glOrtho(-ratio * 8, ratio * 8, -8.0, 8.0, -100.0, 100.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gScene.camera.setupView();
}

// ===================================================================================
// УПРАВЛЕНИЕ СЦЕНОЙ
// ===================================================================================

void changeCamera(double x, double y, double z) {
    gScene.camera.rotate(x, y);
    if (z != 0) gScene.camera.zoom(z);
    rebuildProjectionMatrix();
    glutPostRedisplay();
}

void switchLight(GLenum id, bool val) {
    if (val) glEnable(id); else glDisable(id);
    glutPostRedisplay();
}

void Mouse(int button, int state, int x, int y) {
    if (state == 0) {
        if (button == 3) changeCamera(0, 0, -1.0);
        else if (button == 4) changeCamera(0, 0, 1.0);
    }
}

void switchPerspective(bool val) {
    gScene.usePerspective = val;
    rebuildProjectionMatrix();
    glutPostRedisplay();
}

void switchSmooth(bool val) {
    gScene.smoothShading = val;
    glutPostRedisplay();
}

void switchCarcass(bool val) {
    gScene.showWireframe = val;
    glutPostRedisplay();
}

void switchShowNormals(bool val) {
    gScene.showNormals = val;
    glutPostRedisplay();
}

void switchTexture(int val) {
    gScene.textureIndex = val;
    if (gScene.textureIndex >= (int)gScene.textureIDs.size())
        gScene.textureIndex = -1;
    glutPostRedisplay();
}

// ===================================================================================
// РИСОВАНИЕ
// ===================================================================================

void drawCoordinateNetwork() {
    // Сетка
    glColor3f(0.3f, 0.3f, 0.3f);
    glBegin(GL_LINES);
    for (int i = -20; i <= 20; ++i) {
        glVertex3f(-20.0f, 0.0f, i);
        glVertex3f(20.0f, 0.0f, i);
        glVertex3f(i, 0.0f, -20.0f);
        glVertex3f(i, 0.0f, 20.0f);
    }
    glEnd();

    // Оси
    glLineWidth(2.0f);
    glBegin(GL_LINES);
    glColor3f(1.0f, 0.3f, 0.3f); // X — красная
    glVertex3f(0, 0, 0);
    glVertex3f(5, 0, 0);
    glColor3f(0.3f, 1.0f, 0.3f); // Y — зелёная
    glVertex3f(0, 0, 0);
    glVertex3f(0, 5, 0);
    glColor3f(0.3f, 0.3f, 1.0f); // Z — синяя
    glVertex3f(0, 0, 0);
    glVertex3f(0, 0, 5);
    glEnd();
    glLineWidth(1.0f);
}

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
                n = getNormal(normals, p);
            }
            else {
                n = calcNormal(face.verts[0], face.verts[1], face.verts[2]);
            }

            applyTextureCoordinates(p, gScene.geometry.center, gScene.heightRange);
            glNormal3d(n.x, n.y, n.z);
            glVertex3d(p.x, p.y, p.z);
        }
    }
    glEnd();

    if (gScene.textureIndex != -1) glDisable(GL_TEXTURE_2D);
    glDisable(GL_LIGHTING);
}

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

void Display(void) {
    glClearColor(0.05f, 0.05f, 0.1f, 1.0f); // тёмно-синий фон
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    drawCoordinateNetwork();
    if (gScene.showWireframe) drawCarcass();
    else drawFigure();

    if (gScene.showNormals) drawNormals();

    glutSwapBuffers();
}

void Reshape(GLint newW, GLint newH) {
    gScene.winWidth = newW;
    gScene.winHeight = newH;
    glViewport(0, 0, gScene.winWidth, gScene.winHeight);
    rebuildProjectionMatrix();
}

void Keyboard(unsigned char key, int x, int y) {
    switch (key) {
    case 27: // ESC (более переносимо, чем '\033')
        exit(0);
        break;
    case 'a': changeCamera(5, 0, 0); break; // влево
    case 'd': changeCamera(-5, 0, 0); break; // вправо
    case 'w': changeCamera(0, 5, 0); break; // вверх
    case 's': changeCamera(0, -5, 0); break; // вниз
    case 'z': changeCamera(0, 0, -0.5); break; // приближение
    case 'x': changeCamera(0, 0, 0.5); break;  // отдаление
    case '1': switchLight(GL_LIGHT0, !glIsEnabled(GL_LIGHT0)); break;
    case '2': switchLight(GL_LIGHT1, !glIsEnabled(GL_LIGHT1)); break;
    case '3': switchLight(GL_LIGHT2, !glIsEnabled(GL_LIGHT2)); break;
    case 'q': switchPerspective(!gScene.usePerspective); break;
    case 'e': switchSmooth(!gScene.smoothShading); break;
    case 'r': switchCarcass(!gScene.showWireframe); break;
    case 't': switchShowNormals(!gScene.showNormals); break;
    case 'y': switchTexture(gScene.textureIndex + 1); break;
    case 'i': gScene.printSceneInfo(); break;
    case ' ': gScene.reset(); break;
    default:
        break;
    }
}

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
    case 8: switchTexture(-1); break;
    case 9: switchTexture(0); break;
    case 10: switchTexture(1); break;
    case 11: switchLight(GL_LIGHT0, true); break;
    case 12: switchLight(GL_LIGHT0, false); break;
    case 13: switchLight(GL_LIGHT1, true); break;
    case 14: switchLight(GL_LIGHT1, false); break;
    case 15: switchLight(GL_LIGHT2, true); break;
    case 16: switchLight(GL_LIGHT2, false); break;
    case 17: gScene.reset(); break;
    case 18: gScene.printSceneInfo(); break;
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

    glutDisplayFunc(Display);
    glutReshapeFunc(Reshape);
    glutKeyboardFunc(Keyboard);
    glutMouseFunc(Mouse);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);
    glEnable(GL_LIGHT2);
    initializeLighting();

    GLfloat no_ambient[] = { 0.0f, 0.0f, 0.0f, 1.0f };
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, no_ambient);

    try {
        gScene.reset();
    }
    catch (const std::exception& e) {
        std::cerr << "Ошибка инициализации: " << e.what() << std::endl;
        return false;
    }

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
    setlocale(LC_ALL, "Ru");
    std::cout << "Лабораторная работа №3 - Тиражирование сечений" << std::endl;
    std::cout << "Управление:" << std::endl;
    std::cout << "  WASD - вращение камеры" << std::endl;
    std::cout << "  +/- - приближение/отдаление" << std::endl;
    std::cout << "  1/2/3 - вкл/выкл источники света" << std::endl;
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

    // Меню
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

    int menuLight0 = glutCreateMenu(Menu);
    glutAddMenuEntry("Включить", 11);
    glutAddMenuEntry("Выключить", 12);

    int menuLight1 = glutCreateMenu(Menu);
    glutAddMenuEntry("Включить", 13);
    glutAddMenuEntry("Выключить", 14);

    int menuLight2 = glutCreateMenu(Menu);
    glutAddMenuEntry("Включить", 15);
    glutAddMenuEntry("Выключить", 16);

    int menuSystem = glutCreateMenu(Menu);
    glutAddMenuEntry("Перезагрузить сцену", 17);
    glutAddMenuEntry("Информация о сцене", 18);

    glutCreateMenu(Menu);
    glutAddSubMenu("Проекция", menuProjection);
    glutAddSubMenu("Сглаживание", menuSmooth);
    glutAddSubMenu("Каркас", menuCarcass);
    glutAddSubMenu("Нормали", menuNormal);
    glutAddSubMenu("Текстура", menuTexture);
    glutAddSubMenu("Свет: LIGHT0", menuLight0);
    glutAddSubMenu("Свет: LIGHT1", menuLight1);
    glutAddSubMenu("Свет: LIGHT2", menuLight2);
    glutAddSubMenu("Система", menuSystem);
    glutAttachMenu(GLUT_RIGHT_BUTTON);

    std::cout << "Приложение успешно запущено" << std::endl;
    glutMainLoop();
    return 0;
}