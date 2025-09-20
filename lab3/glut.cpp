#include "stdafx.h"
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdio>
#include <C:\Users\Vladi\Downloads\KG_Lab3_Kostovsky\GL\glut.h>
#include <GL/glu.h>
#include <GL/gl.h>

// Новый стиль: классы и структуры
struct Vector3 {
    double x, y, z;
    bool operator<(const Vector3& o) const {
        if (x != o.x) return x < o.x;
        if (y != o.y) return y < o.y;
        return z < o.z;
    }
    bool operator==(const Vector3& o) const {
        return x == o.x && y == o.y && z == o.z;
    }
};

// Поиск нормали по ключу через ручной перебор map (совместимо со старыми STL)
Vector3 getNormal(const std::map<Vector3, Vector3>& normals, const Vector3& key, const Vector3& def) {
    for (std::map<Vector3, Vector3>::const_iterator it = normals.begin(); it != normals.end(); ++it) {
        if (it->first == key) return it->second;
    }
    return def;
}

struct Face {
    Vector3 verts[3];
};

class Geometry {
public:
    std::vector<Face> meshFaces;
    std::map<Vector3, Vector3> vertexNormals;
    Vector3 center;
    void computeCenter();
    void computeNormals();
};

class FileManager {
public:
    static bool loadSection(const std::string& path, std::vector<double>& out);
    static bool loadTrajectory(const std::string& path, std::vector<double>& out);
    static bool loadScaling(const std::string& path, std::vector<double>& out);
};

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
    Vector3 camAngles = { 45, 45, 10 };
    void reset();
};

class Renderer {
public:
    static void drawAxes();
    static void drawMesh(const Geometry& geom, bool smooth, int texID);
    static void drawWire(const Geometry& geom);
    static void drawNormals(const Geometry& geom);
};

Scene gScene;

// --- Geometry methods ---
Vector3 addVec(const Vector3& a, const Vector3& b) {
    return Vector3{ a.x + b.x, a.y + b.y, a.z + b.z };
}

Vector3 normalizeVec(const Vector3& a) {
    double len = std::sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
    return Vector3{ a.x / len, a.y / len, a.z / len };
}

void Geometry::computeCenter() {
    double sx = 0, sy = 0, sz = 0;
    int cnt = 0;
    for (const auto& f : meshFaces) {
        for (const auto& v : f.verts) {
            sx += v.x; sy += v.y; sz += v.z; ++cnt;
        }
    }
    center = Vector3{ sx / cnt, sy / cnt, sz / cnt };
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

void Geometry::computeNormals() {
    vertexNormals.clear();
    std::map<Vector3, Vector3> tempNormals;
    for (const auto& face : meshFaces) {
        Vector3 n = calcNormal(face.verts[0], face.verts[1], face.verts[2]);
        for (int i = 0; i < 3; ++i) {
            const Vector3& v = face.verts[i];
            if (tempNormals.find(v) == tempNormals.end()) {
                tempNormals[v] = n;
            } else {
                tempNormals[v] = addVec(tempNormals[v], n);
            }
        }
    }
    for (auto& it : tempNormals) {
        vertexNormals[it.first] = normalizeVec(it.second);
    }
}

#define _USE_MATH_DEFINES
#include <cmath>
double degToRad(double deg) {
    return deg * 3.14159265358979323846 / 180.0;
}

// --- FileManager methods ---
bool FileManager::loadSection(const std::string& path, std::vector<double>& out) {
    FILE* f = std::fopen(path.c_str(), "r");
    if (!f) return false;
    double v;
    while (std::fscanf(f, "%lf", &v) == 1) {
        out.push_back(v);
    }
    std::fclose(f);
    return true;
}

bool FileManager::loadTrajectory(const std::string& path, std::vector<double>& out) {
    FILE* f = std::fopen(path.c_str(), "r");
    if (!f) return false;
    double v;
    while (std::fscanf(f, "%lf", &v) == 1) out.push_back(v);
    std::fclose(f);
    return true;
}

bool FileManager::loadScaling(const std::string& path, std::vector<double>& out) {
    FILE* f = std::fopen(path.c_str(), "r");
    if (!f) return false;
    double v;
    while (std::fscanf(f, "%lf", &v) == 1) out.push_back(v);
    std::fclose(f);
    return true;
}

void Scene::reset() {
    std::vector<double> section, trajectory, scaling;
    if (!FileManager::loadSection("files/section.txt", section) ||
        !FileManager::loadTrajectory("files/trajectory.txt", trajectory) ||
        !FileManager::loadScaling("files/scaling.txt", scaling)) {
        std::cerr << "Ошибка открытия файлов" << std::endl;
        exit(1);
    }

    // Проверка данных
    if (section.size() % 2 != 0) {
        std::cerr << "section.txt должен содержать четное число чисел!" << std::endl;
        exit(1);
    }
    if (trajectory.size() % 3 != 0) {
        std::cerr << "trajectory.txt должен содержать кратное 3 число чисел!" << std::endl;
        exit(1);
    }
    if (scaling.size() != trajectory.size() / 3) {
        std::cerr << "scaling.txt должен содержать столько же элементов, сколько точек в trajectory.txt!" << std::endl;
        exit(1);
    }

    geometry.meshFaces.clear();
    const int profilePoints = section.size() / 2; // Для треугольника 3 точки
    int n = scaling.size();

    std::vector< std::vector<Vector3> > rings;
    for (int i = 0; i < n; ++i) {
        double s = scaling[i];
        double tx = trajectory[i * 3 + 0];
        double ty = trajectory[i * 3 + 1];
        double tz = trajectory[i * 3 + 2];
        std::vector<Vector3> ring;
        for (int k = 0; k < profilePoints; ++k) {
            double x = section[2*k] * s + tx;
            double y = section[2*k+1] * s + ty;
            double z = tz;
            ring.push_back(Vector3{x, y, z});
        }
        rings.push_back(ring);
    }

    // Соединяем кольца треугольниками
    for (int i = 0; i < (int)rings.size() - 1; ++i) {
        const std::vector<Vector3>& ring0 = rings[i];
        const std::vector<Vector3>& ring1 = rings[i+1];
        int ringSize = ring0.size();
        for (int j = 0; j < ringSize; ++j) {
            int j1 = (j + 1) % ringSize;
            // Два треугольника на каждую "ячейку"
            geometry.meshFaces.push_back(Face{ ring0[j], ring1[j], ring1[j1] });
            geometry.meshFaces.push_back(Face{ ring0[j], ring1[j1], ring0[j1] });
        }
    }

    geometry.computeCenter();
    geometry.computeNormals();
}

// Инициализация текстуры из .bmp
void initializeTextureFromBmp(const std::string &filename) {
    unsigned int tex;
    unsigned char bytes[54];
    FILE *file = std::fopen(filename.c_str(), "rb");
    if (!file) return;
    int err = std::fread(bytes, 1, 54, file);
    if (err != 54) {
        std::cerr << "Ошибка чтения BMP заголовка" << std::endl;
        fclose(file);
        return;
    }
    
    int offset;
    memcpy(&offset, bytes + 10, sizeof(int));
    int bmpw = *(int*)(bytes + 18);
    int bmph = *(int*)(bytes + 22);
    
    unsigned char *pixels = new unsigned char[3 * bmpw * bmph];
    fseek(file, offset, SEEK_SET);
    err = std::fread(pixels, 3 * bmpw * bmph, 1, file);
    if (err != 1) {
        std::cerr << "Ошибка чтения BMP данных" << std::endl;
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
}

// Создание света и параметров материала
void initializeLighting() {
    // Источник света 1
    float ambientLight0[] = { 0.1f, 0.1f, 0.1f, 1.0f };
    float diffuseLight0[] = { 0.8f, 0.8f, 0.9f, 1.0f };
    float specularLight0[] = { 0.4f, 0.4f, 0.4f, 1.0f };
    float lightPosition0[] = { 10.0f, 10.0f, 10.0f, 1.0f };
    glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight0);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight0);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight0);
    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition0);

    // Установка параметров материала
    float specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
    float shininess[] = { 50.0f };
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shininess);
}

// Пересоздание матрицы проекции на основе режима проекции и позиции камеры
void rebuildProjectionMatrix() {
    double ratio = (double)gScene.winWidth / (double)gScene.winHeight;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if (gScene.usePerspective)
        gluPerspective(90.0, ratio, 0.1, 100.0);
    else
        glOrtho(-ratio, ratio, -ratio, ratio, -100.0, 100.0);
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    double cameraX = gScene.camAngles.z * cos(degToRad(gScene.camAngles.y)) * cos(degToRad(gScene.camAngles.x));
    double cameraY = gScene.camAngles.z * sin(degToRad(gScene.camAngles.y));
    double cameraZ = gScene.camAngles.z * cos(degToRad(gScene.camAngles.y)) * sin(degToRad(gScene.camAngles.x));
    gluLookAt(cameraX, cameraY, cameraZ, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
}

/*                                 */
/* Обработка действий пользователя */
/*                                 */

// Сдвиг камеры
void changeCamera(double x, double y, double z) {
    if (gScene.camAngles.y + y > 90) gScene.camAngles.y = 90;
    else if (gScene.camAngles.y + y < -90) gScene.camAngles.y = -90;
    else gScene.camAngles.y += y;
    
    gScene.camAngles.x = fmod(gScene.camAngles.x + x, 360.0);
    if (gScene.camAngles.x < 0) gScene.camAngles.x += 360.0;
    
    if (gScene.camAngles.z + z < 1) gScene.camAngles.z = 1;
    else gScene.camAngles.z += z;
    
    rebuildProjectionMatrix();
    glutPostRedisplay();
}

// Включить/выключить свет
void switchLight(GLenum id, bool val) {
    if (val)
        glEnable(id);
    else
        glDisable(id);
    glutPostRedisplay();
}

// --- Mouse wheel handler (GLUT) ---
void Mouse(int button, int state, int x, int y) {
    if (state == 0) { // GLUT_DOWN == 0
        if (button == 3) { // wheel up
            changeCamera(0, 0, -1.0); // Zoom in
        }
        else if (button == 4) { // wheel down
            changeCamera(0, 0, 1.0); // Zoom out
        }
    }
}

// Изменить режим проекции
void switchPerspective(bool val) {
    gScene.usePerspective = val;
    rebuildProjectionMatrix();
    glutPostRedisplay();
}

// Переключить глажение
void switchSmooth(bool val) {
    gScene.smoothShading = val;
    glutPostRedisplay();
}

// Переключить показание каркаса
void switchCarcass(bool val) {
    gScene.showWireframe = val;
    glutPostRedisplay();
}

// Переключить показание нормалей
void switchShowNormals(bool val) {
    gScene.showNormals = val;
    glutPostRedisplay();
}

// Переключить на следующую текстуру
void switchTexture(int val) {
    gScene.textureIndex = val;
    if (gScene.textureIndex >= (int)gScene.textureIDs.size()) gScene.textureIndex = -1;
    glutPostRedisplay();
}

/*                   */
/* Функции рисования */
/*                   */

// Рисование сети координат
void drawCoordinateNetwork() {
    // Цвет сетки
    glColor3ub(127, 127, 127);

    // Рисование полосочек по координатам
    glBegin(GL_LINES);
    for (float i = -100; i <= 100; i++) {
        glVertex3f(-100, 0, i);
        glVertex3f(100, 0, i);
        glVertex3f(i, 0, -100);
        glVertex3f(i, 0, 100);
    }

    // Рисование вертикальной полосочки в центре координат
    glVertex3f(0, -100, 0);
    glVertex3f(0, 100, 0);
    glEnd();
}

void drawFigure() {
    glColor3ub(200, 200, 200);
    const auto& faces = gScene.geometry.meshFaces;
    std::map<Vector3, Vector3>& normals = gScene.geometry.vertexNormals;
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
                n = getNormal(normals, p, Vector3{0,0,0});
            } else {
                n = calcNormal(face.verts[0], face.verts[1], face.verts[2]);
            }
            
            // Простые текстурные координаты
            if (j == 0) glTexCoord2f(0.0, 0.0);
            else if (j == 1) glTexCoord2f(1.0, 0.0);
            else glTexCoord2f(0.0, 1.0);
            
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
    std::map<Vector3, Vector3>& normals = gScene.geometry.vertexNormals;
    bool smooth = gScene.smoothShading;
    
    glBegin(GL_LINES);
    for (const auto& face : faces) {
        for (int j = 0; j < 3; j++) {
            const Vector3& p = face.verts[j];
            Vector3 n;
            if (smooth) {
                n = getNormal(normals, p, Vector3{0,0,0});
            } else {
                n = calcNormal(face.verts[0], face.verts[1], face.verts[2]);
            }
            
            glVertex3d(p.x, p.y, p.z);
            glVertex3d(p.x + n.x * 0.1, p.y + n.y * 0.1, p.z + n.z * 0.1);
        }
    }
    glEnd();
}

/*              */
/* Функции GLUT */
/*              */

// Функция рисования
void Display(void) {
    // Чистить экран
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Рисовать все
    drawCoordinateNetwork();
    if (gScene.showWireframe) drawCarcass();
    else drawFigure();
    
    if (gScene.showNormals) drawNormals();

    // Выводить на экран нарисованное
    glutSwapBuffers();
}

// Функция переразмеривания
void Reshape(GLint newW, GLint newH) {
    gScene.winWidth = newW;
    gScene.winHeight = newH;
    glViewport(0, 0, gScene.winWidth, gScene.winHeight);
    rebuildProjectionMatrix();
}

// Нажатие клавиатуры
void Keyboard(unsigned char key, int x, int y) {
    switch (key) {
        case '\033': // ESC
            exit(0); break;
        case 'a':
            changeCamera(5, 0, 0); break;
        case 'd':
            changeCamera(-5, 0, 0); break;
        case 'w':
            changeCamera(0, 5, 0); break;
        case 's':
            changeCamera(0, -5, 0); break;
        case '+': case '=':
            changeCamera(0, 0, -0.5); break;
        case '-': case '_':
            changeCamera(0, 0, 0.5); break;
        case '1':
            switchLight(GL_LIGHT0, !glIsEnabled(GL_LIGHT0));
            break;
        case 'q':
            switchPerspective(!gScene.usePerspective); break;
        case 'e':
            switchSmooth(!gScene.smoothShading); break;
        case 'r':
            switchCarcass(!gScene.showWireframe); break;
        case 't':
            switchShowNormals(!gScene.showNormals); break;
        case 'y':
            switchTexture(gScene.textureIndex + 1); break;
    }
}

// Нажатие меню
void Menu(int option) {
    switch (option) {
        case 0:
            switchPerspective(true);
            break;
        case 1:
            switchPerspective(false);
            break;
        case 2:
            switchSmooth(true);
            break;
        case 3:
            switchSmooth(false);
            break;
        case 4:
            switchCarcass(true);
            break;
        case 5:
            switchCarcass(false);
            break;
        case 6:
            switchShowNormals(true);
            break;
        case 7:
            switchShowNormals(false);
            break;
        case 8:
            switchTexture(-1);
            break;
        case 9:
            switchTexture(0);
            break;
        case 10:
            switchTexture(1);
            break;
        case 11:
            switchLight(GL_LIGHT0, true);
            break;
        case 12:
            switchLight(GL_LIGHT0, false);
            break;
    }
}

// Главняа функция
int main(int argc, char *argv[]) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize(gScene.winWidth, gScene.winHeight);
    glutCreateWindow("Практическое задание 2");

    glutDisplayFunc(Display);
    glutReshapeFunc(Reshape);
    glutKeyboardFunc(Keyboard);
    glutMouseFunc(Mouse); // wheel/кнопки мыши

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    gScene.reset();
    initializeLighting();
    GLfloat no_ambient[] = { 0.0f, 0.0f, 0.0f, 1.0f }; 
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, no_ambient); 
    initializeTextureFromBmp("files/1.bmp");
    initializeTextureFromBmp("files/2.bmp");

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
    glutCreateMenu(Menu);
    glutAddSubMenu("Выбор проекции", menuProjection);
    glutAddSubMenu("Сглаживание", menuSmooth);
    glutAddSubMenu("Отображение каркаса", menuCarcass);
    glutAddSubMenu("Отображение нормалей", menuNormal);
    glutAddSubMenu("Выбор текстуры", menuTexture);
    glutAddSubMenu("Источник света", menuLight1);
    glutAttachMenu(GLUT_RIGHT_BUTTON);

    glutMainLoop();
}