// --- ������ ��������� �� ����� ---
enum GridRenderMode {
    GRID_RENDER_CLASSIC, // �������
    GRID_RENDER_FILL     // ������������ �� �����
};
extern GridRenderMode gridRenderMode;

// --- ����� ��������� ������� ---
enum OutlineMode {
    OUTLINE_GL,
    OUTLINE_BRESENHAM,
    OUTLINE_WU
};

OutlineMode outlineMode = OUTLINE_GL;

#include <C:\Users\Vladi\Downloads\KG_Lab2_Kostovsky_var2\GL\glut.h>
#include <GL/glu.h>
#include <GL/gl.h>


// ���� ��� �������������� ��������� �������� ��������
bool isLoadingTexture = false;

// ������� ����������� ���������, ����� vector ��� ����� ��� struct � ���������� ����������
#include <math.h>
#include <cstdlib>
#include <stdio.h>
#include <locale.h>
#include <windows.h>
#include <commdlg.h> // ��� ������� ������ ����� (������ ���� ����� windows.h)

#include <vector>
using namespace std;

#define _CRT_SECURE_NO_WARNINGS
#define M_PI       3.14159265358979323846   // pi
    
// ������� ���������� ��� ��������� �����
void BresenhamLine(int x0, int y0, int x1, int y1) {
    int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
    int dy = -abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
    int err = dx + dy, e2;
    while (true) {
        glVertex2i(x0, y0);
        if (x0 == x1 && y0 == y1) break;
        e2 = 2 * err;
        if (e2 >= dy) { err += dy; x0 += sx; }
        if (e2 <= dx) { err += dx; y0 += sy; }
    }
}


// --- �������� �� ��� ���������� ����� ---
void plotWu(int x, int y, float c) {
    glColor4f(0.0f, 0.0f, 1.0f, c); // ����� ���� � ������ ��� ��
    glVertex2i(x, y);
}

float rfpart(float x) { return 1.0f - (x - floorf(x)); }
float fpart(float x) { return x - floorf(x); }

void WuLine(int x0, int y0, int x1, int y1) {
    bool steep = abs(y1 - y0) > abs(x1 - x0);
    if (steep) {
        std::swap(x0, y0);
        std::swap(x1, y1);
    }
    if (x0 > x1) {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }
    float dx = x1 - x0;
    float dy = y1 - y0;
    float gradient = dx == 0 ? 1 : dy / dx;

    // ������ �����
    int xend = x0;
    float yend = y0 + gradient * (xend - x0);
    float xgap = rfpart(x0 + 0.5f);
    int xpxl1 = xend;
    int ypxl1 = (int)floorf(yend);
    if (steep) {
        plotWu(ypxl1, xpxl1, rfpart(yend) * xgap);
        plotWu(ypxl1 + 1, xpxl1, fpart(yend) * xgap);
    }
    else {
        plotWu(xpxl1, ypxl1, rfpart(yend) * xgap);
        plotWu(xpxl1, ypxl1 + 1, fpart(yend) * xgap);
    }
    float intery = yend + gradient;

    // ������ �����
    xend = x1;
    yend = y1 + gradient * (xend - x1);
    xgap = fpart(x1 + 0.5f);
    int xpxl2 = xend;
    int ypxl2 = (int)floorf(yend);
    if (steep) {
        plotWu(ypxl2, xpxl2, rfpart(yend) * xgap);
        plotWu(ypxl2 + 1, xpxl2, fpart(yend) * xgap);
    }
    else {
        plotWu(xpxl2, ypxl2, rfpart(yend) * xgap);
        plotWu(xpxl2, ypxl2 + 1, fpart(yend) * xgap);
    }

    // �������� ����
    if (steep) {
        for (int x = xpxl1 + 1; x < xpxl2; ++x) {
            plotWu((int)floorf(intery), x, rfpart(intery));
            plotWu((int)floorf(intery) + 1, x, fpart(intery));
            intery += gradient;
        }
    }
    else {
        for (int x = xpxl1 + 1; x < xpxl2; ++x) {
            plotWu(x, (int)floorf(intery), rfpart(intery));
            plotWu(x, (int)floorf(intery) + 1, fpart(intery));
            intery += gradient;
        }
    }
}

struct Color {
    float r, g, b;
    Color(float r = 1.0f, float g = 0.0f, float b = 0.0f) : r(r), g(g), b(b) {}
};

struct Vertex {
    GLint x, y;
    Vertex(GLint x = 0, GLint y = 0) : x(x), y(y) {}
};

struct Pentagon {
    bool justColorReset;
    Vertex center; // ����� �������������
    float radius; // ������ (���������� �� ������ �� ������)
    Color color;
    bool filled;
    bool textured;
    GLuint textureID;
    float translateX, translateY;
    float rotateAngle;
    float scaleX, scaleY;
    bool selected; // ������ �� ������������

    Pentagon(const Vertex& c, float r, const Color& col = Color())
        : center(c), radius(r), color(col), filled(false), textured(false),
        textureID(0), translateX(0), translateY(0), rotateAngle(0), scaleX(1), scaleY(1),
        selected(false) {}
};




// ���������� ���������� ����������� ������ ����� ���� struct, �� ���� �������
GridRenderMode gridRenderMode = GRID_RENDER_CLASSIC;

// Forward declarations for GLUT callbacks
void Display();
void Reshape(int w, int h);
GLint Width = 800, Height = 600;
int gridRows = 20;
int gridCols = 20;
int currentPentagon = -1; // ������ �������� ����������� �������������
Vertex centerPoint;
bool creatingPentagon = false;
bool dragging = false;
int dragOffsetX = 0, dragOffsetY = 0;
vector<Pentagon> pentagons; // ��� �������������

enum BlendingMode {
    BLENDING_OFF,
    BLENDING_AND,
    BLENDING_NAND
};

BlendingMode currentBlending = BLENDING_OFF;

GLuint loadTexture(const char* filename) {
    FILE* file = fopen(filename, "rb");
    if (!file) {
        printf("�� ������� ������� ���� %s\n", filename);
        return 0;
    }


    GLubyte header[54];
    fread(header, 1, 54, file);

    // ���������, ��� ��� BMP-����
    if (header[0] != 'B' || header[1] != 'M') {
        printf("%s �� �������� BMP-������\n", filename);
        fclose(file);
        return 0;
    }

    // �������� ������� �� ��������� (� little-endian �������)
    int width = *(int*)&header[18];
    int height = *(int*)&header[22];
    int dataSize = width * height * 3;

    // ��������� �������
    if (width <= 0 || height <= 0 || dataSize <= 0) {
        printf("������������ ������� ����������� � ����� %s\n", filename);
        fclose(file);
        return 0;
    }

    // ������ �������
    GLubyte* pixels = new GLubyte[dataSize];
    fread(pixels, 1, dataSize, file);
    fclose(file);

    // ����������� ����������� �� ��������� (BMP ������ ������ ����� �����)
    GLubyte* flippedPixels = new GLubyte[dataSize];
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            for (int c = 0; c < 3; c++) {
                flippedPixels[(y * width + x) * 3 + c] =
                    pixels[((height - 1 - y) * width + x) * 3 + c];
            }
        }
    }

    GLuint texID;
    glGenTextures(1, &texID);
    glBindTexture(GL_TEXTURE_2D, texID);

    // ���������� glTexImage2D ������ gluBuild2DMipmaps ��� ������� ����������
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_BGR_EXT, GL_UNSIGNED_BYTE, flippedPixels);

    // ������������� ��������� ��������
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

    delete[] pixels;
    delete[] flippedPixels;
    return texID;
}

void drawGrid() {
    glColor3f(0.7f, 0.7f, 0.7f); // ���� �����
    glLineWidth(1.0f);
    float cellWidth = (float)Width / gridCols;
    float cellHeight = (float)Height / gridRows;

    glBegin(GL_LINES);
    // ������������ �����
    for (int i = 0; i <= gridCols; ++i) {
        float x = i * cellWidth;
        glVertex2f(x, 0);
        glVertex2f(x, Height);
    }
    // �������������� �����
    for (int j = 0; j <= gridRows; ++j) {
        float y = j * cellHeight;
        glVertex2f(0, y);
        glVertex2f(Width, y);
    }
    glEnd();
}

void drawPentagon(const Pentagon& pent, bool showCenter = false) {
    // ��������� ������� ��������� �������
    glPushMatrix();

    // ��������� ��������������
    glTranslatef(pent.translateX, pent.translateY, 0);
    glRotatef(pent.rotateAngle, 0, 0, 1);
    glScalef(pent.scaleX, pent.scaleY, 1);

    // ������ ����� �������������
    if (showCenter) {
        glColor3f(1.0f, 1.0f, 1.0f);
        glBegin(GL_POINTS);
        glVertex2f(pent.center.x, pent.center.y);
        glEnd();
    }

    // ������������� ����
    glColor3f(pent.color.r, pent.color.g, pent.color.b);

    // �������� �������� ������ � ������������ ������
    bool textureWasEnabled = false;
    if (gridRenderMode == GRID_RENDER_CLASSIC && pent.textured) {
        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, pent.textureID);
        textureWasEnabled = true;
    }

    // ������ ������������
    if (pent.filled) {
        if (gridRenderMode == GRID_RENDER_FILL) {
            // ������������ �� ����� � ���������� (��� �����������)
            float cellWidth = (float)Width / gridCols;
            float cellHeight = (float)Height / gridRows;
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            for (int i = 0; i < gridCols; ++i) {
                for (int j = 0; j < gridRows; ++j) {
                    float cx = (i + 0.5f) * cellWidth;
                    float cy = (j + 0.5f) * cellHeight;
                    // ����� ����� ����� (angle summation)
                    float sumAngles = 0.0f;
                    float angleStep = 2 * M_PI / 5.0f;
                    for (int k = 0; k < 5; ++k) {
                        float x0 = pent.center.x + pent.radius * cos(angleStep * k);
                        float y0 = pent.center.y + pent.radius * sin(angleStep * k);
                        float x1 = pent.center.x + pent.radius * cos(angleStep * ((k+1)%5));
                        float y1 = pent.center.y + pent.radius * sin(angleStep * ((k+1)%5));
                        float ax = x0 - cx, ay = y0 - cy;
                        float bx = x1 - cx, by = y1 - cy;
                        float dot = ax * bx + ay * by;
                        float lenA = sqrt(ax * ax + ay * ay);
                        float lenB = sqrt(bx * bx + by * by);
                        float cosTheta = dot / (lenA * lenB + 1e-6f);
                        // ���������� cosTheta � [-1, 1] ��� acos
                        if (cosTheta < -1.0f) cosTheta = -1.0f;
                        if (cosTheta > 1.0f) cosTheta = 1.0f;
                        float angle = acos(cosTheta);
                        // ���������� ���� ���� ����� ��������� ������������
                        float cross = ax * by - ay * bx;
                        if (cross < 0) angle = -angle;
                        sumAngles += angle;
                    }
                    bool inside = fabs(fabs(sumAngles) - 2 * M_PI) < 0.1f;
                    bool border = false;
                    if (!inside) {
                        float angleStep = 2 * M_PI / 5.0f;
                        for (int k = 0; k < 5; ++k) {
                            float px = pent.center.x + pent.radius * cos(angleStep * k);
                            float py = pent.center.y + pent.radius * sin(angleStep * k);
                            if (px >= i * cellWidth && px < (i+1) * cellWidth &&
                                py >= j * cellHeight && py < (j+1) * cellHeight) {
                                border = true;
                                break;
                            }
                        }
                    }
                    if (inside) {
                        glBegin(GL_QUADS);
                        glColor3f(pent.color.r, pent.color.g, pent.color.b);
                        glVertex2f(i * cellWidth, j * cellHeight);
                        glVertex2f((i+1) * cellWidth, j * cellHeight);
                        glVertex2f((i+1) * cellWidth, (j+1) * cellHeight);
                        glVertex2f(i * cellWidth, (j+1) * cellHeight);
                        glEnd();
                    } else if (border) {
                        glBegin(GL_QUADS);
                        glColor4f(pent.color.r, pent.color.g, pent.color.b, 0.4f);
                        glVertex2f(i * cellWidth, j * cellHeight);
                        glVertex2f((i+1) * cellWidth, j * cellHeight);
                        glVertex2f((i+1) * cellWidth, (j+1) * cellHeight);
                        glVertex2f(i * cellWidth, (j+1) * cellHeight);
                        glEnd();
                    }
                }
            }
            glDisable(GL_BLEND);
        } else {
            // ������� ������� (OpenGL fill)
            if (textureWasEnabled) {
                glColor3f(1.0f, 1.0f, 1.0f); // ���� �� ������ �� ��������, �� ����� � �������� �������
                glBegin(GL_POLYGON);
                for (int j = 0; j < 5; j++) {
                    float angle = 2 * M_PI * j / 5.0f;
                    float x = pent.center.x + pent.radius * cos(angle);
                    float y = pent.center.y + pent.radius * sin(angle);
                    float tx = 0.5f + 0.5f * cos(angle); // [0,1] �� �����
                    float ty = 0.5f + 0.5f * sin(angle);
                    glTexCoord2f(tx, ty);
                    glVertex2f(x, y);
                }
                glEnd();
            } else {
                glColor3f(pent.color.r, pent.color.g, pent.color.b);
                glBegin(GL_POLYGON);
                for (int j = 0; j < 5; j++) {
                    float angle = 2 * M_PI * j / 5.0f;
                    float x = pent.center.x + pent.radius * cos(angle);
                    float y = pent.center.y + pent.radius * sin(angle);
                    glVertex2f(x, y);
                }
                glEnd();
            }
        }
    }
    else {
        // ������ ������ ��������� �������� !!
        if (outlineMode == OUTLINE_GL) {
            glColor3f(0,0,0);
            glBegin(GL_LINE_LOOP);
            for (int j = 0; j < 5; j++) {
                float angle = 2 * M_PI * j / 5.0f;
                float x = pent.center.x + pent.radius * cos(angle);
                float y = pent.center.y + pent.radius * sin(angle);
                glVertex2f(x, y);
            }
            glEnd();
        } else if (outlineMode == OUTLINE_BRESENHAM) {
            glColor3f(1,0,0); // ������� ��� ����������
            glBegin(GL_POINTS);
            for (int j = 0; j < 5; j++) {
                float angle1 = 2 * M_PI * j / 5.0f;
                float angle2 = 2 * M_PI * ((j+1)%5) / 5.0f;
                int x0 = (int)(pent.center.x + pent.radius * cos(angle1) + 0.5f);
                int y0 = (int)(pent.center.y + pent.radius * sin(angle1) + 0.5f);
                int x1 = (int)(pent.center.x + pent.radius * cos(angle2) + 0.5f);
                int y1 = (int)(pent.center.y + pent.radius * sin(angle2) + 0.5f);
                BresenhamLine(x0, y0, x1, y1);
            }
            glEnd();
        } else if (outlineMode == OUTLINE_WU) {
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            glColor3f(0,0,1); // ����� ��� ��
            glBegin(GL_POINTS);
            for (int j = 0; j < 5; j++) {
                float angle1 = 2 * M_PI * j / 5.0f;
                float angle2 = 2 * M_PI * ((j+1)%5) / 5.0f;
                int x0 = (int)(pent.center.x + pent.radius * cos(angle1) + 0.5f);
                int y0 = (int)(pent.center.y + pent.radius * sin(angle1) + 0.5f);
                int x1 = (int)(pent.center.x + pent.radius * cos(angle2) + 0.5f);
                int y1 = (int)(pent.center.y + pent.radius * sin(angle2) + 0.5f);
                WuLine(x0, y0, x1, y1);
            }
            glEnd();
            glDisable(GL_BLEND);
        }
    }

    // ���� ������������ ������, ������ ���������
    if (pent.selected) {
        glColor3f(0.0f, 1.0f, 0.0f);
        glLineWidth(2.0f);
        glBegin(GL_LINE_LOOP);
        for (int j = 0; j < 5; j++) {
            float angle = 2 * M_PI * j / 5.0f;
            float x = pent.center.x + pent.radius * cos(angle);
            float y = pent.center.y + pent.radius * sin(angle);
            glVertex2f(x, y);
        }
        glEnd();
        glLineWidth(1.0f);
    }

    // ��������� ��������, ���� ���� ��������
    if (textureWasEnabled) {
        glDisable(GL_TEXTURE_2D);
    }

    // --- ���������� ������� ���������� ---
    // ���������� (AND/NAND) �������� ��� ������ ������, ����� GRID_RENDER_FILL
    if (gridRenderMode != GRID_RENDER_FILL) {
        glDisable(GL_BLEND); // blending ������ logic op
        if (currentBlending == BLENDING_OFF) {
            glDisable(GL_COLOR_LOGIC_OP);
        } else if (currentBlending == BLENDING_AND) {
            glEnable(GL_COLOR_LOGIC_OP);
            glLogicOp(GL_AND);
        } else if (currentBlending == BLENDING_NAND) {
            glEnable(GL_COLOR_LOGIC_OP);
            glLogicOp(GL_NAND);
        }
    } else {
        glDisable(GL_BLEND);
        glDisable(GL_COLOR_LOGIC_OP);
    }

    glPopMatrix();
}


// ���������, ����� �� ���� � ������������
bool isPointInPentagon(const Pentagon& pent, int x, int y) {
    // ������� ��������: ���� ���������� �� ����� �� ������ ������ �������
    float dx = x - pent.center.x;
    float dy = y - pent.center.y;
    float distance = sqrt(dx * dx + dy * dy);
    return distance < pent.radius * 1.2f; // ��������� ��������� �����
}

void Mouse(int button, int state, int x, int y) {
    y = Height - y; // ����������� ���������� Y

    if (button == GLUT_LEFT_BUTTON) {
        if (state == GLUT_DOWN) {
            if (!creatingPentagon) {
                // ���������, ����� �� ���� � ������������ ������������
                bool clickedOnPentagon = false;
                for (size_t i = 0; i < pentagons.size(); i++) {
                    pentagons[i].selected = false;
                }
                for (size_t i = 0; i < pentagons.size(); i++) {
                    if (isPointInPentagon(pentagons[i], x, y)) {
                        printf("[DEBUG] ���� �� ������������� ������������� #%zu\n", i);
                        currentPentagon = i;
                        pentagons[i].selected = true;
                        pentagons[i].justColorReset = true;
                        clickedOnPentagon = true;
                        // --- Drag&Drop ---
                        dragging = true;
                        dragOffsetX = x - pentagons[i].center.x;
                        dragOffsetY = y - pentagons[i].center.y;
                        break;
                    }
                }
                if (!clickedOnPentagon) {
                    // ������ ����: ��������� �����
                    printf("[DEBUG] ������ ����: ����� (%d, %d), ������ ����� ������������\n", x, y);
                    centerPoint = Vertex(x, y);
                    creatingPentagon = true;
                    currentPentagon = pentagons.size();
                    // �������� ��������� � ����
                    for (size_t i = 0; i < pentagons.size(); i++) pentagons[i].selected = false;
                    Pentagon np(centerPoint, 0, Color(1.0f, 1.0f, 1.0f));
                    np.filled = true;
                    np.selected = true;
                    np.justColorReset = true;
                    pentagons.push_back(np);
                }
            } else {
                // ������ ����: ��������� ������
                float dx = x - centerPoint.x;
                float dy = y - centerPoint.y;
                float radius = sqrt(dx * dx + dy * dy);
                printf("[DEBUG] ������ ����: ������=%.1f, ����������� ������������\n", radius);
                pentagons[currentPentagon].radius = radius;
                creatingPentagon = false;
                currentPentagon = -1;
                glutPostRedisplay();
            }
        } else if (state == GLUT_UP) {
            // ��������� ��� � �������� drag
            dragging = false;
        }
    }
    // ������� ������: ������� ���������� ������������ ������ ��� ����������
    else if (button == GLUT_MIDDLE_BUTTON && state == GLUT_UP) {
        if (currentPentagon != -1 && currentPentagon < pentagons.size()) {
            pentagons.erase(pentagons.begin() + currentPentagon);
            if (pentagons.empty()) {
                currentPentagon = -1;
            }
            else if (currentPentagon >= pentagons.size()) {
                currentPentagon = pentagons.size() - 1;
            }
            glutPostRedisplay();
        }
    }
    // ������ ������: �������� ����
    else if (button == GLUT_RIGHT_BUTTON) {
        glutMenuStatusFunc(nullptr);
        glutPostRedisplay();
    }
}

// --- Drag&Drop ��������� �������� ���� ---
void Motion(int x, int y) {
    y = Height - y;
    if (dragging && currentPentagon != -1 && currentPentagon < pentagons.size()) {
        pentagons[currentPentagon].center.x = x - dragOffsetX;
        pentagons[currentPentagon].center.y = y - dragOffsetY;
        glutPostRedisplay();
    }
}

enum MenuActions {
    MENU_GRID_RENDER_CLASSIC,
    MENU_GRID_RENDER_FILL,
    MENU_COLOR_R,
    MENU_COLOR_G,
    MENU_COLOR_B,
    MENU_COLOR_RESET,
    MENU_TRANSLATE_UP,
    MENU_TRANSLATE_DOWN,
    MENU_TRANSLATE_LEFT,
    MENU_TRANSLATE_RIGHT,
    MENU_ROTATE_LEFT,
    MENU_ROTATE_RIGHT,
    MENU_SCALE_UP,
    MENU_SCALE_DOWN,
    MENU_FILLED,
    MENU_OUTLINE,
    MENU_BLENDING_OFF,
    MENU_BLENDING_AND,
    MENU_BLENDING_NAND,
    MENU_LOAD_TEXTURE,
    MENU_DELETE_CURRENT,
    MENU_DELETE_ALL,
    MENU_SELECT_NEXT,
    MENU_SELECT_PREV,
    MENU_GRID_10,
    MENU_GRID_20,
    MENU_GRID_40,
    MENU_GRID_100,
    MENU_GRID_1000,
    MENU_OUTLINE_GL,
    MENU_OUTLINE_BRESENHAM,
    MENU_OUTLINE_WU,
};

void Menu(int action) {
    switch (action) {
    case MENU_GRID_RENDER_CLASSIC:
        gridRenderMode = GRID_RENDER_CLASSIC;
        break;
    case MENU_GRID_RENDER_FILL:
        gridRenderMode = GRID_RENDER_FILL;
        break;
    case MENU_OUTLINE_GL:
        outlineMode = OUTLINE_GL;
        break;
    case MENU_OUTLINE_BRESENHAM:
        outlineMode = OUTLINE_BRESENHAM;
        break;
    case MENU_OUTLINE_WU:
        outlineMode = OUTLINE_WU;
        break;
    case MENU_COLOR_R:
        if (currentPentagon != -1 && currentPentagon < pentagons.size() && !isLoadingTexture) {
            isLoadingTexture = true;
            pentagons[currentPentagon].color.r = min(pentagons[currentPentagon].color.r + 0.1f, 1.0f);
        }
        break;
    case MENU_COLOR_G:
        if (currentPentagon != -1 && currentPentagon < pentagons.size()) {
            pentagons[currentPentagon].color.g = min(pentagons[currentPentagon].color.g + 0.1f, 1.0f);
        }
        break;
    case MENU_COLOR_B:
        if (currentPentagon != -1 && currentPentagon < pentagons.size()) {
            pentagons[currentPentagon].color.b = min(pentagons[currentPentagon].color.b + 0.1f, 1.0f);
        }
        break;
    case MENU_COLOR_RESET:
        if (currentPentagon != -1 && currentPentagon < pentagons.size()) {
            pentagons[currentPentagon].color = Color(0.0f, 0.0f, 0.0f);
        }
        break;
    case MENU_TRANSLATE_UP:
        if (currentPentagon != -1 && currentPentagon < pentagons.size()) {
            pentagons[currentPentagon].translateY += 5;
        }
        break;
    case MENU_TRANSLATE_DOWN:
        if (currentPentagon != -1 && currentPentagon < pentagons.size()) {
            pentagons[currentPentagon].translateY -= 5;
        }
        break;
    case MENU_TRANSLATE_LEFT:
        if (currentPentagon != -1 && currentPentagon < pentagons.size()) {
            pentagons[currentPentagon].translateX -= 5;
        }
        break;
    case MENU_TRANSLATE_RIGHT:
        if (currentPentagon != -1 && currentPentagon < pentagons.size()) {
            pentagons[currentPentagon].translateX += 5;
        }
        break;
    case MENU_ROTATE_LEFT:
        if (currentPentagon != -1 && currentPentagon < pentagons.size()) {
            pentagons[currentPentagon].rotateAngle -= 10;
        }
        break;
    case MENU_ROTATE_RIGHT:
        if (currentPentagon != -1 && currentPentagon < pentagons.size()) {
            pentagons[currentPentagon].rotateAngle += 10;
        }
        break;
    case MENU_SCALE_UP:
        if (currentPentagon != -1 && currentPentagon < pentagons.size()) {
            pentagons[currentPentagon].scaleX *= 1.1f;
            pentagons[currentPentagon].scaleY *= 1.1f;
        }
        break;
    case MENU_SCALE_DOWN:
        if (currentPentagon != -1 && currentPentagon < pentagons.size()) {
            pentagons[currentPentagon].scaleX /= 1.1f;
            pentagons[currentPentagon].scaleY /= 1.1f;
        }
        break;
    case MENU_FILLED:
        if (currentPentagon != -1 && currentPentagon < pentagons.size()) {
            pentagons[currentPentagon].filled = true;
        }
        break;
    case MENU_OUTLINE:
        if (currentPentagon != -1 && currentPentagon < pentagons.size()) {
            pentagons[currentPentagon].filled = false;
        }
        break;
    case MENU_BLENDING_OFF:
        currentBlending = BLENDING_OFF;
        break;
    case MENU_BLENDING_AND:
        currentBlending = BLENDING_AND;
        break;
    case MENU_BLENDING_NAND:
        currentBlending = BLENDING_NAND;
        break;
    case MENU_LOAD_TEXTURE:
        if (currentPentagon != -1 && currentPentagon < pentagons.size()) {
#ifdef _WIN32
            char fileName[MAX_PATH] = "";
            OPENFILENAMEA ofn = {0};
            ofn.lStructSize = sizeof(ofn);
            ofn.hwndOwner = NULL;
            ofn.lpstrFilter = "BMP Files\0*.bmp\0All Files\0*.*\0";
            ofn.lpstrFile = fileName;
            ofn.nMaxFile = MAX_PATH;
            ofn.Flags = OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST;
            ofn.lpstrTitle = "�������� BMP ��������";
            if (GetOpenFileNameA(&ofn)) {
                GLuint texId = loadTexture(fileName);
                if (texId != 0) {
                    pentagons[currentPentagon].textureID = texId;
                    pentagons[currentPentagon].textured = true;
                } else {
                    MessageBoxA(NULL, "������ �������� BMP-����� ��� ���� ��������!", "������", MB_ICONERROR | MB_OK);
                }
            }
            isLoadingTexture = false;
#else
            // ��� ������ �������� � �������� ������ ���������
            pentagons[currentPentagon].textureID = loadTexture("texture.bmp");
            if (pentagons[currentPentagon].textureID != 0) {
                pentagons[currentPentagon].textured = true;
            }
#endif
        }
        break;
    case MENU_DELETE_CURRENT:
        if (currentPentagon != -1 && currentPentagon < pentagons.size()) {
            pentagons.erase(pentagons.begin() + currentPentagon);
            if (pentagons.empty()) {
                currentPentagon = -1;
            }
            else if (currentPentagon >= pentagons.size()) {
                currentPentagon = pentagons.size() - 1;
            }
        }
        break;
    case MENU_DELETE_ALL:
        pentagons.clear();
        currentPentagon = -1;
        break;
    case MENU_SELECT_NEXT:
        if (pentagons.size() > 0) {
            currentPentagon = (currentPentagon + 1) % pentagons.size();
        }
        break;
    case MENU_SELECT_PREV:
        if (pentagons.size() > 0) {
            currentPentagon = (currentPentagon - 1 + pentagons.size()) % pentagons.size();
        }
        break;
    case MENU_GRID_10:
        gridRows = 10; gridCols = 10;
        break;
    case MENU_GRID_20:
        gridRows = 20; gridCols = 20;
        break;
    case MENU_GRID_40:
        gridRows = 40; gridCols = 40;
        break;
    case MENU_GRID_100:
        gridRows = 100; gridCols = 100;
        break;
    case MENU_GRID_1000:
        gridRows = 1000; gridCols = 1000;
        break;
    }
    glutPostRedisplay();
}

void Keyboard(unsigned char key, int x, int y) {
    if (currentPentagon < 0 || currentPentagon >= (int)pentagons.size()) return;
    Pentagon& pent = pentagons[currentPentagon];

    // �������������
    if (key == 'w') pent.translateY += 5;
    if (key == 's') pent.translateY -= 5;
    if (key == 'a') pent.translateX -= 5;
    if (key == 'd') pent.translateX += 5;
    if (key == 'q') pent.rotateAngle -= 10;
    if (key == 'e') pent.rotateAngle += 10;
    if (key == '+') {
        pent.scaleX *= 1.1f;
        pent.scaleY *= 1.1f;
    }
    if (key == '-') {
        pent.scaleX /= 1.1f;
        pent.scaleY /= 1.1f;
    }


    // ����: ������� ��������� ����� ������, ����� ������� �����������������
    if (key == 'r') {
        if (pent.justColorReset) {
            pent.color.r = 0.0f;
            pent.color.g = 0.0f;
            pent.color.b = 0.0f;
            pent.justColorReset = false;
        }
        pent.color.r += 0.1f;
        if (pent.color.r > 1.0f) pent.color.r = 0.0f;
    }
    if (key == 'g') {
        if (pent.justColorReset) {
            pent.color.r = 0.0f;
            pent.color.g = 0.0f;
            pent.color.b = 0.0f;
            pent.justColorReset = false;
        }
        pent.color.g += 0.1f;
        if (pent.color.g > 1.0f) pent.color.g = 0.0f;
    }
    if (key == 'b') {
        if (pent.justColorReset) {
            pent.color.r = 0.0f;
            pent.color.g = 0.0f;
            pent.color.b = 0.0f;
            pent.justColorReset = false;
        }
        pent.color.b += 0.1f;
        if (pent.color.b > 1.0f) pent.color.b = 0.0f;
    }

    // ����� ����������
    if (key == 'f') pent.filled = true;
    if (key == 'o') pent.filled = false;

    // ������������ ����� ���������������
    if (key == '[') {
        if (pentagons.size() > 0) {
            currentPentagon = (currentPentagon - 1 + pentagons.size()) % pentagons.size();
        }
    }
    if (key == ']') {
        if (pentagons.size() > 0) {
            currentPentagon = (currentPentagon + 1) % pentagons.size();
        }
    }

    glutPostRedisplay();
}

void Display() {
    glClear(GL_COLOR_BUFFER_BIT);
    drawGrid();
    for (const auto& pent : pentagons) {
        drawPentagon(pent);
    }
    glutSwapBuffers();
}

void Reshape(int w, int h) {
    Width = w;
    Height = h;
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, Width, 0, Height, -1.0, 1.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glutPostRedisplay();
}

void main(int argc, char* argv[]) {
    int menu_grid_render = glutCreateMenu(Menu);
    glutAddMenuEntry("������� ������", MENU_GRID_RENDER_CLASSIC);
    glutAddMenuEntry("������������ �� �����", MENU_GRID_RENDER_FILL);
    setlocale(LC_ALL, "Ru");
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE); // �������� ������� �����������
    glutInitWindowSize(Width, Height);
    glutCreateWindow("�������� ���������� ��������������");

    // ������� ����
    int menu_color = glutCreateMenu(Menu);
    glutAddMenuEntry("������� +", MENU_COLOR_R);
    glutAddMenuEntry("������� +", MENU_COLOR_G);
    glutAddMenuEntry("����� +", MENU_COLOR_B);
    glutAddMenuEntry("�������� ����", MENU_COLOR_RESET);

    int menu_transform = glutCreateMenu(Menu);
    glutAddMenuEntry("�������� �����", MENU_TRANSLATE_UP);
    glutAddMenuEntry("�������� ����", MENU_TRANSLATE_DOWN);
    glutAddMenuEntry("�������� �����", MENU_TRANSLATE_LEFT);
    glutAddMenuEntry("�������� ������", MENU_TRANSLATE_RIGHT);
    glutAddMenuEntry("��������� ������ �������", MENU_ROTATE_LEFT);
    glutAddMenuEntry("��������� �� �������", MENU_ROTATE_RIGHT);
    glutAddMenuEntry("��������� �������", MENU_SCALE_UP);
    glutAddMenuEntry("��������� �������", MENU_SCALE_DOWN);

    int menu_fill = glutCreateMenu(Menu);
    glutAddMenuEntry("�����������", MENU_FILLED);
    glutAddMenuEntry("������", MENU_OUTLINE);

    int menu_blend = glutCreateMenu(Menu);
    glutAddMenuEntry("��������� ����������", MENU_BLENDING_OFF);
    glutAddMenuEntry("���������� AND", MENU_BLENDING_AND);
    glutAddMenuEntry("���������� NAND", MENU_BLENDING_NAND);

    int menu_select = glutCreateMenu(Menu);
    glutAddMenuEntry("��������� ������������", MENU_SELECT_NEXT);
    glutAddMenuEntry("���������� ������������", MENU_SELECT_PREV);

    int menu_grid = glutCreateMenu(Menu);
    glutAddMenuEntry("10 x 10", MENU_GRID_10);
    glutAddMenuEntry("20 x 20", MENU_GRID_20);
    glutAddMenuEntry("40 x 40", MENU_GRID_40);
    glutAddMenuEntry("100 x 100", MENU_GRID_100);
    glutAddMenuEntry("1000 x 1000", MENU_GRID_1000);

    int menu_outline = glutCreateMenu(Menu);
    glutAddMenuEntry("GL_LINE_LOOP", MENU_OUTLINE_GL);
    glutAddMenuEntry("���������", MENU_OUTLINE_BRESENHAM);
    glutAddMenuEntry("�� (����������, �����)", MENU_OUTLINE_WU);

    int menu = glutCreateMenu(Menu);
    glutAddSubMenu("����", menu_color);
    glutAddSubMenu("�������������", menu_transform);
    glutAddSubMenu("����� ���������", menu_fill);
    glutAddSubMenu("������ (��������)", menu_outline);
    glutAddSubMenu("����� ������������", menu_grid_render);
    glutAddSubMenu("����������", menu_blend);
    glutAddSubMenu("������ �����", menu_grid);
    glutAddSubMenu("����� �������������", menu_select);
    glutAddMenuEntry("��������� ��������", MENU_LOAD_TEXTURE);
    glutAddMenuEntry("������� ������� ������������", MENU_DELETE_CURRENT);
    glutAddMenuEntry("������� ��� �������������", MENU_DELETE_ALL);
    glutAttachMenu(GLUT_RIGHT_BUTTON);

    glutDisplayFunc(Display);
    glutReshapeFunc(Reshape);
    glutKeyboardFunc(Keyboard);
    glutMouseFunc(Mouse);
    glutMotionFunc(Motion); // Drag&Drop

    // �������������� OpenGL
    glClearColor(0.9f, 0.9f, 0.9f, 1.0f);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, Width, 0, Height, -1.0, 1.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glutMainLoop();
}

