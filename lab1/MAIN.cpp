/*
Отключено консольное (главное) окно:
	Linker ->  Advanced -> Entery Point := "mainCRTStartup"
	Linker ->  System -> SubSystem := "Windows (/SUBSYSTEM:WINDOWS)"
*/


#include <stdlib.h>
#include "glut.h"
#include <vector>
using namespace std;

GLint Width = 512, Height = 512;

struct Color {
	float r, g, b;
	Color(float r = 1.0f, float g = 0.0f, float b = 0.0f) : r(r), g(g), b(b) {}
};

struct Vertex {
	GLint x, y;
	Vertex(GLint x = 0, GLint y = 0) : x(x), y(y) {}
};

struct Primitive {
	vector<Vertex> vertices;
	Color color;
	Primitive(const Color& c = Color()) : color(c) {}
};

vector<Primitive> primitives; // Все примитивы
int currentPrimitive = -1; // Индекс текущего примитива


float NORMAL_POINT_SIZE = 5.0f;
float SELECTED_POINT_SIZE = 10.0f;

Color currentColor(1.0f, 0.0f, 0.0f); // Текущий цвет для новых примитивов

void Display(void)
{
	glClearColor(0.5, 0.5, 0.5, 1); glClear(GL_COLOR_BUFFER_BIT);

	for (size_t i = 0; i < primitives.size(); ++i) {
		const Primitive& prim = primitives[i];
		glColor3f(prim.color.r, prim.color.g, prim.color.b);
		if ((int)i == currentPrimitive) {
			glPointSize(SELECTED_POINT_SIZE);
		}
		else {
			glPointSize(NORMAL_POINT_SIZE);
		}

		// Рисуем точки
		glBegin(GL_POINTS);
		for (const auto& v : prim.vertices) {
			glVertex2i(v.x, v.y);
		}
		glEnd();

		// Рисуем линии (GL_LINES)
		glBegin(GL_LINES);
		for (size_t j = 0; j + 1 < prim.vertices.size(); j += 2) {
			glVertex2i(prim.vertices[j].x, prim.vertices[j].y);
			glVertex2i(prim.vertices[j + 1].x, prim.vertices[j + 1].y);
		}
		glEnd();
	}

	glFinish();
}

/* Функция изменения размеров окна */
void Reshape(GLint w, GLint h)
{
	Width = w;    Height = h;
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, w, 0, h, -1.0, 1.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void Keyboard(unsigned char key, int x, int y)
{
	if (currentPrimitive < 0 || currentPrimitive >= (int)primitives.size()) return;
	Primitive& prim = primitives[currentPrimitive];

	// Изменение цвета текущего примитива
	if (key == 'r') prim.color.r = min(prim.color.r + 0.1f, 1.0f);
	if (key == 'g') prim.color.g = min(prim.color.g + 0.1f, 1.0f);
	if (key == 'b') prim.color.b = min(prim.color.b + 0.1f, 1.0f);

	// Изменение координат всех вершин текущего примитива
	if (key == 'w') for (auto& v : prim.vertices) v.y += 5;
	if (key == 's') for (auto& v : prim.vertices) v.y -= 5;
	if (key == 'a') for (auto& v : prim.vertices) v.x -= 5;
	if (key == 'd') for (auto& v : prim.vertices) v.x += 5;

	// Управление размером точек
	if (key == '+' || key == '=') {
		NORMAL_POINT_SIZE = min(NORMAL_POINT_SIZE + 1.0f, 30.0f);
		SELECTED_POINT_SIZE = min(SELECTED_POINT_SIZE + 2.0f, 40.0f);
	}
	if (key == '-' || key == '_') {
		NORMAL_POINT_SIZE = max(NORMAL_POINT_SIZE - 1.0f, 1.0f);
		SELECTED_POINT_SIZE = max(SELECTED_POINT_SIZE - 2.0f, 2.0f);
	}

	glutPostRedisplay();
}

void Mouse(int button, int state, int x, int y)
{
	if (state != GLUT_DOWN) return;

	if (button == GLUT_LEFT_BUTTON) {
		// Если нет активного примитива, создаём новый
		if (currentPrimitive == -1 || currentPrimitive >= (int)primitives.size()) {
			primitives.push_back(Primitive(currentColor));
			currentPrimitive = (int)primitives.size() - 1;
		}
		primitives[currentPrimitive].vertices.push_back(Vertex(x, Height - y));
	}
	// Удаление последней вершины текущего примитива по средней кнопке мыши
	if (button == GLUT_MIDDLE_BUTTON) {
		if (currentPrimitive != -1 && !primitives[currentPrimitive].vertices.empty()) {
			primitives[currentPrimitive].vertices.pop_back();
			glutPostRedisplay();
		}
	}
	// Завершение текущего примитива и начало нового по правому клику
	if (button == GLUT_RIGHT_BUTTON) {
		if (currentPrimitive != -1 && primitives[currentPrimitive].vertices.size() > 0) {
			primitives.push_back(Primitive(currentColor));
			currentPrimitive = (int)primitives.size() - 1;
		}
	}
	glutPostRedisplay();
}

// Коды для пунктов меню
enum MenuActions {
	MENU_COLOR_R,
	MENU_COLOR_G,
	MENU_COLOR_B,
	MENU_MOVE_UP,
	MENU_MOVE_DOWN,
	MENU_MOVE_LEFT,
	MENU_MOVE_RIGHT,
	MENU_FINISH_PRIMITIVE,
	MENU_DELETE_LAST,
	MENU_DELETE_ALL,
	MENU_DELETE_LAST_VERTEX
};

void Menu(int action)
{
	switch (action) {
	case MENU_COLOR_R: Keyboard('r', 0, 0); break;
	case MENU_COLOR_G: Keyboard('g', 0, 0); break;
	case MENU_COLOR_B: Keyboard('b', 0, 0); break;
	case MENU_MOVE_UP: Keyboard('w', 0, 0); break;
	case MENU_MOVE_DOWN: Keyboard('s', 0, 0); break;
	case MENU_MOVE_LEFT: Keyboard('a', 0, 0); break;
	case MENU_MOVE_RIGHT: Keyboard('d', 0, 0); break;
	case MENU_FINISH_PRIMITIVE:
		if (currentPrimitive != -1 && primitives[currentPrimitive].vertices.size() > 0) {
			primitives.push_back(Primitive(currentColor));
			currentPrimitive = (int)primitives.size() - 1;
			glutPostRedisplay();
		}
		break;
	case MENU_DELETE_LAST:
		if (!primitives.empty()) {
			primitives.pop_back();
			if (primitives.empty()) currentPrimitive = -1;
			else currentPrimitive = (int)primitives.size() - 1;
			glutPostRedisplay();
		}
		break;
	case MENU_DELETE_ALL:
		primitives.clear();
		currentPrimitive = -1;
		glutPostRedisplay();
		break;
	case MENU_DELETE_LAST_VERTEX:
		if (currentPrimitive != -1 && !primitives[currentPrimitive].vertices.empty()) {
			primitives[currentPrimitive].vertices.pop_back();
			glutPostRedisplay();
		}
		break;
	}
}


void main(int argc, char* argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB);
	glutInitWindowSize(Width, Height);
	glutCreateWindow("GL_LINES редактор");

	int menu_color = glutCreateMenu(Menu);
	glutAddMenuEntry("R +", MENU_COLOR_R);
	glutAddMenuEntry("G +", MENU_COLOR_G);
	glutAddMenuEntry("B +", MENU_COLOR_B);

	int menu_move = glutCreateMenu(Menu);
	glutAddMenuEntry("Вверх", MENU_MOVE_UP);
	glutAddMenuEntry("Вниз", MENU_MOVE_DOWN);
	glutAddMenuEntry("Влево", MENU_MOVE_LEFT);
	glutAddMenuEntry("Вправо", MENU_MOVE_RIGHT);

	int menu = glutCreateMenu(Menu);
	glutAddSubMenu("Смена цвета", menu_color);
	glutAddSubMenu("Перемещение", menu_move);
	glutAddMenuEntry("Завершить примитив", MENU_FINISH_PRIMITIVE);
	glutAddMenuEntry("Удалить последнюю вершину", MENU_DELETE_LAST_VERTEX);
	glutAddMenuEntry("Удалить последний примитив", MENU_DELETE_LAST);
	glutAddMenuEntry("Удалить все примитивы", MENU_DELETE_ALL);
	glutAttachMenu(GLUT_RIGHT_BUTTON);

	glutDisplayFunc(Display);
	glutReshapeFunc(Reshape);
	glutKeyboardFunc(Keyboard);
	glutMouseFunc(Mouse);
	glutMainLoop();
}
