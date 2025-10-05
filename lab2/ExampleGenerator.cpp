#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>

using namespace std;
class ExampleGenerator {
private:
    struct Example {
        std::string name;
        std::vector<double> section;
        std::vector<double> trajectory;
        std::vector<double> scaling;
    };

    std::vector<Example> examples;
    int currentExample = 0;

public:
    ExampleGenerator() {
        initializeExamples();
    }

    void initializeExamples() {
        // Пример 1: Простой цилиндр
        examples.push_back({
            "Простой цилиндр",
            // Квадратное сечение
            {-1, -1, 1, -1, 1, 1, -1, 1},
            // Прямая траектория
            {0, 0, 0, 0, 0, 5, 0, 0, 10},
            // Постоянное масштабирование
            {1.0, 1.0, 1.0}
            });

        // Пример 2: Конус
        examples.push_back({
            "Конус",
            // Круглое сечение (8 точек)
            {1,0, 0.707,0.707, 0,1, -0.707,0.707, -1,0, -0.707,-0.707, 0,-1, 0.707,-0.707},
            // Вертикальная траектория
            {0, 0, 0, 0, 0, 3, 0, 0, 6},
            // Уменьшающееся масштабирование
            {1.0, 0.6, 0.2}
            });

        // Пример 3: Спираль
        examples.push_back({
            "Спираль",
            // Треугольное сечение
            {0, 1, -0.866, -0.5, 0.866, -0.5},
            // Спиральная траектория
            {0,0,0, 2,0,1, 0,2,2, -2,0,3, 0,-2,4, 2,0,5, 0,2,6},
            // Переменное масштабирование
            {1.0, 0.8, 1.2, 0.9, 1.1, 0.7, 1.0}
            });

        // Пример 4: Волнообразная труба
        examples.push_back({
            "Волнообразная труба",
            // Шестиугольное сечение
            {1,0, 0.5,0.866, -0.5,0.866, -1,0, -0.5,-0.866, 0.5,-0.866},
            // Синусоидальная траектория
            {0,0,0, 1,0,2, 2,1,4, 3,0,6, 4,-1,8, 5,0,10},
            // Постоянное масштабирование с небольшими вариациями
            {1.0, 1.0, 0.95, 1.0, 1.05, 1.0}
            });

        // Пример 5: Расширяющаяся воронка
        examples.push_back({
            "Расширяющаяся воронка",
            // Квадратное сечение
            {-0.5,-0.5, 0.5,-0.5, 0.5,0.5, -0.5,0.5},
            // Короткая вертикальная траектория
            {0,0,0, 0,0,1, 0,0,2, 0,0,3},
            // Увеличивающееся масштабирование
            {0.5, 1.0, 1.5, 2.0}
            });

        // Пример 6: Зигзагообразный туннель
        examples.push_back({
            "Зигзагообразный туннель",
            // Восьмиугольное сечение
            {1,0, 0.707,0.707, 0,1, -0.707,0.707, -1,0, -0.707,-0.707, 0,-1, 0.707,-0.707},
            // Зигзаг траектория
            {0,0,0, 2,0,2, -2,0,4, 2,0,6, -2,0,8, 0,0,10},
            // Переменное масштабирование
            {1.0, 0.9, 1.1, 0.8, 1.2, 1.0}
            });

        // Пример 7: Сложная пространственная кривая
        examples.push_back({
            "Сложная пространственная кривая",
            // Пятиугольное сечение
            {0,1, 0.951,0.309, 0.588,-0.809, -0.588,-0.809, -0.951,0.309},
            // Сложная 3D траектория
            {0,0,0, 1,1,1, -1,2,2, 2,-1,3, -2,1,4, 1,-2,5, 0,0,6},
            // Случайное масштабирование
            {1.0, 1.3, 0.7, 1.1, 0.9, 1.4, 0.8}
            });
    }

    bool saveExample(int index) {
        if (index < 0 || index >= examples.size()) {
            std::cout << "Неверный индекс примера!" << std::endl;
            return false;
        }

        const Example& ex = examples[index];

        // Создаем папку files если её нет
        std::filesystem::create_directory("files");

        // Сохраняем сечение
        std::ofstream sectionFile("files/section.txt");
        if (!sectionFile) {
            std::cout << "Ошибка создания section.txt" << std::endl;
            return false;
        }
        for (size_t i = 0; i < ex.section.size(); i++) {
            sectionFile << ex.section[i];
            if (i < ex.section.size() - 1) sectionFile << " ";
        }
        sectionFile.close();

        // Сохраняем траекторию
        std::ofstream trajectoryFile("files/trajectory.txt");
        if (!trajectoryFile) {
            std::cout << "Ошибка создания trajectory.txt" << std::endl;
            return false;
        }
        for (size_t i = 0; i < ex.trajectory.size(); i++) {
            trajectoryFile << ex.trajectory[i];
            if (i < ex.trajectory.size() - 1) trajectoryFile << " ";
        }
        trajectoryFile.close();

        // Сохраняем масштабирование
        std::ofstream scalingFile("files/scaling.txt");
        if (!scalingFile) {
            std::cout << "Ошибка создания scaling.txt" << std::endl;
            return false;
        }
        for (size_t i = 0; i < ex.scaling.size(); i++) {
            scalingFile << ex.scaling[i];
            if (i < ex.scaling.size() - 1) scalingFile << " ";
        }
        scalingFile.close();

        currentExample = index;
        std::cout << "Загружен пример: " << ex.name << std::endl;
        std::cout << "Точек сечения: " << ex.section.size() / 2 << std::endl;
        std::cout << "Точек траектории: " << ex.trajectory.size() / 3 << std::endl;
        std::cout << "Коэффициентов масштабирования: " << ex.scaling.size() << std::endl;

        return true;
    }

    void printAvailableExamples() {
        std::cout << "\n=== ДОСТУПНЫЕ ПРИМЕРЫ ===" << std::endl;
        for (size_t i = 0; i < examples.size(); i++) {
            std::cout << i + 1 << ". " << examples[i].name;
            if (i == currentExample) std::cout << " [ТЕКУЩИЙ]";
            std::cout << std::endl;
        }
        std::cout << "=========================" << std::endl;
    }

    void interactiveMenu() {
        int choice;
        do {
            printAvailableExamples();
            std::cout << "\nВыберите пример (1-" << examples.size() << ") или 0 для выхода: ";
            std::cin >> choice;

            if (choice == 0) break;

            if (choice >= 1 && choice <= examples.size()) {
                if (saveExample(choice - 1)) {
                    std::cout << "Пример успешно загружен!" << std::endl;
                }
            }
            else {
                std::cout << "Неверный выбор!" << std::endl;
            }

            std::cout << std::endl;
        } while (choice != 0);
    }

    bool loadNextExample() {
        currentExample = (currentExample + 1) % examples.size();
        return saveExample(currentExample);
    }

    bool loadPreviousExample() {
        currentExample = (currentExample - 1 + examples.size()) % examples.size();
        return saveExample(currentExample);
    }

    int getCurrentExampleIndex() const { return currentExample; }
    int getTotalExamples() const { return examples.size(); }
    std::string getCurrentExampleName() const {
        return examples[currentExample].name;
    }
};

// Глобальный объект для использования в основном приложении
ExampleGenerator gExampleGenerator;

// Функции для интеграции с основным приложением
void initializeExamples() {
    // Автоматически загружаем первый пример при инициализации
    gExampleGenerator.saveExample(0);
}

void switchToExample(int index) {
    gExampleGenerator.saveExample(index);
}

void loadNextExample() {
    gExampleGenerator.loadNextExample();
}

void loadPreviousExample() {
    gExampleGenerator.loadPreviousExample();
}

void showExamplesMenu() {
    gExampleGenerator.interactiveMenu();
}

// Пример использования в основном приложении

int main() {
    std::cout << "Генератор примеров для лабораторной работы №3" << std::endl;
    std::cout << "Создание файлов section.txt, trajectory.txt, scaling.txt" << std::endl;

    ExampleGenerator generator;
    generator.interactiveMenu();

    return 0;
}
