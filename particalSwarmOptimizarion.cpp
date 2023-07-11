#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <algorithm>
#include <limits.h>
#include <GL/glut.h>

using std::vector;
using std::pair;
using std::string;

struct Particle {
    Particle() {}
    Particle(vector<double> velocity, vector<double> position,
    vector<double> best_position, double fitness, 
    double best_fitness
    ) : velocity(velocity), position(position), best_position(best_position),
    fitness(fitness), best_fitness(best_fitness) {}

    vector<double> velocity;
    vector<double> position;
    vector<double> best_position;
    double fitness;
    double best_fitness;
};

class ParticalSwarm {
public:
    ParticalSwarm(
        int number_particles, int number_iterations, int dimensions, 
        double minX, double maxX) : number_particles(number_particles),
        number_iterations(number_iterations),
        dimensions(dimensions), minX(minX), maxX(maxX)
    {
        initializeSwarm();
    }

    void initializeSwarm() {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> disX(minX, maxX);
        std::uniform_real_distribution<double> disY(minV, maxV);

        // Проходимся по каждой частице
        for (int i = 0; i < number_particles; ++i) {
            Particle particle;
            particle.position.resize(dimensions);
            particle.velocity.resize(dimensions);

            // Инициализируем начальные позицию и скорость
            for (int j = 0; j < dimensions; ++j) {
                particle.position[j] = disX(gen);
                particle.velocity[j] = 0.0;
            }

            // Инициализация значения функции в позиции (точке)
            particle.fitness = fitness_function(particle.position);
            particle.best_position = particle.position;
            particle.best_fitness = particle.fitness;

            max_swarm.push_back(particle);
            min_swarm.push_back(particle);

            // После инициализации проверяем текущее значение по глобальным лучшим
            if (min_swarm[i].fitness < bestMinGlobalFitness) {
                bestMinGlobalFitness = min_swarm[i].fitness;
                bestMinGlobalPosition = min_swarm[i].position;
            } else if (max_swarm[i].fitness > bestMaxGlobalFitness) {
                bestMaxGlobalFitness = max_swarm[i].fitness;
                bestMaxGlobalPosition = max_swarm[i].position;
            }
        }
    }

    void set_function(string func) { }
    double fitness_function(vector<double> args) {
        return 3.0 + (args[0] * args[0]) + (args[1] * args[1]);
    }

    void calculate_min() {
        for (int iter = 0; iter < number_iterations; ++iter) {
            for (int i = 0; i < number_particles; ++i) {
                Particle& particle = min_swarm[i];

                // Обновим скорость
                for (int j = 0; j < dimensions; ++j) {
                    double r1 = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
                    double r2 = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);

                    particle.velocity[j] = w * particle.velocity[j]
                                        + c1 * r1 * (particle.best_position[j] - particle.position[j])
                                        + c2 * r2 * (bestMinGlobalPosition[j] - particle.position[j]);
                }

                // Обновим позицию
                for (int j = 0; j < dimensions; ++j) {
                    particle.position[j] = particle.position[j] + particle.velocity[j];

                    // Ограничиваем по max и min
                    particle.position[j] = std::max(minX, std::min(maxX, particle.position[j]));
                }

                particle.fitness = fitness_function(particle.position);

                // Обновим лучший результат для частицы
                if (particle.fitness < particle.best_fitness) {
                    particle.best_position = particle.position;
                    particle.best_fitness = particle.fitness;
                }

                // Обновим лучший глобальный результат
                if (particle.fitness < bestMinGlobalFitness) {
                    bestMinGlobalPosition = particle.position;
                    bestMinGlobalFitness = particle.fitness;
                }
            }
        }
    }

    void calculate_max() {
        for (int iter = 0; iter < number_iterations; ++iter) {
            for (int i = 0; i < number_particles; ++i) {
                Particle& particle = max_swarm[i];

                // Обновим скорость
                for (int j = 0; j < dimensions; ++j) {
                    double r1 = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
                    double r2 = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);

                    particle.velocity[j] = w * particle.velocity[j]
                                        + c1 * r1 * (particle.best_position[j] - particle.position[j])
                                        + c2 * r2 * (bestMaxGlobalPosition[j] - particle.position[j]);
                }

                // Обновим позицию
                for (int j = 0; j < dimensions; ++j) {
                    particle.position[j] = particle.position[j] + particle.velocity[j];

                    // Ограничиваем по max и min
                    particle.position[j] = std::max(minX, std::min(maxX, particle.position[j]));
                }

                particle.fitness = fitness_function(particle.position);

                // Обновим лучший результат для частицы
                if (particle.fitness > particle.best_fitness) {
                    particle.best_position = particle.position;
                    particle.best_fitness = particle.fitness;
                }

                // Обновим лучший глобальный результат
                if (particle.fitness > bestMaxGlobalFitness) {
                    bestMaxGlobalPosition = particle.position;
                    bestMaxGlobalFitness = particle.fitness;
                }
            }
        }
    }

    pair<vector<double>, double> get_min_position_and_fitness() {
        return {bestMinGlobalPosition, bestMinGlobalFitness};
    }

    pair<vector<double>, double> get_max_position_and_fitness() {
        return {bestMaxGlobalPosition, bestMaxGlobalFitness};
    }

    void make_graphic() {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        double minX = this->minX;
        double minY = this->minX;
        double maxX = this->maxX;
        double maxY = this->maxX;
        glOrtho(minX, maxX, minY, maxY, -1.0, 1.0);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        // Рисуем трехмерный график функции
        glBegin(GL_POINTS);
        double step = 0.1;
        for (double x = minX; x <= maxX; x += step) {
            for (double y = minY; y <= maxY; y += step) {
                double z = fitness_function({x, y});
                glVertex3d(x, y, z);
            }
        }
        glEnd();

        // Отмечаем минимум
        double min_x = bestMinGlobalPosition[0];
        double min_y = bestMinGlobalPosition[1];
        double min_z = fitness_function({min_x, min_y});
        glColor3d(1.0, 0.0, 0.0);
        glPointSize(5.0);
        glBegin(GL_POINTS);
        glVertex3d(min_x, min_y, min_z);
        glEnd();

        glutSwapBuffers();
    }

    static ParticalSwarm* instance;

    static void make_graphic_static() {
        instance->make_graphic();
    }

private:
    int number_particles;
    int number_iterations;
    int dimensions;
    
    double minX;
    double maxX;
    double minV;
    double maxV;

    double w = 0.729; // весовая доля инерции
    double c1 = 1.49445; // когнитивная весовая доля
    double c2 = 1.49445; // социальная весовая доля

    double bestMinGlobalFitness = std::numeric_limits<double>::max();
    vector<double> bestMinGlobalPosition;

    double bestMaxGlobalFitness = std::numeric_limits<double>::min();
    vector<double> bestMaxGlobalPosition;

    vector<Particle> min_swarm, max_swarm;

    friend Particle;
};

ParticalSwarm* ParticalSwarm::instance = nullptr;

int main(int argc, char** argv) {
    // Задаём константы: количество итераций, количество частиц, число пространств и рамки для координаты minX - maxX
    ParticalSwarm swarm(1000, 1000, 2, -100, 100);
    swarm.calculate_max();
    swarm.calculate_min();

    pair<vector<double>, double> min_result = swarm.get_min_position_and_fitness(); 
    pair<vector<double>, double> max_result = swarm.get_max_position_and_fitness(); 

    std::cout << "Минимальная позиция и значение функции:\nТочка: ";
    for (int i = 0; i < min_result.first.size(); ++i) {
        std::cout << min_result.first[i] << " ";
    }

    std::cout << "\nЗначение: " << min_result.second << std::endl;

    std::cout << "Максимальная позиция и значение функции:\nТочка: ";
    for (int i = 0; i < max_result.first.size(); ++i) {
        std::cout << max_result.first[i] << " ";
    }

    std::cout << "\nЗначение: " << max_result.second << std::endl;

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(800, 600);
    glutCreateWindow("Fitness Function");

    glEnable(GL_DEPTH_TEST);
    ParticalSwarm::instance = &swarm;
    glutDisplayFunc(ParticalSwarm::make_graphic_static);
    glutMainLoop();
    
    return 0;
}
