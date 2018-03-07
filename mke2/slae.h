/*slae.h*/
#pragma once
#include "tests.h"
#define _USE_MATH_DEFINES
#include <math.h>
using namespace matrix;
using namespace basis;
using namespace integration;
using namespace tests;

namespace slae
{
	class SLAE : private Basis, private GaussIntegration
	{
		// Размерность задачи
		int n;
		// Максимальное количество итераций в решателе
		int maxiter = 10000;
		// Точность решения СЛАУ
		const double eps = 1e-10;

		// Значение мю в вакууме(4 * пи * 10 ^ -7)
		double mu0 = 4 * M_PI * 0.0000001;
		double Jz = 1e7;
		bool liner = true;
		vector<pair<double, double>> tableMu;
		vector<double> tableDetMu;

		// Сетка
		Grid grid;
		// Хранилище тестовых функций
		Tests tests;
		// Глобальная матрица
		Matrix A;
		// Локальные матрицы
		// Матрица жёсткости
		array<array<double, 4>, 4> G;
		// Матрица массы
		array<array<double, 4>, 4> M;
		array<array<double, 4>, 4> Mg;
		// Локальный вектор правой части
		array <double, 4> locF;
		// Глобальный вектор правой части
		vector <double> F;
		array <double, 4> newF;
		// Вектор приближенного решения на предыдущей итерации
		// по нелинейности
		vector <double> u_n;
		// Вектор приближенного решения
		vector <double> u;
		// Норма вектора правой части
		double normF;

		// Сборка локальных матриц жёсткости
		void CalculateG(int elementNumber);
		// Сборка локальных правых частей
		void CalculateLocalF(int elementNumber);
		// Добавка локального элемента в глобальный
		void AddElementToGlobalMatrix(Matrix &B, int i, int j, double element);
		// Сборка локальных матриц(векторов) и добавление в глобальные
		void CalculateLocals(int elementNumber);

		int detectIntervalB(double B);

		double CalculateMu(double x, double y, int elementNumber);
		void CalculateDerivatives();
		double CalculateB(double x, double y, int elementNumber, double *Bx, double *By);
		double CalculateBInPoint(double x, double y, double *Bx, double *By);
		double SplineInterpolation(int interval_B_number, double B);

		// Вектор праввой части для первого краевого 
		array<double, 3> g;
		// Учёт первого краевого условия
		void CalculateBoundaries1();

		// Компоненты матрицы с факторизацией
		vector <double> L;
		vector <double> D;
		vector <double> U;

		// Вектор невязки
		vector <double> r;
		// Вектор спуска
		vector <double> z;

		// Вычисление нормы вектора
		double Norm(const vector<double>& x);
		// Скалярное произведение векторов
		double Scalar(const vector<double>& x, const vector<double>& y);

		// Генерация СЛАУ на i-ой итерации по времени
		void GenerateSLAE();
		// LU-факторизация
		void LU();
		// Вспомогательные функции для решателя
		void LYF(const vector<double>& C, vector<double>& yl);
		void UXY(const vector<double>& C, vector<double>& yu);
		double Rel_Discrepancy();
		// Решатель ЛОС с LU-факторизацией
		void LULOS();

	public:
		SLAE();

		void Solve();

		~SLAE() {};
	};
}
