/*matrix.h*/
#pragma once
#include <algorithm>
#include "basis.h"

namespace matrix
{
	struct Matrix
	{
		//Матрица СЛАУ
		//Указатели начала строк(столбцов)
		vector <int> ig;
		//Номера столбцов внедиагональных элементов 
		vector <int> jg;
		//Диагональные элементы матрицы
		vector <double> di;
		//Внедиагональные элементы нижнего треугольника матрицы
		vector <double> ggl;
		//Внедиагональные элементы верхнего треугольника матрицы
		vector <double> ggu;

		//Размерность матрицы
		int n;
		//Размер векторов ggl, ggu, jg
		int size;

		//Генерация портрета матрицы
		void CreatePortret(int slaeSize, Grid grid);
		//Умножение матрицы на вектор
		void MultiplyAx(const vector <double> a, vector<double>& result);
		void MultiplyATx(vector<double> a, vector<double>& result);
		void Sum(vector<double>& result);
	};
}
