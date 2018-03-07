/*grid.h*/
#define _CRT_SECURE_NO_WARNINGS
#pragma once
#include <stdio.h>
#include <vector>
#include <fstream>
#include <functional>
#include <algorithm>
#include <array>
using namespace std;

namespace grid
{
	//Конечный элемент
	struct Element
	{
		//Узлы
		int nodes[4];
		int numberOfMaterial;

		Element& operator=(Element element)
		{
			for (int i = 0; i < 4; i++)
				nodes[i] = element.nodes[i];
			numberOfMaterial = element.numberOfMaterial;

			return *this;
		}
		//Поиск глобального номера nodesNumber в элементе
		bool SearchNode(int nodesNumber);
	};

	struct Point
	{
		double x;
		double y;
		Point() {};
		~Point() {};
		Point(double xx, double yy)
		{
			x = xx;
			y = yy;
		}

		bool operator==(Point point)
		{
			if (point.x == x && point.y == y)
				return true;
			else
				return false;
		}

	};

	//Разбиение общей области
	struct Grid
	{
		Grid();
		//Массив конечных элементов
		vector <Element> elements;
		//Массив узлов
		vector <Point> nodes;
		//Массив узлов с первыми краевыми условиями
		vector <int> ku;
		//Формирование списка элементов, содержащих глобальный номер б.ф.
		//равный nodesNumber
		void SearchElements(int nodesNumber, vector<int>& elList);

		~Grid();
	};
}