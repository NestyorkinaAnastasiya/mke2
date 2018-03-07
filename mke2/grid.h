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
	//�������� �������
	struct Element
	{
		//����
		int nodes[4];
		int numberOfMaterial;

		Element& operator=(Element element)
		{
			for (int i = 0; i < 4; i++)
				nodes[i] = element.nodes[i];
			numberOfMaterial = element.numberOfMaterial;

			return *this;
		}
		//����� ����������� ������ nodesNumber � ��������
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

	//��������� ����� �������
	struct Grid
	{
		Grid();
		//������ �������� ���������
		vector <Element> elements;
		//������ �����
		vector <Point> nodes;
		//������ ����� � ������� �������� ���������
		vector <int> ku;
		//������������ ������ ���������, ���������� ���������� ����� �.�.
		//������ nodesNumber
		void SearchElements(int nodesNumber, vector<int>& elList);

		~Grid();
	};
}