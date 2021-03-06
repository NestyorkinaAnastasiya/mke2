/*grid.cpp*/
#include "grid.h"

namespace grid
{
	// ����� ����������� ������ �.�. � ��������
	bool Element::SearchNode(int nodesNumber)
	{
		for (int i = 0; i < 4; i++)
			if (nodesNumber == nodes[i]) return true;

		return false;
	}

	// �������� ����������� ������
	/*
	# inf2tr.dat (��������� ����):
	kuzlov - ����� ����� �����
	ktr    - ����� �������� ��������� (���������������)
	kt1    - ����� ����� � ������� �������� ���������
	# rz.dat (�������� ����):
	��������� i-� ������d (i=1..kuzlov):
	2*double (x,y),	��� x,y - (x,y)-���������� i-� �������
	# nvtr.dat (�������� ����):
	��������� i-� ������ (i=1..ktr):
	6*long (i1,i2,i3,i4,0,1), ��� i1,i2,i3,i4 - ���������� ������
	������ i-�� �������������� (����� �������, ������ �������,
	����� ������, ������ ������)
	# nvkat2d.dat (�������� ����):
	��������� i-� ������ (i=1..ktr):
	1*long (m), ��� m - ����� ��������� i-�� �������������� �
	������������ � ������� sreda � mu
	# l1.dat (�������� ����):
	��������� i-� ������ (i=1..kt1):
	1*long (k), ��� k - ���������� ����� i-� ������� � ������
	�������� ������� ��������
	*/

	// ���������� �����
	Grid::Grid()
	{
		int countOfNodes, countOfElements, countOfBC, l;
		FILE *fNodes, *fElements, *l1, *fp;

		/* ������� fseek ���������� ��������� ������� � ������.
		������������� ���������� ��������� ��������� � �����,
		� ����� �������, ������� ������������ ����� ����������
		�������� � ��������� ���������.

		int fseek( FILE * filestream, long int offset, int origin );
		- filestream	��������� �� ������ ���� FILE, ���������������� �����.

		- offset		���������� ���� ��� ��������, ������������ ���������� ��������� ���������.

		- origin		������� ���������, ������������ ������� ����� ����������� ��������.
		����� ������� ������� ����� �� ��������� ��������, ����������� �
		������������ ����� <cstdio>:

		SEEK_SET	������ �����
		SEEK_CUR	������� ��������� �����
		SEEK_END	����� �����*/


		fopen_s(&fp, "inf2tr.dat", "r");
		char c = ' ';
		fseek(fp, 56, SEEK_SET);
		fscanf(fp, "%d", &countOfNodes);
		while (c != '=')
			fscanf(fp, "%c", &c);
		fscanf(fp, "%d", &countOfElements);
		c = ' ';
		while (c != '=')
			fscanf(fp, "%c", &c);
		fscanf(fp, "%d", &countOfBC);
		fclose(fp);
		nodes.resize(countOfNodes);
		elements.resize(countOfElements);
		ku.resize(countOfBC);

		fopen_s(&fNodes, "rz.dat", "rb");
		for (int i = 0; i < countOfNodes; i++)
		{
			double xy[2];
			fread(xy, sizeof(double), 2, fNodes);
			nodes[i].x = xy[0];
			nodes[i].y = xy[1];
		}
		fclose(fNodes);

		fopen_s(&fElements, "nvtr.dat", "rb");
		fp = fopen("nvkat2d.dat", "rb");
		double g = 0;
		for (int i = 0; i < countOfElements; i++)
		{
			int nodes[6], material;
			fread(nodes, sizeof(int), 6, fElements);
			fread(&material, sizeof(int), 1, fp);


			if (material == 3 || material == 4)
				g++;

			// ��������� � ����� ���������� � 1
			for (int j = 0; j < 4; j++)
				elements[i].nodes[j] = nodes[j] - 1;

			sort(elements[i].nodes, elements[i].nodes + 4);
			elements[i].numberOfMaterial = material;
		}
		fclose(fp);
		fclose(fElements);

		// ���������� ����� � �������� ���������
		fopen_s(&l1, "l1.dat", "rb");
		for (int i = 0; i < countOfBC; i++)
		{
			fread(&l, sizeof(int), 1, l1);
			ku[i] = l - 1;
		}

		sort(ku.begin(), ku.end());
		fclose(l1);
	}

	Grid::~Grid() {}

	// ������������ ������ ������� ���������, ���������� ���������� ����� �.�.
	// ������ nodesNumber
	void Grid::SearchElements(int nodesNumber, vector <int> &elList)
	{
		int count;
		int size = elements.size();
		elList.reserve(4);

		count = 0;
		for (int i = 0; i < size && count < 4; i++)
		{
			if (elements[i].SearchNode(nodesNumber))
			{
				count++;
				elList.push_back(i);
			}
		}
	}
}