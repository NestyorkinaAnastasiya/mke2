/*grid.cpp*/
#include "grid.h"

namespace grid
{
	// Поиск глобального номера б.ф. в элементе
	bool Element::SearchNode(int nodesNumber)
	{
		for (int i = 0; i < 4; i++)
			if (nodesNumber == nodes[i]) return true;

		return false;
	}

	// Описание содержимого файлов
	/*
	# inf2tr.dat (текстовый файл):
		kuzlov - число узлов сетки
		ktr    - число конечных элементов (прямоугольников)
		kt1    - число узлов с первыми краевыми условиями
	# rz.dat (двоичный файл):
		структура i-й записи (i=1..kuzlov):
			2*double (x,y),	где x,y - (x,y)-координаты i-й вершины
	# nvtr.dat (двоичный файл):
		структура i-й записи (i=1..ktr): 
			6*long (i1,i2,i3,i4,0,1), где i1,i2,i3,i4 - глобальные номера 
			вершин i-го прямоугольника (левый верхний, правый верхний,
			левый нижний, правый нижний)
	# nvkat2d.dat (двоичный файл):
		структура i-й записи (i=1..ktr):
			1*long (m), где m - номер материала i-го прямоугольника в 
			соответствии с файлами sreda и mu
	# l1.dat (двоичный файл):
		структура i-й записи (i=1..kt1):
			1*long (k), где k - глобальный номер i-й вершины с первым
			нулевыми краевым условием
	*/
	
	// Заполнение сетки
	Grid::Grid()
	{
		int countOfNodes, countOfElements, l;
		FILE *fNodes, *fElements, *l1, *fp ;
		
		fNodes = fopen("rz.dat", "rb");
		fElements = fopen("nvtr.dat", "rb");
		fp = fopen("inf2tr.dat", "r");		
		l1 = fopen("l1.dat", "rb");

		/* Функция fseek перемещает указатель позиции в потоке. 
		Устанавливает внутренний указатель положения в файле, 
		в новую позицию, которая определяются путем добавления 
		смещения к исходному положению.

		int fseek( FILE * filestream, long int offset, int origin );
		- filestream	Указатель на объект типа FILE, идентифицируемый поток.

		- offset		Количество байт для смещения, относительно некоторого положения указателя.

		- origin		Позиция указателя, относительно которой будет выполняться смещение. 
						Такая позиция задаётся одной из следующих констант, определённых в
						заголовочном файле <cstdio>:

		SEEK_SET	Начало файла
		SEEK_CUR	Текущее положение файла
		SEEK_END	Конец файла*/

		fseek(fp, 56, SEEK_SET);
		fscanf(fp, "%d", &countOfNodes);
		fseek(fp, 8, SEEK_CUR);
		fscanf(fp, "%d", &countOfElements);
		fclose(fp);

		nodes.resize(countOfNodes);
		elements.resize(countOfElements);

		for (int i = 0; i < countOfNodes; i++)
		{
			double xy[2];
			fread(xy, sizeof(double), 2, fNodes);
			nodes[i].x = xy[0];
			nodes[i].y = xy[1];
		}

		fp = fopen("nvkat2d.dat", "rb");
		for (int i = 0; i < countOfElements; i++)
		{
			int nodes[6], material;
			fread(nodes, sizeof(int), 6, fElements);
			fread(&material, sizeof(int), 1, fp);

			// Нумерация в файле начинается с 1
			for (int j = 0; j < 4; j++)
				elements[i].nodes[j] = nodes[j] - 1;

			sort(elements[i].nodes, elements[i].nodes + 4);
			elements[i].numberOfMaterial = material;
		}
		fclose(fp);

		// Заполнение узлов с краевыми условиями
		while (fread(&l, sizeof(int), 1, l1) == 1)
		{
			ku.push_back(l - 1);
		}
	}

	Grid::~Grid() {}
	
	// Формирование списка элементов, содержащих глобальный номер б.ф.
	// равный nodesNumber
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