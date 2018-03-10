/*matrix.cpp*/
#include "matrix.h"

namespace matrix
{
	// Генерация портрета матрицы
	void Matrix::CreatePortret(int slaeSize, Grid grid)
	{
		vector <int> elList, unzeroNumbersList;

		vector<vector <int>> list;
		list.reserve(slaeSize);

		for (int i = 0; i < slaeSize; i++)
		{
			// Нашли,в каких элементах присутствует данный i-ый узел
			grid.SearchElements(i, elList);
			// Собрать список узлов этих элементов, меньших по номеру
			for (unsigned int j = 0; j < elList.size(); j++)
			{
				for (int k = 0; k < 4; k++)
				{
					//если такого номера ещё нет и он меньше по номеру,то добавляем
					if (find(unzeroNumbersList.begin(), unzeroNumbersList.end(), grid.elements[elList[j]].nodes[k])
						== unzeroNumbersList.end() && grid.elements[elList[j]].nodes[k] < i)
						unzeroNumbersList.push_back(grid.elements[elList[j]].nodes[k]);
				}
			}

			sort(unzeroNumbersList.begin(), unzeroNumbersList.end());
			list.push_back(unzeroNumbersList);
			unzeroNumbersList.clear();
			elList.clear();
		}

		//вычисляем размерность ggl,ggu
		int gg_size = 0;
		for (int i = 0; i < slaeSize; i++)
		{
			if (!list[i].empty())
				gg_size += list[i].size();
		}

		//инициализируем матрицу и собираем портрет
		n = slaeSize; size = gg_size;

		ggl.resize(size);
		ggu.resize(size);
		di.resize(n);
		ig.resize(n + 1);
		jg.resize(size);

		ig[0] = 0;

		for (int i = 0; i < n; i++)
		{
			if (!list[i].empty())
				ig[i + 1] = ig[i] + list[i].size();
			else
				ig[i + 1] = ig[i];
		}

		int k = 0;
		for (int i = 0; i < n; i++)
		{
			if (!list[i].empty())
			{
				for (unsigned int j = 0; j < list[i].size(); j++)
				{
					jg[k] = list[i][j];
					k++;
				}
				list[i].clear();
			}
		}

		list.clear();
	}

	// Умножение матрицы на вектор
	void Matrix::MultiplyAx(const vector <double> a, vector <double> &result)
	{
		int i, j, l, ik, iend, k;

		for (i = 0; i < n; i++)
		{
			// начало i-ой строки(столбца)
			l = ig[i];
			// начало (i+1)-ой строки(столбца)
			iend = ig[i + 1];
			// количество элементов в i строке(столбце)
			ik = iend - l;

			result[i] = di[i] * a[i];

			// проходим по всем элементам i строки (столбца)
			for (k = 0; k < ik; k++, l++)
			{
				j = jg[l];
				result[i] += ggl[l] * a[j];
				result[j] += ggu[l] * a[i];
			}
		}
	}
	void Matrix::MultiplyATx(vector <double> a, vector <double> &result)
	{
		int i, j, l, ik, iend, k;
		for (i = 0; i < n; i++)
		{
			//начало i-ой строки(столбца)
			l = ig[i];
			//начало (i+1)-ой строки(столбца)
			iend = ig[i + 1];
			//количество элементов в i строке(столбце)
			ik = iend - l;

			result[i] = di[i] * a[i];

			//проходим по всем элементам i строки (столбца)
			for (k = 0; k < ik; k++, l++)
			{
				j = jg[l];
				result[i] += ggu[l] * a[j];
				result[j] += ggl[l] * a[i];
			}
		}
	}
}

