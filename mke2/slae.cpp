/*slae.cpp*/
#include "slae.h"

namespace slae
{
	/*
	# Файл mu.001 (зависимость mu с крышечкой от вектора индукции B магнитного поля) содержит в себе:
	1. общее количество мюшек
	2. пары чисел: Bx и By, наверное
	Общее mu высчитывается следующим образом: mu = mu0 * mu с крышечкой
	*/

	SLAE::SLAE()
	{
		// Размерность задачи соответствует общему числу базисных функций
		n = grid.nodes.size();

		F.resize(n);
		u.resize(n);
		for (int i = 0; i < n; i++)
			u[i] = 10;
		u_n.resize(n);
		r.resize(n);
		z.resize(n);
		// Генерация портрета матрицы и её инициализация
		A.CreatePortret(n, grid);

		if (!liner)
		{
			// Считываем табличные значения мю
			FILE *fp;
			fopen_s(&fp, "mu.001", "r");
			int size;
			fscanf_s(fp, "%d", &size);
			tableMu.resize(size);
			for (int i = 0; i < size; i++)
			{
				fscanf_s(fp, "%lf %lf", &tableMu[i].second, &tableMu[i].first);
				tableMu[i].second *= mu0;
			}
			fclose(fp);

			// Вычисление значений производных полиномов Лагранжа в точках 
			// (tableDetMu)
			CalculateDerivatives();
		}
	}

#pragma region matrix
	// Сборка локальных матриц жёсткости
	void SLAE::CalculateG(int elementNumber)
	{
		Element element = grid.elements[elementNumber];
		double g1, g2, ksi, etta, x_, y_, lambda,
			x1 = grid.nodes[element.nodes[0]].x, x3 = grid.nodes[element.nodes[1]].x,
			y1 = grid.nodes[element.nodes[0]].y, y3 = grid.nodes[element.nodes[2]].y,
			hx = x3 - x1, hy = y3 - y1,
			hx2 = hx * hx, hy2 = hy * hy,
			jacobian = hx * hy / 4.0;

		for (int i = 0; i < 4; i++)
		{
			for (int j = i; j < 4; j++)
			{
				g1 = 0; g2 = 0;
				for (int k = 0; k < 25; k++)
				{
					ksi = 0.5 + 0.5 * gaussPoints[0][k]; etta = 0.5 + 0.5 * gaussPoints[1][k];
					x_ = x1 + ksi * hx; y_ = y1 + etta * hy;

					// Номер материала соответствует железу
					if (element.numberOfMaterial == 2)
					{
						if (!liner)
							lambda = 1 / CalculateMu(x_, y_, elementNumber);
						else lambda = 1 / (1000 * mu0);
					}
					else // В остальных случаях нет зависмости от mu
					{
						lambda = 1.0 / mu0;
					}
					g1 += gaussWeights[k] * dphiksi[i](ksi, etta) * dphiksi[j](ksi, etta) * lambda;
					g2 += gaussWeights[k] * dphietta[i](ksi, etta) * dphietta[j](ksi, etta) * lambda;
				}
				G[i][j] = g1 * jacobian / hx2 + g2 * jacobian / hy2;
			}
		}
		// матрица симметричная, заполняем нижний треугольник
		for (int i = 1; i < 4; i++)
			for (int j = 0; j < i; j++)
				G[i][j] = G[j][i];

	}

	// Сборка локальных правых частей
	void SLAE::CalculateLocalF(int elementNumber)
	{
		Element element = grid.elements[elementNumber];
		double ksi, etta, x_, y_,
			x1 = grid.nodes[element.nodes[0]].x, x3 = grid.nodes[element.nodes[1]].x,
			y1 = grid.nodes[element.nodes[0]].y, y3 = grid.nodes[element.nodes[2]].y,
			hx = x3 - x1, hy = y3 - y1,
			jacobian = hx * hy / 4.0;

		// интегрирование(Гаусс 5)
		for (int i = 0; i < 4; i++)
		{
			locF[i] = 0;
			if (element.numberOfMaterial == 3 || element.numberOfMaterial == 4)
			{
				for (int k = 0; k < 25; k++)
				{
					ksi = 0.5 + 0.5 * gaussPoints[0][k]; etta = 0.5 + 0.5 * gaussPoints[1][k];
					x_ = x1 + ksi * hx; y_ = y1 + etta * hy;
					if (element.numberOfMaterial == 3)
						locF[i] += Jz * gaussWeights[k] * phi[i](ksi, etta);
					else if (element.numberOfMaterial == 4)
						locF[i] -= Jz * gaussWeights[k] * phi[i](ksi, etta);
				}
				locF[i] *= jacobian;
			}
		}
		/*double fi;
		Element element = grid.elements[elementNumber];
		double x1 = grid.nodes[element.nodes[0]].x, x3 = grid.nodes[element.nodes[1]].x;
		double y1 = grid.nodes[element.nodes[0]].y, y3 = grid.nodes[element.nodes[2]].y;

		double hx = x3 - x1;
		double hy = y3 - y1;

		double f1, f2, f3, f4;
		if (element.numberOfMaterial == 3)
			f1 = f2 = f3 = f4 = Jz;
		else
		{
			if (element.numberOfMaterial == 4)
				f1 = f2 = f3 = f4 = -Jz;
			else f1 = f2 = f3 = f4 = 0;
		}

		locF[0] = hy * hx*(4 * f1 + 2 * f2 + 2 * f3 + f4) / 36;
		locF[1] = hy * hx*(2 * f1 + 4 * f2 + f3 + 2 * f4) / 36;
		locF[2] = hy * hx*(2 * f1 + f2 + 4 * f3 + 2 * f4) / 36;
		locF[3] = hy * hx*(f1 + 2 * f2 + 2 * f3 + 4 * f4) / 36;*/
	}

	// Добавка локального элемента в глобальный
	void SLAE::AddElementToGlobalMatrix(Matrix &B, int i, int j, double element)
	{
		int id;
		bool flag;

		if (i == j)
			B.di[i] += element;
		else
		{
			if (i < j)
			{
				flag = false;
				for (id = B.ig[j]; !flag && id < B.ig[j + 1]; id++)
					if (B.jg[id] == i) flag = true;

				if (flag) B.ggu[id - 1] += element;
			}
			else
			{
				flag = false;
				for (id = B.ig[i]; !flag && id < B.ig[i + 1]; id++)
					if (B.jg[id] == j) flag = true;

				if (flag) B.ggl[id - 1] += element;
			}
		}
	}

	// Сборка локальных матриц(векторов) и добавление в глобальные
	void SLAE::CalculateLocals(int elementNumber)
	{
		Element element = grid.elements[elementNumber];
		int ki, kj;

		// вычисление локальных матриц
		CalculateG(elementNumber);
		CalculateLocalF(elementNumber);

		for (int i = 0; i < 4; i++)
		{
			ki = element.nodes[i];
			for (int j = 0; j < 4; j++)
			{
				kj = element.nodes[j];
				// добавка в глобальную матрицу А
				AddElementToGlobalMatrix(A, ki, kj, G[i][j]);

			}
			// добавка в глобальную правую часть
			F[ki] = locF[i];
		}
	}

	// Вычисление 1ого краевого условия для одного узла
	void SLAE::CalculateBoundaries1()
	{
		int id;
		int nNodes = grid.ku.size();
		int matrixSize = A.n;
		int numberOfNode;
		bool flag;

		for (int i = 0; i < nNodes; i++)
		{
			numberOfNode = grid.ku[i];

			F[numberOfNode] = 0.0;
			A.di[numberOfNode] = 1e307;
			//u[numberOfNode] = 0;

			/*	for (int j = 0; j < matrixSize; j++)
			if (numberOfNode < j)
			{
			flag = false;
			for (id = A.ig[j]; !flag && id <= A.ig[j + 1] - 1; id++)
			if (A.jg[id] == numberOfNode) flag = true;
			if (flag) A.ggu[id - 1] = 0.0;
			}
			else
			{
			flag = false;
			for (id = A.ig[numberOfNode]; !flag && id <= A.ig[numberOfNode + 1] - 1; id++)
			if (A.jg[id] == j) flag = true;
			if (flag) A.ggl[id - 1] = 0.0;
			}*/
		}
	}
#pragma endregion

#pragma region B, mu
	// Определение интервала в котором находится B
	int SLAE::detectIntervalB(double B)
	{
		int i;
		int size = tableMu.size();
		for (i = 0; i < size; i++)
			if (B < tableMu[i].first)
				// если B меньше минимально возможного,то будет возвращаться -1
				// если B будет меньше какого-либо внутреннего i-ого B из таблицы,
				//	то будет возвращаться номер B из таблицы, который меньше рассматриваемого B
				return i - 1;
		// если B больше максимально возможного, то будет возвращаться -2 
		return -2;
	}

	// Вычисление мю в зависимости от B
	double SLAE::CalculateMu(double x, double y, int elementNumber)
	{
		double Bx, By;
		double B = CalculateB(x, y, elementNumber, &Bx, &By);
		int intervalBNumber = detectIntervalB(B);
		// B меньше минимально возможного, берём минимально возможный в качестве значения
		if (intervalBNumber == -1)
		{
			return tableMu[0].second;
		}
		else // B больше максимально возможного
			if (intervalBNumber == -2)
			{
				// Вычисление mu по формуле (2.12) из методы
				return 1.0 / ((tableMu.end() - 1)->first / B * (1.0 / (tableMu.end() - 1)->second - 1) + 1);
			}
		// если B внутренняя точка - построение сплайна
		return SplineInterpolation(intervalBNumber, B);
	}

	void SLAE::CalculateDerivatives()
	{
		// Шаги, i-тый и (i-1)-ый
		double hi, hi1;
		double temp1, temp2, temp3;
		int size = tableMu.size();
		tableDetMu.resize(size);
		// Формирование fi' = dPi^L(x)/dx, f1' и fn' считаются отдельно
		// Формула (4.91)
		for (int i = 1; i < size - 2; i++)
		{
			hi = tableMu[i + 1].first - tableMu[i].first; // hi
			hi1 = tableMu[i].first - tableMu[i - 1].first; // hi-1

			temp1 = -hi / (hi1 * (hi1 + hi));
			temp1 *= tableMu[i - 1].second;

			temp2 = (hi - hi1) / (hi1 * hi);
			temp2 *= tableMu[i].second;

			temp3 = hi1 / (hi * (hi1 + hi));
			temp3 *= tableMu[i + 1].second;

			tableDetMu[i] = temp1 + temp2 + temp3;
		}

		// Формирование f1' = dP2^L(x)/dx
		// Формула (4.92)
		hi1 = tableMu[1].first - tableMu[0].first;	// h1
		hi = tableMu[2].first - tableMu[1].first;		// h2

		temp1 = -(2 * hi1 + hi) / (hi1 * (hi1 + hi));
		temp1 *= tableMu[0].second;

		temp2 = (hi1 + hi) / (hi1 * hi);
		temp2 *= tableMu[1].second;

		temp3 = -hi1 / (hi * (hi1 + hi));
		temp3 *= tableMu[2].second;

		tableDetMu[0] = temp1 + temp2 + temp3;

		// Формирование fn' = dPn-1^L(x)/dx
		// Формула (4.93)
		hi = tableMu[size - 1].first - tableMu[size - 2].first; // hn-1
		hi1 = tableMu[size - 2].first - tableMu[size - 3].first; // hn-2

		temp1 = hi / (hi1 * (hi1 + hi));
		temp1 *= tableMu[size - 3].second;

		temp2 = -(hi1 + hi) / (hi1 * hi);
		temp2 *= tableMu[size - 2].second;

		temp3 = (2 * hi + hi1) / (hi * (hi1 + hi));
		temp3 *= tableMu[size - 1].second;

		tableDetMu[size - 1] = temp1 + temp2 + temp3;
	}

	// Вычисление B по формуле B = sqrt((dAz/dy)^2+(dAz/dx))
	double SLAE::CalculateB(double x, double y, int elementNumber, double *Bx, double *By)
	{
		Element element = grid.elements[elementNumber];

		double qi[4], dAx, dAy, hx, hy, dksi, detta, ksi, etta;

		hx = grid.nodes[element.nodes[1]].x - grid.nodes[element.nodes[0]].x;
		hy = grid.nodes[element.nodes[2]].y - grid.nodes[element.nodes[0]].y;
		ksi = (x - grid.nodes[element.nodes[0]].x) / hx;
		etta = (y - grid.nodes[element.nodes[0]].y) / hy;

		//заполняем соответствующими весами из предыдущего решения по нелинейности
		if (liner)
			for (int i = 0; i < 4; i++)
				qi[i] = u[element.nodes[i]];
		else
			for (int i = 0; i < 4; i++)
				qi[i] = u_n[element.nodes[i]];
		dAx = 0;
		dAy = 0;
		for (int i = 0; i < 4; i++)
		{
			dAx += dphiksi[i](ksi, etta) * qi[i] / hx;
			dAy += dphietta[i](ksi, etta) * qi[i] / hy;
		}
		*Bx = dAy;
		*By = -dAx;

		//return sqrt(B2 / S);
		return sqrt(dAx * dAx + dAy * dAy);
	}

	// Определение интервала в котором находится B
	double SLAE::CalculateAzInPoint(double x, double y)
	{
		int size = grid.elements.size();
		int numberOfElement = -1;
		// Пока не прошли все элементы и не нашли элемента, в котором лежит данная точка
		for (int i = 0; i < size && numberOfElement == -1; i++)
			if (grid.nodes[grid.elements[i].nodes[0]].x <= x &&
				x <= grid.nodes[grid.elements[i].nodes[1]].x &&
				grid.nodes[grid.elements[i].nodes[0]].y <= y &&
				y <= grid.nodes[grid.elements[i].nodes[2]].y)
				numberOfElement = i;
		if (numberOfElement != -1)
			return CalculateAz(x, y, numberOfElement);
		else return -1;
	}

	double SLAE::CalculateAz(double x, double y, int elementNumber)
	{
		Element element = grid.elements[elementNumber];

		double qi[4], Az, hx, hy, dksi, detta, ksi, etta;

		hx = grid.nodes[element.nodes[1]].x - grid.nodes[element.nodes[0]].x;
		hy = grid.nodes[element.nodes[2]].y - grid.nodes[element.nodes[0]].y;
		ksi = (x - grid.nodes[element.nodes[0]].x) / hx;
		etta = (y - grid.nodes[element.nodes[0]].y) / hy;

		//заполняем соответствующими весами из предыдущего решения по нелинейности
		if (liner)
			for (int i = 0; i < 4; i++)
				qi[i] = u[element.nodes[i]];
		else
			for (int i = 0; i < 4; i++)
				qi[i] = u_n[element.nodes[i]];
		Az = 0;
		for (int i = 0; i < 4; i++)
		{
			Az += qi[i] * phi[i](ksi, etta);
		}

		return Az;
	}

	// Получение значения B в точке (x,y)
	double SLAE::CalculateBInPoint(double x, double y, double *Bx, double *By)
	{
		int size = grid.elements.size();
		int numberOfElement = -1;
		// Пока не прошли все элементы и не нашли элемента, в котором лежит данная точка
		for (int i = 0; i < size && numberOfElement == -1; i++)
			if (grid.nodes[grid.elements[i].nodes[0]].x <= x &&
				x <= grid.nodes[grid.elements[i].nodes[1]].x &&
				grid.nodes[grid.elements[i].nodes[0]].y <= y &&
				y <= grid.nodes[grid.elements[i].nodes[2]].y)
				numberOfElement = i;
		if (numberOfElement != -1)
			return CalculateB(x, y, numberOfElement, Bx, By);
		else return -1;
	}

	// Построение кубического интерполяционного сплайна
	double SLAE::SplineInterpolation(int interval_B_number, double B)
	{
		double psi1, psi2, psi3, psi4, h, ksi;

		// Билинейные эрмитовы элементы третьего порядка
		// Кирпич страница 151, формула (4.33)
		h = (tableMu[interval_B_number + 1].first - tableMu[interval_B_number].first);	// hk
		ksi = (B - tableMu[interval_B_number].first) / h;								// ksi = (x-xk)/hk, xk = tableMu[interval_B_number].first

		psi1 = (1 - 3 * ksi * ksi + 2 * ksi * ksi * ksi);
		// Для того, чтобы локальная базисная функция имела равную единице производную, 
		// домножаем на шаг
		psi2 = (h * (ksi - 2 * ksi * ksi + ksi * ksi * ksi));
		psi3 = (3 * ksi * ksi - 2 * ksi * ksi * ksi);
		psi4 = (h * (-ksi * ksi + ksi * ksi * ksi));

		// Кирпич, стр 204, формула (4.90)
		return tableMu[interval_B_number].second * psi1 +
			tableDetMu[interval_B_number] * psi2 +
			tableMu[interval_B_number + 1].second * psi3 +
			tableDetMu[interval_B_number + 1] * psi4;
	}
#pragma endregion
	// Генерация СЛАУ
	void SLAE::GenerateSLAE()
	{
		// Высчитывание локальных матриц(векторов) и добавление в глобальные
		for (int i = 0; i < grid.elements.size(); i++)
			CalculateLocals(i);

		CalculateBoundaries1();

		normF = Norm(F);
		/*F.resize(3);
		F = { 9,3,5 };
		A.di.resize(3);
		A.di = { 1,1,1 };
		A.n = 3;
		n = 3;
		A.ig.resize(4);
		A.ig = { 0,0,1,2 };
		A.jg.resize(2);
		A.jg = {0,0};
		A.ggl.resize(2);
		A.ggu.resize(2);
		A.ggl = A.ggu = { 2,3 };
		u.resize(3);
		r.resize(3);
		z.resize(3);*/

	}


#pragma region LU_LOS
	// Вычисление нормы вектора
	double SLAE::Norm(const vector<double> &x)
	{
		double norm = 0;
		int size = x.size();

		for (int i = 0; i < size; i++)
			norm += x[i] * x[i];

		return sqrt(norm);
	}

	// Скалярное произведение векторов
	double SLAE::Scalar(const vector<double> &x, const vector<double> &y)
	{
		double sum = 0;
		int size = x.size();
		for (int i = 0; i < size; i++)
			sum += x[i] * y[i];

		return sum;
	}

	// LU-факторизация
	void SLAE::LU()
	{
		int i, i0, j0, iend, num, ki, kj, jend;
		double suml, sumu, sumdg;

		L.resize(A.ggl.size());
		L = A.ggl;
		U.resize(A.ggu.size());
		U = A.ggu;
		D.resize(A.di.size());
		D = A.di;

		for (i = 0; i < n; i++)
		{
			i0 = A.ig[i];
			iend = A.ig[i + 1];

			for (num = i0, sumdg = 0; num < iend; num++)
			{
				j0 = A.ig[A.jg[num]];
				jend = A.ig[A.jg[num] + 1];
				ki = i0;
				kj = j0;
				// для num учитываются все предыдущие элементы
				for (suml = 0, sumu = 0, ki = i0; ki < num; ki++)
				{
					for (int m = kj; m < jend; m++)
						// ищем соответствующие ненулевые элементы для умножения
						if (A.jg[ki] == A.jg[m])
						{
							suml += L[ki] * U[m];
							sumu += L[m] * U[ki];
						}
				}
				L[num] -= suml;
				U[num] = (U[num] - sumu) / D[A.jg[num]];
				// умножаются симметричные элементы
				sumdg += L[num] * U[num];
			}
			D[i] -= sumdg;
		}
	}

	void SLAE::LYF(const vector <double> &C, vector <double> &yl)
	{
		int i, i0, iend; // i0-адрес начала строки, iend-адрес конца строки
		double sum;
		for (i = 0; i < n; i++)
		{
			i0 = A.ig[i]; iend = A.ig[i + 1];

			yl[i] = C[i];

			for (i0, sum = 0; i0 < iend; i0++)
				yl[i] -= yl[A.jg[i0]] * L[i0];
			yl[i] /= D[i];
		}
	}

	void SLAE::UXY(const vector <double> &C, vector <double> &yu)
	{
		int i, i0, iend;

		for (i = 0; i < n; i++)
			yu[i] = 0.0;

		for (i = n - 1; i >= 0; i--)// проход по столбцам с конца
		{
			yu[i] += C[i];

			i0 = A.ig[i]; iend = A.ig[i + 1]; iend--;

			for (; iend >= i0; iend--)// идём по столбцу с конца
				yu[A.jg[iend]] -= yu[i] * U[iend];
		}
	}

	double SLAE::Rel_Discrepancy()
	{
		double dis1, dis2;
		dis1 = Scalar(r, r);
		dis2 = Scalar(F, F);
		double dis = dis1 / dis2;
		return sqrt(dis);
	}

	void SLAE::LULOS()
	{
		double a, b, pp, dis, rr;
		int i, k;
		vector <double> Ax(n), C(n), p(n);

		LU();
		// Ax0
		A.MultiplyAx(u, Ax);
		// f-Ax0
		for (i = 0; i < n; i++)
			r[i] = F[i] - Ax[i];

		// r0=L^(-1)(f-Ax0)
		LYF(r, r);

		// z0=U^(-1)r0->r0=Uz0
		UXY(r, z);

		// p0=L^(-1)Az0
		A.MultiplyAx(z, Ax);//Az0
		LYF(Ax, p);

		rr = Scalar(r, r);
		dis = Scalar(r, r) / rr;
		dis = sqrt(dis);
		k = 0;

		for (k = 1; dis > eps && k <= 10000; k++)
		{
			// Аk
			pp = Scalar(p, p);
			a = Scalar(p, r) / pp;

			// Xk, Rk
			for (i = 0; i < n; i++)
			{
				u[i] = u[i] + a * z[i];
				r[i] = r[i] - a * p[i];
			}

			// UY=rk->Y=U^(-1)rk
			UXY(r, C);
			// AU^(-1)rk=Ax
			A.MultiplyAx(C, Ax);
			// L^(-1)AU^(-1)rk=Y2->L^(-1)B=Y2->LY2=B->Y2=L^(-1)AU^(-1)rk
			LYF(Ax, Ax);

			// bk
			b = -Scalar(p, Ax) / pp;

			// zk=U^(-1)rk+bkz[k-1]
			// pk
			for (i = 0; i < n; i++)
			{
				z[i] = C[i] + b * z[i];
				p[i] = Ax[i] + b * p[i];
			}
			dis = Scalar(r, r) / rr;
			dis = sqrt(dis);
		}

		dis = Rel_Discrepancy();
	}
#pragma endregion


	void SLAE::Solve()
	{
		GenerateSLAE();
		LU_MSG();

		FILE *fo;
		fopen_s(&fo, "result.txt", "w");
		fprintf(fo, "Az:\n");
		for (int i = 0; i < n; i++)
			fprintf(fo, "%.14lf\n", u[i]);
		fprintf(fo, "\n\n");
		FILE *fp;
		fopen_s(&fp, "points.txt", "r");
		int n_;
		fscanf_s(fp, "%d", &n);
		double x, y, Bx, By, B, Az;
		for (int i = 0; i < n; i++)
		{
			fscanf_s(fp, "%lf%lf", &x, &y);
			B = CalculateBInPoint(x, y, &Bx, &By);
			Az = CalculateAzInPoint(x, y);
			fprintf(fo, "X= %.3le  Y= %.3le  Bx= %.4le\tBy= %.4le\t||B||= %.4le\tAz =  %lf\n", x, y, Bx, By, B, Az);
		}
		fclose(fp);
		fclose(fo);
	}
	void SLAE::LFx(const vector<double> &b, vector <double> &result)
	{
		int i, k;
		int i0;//адрес начала строки
		int iend;//адрес конца строки
		double sum;
		for (i = 0; i < n; i++)
		{
			i0 = A.ig[i];
			iend = A.ig[i + 1];

			result[i] = b[i];

			for (k = i0, sum = 0; k < iend; k++)
				result[i] -= result[A.jg[k]] * L[k];

			result[i] /= D[i];
		}
	}

	void SLAE::LTFx(const vector<double> &b, vector <double> &result)
	{
		int i, k;
		int i0;
		int iend;
		double sum;

		for (i = 0; i < n; i++)
			result[i] = b[i];

		for (i = n - 1; i >= 0; i--)
		{

			i0 = A.ig[i];
			iend = A.ig[i + 1];

			result[i] /= D[i];

			for (k = iend - 1; k >= i0; k--)
				result[A.jg[k]] -= result[i] * L[k];
		}
	}

	void SLAE::UFx(const vector<double> &b, vector <double> &result)
	{
		int i, k;
		int i0;
		int iend;

		for (i = 0; i < n; i++)
			result[i] = b[i];

		for (i = n - 1; i >= 0; i--)
		{
			i0 = A.ig[i];
			iend = A.ig[i + 1];

			for (k = iend - 1; k >= i0; k--)
				result[A.jg[k]] -= result[i] * U[k];
		}
	}

	void SLAE::UTFx(const vector<double> &b, vector <double> &result)
	{
		int i, k;
		int i0;
		int iend;
		for (i = 0; i < n; i++)
		{
			i0 = A.ig[i];
			iend = A.ig[i + 1];

			result[i] = b[i];

			for (k = i0; k < iend; k++)
				result[i] -= result[A.jg[k]] * U[k];
		}
	}

	double SLAE::IterMSG(const vector<double> &Az)
	{
		double ak, bk, scr;

		scr = Scalar(r, r);
		ak = scr / Scalar(Az, z);

		for (int i = 0; i < n; i++)
		{
			r[i] -= ak * Az[i];
			u[i] += ak * z[i];
		}
		bk = Scalar(r, r) / scr;

		for (int i = 0; i < n; i++)
			z[i] = r[i] + bk * z[i];

		return Norm(r) / normF;
	}

	void SLAE::MultiplyUx(const vector<double> &a, vector <double> &result)
	{
		int i, j, l, ik, iend, k;
		for (i = 0; i < n; i++)
		{
			//начало i-ой строки(столбца)
			l = A.ig[i];
			//начало (i+1)-ой строки(столбца)
			iend = A.ig[i + 1];
			//количество элементов в i строке(столбце)
			ik = iend - l;

			result[i] =/* A.di[i] **/ a[i];

			//проходим по всем элементам i строки (столбца)
			for (k = 0; k < ik; k++, l++)
			{
				j = A.jg[l];
				result[j] += U[l] * a[i];
			}
		}
	}
	void SLAE::LU_MSG()
	{
		int iter;
		double res;
		vector <double> Az;
		Az.resize(n);
		
		for (int i = 0; i < n; i++)
			u[i] = 0;

		LU();

		A.MultiplyAx(u, r);
		LFx(r, r);
		LTFx(r, r);
		A.MultiplyATx(r, r);
		UTFx(r, r);

		LFx(F, z);
		LTFx(z, z);
		A.MultiplyATx(z, z);
		UTFx(z, z);

		normF = Norm(z);
		for (int i = 0; i < n; i++)
			r[i] = z[i] - r[i];

		z = r;
		MultiplyUx(u, u);

		res = Norm(r) / normF;
		for (iter = 0; res >= eps && iter < 10000; iter++)
		{
			UFx(z, Az);
			A.MultiplyAx(Az, Az);
			LFx(Az, Az);
			LTFx(Az, Az);
			A.MultiplyATx(Az, Az);
			UTFx(Az, Az);

			res = IterMSG(Az);
		}
		UFx(u, u);
		printf("res = %e\t iter = %d", res, iter);
		getchar();

		/*for (int i = 0; i < n; i++)
			fprintf(fo, "%.14lf\n", x[i]);
		fprintf(fo, "%E\t", res);
		fprintf(fo, "%d\t", iter);*/
	}
	void SLAE::MSG()
	{
		int iter;
		double res;
		vector <double> Az;
		Az.resize(n);
		for (int i = 0; i < n; i++)
			u[i] = 1;

		A.MultiplyAx(u, r);
		A.MultiplyATx(r, r);

		A.MultiplyATx(F, z);
		normF = Norm(z);

		for (int i = 0; i < n; i++)
			r[i] = z[i] - r[i];

		z = r;
		res = Norm(r) / normF;

		for (iter = 0; res >= eps && iter < maxiter; iter++)
		{
			A.MultiplyAx(z, Az);
			A.MultiplyATx(Az, Az);

			res = IterMSG(Az);
		}

		/*for (int i = 0; i < n; i++)
			fprintf(fo, "%.14lf\n", x[i]);
		fprintf(fo, "%E\t", res);
		fprintf(fo, "%d\t", iter);*/
	}

#pragma region MSG Sim

	//вычисление значения вектора
	void SLAE::CulcVect(vector<double>& res, const vector<double>&v1, double var, vector<double>& v2)
	{
		for (int i = 0; i < n; i++)
			res[i] = v1[i] + var * v2[i];
	}

	//равенство векторов
	void SLAE::Equal(vector<double>& vo, vector<double>& vi)
	{
		// vo = vi;
		for (int i = 0; i < n; i++)
			vo[i] = vi[i];
	}

	//невязка
	double SLAE::Residual()
	{
		double resid = 0.;
		resid = Norm(r) / Norm(F);
		return resid;
	}

	//LL^t факторизация
	void SLAE::LLt_factorization()
	{
		double sum = 0., sum2 = 0.;
		L.resize(A.ggl.size());
		//U.resize(A.ggu.size());
		D.resize(A.di.size());

		for (int i = 0; i < n; i++)
		{
			int i0 = A.ig[i];      //индекс первого элемента строки
			int i1 = A.ig[i + 1];  //индекс первого элемента следующей строки
			sum2 = 0.;
			for (int p = i0; p < i1; p++)
			{
				sum = 0.;
				int j = A.jg[p];
				int iline = A.ig[i], jline = A.ig[j];
				int imax = A.ig[i + 1], jmax = A.ig[j + 1];
				for (; iline < imax && jline < jmax;)
				{
					int ki = A.jg[iline], kj = A.jg[jline];
					if (ki < kj) iline++;
					else
						if (kj < ki) jline++;
						else
						{
							sum += L[iline] * L[jline];
							iline++; jline++;
						}
				}
				L[p] = (A.ggl[p] - sum) / D[j];
				sum2 += L[p] * L[p];
			}
			D[i] = sqrt(A.di[i] - sum2);
		}
	}

	//обратный метод нахождения x из L^t * x = y
	void SLAE::ReverseLLt(vector<double>& res, vector<double>& y)
	{
		for (int i = n - 1; i >= 0; i--)
		{
			res[i] /= D[i];
			int p = A.ig[i], i1 = A.ig[i + 1];
			for (; p < i1; p++)
				res[A.jg[p]] -= L[p] * res[i];
		}
	}

	//прямой метод нахождения y из L * y = f
	void SLAE::DirectlyLLt(vector<double>& res, vector<double>& f)
	{
		double sum;
		for (int i = 0; i < n; i++)//проход по всем строкам
		{
			sum = 0.;
			int p = A.ig[i], i1 = A.ig[i + 1];
			for (; p < i1; p++)
				//суммирование произведения элементов l(ij)*y(i)
				sum += L[p] * res[A.jg[p]];
			res[i] = (f[i] - sum) / D[i];
		}
	}

	void SLAE::LLt_MSG()
	{
		double alpha, betta, a, b;
		vector<double>Av(n), q(n);
		LLt_factorization();//LL^t факторизация
		A.MultiplyAx(u, Av);
		for (int i = 0; i < n; i++)
			r[i] = F[i] - Av[i];
		Equal(z, r);
		DirectlyLLt(z, z);
		ReverseLLt(z, z);
		double residual = Residual();
		int k = 0;
		//printf("Количество итераций: %i Невязка: %.15e\r", k, residual);
		for (k = 1; k <= maxiter && residual >= eps; k++)
		{
			A.MultiplyAx(Av, z);
			//ak
			DirectlyLLt(q, r);
			ReverseLLt(q, q);
			a = Scalar(q, r);
			b = Scalar(Av, z);
			alpha = a / b;
			//xk = x[k-1] + ak * z[k-1]
			CulcVect(u, u, alpha, z);
			//rk = r[k-1] - ak * Av[k-1]
			CulcVect(r, r, -alpha, Av);
			residual = Residual();
			//bk
			DirectlyLLt(q, r);
			ReverseLLt(q, q);
			b = Scalar(q, r);
			betta = b / a;
			//zk = r[k] + bk * z[k-1]
			CulcVect(z, q, betta, z);
			//printf("Количество итераций: %i Невязка: %.15e\r", k, residual);
		}
	}

#pragma endregion

}