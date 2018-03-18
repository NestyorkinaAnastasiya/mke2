#define _CRT_SECURE_NO_WARNINGS 
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
			u[i] = 0;
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
			x1 = grid.nodes[element.nodes[0]].x,
			x3 = grid.nodes[element.nodes[1]].x,
			y1 = grid.nodes[element.nodes[0]].y,
			y3 = grid.nodes[element.nodes[2]].y,
			hx = x3 - x1,
			hy = y3 - y1,
			hx2 = hx * hx,
			hy2 = hy * hy,
			jacobian = hx * hy / 4.0,
			lambda_[4];
		/* lambda_[0] = 1 / CalculateMu(x1, y1, elementNumber);
		 lambda_[1] = 1 / CalculateMu(x3, y1, elementNumber);
		 lambda_[2] = 1 / CalculateMu(x1, y3, elementNumber);
		 lambda_[3] = 1 / CalculateMu(x3, y3, elementNumber);*/
		for (int i = 0; i < 4; i++)
		{
			for (int j = i; j < 4; j++)
			{
				g1 = 0; g2 = 0;
				for (int k = 0; k < 25; k++)
				{
					ksi = 0.5 + 0.5 * gaussPoints[0][k];
					etta = 0.5 + 0.5 * gaussPoints[1][k];
					x_ = x1 + ksi * hx;
					y_ = y1 + etta * hy;

					// Номер материала соответствует железу
					if (element.numberOfMaterial == 2)
					{
						if (!liner)
							lambda = 1 / CalculateMu(x_, y_, elementNumber);
							//lambda = lambda_[j];
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
		/*double Gx11 = 1 / hx, Gx12 = -1 / hx, Gx21 = -1 / hx, Gx22 = 1 / hx,
			Gy11 = 1 / hy, Gy12 = -1 / hy, Gy21 = -1 / hy, Gy22 = 1 / hy,
			Mx11 = hx / 3, Mx12 = hx / 6, Mx21 = hx / 6, Mx22 = hx / 3,
			My11 = hy / 3, My12 = hy / 6, My21 = hy / 6, My22 = hy / 3;
			if (element.numberOfMaterial == 2)
		{
			//if (!liner)
				//lambda = 1 / CalculateMu(x_, y_, elementNumber);
			/*else //lambda = 1 / (1000 * mu0);
		}
		else // В остальных случаях нет зависмости от mu
		{
			lambda = 1.0 / mu0;
		}
		G[0][0] = lambda*(Gx11 * My11 + Mx11 * Gy11);
		G[0][1] = lambda*(Gx12 * My11 + Mx12 * Gy11);
		G[0][2] = lambda*(Gx11 * My12 + Mx11 * Gy12);
		G[0][3] = lambda*(Gx12 * My12 + Mx12 * Gy12);

		G[1][1] = lambda*(Gx22 * My11 + Mx22 * Gy11);
		G[1][2] = lambda*(Gx21 * My12 + Mx21 * Gy12);
		G[1][3] = lambda*(Gx22 * My12 + Mx22 * Gy12);

		G[2][2] = lambda*(Gx11 * My22 + Mx11 * Gy22);
		G[2][3] = lambda*(Gx12 * My22 + Mx12 * Gy22);

		G[3][3] = lambda*(Gx22 * My22 + Mx22 * Gy22);

		// матрица симметричная, заполняем нижний треугольник
		for (int i = 1; i < 4; i++)
			for (int j = 0; j < i; j++)
				G[i][j] = G[j][i];*/
	}

	// Сборка локальных правых частей
	void SLAE::CalculateLocalF(int elementNumber)
	{
		Element element = grid.elements[elementNumber];
		double ksi, etta, x_, y_, s,
			x1 = grid.nodes[element.nodes[0]].x, x3 = grid.nodes[element.nodes[1]].x,
			y1 = grid.nodes[element.nodes[0]].y, y3 = grid.nodes[element.nodes[2]].y,
			hx = x3 - x1, hy = y3 - y1,
			jacobian = hx * hy / 4;

		// интегрирование(Гаусс 5)
		for (int i = 0; i < 4; i++)
		{
			locF[i] = 0;
			if (element.numberOfMaterial == 3 || element.numberOfMaterial == 4)
			{
				for (int k = 0; k < 25; k++)
				{
					ksi = 0.5 + 0.5 * gaussPoints[0][k];
					etta = 0.5 + 0.5 * gaussPoints[1][k];
					x_ = x1 + ksi * hx;
					y_ = y1 + etta * hy;
					if (element.numberOfMaterial == 3)
						locF[i] += 4 * Jz* gaussWeights[k] * phi[i](ksi, etta);
					else if (element.numberOfMaterial == 4)
						locF[i] -= 4 * Jz* gaussWeights[k] * phi[i](ksi, etta);
				}
				locF[i] *= jacobian;
			}
		}/**/
		/*
		double fi;
		Element element = grid.elements[elementNumber];
		double x1 = grid.nodes[element.nodes[0]].x, x3 = grid.nodes[element.nodes[1]].x;
		double y1 = grid.nodes[element.nodes[0]].y, y3 = grid.nodes[element.nodes[2]].y;

		double hx = x3 - x1;
		double hy = y3 - y1;

		double f1, f2, f3, f4;
		if (element.numberOfMaterial == 3)
		f1 = f2 = f3 = f4 = Jz * 4;
		else
		{
		if (element.numberOfMaterial == 4)
		f1 = f2 = f3 = f4 = -Jz*4;
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
			A.di[numberOfNode] = 1e37;
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
		for (int i = 0; i < n; i++)
			A.di[i] = F[i] = 0;
		for (int i = 0; i < A.ggl.size(); i++)
			A.ggl[i] = A.ggu[i] = 0;
		// Высчитывание локальных матриц(векторов) и добавление в глобальные
		for (int i = 0; i < grid.elements.size(); i++)
			CalculateLocals(i);

		CalculateBoundaries1();

		normF = Norm(F);
	}

	double SLAE::StopIteration()
	{
		vector <double> b(n), z(n);

		//генерируем СЛАУ с текущим решением
		GenerateSLAE();
		//умножаем матрицу на текущее решение
		A.MultiplyAx(u, z);
		//и находим абсолютную невязку
		for (int i = 0; i < n; i++)
			b[i] = z[i] - F[i];

		return Norm(b) / Norm(F);
	}

	void SLAE::Solve()
	{
		int k;
		double exit_condition;
		if (liner)
		{
			GenerateSLAE();
			LU_MSG();
		}
		else
		{
			liner = true;
			GenerateSLAE();
			LU_MSG();
			liner = false;
			//цикл по нелинейности
			for (k = 1; k < 50; k++)
			{
				//сохраняем предыдущее решение по нелинейности
				u_n = u;
				//находим новое решение
				GenerateSLAE();
				vector<double>r(n);
				A.Sum(r);

				LU_MSG();
				//проверяем условие останова 
				//+ генерация СЛАУ с текущим решением
				exit_condition = StopIteration();
				printf("%e\r", exit_condition);
				//если оно не выполнено, то, используя параметр
				//релаксации, находим новое решение
				if (exit_condition > eps)
				{
					for (int i = 0; i < n; i++)
						u[i] = w * u[i] + (1 - w)*u_n[i];
				}
				//иначе выходим из цикла
				else break;

			}
		}

		FILE *fo;
		fopen_s(&fo, "result.txt", "w");
		fprintf(fo, "Az:\n");
		for (int i = 0; i < n; i++)
			fprintf(fo, "%.14lf\n", u[i]);
		fprintf(fo, "\n\n");
		FILE *fp;
		fopen_s(&fp, "points.txt", "r");
		int n_;
		fscanf_s(fp, "%d", &n_);
		double x, y, Bx, By, B, Az;
		for (int i = 0; i < n_; i++)
		{
			fscanf_s(fp, "%lf%lf", &x, &y);
			B = CalculateBInPoint(x, y, &Bx, &By);
			Az = CalculateAzInPoint(x, y);
			fprintf(fo, "X= %.3le  Y= %.3le  Bx= %.4le\tBy= %.4le\tAz =  %le\t||B||= %.4le\n", x, y, Bx, By, Az, B);
		}
		fclose(fp);
		fclose(fo);
		fo = fopen("v2.dat", "wb");
		for (int i = 0; i < u.size(); i++)
			fwrite(&u[i], sizeof(double), 1, fo);
		fclose(fo);
	}

#pragma region solver
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
		/*FILE *fn = fopen("v2.dat", "rb");
		for (int i = 0; i < u.size(); i++)
		fread(&u[i], sizeof(double), 1, fn);
		fclose(fn);*/
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
		for (iter = 0; res >= eps && iter < maxiter; iter++)
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
		
		//getchar();

		/*for (int i = 0; i < n; i++)
		fprintf(fo, "%.14lf\n", x[i]);
		fprintf(fo, "%E\t", res);
		fprintf(fo, "%d\t", iter);*/
	}
#pragma endregion

}