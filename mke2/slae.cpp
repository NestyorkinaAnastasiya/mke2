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
		u_n.resize(n);
		r.resize(n);
		z.resize(n);
		// Генерация портрета матрицы и её инициализация
		A.CreatePortret(n, grid);

		// Считываем табличные значения мю
		FILE *fp = fopen("mu.001", "r");
		int size;
		fscanf(fp, "%d", &size);
		tableMu.resize(size);
		for (int i = 0; i < size; i++)
		{
			fscanf(fp, "%lf %lf", &tableMu[i].second, &tableMu[i].first);
			tableMu[i].second *= mu0;
		}
		fclose(fp);

		// Вычисление значений производных полиномов Лагранжа в точках 
		// (tableDetMu)
		CalculateDerivatives();

	}

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
					x_ = x1 + ksi*hx; y_ = y1 + etta*hy;
					
					// Номер материала соответствует железу
					if (element.numberOfMaterial == 2)
					{
						lambda = 1 / CalculateMu(x_, y_, elementNumber);
					}
					else // В остальных случаях нет зависмости от mu
					{
						lambda = 1.0 / mu0;
					}
					lambda = tests.Lambda(x_, y_);

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

	// Сборка локальных матриц масс
	void SLAE::CalculateM(int elementNumber)
	{
		Element element = grid.elements[elementNumber];
		double  g, ksi, etta, sigma, x_, y_, 
			x1 = grid.nodes[element.nodes[0]].x, x3 = grid.nodes[element.nodes[1]].x,
			y1 = grid.nodes[element.nodes[0]].y, y3 = grid.nodes[element.nodes[2]].y,
			hx = x3 - x1, hy = y3 - y1,
			jacobian = hx * hy / 4.0;

		for (int i = 0; i < 4; i++)
		{	
			for (int j = i; j < 4; j++)
			{	
				g = 0;
				for (int k = 0; k < 25; k++)
				{
					ksi = 0.5 + 0.5 * gaussPoints[0][k]; etta = 0.5 + 0.5 * gaussPoints[1][k];
					x_ = x1 + ksi*hx; y_ = y1 + etta*hy;

					sigma = tests.Sigma(x_, y_);
					g +=  sigma *  gaussWeights[k] * phi[i](ksi, etta) * phi[j](ksi, etta);
				}
				M[i][j] = g * jacobian;
			}
		}
		// матрица симметричная, заполняем нижний треугольник
		for (int i = 1; i < 4; i++)
			for (int j = 0; j < i; j++)
				M[i][j] = M[j][i];
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
			
				for (int k = 0; k < 25; k++)
				{
					ksi = 0.5 + 0.5 * gaussPoints[0][k]; etta = 0.5 + 0.5 * gaussPoints[1][k];
					x_ = x1 + ksi*hx; y_ = y1 + etta*hy;

					locF[i] += tests.Fi(x_, y_) * gaussWeights[k] * phi[i](ksi, etta);
				}
			locF[i] *= jacobian;
		}
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
		double tmpM;

		// вычисление локальных матриц
		CalculateG(elementNumber);
		CalculateM(elementNumber);
		CalculateLocalF(elementNumber);

		for (int i = 0; i < 4; i++)
		{
			ki = element.nodes[i];
			for (int j = 0; j < 4; j++)
			{
				kj = element.nodes[j];
				// добавка в глобальную матрицу А
				AddElementToGlobalMatrix(A, ki, kj, G[i][j] + M[i][j]);
				
			}
			// добавка в глобальную правую часть
			F[ki] += locF[i];
		}
	}

	//чо-то не понятное
	int SLAE::detectIntervalB(double B)
	{
		int i;
		int size = tableMu.size();
		for (i = 0; i < size; i++)
			if (B < tableMu[i].first)
				return i - 1; //если левее,то будет возвращаться -1
		return -2;
	}

	//чо-то не понятное
	double SLAE::CalculateMu(double x, double y, int elementNumber)
	{
		double B = CalculateB(x, y, elementNumber);
		int intervalBNumber = detectIntervalB(B);

		if (intervalBNumber == -1) //левее по В
		{
			return tableMu[0].second;
		}
		else
			if (intervalBNumber == -2)//правее по В
			{
				// Вычисление mu по формуле (2.12) из методы
				return 1.0 / ((tableMu.end() - 1)->first / B * (1.0 / (tableMu.end() - 1)->second - 1) + 1);
			}

		return Spline(intervalBNumber, B);
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
	double SLAE::CalculateB(double x, double y, int elementNumber)
	{
		Element element = grid.elements[elementNumber];
		CalculateG(elementNumber);
				
		double qi[4], dAx, dAy, hx, hy, dksi, detta, ksi, etta;

		hx = grid.nodes[element.nodes[1]].x - grid.nodes[element.nodes[0]].x;
		hy = grid.nodes[element.nodes[2]].y - grid.nodes[element.nodes[0]].y;
		ksi = (x - grid.nodes[element.nodes[0]].x) / hx;
		etta = (y - grid.nodes[element.nodes[0]].y) / hy;

		//заполняем соответствующими весами из предыдущего решения по нелинейности
		for (int i = 0; i < 4; i++)
			qi[i] = u_n[element.nodes[i]];

		dAx = 0;
		dAy = 0;
		for (int i = 0; i < 4; i++)
		{
			dAx += dphiksi[i](ksi, etta) * qi[i];
			dAy += dphietta[i](ksi, etta) * qi[i];
		}

		//return sqrt(B2 / S);
		return sqrt(dAx * dAx + dAy * dAy);
	}

	// Построение кубического интерполяционного сплайна
	double SLAE::Spline(int interval_B_number, double B)
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

	// Генерация СЛАУ
	void SLAE::GenerateSLAE()
	{
		// Высчитывание локальных матриц(векторов) и добавление в глобальные
		for (int i = 0; i < grid.elements.size(); i++)
			CalculateLocals(i);

		vector<double> b;
		b.resize(n);
		CalculateBoundaries1();

		normF = Norm(F);
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
			A.di[numberOfNode] = 1.0;

			for (int j = 0; j < matrixSize; j++)
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
				}
		}
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
		int i,k;
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

		for (k = 1; dis > eps && k <= maxiter; k++)
		{
			// Аk
			pp = Scalar(p, p);
			a = Scalar(p, r) / pp;

			// Xk, Rk
			for (i = 0; i < n; i++)
			{
				u[i] = u[i] + a*z[i];
				r[i] = r[i] - a*p[i];
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
		FILE *fo;
		fopen_s(&fo, "result.txt", "w");
		//printf("Calculating......([t] = %f)\n", t);
		// находим новое решение
		GenerateSLAE();
		LULOS();

		for (int i = 0; i < n; i++)
			fprintf(fo, "%.14lf\n", u[i]);
		//fprintf(fo, "%.14lf\t", t);
		fprintf(fo, "\n\n");

		printf("Vector:");
		for (int i = 0; i < n; i++)
			printf("\t\t%.14lf\n", u[i]);
		printf("\n\n");
	}
}