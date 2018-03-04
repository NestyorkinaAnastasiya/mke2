/*slae.h*/
#pragma once
#include "tests.h"
using namespace matrix;
using namespace basis;
using namespace integration;
using namespace tests;

namespace slae
{
	class SLAE : private Basis, private GaussIntegration
	{
		// ����������� ������
		int n;
		// ������������ ���������� �������� � ��������
		int maxiter = 10000;
		// �������� ������� ����
		const double eps = 1e-10;
		// �����
		Grid grid;
		// ��������� �������� �������
		Tests tests;
		// ���������� �������
		Matrix A;
		// ��������� �������
		// ������� ��������
		array<array<double, 4>, 4> G;
		// ������� �����
		array<array<double, 4>, 4> M;
		array<array<double, 4>, 4> Mg;
		// ��������� ������ ������ �����
		array <double, 4> locF;
		// ���������� ������ ������ �����
		vector <double> F;
		array <double, 4> newF;
		// ������ ������������� ������� �� ���������� ��������
		// �� ������������
		vector <double> u_n;
		// ������ ������������� �������
		vector <double> u;
		// ����� ������� ������ �����
		double normF;

		// ������ ��������� ������ ��������
		void CalculateG(int elementNumber);
		// ������ ��������� ������ ����
		void CalculateM(int elementNumber);
		// ������ ��������� ������ ������
		void CalculateLocalF(int elementNumber);
		// ������� ���������� �������� � ����������
		void AddElementToGlobalMatrix(Matrix &B, int i, int j, double element);
		// ������ ��������� ������(��������) � ���������� � ����������
		void CalculateLocals(int elementNumber);


		// ������ ������� ����� ��� ������� �������� 
		array<double, 3> g;
		// ���� ������� �������� �������
		void CalculateBoundaries1();

		// ���������� ������� � �������������
		vector <double> L;
		vector <double> D;
		vector <double> U;

		// ������ �������
		vector <double> r;
		// ������ ������
		vector <double> z;

		// ���������� ����� �������
		double Norm(const vector<double>& x);
		// ��������� ������������ ��������
		double Scalar(const vector<double>& x, const vector<double>& y);

		// ��������� ���� �� i-�� �������� �� �������
		void GenerateSLAE();
		// LU-������������
		void LU();
		// ��������������� ������� ��� ��������
		void LYF(const vector<double>& C, vector<double>& yl);
		void UXY(const vector<double>& C, vector<double>& yu);
		double Rel_Discrepancy();
		// �������� ��� � LU-�������������
		void LULOS();

	public:
		SLAE();

		void Solve();

		~SLAE() {};
	};
}
