/*slae.h*/
#pragma once
#include "tests.h"
#define _CRT_SECURE_NO_WARNINGS 
#define _USE_MATH_DEFINES
#include <math.h>
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

		// �������� �� � �������(4 * �� * 10 ^ -7)
		double mu0 = 4 * M_PI * 1e-7;
		double Jz = 10e6;
		bool liner = false;
		vector<pair<double, double>> tableMu;
		vector<double> tableDetMu;

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
		double w = 1;

		// ������ ��������� ������ ��������
		void CalculateG(int elementNumber);
		// ������ ��������� ������ ������
		void CalculateLocalF(int elementNumber);
		// ������� ���������� �������� � ����������
		void AddElementToGlobalMatrix(Matrix &B, int i, int j, double element);
		// ������ ��������� ������(��������) � ���������� � ����������
		void CalculateLocals(int elementNumber);

		int detectIntervalB(double B);

		double CalculateMu(double x, double y, int elementNumber);
		void CalculateDerivatives();
		double CalculateB(double x, double y, int elementNumber, double *Bx, double *By);

		double CalculateAzInPoint(double x, double y);
		double CalculateAz(double x, double y, int elementNumber);
		double CalculateBInPoint(double x, double y, double *Bx, double *By);
		double SplineInterpolation(int interval_B_number, double B);

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
		double StopIteration();
		// LU-������������
		void LU();

	public:
		SLAE();

		void calculate_M(int elementNumber);

		void Solve();

		void LFx(const vector<double>& b, vector<double>& result);

		void LTFx(const vector<double>& b, vector<double>& result);

		void UFx(const vector<double>& b, vector<double>& result);

		void UTFx(const vector<double>& b, vector<double>& result);

		double IterMSG(const vector<double>& Az);

		void MultiplyUx(const vector<double>& a, vector<double>& result);

		void LU_MSG();

		void MSG();

		void CulcVect(vector<double>& res, const vector<double>& v1, double var, vector<double>& v2);

		void Equal(vector<double>& vo, vector<double>& vi);

		double Residual();

		void LLt_factorization();

		void ReverseLLt(vector<double>& res, vector<double>& y);

		void DirectlyLLt(vector<double>& res, vector<double>& f);

		void LLt_MSG();

		~SLAE() {};
	};
}
