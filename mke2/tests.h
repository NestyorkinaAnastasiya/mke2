/*tests.h*/
#include "integration.h"

namespace tests
{
	struct Tests
	{
		int test;
		double Lambda(double x, double y);
		double Sigma(double x, double y);
		double Gamma(double x, double y);
		double Ug(double x, double y);
		double Betta(int formNumber, double x, double y);
		double Ubetta(int formNumber, double x, double y);
		double Tetta(int formNumber, double x, double y);
		double Fi(double x, double y);
		Tests();
	};
}
