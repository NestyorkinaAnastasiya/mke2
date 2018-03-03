/*integration.h*/
#pragma once
#include "matrix.h"
namespace integration
{
	const int n_ip = 25;
	const int n_ip1D = 5;
	struct GaussIntegration
	{
		//Точки гаусса
		array<array<double, n_ip>, 2> gaussPoints;
		//Веса гаусса
		array<double, n_ip> gaussWeights;
		//Точки гаусса
		array<double, n_ip1D> gaussPoints1;
		//Веса гаусса
		array<double, n_ip1D> gaussWeights1;

		GaussIntegration();
	};
}