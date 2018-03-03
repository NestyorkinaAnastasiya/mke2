/*integration.cpp*/
#include "integration.h"

namespace integration
{
	GaussIntegration::GaussIntegration()
	{
		gaussPoints1[0] = 0;
		gaussPoints1[1] = (1.0 / 3.0) * sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0));
		gaussPoints1[2] = -gaussPoints1[1];
		gaussPoints1[3] = (1.0 / 3.0) * sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0));
		gaussPoints1[4] = -gaussPoints1[3];

		gaussWeights1[0] = 128.0 / 225.0;
		gaussWeights1[1] = (322.0 + 13.0 * sqrt(70.0)) / 900.0;
		gaussWeights1[2] = gaussWeights1[1];
		gaussWeights1[3] = (322.0 - 13.0 * sqrt(70.0)) / 900.0;
		gaussWeights1[4] = gaussWeights1[3];

		for (size_t i = 0; i < n_ip1D; i++)
		{
			for (size_t j = 0; j < n_ip1D; j++)
			{
				gaussPoints[0][i * n_ip1D + j] = gaussPoints1[j];
				gaussPoints[1][i * n_ip1D + j] = gaussPoints1[i];
				gaussWeights[i * n_ip1D + j] = gaussWeights1[i] * gaussWeights1[j];
			}
		}
	}
}
