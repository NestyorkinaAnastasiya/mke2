/*basis.h*/
#pragma once
#include "grid.h"

using namespace grid;
namespace basis
{
	//число билинейных функций
	const int nFunc = 4;
	//число линейных функций
	const int nFunc1D = 2;

	class Basis
	{
	protected:
		//============Билиненйые базисные функции============//
		//указатели на функции вычисления базисных функций в точке
		array <function <double(double, double)>, nFunc> phi;
		//указатели на функции вычисления d/dksi базисных функций в точке
		array <function <double(double, double)>, nFunc> dphiksi;
		//указатели на функции вычисления d/detta базисных функций в точке
		array <function <double(double, double)>, nFunc> dphietta;

		//конструктор(вычисление базисных функций)
		Basis();
		//деструктор
		~Basis();
	};
}