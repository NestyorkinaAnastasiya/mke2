/*basis.cpp*/
#include "basis.h"

namespace basis
{
	Basis::Basis()
	{
		//��������� �� ������� ���������� ���������� �������� ������� � �����
		array <function <double(double)>, nFunc1D> _phi_;
		//��������� �� ������� ���������� d/dksi ���������� �������� ������� � �����
		array <function <double(double)>, nFunc1D> dphi_ksi;

		_phi_[0] = [](double ksi) { return 1 - ksi; };
		_phi_[1] = [](double ksi) { return ksi; };
		dphi_ksi[0] = [](double ksi) { return -1; };
		dphi_ksi[1] = [](double ksi) { return  1; };

		//��������� �� ������� ���������� �������� ������� � �����
		phi[0] = [_phi_](double ksi, double etta) { return _phi_[0](ksi) * _phi_[0](etta); };
		phi[1] = [_phi_](double ksi, double etta) { return _phi_[1](ksi) * _phi_[0](etta); };
		phi[2] = [_phi_](double ksi, double etta) { return _phi_[0](ksi) * _phi_[1](etta); };
		phi[3] = [_phi_](double ksi, double etta) { return _phi_[1](ksi) * _phi_[1](etta); };
		//��������� �� ������� ���������� d/dksi �������� ������� � �����
		dphiksi[0] = [_phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[0](ksi) * _phi_[0](etta); };
		dphiksi[1] = [_phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[1](ksi) * _phi_[0](etta); };
		dphiksi[2] = [_phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[0](ksi) * _phi_[1](etta); };
		dphiksi[3] = [_phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[1](ksi) * _phi_[1](etta); };
		//��������� �� ������� ���������� d/detta �������� ������� � �����
		dphietta[0] = [_phi_, dphi_ksi](double ksi, double etta) { return _phi_[0](ksi) * dphi_ksi[0](etta); };
		dphietta[1] = [_phi_, dphi_ksi](double ksi, double etta) { return _phi_[1](ksi) * dphi_ksi[0](etta); };
		dphietta[2] = [_phi_, dphi_ksi](double ksi, double etta) { return _phi_[0](ksi) * dphi_ksi[1](etta); };
		dphietta[3] = [_phi_, dphi_ksi](double ksi, double etta) { return _phi_[1](ksi) * dphi_ksi[1](etta); };
	}
}
