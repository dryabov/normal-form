#pragma once

#include <complex>
#include "normalform/monom.h"

namespace normalform {

	using std::complex;

	template<size_t N,class Tfloat> class CPolynom;

	//friend function
	template<size_t N,class Tfloat> class CMonomCoeff;
	template<size_t N,class Tfloat> CMonomCoeff<N,Tfloat> operator*(CMonomCoeff<N,Tfloat>, const CMonomCoeff<N,Tfloat>&);

	template<size_t N,class Tfloat=double>
	class CMonomCoeff
	{
	public:
		complex<Tfloat> coeff;
		CMonom<N> monom;

		CMonomCoeff()
		{
			coeff = (Tfloat)0;
		};
		CMonomCoeff(const complex<Tfloat> c, const IntPower powers[2*N])
		{
			coeff = c;
			for(size_t i=0; i<2*N; i++)
				monom[i] = powers[i];
		};
		CMonomCoeff(const typename CPolynom<N,Tfloat>::CMonomMap::const_iterator it)
		{
			coeff = it->second;
			monom = it->first;
		};

		friend CMonomCoeff<N,Tfloat> operator*<>(CMonomCoeff<N,Tfloat> lhs, const CMonomCoeff<N,Tfloat>& rhs);
	};

	template<size_t N,class Tfloat>
	inline CMonomCoeff<N,Tfloat> operator*(CMonomCoeff<N,Tfloat> lhs, const CMonomCoeff<N,Tfloat>& rhs)
	{
		lhs.monom *= rhs.monom;
		lhs.coeff *= rhs.coeff;
		return lhs;
	}

} // namespace normalform