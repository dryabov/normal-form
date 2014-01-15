#pragma once

#include <utility>
#include <complex>

#include <boost/unordered_map.hpp>

#include "normalform/monom.h"
#include "normalform/monomcoeff.h"

namespace normalform {

	using std::complex;

	template<class Tfloat>
	inline bool isZero(const complex<Tfloat> x)
	{
		return abs(x.real()) < (Tfloat)1e-8 && abs(x.imag()) < (Tfloat)1e-8;
	}


	// friend functions
	template<unsigned int N,class Tfloat> class CPolynom;
	template<unsigned int N,class Tfloat> CPolynom<N,Tfloat> operator +(const CPolynom<N,Tfloat>&, const CPolynom<N,Tfloat>&);
	template<unsigned int N,class Tfloat> CPolynom<N,Tfloat> operator +(const CPolynom<N,Tfloat>&, const CMonomCoeff<N,Tfloat>&);
	template<unsigned int N,class Tfloat> CPolynom<N,Tfloat> operator -(const CPolynom<N,Tfloat>&, const CPolynom<N,Tfloat>&);
	template<unsigned int N,class Tfloat> CPolynom<N,Tfloat> operator -(const CPolynom<N,Tfloat>&, const CMonomCoeff<N,Tfloat>&);
	template<unsigned int N,class Tfloat> CPolynom<N,Tfloat> operator *(const CPolynom<N,Tfloat>&, const CPolynom<N,Tfloat>&);
	template<unsigned int N,class Tfloat> CPolynom<N,Tfloat> operator *(const CPolynom<N,Tfloat>&, const complex<Tfloat>&);
	template<unsigned int N,class Tfloat> CPolynom<N,Tfloat> operator *(const complex<Tfloat>&, const CPolynom<N,Tfloat>&);
	template<unsigned int N,class Tfloat> CPolynom<N,Tfloat> operator ^(const CPolynom<N,Tfloat>&, const CPolynom<N,Tfloat>&);
	template<unsigned int N,class Tfloat> CPolynom<N,Tfloat> operator -(const CPolynom<N,Tfloat>&);

	template<unsigned int N,class Tfloat=double>
	class CPolynom
	{
	public:
		typedef boost::unordered_map<CMonom<N>,complex<Tfloat>> CMonomMap;
		CMonomMap list;

		CPolynom<N,Tfloat>()
		{};
		CPolynom<N,Tfloat>(const CPolynom<N,Tfloat>& p)
		{
			list = p.list;
		};
		CPolynom<N,Tfloat>& operator =(const CPolynom<N,Tfloat>& p)
		{
			list = p.list;
			return *this;
		};
		void Clear()
		{
			list.clear();
		};
		void Simplify();
		friend CPolynom<N,Tfloat> operator +<>(const CPolynom<N,Tfloat>& p1, const CPolynom<N,Tfloat>& p2);
		friend CPolynom<N,Tfloat> operator +<>(const CPolynom<N,Tfloat>& p, const CMonomCoeff<N,Tfloat>& m);
		friend CPolynom<N,Tfloat> operator -<>(const CPolynom<N,Tfloat>& p1, const CPolynom<N,Tfloat>& p2);
		friend CPolynom<N,Tfloat> operator -<>(const CPolynom<N,Tfloat>& p, const CMonomCoeff<N,Tfloat>& m);
		friend CPolynom<N,Tfloat> operator *<>(const CPolynom<N,Tfloat>& p1, const CPolynom<N,Tfloat>& p2);
		friend CPolynom<N,Tfloat> operator *<>(const CPolynom<N,Tfloat>& p, const complex<Tfloat>& r);
		friend CPolynom<N,Tfloat> operator *<>(const complex<Tfloat>& r, const CPolynom<N,Tfloat>& p);
		friend CPolynom<N,Tfloat> operator ^<>(const CPolynom<N,Tfloat>& p1, const CPolynom<N,Tfloat>& p2);
		friend CPolynom<N,Tfloat> operator -<>(const CPolynom<N,Tfloat>& p);
		CPolynom<N,Tfloat>& operator +=(const CPolynom<N,Tfloat>& p);
		CPolynom<N,Tfloat>& operator -=(const CPolynom<N,Tfloat>& p);
		CPolynom<N,Tfloat>& operator +=(const CMonomCoeff<N,Tfloat>& m);
		CPolynom<N,Tfloat>& operator -=(const CMonomCoeff<N,Tfloat>& m);
		CPolynom<N,Tfloat>& operator *=(const complex<Tfloat>& r);
	};

	template<unsigned int N,class Tfloat>
	inline void CPolynom<N,Tfloat>::Simplify()
	{
		if(list.empty())
			return;

		for(typename CMonomMap::const_iterator it = list.begin(); it != list.end();)
		{
			if(isZero(it->second))
				list.erase(it++);
			else
				++it;
		}
	}

	template<unsigned int N,class Tfloat>
	inline CPolynom<N,Tfloat> operator -(const CPolynom<N,Tfloat>& p)
	{
		CPolynom<N,Tfloat> pm(p);
		for(typename CPolynom<N,Tfloat>::CMonomMap::const_iterator it = pm.list.begin(); it != pm.list.end(); ++it)
			pm.list[it->first] = - pm.list[it->first];
		return pm;
	}

	template<unsigned int N,class Tfloat>
	inline CPolynom<N,Tfloat> operator +(const CPolynom<N,Tfloat>& p1, const CPolynom<N,Tfloat>& p2)
	{
		CPolynom<N,Tfloat> p(p1);
		for(typename CPolynom<N,Tfloat>::CMonomMap::const_iterator it = p2.list.begin(); it != p2.list.end(); ++it)
			p.list[it->first] += it->second;
		return p;
	}

	template<unsigned int N,class Tfloat>
	inline CPolynom<N,Tfloat>& CPolynom<N,Tfloat>::operator +=(const CPolynom<N,Tfloat>& p)
	{
		for(typename CPolynom<N,Tfloat>::CMonomMap::const_iterator it = p.list.begin(); it != p.list.end(); ++it)
			list[it->first] += it->second;
		return *this;
	}

	template<unsigned int N,class Tfloat>
	inline CPolynom<N,Tfloat> operator +(const CPolynom<N,Tfloat>& p, const CMonomCoeff<N,Tfloat>& mc)
	{
		CPolynom<N,Tfloat> p1(p);
		p1.list[mc.monom] += mc.coeff;
		return p1;
	}

	template<unsigned int N,class Tfloat>
	CPolynom<N,Tfloat>& CPolynom<N,Tfloat>::operator +=(const CMonomCoeff<N,Tfloat>& mc)
	{
		list[mc.monom] += mc.coeff;
		return *this;
	}

	template<unsigned int N,class Tfloat>
	inline CPolynom<N,Tfloat> operator -(const CPolynom<N,Tfloat>& p1, const CPolynom<N,Tfloat>& p2)
	{
		CPolynom<N,Tfloat> p(p1);
		for(typename CPolynom<N,Tfloat>::CMonomMap::const_iterator it = p2.list.begin(); it != p2.list.end(); ++it)
			p.list[it->first] -= it->second;
		return p;
	}

	template<unsigned int N,class Tfloat>
	inline CPolynom<N,Tfloat>& CPolynom<N,Tfloat>::operator -=(const CPolynom<N,Tfloat>& p)
	{
		for(typename CPolynom<N,Tfloat>::CMonomMap::const_iterator it = p.list.begin(); it != p.list.end(); ++it)
			list[it->first] -= it->second;
		return *this;
	}

	template<unsigned int N,class Tfloat>
	inline CPolynom<N,Tfloat> operator -(const CPolynom<N,Tfloat>& p, const CMonomCoeff<N,Tfloat>& mc)
	{
		CPolynom<N,Tfloat> p1(p);
		p1.list[mc.monom] -= mc.coeff;
		return p1;
	}

	template<unsigned int N,class Tfloat>
	inline CPolynom<N,Tfloat>& CPolynom<N,Tfloat>::operator -=(const CMonomCoeff<N,Tfloat>& mc)
	{
		list[mc.monom] -= mc.coeff;
		return *this;
	}

	template<unsigned int N,class Tfloat>
	inline CPolynom<N,Tfloat> operator *(const CPolynom<N,Tfloat>& p1, const CPolynom<N,Tfloat>& p2)
	{
		CPolynom<N,Tfloat> p;
		if(p1.list.empty() || p2.list.empty())
			return p;

		for(typename CPolynom<N,Tfloat>::CMonomMap::const_iterator it1 = p1.list.begin(); it1 != p1.list.end(); ++it1)
		{
			CMonomCoeff<N,Tfloat> mc1(it1);
			for(typename CPolynom<N,Tfloat>::CMonomMap::const_iterator it2 = p2.list.begin(); it2 != p2.list.end(); ++it2)
			{
				CMonomCoeff<N,Tfloat> mc2(it2);
				p += mc1 * mc2;
			}
		}

		return p;
	}

	template<unsigned int N,class Tfloat>
	inline CPolynom<N,Tfloat> operator *(const CPolynom<N,Tfloat>& p, const complex<Tfloat>& r)
	{
		if(isZero(r))
			return CPolynom<N,Tfloat>();

		CPolynom<N,Tfloat> p1(p);
		if(isZero(r - (complex<Tfloat>)1))
			return p1;

		for(typename CPolynom<N,Tfloat>::CMonomMap::const_iterator it = p.list.begin(); it != p.list.end(); ++it)
			p1.list[it->first] = it->second * r;

		return p1;
	}

	template<unsigned int N,class Tfloat>
	inline CPolynom<N,Tfloat> operator *(const complex<Tfloat>& r, const CPolynom<N,Tfloat>& p)
	{
		if(isZero(r))
			return CPolynom<N,Tfloat>();

		CPolynom<N,Tfloat> p1(p);
		if(isZero(r - (complex<Tfloat>)1))
			return p1;

		for(typename CPolynom<N,Tfloat>::CMonomMap::const_iterator it = p.list.begin(); it != p.list.end(); ++it)
			p1.list[it->first] = it->second * r;

		return p1;
	}

	template<unsigned int N,class Tfloat>
	inline CPolynom<N,Tfloat>& CPolynom<N,Tfloat>::operator *=(const complex<Tfloat>& r)
	{
		if(isZero(r))
		{
			Clear();
			return *this;
		}

		if(isZero(r-(complex<Tfloat>)1))
			return *this;

		for(typename CPolynom<N,Tfloat>::CMonomMap::const_iterator it = list.begin(); it != list.end(); ++it)
			list[it->first] = it->second * r;

		return *this;
	}

	template<unsigned int N,class Tfloat>
	inline CPolynom<N,Tfloat> diff(const CPolynom<N,Tfloat>& p, const unsigned int j)
	{
		CPolynom<N,Tfloat> D;
		for(typename CPolynom<N,Tfloat>::CMonomMap::const_iterator it = p.list.begin(); it != p.list.end(); ++it)
			if(it->first[j])
			{
				CMonomCoeff<N,Tfloat> d(it);
				d.coeff *= (complex<Tfloat>)d.monom[j];
				d.monom[j]--;
				D += d;
			}
			return D;
	}

	template<unsigned int N,class Tfloat>
	inline CPolynom<N,Tfloat> operator ^(const CPolynom<N,Tfloat>& F, const CPolynom<N,Tfloat>& G)
	{
		//	CPolynom<N,Tfloat> C;
		//	for(unsigned int j = 0; j < N; j++)
		//		C += diff(F, j) * diff(G, j + N) - diff(G, j) * diff(F, j + N);
		//	return C;

		CPolynom<N,Tfloat> C1, C2;

#pragma omp parallel sections
		{
#pragma omp section
			{
				for(unsigned int j = 0; j < N; j++)
					C1 += diff(F, j) * diff(G, j + N);
			}
#pragma omp section
			{
				for(unsigned int j = 0; j < N; j++)
					C2 += diff(G, j) * diff(F, j + N);
			}
		}
		C1 -= C2;
		C1.Simplify();
		return C1;
	}

} // namespace normalform