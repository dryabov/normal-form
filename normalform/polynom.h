#pragma once

#include <utility>
#include <complex>
#include <cmath>

#include <boost/unordered_map.hpp>

#include "normalform/monom.h"
#include "normalform/monomcoeff.h"

#ifdef _OPENMP
#include <omp.h>
#endif


namespace normalform {

	using std::complex;

	template<class Tfloat>
	inline bool isZero(const complex<Tfloat> x)
	{
		return (std::abs(x.real()) < (Tfloat)1e-8) && (std::abs(x.imag()) < (Tfloat)1e-8);
	}


	// friend functions
	template<size_t,class Tfloat> class CPolynom;
	template<size_t N,class Tfloat> CPolynom<N,Tfloat> operator +(const CPolynom<N,Tfloat>&, const CPolynom<N,Tfloat>&);
	template<size_t N,class Tfloat> CPolynom<N,Tfloat> operator +(const CPolynom<N,Tfloat>&, const CMonomCoeff<N,Tfloat>&);
	template<size_t N,class Tfloat> CPolynom<N,Tfloat> operator +(const CMonomCoeff<N,Tfloat>&, const CMonomCoeff<N,Tfloat>&);
	template<size_t N,class Tfloat> CPolynom<N,Tfloat> operator -(const CPolynom<N,Tfloat>&, const CPolynom<N,Tfloat>&);
	template<size_t N,class Tfloat> CPolynom<N,Tfloat> operator -(const CPolynom<N,Tfloat>&, const CMonomCoeff<N,Tfloat>&);
	template<size_t N,class Tfloat> CPolynom<N,Tfloat> operator -(const CMonomCoeff<N,Tfloat>&, const CMonomCoeff<N,Tfloat>&);
	template<size_t N,class Tfloat> CPolynom<N,Tfloat> operator *(const CPolynom<N,Tfloat>&, const CPolynom<N,Tfloat>&);
	template<size_t N,class Tfloat> CPolynom<N,Tfloat> operator *(const CPolynom<N,Tfloat>&, const complex<Tfloat>&);
	template<size_t N,class Tfloat> CPolynom<N,Tfloat> operator *(const complex<Tfloat>&, const CPolynom<N,Tfloat>&);
	template<size_t N,class Tfloat> CPolynom<N,Tfloat> operator ^(const CPolynom<N,Tfloat>&, const CPolynom<N,Tfloat>&);
	template<size_t N,class Tfloat> CPolynom<N,Tfloat> operator -(const CPolynom<N,Tfloat>&);

	template<size_t N,class Tfloat=double>
	class CPolynom
	{
	public:
		typedef boost::unordered_map<CMonom<N>,complex<Tfloat> > CMonomMap;
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
		friend CPolynom<N,Tfloat> operator +<>(const CMonomCoeff<N,Tfloat>& m1, const CMonomCoeff<N,Tfloat>& m2);
		friend CPolynom<N,Tfloat> operator -<>(const CPolynom<N,Tfloat>& p1, const CPolynom<N,Tfloat>& p2);
		friend CPolynom<N,Tfloat> operator -<>(const CPolynom<N,Tfloat>& p, const CMonomCoeff<N,Tfloat>& m);
		friend CPolynom<N,Tfloat> operator -<>(const CMonomCoeff<N,Tfloat>& m1, const CMonomCoeff<N,Tfloat>& m2);
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

	template<size_t N,class Tfloat>
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

	template<size_t N,class Tfloat>
	inline CPolynom<N,Tfloat> operator -(const CPolynom<N,Tfloat>& p)
	{
		CPolynom<N,Tfloat> pm(p);
		for(typename CPolynom<N,Tfloat>::CMonomMap::const_iterator it = pm.list.begin(); it != pm.list.end(); ++it)
			pm.list[it->first] = - pm.list[it->first];
		return pm;
	}

	template<size_t N,class Tfloat>
	inline CPolynom<N,Tfloat> operator +(const CPolynom<N,Tfloat>& p1, const CPolynom<N,Tfloat>& p2)
	{
		CPolynom<N,Tfloat> p(p1);
		for(typename CPolynom<N,Tfloat>::CMonomMap::const_iterator it = p2.list.begin(); it != p2.list.end(); ++it)
			p.list[it->first] += it->second;
		return p;
	}

	template<size_t N,class Tfloat>
	inline CPolynom<N,Tfloat>& CPolynom<N,Tfloat>::operator +=(const CPolynom<N,Tfloat>& p)
	{
		for(typename CPolynom<N,Tfloat>::CMonomMap::const_iterator it = p.list.begin(); it != p.list.end(); ++it)
			list[it->first] += it->second;
		return *this;
	}

	template<size_t N,class Tfloat>
	inline CPolynom<N,Tfloat> operator +(const CPolynom<N,Tfloat>& p, const CMonomCoeff<N,Tfloat>& mc)
	{
		CPolynom<N,Tfloat> p1(p);
		p1.list[mc.monom] += mc.coeff;
		return p1;
	}

	template<size_t N,class Tfloat>
	CPolynom<N,Tfloat>& CPolynom<N,Tfloat>::operator +=(const CMonomCoeff<N,Tfloat>& mc)
	{
		list[mc.monom] += mc.coeff;
		return *this;
	}

	template<size_t N,class Tfloat>
	inline CPolynom<N,Tfloat> operator +(const CMonomCoeff<N,Tfloat>& m1, const CMonomCoeff<N,Tfloat>& m2)
	{
		CPolynom<N,Tfloat> p;
		p.list[m1.monom] += m1.coeff;
		p.list[m2.monom] += m2.coeff;
		return p;
	}

	template<size_t N,class Tfloat>
	inline CPolynom<N,Tfloat> operator -(const CPolynom<N,Tfloat>& p1, const CPolynom<N,Tfloat>& p2)
	{
		CPolynom<N,Tfloat> p(p1);
		for(typename CPolynom<N,Tfloat>::CMonomMap::const_iterator it = p2.list.begin(); it != p2.list.end(); ++it)
			p.list[it->first] -= it->second;
		return p;
	}

	template<size_t N,class Tfloat>
	inline CPolynom<N,Tfloat>& CPolynom<N,Tfloat>::operator -=(const CPolynom<N,Tfloat>& p)
	{
		for(typename CPolynom<N,Tfloat>::CMonomMap::const_iterator it = p.list.begin(); it != p.list.end(); ++it)
			list[it->first] -= it->second;
		return *this;
	}

	template<size_t N,class Tfloat>
	inline CPolynom<N,Tfloat> operator -(const CPolynom<N,Tfloat>& p, const CMonomCoeff<N,Tfloat>& mc)
	{
		CPolynom<N,Tfloat> p1(p);
		p1.list[mc.monom] -= mc.coeff;
		return p1;
	}

	template<size_t N,class Tfloat>
	inline CPolynom<N,Tfloat>& CPolynom<N,Tfloat>::operator -=(const CMonomCoeff<N,Tfloat>& mc)
	{
		list[mc.monom] -= mc.coeff;
		return *this;
	}

	template<size_t N,class Tfloat>
	inline CPolynom<N,Tfloat> operator -(const CMonomCoeff<N,Tfloat>& m1, const CMonomCoeff<N,Tfloat>& m2)
	{
		CPolynom<N,Tfloat> p;
		p.list[m1.monom] -= m1.coeff;
		p.list[m2.monom] -= m2.coeff;
		return p;
	}

	template<size_t N,class Tfloat>
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

	template<size_t N,class Tfloat>
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

	template<size_t N,class Tfloat>
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

	template<size_t N,class Tfloat>
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
/*
	template<size_t N,class Tfloat>
	inline CPolynom<N,Tfloat> diff(const CPolynom<N,Tfloat>& p, const size_t j)
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
*/
	template<size_t N,class Tfloat>
	inline CPolynom<N,Tfloat> operator ^(const CPolynom<N,Tfloat>& F, const CPolynom<N,Tfloat>& G)
	{
		//	CPolynom<N,Tfloat> C;
		//	for(size_t j = 0; j < N; j++)
		//		C += diff(F, j) * diff(G, j + N) - diff(G, j) * diff(F, j + N);
		//	return C;

#ifdef _OPENMP
		if(F.list.size() < G.list.size())
		{
			return -(G^F);
		}
#endif

		CPolynom<N,Tfloat> C;
		#pragma omp parallel shared(C)
		{
#ifdef _OPENMP
			int thread_count = omp_get_num_threads();
			int thread_num = omp_get_thread_num();
			size_t chunk_size = F.list.size() / thread_count;

			typename CPolynom<N,Tfloat>::CMonomMap::const_iterator begin = F.list.begin();
			std::advance(begin, thread_num * chunk_size);

			typename CPolynom<N,Tfloat>::CMonomMap::const_iterator end = begin;
			if(thread_num == thread_count - 1)
				end = F.list.end();
			else
				std::advance(end, chunk_size);
#else
			typename CPolynom<N,Tfloat>::CMonomMap::const_iterator begin = F.list.begin();
			typename CPolynom<N,Tfloat>::CMonomMap::const_iterator end = F.list.end();
#endif
			#pragma omp barrier
			CPolynom<N,Tfloat> Cchunk;
			for(typename CPolynom<N,Tfloat>::CMonomMap::const_iterator itF = begin; itF != end; ++itF)
			{
				CMonomCoeff<N,Tfloat> mcF(itF);
				for(typename CPolynom<N,Tfloat>::CMonomMap::const_iterator itG = G.list.begin(); itG != G.list.end(); ++itG)
				{
					CMonomCoeff<N,Tfloat> mcG(itG);
					CMonomCoeff<N,Tfloat> mcFG = mcF*mcG;
					for(size_t j = 0; j < N; j++)
					{
						int diff = mcF.monom[j] * mcG.monom[j+N] - mcG.monom[j] * mcF.monom[j+N];
						if(diff)
						{
							CMonomCoeff<N,Tfloat> mc(mcFG);
							mc.coeff *= complex<Tfloat>(diff);
							mc.monom[j]--;
							mc.monom[j+N]--;
							Cchunk += mc;
						}
					}
				}
			}
			#pragma omp critical
			C += Cchunk;
		}
		C.Simplify();
		return C;
	}

} // namespace normalform