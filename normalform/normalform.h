#pragma once

#include <complex>
#include <boost/array.hpp>

#include "normalform/monom.h"
#include "normalform/monomcoeff.h"
#include "normalform/polynom.h"

#ifdef NF_LOGGING
#include <iostream>
#endif

namespace normalform {

	using std::complex;
	using boost::array;

	// binom coeffs
	size_t C(const size_t n, size_t k)
	{
		if (0 == k || n == k)
			return 1;
		if (k > n - k)
			k = n - k;
		if (1 == k)
			return n;

		size_t b = n + 1 - k;

		for (size_t i = 2; i <= k; ++i)
		{
			b *= n + i - k;
			b /= i;
		}
		return b;
	}

	template<size_t N,size_t order,class Tfloat=double>
	class NormalForm {
	public:
		typedef array<CPolynom<N,Tfloat>,order> serie;

		serie H, K, S;

		NormalForm(CPolynom<N,Tfloat> p)
		{
			for(typename CPolynom<N,Tfloat>::CMonomMap::const_iterator it = p.list.begin(); it != p.list.end(); ++it)
			{
				int m_order = -2;
				for(size_t i = 0; i < 2*N; i++)
					m_order += it->first[i];
				if(m_order >= 0 && m_order < order)
					H[m_order] += CMonomCoeff<N,Tfloat>(it);
			}
		}

		void normalize()
		{
			CPolynom<N,Tfloat>& H0 = H[0];
			array<complex<Tfloat>,N> lambda;
			array<serie,order> L;

			//Get linear part
#ifdef NF_LOGGING
			std::cout << "Get frequencies...\n";
#endif
			for(typename CPolynom<N,Tfloat>::CMonomMap::const_iterator it = H0.list.begin(); it != H0.list.end(); ++it)
			{
				CMonomCoeff<N,Tfloat> m(it);
				for(size_t i = 0; i < N; i++)
					if(m.monom[i])
					{
						lambda[i] = m.coeff;
						break;
					}
			}

			//Normalization
#ifdef NF_LOGGING
			std::cout << "Normalization...\n";
#endif
			L[0][0] = H[0];
			K[0] = H[0];
			for(size_t n = 1; n < order; n++)
			{
#ifdef NF_LOGGING
				std::cout << "\n" << n << "-th order (" << (n+2) << "-th in H)";
#endif
				L[n][0] = H[n];
				for(size_t i = 1; i <= n; i++)
				{
					L[n][i] = L[n][i-1];
					for(size_t k = 0; k <= n-i; k++)
					{
						//L[n][i] += (complex<Tfloat>)C(n-i,k) * (L[n-1-k][i-1] ^ S[k]);
						CPolynom<N,Tfloat> LS = L[n-1-k][i-1] ^ S[k];
						LS *= (complex<Tfloat>)C(n-i,k);
						L[n][i] += LS;
					}
				}
				K[n] = L[n][n];
				for(typename CPolynom<N,Tfloat>::CMonomMap::const_iterator it = K[n].list.begin(); it != K[n].list.end(); ++it)
				{
					CMonomCoeff<N,Tfloat> f(it);

					complex<Tfloat> res = 0;
					for(size_t i = 0; i < N; i++)
						res += lambda[i] * (complex<Tfloat>)(f.monom[i] - f.monom[i+N]);

					if(!isZero(res))
					{
						f.coeff /= res;
						S[n-1] += f;
					}
				}
				CPolynom<N,Tfloat> dL = H[0] ^ S[n-1];
				for(size_t i = 1; i <= n; i++)
				{
					L[n][i] += dL;
					L[n][i].Simplify();
				}
				K[n] = L[n][n];
			}
		}

		serie getForwardTransform(const size_t i)
		{
			serie X;
			CMonomCoeff<N,Tfloat> mc;
			mc.coeff = 1;
			mc.monom[i]++;

			array<serie,order> Xnj;
			X[0] += mc;
			Xnj[0][0] = X[0];

			for(size_t n = 1; n < order; n++)
			{
#ifdef NF_LOGGING
				std::cout << ".";
#endif
				for(size_t j = 1; j <= n; j++)
				{
					Xnj[n][j] = Xnj[n][j-1];
					for(size_t k = 0; k <= n-j; k++)
					{
						//Xnj[n][j] += (complex<Tfloat>)C(n-j,k) * (Xnj[j+k-1][j-1] ^ S[n-(j+k)]);
						CPolynom<N,Tfloat> XS = Xnj[j+k-1][j-1] ^ S[n-(j+k)];
						XS *= (complex<Tfloat>)C(n-j,k);
						Xnj[n][j] += XS;
					}
				}
				X[n] = Xnj[n][n];
				X[n].Simplify();
			}
			return X;
		}

		serie getBackwardTransform(const size_t i)
		{
			serie Y;
			CMonomCoeff<N,Tfloat> mc;
			mc.coeff = 1;
			mc.monom[i]++;

			array<serie,order> Ynj;
			Y[0] += mc;
			Ynj[0][0] = Y[0];
			for(size_t n = 1; n < order; n++)
			{
#ifdef NF_LOGGING
				std::cout << ".";
#endif
				for(size_t j = n; j > 0; j--)
				{
					Ynj[n][j-1] = Ynj[n][j];
					for(size_t k = 0; k <= n-j; k++)
					{
						//Ynj[n][j-1] -= (complex<Tfloat>)C(n-j,k) * (Ynj[n-k-1][j-1] ^ S[k]);
						CPolynom<N,Tfloat> YS = Ynj[n-k-1][j-1] ^ S[k];
						YS *= (complex<Tfloat>)C(n-j,k);
						Ynj[n][j-1] -= YS;
					}
				}
				Y[n] = Ynj[n][0];
				Y[n].Simplify();
			}
			return Y;
		}
	};

} // namespace normalform