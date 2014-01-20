#pragma once

#include <fstream>
#include <map>
#include "normalform/normalform.h"

namespace normalform {

	template<size_t N>
	std::ostream& operator <<(std::ostream& stream, const CMonom<N>& m)
	{
		for(size_t i = 0; i < 2*N; i++)
		{
			if(m[i])
			{
				if(i < N)
					stream << "q" << (i+1);
				else
					stream << "p" << (i+1-N);

				if(m[i] > 1)
					stream << "^" << (size_t)m[i];

				stream << " ";
			}
		}
		return stream;
	}

	template<size_t N,class Tfloat>
	std::ostream& operator <<(std::ostream& stream, const CPolynom<N,Tfloat>& p)
	{
		typedef std::map<CMonom<N>,complex<Tfloat> > ordered_serie;
		ordered_serie serie(p.list.begin(), p.list.end());

		for(typename ordered_serie::const_iterator it = serie.begin(); it != serie.end(); ++it)
		{
			if(it->second.imag() == 0)
				stream << std::showpos << it->second.real() << std::noshowpos;
			else if(it->second.real() == 0)
				stream << std::showpos << it->second.imag() << std::noshowpos << "*I";
			else
				stream << "+(" << it->second.real() << std::showpos << it->second.imag() << std::noshowpos << "*I";

			stream << " " << it->first;
		}
		return stream;
	}

	template<size_t N,size_t order,class Tfloat>
	std::ostream& operator <<(std::ostream& stream, const boost::array<CPolynom<N,Tfloat>,order>& H)
	{
		std::ios::fmtflags savedFlags(stream.flags());
		stream << std::noshowpos;

		Tfloat factorial = (Tfloat)1;
		for(size_t i = 0; i < order; i++)
		{
			if(i > 0)
				factorial *= i;

			if(H[i].list.empty())
				continue;

			if(i)
			{
				stream << " + eps";
				if(i > 1)
					stream << "^" << i;

				if(factorial > (Tfloat)1.5)
					stream << "/" << factorial;

				stream << " ";
			}

			stream << "(" << H[i] << ")";
		}

		stream.flags(savedFlags);

		return stream;
	}

} // namespace normalform