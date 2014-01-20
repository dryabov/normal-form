#pragma once

#include <fstream>
#include <map>
#include "normalform/normalform.h"

namespace normalform {

	template<size_t N,size_t order,class Tfloat>
	std::ostream& operator <<(std::ostream& stream, const boost::array<CPolynom<N,Tfloat>,order>& H)
	{
		ios::fmtflags savedFlags(stream.flags());
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
				stream << "+eps";
				if(i > 1)
					stream << "^" << i;
				stream << " ";
			}

			stream << "(";

			// @TODO: convert boost::unordered_map to boost::map for monom ordering
			typedef std::map<CMonom<N>,complex<Tfloat>> ordered_serie;
			ordered_serie serie(H[i].list.begin(), H[i].list.end());

			for(ordered_serie::const_iterator it = serie.begin(); it != serie.end(); ++it)
			{
				if(it->second.imag() == 0)
					stream << std::showpos << it->second.real()/factorial << std::noshowpos;
				else if(it->second.real() == 0)
					stream << std::showpos << it->second.imag()/factorial << std::noshowpos << "*I";
				else
					stream << "+(" << it->second.real()/factorial << std::showpos << it->second.imag()/factorial << std::noshowpos << "*I";

				stream << " ";
				for(size_t i = 0; i < 2*N; i++)
				{
					if(it->first[i])
					{
						if(i < N)
							stream << "q" << (i+1);
						else
							stream << "p" << (i+1-N);

						if(it->first[i] > 1)
							stream << "^" << (size_t)it->first[i];

						stream << " ";
					}
				}
			}

			stream << ")";
		}

		stream.flags(savedFlags);

		return stream;
	}

} // namespace normalform