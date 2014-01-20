#pragma once

#include "normalform/monom.h"
#include "normalform/monomcoeff.h"
#include "normalform/polynom.h"

#include <boost/serialization/array.hpp>
#include <boost/serialization/complex.hpp>
#include "normalform/unordered_map_serialization.h"

namespace boost {
	namespace serialization {

		using namespace normalform;

		template<class Archive,size_t N>
		void serialize(Archive &ar, CMonom<N> &m, const unsigned int version)
		{
			ar & m.powers;
		}

		template<class Archive,size_t N,class Tfloat>
		void serialize(Archive &ar, CMonomCoeff<N,Tfloat> &mc, const unsigned int version)
		{
			ar & mc.coeff & mc.monom;
		}

		template<class Archive,size_t N,class Tfloat>
		void serialize(Archive &ar, CPolynom<N,Tfloat> &p, const unsigned int version)
		{
			ar & p.list;
		}

	}
}