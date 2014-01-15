#pragma once

#include <boost/array.hpp>

namespace normalform {

	typedef unsigned char IntPower;

	//friend function
	template<unsigned int N> class CMonom;
	template<unsigned int N> CMonom<N> operator*(CMonom<N>, const CMonom<N>&);

	template<unsigned int N>
	class CMonom
	{
	public:
		boost::array<IntPower,2*N> powers;

		CMonom<N>()
		{
			powers.assign(0);
		};

		CMonom<N>(const CMonom<N>& p)
		{
			powers = p.powers;
		};

		CMonom<N>& operator=(const CMonom<N>& rhs)
		{
			powers = rhs.powers;
			return *this;
		};

		bool operator<(const CMonom<N>& rhs) const
		{
			for(size_t i = 0; i < 2*N; i++)
				if(powers[i] != rhs.powers[i])
					return powers[i] > rhs.powers[i];
			return false;
		};
		bool operator==(const CMonom<N>& rhs) const
		{
			return powers == rhs.powers;
		};
		bool operator!=(const CMonom<N>& rhs) const{return !(*this==rhs);};
		bool operator> (const CMonom<N>& rhs) const{return (rhs<*this);};
		bool operator<=(const CMonom<N>& rhs) const{return !(*this>rhs);};
		bool operator>=(const CMonom<N>& rhs) const{return !(*this<rhs);};

		IntPower& operator[](const unsigned int j)
		{
			return powers[j];
		};
		const IntPower& operator[](const unsigned int j) const
		{
			return powers[j];
		};

		CMonom<N>& operator*=(const CMonom<N>& rhs)
		{
			for(size_t i = 0; i < 2*N; i++)
				powers[i] += rhs.powers[i];
			return *this;
		};

		friend CMonom<N> operator* <>(CMonom<N> lhs, const CMonom<N>& rhs);
	};

	template<unsigned int N>
	inline CMonom<N> operator*(CMonom<N> lhs, const CMonom<N>& rhs)
	{
		lhs *= rhs;
		return lhs;
	}

	template<unsigned int N>
	inline size_t hash_value(const CMonom<N>& m)
	{
		size_t seed = 0;

		for(unsigned int i = 0; i < 2*N; i++)
			boost::hash_combine(seed, m[i]);

		return seed;
	}

} // namespace normalform