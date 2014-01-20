// Example for NormalForm library
//
// Normal form for Henon-Heiles hamiltonian

#include "stdafx.h"

#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>

// enable some logging output
#define NF_LOGGING

#include "normalform/normalform.h"
#include "normalform/prettyprint.h"


using namespace normalform;
using namespace std;


#ifdef MSVC
int _tmain(int argc, _TCHAR* argv[])
#else
int main(int argc, char* argv[])
#endif
{
	const int N = 2;     // number of degree of freedom
	const int order = 5; // normalization order

	typedef CPolynom<N>     Polynom;
	typedef CMonomCoeff<N>  Term;
	typedef complex<double> Complex;

	// construct hamiltonian
	IntPower m1000[] = {1,0,0,0}; // Q1 monom powers
	IntPower m0100[] = {0,1,0,0}; // Q2 monom powers
	IntPower m0010[] = {0,0,1,0}; // P1 monom powers
	IntPower m0001[] = {0,0,0,1}; // P2 monom powers
	// complex coordinates to diagonalize linear part
	// q1 = (Q1+i*P1)/sqrt(2),  p1 = (i*Q1+P1)/sqrt(2)
	// q2 = (Q2+i*P2)/sqrt(2),  p2 = (i*Q2+P2)/sqrt(2)
	Polynom q1 = Polynom() + Term(Complex(sqrt(0.5),0), m1000) + Term(Complex(0,sqrt(0.5)), m0010);
	Polynom p1 = Polynom() + Term(Complex(0,sqrt(0.5)), m1000) + Term(Complex(sqrt(0.5),0), m0010);
	Polynom q2 = Polynom() + Term(Complex(sqrt(0.5),0), m0100) + Term(Complex(0,sqrt(0.5)), m0001);
	Polynom p2 = Polynom() + Term(Complex(0,sqrt(0.5)), m0100) + Term(Complex(sqrt(0.5),0), m0001);

	// H = (p1^2+p2^2)/2 + (q1^2+q2^2)/2 + q1^2*q2 - q2^3/3
	Polynom H = Complex(0.5)*(p1*p1 + p2*p2 + q1*q1 + q2*q2) + q1*q1*q2 - Complex(1.0/3.0) * q2*q2*q2;

	H.Simplify();

	// set output precision (15 is maximal for double type)
	cout << std::setprecision(15);

	// initialize normal form
	NormalForm<N,order> NF(H);
	cout << "H=\n" << NF.H << "\n\n";

	// construct normal Form
	NF.normalize();

	// save normal form and generating function
	cout << "K=\n" << NF.K << "\n\n";
	cout << "S=\n" << NF.S << "\n\n";

	// Forward transformation
	for(int i = 0; i < 2 * N; i++)
		cout << "X" << i << "=\n" << NF.getForwardTransform(i) << "\n\n";

	// Backward transformation
	for(int i = 0; i < 2 * N; i++)
		cout << "Y" << i << "=\n" << NF.getBackwardTransform(i) << "\n\n";

	return 0;
}
