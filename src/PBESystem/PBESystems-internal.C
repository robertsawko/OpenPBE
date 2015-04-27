#include "PBESystems-internal.H"

namespace Foam
{

scalar PBESystems::internal::gamma(const scalar& x)
{
	scalar e = constant::mathematical::e;
	scalar pi = constant::mathematical::pi;

	scalar g = 7.0;
	scalar p[] = 
	{
		0.99999999999980993,
		676.5203681218851,
		-1259.1392167224028,
		771.32342877765313,
		-176.61502916214059,
		12.507343278686905,
		-0.13857109526572012,
		9.9843695780195716e-6,
		1.5056327351493116e-7
	};

	scalar value, t, z;

	value = p[0];
	z = x - 1;
	for (label i = 1; i < g + 2; i++)
        {
		value += p[i] / (z + i);
	}
	t = z + g + 0.5;
	return sqrt(2 * pi) * pow(t, z + 0.5) * pow(e, -t) * value;
}

}
