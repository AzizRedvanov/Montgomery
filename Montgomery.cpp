#include <iostream>
#include <cmath>
#include <ctime>
#include "NTL/ZZ.h"

NTL_CLIENT
using namespace std;

ZZ mod(ZZ a, ZZ& b) {
	return (a % b < 0 ? (a % b) + b : a % b);
}
struct THREE {
	ZZ x, y, z;
}Three[3];

void egcd(ZZ a1, ZZ b1) {
	ZZ s{ conv<ZZ>(1) }, t{ conv<ZZ>(0) }, x{ conv<ZZ>(0) }, y{ conv<ZZ>(1) };
	ZZ q = a1 / b1, r, old_x, old_y;
	r = mod(a1, b1);
	while (r > 0) {
		a1 = b1;
		b1 = r;
		old_x = x;
		old_y = y;
		x = s - q * x;
		y = t - q * y;
		s = old_x;
		t = old_y;
		q = a1 / b1;
		r = mod(a1, b1);
	}
	Three[2].x = b1;
	Three[2].y = x;
	Three[2].z = y;
}

ZZ modinverse(ZZ c, ZZ d) {
	egcd(c, d);
	ZZ e, x1, y1;
	e = Three[2].x;
	x1 = Three[2].y;
	y1 = Three[2].z;
	if (e > conv<ZZ>(1))
		cout << "No inverse";
	else
		return mod(x1, d);
}

void Montgomery(ZZ& a, ZZ& b);
ZZ p{ power(conv<ZZ>(2), 384) - power(conv<ZZ>(2), 128) - power(conv<ZZ>(2), 96) + power(conv<ZZ>(2), 32) - conv<ZZ>(1) },
	r{ power(conv<ZZ>(2), 384) };
int BitsR = NumBits(r) - 1;
ZZ r1 = modinverse(r, p),
	p1 = ((operator<<(r1, BitsR)) - 1) / p;

int main()
{
	ZZ a{ conv<ZZ>("27580193559959705877849011840389048093056905856361568521428707301988689241309860865136260764883745107765439761230575") },
		b{ conv<ZZ>("27580193559959705877849011840389048093056905856361568521428707301988689241309860865136260764883745107765439761230575") };
	
	ZZ A{ conv<ZZ>("27580193559959705877849011840389048093056905856361568521428707301988689241309860865136260764883745107765439761230575") },
		B;

	time_t start1, end1;
	start1 = time(NULL);
	PowerMod(B, A, 1000000, p);
	end1 = time(NULL);

	time_t start2, end2;
	start2 = time(NULL);
	for (int i = 1; i < 1000000; i++)
		Montgomery(a, b);
	end2 = time(NULL);

	cout << "Classic\tMontgomery" << endl;
	cout << difftime(end1,start1) << "\t" << difftime(end2, start2) << endl;
	return 0;
}

void Montgomery(ZZ& a, ZZ& b) {
	ZZ a1, T;
	a1 = operator<<(a, BitsR) % p;

	T = a1 * b;
	ZZ sh = T * p1;
	ZZ sh1 = operator>>(sh, BitsR);
	operator<<=(sh1, BitsR);
	sh -= sh1;
	a = operator>>((T + (sh * p)), BitsR);

	if (a >= p)
		a -= p;
}