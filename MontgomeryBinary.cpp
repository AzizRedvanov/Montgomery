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
void LeftShiftMod(ZZ& x, ZZ& a, int k, ZZ& p);
void Montgomery(ZZ& x, ZZ a, ZZ& b, ZZ& p);
void PowerMod_Montgomery(ZZ& x, ZZ& a, ZZ& e, ZZ& p);
ZZ r{ power(conv<ZZ>(2), 384) };
ZZ r1, p1;
int BitsR = NumBits(r) - 1;

int main()
{
	ZZ a{ conv<ZZ>("27580193559959705877849011840389048093056905856361568521428707301988689241309860865136260764883745107765439761230575") },
		b;

	ZZ A{ conv<ZZ>("27580193559959705877849011840389048093056905856361568521428707301988689241309860865136260764883745107765439761230575") },
		B;
	ZZ p{ power(conv<ZZ>(2), 384) - power(conv<ZZ>(2), 128) - power(conv<ZZ>(2), 96) + power(conv<ZZ>(2), 32) - conv<ZZ>(1) };
	r1 = modinverse(r, p);
	p1 = ((operator<<(r1, BitsR)) - 1) / p;

	ZZ e1{ conv<ZZ>("8904809305690585636156852142870730198868924130986086513626076488374165465474984106541006514105454118461") };
	ZZ e2{ conv<ZZ>("8904809305690585636156852142870730198868924130986086513626076488374165465474984106541006514105454118461") };

	time_t start1, end1;
	start1 = clock();
	B = PowerMod(A, e1, p);
	end1 = clock();

	LeftShiftMod(b, a, BitsR, p);
	ZZ x{ conv<ZZ>("1") };
	time_t start2, end2;
	start2 = clock();
	PowerMod_Montgomery(x, b, e2, p);
	end2 = clock();

	cout << "Classic\tMontgomery" << endl;
	cout << (end1 - start1) / 1000. << "\t" << (end2 - start2) / 1000. << endl;
	cout << endl << B << "\n" << x << endl;
	return 0;
}

void LeftShiftMod(ZZ& x, ZZ& a, int k, ZZ& p) {
	int BitsP = NumBits(p),
		K = BitsP - NumBits(a) + 1;
	x = a;
	while (k>K) {
		operator<<=(x, K);
		operator-=(x, p);
		if (operator>=(x, p))
			operator-=(x, p);
		k -= K;
		K = BitsP - NumBits(x) + 1;
	}
	operator<<=(x, k);
	while (operator>=(x, p))
		operator-=(x, p);
}
void Montgomery(ZZ& x, ZZ a, ZZ& b, ZZ& p) {
	ZZ T,sh,sh1;

	mul(T, a, b);
	mul(sh, T, p1);
	RightShift(sh1, sh, BitsR);
	operator<<=(sh1, BitsR);
	operator-=(sh, sh1);
	MulAddTo(T, sh, p);
	RightShift(x, T, BitsR);

	if (operator>=(x, p))
		operator-=(x, p);
}

void PowerMod_Montgomery(ZZ &x, ZZ& a, ZZ& e, ZZ& p) {

	while (NumBits(e)) {
		if (bit(e, 0))
			Montgomery(x, x, a, p);
		ZZ b{ a };
		Montgomery(a, a, b, p);
		operator>>=(e, 1);
	}
}