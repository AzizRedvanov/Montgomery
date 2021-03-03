#include <iostream>
#include <cmath>
#include <ctime>
#include "NTL/ZZ.h"
#include "ExpandingBase.h"

NTL_CLIENT
using namespace std;

void C(ZZ* m, ZZ* r, ZZ& d);
void smenaRNS(ZZ* x, ZZ* s, ZZ* m, ZZ* r, ZZ** Tau, ZZ** c);
void SD(ZZ& x, ZZ* s, ZZ* m, ZZ& r, ZZ** Tau);
void Montgomery(ZZ* xm, ZZ* xr, ZZ* um, ZZ* ur, ZZ* vm, ZZ* vr, ZZ* m, ZZ* r, ZZ& D, ZZ* p, ZZ& M);
void PowerMod_Montgomery(ZZ* xm, ZZ* xr, ZZ* am, ZZ* ar, ZZ& e, ZZ* m, ZZ* r, ZZ& d, ZZ* p, ZZ& M);
const int k{ 4 }, l{ 4 };
ZZ* alpha = new ZZ[k], * gamma = new ZZ[k + 1];
ZZ U{ conv<ZZ>("4049") },
   V{ conv<ZZ>("276") };

struct SrtuctC {
	ZZ** cm = new ZZ * [k];
	ZZ** cr = new ZZ * [k];
	ZZ* cd = new ZZ [k];
}structC;

int main()
{
	ZZ P{ conv<ZZ>("10463") }, D{ conv<ZZ>("2") }, M, R;
	ZZ* p = new ZZ[k + 1],
		* m = new ZZ[k]{ conv<ZZ>("5"),
						 conv<ZZ>("7"),
						 conv<ZZ>("13"),
						 conv<ZZ>("23") },
		* r = new ZZ[k]{ conv<ZZ>("3"),
						 conv<ZZ>("11"),
						 conv<ZZ>("17"),
						 conv<ZZ>("19") },
		* um = new ZZ[k + 1], * ur = new ZZ[k],
		* vm = new ZZ[k + 1], * vr = new ZZ[k];
	set(M);set(R);
	p[0] = mod(P, D);
	for (int i = 0; i < k; i++) {
		operator*=(M, m[i]);
		operator*=(R, r[i]);
		p[i+1] = mod(P, m[i]);
	}

	cout << "R-1=" << modinverse(R, P) << endl;
	cout << "UVR-1 mod P" << endl<<"\t";
	ZZ U1 = MulMod(U, modinverse(R, P), P);
	ZZ proverka = MulMod(U1, V,P);
	for (int i = 0; i < k; i++)
		cout << proverka % m[i] << "\t";
	cout << endl;

	gamma[0] = modinverse(R, D);
	um[0] = mod(U, D);
	vm[0] = mod(V, D);
	for (int i = 0; i < k; i++) {
		alpha[i] = r[i]-modinverse(P, r[i]);
		gamma[i+1] = modinverse(R, m[i]);
		um[i + 1] = mod(U, m[i]);
		ur[i] = mod(U, r[i]);
		vm[i + 1] = mod(V, m[i]);
		vr[i] = mod(V, r[i]);
	}
	C(m, r, D);
	TAU(m,r);
	ZZ* xm = new ZZ[k + 1], * xr = new ZZ[k];
	time_t start, end;
	start = clock();
	Montgomery(xm, xr, um, ur, vm,vr, m, r, D, p, M);
	end = clock();

	cout << "Xm\n";
	for (int i{ 0 }; i <= k; i++)
		cout << xm[i] << "\t";
	cout << "\n\nXr\n";
	for (int i{ 0 }; i < k; i++)
		cout << xr[i] << "\t";
	cout << endl<<endl<< (end - start) / 1000.<<endl;

	delete[] p, m, r, xm,xr, alpha, gamma, um, ur, vm, vr, structC.cd;
	for (int i{ 0 }; i < n; i++)
		delete[] tau.m[i], tau.r[i],
				 structC.cm[i], structC.cr[i];
	return 0;
}

void C(ZZ* m, ZZ* r, ZZ& d) {
	structC.cm[0] = new ZZ[l]; structC.cr[0] = new ZZ[l];
	ZZ mm{ conv<ZZ>("1") }, mr{ conv<ZZ>("1") };
	set(structC.cm[0][0]); set(structC.cr[0][0]); set(structC.cd[0]);
	for (int i{ 1 }; i < k; i++) {
		structC.cm[i] = new ZZ[l]; structC.cr[i] = new ZZ[l];
		set(structC.cm[0][i]); set(structC.cr[0][i]);
		mm *= m[i-1];
		mr *= r[i-1];
		for (int j{ 0 }; j < l; j++) {
			structC.cm[i][j] = mod(mm, r[j]);
			structC.cr[i][j] = mod(mr, m[j]);
		}
		structC.cd[i] = mod(mm, d);
	}
}
void smenaRNS(ZZ* x, ZZ* s, ZZ* m, ZZ* r, ZZ** Tau, ZZ** c) {
	ZZ* b = new ZZ[k];
	RNS_MRNS(b, s, Tau, m);
	ZZ sum;
	for (int j = 0; j < l; j++) {
		clear(sum);
		for (int i = 0; i < k; i++)
			MulAddTo(sum, b[i], c[i][j]);
		x[j] = mod(sum, r[j]);
	}
	delete[] b;
}
void SD(ZZ& x, ZZ* s, ZZ* m, ZZ& d, ZZ** Tau) {
	ZZ* b = new ZZ[k];
	ZZ* sm = new ZZ[k]{ s[1],s[2],s[3],s[4] };
	RNS_MRNS(b, sm, Tau, m);
	ZZ sum{ conv<ZZ>("0") };
	for (int i = 0; i < k; i++)
		MulAddTo(sum, b[i], structC.cd[i]);
	x = mod(sum, d);
	delete[] b, sm;
}

void Montgomery(ZZ* xm, ZZ* xr, ZZ* um, ZZ* ur, ZZ* vm, ZZ* vr, ZZ* m, ZZ* r, ZZ& D, ZZ* p, ZZ& M) {
	ZZ	*tr = new ZZ[k],
		*tm = new ZZ[k+1],
		*qr = new ZZ[k],
		*qm = new ZZ[k+1];

	MulMod(tm[0], um[0], vm[0], D);
	for (int i = 0; i < k; i++) {
		MulMod(tr[i], ur[i], vr[i], r[i]);
		MulMod(qr[i], alpha[i], tr[i], r[i]);
		MulMod(tm[i + 1], um[i + 1], vm[i + 1], m[i]);
	}
	ZZ* Q = new ZZ[k];
	smenaRNS(Q, qr, r, m, tau.r, structC.cr);
	ExpandingBase(qm, Q, m, D, tau.r);
	xm[0] = bit(operator*(operator+(tm[0], operator*(qm[0], p[0])), gamma[0]),0);
	for (int i = 1; i <= k; i++)
		MulMod(xm[i], operator+(tm[i], operator*(qm[i],p[i])), gamma[i], m[i-1]);
	ZZ sd;
	SD(sd, xm, m, D, tau.m);
	ZZ* X = new ZZ[k];
	ZZ XM;
	if ((xm[0] - sd) == bit(M,0)) {
		sub(XM, xm[0], p[0]);
		xm[0] = mod(XM, D);
		for (int i = 1; i <= k; i++) {
			sub(XM, xm[i], p[i]);
			xm[i] = mod(XM, m[i - 1]);
			X[i - 1] = xm[i];
		}
	}
	smenaRNS(xr, X, m, r, tau.m, structC.cm);
	delete[] tr, tm, qr, qm, Q, X;
}

void PowerMod_Montgomery(ZZ* xm, ZZ* xr, ZZ* am, ZZ* ar, ZZ& e, ZZ* m, ZZ* r, ZZ& d, ZZ* p, ZZ& M) {

	while (NumBits(e)) {
		ZZ* bm = am, * br = ar;
		if (bit(e, 0)) {
			Montgomery(xm, xr, am, ar, bm, br, m, r, d, p, M);
		}
		ZZ* Am = new ZZ[k+1], * Ar = new ZZ[k];
		Montgomery(Am, Ar, am, ar, bm, br, m, r, d, p, M);
		am = Am; ar = Ar;
		//Montgomery(a, a, b, p);
		operator>>=(e, 1);
	}
}