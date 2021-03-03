#pragma once

#include<iostream>
#include <stdio.h>
#include "NTL/ZZ.h"

NTL_CLIENT

const int n = 4;

struct MySrtuct {
    ZZ** m = new ZZ * [n];
    ZZ** r = new ZZ * [n];
}tau;

ZZ mod(ZZ a, ZZ& b) {
    ZZ c = a % b < 0 ? (a % b) + b : a % b;
    return c;
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

void TAU(ZZ* m, ZZ* r) {
    for (int i{ 0 }; i < n; i++) {
        tau.m[i] = new ZZ[n];
        tau.r[i] = new ZZ[n];
        for (int j = i + 1; j < n; j++) {
            tau.m[i][j] = modinverse(m[i], m[j]);
            tau.r[i][j] = modinverse(r[i], r[j]);
        }
    }
}

void RNS_MRNS(ZZ* x, ZZ* a, ZZ** Tau, ZZ* p) {
    x[0] = a[0];
    for (int i{ 1 }; i < n; i++) {
        x[i] = a[i];
        for (int j{ 0 }; j < i; j++) {
            x[i] -= x[j];
            x[i] *= Tau[j][i];
            x[i] = mod(x[i], p[i]);
        }
    }
}

void W(ZZ* p, ZZ q, ZZ* w) {
    ZZ pp{ conv<ZZ>("1") };
    set(w[0]);
    for (int j{ 1 }; j < n; j++) {
        pp *= p[j - 1];
        w[j] = mod(pp, q);
    }
}

void Garner_method(ZZ* p, ZZ q, ZZ* w, ZZ* b, ZZ* x) {
    ZZ summ{ conv<ZZ>("0") };
    for (int j{ 0 }; j < n; j++)
        summ += mod((b[j] * w[j]), q);
    x[0] = mod(summ, q);
}

void ExpandingBase(ZZ* x, ZZ* a, ZZ* p, ZZ q, ZZ** Tau) {

    ZZ* b = new ZZ[n];
    RNS_MRNS(b, a, Tau, p);

    ZZ* w = new ZZ[n];
    W(p, q, w);

    Garner_method(p, q, w, b, x);

    for (int i{ 1 }; i <= n; i++)
        x[i] = a[i-1];
    delete[] b, w;
}