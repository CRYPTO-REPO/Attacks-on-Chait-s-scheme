#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <chrono>
#include <vector>
#include <string>
#include <algorithm>
#include <string.h>
#include <stdio.h>
#include <cmath>
#include <chrono>

#include <openssl/evp.h>
#include <openssl/sha.h>
#include <openssl/bn.h>
#include <openssl/err.h>

using namespace std;

enum class secuParam{
    lambda_1024 = 1024,
    lambda_2048 = 2048,
    lambda_3072 = 3072,
};

enum class numUser{
    user_2 = 2,
    user_3 = 3,
    user_50 = 50,
    user_100 = 100,
    user_250 = 250,
};

enum class numThreshold{
    th_13 = 13,
    th_29 = 29,
    th_109 = 109,
    th_251 = 251,
    th_503 = 503,
    th_1009 = 1009,
};

#define LAMBDA static_cast<int>(secuParam::lambda_2048)
#define threshold static_cast<int>(numThreshold::th_13)
#define Lsize threshold*2
#define NUM_USER threshold

#define msgLen 256 // bit

#define safe_prime true // true is safe for N
#define e_rand false // true is e random, else small random
#define bits_small_e 40

struct publicKey_TTP{
    BIGNUM* N;
    BIGNUM* t;
    BIGNUM* e_r;
    BIGNUM* e_user[NUM_USER];
};

struct secretKey_TTP{
    BIGNUM* p;
    BIGNUM* q;
    BIGNUM* r;
    BIGNUM* d_r;
    BIGNUM* d_user[NUM_USER];
};

void initParamTTP(publicKey_TTP& pk, secretKey_TTP& sk, int n);
void clearParamTTP(publicKey_TTP& pk, secretKey_TTP& sk, int n);
void freeParamTTP(publicKey_TTP& pk, secretKey_TTP& sk, int n);

void initList(BIGNUM* list[], int n);
void clearList(BIGNUM* list[], int n);
void freeList(BIGNUM* list[], int n);

void euler_phi(BIGNUM* phi, BIGNUM* N, BIGNUM* p, BIGNUM* q);

void BN_rand_coprime(BIGNUM* res, BIGNUM* eulerPhi);
void BN_rand_coprime(BIGNUM* res, BIGNUM* eulerPhi, int bits);

void get_target_euler_phi(BIGNUM* res, publicKey_TTP& pkTTP, secretKey_TTP& skTTP, BIGNUM* eulerPhi);

bool keyVerify(BIGNUM* msg, BIGNUM* sig, publicKey_TTP& pkTTP, BIGNUM* pk_target, BIGNUM* eulerPhi);

void handleErrors();

void sha3_256(BIGNUM* res, const BIGNUM *src);
void sha3_256(BIGNUM* res, const unsigned char* src, int sLen);

bool AGSVerify(BIGNUM* msg, BIGNUM* pk_target, BIGNUM* sig_total, unsigned char* L, BIGNUM* SL, publicKey_TTP& pkTTP);

#endif
