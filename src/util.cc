#include "util.h"



void initParamTTP(publicKey_TTP& pk, secretKey_TTP& sk, int n){
    pk.e_r = BN_new();
    pk.N = BN_new();
    pk.t = BN_new();
    for(int i = 0; i < n; i++){
        pk.e_user[i] = BN_new();
    }

    sk.d_r = BN_new();
    sk.p = BN_new();
    sk.q = BN_new();
    sk.r = BN_new();
    for(int i = 0; i < n; i++){
        sk.d_user[i] = BN_new();
    }
}

void clearParamTTP(publicKey_TTP& pk, secretKey_TTP& sk, int n){
    BN_clear(pk.e_r);
    BN_clear(pk.N);
    BN_clear(pk.t);
    for(int i = 0; i < n; i++){
        BN_clear(pk.e_user[i]);
    }

    BN_clear(sk.d_r);
    BN_clear(sk.p);
    BN_clear(sk.q);
    BN_clear(sk.r);
    for(int i = 0; i < n; i++){
        BN_clear(sk.d_user[i]);
    }
}

void freeParamTTP(publicKey_TTP& pk, secretKey_TTP& sk, int n){
    BN_free(pk.e_r);
    BN_free(pk.N);
    BN_free(pk.t);
    for(int i = 0; i < n; i++){
        BN_free(pk.e_user[i]);
    }

    BN_free(sk.d_r);
    BN_free(sk.p);
    BN_free(sk.q);
    BN_free(sk.r);
    for(int i = 0; i < n; i++){
        BN_free(sk.d_user[i]);
    }
}


void initList(BIGNUM* list[], int n){
    for(int i = 0; i < n; i++){
        list[i] = BN_new();
    }
}
void clearList(BIGNUM* list[], int n){
    for(int i = 0; i < n; i++){
        BN_clear(list[i]);
    }
}
void freeList(BIGNUM* list[], int n){
    for(int i = 0; i < n; i++){
        BN_free(list[i]);
    }
}


void euler_phi(BIGNUM* phi, BIGNUM* N, BIGNUM* p, BIGNUM* q){
    BN_sub(phi, N, p);
    BN_sub(phi, phi, q);
    BN_add(phi, phi, BN_value_one());
}

void BN_rand_coprime(BIGNUM* res, BIGNUM* eulerPhi){
    BIGNUM* zero = BN_new();
    BN_zero(zero);

    BIGNUM* gcd = BN_new();
    BN_CTX* c = BN_CTX_new();

    do{
        do{
            BN_rand(res, BN_num_bits(eulerPhi), BN_RAND_TOP_ANY, BN_RAND_BOTTOM_ODD);
        }while(BN_is_zero(res) || BN_is_one(res) || BN_cmp(eulerPhi, res) != 1);
        BN_gcd(gcd, eulerPhi, res, c);
    }while(!BN_is_one(gcd));

    BN_free(gcd);
    BN_free(zero);
    BN_CTX_free(c);
}

void BN_rand_coprime(BIGNUM* res, BIGNUM* eulerPhi, int bits){
    BIGNUM* zero = BN_new();
    BN_zero(zero);

    BIGNUM* gcd = BN_new();
    BN_CTX* c = BN_CTX_new();

    do{
        do{
            BN_rand(res, bits, BN_RAND_TOP_ANY, BN_RAND_BOTTOM_ODD);
        }while(BN_is_zero(res) || BN_is_one(res) || BN_cmp(eulerPhi, res) != 1);
        BN_gcd(gcd, eulerPhi, res, c);
    }while(!BN_is_one(gcd));

    BN_free(gcd);
    BN_free(zero);
    BN_CTX_free(c);
}

void get_target_euler_phi(BIGNUM* res, publicKey_TTP& pkTTP, secretKey_TTP& skTTP, BIGNUM* eulerPhi){
    BIGNUM* res1 = BN_new();
    BIGNUM* res2 = BN_new();

    BN_CTX* c = BN_CTX_new();
    
    BN_mul(res1, pkTTP.e_user[1], pkTTP.e_user[2], c);
    BN_mul(res1, res1, skTTP.d_user[1], c);
    
    BN_mul(res2, pkTTP.e_user[1], pkTTP.e_user[2], c);
    BN_mul(res2, res2, skTTP.d_user[2], c);

    BN_sub(res, res1, res2);
    BN_add(res, res, pkTTP.e_user[1]);
    BN_sub(res, res, pkTTP.e_user[2]);

    if(BN_is_negative(res)){
        BN_set_negative(res, 0);
    }

    BN_free(res1);
    BN_free(res2);
    
    BN_CTX_free(c);
}

bool keyVerify(BIGNUM* msg, BIGNUM* sig, publicKey_TTP& pkTTP, BIGNUM* e, BIGNUM* eulerPhi){
    BIGNUM* res_msg = BN_new();
    BIGNUM* res_sig = BN_new();
    BIGNUM* vexp_msg = BN_new();
    BIGNUM* vexp_sig = BN_new();

    BN_CTX* c = BN_CTX_new();

    // compute e_r*t
    BN_mod_mul(vexp_msg, pkTTP.e_r, pkTTP.t, eulerPhi, c);

    // compute e_r*t*e_i & res_sig
    BN_mod_mul(vexp_sig, vexp_msg, e, eulerPhi, c);
    BN_mod_exp(res_sig, sig, vexp_sig, pkTTP.N, c);

    // compute e_r*t+e_i & res_m
    BN_mod_add(vexp_msg, vexp_msg, e, eulerPhi, c);
    BN_mod_exp(res_msg, msg, vexp_msg, pkTTP.N, c);

    if(BN_cmp(res_msg, res_sig) == 0){
        return true;
    } else{
        return false;
    }

    BN_free(res_msg);
    BN_free(res_sig);
    BN_free(vexp_msg);
    BN_free(vexp_sig); 

    BN_CTX_free(c);
}

void handleErrors()
{
    ERR_print_errors_fp(stderr);
    abort();
}

void sha3_256(BIGNUM* res, const BIGNUM* src)
{
    const EVP_MD *md = NULL;
    EVP_MD_CTX *mdctx;
    unsigned dlen;
    unsigned char* dest;

    md = EVP_sha3_256();
    dlen = SHA256_DIGEST_LENGTH;

    if ((mdctx = EVP_MD_CTX_create()) == NULL)
    {
        handleErrors();
    }

    if (EVP_DigestInit_ex(mdctx, md, NULL) != 1)
    { // returns 1 if successful
        handleErrors();
    }

    unsigned char *bytes = NULL;
    size_t size = BN_num_bytes(src);
    bytes = (unsigned char*)malloc(size);
    BN_bn2bin(src, bytes);
    EVP_DigestUpdate(mdctx, bytes, size);
    free(bytes);

    if ((dest = (unsigned char *)OPENSSL_malloc(dlen)) == NULL)
    {
        handleErrors();
    }

    memset(dest, 0x00, dlen);
    if (EVP_DigestFinal_ex(mdctx, dest, &dlen) != 1)
    { // returns 1 if successful
        OPENSSL_free(dest);
        handleErrors();
    }
    BN_bin2bn(dest, dlen, res);

    EVP_MD_CTX_destroy(mdctx);
}

void sha3_256(BIGNUM* res, const unsigned char* src, int sLen)
{
    const EVP_MD *md = NULL;
    EVP_MD_CTX *mdctx;
    unsigned dlen;
    unsigned char* dest;

    md = EVP_sha3_256();
    dlen = SHA256_DIGEST_LENGTH;

    if ((mdctx = EVP_MD_CTX_create()) == NULL)
    {
        handleErrors();
    }

    if (EVP_DigestInit_ex(mdctx, md, NULL) != 1)
    { // returns 1 if successful
        handleErrors();
    }

    EVP_DigestUpdate(mdctx, src, sLen);

    if ((dest = (unsigned char *)OPENSSL_malloc(dlen)) == NULL)
    {
        handleErrors();
    }

    memset(dest, 0x00, dlen);
    if (EVP_DigestFinal_ex(mdctx, dest, &dlen) != 1)
    { // returns 1 if successful
        OPENSSL_free(dest);
        handleErrors();
    }
    BN_bin2bn(dest, dlen, res);

    EVP_MD_CTX_destroy(mdctx);
}

bool AGSVerify(BIGNUM* msg, BIGNUM* pk_target, BIGNUM* sig_total, unsigned char* L, BIGNUM* SL, publicKey_TTP& pkTTP){
    BIGNUM* h = BN_new(); BIGNUM* h_ = BN_new(); BIGNUM* hl = BN_new();
    BIGNUM* y = BN_new(); BIGNUM* y_ = BN_new();
    BIGNUM* s = BN_new(); BIGNUM* s_ = BN_new();
    BIGNUM* E = BN_new(); BIGNUM* E_ = BN_new(); BIGNUM* EP = BN_new();
    BIGNUM* tmp = BN_new();
    
    BN_CTX* c = BN_CTX_new();

    sha3_256(h, msg);
    sha3_256(hl, L, Lsize);
    
    BN_mul(tmp, pkTTP.e_r, pkTTP.t, c);
    BN_add(tmp, tmp, pk_target);
    BN_mod_exp(y, hl, tmp, pkTTP.N, c);
    BN_clear(tmp);

    BN_mul(tmp, pkTTP.e_r, pkTTP.t, c);
    BN_mul(tmp, tmp, pk_target, c);
    BN_mod_exp(y_, SL, tmp, pkTTP.N, c);
    BN_clear(tmp);

    if(BN_cmp(y, y_) == 0){
        BN_mod_exp(h_, h, pkTTP.e_r, pkTTP.N, c);

        BN_one(E_);
        for(int i = 0; i < threshold; i++){
            BN_mul(E_, E_, pkTTP.e_user[i], c);
        }

        BN_zero(E);
        for(int i = 0; i < threshold; i++){
            BN_div(tmp, NULL, E_, pkTTP.e_user[i], c);
            BN_add(E, E, tmp);
            BN_clear(tmp);
        }

        BN_mod_exp(s, h_, E, pkTTP.N, c);
        
        BN_mod_exp(s_, sig_total, pkTTP.e_r, pkTTP.N, c);
        BN_mod_inverse(tmp, h, pkTTP.N, c);
        BN_mod_mul(s_, s_, tmp, pkTTP.N, c);
        BN_mod_exp(s_, s_, E_, pkTTP.N, c);
        BN_clear(tmp);

        if(BN_cmp(s, s_) == 0){
            return true;
        } else{
            return false;
        }
    } else{
        return false;
    }

    BN_free(h); BN_free(h_); BN_free(hl);
    BN_free(y); BN_free(y_);
    BN_free(s); BN_free(s_);
    BN_free(E); BN_free(E_); BN_free(EP);
    BN_free(tmp);
    
    BN_CTX_free(c);
}

void printhex(unsigned char* hex, int size){
    for(int i = 0; i < size; i++){
        printf("%02x", hex[i]);
    } printf("\n");
}
