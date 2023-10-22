#include "util.h"

using namespace std;

int main(){
    typedef chrono::high_resolution_clock elapse_time;
    auto algo1DurationTotal = 0, algo2DurationTotal = 0, indSigDurationTotal = 0, aggDurationTotal = 0, verifyDurationTotal = 0;
    
    publicKey_TTP pkTTP; secretKey_TTP skTTP;
    BIGNUM* pk_target = BN_new(); BIGNUM* sk_target = BN_new(); BIGNUM* r_target = BN_new(); BIGNUM* sig_target = BN_new();
    BIGNUM* sig_list[NUM_USER]; BIGNUM* sig_total = BN_new(); BIGNUM* SL = BN_new();
    BIGNUM* msg = BN_new(); BIGNUM* hash = BN_new(); BIGNUM* gcd = BN_new();
    BIGNUM* eulerPhi = BN_new(); BIGNUM* target_eulerPhi = BN_new();
    BIGNUM* zero = BN_new(); BN_zero(zero);
    BN_CTX* c = BN_CTX_new();

    BIGNUM* res1 = BN_new(); BIGNUM* res2 = BN_new();

    unsigned char* L = (unsigned char*)malloc(Lsize);
    size_t size = 1;
    for(int i = 1; i <= threshold; i++){
        if(i > 255){
            memset(L + i*2-2, i/256, size);
        }
        memset(L + i*2-1, i, size);
    }

    bool flag = false;
    int iteration = 100;
    int cnt = 0;

    bool result_key = false, result_sigverify = false;

    int algo1_elapsed = 0, algo2_elapsed = 0, indSig_elapsed = 0, agg_elapsed = 0, verify_elapsed = 0;

    initParamTTP(pkTTP, skTTP, NUM_USER); initList(sig_list, NUM_USER);

    for(int i = 0; i < iteration; i++){
        cout << cnt << endl;
        // key gen TTP
        while(BN_num_bits(pkTTP.N) != LAMBDA){ // safe prime & non 모두 수행
            BN_generate_prime_ex(skTTP.p, LAMBDA/2, safe_prime, NULL, NULL, NULL);
            BN_generate_prime_ex(skTTP.q, LAMBDA/2, safe_prime, NULL, NULL, NULL);
            BN_mul(pkTTP.N, skTTP.p, skTTP.q, c);
        }
        euler_phi(eulerPhi, pkTTP.N, skTTP.p, skTTP.q);

        // set t
        BN_set_word(pkTTP.t, threshold);

        BN_gcd(gcd, eulerPhi, pkTTP.t, c);
        if(!BN_is_one(gcd)){
            cout << "gcd(eulerPhi, pkTTP.t) is not 1" << endl;
            cout << "gcd : " << BN_bn2dec(gcd) << endl;
            return 0;
        }
        BN_clear(gcd);

        // choose r
        BN_rand_coprime(skTTP.r, eulerPhi);
        // compute d_r
        BN_mod_mul(skTTP.d_r, pkTTP.t, skTTP.r, eulerPhi, c);
        // compute e_r
        BN_mod_inverse(pkTTP.e_r, skTTP.d_r, eulerPhi, c);


        // compute user key
        for(int i = 0; i < NUM_USER; i++){
            // choose e_i, no-duplicate
            do{
                if(e_rand){ // if true : e is random with large size
                    BN_rand_coprime(pkTTP.e_user[i], eulerPhi);
                } else{ // if false : e is small random value
                    BN_rand_coprime(pkTTP.e_user[i], eulerPhi, bits_small_e);
                }
                for(int j = 0; j < i; j++){
                    if(BN_cmp(pkTTP.e_user[j], pkTTP.e_user[i]) == 0){
                        flag = true;
                        break;
                    } else{
                        flag = false;
                    }
                }
            }while(flag); // if user's pk is equal, re-choose

            // compute d_i
            BN_mod_inverse(skTTP.d_user[i], pkTTP.e_user[i], eulerPhi, c);
            // compute d'_i
            BN_mod_add(skTTP.d_user[i], skTTP.d_user[i], skTTP.r, eulerPhi, c);
        }

        //////// algo1 start
        auto algo1Start = elapse_time::now();
        // get_target_euler_phi(target_eulerPhi, pkTTP, skTTP, eulerPhi);

        BN_mul(res1, pkTTP.e_user[1], pkTTP.e_user[2], c);
        BN_mul(res1, res1, skTTP.d_user[1], c);
        
        BN_mul(res2, pkTTP.e_user[1], pkTTP.e_user[2], c);
        BN_mul(res2, res2, skTTP.d_user[2], c);

        BN_sub(target_eulerPhi, res1, res2);
        BN_add(target_eulerPhi, target_eulerPhi, pkTTP.e_user[1]);
        BN_sub(target_eulerPhi, target_eulerPhi, pkTTP.e_user[2]);

        if(BN_is_negative(target_eulerPhi)){
            BN_set_negative(target_eulerPhi, 0);
        }

        auto algo1End = elapse_time::now();
        //////// algo1 end

        // choose target e over target phiN
        // this is given
        BN_rand_coprime(pk_target, eulerPhi);

        //////// algo2 start
        auto algo2Start = elapse_time::now();
        // if gcd(target_eulerPhi, pk_target) != 1, divide target_eulerPhi / gcd
        // do the same with e_r, t
        do{
            BN_clear(gcd);
            BN_gcd(gcd, target_eulerPhi, pk_target, c);
            BN_div(target_eulerPhi, NULL, target_eulerPhi, gcd, c);
        }while(!BN_is_one(gcd));
        BN_clear(gcd);
        do{
            BN_clear(gcd);
            BN_gcd(gcd, target_eulerPhi, pkTTP.e_r, c);
            BN_div(target_eulerPhi, NULL, target_eulerPhi, gcd, c);
        }while(!BN_is_one(gcd));
        BN_clear(gcd);
        do{
            BN_clear(gcd);
            BN_gcd(gcd, target_eulerPhi, pkTTP.t, c);
            BN_div(target_eulerPhi, NULL, target_eulerPhi, gcd, c);
        }while(!BN_is_one(gcd));
        BN_clear(gcd);

        // compute r_k
        BN_mod_mul(r_target, pkTTP.e_r, pkTTP.t, target_eulerPhi, c);
        BN_mod_inverse(r_target, r_target, target_eulerPhi, c);

        // compute d_k
        BN_mod_inverse(sk_target, pk_target, target_eulerPhi, c);
        // compute d'_k
        BN_mod_add(sk_target, sk_target, r_target, target_eulerPhi, c);

        auto algo2End = elapse_time::now();
        //////// algo2 end


        //////// indSig partial start
        auto indSigStart = elapse_time::now();

        // select msg and hash msg
        BN_rand(msg, msgLen, BN_RAND_TOP_ANY, BN_RAND_BOTTOM_ANY);
        sha3_256(hash, msg);

        // compute sig = m^d'_k mod N
        BN_mod_exp(sig_target, hash, sk_target, pkTTP.N, c);

        auto indSigEnd = elapse_time::now();
        //////// indSig partial end

        // key check
        result_key = keyVerify(hash, sig_target, pkTTP, pk_target, eulerPhi);
        if(result_key){
            cout << "key success" << endl;
        } else{
            cout << "key fail" << endl;
        }
///////////////////////////////////
        // forge key pair
        BN_copy(skTTP.d_user[0], sk_target);
        BN_copy(pkTTP.e_user[0], pk_target);

        // forge a signature
        //////// signature aggregate start
        auto aggStart = elapse_time::now();

        BN_one(sig_total);
        for(int i = 0; i < threshold; i++){
            if(i == 0){
                BN_copy(sig_list[i], sig_target); // forge a signature
            } else{
                BN_mod_exp(sig_list[i], hash, skTTP.d_user[i], pkTTP.N, c); // make sig per user without forge
            }
            BN_mod_mul(sig_total, sig_total, sig_list[i], pkTTP.N, c); // product all user's signature and forged signature
        }

        sha3_256(SL, L, Lsize);
        BN_mod_exp(SL, SL, skTTP.d_user[0], pkTTP.N, c); // use forged sk

        auto aggEnd = elapse_time::now();
        //////// signature aggregate end

        //////// verify start
        auto verifyStart = elapse_time::now();
        result_sigverify = AGSVerify(msg, pkTTP.e_user[0], sig_total, L, SL, pkTTP); // use forged pk
        auto verifyEnd = elapse_time::now();
        //////// verify end

        if(result_sigverify){
            cout << "sig verify success" << endl;
        } else{
            cout << "sig verify fail" << endl;
        }

///////////////////////////////////

        algo1_elapsed = chrono::duration_cast<std::chrono::nanoseconds>(algo1End - algo1Start).count();
        algo2_elapsed = chrono::duration_cast<std::chrono::nanoseconds>(algo2End - algo2Start).count();
        indSig_elapsed = chrono::duration_cast<std::chrono::nanoseconds>(indSigEnd - indSigStart).count();
        agg_elapsed = chrono::duration_cast<std::chrono::microseconds>(aggEnd - aggStart).count();
        verify_elapsed = chrono::duration_cast<std::chrono::microseconds>(verifyEnd - verifyStart).count();

        algo1DurationTotal += algo1_elapsed;
        algo2DurationTotal += algo2_elapsed;
        indSigDurationTotal += indSig_elapsed;
        aggDurationTotal += agg_elapsed;
        verifyDurationTotal += verify_elapsed;
        
        cnt++;
        flag = false;

        clearParamTTP(pkTTP, skTTP, NUM_USER); clearList(sig_list, NUM_USER);
        BN_clear(msg); BN_clear(hash); BN_clear(gcd);
        BN_clear(eulerPhi); BN_clear(target_eulerPhi);
        BN_clear(pk_target); BN_clear(sk_target); BN_clear(r_target); BN_clear(sig_target);
        BN_clear(sig_total); BN_clear(SL);
        BN_clear(res1); BN_clear(res2);
    }

    std::cout << "algo1 duration average time(nanoseconds) : " << algo1DurationTotal / iteration << endl;
    std::cout << "algo2 duration average time(nanoseconds) : " << algo2DurationTotal / iteration << endl;
    std::cout << "indSig duration average time(nanoseconds) : " << indSigDurationTotal / iteration << endl;
    std::cout << "agg duration average time(microseconds) : " << aggDurationTotal / iteration << endl;
    std::cout << "verify duration average time(microseconds) : " << verifyDurationTotal / iteration << endl;


    free(L);
    BN_free(zero);

    freeParamTTP(pkTTP, skTTP, NUM_USER); freeList(sig_list, NUM_USER);
    BN_free(msg); BN_free(hash); BN_free(gcd);
    BN_free(eulerPhi); BN_free(target_eulerPhi);
    BN_free(pk_target); BN_free(sk_target); BN_free(r_target); BN_free(sig_target);
    BN_free(sig_total); BN_free(SL);

    BN_CTX_free(c);

    BN_free(res1); BN_free(res2);

    return 0;
}
