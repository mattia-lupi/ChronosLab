#if defined(_MSC_VER)
    #define RESTRICT __restrict
#elif defined(__GNUC__) || defined(__clang__)
    #define RESTRICT __restrict__
#else
    #define RESTRICT
#endif

void spmatv_blk(const int np, const int nblk, const int* RESTRICT pt_blk,
                const int* RESTRICT iat, const int* RESTRICT ja,
                const double* RESTRICT coef, const double* RESTRICT v_in,
                double* RESTRICT v_out);
