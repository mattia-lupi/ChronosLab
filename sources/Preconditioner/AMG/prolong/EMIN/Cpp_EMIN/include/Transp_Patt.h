#if defined(_MSC_VER)
    #define RESTRICT __restrict
#elif defined(__GNUC__) || defined(__clang__)
    #define RESTRICT __restrict__
#else
    #define RESTRICT
#endif

int Transp_Patt(const int nthreads, const int nrows, const int ncols, const int nterm,
                const int* RESTRICT iat, const int* RESTRICT ja,
                int* RESTRICT iat_T, int* RESTRICT ja_T, int* RESTRICT perm,
                int* RESTRICT iperm);
