#if defined(_MSC_VER)
    #define RESTRICT __restrict
#elif defined(__GNUC__) || defined(__clang__)
    #define RESTRICT __restrict__
#else
    #define RESTRICT
#endif

void mvjcol( const int firstrow, const int nrows, const int* RESTRICT iat,
             const int* RESTRICT ja, int* RESTRICT ja_T, int* RESTRICT punt,
             int* RESTRICT perm);
