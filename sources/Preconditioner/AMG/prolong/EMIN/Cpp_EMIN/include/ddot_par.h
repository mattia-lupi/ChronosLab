#if defined(_MSC_VER)
    #define RESTRICT __restrict
#elif defined(__GNUC__) || defined(__clang__)
    #define RESTRICT __restrict__
#else
    #define RESTRICT
#endif

double ddot_par(const int np, const int nn, const double* RESTRICT x,
                const double* RESTRICT y, double* RESTRICT reduc);
