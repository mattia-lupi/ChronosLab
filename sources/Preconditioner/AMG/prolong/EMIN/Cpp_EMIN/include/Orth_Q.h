#if defined(_MSC_VER)
    #define RESTRICT __restrict
#elif defined(__GNUC__) || defined(__clang__)
    #define RESTRICT __restrict__
#else
    #define RESTRICT
#endif

void Orth_Q(const int np, const int nn, const int ntv, const int* RESTRICT iat_patt,
            const double* RESTRICT mat_Q, const double* RESTRICT v_in,
            double* RESTRICT v_ntv, double* RESTRICT v_out);
