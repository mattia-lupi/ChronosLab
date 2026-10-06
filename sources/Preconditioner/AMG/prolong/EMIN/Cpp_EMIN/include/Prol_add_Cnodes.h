#if defined(_MSC_VER)
    #define RESTRICT __restrict
#elif defined(__GNUC__) || defined(__clang__)
    #define RESTRICT __restrict__
#else
    #define RESTRICT
#endif

int Prol_add_Cnodes(const int np, const int nn, const int nn_C,
                    const int* RESTRICT fcnode, const int* RESTRICT iat_in,
                    const int* RESTRICT ja_in, const double* RESTRICT coef_in,
                    int *&iat_out, int *&ja_out, double *&coef_out);
