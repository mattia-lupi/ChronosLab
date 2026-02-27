int solve_RRQR( int const istart_patt, int const iend_patt, int const * const ja_patt,
                int const * const c2glo, int const ntv, double const * const * const TV,
                int const nr_BB_loc, int const l_mm, int const l_nn, double * const BBT,
                lapack_int * const JPVT, int const ind_tau, double * const tau, int const lwork,
                double * const work, double const condmax, int const ind_g,
                double * const g_scr, double * const Rmat, double * const RRT,
                double * const delta, int const ind_BB, double * const BB_scr,
                double * const coef_P0 );
