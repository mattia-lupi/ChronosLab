#pragma once

int count_merged_row_patt(const int len_1, const int *const ja_1,
                          const int len_2, const int *const ja_2);

void merge_row_patt(const int len_1, const int *const ja_1, const double *const coef_1,
                    const int len_2, const int *const ja_2, const double *const coef_2,
                    int &len_3, int *const ja_3, double *const coef_3);
