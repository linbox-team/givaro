#include <immintrin.h>
int main() {
	__m256d P ;
	double p = 0;
	P   = _mm256_set1_pd(p);
	P = _mm256_fnmadd_pd(P,P,P);
	return 0;
}
