#include <immintrin.h>

int main() {
	// SSE 2
	__m128d P ;
	double p = 0;
	P   = _mm_set1_pd(p);
	P = _mm_mul_pd(P,P);
	return 0;
}
