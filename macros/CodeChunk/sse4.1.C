#include <immintrin.h>

int main() {
	__m128d P ;
	double p = 0;
	P   = _mm_set1_pd(p);
	P = _mm_floor_pd(P);
	return 0;
}
