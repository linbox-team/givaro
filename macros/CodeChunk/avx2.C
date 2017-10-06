#include <immintrin.h>
int main() {
	__m256i P ;
	int32_t p = 0;
	P   = _mm256_set1_epi32(p);
	P = _mm256_mul_epi32(P,P);
	return 0;
}
