#include <immintrin.h>

int main() {
	__m128i P ;
	int16_t p = 0;
	P   = _mm_set1_epi16(p);
	P = _mm_mulhrs_epi16(P,P);
	return 0;
}
