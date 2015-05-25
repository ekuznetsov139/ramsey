#pragma once

#include <emmintrin.h>
#include <immintrin.h>

#define MAX_DEPTH 12
#define MAX_LEN 4
#ifdef WIN32
typedef unsigned __int64 uint64_t;
typedef signed __int64 int64_t;
#else
#include <inttypes.h>
#endif
#ifdef WIN32
#define One 1i64
#else
#define One 1ll
#endif

#ifndef WIN32
#define _BitScanForward(a, b) *(a)=__builtin_ctz(b)
#define _BitScanForward64(a, b) *(a)=__builtin_ctzll(b)
#define _BitScanReverse64(a, b) *(a)=63-__builtin_clzll(b)
#define _popcnt32 __builtin_popcount
#define __popcnt64 __builtin_popcountll
#define _byteswap_uint64 __builtin_bswap64

static uint64_t __rdtsc () 
{
  uint64_t ret;
  uint32_t v1, v2;
  __asm__ __volatile__("rdtsc\n"
 //: "=A" (ret) : :);
 : "=a" (v1), "=d" (v2) : :);
// printf("%d %d\n", v1, v2);
 ret = ((uint64_t)v2)<<32;
 ret |= v1;
  return ret;
}

#endif


inline __m128i bit_reverse(__m128i x)
{
	x = _mm_shuffle_epi8(x, _mm_set_epi8(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15));
	__m128i mm_0f = _mm_set1_epi8(0xf);
	__m128i x0 = _mm_slli_epi16(_mm_and_si128(x, mm_0f), 4);
	__m128i x1 = _mm_and_si128(_mm_srli_epi16(x, 4), mm_0f);
	x = _mm_or_si128(x0, x1);
	__m128i mm_33 = _mm_set1_epi8(0x33);
	x0 = _mm_slli_epi16(_mm_and_si128(x, mm_33), 2);
	x1 = _mm_and_si128(_mm_srli_epi16(x, 2), mm_33);
	x = _mm_or_si128(x0, x1);
	__m128i mm_55 = _mm_set1_epi8(0x55);
	x0 = _mm_slli_epi16(_mm_and_si128(x, mm_55), 1);
	x1 = _mm_and_si128(_mm_srli_epi16(x, 1), mm_55);
	x = _mm_or_si128(x0, x1);
	return x;
}

inline __m128i shift_right(__m128i x, int y)
{
	if(y<64)
	{
		__m128i t = _mm_srl_epi64(x, _mm_set_epi64x(0, y));
		__m128i rem = _mm_sll_epi64(x, _mm_set_epi64x(0, 64-y));
		return  _mm_or_si128(t, _mm_srli_si128(rem, 8));
	}
	else
	{
		__m128i t = _mm_srl_epi64(x, _mm_set_epi64x(0, y-64));
		return _mm_srli_si128(t, 8);
	}
}

inline uint64_t bit_reverse(uint64_t x)
{
	x = _byteswap_uint64(x);
	x = ((x&0x0F0F0F0F0f0f0f0f)<<4) | ((x&0xF0F0F0F0F0F0F0F0)>>4);
	x = ((x&0x3333333333333333)<<2) | ((x&0xcccccccccccccccc)>>2);
	x = ((x&0x5555555555555555)<<1) | ((x&0xaaaaaaaaaaaaaaaa)>>1);
	return x;
}

inline __m128i shift_left(__m128i x, int y)
{
	if(y<64)
	{
		__m128i t = _mm_sll_epi64(x, _mm_set_epi64x(0, y));
		__m128i rem = _mm_srl_epi64(x, _mm_set_epi64x(0, 64-y));
		return  _mm_or_si128(t, _mm_slli_si128(rem, 8));
	}
	else
	{
		__m128i t = _mm_sll_epi64(x, _mm_set_epi64x(0, y-64));
		return _mm_slli_si128(t, 8);
	}
}


struct bare_bigint
{
	uint64_t n[MAX_LEN];

	bare_bigint() 
	{
//		for(int k=0; k<MAX_LEN; k++)
//			n[k]=0;
	}
	void set(int x) { n[x>>6] |= One<<(x&63); }
	void unset(int x) { n[x>>6] &= ~(One<<(x&63)); }

	void unset_from(int x)
	{
		if(x>=MAX_LEN*64)
			return;
		int qw = x >> 6;
		int fract = x & 63;
		n[qw] &= (One<<fract)-One;
		for(int i=qw+1; i<MAX_LEN; i++)
			n[i]=0;
	}

	int bit(int x) const { return (n[x>>6] >> (x&63)) & 1; }
	int word(int x) const
	{
		if((x&63)<=48)
			return (n[x>>6]>>(x&63)) & 0xFFFF;
		else
		{
			uint64_t lo_part = (n[x>>6]>>(x&63));
			uint64_t hi_part = (n[(x>>6)+1]) << (64-(x&63));
			return int(lo_part|hi_part) & 0xFFFF;
		}
	}

	bool zero() const 
	{
#if MAX_LEN==2
		return ((n[0]|n[1])==0);
#else
		for(int k=0; k<MAX_LEN; k++)
			if(n[k]!=0)
				return false;
		return true;
#endif
	}

	int leading_bit() const
	{
		int i=0;
		for(i=0; i<MAX_LEN; i++)
		{
			if(n[i]!=0)
			{
				unsigned long x;
				_BitScanForward64(&x, n[i]);
				return x+i*64;
			}
		}
		return 0;
	}

	int leading_zero() const
	{
		int i=0;
		for(i=0; i<MAX_LEN; i++)
		{
#ifdef WIN32
			if(n[i]!=0xFFFFFFFFFFFFFFFFui64)
#else
			if(n[i]!=0xFFFFFFFFFFFFFFFFull)
#endif
			{
				unsigned long x;
				_BitScanForward64(&x, ~n[i]);
				return x+i*64;
			}
		}
		return 0;
	}

	void operator>>=(int x)
	{
#if MAX_LEN==2
		n[0] = (n[0]>>x) | (n[1]<<(64-x));
		n[1] >>= x;
#else
		for(int i=0; i<MAX_LEN; i++)
		{
			if(i>0)
				n[i-1] |= n[i]<<(64-x);
			n[i] >>= x;
		}
#endif	
	}
	void operator|=(const bare_bigint& x)
	{
#if MAX_LEN==2
		n[0]|=x.n[0];
		n[1]|=x.n[1];
#else
		for(int i=0; i<MAX_LEN; i++)
			n[i]|=x.n[i];
#endif
	}
	void operator&=(const bare_bigint& x)
	{
#if MAX_LEN==2
		n[0]&=x.n[0];
		n[1]&=x.n[1];
#else
		for(int i=0; i<MAX_LEN; i++)
			n[i]&=x.n[i];
#endif
	}

	int trailing_bit() const
	{
#if MAX_LEN==2
		int byte = (n[1]!=0 ? 1 : 0);
		unsigned long x;
		_BitScanReverse64(&x, n[byte]);
		return x+byte*64;
#else
		int i=0;
		for(i=MAX_LEN-1; i>=0; i--)
		{
			if(n[i]!=0)
			{
				unsigned long x;
				_BitScanReverse64(&x, n[i]);
				return x+i*64;
			}
		}
#endif
		return 0;
	}
	void neg()
	{
		int i;
		for(i=0; i<MAX_LEN; i++)
			n[i]^=(uint64_t)-1;
	}

	void left_shift()
	{
		int i;
		for(i=MAX_LEN-1; i>0; i--)
			n[i] = (n[i]<<1) | (n[i-1]>>63);
		n[0]<<=1;
	}
	void right_shift()
	{
		int i;
		for(i=0; i<MAX_LEN-1; i++)
			n[i] = (n[i]>>1) | (n[i+1]<<63);
		n[MAX_LEN-1]>>=1;
	}
	void clear()
	{
		for(int k=0; k<MAX_LEN; k++)
			n[k]=0;
	}

	void printbits() const
	{
		for(int k=0; k<MAX_LEN*64; k++)
			if(bit(k))
				printf("%d ", k);
		printf("\n");
	}


};

struct bigint: bare_bigint
{
	int len;
	int popcnt;
#if !(MAX_LEN & 1)
	uint64_t padding;
#endif
	void set_len(int x) { len=x; }
	bigint()
	{
//		popcnt=0;
	}
	bigint(const bare_bigint& x)
	{
		for(int i=0; i<MAX_LEN; i++)
			n[i]=x.n[i];
	}
	void print() const
	{
		printf("%d ", len);
		if(len>=0)
		{
			for(int k=0; k<len; k++)
				printf("%d", bit(k));
		}
		printf("\n");
	}
	void printbits() const
	{
		printf("(%d/%d) ", popcnt, len);
		for(int k=0; k<=len; k++)
			if(bit(k))
				printf("%d ", k);
		printf("\n");
	}

	bool operator==(const bigint& x) const
	{
		if(len!=x.len)
			return false;
		for(int k=0; k<MAX_LEN; k++)
			if(n[k]!=x.n[k])
				return false;
		return true;
	}
	bool operator<(const bigint& x) const
	{
		if(n[0]!=x.n[0])
			return n[0]<x.n[0];
		return n[1]<x.n[1];
	}
	void clear()
	{
		len=popcnt=0;
		for(int k=0; k<MAX_LEN; k++)
			n[k]=0;
	}
};

inline void store_and(bare_bigint& dst, const bare_bigint& src_a, const bare_bigint& src_b)
{
#if MAX_LEN==2
	__m128i a = _mm_load_si128((const __m128i*)&src_a.n[0]);
	__m128i b = _mm_load_si128((const __m128i*)&src_b.n[0]);
	_mm_store_si128((__m128i*)(&dst.n[0]), _mm_and_si128(a, b));
#else
	for(int i=0; i<MAX_LEN; i++)
		dst.n[i] = src_a.n[i] & src_b.n[i];
#endif
}

inline bool bit_equal(const bigint& a, const bigint& b)
{
#if MAX_LEN==2
	return (a.n[0]==b.n[0] && a.n[1]==b.n[1]);
#else
	int i;
	for(i=0; i<MAX_LEN; i++)
		if(a.n[i]!=b.n[i])
			return false;
	return true;
#endif
}

typedef pair<bigint,bigint> bigint2;


inline uint64_t popcnt(const bare_bigint& a)
{
	uint64_t sum=0;
	for(int k=0; k<MAX_LEN; k++)
		sum += __popcnt64(a.n[k]);
	return sum;
}

inline bare_bigint andnot(const bare_bigint& a, const bare_bigint& b)
{
	bare_bigint x;
	for(int k=0; k<MAX_LEN; k++)
		x.n[k] = a.n[k] & ~b.n[k];
	return x;
}

inline uint64_t popcnt(const bare_bigint& a, const bare_bigint& b)
{
	uint64_t sum=0;
#if MAX_LEN==1
	sum = __popcnt64(a.n[0] & b.n[0]);
#elif MAX_LEN==2
	sum = __popcnt64(a.n[0] & b.n[0])
		+ __popcnt64(a.n[1] & b.n[1]);
#elif MAX_LEN==3
	sum = __popcnt64(a.n[0] & b.n[0])
		+ __popcnt64(a.n[1] & b.n[1])
		+ __popcnt64(a.n[2] & b.n[2]);
#else
	for(int k=0; k<MAX_LEN; k++)
		sum += __popcnt64(a.n[k] & b.n[k]);
#endif
	return sum;
}



inline bool pattern_match(const bigint& a, const bigint& b)
{
#if MAX_LEN==2
	//return ((a.n[0] & b.n[0]) == a.n[0])
	//	&& ((a.n[1] & b.n[1]) == a.n[1]);
	
	__m128i m0 = _mm_loadu_si128((const __m128i*)&a.n[0]);
	__m128i m1 = _mm_loadu_si128((const __m128i*)&b.n[0]);
	__m128i m2 = _mm_cmpeq_epi64(m0, _mm_and_si128(m0,m1));
	return (_mm_movemask_epi8(m2) == 65535);
	
#elif MAX_LEN==4
	__m128i m0 = _mm_andnot_si128(_mm_loadu_si128((const __m128i*)&b.n[0]), _mm_loadu_si128((const __m128i*)&a.n[0]));
	m0 = _mm_or_si128(m0, _mm_andnot_si128(_mm_loadu_si128((const __m128i*)&b.n[2]), _mm_loadu_si128((const __m128i*)&a.n[2])));
	m0 = _mm_cmpeq_epi16(m0, _mm_setzero_si128());
	return (_mm_movemask_epi8(m0) == 65535);
		
#else
	bool ok = true;
	uint64_t mismatch = 0;
	for(int i=0; i<MAX_LEN; i++)
	{
//		if((a.n[i] & b.n[i]) != a.n[i])
//			return false;
		mismatch |= a.n[i] & ~b.n[i];
	}
	return (mismatch==0);
#endif
}

inline bigint operator|(bigint a, bigint b)
{
#if MAX_LEN==2
	a.n[0]|=b.n[0];
	a.n[1]|=b.n[1];
#else
	for(int i=0; i<MAX_LEN; i++)
		a.n[i] |= b.n[i];
#endif
	return a;
}

inline bigint operator&(bigint a, bigint b)
{
#if MAX_LEN==2
	a.n[0]&=b.n[0];
	a.n[1]&=b.n[1];
#else
	for(int i=0; i<MAX_LEN; i++)
		a.n[i] &= b.n[i];
#endif
	return a;
}

inline bigint operator<<(bigint x, int n)
{
	bigint retval;
#if MAX_LEN==2
	retval.set_len(x.len);
	retval.popcnt = x.popcnt;
	__m128i m = _mm_loadu_si128((const __m128i*)&x.n[0]);
	m = shift_left(m, n);
	_mm_storeu_si128((__m128i*)&retval.n[0], m);

#else
	int i;
	retval.clear();
	retval.set_len(x.len);
	retval.popcnt = x.popcnt;
	int qw_n = n >> 6;
	int fract_n = n & 63;
	if(fract_n != 0)
	{
		for(i=0; i+qw_n<MAX_LEN; i++)
		{
			retval.n[i+qw_n] |= x.n[i] << fract_n;
			if(i+qw_n+1<MAX_LEN)
				retval.n[i+qw_n+1] |= x.n[i] >> (64-fract_n);
		}
	}
	else
	{
		for(i=0; i+qw_n<MAX_LEN; i++)
			retval.n[i+qw_n] = x.n[i];
	}
#endif
	return retval;
}

inline bigint operator>>(bigint x, int n)
{
	bigint retval;
	int i;
	retval.clear();
	retval.set_len(x.len);
	retval.popcnt = x.popcnt;
	int qw_n = n >> 6;
	int fract_n = n & 63;
	if(fract_n != 0)
	{
		for(i=0; i+qw_n<MAX_LEN; i++)
		{
			retval.n[i] |= x.n[i+qw_n] >> fract_n;
			if(i+qw_n+1<MAX_LEN)
				retval.n[i] |= x.n[i+qw_n+1] << (64-fract_n);
		}
	}
	else
	{
		for(i=0; i+qw_n<MAX_LEN; i++)
			retval.n[i] = x.n[i+qw_n];
	}
	return retval;
}


inline bool pattern_match(const bigint& a, const bigint& b, int shift)
{ 
#if MAX_LEN==2
	__m128i m0 = _mm_loadu_si128((const __m128i*)&a.n[0]);
	__m128i m1 = _mm_loadu_si128((const __m128i*)&b.n[0]);
	int shift_save = shift;
	if(shift >= 64)
	{
		m0 = _mm_unpacklo_epi64(_mm_setzero_si128(), m0);
		shift &= 63;
	}
	if(shift != 0)
	{
		m0 = _mm_or_si128(_mm_slli_epi64(m0, shift), _mm_unpacklo_epi64(_mm_setzero_si128(), _mm_srli_epi64(m0, 64-shift)));
	}
	__m128i m2 = _mm_cmpeq_epi64(m0, _mm_and_si128(m0,m1));
	bool retval = (_mm_movemask_epi8(m2) == 65535);
	return retval;
#else
	return pattern_match(a<<shift, b);
#endif
}

inline int inv_batch_pattern_match(const bigint& a, const bigint* pb)
{
#if MAX_LEN==2
	__m128i m0 = _mm_loadu_si128((const __m128i*)&a.n[0]);
	
	__m128i r0 = _mm_loadu_si128((const __m128i*)&pb[0].n[0]);
	__m128i r1 = _mm_loadu_si128((const __m128i*)&pb[1].n[0]);
	__m128i r2 = _mm_loadu_si128((const __m128i*)&pb[2].n[0]);
	__m128i r3 = _mm_loadu_si128((const __m128i*)&pb[3].n[0]);
	r0 = _mm_cmpeq_epi64(r0, _mm_and_si128(m0,r0));
	r1 = _mm_cmpeq_epi64(r1, _mm_and_si128(m0,r1));
	r2 = _mm_cmpeq_epi64(r2, _mm_and_si128(m0,r2));
	r3 = _mm_cmpeq_epi64(r3, _mm_and_si128(m0,r3));
	r0 = _mm_packs_epi32(r0, r1);
	r2 = _mm_packs_epi32(r2, r3);

	__m128i r4 = _mm_loadu_si128((const __m128i*)&pb[4].n[0]);
	__m128i r5 = _mm_loadu_si128((const __m128i*)&pb[5].n[0]);
	__m128i r6 = _mm_loadu_si128((const __m128i*)&pb[6].n[0]);
	__m128i r7 = _mm_loadu_si128((const __m128i*)&pb[7].n[0]);
	r4 = _mm_cmpeq_epi64(r4, _mm_and_si128(m0,r4));
	r5 = _mm_cmpeq_epi64(r5, _mm_and_si128(m0,r5));
	r6 = _mm_cmpeq_epi64(r6, _mm_and_si128(m0,r6));
	r7 = _mm_cmpeq_epi64(r7, _mm_and_si128(m0,r7));
	r4 = _mm_packs_epi32(r4, r5);
	r6 = _mm_packs_epi32(r6, r7);

	r0 = _mm_packs_epi32(r0, r2);
	r4 = _mm_packs_epi32(r4, r6);
	r0 = _mm_packs_epi16(r0, r4);

	__m128i s0 = r0;
	pb+=8;
	r0 = _mm_loadu_si128((const __m128i*)&pb[0].n[0]);
	r1 = _mm_loadu_si128((const __m128i*)&pb[1].n[0]);
	r2 = _mm_loadu_si128((const __m128i*)&pb[2].n[0]);
	r3 = _mm_loadu_si128((const __m128i*)&pb[3].n[0]);
	r0 = _mm_cmpeq_epi64(r0, _mm_and_si128(m0,r0));
	r1 = _mm_cmpeq_epi64(r1, _mm_and_si128(m0,r1));
	r2 = _mm_cmpeq_epi64(r2, _mm_and_si128(m0,r2));
	r3 = _mm_cmpeq_epi64(r3, _mm_and_si128(m0,r3));
	r0 = _mm_packs_epi32(r0, r1);
	r2 = _mm_packs_epi32(r2, r3);

	r4 = _mm_loadu_si128((const __m128i*)&pb[4].n[0]);
	r5 = _mm_loadu_si128((const __m128i*)&pb[5].n[0]);
	r6 = _mm_loadu_si128((const __m128i*)&pb[6].n[0]);
	r7 = _mm_loadu_si128((const __m128i*)&pb[7].n[0]);
	r4 = _mm_cmpeq_epi64(r4, _mm_and_si128(m0,r4));
	r5 = _mm_cmpeq_epi64(r5, _mm_and_si128(m0,r5));
	r6 = _mm_cmpeq_epi64(r6, _mm_and_si128(m0,r6));
	r7 = _mm_cmpeq_epi64(r7, _mm_and_si128(m0,r7));
	r4 = _mm_packs_epi32(r4, r5);
	r6 = _mm_packs_epi32(r6, r7);

	r0 = _mm_packs_epi32(r0, r2);
	r4 = _mm_packs_epi32(r4, r6);
	r0 = _mm_packs_epi16(r0, r4);
	r0 = _mm_packs_epi16(s0, r0);
	__m128i rsave = r0;
	__m128i mm_ones=_mm_cmpeq_epi8(r0, r0);
	r0 = _mm_cmpeq_epi8(r0, mm_ones);
	int retval = _mm_movemask_epi8(r0);
//	int retval = _mm_movemask_epi8(r0);
	return retval;
#else
	int retval=0;
	for(int i=0; i<16; i++)
		retval |= (pattern_match(pb[i],a) << i);
	return retval;
#endif
}

inline int batch_pattern_match(const bigint& a, const bigint* pb)
{
#if MAX_LEN==2
	__m128i m0 = _mm_loadu_si128((const __m128i*)&a.n[0]);
	
	__m128i r0 = _mm_loadu_si128((const __m128i*)&pb[0].n[0]);
	__m128i r1 = _mm_loadu_si128((const __m128i*)&pb[1].n[0]);
	__m128i r2 = _mm_loadu_si128((const __m128i*)&pb[2].n[0]);
	__m128i r3 = _mm_loadu_si128((const __m128i*)&pb[3].n[0]);
	r0 = _mm_cmpeq_epi64(m0, _mm_and_si128(m0,r0));
	r1 = _mm_cmpeq_epi64(m0, _mm_and_si128(m0,r1));
	r2 = _mm_cmpeq_epi64(m0, _mm_and_si128(m0,r2));
	r3 = _mm_cmpeq_epi64(m0, _mm_and_si128(m0,r3));
	r0 = _mm_packs_epi32(r0, r1);
	r2 = _mm_packs_epi32(r2, r3);

	__m128i r4 = _mm_loadu_si128((const __m128i*)&pb[4].n[0]);
	__m128i r5 = _mm_loadu_si128((const __m128i*)&pb[5].n[0]);
	__m128i r6 = _mm_loadu_si128((const __m128i*)&pb[6].n[0]);
	__m128i r7 = _mm_loadu_si128((const __m128i*)&pb[7].n[0]);
	r4 = _mm_cmpeq_epi64(m0, _mm_and_si128(m0,r4));
	r5 = _mm_cmpeq_epi64(m0, _mm_and_si128(m0,r5));
	r6 = _mm_cmpeq_epi64(m0, _mm_and_si128(m0,r6));
	r7 = _mm_cmpeq_epi64(m0, _mm_and_si128(m0,r7));
	r4 = _mm_packs_epi32(r4, r5);
	r6 = _mm_packs_epi32(r6, r7);

	r0 = _mm_packs_epi32(r0, r2);
	r4 = _mm_packs_epi32(r4, r6);
	r0 = _mm_packs_epi16(r0, r4);

	__m128i s0 = r0;
	pb+=8;
	r0 = _mm_loadu_si128((const __m128i*)&pb[0].n[0]);
	r1 = _mm_loadu_si128((const __m128i*)&pb[1].n[0]);
	r2 = _mm_loadu_si128((const __m128i*)&pb[2].n[0]);
	r3 = _mm_loadu_si128((const __m128i*)&pb[3].n[0]);
	r0 = _mm_cmpeq_epi64(m0, _mm_and_si128(m0,r0));
	r1 = _mm_cmpeq_epi64(m0, _mm_and_si128(m0,r1));
	r2 = _mm_cmpeq_epi64(m0, _mm_and_si128(m0,r2));
	r3 = _mm_cmpeq_epi64(m0, _mm_and_si128(m0,r3));
	r0 = _mm_packs_epi32(r0, r1);
	r2 = _mm_packs_epi32(r2, r3);

	r4 = _mm_loadu_si128((const __m128i*)&pb[4].n[0]);
	r5 = _mm_loadu_si128((const __m128i*)&pb[5].n[0]);
	r6 = _mm_loadu_si128((const __m128i*)&pb[6].n[0]);
	r7 = _mm_loadu_si128((const __m128i*)&pb[7].n[0]);
	r4 = _mm_cmpeq_epi64(m0, _mm_and_si128(m0,r4));
	r5 = _mm_cmpeq_epi64(m0, _mm_and_si128(m0,r5));
	r6 = _mm_cmpeq_epi64(m0, _mm_and_si128(m0,r6));
	r7 = _mm_cmpeq_epi64(m0, _mm_and_si128(m0,r7));
	r4 = _mm_packs_epi32(r4, r5);
	r6 = _mm_packs_epi32(r6, r7);

	r0 = _mm_packs_epi32(r0, r2);
	r4 = _mm_packs_epi32(r4, r6);
	r0 = _mm_packs_epi16(r0, r4);
	r0 = _mm_packs_epi16(s0, r0);
	__m128i rsave = r0;
	__m128i mm_ones=_mm_cmpeq_epi8(r0, r0);
	r0 = _mm_cmpeq_epi8(r0, mm_ones);
	int retval = _mm_movemask_epi8(r0);
//	int retval = _mm_movemask_epi8(r0);
	return retval;
	/*
	pb-=8;
	int check=0;
	for(int i=0; i<16; i++)
		check |= (pattern_match(a,pb[i]) << i);
	if(check != retval)
	{
		printf("%x %x\n", check, retval);
		print(rsave);
		print(r0);

		exit(-1);
	}
	return retval;
*/
#else
	int retval=0;
	for(int i=0; i<16; i++)
		retval |= (pattern_match(a,pb[i]) << i);
	return retval;
#endif
}

inline int batch_pattern_match_2(const bigint& a, const bigint* pb)
{
#if MAX_LEN==2
	__m128i m0 = _mm_loadu_si128((const __m128i*)&a.n[0]);
	//m0 = _mm_srli_si128(m0, shift);
//	m0 = _mm_alignr_epi8(_mm_setzero_si128(), m0, shift);
	__m128i m0_ = _mm_srli_si128(m0,1);

//	__m128i m04 = _mm_or_si128(_mm_srli_epi64(m0,4), _mm_srli_si128(_mm_slli_epi64(m0, 60), 8));
//	__m128i m04_ = _mm_or_si128(_mm_srli_epi64(m0_,4), _mm_srli_si128(_mm_slli_epi64(m0_, 60), 8));


	__m128i r0 = _mm_loadu_si128((const __m128i*)&pb[0].n[0]);
	__m128i r1 = _mm_loadu_si128((const __m128i*)&pb[1].n[0]);
	__m128i r2 = _mm_loadu_si128((const __m128i*)&pb[2].n[0]);
	__m128i r3 = _mm_loadu_si128((const __m128i*)&pb[3].n[0]);

	__m128i b0 = _mm_cmpeq_epi64(m0, _mm_and_si128(m0,r0));
	__m128i b0_ = _mm_cmpeq_epi64(m0_, _mm_and_si128(m0_, r0));
	__m128i b1 = _mm_cmpeq_epi64(m0, _mm_and_si128(m0,r1));
	__m128i b1_ = _mm_cmpeq_epi64(m0_, _mm_and_si128(m0_, r1));
	__m128i b2 = _mm_cmpeq_epi64(m0, _mm_and_si128(m0,r2));
	__m128i b2_ = _mm_cmpeq_epi64(m0_, _mm_and_si128(m0_, r2));
	__m128i b3 = _mm_cmpeq_epi64(m0, _mm_and_si128(m0,r3));
	__m128i b3_ = _mm_cmpeq_epi64(m0_, _mm_and_si128(m0_, r3));
	b0 = _mm_packs_epi32(b0, b1);
	b2 = _mm_packs_epi32(b2, b3);
	b0_ = _mm_packs_epi32(b0_, b1_);
	b2_ = _mm_packs_epi32(b2_, b3_);
	__m128i dword0 = _mm_packs_epi32(b0, b2);
	__m128i dword2 = _mm_packs_epi32(b0_, b2_);


//	m0 = _mm_or_si128(_mm_srli_epi64(m0,4), _mm_srli_si128(_mm_slli_epi64(m0, 60), 8));
//	m0_ = _mm_or_si128(_mm_srli_epi64(m0_,4), _mm_srli_si128(_mm_slli_epi64(m0_, 60), 8));

	r0 = _mm_loadu_si128((const __m128i*)&pb[4].n[0]);
	r1 = _mm_loadu_si128((const __m128i*)&pb[5].n[0]);
	r2 = _mm_loadu_si128((const __m128i*)&pb[6].n[0]);
	r3 = _mm_loadu_si128((const __m128i*)&pb[7].n[0]);
	
	b0 = _mm_cmpeq_epi64(m0, _mm_and_si128(m0,r0));
	b0_ = _mm_cmpeq_epi64(m0_, _mm_and_si128(m0_, r0));
	b1 = _mm_cmpeq_epi64(m0, _mm_and_si128(m0,r1));
	b1_ = _mm_cmpeq_epi64(m0_, _mm_and_si128(m0_, r1));
	b2 = _mm_cmpeq_epi64(m0, _mm_and_si128(m0,r2));
	b2_ = _mm_cmpeq_epi64(m0_, _mm_and_si128(m0_, r2));
	b3 = _mm_cmpeq_epi64(m0, _mm_and_si128(m0,r3));
	b3_ = _mm_cmpeq_epi64(m0_, _mm_and_si128(m0_, r3));
	b0 = _mm_packs_epi32(b0, b1);
	b2 = _mm_packs_epi32(b2, b3);
	b0_ = _mm_packs_epi32(b0_, b1_);
	b2_ = _mm_packs_epi32(b2_, b3_);
	__m128i dword1 = _mm_packs_epi32(b0, b2);
	__m128i dword3 = _mm_packs_epi32(b0_, b2_);

	__m128i qword0 = _mm_packs_epi16(dword0, dword1);
	__m128i qword1 = _mm_packs_epi16(dword2, dword3);
	__m128i r = _mm_packs_epi16(qword0, qword1);

	__m128i mm_ones=_mm_cmpeq_epi8(r, r);
	r = _mm_cmpeq_epi8(r, mm_ones);
	int retval = _mm_movemask_epi8(r);
//	int retval = _mm_movemask_epi8(r0);
	return retval;
#else
	int retval=0;
	for(int i=0; i<16; i++)
		retval |= (pattern_match(a,pb[i]) << i);
	return retval;
#endif
}

inline void invert(bigint& out, const bigint& in)
{
	out.len = in.len;
	int k;
	for(k=0; k<MAX_LEN; k++)
		out.n[k]=0;
	int shift = in.len & 63;
	int qw_len = (in.len+63) >> 6;
	
	if(shift != 0)
	{
		for(k=0; k<qw_len-1; k++)
		{
			uint64_t x = bit_reverse(in.n[k]);
			out.n[qw_len-2-k] |= x << shift;
			out.n[qw_len-1-k] |= x >> (64-shift);
			// bits k*64 to k*64+63
			// go into n-(k+1)*64 to n-1-k*64
		}
		uint64_t x = bit_reverse(in.n[qw_len-1]);
		out.n[0] |= x >> (64-shift);
	}
	else
	{
		for(k=0; k<qw_len; k++)
		{
			uint64_t x = bit_reverse(in.n[k]);
			out.n[qw_len-1-k] = x;
		}
	}
	out.popcnt = 0;//opcnt(out);
	for(k=0; k<MAX_LEN; k++)
		out.popcnt += __popcnt64(out.n[k]);
	/*
	if(popcnt(out)!=out.popcnt)
	{
		printf("Error: popcnt mismatch (%d %d)\n", popcnt(out), out.popcnt);
		out.print();
		exit(-1);
	}
	*/
}

/**
0 1  2      n-1 n  n+1 n+2
n n-1n-2 ...1   -  1   2
n-1n-2n-3...0   -  0
***/
inline void invert(bare_bigint& out, const bare_bigint& in, int len)
{
#if MAX_LEN==2
	if(len==0)
	{
		out.n[0]=0;
		out.n[1]=0;
		return;
	}

	__m128i inval = _mm_loadu_si128((const __m128i*)&in.n[0]);
	__m128i rev = bit_reverse(inval);
//	__m128i a = shift_right(rev, 128-len);
//	__m128i b = _mm_or_si128(a, shift_left(inval, len+1));
	
	int y = 128-len;
	__m128i a, b;
	if(128-len<=64)
	{
		int y = 128-len;
		__m128i t = _mm_srl_epi64(rev, _mm_set_epi64x(0, y));
		__m128i rem = _mm_sll_epi64(rev, _mm_set_epi64x(0, 64-y));
		a = _mm_or_si128(t, _mm_srli_si128(rem, 8));
		y = len+1;
		t = _mm_sll_epi64(inval, _mm_set_epi64x(0, y-64));
		b = _mm_slli_si128(t, 8);
	}
	else
	{
		int y = 128-len;
		__m128i t = _mm_srl_epi64(rev, _mm_set_epi64x(0, y-64));
		a = _mm_srli_si128(t, 8);
		y = len+1;
		t = _mm_sll_epi64(inval, _mm_set_epi64x(0, y));
		__m128i rem = _mm_srl_epi64(inval, _mm_set_epi64x(0, 64-y));
		b = _mm_or_si128(t, _mm_slli_si128(rem, 8));

	}
	/*
	y = len+1;
	if(y<64)
	{
		__m128i t = _mm_sll_epi64(inval, _mm_set_epi64x(0, y));
		__m128i rem = _mm_srl_epi64(inval, _mm_set_epi64x(0, 64-y));
		b = _mm_or_si128(t, _mm_slli_si128(rem, 8));
	}
	else
	{
		__m128i t = _mm_sll_epi64(inval, _mm_set_epi64x(0, y-64));
		b = _mm_slli_si128(t, 8);
	}
	*/
	b = _mm_or_si128(a, b);
	
	//uint64_t temp[2];
	_mm_storeu_si128((__m128i*)&out.n[0], b);
	/*
	_mm_storeu_si128((__m128i*)temp, b);

	int shift = len & 63;
	int qw_len = (len+63) >> 6;
	if(len > 64)
	{
		int shift = len - 64;
		out.n[0] = bit_reverse(in.n[1]); // 0 -> 127, 63 -> 64
		out.n[1] = bit_reverse(in.n[0]); // 64 -> 63, 127 -> 0
		out.n[0] = (out.n[0]>>(64-shift)) | (out.n[1] << shift);
		out.n[1] >>= 64-shift;
		//out.n[1] |= in.n[0] << shift;
	}
	else
	{
		out.n[0] = bit_reverse(in.n[0]); // 0 -> 63
		out.n[0] >>= 64-len;
		out.n[1] = 0;
		//out.n[0] |= in.n[0] << (len+1);
		//out.n[1] = in.n[0] << (len+1);
	}

	shift = (len+1) & 63;
	qw_len = (len+1) >> 6;

	if(qw_len==0)
	{
		out.n[0] |= in.n[0] << shift;
		out.n[1] |= in.n[1] << shift;
		out.n[1] |= in.n[0] >> (64-shift);
	}
	else if(qw_len==1)
	{
		out.n[1] |= in.n[0]<<shift;
	}
	if(out.n[0]!=temp[0] || out.n[1]!=temp[1])
	{
		printf("Mismatch: %I64x %I64x  vs %I64x %I64x (len %d)\n",
			out.n[0], out.n[1], temp[0], temp[1], len);
		exit(-1);
	}
	*/
#else
	int k;
	for(k=0; k<MAX_LEN; k++)
		out.n[k]=0;
	if(len==0)
		return;
	int shift = len & 63;
	int qw_len = (len+63) >> 6;
	if(shift != 0)
	{
		for(k=0; k<qw_len-1; k++)
		{
			uint64_t x = bit_reverse(in.n[k]);
			out.n[qw_len-2-k] |= x << shift;
			out.n[qw_len-1-k] |= x >> (64-shift);
			// bits k*64 to k*64+63
			// go into n-(k+1)*64 to n-1-k*64
		}
		uint64_t x = bit_reverse(in.n[qw_len-1]);
		out.n[0] |= x >> (64-shift);
	}
	else
	{
		for(k=0; k<qw_len; k++)
		{
			uint64_t x = bit_reverse(in.n[k]);
			out.n[qw_len-1-k] = x;
		}
	}
	shift = (len+1) & 63;
	qw_len = (len+1) >> 6;

	if(shift != 0)
	{
		int i;
		for(i=0; i+qw_len<MAX_LEN; i++)
		{
			out.n[i+qw_len] |= in.n[i] << shift;
			if(i+qw_len+1<MAX_LEN)
				out.n[i+qw_len+1] |= in.n[i] >> (64-shift);
		}
	}
	else
	{
		int i;
		for(i=0; i+qw_len<MAX_LEN; i++)
			out.n[i+qw_len] |= in.n[i];
	}
#endif
}

inline void invert(bigint& out, const bigint& in, int len)
{
	out.len = in.len;
	invert((bare_bigint&)out, (const bare_bigint&)in, len);
}


inline void invert2(bigint& out, const bigint& in, int len)
{
	out.len = in.len;
	int k;
	for(k=0; k<MAX_LEN; k++)
		out.n[k]=0;
	int shift = len & 63;
	int qw_len = (len+63) >> 6;
	
#if MAX_LEN==2
	if(len > 64)
	{
		int shift = len - 64;
		out.n[0] = bit_reverse(in.n[1]); // 0 -> 127, 63 -> 64
		out.n[1] = bit_reverse(in.n[0]); // 64 -> 63, 127 -> 0
		out.n[0] = (out.n[0]>>(64-shift)) | (out.n[1] << shift);
		out.n[1] >>= 64-shift;
	}
	else
	{
		out.n[0] = bit_reverse(in.n[0]); // 0 -> 63
		out.n[0] >>= 64-len;
	}
#else
	bigint out2;
	if(shift != 0)
	{
		for(k=0; k<qw_len-1; k++)
		{
			uint64_t x = bit_reverse(in.n[k]);
			out2.n[qw_len-2-k] |= x << shift;
			out2.n[qw_len-1-k] |= x >> (64-shift);
			// bits k*64 to k*64+63
			// go into n-(k+1)*64 to n-1-k*64
		}
		uint64_t x = bit_reverse(in.n[qw_len-1]);
		out2.n[0] |= x >> (64-shift);
	}
	else
	{
		for(k=0; k<qw_len; k++)
		{
			uint64_t x = bit_reverse(in.n[k]);
			out2.n[qw_len-1-k] = x;
		}
	}
	for(int i=0; i<len; i++)
		if(out.bit(i) != out2.bit(i))
		{
			printf("Mismatch in bit %d, len %d\n", i, len);
			exit(-1);
		}

#endif
}

inline void invert3(bigint& out, const bigint& in, int len)
{
	out.len = in.len;
	int k;
	for(k=0; k<MAX_LEN; k++)
		out.n[k]=0;
	if(len==0)
		return;
	int shift = len & 63;
	int qw_len = (len+63) >> 6;
#if MAX_LEN==2

	if(len > 64)
	{
		int shift = len - 64;
		out.n[0] = bit_reverse(in.n[1]); // 0 -> 127, 63 -> 64
		out.n[1] = bit_reverse(in.n[0]); // 64 -> 63, 127 -> 0
		out.n[0] = (out.n[0]>>(64-shift)) | (out.n[1] << shift);
		out.n[1] >>= 64-shift;
		//out.n[1] |= in.n[0] << shift;
	}
	else
	{
		out.n[0] = bit_reverse(in.n[0]); // 0 -> 63
		out.n[0] >>= 64-len;
		out.n[1] = 0;
		//out.n[0] |= in.n[0] << (len+1);
		//out.n[1] = in.n[0] << (len+1);
	}

	shift = (len+1) & 63;
	qw_len = (len+1) >> 6;

	if(qw_len==0)
	{
		out.n[0] |= in.n[0] >> shift;
		out.n[1] |= in.n[1] >> shift;
		out.n[0] |= in.n[1] << (64-shift);
	}
	else if(qw_len==1)
	{
		out.n[0] |= in.n[1]>>shift;
	}
#else
//	...not implemented ...
	if(shift != 0)
	{
		for(k=0; k<qw_len-1; k++)
		{
			uint64_t x = bit_reverse(in.n[k]);
			out.n[qw_len-2-k] |= x << shift;
			out.n[qw_len-1-k] |= x >> (64-shift);
			// bits k*64 to k*64+63
			// go into n-(k+1)*64 to n-1-k*64
		}
		uint64_t x = bit_reverse(in.n[qw_len-1]);
		out.n[0] |= x >> (64-shift);
	}
	else
	{
		for(k=0; k<qw_len; k++)
		{
			uint64_t x = bit_reverse(in.n[k]);
			out.n[qw_len-1-k] = x;
		}
	}
	shift = (len+1) & 63;
	qw_len = (len+1) >> 6;

	if(shift != 0)
	{
		int i;
		for(i=0; i+qw_len<MAX_LEN; i++)
		{
			out.n[i] |= in.n[i+qw_len] >> shift;
			if(i+qw_len+1<MAX_LEN)
				out.n[i] |= in.n[i+qw_len+1] << (64-shift);
		}
	}
	else
	{
		int i;
		for(i=0; i+qw_len<MAX_LEN; i++)
			out.n[i] |= in.n[i+qw_len];
	}
#endif
}

inline uint64_t rand64()
{
#ifndef WIN32
	uint64_t retval = (rand() ^ (rand()<<14) ^ (uint64_t(rand())<<28)) ^ (uint64_t(rand())<<42) ^ (uint64_t(rand())<<56);
#else
	uint64_t retval;
	while(_rdrand64_step(&retval) == 0);
#endif
	return retval;
}

