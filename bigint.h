#pragma once

#include <emmintrin.h>
#include <immintrin.h>

#define MAX_DEPTH 16
// Set to maximum graph length, in units of 64 (e.g. 6 allows graphs up to order 384).
// Must be an even number.
// Operations on shorter graphs are faster (but the code will automatically use 128-bit or 256-bit ops when possible, 
// so performance impact is relatively low)
#define GRAPH_MAX_LEN 4


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

template <int MAX_LEN>
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
	void flip(int x) { if(bit(x)) unset(x); else set(x); }

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
	void unset_to(int x)
	{
		int qw = x >> 6;
		int fract = x & 63;
		for(int i=0; i<qw; i++)
			n[i]=0;
		n[qw] &= ~((One<<fract)-One);
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

	bool zero() const;

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
		if(MAX_LEN==2)
		{
			n[0] = (n[0]>>x) | (n[1]<<(64-x));
			n[1] >>= x;
		}
		else
		{
			for(int i=0; i<MAX_LEN; i++)
			{
				if(i>0)
					n[i-1] |= n[i]<<(64-x);
				n[i] >>= x;
			}
		}
	}
	void operator|=(const bare_bigint& x)
	{
		if(MAX_LEN==2)
		{
			n[0]|=x.n[0];
			n[1]|=x.n[1];
		}
		else
		{
			for(int i=0; i<MAX_LEN; i++)
				n[i]|=x.n[i];
		}
	}
	void operator&=(const bare_bigint& x)
	{
		if(MAX_LEN==2)
		{
			n[0]&=x.n[0];
			n[1]&=x.n[1];
		}
		else
		{
			for(int i=0; i<MAX_LEN; i++)
				n[i]&=x.n[i];
		}
	}

	int trailing_bit() const
	{
		if(MAX_LEN==2)
		{
			int byte = (n[1]!=0 ? 1 : 0);
			unsigned long x;
			_BitScanReverse64(&x, n[byte]);
			return x+byte*64;
		}
		else
		{
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
		}
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


	template <int Y>
	bare_bigint& operator=(const bare_bigint<Y>& y)
	{
		int i;
		for(i=0; i<MAX_LEN && i<Y; i++)
			n[i] = y.n[i];
		for(; i<MAX_LEN; i++)
			n[i]=0;
		return *this;
	}
};


template <>
bool bare_bigint<2>::zero() const
{
	return ((n[0]|n[1])==0);
}

template <>
bool bare_bigint<4>::zero() const
{
	return ((n[0]|n[1]|n[2]|n[3])==0);
}

template<int MAX_LEN>
bool bare_bigint<MAX_LEN>::zero() const
{
	for(int k=0; k<MAX_LEN; k++)
		if(n[k]!=0)
			return false;
	return true;
}


template <int MAX_LEN>
struct ext_bigint: bare_bigint<MAX_LEN>
{
	uint64_t padding;
	int len;
	int popcnt;
	void set_len(int x) { len=x; }
	ext_bigint()
	{
//		popcnt=0;
	}
	ext_bigint(const bare_bigint& x)
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
		for(int k=0; k<MAX_LEN*64; k++)
			if(bit(k))
				printf("%d ", k);
		printf("\n");
	}

	bool operator==(const ext_bigint& x) const
	{
		if(len!=x.len)
			return false;
		for(int k=0; k<MAX_LEN; k++)
			if(n[k]!=x.n[k])
				return false;
		return true;
	}
	bool operator<(const ext_bigint& x) const;
	void clear()
	{
		len=popcnt=0;
		for(int k=0; k<MAX_LEN; k++)
			n[k]=0;
	}

	template <int Y>
	ext_bigint& operator=(const ext_bigint<Y>& y)
	{

		clear();
		set_len(y.len);
		popcnt = y.popcnt;
		for(int i=0; i<MAX_LEN && i<Y; i++)
			n[i] = y.n[i];
		return *this;
	}
	void flip2(int x) 
	{
		flip(x);
		if(len-1-x!=x)
			flip(len-1-x);
	}
};

template <int T>
bool ext_bigint<T>::operator<(const ext_bigint<T>& x) const
{
	if(T==6)
	{
		if(n[0]!=x.n[0])
			return n[0]<x.n[0];
		if(n[1]!=x.n[1])
			return n[1]<x.n[1];
		if(n[2]!=x.n[2])
			return n[2]<x.n[2];
		if(n[3]!=x.n[3])
			return n[3]<x.n[3];
		if(n[4]!=x.n[4])
			return n[4]<x.n[4];
		return n[5]<x.n[5];
	}
	else
	{
		exit(-1);
		return false;
	}
}

template <>
bool ext_bigint<2>::operator<(const ext_bigint<2>& x) const
{
	if(n[0]!=x.n[0])
		return n[0]<x.n[0];
	return n[1]<x.n[1];
}

template <>
bool ext_bigint<4>::operator<(const ext_bigint<4>& x) const
{
	if(n[0]!=x.n[0])
		return n[0]<x.n[0];
	if(n[1]!=x.n[1])
		return n[1]<x.n[1];
	if(n[2]!=x.n[2])
		return n[2]<x.n[2];
	return n[3]<x.n[3];
}

template<int T>
inline void store_and(bare_bigint<T>& dst, const bare_bigint<T>& src_a, const bare_bigint<T>& src_b)
{
	for(int i=0; i<T; i++)
		dst.n[i] = src_a.n[i] & src_b.n[i];
}

template<>
inline void store_and(bare_bigint<2>& dst, const bare_bigint<2>& src_a, const bare_bigint<2>& src_b)
{
	__m128i a = _mm_load_si128((const __m128i*)&src_a.n[0]);
	__m128i b = _mm_load_si128((const __m128i*)&src_b.n[0]);
	_mm_store_si128((__m128i*)(&dst.n[0]), _mm_and_si128(a, b));
}


template <int T>
inline bool bit_equal(const bare_bigint<T>& a, const bare_bigint<T>& b)
{
	int i;
	for(i=0; i<T; i++)
		if(a.n[i]!=b.n[i])
			return false;
	return true;
}

template <>
inline bool bit_equal(const bare_bigint<2>& a, const bare_bigint<2>& b)
{
	return (a.n[0]==b.n[0] && a.n[1]==b.n[1]);
}


template <int T>
inline bool bit_equal(const ext_bigint<T>& a, const ext_bigint<T>& b)
{
	int i;
	for(i=0; i<T; i++)
		if(a.n[i]!=b.n[i])
			return false;
	return true;
}

template <>
inline bool bit_equal(const ext_bigint<2>& a, const ext_bigint<2>& b)
{
	return (a.n[0]==b.n[0] && a.n[1]==b.n[1]);
}


typedef ext_bigint<GRAPH_MAX_LEN> bigint;

typedef pair<bigint,bigint> bigint2;

typedef ext_bigint<2> b128;


template<int T>
inline uint64_t popcnt(const bare_bigint<T>& a)
{
	uint64_t sum=0;
	for(int k=0; k<T; k++)
		sum += __popcnt64(a.n[k]);
	return sum;
}

template <int T>
inline bare_bigint<T> andnot(const bare_bigint<T>& a, const bare_bigint<T>& b)
{
	bare_bigint<T> x;
	for(int k=0; k<T; k++)
		x.n[k] = a.n[k] & ~b.n[k];
	return x;
}

template <int T>
inline uint64_t popcnt(const bare_bigint<T>& a, const bare_bigint<T>& b)
{
	uint64_t sum=0;
	if(T==1)
		sum = __popcnt64(a.n[0] & b.n[0]);
	else if(T==2)
		sum = __popcnt64(a.n[0] & b.n[0])
			+ __popcnt64(a.n[1] & b.n[1]);
	else if(T==3)
		sum = __popcnt64(a.n[0] & b.n[0])
			+ __popcnt64(a.n[1] & b.n[1])
			+ __popcnt64(a.n[2] & b.n[2]);
	else
		for(int k=0; k<T; k++)
			sum += __popcnt64(a.n[k] & b.n[k]);
	return sum;
}


template<int MAX_LEN>
inline ext_bigint<MAX_LEN> operator|(ext_bigint<MAX_LEN> a, ext_bigint<MAX_LEN> b)
{
	for(int i=0; i<MAX_LEN; i++)
		a.n[i] |= b.n[i];
	return a;
}

template<int MAX_LEN>
inline ext_bigint<MAX_LEN> operator^(ext_bigint<MAX_LEN> a, ext_bigint<MAX_LEN> b)
{
	for(int i=0; i<MAX_LEN; i++)
		a.n[i] ^= b.n[i];
	return a;
}


template<int MAX_LEN>
inline ext_bigint<MAX_LEN> operator&(ext_bigint<MAX_LEN> a, ext_bigint<MAX_LEN> b)
{
	for(int i=0; i<MAX_LEN; i++)
		a.n[i] &= b.n[i];
	return a;
}

template<>
inline ext_bigint<2> operator|(ext_bigint<2> a, ext_bigint<2> b)
{
	a.n[0]|=b.n[0];
	a.n[1]|=b.n[1];
	return a;
}

template<>
inline ext_bigint<2> operator^(ext_bigint<2> a, ext_bigint<2> b)
{
	a.n[0]^=b.n[0];
	a.n[1]^=b.n[1];
	return a;
}

template<>
inline ext_bigint<2> operator&(ext_bigint<2> a, ext_bigint<2> b)
{
	a.n[0]&=b.n[0];
	a.n[1]&=b.n[1];
	return a;
}

template<int MAX_LEN>
inline ext_bigint<MAX_LEN> operator<<(ext_bigint<MAX_LEN> x, int n)
{
	ext_bigint<MAX_LEN> retval;
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
	return retval;
}


template<>
inline ext_bigint<2> operator<<(ext_bigint<2> x, int n)
{
	ext_bigint<2> retval;
	retval.set_len(x.len);
	retval.popcnt = x.popcnt;
	__m128i m = _mm_loadu_si128((const __m128i*)&x.n[0]);
	m = shift_left(m, n);
	_mm_storeu_si128((__m128i*)&retval.n[0], m);
	return retval;
}


template<int MAX_LEN>
inline ext_bigint<MAX_LEN> operator>>(ext_bigint<MAX_LEN> x, int n)
{
	ext_bigint<MAX_LEN> retval;
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

template<int MAX_LEN1, int MAX_LEN2>
inline void invert(ext_bigint<MAX_LEN1>& out, const ext_bigint<MAX_LEN2>& in)
{
	out.len = in.len;
	int k;
	for(k=0; k<MAX_LEN1; k++)
		out.n[k]=0;
	int shift = in.len & 63;
	int qw_len = (in.len+63) >> 6;
	if(qw_len>MAX_LEN1)
	{
		printf("Error: invert() overflow\n");
		exit(-1);
	}
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
	for(k=0; k<MAX_LEN1; k++)
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
template<int MAX_LEN1, int MAX_LEN2>
inline void invert(bare_bigint<MAX_LEN1>& out, const bare_bigint<MAX_LEN2>& in, int len)
{
	if(MAX_LEN1==2)
	{
		if(len==0)
		{
			//out = in;
			//out >>= 1;
			out.clear();
			return;
		}
		if(len>128)
		{
			printf("invert() overflow\n");
			exit(-1);
		}

		__m128i inval = _mm_loadu_si128((const __m128i*)&in.n[0]);
		__m128i rev = bit_reverse(inval);
	
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
		b = _mm_or_si128(a, b);
		_mm_storeu_si128((__m128i*)&out.n[0], b);
	}
	else
	{
		out.clear();
		if(len==0)
		{
			//out = in;
			//out >>= 1;
			return;
		}
		if(len<=128)
		{
			__m128i inval = _mm_loadu_si128((const __m128i*)&in.n[0]);
			__m128i rev = bit_reverse(inval);

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
			b = _mm_or_si128(a, b);
			_mm_storeu_si128((__m128i*)&out.n[0], b);
			return;
		}

		int k;
		for(k=0; k<MAX_LEN1; k++)
			out.n[k]=0;
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
			for(i=0; i+qw_len<MAX_LEN1; i++)
			{
				out.n[i+qw_len] |= in.n[i] << shift;
				if(i+qw_len+1<MAX_LEN1)
					out.n[i+qw_len+1] |= in.n[i] >> (64-shift);
			}
		}
		else
		{
			int i;
			for(i=0; i+qw_len<MAX_LEN1; i++)
				out.n[i+qw_len] |= in.n[i];
		}
	}
}

template<int T>
inline void invert(ext_bigint<T>& out, const ext_bigint<T>& in, int len)
{
	out.len = in.len;
	invert((bare_bigint<T>&)out, (const bare_bigint<T>&)in, len);
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

