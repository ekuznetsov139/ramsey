#pragma once
#include "bigint.h"

template <class T>
struct buffer
{
	T* ptr;
	size_t dim, cap;
	buffer() { ptr=(T*)malloc(1024*sizeof(T)); dim=0; cap=1024; }
	void resize(size_t n, bool exact=false) 
	{
		if(n>cap) 
		{
			size_t new_cap=exact ? n : (n+n/2);
			ptr=(T*)realloc(ptr,new_cap*sizeof(T));
			if(ptr==0)
			{
				printf("ERROR: realloc(size=%I64d*%d) failure\n", new_cap, sizeof(T));
				exit(-1);
			}
			cap=new_cap;
		}
		dim=n;
	}
	void push_back(const T& x)
	{
		resize(dim+1);
		ptr[dim-1]=x;
	}
	T& operator[](size_t x) { return ptr[x]; }
	const T& operator[](size_t x) const { return ptr[x]; }
	size_t size() const { return dim; }
	~buffer() { free(ptr); }
	void swap(buffer<T>& x)
	{
		std::swap(ptr, x.ptr);
		std::swap(dim, x.dim);
		std::swap(cap, x.cap);
	}
	void alloc(size_t x)
	{
		if(x < cap)
			return;
		free(ptr);
		ptr = (T*)malloc(x*sizeof(T));
		cap=x;
	}
	void clear() { dim=0; }
	void operator=(const buffer& buf) 
	{
		resize(buf.size());
		memcpy(ptr, buf.ptr, buf.size()*sizeof(T)); 
	}
private:
	buffer(const buffer& buf) { printf("Copy construction not allowed\n"); }
};



typedef buffer<bigint> lite_vector;
typedef buffer<bigint2> lite_vector2;

struct graph
{
	int n;
	bigint mask0;
	bigint mask;
	bool dirty;
	lite_vector cliques[8];
	int clique_counts[2][MAX_DEPTH];
	const graph* parent_ext;
	buffer<int> new_marks[2];
	buffer<pair<int,int> > new_bits;
	buffer<int> vmatch;
#if MAX_LEN==2
	bigint extensions[4][MAX_LEN*64];
#endif
	bool set_extensions;
	int max_set_ext;
	bigint combo_mask[2][2*MAX_LEN*64];
	bigint* combo_ptrs[2];

	void merge();

	float save_rates[64];
	int save_base;

	int best_next_bit;
	uint64_t v_a0, v_b0, v_a1, v_b1;
	int second_best_next_bit;

	graph() : n(0), set_extensions(false), parent_ext(0), dirty(false) {}

	buffer<const lite_vector*> parent_cliques[8];
	void reset();
	void print() const;
	void print_short() const;

	int search_range;
	int first_unset_bit;
	void check_first_unset_bit();
	void construct_shifted_masks(int target);
	void update_shifted_masks(int color, int target, int n);
	void forced_bit_search(int kmax[2], int target, buffer<int> new_marks[2], bool use_parent_cliques=true, int color_mask=3, int minpos0=0, int minpos1=0, bool debug=false);
	void forced_bit_search2(int kmax[2], int target, int color, int new_bit, buffer<int> new_marks[2], int sizePrev, bool debug=false);
	int new_clique_search(int kmax[2], int target, buffer<int> new_marks[2], bigint* badmask=0);
	bool set_bit(int color, int new_bit, int kmax[2], bool exact=false, bigint* badmask=0);
	bool set_bits(int color, buffer<int>& new_bits, int kmax[2], int target, bigint* badmask=0);
private:
	graph(const graph& g) {}
public:
	void operator=(graph& g);
};

bool build_graph2(graph& gOut, const bigint2& mask, int kmax[2]);
bool build_graph3(graph& gOut, const bigint2& mask, int kmax[2]);

bool extend(graph& g, int color, int kmax[2], int target, int new_bit);
bool extend(graph& gOut, const graph& g, int color, int kmax[2], int target, int new_bit, int search_mode=1, bool debug = false, int force_range = 0);

int max_extensibility(graph& g, int kmax[2], int target, bool debug=false);

bool link_pair_search(int target, int kmax[2], bigint* combo_ptrs[2], graph& g, buffer<pair<int,int> >& new_bits, bool debug=false, bool use_parents=false);


inline void bigint_to_positions(int* positions, int& c, int len, const bigint& x)
{
	if(len>MAX_LEN*64)
		len=MAX_LEN*64;
	memset(positions, 0, MAX_DEPTH*sizeof(int));
	int m;
	c=0;
	for(m=0; m<len; m++)
	{
		if(x.bit(m))
		{
			positions[c]=m+1;
			c++;
		}
	}
}

bool verify_clique(const bigint& x, int target, const bigint& mask);
void verify_extension(const bigint& x, int n, int target, const bigint& mask);
void verify_extension(const bigint& x, int gap, int n, int target, const bigint& mask);
void verify_extension_2x(const bigint& x, int n1, int n2, int target, const bigint& mask);
void verify_extension_2x(const bigint& x, int idx1, int idx2, int n1, int n2, int target, const bigint& mask);

void validate_graph_build_native(const graph& g, const bigint2& mask);


inline void validate_graph_build(const graph& g, const bigint2& mask)
{
	return;
	//validate_graph_build_native(g, mask);
}

inline void validate_graph_build(const graph& g)
{
	validate_graph_build(g, make_pair(g.mask0, g.mask));
}