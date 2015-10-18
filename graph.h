#pragma once
#include "bigint.h"


/**
Basically a stripped down std::vector (substantially faster than the original in some situations)
**/
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
	buffer<int> vmatch;

	bigint implied_mask[2];
	bool implied_mask_set;
	/*
	bigint combo_mask[2][2*MAX_LEN*64];
	bigint* combo_ptrs[2];
	*/
	bigint explored[2];

	bool known_shifted_ref_masks;
	int shifted_ref_mask_len;
	bigint shifted_ref_masks[GRAPH_MAX_LEN*64+5];

	void merge();
	int best_next_bit;

	graph() : n(0), parent_ext(0), dirty(false), known_shifted_ref_masks(false), shifted_ref_mask_len(0), implied_mask_set(false){}

	buffer<const lite_vector*> parent_cliques[8];
	void reset();
	void print() const;
	void print_short() const;

	int search_range;
	int first_unset_bit;
	void check_first_unset_bit();
	void construct_shifted_masks(int target);
	void update_shifted_masks(int color, int target, int n);
	void forced_bit_search(int kmax[2], int target, bigint new_marks[2], bool use_parent_cliques=true, int color_mask=3, int minpos0=0, int minpos1=0);
	void forced_bit_search2(int kmax[2], int target, int color, int new_bit, bigint& new_mask, int sizePrev);
	bool new_clique_search(int kmax[2], int color, int target, bigint new_mark);
	bool set_bit(int color, int new_bit, int kmax[2], int target, bool test_only=false);

	void forced_bit_search_oneside(int kmax, int target, int color, bigint& new_mask, int minpos);
private:
	graph(const graph& g) {}
public:
	void operator=(graph& g);
};

bool build_graph_old(graph& gOut, const bigint2& mask, int kmax[2]);
bool build_graph(graph& gOut, const bigint2& mask, int kmax[2], bool do_minus2=true);

bool extend(graph& g, int color, int kmax[2], int target, int new_bit, bool full_reconst=true);
bool extend(graph& gOut, const graph& g, int color, int kmax[2], int target, int new_bit, int search_mode=1);

bool max_extensibility(graph& g, int kmax[2], int target, bool debug=false);

bool link_pair_search(int target, int kmax[2], graph& g, bigint new_masks[2], bool debug, bool use_parents);


inline void bigint_to_positions(int* positions, int& c, int len, const bigint& x)
{
	if(len>GRAPH_MAX_LEN*64)
		len=GRAPH_MAX_LEN*64;
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
