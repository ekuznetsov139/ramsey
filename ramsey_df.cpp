#include "stdafx.h"
#include "bigint.h"
#include "graph.h"
#include <omp.h>
#define USE_TBB
#ifdef WIN32
typedef unsigned char uint8_t;
#include <Windows.h> // InterlockedIncrement
#endif
#include <set>
#include <map>
#include <tbb/parallel_sort.h>
#include <concurrent_vector.h>
using namespace tbb;
using namespace Concurrency;


#ifndef WIN32
#define _ftelli64 ftell
#define _unlink unlink
#endif

int g_kmax[2]={5,5};

double g_fFreq;
int nthreads;
const int max_nthreads = 48;
#ifdef WIN32
#define PATH "c:\\temp\\"
#else
#define PATH
#endif

#if GRAPH_MAX_LEN==2
#define graph_pool_dim 160
#elif GRAPH_MAX_LEN==4
#define graph_pool_dim 280
#else
#define graph_pool_dim 400
#endif

#define pop_acc_dim (GRAPH_MAX_LEN*64+2)

/*** Main top-level methods ***/

// Top-level entry point of the signature-based search. 
// Searches for R(k0+2,k1+2) 's'-bit signatures extending to colorings of order 'd', while automatically locating and recording
// cyclic colorings.
// Infinite loop combining nearest neighbor search and relabeling. Simply pick parameters, provide some starting data 
// (anything either in the hit file or in the sym file), and leave it running. Terminates when it runs out of things to do.
// General guidelines:
// 'tests_per_pass' controls the duration of each pass. Keep high enough to avoid wasting CPU time (aim for at least 1 min / pass)
// 'max_ops' controls depth search timeout (see paper). Ideal range is 10^4 to 10^5.
// If in doubt, look at "Self check 4" output. If it contains "1000 -> x" with x<980, raise max_ops. If x>995, consider lowering.
// If max_ops<10^4 and x>995, lower 's'.
// If max_ops>10^5 and x<980, increase 's'.
// It'll still run outside these boundaries, but may be inefficient or miss valid signatures.
// Current implementation does not allow 's' to exceed 127. 'd' can be arbitrarily large, but you must set GRAPH_MAX_LEN in bigint.h so 
// that d<GRAPH_MAX_LEN*64, and, due to the limit on 's', it will likely become quite inefficient if d>300.
// If 'self_checks' is 'true', the function will spot check data in list files (or, equivalently, check that it works correctly on known 
// data) before getting to work. Generally safe to leave it on.
void unified_search(int k0, int k1, int s, int d, int tests_per_pass, int max_ops, bool self_checks=true);

// Top-level entry point of the cyclic-based search.
// if 'fwd_only' is false, explores neighbors of cyclic colorings of order 'd' or above in the sym list by doing up to 2 flips or 1 flip + reflection.
// if 'fwd_only' is true, explores longer neighbors of colorings of order 'd' by doing reflection and up to 1 flip.
// Keeps track of previously processed colorings.
void cyclic_network_search(int k0, int k1, int d, bool fwd_only=false, int min_log=0, int max_passes=INT_MAX);

// Cyclic-based ladder search.
// Starting with 1 or more short R(k0+2,k1+2) cyclic colorings, tries to find colorings of order 'min_log_len' or larger.
// Will write any colorings of sufficient size into the corresponding sym file.
// Typical run time ~1..20 minutes for each cluster in init_sym
void cyclic_ladder_search(int k0, int k1, int start_length, vector<bigint>& init_sym, int min_log_len);

// Wrapper for the function above that combines a Monte Carlo run (to generate initial cyclic colorings) and an actual ladder search.
// Parameters (dim, start_size, mc_test_count) need to be set on case-by-case basis. The goal is to get a reasonable number of initial cyclic clusters
// (more than 1 and fewer than 50) in reasonable amount of time (ladder search is relatively insensitive to the start level, so the goal is
// for MC time to be <10% of total run time.)
// Typically start_size is in the low 100s (e.g. 120 for (7,8), 140 for (8,8), 140 for (6,12), 150 for (7,10))
// Set dim=start_size/2 or a bit lower
// mc_test_count in 10^4 to 10^6 range 
// As a general rule, the number of clusters generated per unit of time declines by a factor of 5 every time start_size is incremented by 5.
// If MC generates some signatures but not enough cyclics, it will keep running nearest neighbor search till it gets at least 5 cyclics (generally 
// undesirable since resulting cyclics likely won't be independent.)
// If mc_test_count is set too low or 'start_size' is set too high, MC won't find any hits at all, in which case the function returns right away.
void cyclic_ladder_search(int dim, int k0, int k1, int start_size, int mc_test_count, int min_log_len);


/*** General maintenance ***/

//List reading/writing

// List file name generators. Return values depend on current values of g_kmax[0]/g_kmax[1]
// E.g., in Windows, for R(8,8), s=80, d=180 (that is, 80-bit signatures corresponding to order 81 graphs), file names would be
// 'c:\\temp\\88-81-180.txt' 
// 'c:\\temp\\cover88-81-180.txt'
// 'c:\\temp\\sym88.txt'
string format_name(int len, int target, string prefix=string()); // Hit file (contains 'len' bit sigs which are known to extend to 'target')
string format_sym_name(string prefix=""); // Sym file (contains known cyclic colorings)
string format_cover_name(int len, int target); // Cover file (contains sigs which were checked at some prior point. Would normally contain the data in the hit file as a subset.

// List file readers/writers.
// For historical reasons, file name generators produce strings with '.txt' extensions, and functions below can still read txt files,
// but, except for the sym file, they've been converted to use binary files and they append '.bin' to file names internally,
// since reading binary is substantially faster than parsing text, and some jobs result in extremely large files
// (e.g. cover files exceeding 1 GB are common and 5 GB was exceeded on a few occasions)
// 'len' argument is now obsolete but should be set equal to 'len' passed to format_name()/format_cover_name()
template <int T> void readfile(string fn, vector<ext_bigint<T> >& sigs, int len);
template <int T> void readfile(const char* fn, vector<ext_bigint<T> >& sigs, int len);
template <int T> void writefile(string str, vector<ext_bigint<T> >& sigs, int len, bool append=false);
template<int T> void writefile(const char* fn, vector<ext_bigint<T> >& sigs, int len, bool append=false);
void readsym(string s, vector<bigint>& v, int target=0);
template <class T> void writesym(string s, T& v);
template <class T> void writesym(T it0, T it1, const char* filename, int min_len=0);

// Cleans the sym file by removing duplicates, removing relabelings, optionally checking for errors, and optionally merging in a secondary sym file 
void cleansym(int k0, int k1, bool check_errors=false, const char* secondary = 0, int sec_min_len=0);

// 'Fixes' the corresponding signature list pair, by removing duplicates and ensuring that all 'hit' entries are present in the 'cover' file
void fix_list(int d, int k0, int k1, int target);

/**
Core of the signature based search. Walks through the graph tree rooted at 'sig' attempting to construct a complete valid graph of order 'target'.
'detect_only' should nearly always be set to 'true', and 'stop_point' to target+1.
Normally won't be called directly - use JobArray instead.

Return code:
* 0: no valid graphs were found (the tree was exhausted or the search timed out because g_max_ops was nonzero and the number of ops was exceeded)
* 1: at least one valid graph was found. Observed graphs are inserted into 'pHits'. If 'detect_only' is set, the call returns 1 with a single entry in 'pHits' as soon as the first valid graph is found.
* 2: bad signature, no need to run the search.
**/
int depth_search(ext_bigint<2> sig, int kmax[2], int target, int stop_point, bool detect_only, bool quiet, int* populations, graph* graph_pool, int thr, buffer<bigint>* pHits=0);

/**
Parallelizer for depth_search calls.
Create an instance and call submit(), passing an array of ext_bigint<2>'s. 
Reports results via callbacks. Callbacks calls are controlled by a mutex (it's safe to e.g. update global container objects),
but they should be lightweight or they'll interfere with parallelization.
It's possible but generally not recommended for callbacks to write into files (it's especially "wrong" for a callback to call fopen/fclose each time,
since that can drastically hurt performance.)
**/
typedef void hit_detect_callback_fun(uint64_t pos, ext_bigint<2> sig, buffer<bigint>* pHits, int max_len);
typedef void base_cover_callback_fun(ext_bigint<2> sig, uint64_t max_len);

struct depth_search_state
{
	depth_search_state() : pHitMiss(0) {}
	int* populations;
	graph* graph_pool;
	uint64_t base;
	uint64_t len;
	uint64_t fractions;
	bool detect_only;
	int target;
	int stop_size;
	int thr;
	omp_lock_t* pCritSect;
	unsigned int* pNextPos;
	uint64_t cputime;
	const b128* bases;
	uint64_t* base_masks;
	uint64_t nbases;
	int nthr;
	uint64_t t0;
	hit_detect_callback_fun* pHitDetect;
	hit_detect_callback_fun* pHitMiss;
	base_cover_callback_fun* pBaseCover;
};

struct JobArray
{
	omp_lock_t critSect;
	depth_search_state jobs[max_nthreads];
	int* populations;
	uint64_t sum_populations[pop_acc_dim];
	graph* graph_pool;
	uint64_t cputime;
	uint64_t realtime;

	JobArray() : populations(0), graph_pool(0), cputime(0), realtime(0)
	{
		omp_init_lock(&critSect);
		populations=new int[nthreads*pop_acc_dim];
	}
	~JobArray()
	{
		omp_destroy_lock(&critSect);
		delete[] populations;
	}
	void submit(const ext_bigint<2>* p, size_t n, int len, int target, int stop_size, bool detect_only, int fractions=1, 
		base_cover_callback_fun* pBaseCover=0, hit_detect_callback_fun* pHitDetect=0, 
		hit_detect_callback_fun* pHitFail=0);
	void print_populations();
	void print_cputime();
};

/** Signature based search - higher level operations **/

// Loads the sym list into memory. Should be called before any calls to lengthen/narrow/in_memory_search to prevent them
// from creating duplicates in the sym file.
void load_known_sym_graphs(int k0, int k1, int d);


// Merges a (dPrev,targetPrev) list into a (d,target) list (with d<=dPrev and/or target<=targetPrev). 
void merge(int k0, int k1, int dPrev, int d, int targetPrev, int target);

// Full relabeling. Extremely slow.
// 'initial'=false: only loads cyclics from the sym-file, generates signatures from them, and writes them into the hit file.
// 'initial'=true: also does the actual work (goes through all signatures looking for cyclics).
void do_relabeling(int d, int target, int k0, int k1, bool initial, double fraction=1.0);

// Shortcut relabeling. Fast. Runs through all elements in the list, imputes possible cyclics by inspecting
// 'target' order colorings returned by depth_search(), and tests them. Does not update the sig file.
// Typically unnecessary because most search methods do this implicitly.
void do_relabeling_new(int d, int target, int k0, int k1, double fraction=1.0);

// Full enumeration of R(k0+2,k1+2) colorings.
// 'dim' controls search space segmentation (internally, this will be executed as 2^dim calls to depth_search().)
// If the projected run time is <1 min, set dim=10.
// For 1..100 min, dim=15.
// For 100+ min, dim=20.
void do_full_search(int k0, int k1, int dim, int target);

// "Lengthen": take a (s,d) list and extend it into a (s',d') list with s'>=s and/or d'>=d.
// This entails creating 2^(s'-s) signatures for each list entry and checking each one for extensibility to d'.
// File-based versions don't carry over existing (s,d) cover since that would be prohibitively expensive in terms of storage space.
// File-based versions will merge the output into the existing (s',d') list if it's there, and exclude entries in the existing (s',d') cover from testing.
// Don't use when s=s' (use "narrow()" instead)
int lengthen(const char* fn, int k0, int k1, int len, int lenNext, int target, double fraction = 1.0);
int lengthen(int k0, int k1, int len, int lenNext, int targetPrev, int target, double fraction = 1.0);
int lengthen_memory(const vector<b128>& v, vector<b128>& out, int k0, int k1, int len, int lenNext, int target);

// An alternate implementation of 'lengthen'. Instead of checking 2^(s'-s) (s',d) signatures per entry, checks only one (s,d) signature and pulls extra bits
// from the bigint returned by depth_search.
// Consequently, g_max_ops should be set as if this were a (s,d') run.
// If 'detect_only' is 'true', it will report one hit per input list entry (occasionally even slightly less if some runs timeout)
// If 'detect_only' is 'false', the results are, in theory, identical to normal "lengthen" (slightly more than 1 hit/entry).
// In practice, with 'detect_only'='false', this function is so slow that it's not evident that it's better than normal "lengthen".
int lengthen_fast(const char* fn, int k0, int k1, int len, int lenNext, int target, double fraction = 1.0, bool detect_only=false);
int lengthen_fast(int k0, int k1, int len, int lenNext, int targetPrev, int target, double fraction = 1.0, bool detect_only=false);

// Take a (s,d) list and extend it into a (s,d') list with d'>d. Carries over and merges the (s,d) cover into the (s,d') cover.
// file-based
int narrow(int len, int targetPrev, int target, double fraction = 1.0);
// memory-based
int narrow_memory(vector<b128>& vsig, int len, int targetPrev, int target, vector<b128>& results, vector<bigint>& sym, int min_save);

// In-memory nearest-neighbor search. 
// 'sparse'=false: clusterize the input list and keep running nearest-neighbor search until all clusters are either complete (with no unexplored neighbors) or are over 'cap' elements.
// 'sparse'=true: clusterize the input list and run nearest-neighbor search starting from these clusters, forming "trails" quickly leading "away" from original clusters.
void in_memory_search(int d, int k0, int k1, int target0, int target1, vector<b128>& v, vector<b128>& cover, vector<bigint>& sym, int min_save, bool sparse = true, int cap = 80);



/** Cyclic based search **/

// Thin wrappers for build_graph() 
bool test_build(graph& g, bigint b, int k0, int k1);
bool test_build(graph& g, const graph& gRef, bigint b, int k0, int k1);

// Obsolete. Searches for cyclic colorings differing by up to 2 bits from colorings in the known sym list.
void horizontal_search(int k0, int k1, int d, int dmax=INT_MAX, int max_tests=INT_MAX);

/**************** Implementation **************************/

// Graphs of all levels and threads are pre-allocated and kept in memory simultaneously.
// Underlying object never releases any memory (realloc()'s coming from multiple threads at the same time
// are very bad for performance)
graph* g_graph_pool = 0;

#define PARALLEL



/***
Controls direct vs. indirect pre-vet mode in cyclic search.

Pre-vet involves listing a subset of 1-incomplete maximal cliques in a cyclic coloring and quickly checking them against each candidate generated from this coloring.
If any of them is complete in the new coloring, it can be discarded and a build_graph() can be skipped.
Quite often it's possible to pick a relatively small number (10^3 to 10^4) of 1-incompletes (fewer than 10% of all 1-incompletes) and effect a 90-95% reduction in the 
number of candidates that reach build_graph() stage.

In direct mode, build_graph() is called for each relabeling and generated cliques are used directly to make pre-vet masks (see ::do_minus_ones())
In indirect mode, build_graph() is called for the original cyclic only, and generated cliques go through relabeling.

Indirect mode might seem like a no-brainer (all relabelings contain the same cliques, just node labels are different, right?)
but it's quite a bit trickier than that.
* build_graph() creates a "nicely" ordered clique list, with ones most valuable as pre-vet masks naturally occurring near the front. Attempting to relabel shuffles cliques.
* Getting the right logic to pick pre-vet masks out of all 1-incomplete cliques is easier if cliques are nicely ordered.
* Relabeling isn't all that cheap, and the number of submaximal cliques can be huge. It's sometimes cheaper to do build_graph() than to relabel all of them.

As implemented, indirect is _generally_ somewhat faster, with exceptions (this varies from 3x speedup for (6,12) and (5,13), to about even for (8,8), to 50% slower for (7,8).)

TL;DR: keep it off.
****/
//#define DIRECT_PRE_VET_MODE

template <class X>
inline void append(X& v, const X& w)
{
	if(w.size()==0)
		return;
	size_t n = v.size();
	v.resize(n+w.size());
	memcpy(&v[n], &w[0], w.size()*sizeof(w[0]));
}

template <int T>
inline void neg(vector<ext_bigint<T> >& v)
{
	for (size_t i = 0; i < v.size(); i++)
	{
		v[i].neg();
		v[i].unset_from(v[i].len);
	}
}


template <class T>
void refilter(vector<T>& v, int d=2)
{
	vector<T> v2;
	for (int i = 0; i < v.size(); i++)
		if (v[i].len + 1 >= d)
			v2.push_back(v[i]);
	v.swap(v2);
}


template <class T>
void shuffle(vector<T>& v)
{
	if(v.size()<nthreads*4)	
		return;
	size_t dim = v.size() / nthreads;
#pragma omp parallel for
	for(int thr=0; thr<nthreads; thr++)
	{
		for(size_t i=0; i+1<dim; i+=2)
		{
			size_t r = rand64();
			size_t p = (r & 0xFFFFFFFF) % (dim-1);
			if(p==i)
				p++;
			swap(v[nthreads*i+thr], v[nthreads*p+thr]);

			p = (r>>32) % (dim-1);
			if(p==i+1)
				p++;
			swap(v[nthreads*(i+1)+thr], v[nthreads*p+thr]);
		}
	}
}

bool symmetry_check(graph& g, int kmax[2], int target, bool debug);

inline bool operator<(const pair<uint64_t,int>& a, const pair<uint64_t,int>& b) { return a<b; }

template <class T>
inline size_t rem_dup(vector<T>& v)
{
	if(v.empty())
		return 0;
#ifdef USE_TBB
	tbb::parallel_sort(v.begin(), v.end());
#else
	sort(v.begin(), v.end());
#endif
	size_t sizePrev = v.size();
	v.resize(unique(v.begin(), v.end())-v.begin());
	return sizePrev-v.size();
}

void prealloc_graph_pool(graph*& graph_pool, int len, int min_entries=0)
{
	if(g_graph_pool!=0)
	{
		graph_pool = g_graph_pool;
		return;
	}
	min_entries=max(min_entries, nthreads*graph_pool_dim*2);
	g_graph_pool=new graph[min_entries];
	graph_pool = g_graph_pool;
//	printf("Attempting to alloc %
	for(int i=len-1; i<graph_pool_dim; i++)
	{
		int expect_size0 = 2000*i/50;
		if(expect_size0<1000)
			expect_size0=1000;
		int expect_size1 = expect_size0;
		int def_size = expect_size0;
		if(g_kmax[0]>=7)
			expect_size0=100000;
		if(g_kmax[1]>=7)
			expect_size1=100000;
		if(i>300)
			expect_size0=expect_size1=1000;
		if(g_kmax[0]+g_kmax[1]>=11 && i<300)
		{
			if(g_kmax[0]>g_kmax[1])
			{
				expect_size0=(g_kmax[0]>=8 ? 800*i : 500*i);
				expect_size1=300*i;
			}
			else
			{
				expect_size0=300*i;
				expect_size1=(g_kmax[1]>=8 ? 800*i : 500*i);
			}
		}

		
//		printf("%d\t%d\t%d\n", i, expect_size, expect_size*nthreads);
		for(int thr=0; thr<nthreads; thr++)
		{
			for(int cl=0; cl<6; cl++)
			{
				graph_pool[thr*graph_pool_dim*2+i].cliques[cl].alloc(cl>=2 ? def_size : ((cl&1) ? expect_size1 : expect_size0));
				if(cl<2)
					graph_pool[thr*graph_pool_dim*2+i+graph_pool_dim].cliques[cl].alloc(cl>=2 ? def_size : ((cl&1) ? expect_size1 : expect_size0));
			}
		}
	}
}

inline bool extend_from(graph& g, const graph& gmin, bigint2 m, int kmax[2], int target, bool full_reconst = true)
{
	int last_bit=max(m.first.trailing_bit(), m.second.trailing_bit());
	bool first_bit=true;
	int last_known_bit=max(gmin.mask.trailing_bit(), gmin.mask0.trailing_bit());
	if(last_known_bit >= last_bit)
	{
		printf("%d %d\n", last_known_bit, last_bit);
		exit(-1);
	}
	if (!andnot(gmin.mask0, m.first).zero())
	{
		printf("extend_from internal error\n");
		exit(-1);
	}
	if (!andnot(gmin.mask, m.second).zero())
	{
		printf("extend_from internal error\n");
		exit(-1);
	}
	for (int i = last_known_bit + 1; i <= last_bit; i++)
	{
		int color;
		if(m.first.bit(i))
			color=0;
		else if(m.second.bit(i))
			color=1;
		else
			continue;
		if(first_bit)
		{
			first_bit=false;
			if(!extend(g, gmin, color, kmax, target, i, full_reconst ? 0 : 3))
				return false;
		}
		else
		{
			if (!extend(g, color, kmax, target, i, full_reconst))
				return false;
		}
	}
	if(first_bit)
	{
		printf("extend_from error\n");
		exit(-1);
	}
	return true;
}


bool validate_solution(graph& g, int kmax[2], int target, graph* gmin, int pos)
{
	bigint2 m = make_pair(g.mask0, g.mask);
	int first_unset=0;
	while((g.mask0.bit(first_unset) || g.mask.bit(first_unset)) && first_unset<=GRAPH_MAX_LEN*64)
		first_unset++;

	if(first_unset > target+3)
		return true;
	/*
	for(int i=first_unset; i<MAX_LEN*64; i++)
	{
		m.first.unset(i);
		m.second.unset(i);
	}
	*/
	int last_known=pos;
	while(gmin[last_known].n==0 && last_known>=0)
		last_known--;
	if(last_known<0)
	{
		printf("Error: last_known<0\n");
		exit(-1);
	}
	/*
	if(popcnt(gmin[last_known].mask)+popcnt(gmin[last_known].mask0)+20
		< popcnt(m.first)+popcnt(m.second))
	{
		printf("Warning: extending to %d from %d: %d vs %d bits set\n",
			target, last_known, 
			popcnt(gmin[last_known].mask)+popcnt(gmin[last_known].mask0),
			popcnt(m.first)+popcnt(m.second)
			);
	}
	*/
	if(popcnt(gmin[last_known].mask)+popcnt(gmin[last_known].mask0)
		< popcnt(m.first)+popcnt(m.second))
	{
		uint64_t t0, t1;
		gmin[0].reset();
		t0 = __rdtsc();
		bool ok = extend_from(gmin[0], gmin[last_known], m, kmax, target + 2);
		t1 = __rdtsc();
		if (!ok)
			return false;
		if(gmin[0].clique_counts[0][kmax[0]]>0 || gmin[0].clique_counts[1][kmax[1]]>0)
			return false;
	}

	return true;
}

vector<int> op_counts[3];
omp_lock_t op_count_lock;

int g_max_ops=0;
int g_max_ops_2=0;
bool g_search_trace = false;

void strip(graph& g, int* kmax)
{
	for(int i=0; i<2; i++)
	{
		int readpos=0, writepos=0;
		while(readpos<g.cliques[i].size())
		{
			if(g.cliques[i][readpos].popcnt < kmax[i]-1 && andnot(g.cliques[i][readpos], g.explored[i]).zero())
			{
				readpos++;
				continue;
			}
			if(readpos!=writepos)
			{
				g.cliques[writepos]=g.cliques[readpos];
			}
			readpos++;
			writepos++;
		}
		g.cliques[i].resize(writepos);
	}
}

bool construct_incomplete_cliques_simple(graph& g, int target, int kmax[2], bool use_parents, bool min_only);

int depth_search(ext_bigint<2> sig, int kmax[2], int target, int stop_point, bool detect_only, bool quiet, int* populations, graph* graph_pool, int thr, buffer<bigint>* pHits)
{
	if(!sig.zero() && sig.trailing_bit()>=sig.len)
	{
		printf("ERROR: sig has length %d, trailing bit %d\n", sig.len, sig.trailing_bit());
		exit(-1);
	}
	int retval=0;

	target -= 2;
	stop_point -= 2;

	uint64_t tStart = __rdtsc();
	int tIntReportPrev = 0;
	int popPrev = populations[target];

	int ext_target = target+2;
	if(stop_point>GRAPH_MAX_LEN*64+1)
		stop_point=GRAPH_MAX_LEN*64+1;
	int bit_positions[graph_pool_dim];
	int bit_orders[graph_pool_dim];
	
	graph* g, *gmin;
	g = graph_pool;
	gmin = graph_pool+graph_pool_dim;
	bigint2 base_mask;
	base_mask.first=sig;
	base_mask.second=sig;
	base_mask.first.neg();
	base_mask.first.unset_from(sig.len);
//	base_mask.first.printbits();
//	base_mask.second.printbits();
//	exit(-1);
	uint64_t t0=__rdtsc();
	uint64_t op_count=0;
	int len = sig.len;
	bool ok = build_graph(g[len-1], base_mask, kmax);
	if (!ok)
	{
		//omp_set_lock(&op_count_lock);
		//op_counts[2].push_back(op_count);
		//omp_unset_lock(&op_count_lock);
		return 2;
	}
	g[len-1].check_first_unset_bit();
	g[len-1].n=len;
	g[len-1].construct_shifted_masks(target);
	construct_incomplete_cliques_simple(g[len-1], target, kmax, true, true);
	gmin[len-1]=g[len-1];
	
//	build_graph(gmin[len-1], base_mask, kmax);
	gmin[len-1].n=len;
	gmin[len-1].merge();

	int pos = len, bit=0;
	int max_len=0;
	int sym_len=0;
	int sym_known_pos=0;
	int last_bit=0;
	graph& gc=g[len-1];
	gc.n=len;
	gc.check_first_unset_bit();
	ok = max_extensibility(gc, kmax, ext_target, false);
	if(!ok)
		goto done;


	last_bit=max(gc.mask.trailing_bit(), gc.mask0.trailing_bit());
	if(last_bit>len)
		gc.dirty=true;

	for(int i=len; i<graph_pool_dim; i++)
	{
		g[i].reset();
		gmin[i].reset();
	}


	while(true)
	{
		graph& gc = g[pos];
		graph& gp = g[pos-1];

		int delta = pos-len+1;
		int set_count = pos;
		if(delta >= 8)
			set_count -= (delta-2)/6;
		//gc.set_extensions = false;

		int next_bit=set_count;
		bool dup = false;
		if (delta >= 6)
		{
			if (((delta % 6) == 1)
				&& gp.best_next_bit > set_count 
				&& gp.best_next_bit < target - 20
				&& !gp.mask0.bit(gp.best_next_bit) 
				&& !gp.mask.bit(gp.best_next_bit)
				)
			{
				next_bit = gp.best_next_bit;
				dup = true;
			}

			int out_of_order_range = 3; //(g_kmax[0]>=10 || g_kmax[1]>=10) ? 3 : 3;
			if(!dup && (delta%6)<=out_of_order_range)
			{
				int first_unset = 0;
				while ((gp.mask0.bit(first_unset) || gp.mask.bit(first_unset)) && first_unset <= GRAPH_MAX_LEN * 64 && first_unset<=target)
					first_unset++;
				if (first_unset != next_bit && first_unset<=target)
				{
					next_bit = first_unset;
					dup = true;
				}
			}
		}
		int bit_val = bit;
		
		if(g_kmax[0] > g_kmax[1])
			bit_val=1-bit;
		
//		if(delta>=6 && (delta % 6)==1)

		if(gp.mask0.bit(next_bit) || gp.mask.bit(next_bit))
			bit_positions[pos]=-1;
		else
			bit_positions[pos]=next_bit;
		if(bit_positions[pos]>=0)
			bit_orders[pos]=(bit==bit_val);

		bool debug = false;
		bool finish_mode = (kmax[0] >= 10 || kmax[1] >= 10) && (popcnt(gp.mask) + popcnt(gp.mask0) >= target - 30);
		bool detect_only_mode = (kmax[bit_val] >= 10 && finish_mode);

		if(gp.mask0.bit(next_bit) && bit_val==1)
			goto next;
		if(gp.mask.bit(next_bit) && bit_val==0)
			goto next;

		if(g_search_trace)
		{
			bigint testval;
			testval.n[0] = 0x61a39b78f6eede ;
			testval.n[1] = 0x3070ee85dc3830 ;
			testval.n[2] = 0x1edddbc7b671618 ;
			testval.n[3] = 0;
			bigint dual_testval = testval;
			dual_testval.neg();
			dual_testval.unset_from(185);
			debug = bit_equal(gp.mask, gp.mask & testval) && bit_equal(gp.mask0, gp.mask0 & dual_testval) && (testval.bit(next_bit) == bit_val);
			if (debug)
			{
				printf("%d ", op_count);
				for(int i=sig.len; i<target; i++)
				{
					if(gp.mask0.bit(i))
						printf("0");
					else if(gp.mask.bit(i))
						printf("1");
					else
						printf("_");
				}
				printf(" %d<=%d", next_bit, bit_val);
				printf(" (%d)\n", gp.clique_counts[0][kmax[0]-1] + gp.clique_counts[1][kmax[1]-1]);
			}
		}

		{
			uint64_t t = __rdtsc();
			int tInt = (int)floor(((t-tStart)/g_fFreq)/900.);
			if(t>=tStart && tInt>tIntReportPrev && set_count>=len+10)
			{
				double progress = 0.0;
				for(int i=len; i<len+8; i++)
					if(gp.mask.bit(i))
						progress += pow(2.0, -(i-len)-1);
				printf("...thread %d sig 0x%llx_0x%llx: %.1f%% done, max %d, %I64d ops\n", thr, 
					sig.n[0], sig.n[1], progress*100.0, max_len+1, op_count);
				tIntReportPrev = tInt;
			}
		}
		op_count++;
		if(g_max_ops!=0 && op_count > g_max_ops)
		{
			if(g_search_trace)
			printf("=> timeout\n");
			if(detect_only)
				return 0;
			else if(g_max_ops_2==0)
				return (populations[target]!=popPrev) ? 1 : 0;
			else
			{
				if(pHits->size()==0)
					return 0;
				if(op_count>g_max_ops_2)
					return 1;
			}
		}
		max_len = max(max_len, set_count);

		if(set_count<stop_point && extend(gc, gp, bit_val, kmax, ext_target, next_bit, detect_only_mode ? 3 : 2))
		{
			int first_unset=0;
			while((gc.mask0.bit(first_unset) || gc.mask.bit(first_unset)) && first_unset<=GRAPH_MAX_LEN*64)
				first_unset++;
			int last_bit=max(gc.mask.trailing_bit(), gc.mask0.trailing_bit());
			gc.dirty = gp.dirty || (last_bit != first_unset-1);
			gc.n = gp.n+1;

			if(gc.clique_counts[0][kmax[0]]>0  || gc.clique_counts[1][kmax[1]]>0)
				goto next;

			if(delta >= 6 && !(delta % 6) && set_count<target)
			{
				if(first_unset <= set_count)
				{
					printf("ERROR: assume first_unset>set_count (%d, %d)\n", first_unset, set_count);
					exit(-1);
				}
				if(!gc.dirty)
				{
					gmin[pos] = gc;
					gmin[pos].merge();
				}
				else
				{
					if (g_search_trace && debug)
						printf("dirty rebuild\n");
					bigint mask0 = gc.mask0;
					bigint mask = gc.mask;
					graph& gprev = (g[pos - 6].dirty) ? gmin[pos - 6] : g[pos - 6];
					if (gprev.n == 0)
					{
						printf("Error: null gprev\n");
						exit(-1);
					}
					bigint prev_unset_mask = gprev.mask|gprev.mask0;
					prev_unset_mask.neg();
					int prev_first_unset=prev_unset_mask.leading_bit();
					if(mask.bit(prev_first_unset) || mask0.bit(prev_first_unset))
					{
						gmin[pos].reset();
#if 1
						if(!extend(gmin[pos], gprev, mask.bit(prev_first_unset), kmax, ext_target, prev_first_unset, 0))
//						if(!extend(gmin[pos], gprev, mask.bit(set_count-5), kmax, ext_target, set_count-5, 0))
						{
							if (g_search_trace && debug)
								printf("Fail 1\n");
							goto next;
						}

						for (int i = 0; i<=target; i++)
						{
							if((gprev.mask.bit(i) && !mask.bit(i)) 
								|| (gprev.mask0.bit(i) && !mask0.bit(i)))
							{
								printf("ERROR: unexpected gprev state\n");
								gprev.mask0.print();
								mask0.print();
								gprev.mask.print();
								mask.print();
								exit(-1);
							}								
						}

						for (int i = prev_first_unset+1; i<=target; i++)
						{
							if (!(mask0.bit(i) || mask.bit(i)))
								break;
							if (!extend(gmin[pos], mask.bit(i), kmax, ext_target, i, !finish_mode))
							{
								if(debug)
								{
									printf("extend(gmin[pos]), set_count=%d, i=%d: fail\n", set_count, i);
								}
								goto next;
							}
						}
#else
						if(!extend(gmin[pos], gprev, mask.bit(set_count-5), kmax, ext_target, set_count-5, 0))
							goto next;
						for (int i = 1; ; i++)
						{
							int b = set_count - 5 + i;
							if (b > target)
								break;
							if (!(mask0.bit(b) || mask.bit(b)))
								break;
							if (!extend(gmin[pos], mask.bit(b), kmax, ext_target, b, !finish_mode))
								goto next;
						}
#endif
					}
					else
					{
						gmin[pos]=gprev;
					}


					int last_known_bit = max(gmin[pos].mask.trailing_bit(), gmin[pos].mask0.trailing_bit());
					if(last_known_bit+1 != popcnt(gmin[pos].mask)+popcnt(gmin[pos].mask0))
					{
						printf("ERROR: invalid gmin\n");
						printf("pos=%d\n", pos);
						printf("delta=%d\n", delta);
						printf("last_known_bit=%d\n", last_known_bit);
						printf("popcnt %d\n", popcnt(gmin[pos].mask)+popcnt(gmin[pos].mask0));
						gmin[pos].mask0.print();
						gmin[pos].mask.print();
						exit(-1);
					}


					if(gmin[pos].parent_cliques[0].size()>=8)
						gmin[pos].merge();
					strip(gmin[pos], kmax);

					gmin[pos].check_first_unset_bit();
					if(target >= gmin[pos].first_unset_bit + 10)
					{
						gmin[pos].construct_shifted_masks(target);
						construct_incomplete_cliques_simple(gmin[pos], target, kmax, true, true);
					}

					int last_bit = max(gc.mask.trailing_bit(), gc.mask0.trailing_bit());
					if(last_bit == last_known_bit)
					{
						gc = gmin[pos];
						gc.dirty = false;
					}
					else
					{
						if(!extend_from(gc, gmin[pos], make_pair(gc.mask0,gc.mask), kmax, ext_target))
							goto next;
//						if(last_bit+1 == popcnt(gc.mask0)+popcnt(gc.mask))
//							gc.dirty = false;
//						else
							gc.dirty = true;
					}

					if (g_search_trace && debug)
						printf("dirty rebuild done (dirty=%d)\n", (int)gc.dirty);
				}
				//else
				if(gc.clique_counts[0][kmax[0]]>0  || gc.clique_counts[1][kmax[1]]>0)
				{
					if (g_search_trace && debug)
						printf("%d %d\n", gc.clique_counts[0][kmax[0]], gc.clique_counts[1][kmax[1]]);
					goto next;
				}

				int min_unset_bit = 35;
				if(kmax[0]<5 || kmax[1]<5)
					min_unset_bit=10;

				if((set_count > min_unset_bit) && (set_count+20 < target) && !finish_mode)
				{
					gc.n=pos;
					gc.check_first_unset_bit();
					bool ok = max_extensibility(gc, kmax, ext_target, debug);
					if(!ok)
					{
						goto next;
					}
					else
					{
//						int last_bit=max(gc.mask.trailing_bit(), gc.mask0.trailing_bit());
//						if(last_bit>pos)
//							gc.dirty=true;
						int last_known_bit = max(gc.mask.trailing_bit(), gc.mask0.trailing_bit());
						if(last_known_bit+1 != popcnt(gc.mask)+popcnt(gc.mask0))
							gc.dirty=true;
					}
				}
			}
			
			if(set_count == target)
			{
				int first_unset=0;
				while((gc.mask0.bit(first_unset) || gc.mask.bit(first_unset)) && first_unset<=GRAPH_MAX_LEN*64)
					first_unset++;
				if(first_unset<target)
				{
					printf("thread %d, sig 0x%llx: ERROR: set_count %d, but bit %d is not set\n", thr, sig.n[0], set_count, first_unset);
					exit(-1);
				}
				if(!validate_solution(gc, kmax, target, gmin, pos))
				{
//					if(g_search_trace)
//						printf(" - validate_solution fail\n");
					goto next;
				}
				else
				{
//					if(g_search_trace)
//						printf(" - hit\n");
				}
				if(detect_only)
				{
					populations[set_count]++;

					if(pHits!=0)
					{
						bigint x = gc.mask;
						x.set_len(set_count+1);
						pHits->push_back(x);
					}

					retval=1;
					goto done;
				}
			}
			if(!dup)
			{
				if(set_count>=target && pHits!=0)
				{
					if(!validate_solution(gc, kmax, target, gmin, pos))
					{
//						if(g_search_trace)
//							printf(" - validate_solution fail\n");
						goto next;
					}
					else
					{
//						if(g_search_trace)
//							printf(" - hit\n");
					}
					populations[set_count]++;

					int first_unset=0;
					while((gc.mask0.bit(first_unset) || gc.mask.bit(first_unset)) && first_unset<=GRAPH_MAX_LEN*64)
						first_unset++;
					if(first_unset<set_count)
					{
						printf("ERROR: set_count %d, but bit %d is not set\n", set_count, first_unset);
						graph gg;
						bigint2 m = make_pair(gc.mask0, gc.mask);
						for(int i=0; i<4; i++)
							printf("%llx %llx %llx\n", gc.mask0.n[i], gc.mask.n[i], gc.mask0.n[i] | gc.mask.n[i]);
						bool ok = build_graph(gg, m, kmax);
						gg.print_short();
						exit(-1);
					}

					bigint x = gc.mask;
					x.set_len(set_count+1);
					pHits->push_back(x);
					retval=1;
					
					bigint inv;
					invert(inv, x);
					if(bit_equal(x, inv))
					{
//						printf("%d\t0x%llx\t0x%llx\t0x%llx\t0x%llx\n", 
//							x.len+1, x.n[0], x.n[1], x.n[2], x.n[3]);
					}
					else
					{
						//printf("%d\t0x%llx\t0x%llx\t0x%llx\t0x%llx\n", 
						//	x.len+1, x.n[0], x.n[1], x.n[2], x.n[3]);
					//	x.print();
					}
					
					if(detect_only)
					{
						goto done;
					}
					else if(set_count>=stop_point)
					{
						goto next;
					}
				}
				else
				{
					populations[set_count]++;
				}

				gc.check_first_unset_bit();
			
				if(detect_only && set_count < target && gc.first_unset_bit >= target+2)
				{
					if(!validate_solution(gc, kmax, target, gmin, pos))
						goto next;
					for(int i=set_count+1; i<=target; i++)
						populations[i]++;
					if(pHits!=0)
					{
						bigint x = gc.mask;
						x.set_len(gc.first_unset_bit);
						pHits->push_back(x);
					}
					retval=1;
					goto done;
				}
				if((!detect_only) && set_count < stop_point && gc.first_unset_bit >= stop_point+2)
				{
					if(!validate_solution(gc, kmax, stop_point, gmin, pos))
						goto next;
					for(int i=set_count+1; i<=stop_point; i++)
						populations[i]++;
					if(pHits!=0)
					{
						bigint x = gc.mask;
						x.set_len(gc.first_unset_bit);
						pHits->push_back(x);
					}
					retval=1;
					goto next;
				}
			}
			pos++;
//			if((delta % 6)!=1)
//				set_count++;
			bit=0;
			//printf("=> level down\n");
			continue;
		}

next:
		if(bit==0)
		{
			bit=1;
		}
		else
		{
			pos--;
			while(true)
			{
				if(pos>=len && (bit_positions[pos]<0 || 
					(bit_orders[pos] 
						? g[pos].mask.bit(bit_positions[pos])
						: g[pos].mask0.bit(bit_positions[pos]))
						))
					pos--;
				else
					break;
			}
			if(pos<len)
				break;
			bit=1;
		}
	}
	if(!quiet)
	{
		uint64_t t1=__rdtsc();
		printf("<%.3f>\n", (t1-t0)/g_fFreq);
		for(int i=30; i<pop_acc_dim; i++)
		{
			printf("%d\t%d\n", i+2, populations[i]);
		}
	}

done:
//	if(g_search_trace)
//		printf("=> %d (%d, %d)\n", retval, (*pHits).size(), populations[target]-popPrev);

	if(retval==2)
	{
		sig.print();
		printf("retval==2?\n");
		return 2;
	}
	if((retval==1) || (populations[target]!=popPrev))
	{
		//printf("Sig 0x%llx 0x%llx success\n", sig.n[0], sig.n[1]);

		omp_set_lock(&op_count_lock);
		op_counts[1].push_back(op_count);
		omp_unset_lock(&op_count_lock);
		return 1;
	}
	else
	{
//		printf("Sig 0x%llx 0x%llx fail\n", sig.n[0], sig.n[1]);
		omp_set_lock(&op_count_lock);
		op_counts[0].push_back(op_count);
		omp_unset_lock(&op_count_lock);
		return 0;
	}
}

int depth_search(uint64_t sig, int len, int kmax[2], int target, int stop_point, bool detect_only, bool quiet, int* populations, graph* graph_pool, int thr, buffer<bigint>* pHits=0)
{
	ext_bigint<2> b;	
	b.clear();
	b.n[0]=sig;
	b.set_len(len);
	return depth_search(b, kmax, target, stop_point, detect_only, quiet, populations, graph_pool, thr, pHits);
}




bool quiet=false;

int undivided_depth_search_job(void* arg)
{
	depth_search_state* pState =(depth_search_state*)arg;
	uint64_t t0 = __rdtsc();
	uint64_t pos = pState->thr;
	int thr = pState->thr;
	buffer<bigint> hits;
	int local_populations[pop_acc_dim];
	buffer<b128> cover_list;
	cover_list.alloc(1.5*pState->nbases/pState->nthr);
	while(true)
	{
		uint64_t share = (pos+pState->nbases/100 - 1) * 100 / pState->nbases;
		if (pos == pState->nbases*share / 100 && !quiet)
			printf("%d%%...\n", share);

		b128 val = pState->bases[pos];
		hits.clear();
		memset(local_populations, 0, sizeof(local_populations));
		int ok = depth_search(val, g_kmax, pState->target, pState->stop_size, pState->detect_only, true, local_populations, pState->graph_pool+thr*graph_pool_dim*2, thr, &hits);
		if(ok!=2)
			for(int i=0; i<pop_acc_dim; i++)
				pState->populations[thr*pop_acc_dim+i]+=local_populations[i];
		int max_len=pState->len+1;
		while(max_len<pop_acc_dim && local_populations[max_len]!=0)
			max_len++;
		max_len--;
		max_len+=2;
		uint64_t prevpos = pos;
		unsigned int* x = pState->pNextPos;
#ifdef KNIGHTS_CORNER
#pragma omp atomic capture
		{pos = *x; (*x)++;}
#elif !defined(WIN32)
		pos=__sync_fetch_and_add(x, 1);//, __ATOMIC_SEQ_CST);
#else
		pos=InterlockedIncrement(pState->pNextPos)-1;
#endif
		if((ok==1 && pState->pHitDetect)
			 || (ok!=1 && pState->pHitMiss)
			 || (pState->len<49 && pState->pBaseCover))
		{
			omp_set_lock(pState->pCritSect);
			if(ok==1)
			{
				if(pState->pHitDetect)
					pState->pHitDetect(prevpos, val, &hits, max_len);
			}
			else
			{
				if(pState->pHitMiss)
					pState->pHitMiss(prevpos, val, &hits, max_len);
			}

			if(pState->len<49 && pState->pBaseCover)
				pState->pBaseCover(val, 0);

			omp_unset_lock(pState->pCritSect);
		}
		if(pState->len>=49)
			cover_list.push_back(val);

		if(pos >= pState->nbases)
			break;
	}
	if(pState->len>=49)
	{
		if(pState->pBaseCover)
		{
			omp_set_lock(pState->pCritSect);
			for(int i=0; i<cover_list.size(); i++)
				pState->pBaseCover(cover_list[i], 0);
			omp_unset_lock(pState->pCritSect);
		}
	}
	uint64_t t1 = __rdtsc();
	pState->cputime = t1-t0;
	return 0;
}

void JobArray::print_populations()
{
	int stop_pos = pop_acc_dim-1;
	while(sum_populations[stop_pos]==0 && stop_pos>0)
		stop_pos--;
	for(int i=2; i<=stop_pos; i++)
		printf("%d %d\n", i, sum_populations[i]);
}

void JobArray::print_cputime()
{
	printf("<%.3f><%.3f>\n", realtime/g_fFreq, cputime/g_fFreq);
}

void JobArray::submit(const ext_bigint<2>* p, size_t n, int len, int target, int stop_size, bool detect_only, int fractions, base_cover_callback_fun* pBaseCover, hit_detect_callback_fun* pHitDetect, hit_detect_callback_fun* pHitFail)
{
	memset(populations, 0, nthreads*pop_acc_dim*sizeof(int));
	memset(sum_populations, 0, sizeof(sum_populations));
	prealloc_graph_pool(graph_pool, len);
	uint64_t* base_masks=new uint64_t[n];
	for(int i=0; i<n; i++)
		base_masks[i]=0;
	unsigned int next_pos=nthreads;
	uint64_t t0 = __rdtsc();
	for(int thr=0; thr<nthreads; thr++)
	{
		jobs[thr].cputime = 0;
		if(thr>=n)
			continue;
		jobs[thr].bases = p;
		jobs[thr].nbases = n;
		jobs[thr].base_masks = base_masks;
		jobs[thr].detect_only = detect_only;
		jobs[thr].graph_pool = graph_pool;
		jobs[thr].len = len;
		jobs[thr].fractions = 0;
		jobs[thr].pCritSect = &critSect;
		jobs[thr].pNextPos = &next_pos;
		jobs[thr].populations = populations;
		jobs[thr].target = target;
		jobs[thr].stop_size = stop_size;
		jobs[thr].thr = thr;
		jobs[thr].nthr = nthreads;
		jobs[thr].t0 = t0;
		jobs[thr].pBaseCover=pBaseCover;
		jobs[thr].pHitDetect=pHitDetect;
		jobs[thr].pHitMiss=pHitFail;
	}

	int nthr = min(nthreads, (int)n);

#ifdef PARALLEL
#pragma omp parallel for
	for(int thr=0; thr<nthr; thr++)
		undivided_depth_search_job((void*)(jobs+thr));

#else
	for(int thr=0; thr<nthr; thr++)
		undivided_depth_search_job((void*)(jobs+thr));
#endif

	for(int thr=0; thr<nthreads; thr++)
		for(int i=0; i<pop_acc_dim-2; i++)
			sum_populations[i+2]+=populations[thr*pop_acc_dim+i];

	cputime=0;
	for(int thr=0; thr<nthreads; thr++)
		cputime+=jobs[thr].cputime;

	uint64_t t1 = __rdtsc();
	realtime = t1-t0;
//	printf("<%.3f>\n", (t1-t0)/g_fFreq);
	delete[] base_masks;
}



vector<ext_bigint<2> > mc_memory_hits;
void mc_memory_hit_detect(uint64_t pos, ext_bigint<2> sig, buffer<bigint>* p, int max_len)
{
	mc_memory_hits.push_back(sig);
}



string g_hit_detect_name;
string g_base_cover_name;

string format_name(int len, int target, string prefix)
{
	char temp[128];
	sprintf(temp, PATH "%s%d%d-%d-%d.txt", prefix.c_str(), g_kmax[0]+2, g_kmax[1]+2, len+1, target);
	return string(temp);
}

string format_symcover_name(int len, int target)
{
	char temp[128];
#if GRAPH_MAX_LEN==6
	sprintf(temp, PATH "symcover6_%d%d-%d.txt", g_kmax[0]+2, g_kmax[1]+2, len+1);
#else
	sprintf(temp, PATH "symcover%d%d-%d.txt", g_kmax[0]+2, g_kmax[1]+2, len+1);
#endif
	return string(temp);
}

string format_sym_name(string prefix)
{
	char temp[128];
	sprintf(temp, PATH "sym%s%d%d.txt", prefix.c_str(), g_kmax[0]+2, g_kmax[1]+2);
	return string(temp);
}

string format_cover_name(int len, int target)
{
	char temp[128];
	sprintf(temp, PATH "cover%d%d-%d-%d.txt", g_kmax[0]+2, g_kmax[1]+2, len+1, target);
	return string(temp);
}

string format_cover_name_offby2(int len, int target)
{
	char temp[128];
	sprintf(temp, PATH "cover%d%d-%d-%d-offby2.txt", g_kmax[0]+2, g_kmax[1]+2, len+1, target);
	return string(temp);
}

string format_review_name(int len, int target)
{
	char temp[128];
	sprintf(temp, PATH "review%d%d-%d-%d.txt", g_kmax[0]+2, g_kmax[1]+2, len+1, target);
	return string(temp);
}


void readi64(FILE* f, uint8_t* p, size_t sz)
{
	while(sz >= 1000000000)
	{
		fread(p, sz, 1, f);
		sz-=1000000000;
		p+=1000000000;
	}
	if(sz!=0)
		fread(p, sz, 1, f);
}

template <int T>
void readfile(const char* fn, vector<ext_bigint<T> >& sigs, int len)
{
	string s = string(fn)+".bin";
	FILE* f;

	if(!strstr(fn,"sym"))
	{
		f=fopen(s.c_str(), "rb");
		if(f)
		{
			if(len>128)
			{
				printf("ERROR: readfile() can't work with len>128\n");
				exit(-1);
			}
			_fseeki64(f, 0, SEEK_END);
			uint64_t pos = _ftelli64(f);
			vector<uint64_t> temp;
			temp.resize(pos/8);
			rewind(f);
			readi64(f, (uint8_t*)&temp[0], pos);
			fclose(f);
			size_t sizePrev = sigs.size();
			sigs.resize(sizePrev+temp.size()/2);
			for(int i=0; i<temp.size(); i+=2)
			{
				ext_bigint<T>& x = sigs[sizePrev+i/2];
				x.clear();
				x.set_len(len);
				x.n[0]=temp[i];
				x.n[1]=temp[i+1];
			}

			return;
		}
	}
	f=fopen(fn, "r");
	if(!f)
		return;
	
	size_t sizePrev = sigs.size();
	while(!feof(f))
	{
		char line[1024];
		line[0]=0;
		fscanf(f, "%[^\n]s", line);
		fscanf(f, "\n");
		uint64_t v1=0, v2=0;
		if(len<=64)
		{
			sscanf(line, "0x%llx", &v1);
			if(v1!=0)
			{
				ext_bigint<T> b;
				b.clear();
				b.set_len(len);
				b.n[0]=v1;
				sigs.push_back(b);
			}
		}
		else
		{
			sscanf(line, "0x%llx\t0x%llx", &v1, &v2);
			if(v1!=0)
			{
				ext_bigint<T> b;
				b.clear();
				b.set_len(len);
				b.n[0]=v1;
				b.n[1]=v2;
				sigs.push_back(b);
			}
		}
	}
	printf("%s: %d graphs\n", fn, sigs.size()-sizePrev);
	fclose(f);
}

template<int T>
void writefile(const char* fn, vector<ext_bigint<T> >& sigs, int len, bool append)
{
	if(append==false && sigs.empty())
	{
		_unlink(fn);
		return;
	}

	if(!strstr(fn,"sym"))
	{
		string s = string(fn)+".bin";
		if(append)
		{
			FILE* f=fopen(s.c_str(), "rb");
			
			if(f==0)
			{
				vector<ext_bigint<T> > x;
				readfile(fn,x,len);
				f=fopen(s.c_str(), "wb");
//				fwrite(&x[0], 1, x.size()*sizeof(bigint), f);
				vector<uint64_t> temp;
				temp.resize(x.size()*2);
				for(uint64_t i=0; i<x.size(); i++)
				{
					temp[2*i]=x[i].n[0];
					temp[2*i+1]=x[i].n[1];
				}
				uint64_t sz = temp.size()*8;
				uint8_t* p = (uint8_t*)(&temp[0]);
				while(sz > 1000000000)
				{
					fwrite(p, 1000000000, 1, f);
					p+=1000000000;
					sz-=1000000000;
				}
				if(sz!=0)
					fwrite(p, sz, 1, f);
				fclose(f);
			}
			else
			{
				fclose(f);
			}
		}

		FILE* f=fopen(s.c_str(), append ? "ab" : "wb");
//		fwrite(&sigs[0], 1, sigs.size()*sizeof(bigint), f);
		vector<uint64_t> temp;
		temp.resize(sigs.size()*2);
		for(uint64_t i=0; i<sigs.size(); i++)
		{
			temp[2*i]=sigs[i].n[0];
			temp[2*i+1]=sigs[i].n[1];
		}
		uint64_t sz = temp.size()*8;
		uint8_t* p = (uint8_t*)(&temp[0]);
		while(sz > 1000000000)
		{
			fwrite(p, 1000000000, 1, f);
			p+=1000000000;
			sz-=1000000000;
		}
		if(sz!=0)
			fwrite(p, sz, 1, f);
		fclose(f);
		return;
	}

	FILE* f=fopen(fn, append ? "a" : "w");
	if(!f)
	{
		printf("ERROR: failed to open %s for writing\n", fn);
		exit(-1);
		return;
	}
	
	if(len<=64)
	{
		for(int i=0; i<sigs.size(); i++)
			fprintf(f, "0x%llx\n", sigs[i].n[0]);
	}
	else
	{
		for(int i=0; i<sigs.size(); i++)
			fprintf(f, "0x%llx\t0x%llx\n", sigs[i].n[0], sigs[i].n[1]);
	}
	fclose(f);
}

template <int T>
void writefile(string str, vector<ext_bigint<T> >& sigs, int len, bool append)
{
	writefile(str.c_str(), sigs, len, append);
}

template <int T>
void readfile(string fn, vector<ext_bigint<T> >& sigs, int len)
{
	readfile(fn.c_str(), sigs, len);
}


template <class T>
void writesym(T it0, T it1, const char* filename, int min_len)
{
	FILE* f = 0;
	while(it0!=it1)
	{
		bigint x = *it0;
		if(x.len+1<min_len)
		{
			it0++;
			continue;
		}
		if(f==0)
			 f=fopen(filename, "a");
#if GRAPH_MAX_LEN>=6
		fprintf(f, "%d\t0x%llx\t0x%llx\t0x%llx\t0x%llx\t0x%llx\t0x%llx\n", x.len + 1, x.n[0], x.n[1], x.n[2], x.n[3], x.n[4], x.n[5]);
#elif GRAPH_MAX_LEN>=4
		fprintf(f, "%d\t0x%llx\t0x%llx\t0x%llx\t0x%llx\t0x%llx\t0x%llx\n", x.len + 1, x.n[0], x.n[1], x.n[2], x.n[3], 0, 0);
#else
		fprintf(f, "%d\t0x%llx\t0x%llx\t0x%llx\t0x%llx\t0x%llx\t0x%llx\n", x.len + 1, x.n[0], x.n[1], 0, 0, 0, 0);
#endif
		it0++;
	}
	if(f!=0)
		fclose(f);
}

template <int T>
size_t exclude(vector<ext_bigint<T> >& vsig, vector<ext_bigint<T> >& vsigExcl, int len=0)
{
	//printf("exclude(%d, %d, %d)\n", vsig.size(), vsigExcl.size(), len);
	if(vsigExcl.size()==0)
		return 0;
	size_t i;
	if(len!=0)
	{
		for(i=0; i<vsigExcl.size(); i++)
			vsigExcl[i].unset_from(len);
	}
	rem_dup(vsigExcl);

	int64_t sizeInit = vsig.size();

	if(sizeInit < 10000)
	{
		int64_t sizePrev = sizeInit;
		for(int64_t i=0; i<sizeInit; i++)
		{
			ext_bigint<T> q = vsig[i];
			if(len!=0)
				q.unset_from(len);
			bool dup=binary_search(vsigExcl.begin(), vsigExcl.end(), q);
			if(dup)
			{
				vsig[i] = vsig.back();
				vsig.pop_back();
				i--;
				sizeInit--;
			}
		}

		return sizePrev-vsig.size();
	}

	int64_t endpos[max_nthreads];

#pragma omp parallel for
	for(int64_t thr=0; thr<nthreads; thr++)
	{
		int64_t low=sizeInit*thr/nthreads;
		int64_t high=sizeInit*(thr+1)/nthreads;
		
		for(int64_t i=low; i<high; i++)
		{
			ext_bigint<T> q = vsig[i];
			if(len!=0)
				q.unset_from(len);
			bool dup=binary_search(vsigExcl.begin(), vsigExcl.end(), q);
			if(dup)
			{
				vsig[i] = vsig[high-1];
				high--;
				i--;
			}
		}
		endpos[thr]=high;
	}
//	for(int64_t thr=0; thr<16; thr++)
//		printf("%d ", vsig2[thr].size());
//	printf("\n");

	vsig.clear();
	size_t total_size=endpos[0];
	for(int64_t thr=1; thr<nthreads; thr++)
	{
		int64_t low=sizeInit*thr/nthreads;
		int64_t high = endpos[thr];
		//size_t sizePrev = total_size;
		//total_size+=vsig2[thr].size();
		//vsig.resize(total_size);
		memmove(&vsig[total_size], &vsig[low], (high-low)*sizeof(vsig[0]));
		total_size += high-low;
	}
	vsig.resize(total_size);
	return sizeInit-vsig.size();
}

bool mc_detect_only=true;

ext_bigint<2> random_bigint(int len)
{
	ext_bigint<2> x;
	x.clear();
	x.set_len(len);
	if(len<=64)
	{
		x.n[0]=rand64() & ((One<<len)-2);
	}
	else
	{
		x.n[0]=rand64() & (-2ll);
		x.n[1]=rand64() & ((One<<(len-64))-2);
	}
	return x;
}

template <int T>
void clusterize(vector<vector<ext_bigint<T> > >& vv, const vector<ext_bigint<T> >& b, bool sym=false);

bool ns_detect_only=true;

template <int T>
bool is_incomplete(ext_bigint<T>  x, const vector<ext_bigint<T> >& cover, int d)
{
	for (int j = 0; j < d - 1; j++)
	{
		ext_bigint<T>  y = x;
		y.flip(j + 1);
		bool dup = binary_search(cover.begin(), cover.end(), y);
		if (!dup)
			return true;
	}
	return false;
}


void narrow_hit_detect(uint64_t pos, b128 sig, buffer<bigint>* p, int size)
{
	FILE* f=fopen((g_hit_detect_name+".bin").c_str(), "ab");
	fwrite(&sig.n[0],16,1,f);
	fclose(f);
}

void finish_cyclic_search(int k0, int k1, const vector<b128>& hits, vector<bigint>& sym, int min_save);

template <int T>
int dist(ext_bigint<T>  a, ext_bigint<T>  b)
{
	return __popcnt64(a.n[0]^b.n[0]) + __popcnt64(a.n[1]^b.n[1]);
}

template <int T>
int dist(ext_bigint<T>  a, const vector<ext_bigint<T>  >& vb)
{
	int min_dist=100;
	for(int i=0; i<vb.size(); i++)
	{
		int d=dist(a,vb[i]);
		if(min_dist>d)
			min_dist=d;
	}
	return min_dist;
}

template <int T>
int dist2(ext_bigint<T>  a, const vector<ext_bigint<T> >& vb)
{
	int min_dist=100;
	for(int i=0; i<vb.size(); i++)
	{
		min_dist=min(min_dist, dist(a,vb[i]));
	}
	return min_dist;
}

template <int T>
int dist(const vector<ext_bigint<T> >& va, const vector<ext_bigint<T> >& vb)
{
	int min_dist=100;
	for(int i=0; i<va.size(); i++)
		min_dist=min(min_dist, dist(va[i],vb));
	return min_dist;
}

template <int T>
void clusterize(vector<vector<ext_bigint<T> > >& vv, const vector<ext_bigint<T> >& v, bool sym)
{
	uint64_t t1, t2;
	t1=__rdtsc();
	if(v.size()==0)
	{
		vv.clear();
		return;
	}

	vector<ext_bigint<T> > v2 = v;
	rem_dup(v2);

	vv.clear();
	while(v2.size()>0)
	{
		vector<ext_bigint<T> > cl;
		cl.push_back(v2.back());
		v2.pop_back();
		int pos1=0, pos2=1;
		while(true)
		{
			vector<int> positions;
			for(int i=pos1; i<pos2; i++)
			{
				int d = cl[i].len;
				int range = sym ? (d + 1) / 2 : d;
				if(sym)
				{
					for(int j=0; j<range; j++)
					{
						for(int k=j; k<range; k++)
						{
							ext_bigint<T>  x = cl[i];
							x.flip(j);
							if (sym && d - j - 1 != j)
								x.flip(d - j - 1);
							if(k!=j)
							{
								x.flip(k);
								if (sym && d - k - 1 != k)
									x.flip(d - k - 1);
							}
							int pos = lower_bound(v2.begin(), v2.end(), x) - v2.begin();
							if(pos>=0 && pos<v2.size() && bit_equal(v2[pos], x) && v2[pos].len>0)
							{
								positions.push_back(pos);
		//						v2[pos].len=0;
							}
						}
					}
				}
				else
				{
					for(int j=0; j<range; j++)
					{
						ext_bigint<T>  x = cl[i];
						x.flip(j);
						if (sym && d - j - 1 != j)
							x.flip(d - j - 1);
						int pos = lower_bound(v2.begin(), v2.end(), x) - v2.begin();
						if(pos>=0 && pos<v2.size() && bit_equal(v2[pos], x) && v2[pos].len>0)
						{
							positions.push_back(pos);
	//						v2[pos].len=0;
						}
					}
				}
			}
			rem_dup(positions);
			for(int i=0; i<positions.size(); i++)
			{
				cl.push_back(v2[positions[i]]);
				v2[positions[i]].len=0;
			}
			pos1=pos2;
			pos2=cl.size(); 
			if(pos1==pos2)
				break;
		}
		vv.push_back(cl);
		vector<ext_bigint<T> > v3;
		v3.resize(v2.size());
		int pos=0;
		for(int i=0; i<v2.size(); i++)
		{
			if(v2[i].len>0)
				v3[pos++]=v2[i];
		}
		v3.resize(pos);
		v2.swap(v3);

	}
	t2=__rdtsc();
	printf("%d clusters in %.3f s\n", vv.size(), (t2-t1)/g_fFreq);
}

void cover_rebuild(int len, int k0, int k1, int target, const char* secondary=0)
{
	g_kmax[0]=k0;
	g_kmax[1]=k1;
	string s=format_name(len, target);
	string sCover=format_cover_name(len, target);
	vector<b128> v, vsig;
	readfile(s.c_str(), v, len);

	vsig.resize((len-1)*v.size());
	for(int i=0; i<v.size(); i++)
		for(int j=0; j<len-1; j++)
		{
			vsig[i*(len-1)+j]=v[i];
			if(v[i].bit(j+1))
				vsig[i*(len-1)+j].unset(j+1);
			else
				vsig[i*(len-1)+j].set(j+1);
		}
	append(vsig, v);
	if(secondary!=0)
	{
		v.clear();
		readfile(secondary, v, len);
		append(vsig, v);
	}
	rem_dup(vsig);
	writefile(sCover.c_str(), vsig, len);
}

string fs_hit_name;
void fs_hit_detect(uint64_t pos, b128 sig, buffer<bigint>* p, int max_len)
{
	//FILE* f = fopen(PATH "413-21-128.txt", "a");
	FILE* f = fopen(fs_hit_name.c_str(), "a");
	if(sig.len>64)
		fprintf(f, "%d\t%d\t0x%llx\t0x%llx\n", g_kmax[0], g_kmax[1], sig.n[0], sig.n[1]);
	else
	//	fprintf(f, "%d\t%d\t0x%llx\n", g_kmax[0], g_kmax[1], sig.n[0]);
		fprintf(f, "0x%llx\n", sig.n[0]);
	fclose(f);
	for (int i = 0; i<p->size(); i++)
	{
		bigint& x = (*p)[i];
		int d = x.len + 1;
		if (d == max_len)
		{
			bigint inv;
			invert(inv, x);
			if (bit_equal(x, inv))
			{
				printf("%d\t0x%llx\t0x%llx\t0x%llx\t0x%llx\n", d, x.n[0], x.n[1], x.n[2], x.n[3]);
			}
		}
	}
}

void do_full_search(int k0, int k1, int dim, int target)
{
	g_max_ops=0;
	g_kmax[0]=k0;
	g_kmax[1]=k1;
	fs_hit_name = format_name(dim, target);
	vector<b128> v;
	size_t nOpts = 1<<dim;
	int covered = 0;//1024*3;
	v.resize(nOpts-covered);
	for(size_t i=covered; i<nOpts; i++)
	{
		v[i-covered].clear();
		v[i-covered].set_len(dim);
		v[i-covered].n[0]=i;
	}
	vector<b128> vExcl;
	printf("%d left\n", v.size());
	random_shuffle(v.begin(), v.end());
	JobArray array;
	array.submit(&v[0], v.size(), dim, target, GRAPH_MAX_LEN*64+1, false, 0, 0, fs_hit_detect);
	printf("%d -> %d\n", v.size(), array.sum_populations[target]);
	array.print_populations();
	array.print_cputime();
}


string g_symfn, g_symcoverfn;
vector<bigint> g_sym_graphs, g_sym_graphs_ref;

bool do_early_enum=true;
int max_early_enum=7;
int max_errors_detect_only = 10;
int max_errors=6;

vector<bigint> sym_potential_graphs;
bool test_build(graph& g, bigint b, int k0, int k1);
int g_sym_target = 0;
void sym_hit_detect(uint64_t pos, b128 sig, buffer<bigint>* p, int max_len)
{
	int i;
	for (i = 0; i<p->size(); i++)
	{
		bigint& xx = (*p)[i];
//		int d = x.len + 1;
		for(int d=xx.len+1; d>=g_sym_target; d--)
		{
			bigint x = xx;
			x.set_len(d-1);
			x.unset_from(d-1);
			bigint inv;
			invert(inv, x);
			if (bit_equal(x, inv) && !binary_search(g_sym_graphs_ref.begin(), g_sym_graphs_ref.end(), x))
			{
				FILE* f = fopen(g_symfn.c_str(), "a");
#if GRAPH_MAX_LEN==4
				fprintf(f, "%d\t0x%llx\t0x%llx\t0x%llx\t0x%llx\t0x0\t0x0\n",
					d, x.n[0], x.n[1], x.n[2], x.n[3]);
#elif GRAPH_MAX_LEN==6
				fprintf(f, "%d\t0x%llx\t0x%llx\t0x%llx\t0x%llx\t0x%llx\t0x%llx\n",
					d, x.n[0], x.n[1], x.n[2], x.n[3], x.n[4], x.n[5]);
#else
				fprintf(f, "%d\t0x%llx\t0x%llx\n",
					d, x.n[0], x.n[1]);
#endif
				fclose(f);
				g_sym_graphs.push_back(x);
			}
		}
	}
	FILE* f=fopen(g_symcoverfn.c_str(), "a");
	if(sig.len<=64)
		fprintf(f, "0x%llx\n", sig.n[0]);
	else
		fprintf(f, "0x%llx\t0x%llx\n", sig.n[0], sig.n[1]);
	fclose(f);
}

void sym_hit_miss(uint64_t pos, b128 sig, buffer<bigint>* p, int max_len)
{
	if (!g_symcoverfn.empty())
	{
		FILE* f = fopen(g_symcoverfn.c_str(), "a");
		if (sig.len <= 64)
			fprintf(f, "0x%llx\n", sig.n[0]);
		else
			fprintf(f, "0x%llx\t0x%llx\n", sig.n[0], sig.n[1]);
		fclose(f);
	}
}

void sym_base_cover(b128 sig, uint64_t max_len)
{
	if (!g_symcoverfn.empty())
	{
		FILE* f = fopen(g_symcoverfn.c_str(), "a");
		if (sig.len <= 64)
			fprintf(f, "0x%llx\n", sig.n[0]);
		else
			fprintf(f, "0x%llx\t0x%llx\n", sig.n[0], sig.n[1]);
		fclose(f);
	}
}

int gcd(int x, int y)
{
	if (x < y)
		swap(x, y);
	while (true)
	{
		x = x%y;
		if (x == 0)
			return y;
		swap(x, y);
	}
}



int g_phi[1000];

void init_phi()
{
	for (int len = 2; len<1000; len++)
	{
		int phi = 1;
		for (int j = 2; j * 2 <= len; j++)
			//for (int j = 2; j <= 2; j++)
		{
			if (gcd(j, len + 1) != 1)
				continue;
			phi++;
		}
		g_phi[len] = phi;
	}
}

void relabel(vector<b128>& vnew, int d, const vector<bigint>& v, bool diagonal)
{
	vector<b128> vnew_thr[nthreads];
#pragma omp parallel for
	for(int thr=0; thr<nthreads; thr++)
	{
		for (int i = thr; i < v.size(); i+=nthreads)
		{
			bigint x = v[i];
			for (int j = 1; j*2 <= x.len; j++)
			{
				if (gcd(j, x.len + 1) != 1)
					continue;
				b128 y;
				y.clear();
				y.set_len(x.len);
				for (int k = 0; k < d; k++)
					if (x.bit((((k + 1)*j) % (x.len + 1)) - 1))
						y.set(k);
				y.unset_from(d);
				y.set_len(d);
				if (diagonal && (y.n[0] & 1))
				{
					y.neg();
					y.unset_from(d);
				}
				vnew_thr[thr].push_back(y);
			}
		}
	}
	for(int thr=0; thr<nthreads; thr++)
	{
		append(vnew, vnew_thr[thr]);
		vnew_thr[thr].clear();
	}
	rem_dup(vnew);
}


void relabel(vector<bigint>& vnew, bigint b, bool diagonal)
{
	vnew.resize(g_phi[b.len]);
	int pos = 0;

	vnew[pos++] = b;
	bigint x = b;
	for (int j = 2; j * 2 <= x.len; j++)
		//for (int j = 2; j <= 2; j++)
	{
		if (gcd(j, x.len + 1) != 1)
			continue;
		bigint y;
		y.clear();
		y.set_len(x.len);
		int p = j;
		for (int k = 0; k < x.len; k++)
		{
		//	if (x.bit((((k + 1)*j) % (x.len + 1)) - 1))
			if(x.bit(p-1))
				y.set(k);
			p+=j;
			if(p>=x.len+1)
				p-=x.len+1;
		}
		if (diagonal && y.bit(0))
		{
			y.neg();
			y.unset_from(y.len);
		}
		vnew[pos++] = y;
	}
}


void relabel(bigint& y, bigint x, int mul, bool diagonal)
{
	if (gcd(mul, x.len + 1) != 1)
	{
		printf("ERROR: relabel() with bad mul called (mul %d, len %d)\n", mul, x.len);
		exit(-1);
	}
	y.clear();
	y.set_len(x.len);
	int p = mul;
	for (int k = 0; k < x.len; k++)
	{
	//	if (x.bit((((k + 1)*j) % (x.len + 1)) - 1))
		if(x.bit(p-1))
			y.set(k);
		p+=mul;
		if(p>=x.len+1)
			p-=x.len+1;
	}
	if (diagonal && y.bit(0))
	{
		y.neg();
		y.unset_from(y.len);
	}
}

void relabel(vector<bigint>& vnew, vector<int>& vnew_idx, bigint b, bool diagonal)
{
	vnew.resize(g_phi[b.len]);
	vnew_idx.resize(g_phi[b.len]);
	int pos = 0;

	vnew[pos] = b;
	vnew_idx[pos] = 1;
	pos++;
	bigint x = b;
	for (int j = 2; j * 2 <= x.len; j++)
		//for (int j = 2; j <= 2; j++)
	{
		if (gcd(j, x.len + 1) != 1)
			continue;
		bigint y;
		y.clear();
		y.set_len(x.len);
		int p = j;
		for (int k = 0; k < x.len; k++)
		{
		//	if (x.bit((((k + 1)*j) % (x.len + 1)) - 1))
			if(x.bit(p-1))
				y.set(k);
			p+=j;
			if(p>=x.len+1)
				p-=x.len+1;
		}
		if (diagonal && y.bit(0))
		{
			y.neg();
			y.unset_from(y.len);
		}
		vnew[pos] = y;
		vnew_idx[pos] = j;
		pos++;
	}
	if(pos!=vnew.size())
	{
		printf("Error: relabel mismatch (%d != %d, b.len %d)\n", pos, vnew.size(), b.len);
		exit(-1);
	}
}

void relabel_old(vector<bigint>& vnew, bigint b, bool diagonal)
{
	vnew.resize(g_phi[b.len]);
	int pos = 0;

	vnew[pos++] = b;
	bigint x = b;
	for (int j = 2; j * 2 <= x.len; j++)
		//for (int j = 2; j <= 2; j++)
	{
		if (gcd(j, x.len + 1) != 1)
			continue;
		bigint y;
		y.clear();
		y.set_len(x.len);
		int p = j;
		for (int k = 0; k < x.len; k++)
		{
			if (x.bit((((k + 1)*j) % (x.len + 1)) - 1))
		//	if(x.bit(p-1))
				y.set(k);
		//	p+=j;
		//	if(p>=x.len+1)
		//		p-=x.len+1;
		}
		if (diagonal && y.bit(0))
		{
			y.neg();
			y.unset_from(y.len);
		}
		vnew[pos++] = y;
	}
}

void relabel(vector<bigint>& vnew, const vector<bigint>& vb, bool diagonal)
{
	vector<bigint> vnew_thr[nthreads];
#pragma omp parallel for
	for(int thr=0; thr<nthreads; thr++)
	{
		for (int i = thr; i < vb.size(); i+=nthreads)
		{
			bigint x = vb[i];
			if (diagonal && x.bit(0))
			{
				x.neg();
				x.unset_from(x.len);
			}
			vnew_thr[thr].push_back(x);
			for (int j = 2; j * 2 <= x.len; j++)
				//for (int j = 2; j <= 2; j++)
			{
				if (gcd(j, x.len + 1) != 1)
					continue;
				bigint y;
				y.clear();
				y.set_len(x.len);
				int p = j;
				for (int k = 0; k < x.len; k++)
				{
//					if (x.bit((((k + 1)*j) % (x.len + 1)) - 1))
//						y.set(k);
					if(x.bit(p-1))
						y.set(k);
					p+=j;
					if(p>=x.len+1)
						p-=x.len+1;
				}
				if (diagonal && y.bit(0))
				{
					y.neg();
					y.unset_from(y.len);
				}
				vnew_thr[thr].push_back(y);
			}
		}
	}

	for(int thr=0; thr<nthreads; thr++)
	{
		append(vnew, vnew_thr[thr]);
		vnew_thr[thr].clear();
	}
	rem_dup(vnew);
}

void relabel(vector<bigint>& vnew, const bigint* vb, size_t size, bool diagonal)
{
	vector<bigint> vnew_thr[nthreads];
#pragma omp parallel for
	for(int thr=0; thr<nthreads; thr++)
	{
		for (int i = thr; i < size; i+=nthreads)
		{
			bigint x = vb[i];
			if (diagonal && x.bit(0))
			{
				x.neg();
				x.unset_from(x.len);
			}
			vnew_thr[thr].push_back(x);
			for (int j = 2; j * 2 <= x.len; j++)
				//for (int j = 2; j <= 2; j++)
			{
				if (gcd(j, x.len + 1) != 1)
					continue;
				bigint y;
				y.clear();
				y.set_len(x.len);
				int p = j;
				for (int k = 0; k < x.len; k++)
				{
//					if (x.bit((((k + 1)*j) % (x.len + 1)) - 1))
//						y.set(k);
					if(x.bit(p-1))
						y.set(k);
					p+=j;
					if(p>=x.len+1)
						p-=x.len+1;
				}
				if (diagonal && y.bit(0))
				{
					y.neg();
					y.unset_from(y.len);
				}
				vnew_thr[thr].push_back(y);
			}
		}
	}

	for(int thr=0; thr<nthreads; thr++)
	{
		append(vnew, vnew_thr[thr]);
		vnew_thr[thr].clear();
	}
	rem_dup(vnew);
}

template <class T>
void writesym(string s, T& v)
{
	if (v.empty())
	{
		_unlink(s.c_str());
		return;
	}		
	FILE* f = fopen(s.c_str(), "w");
	if (f != 0)
	{
		for(auto i=v.begin(); i!=v.end(); i++)
		{
#if GRAPH_MAX_LEN==2
			fprintf(f, "%d\t0x%llx\t0x%llx\t0x%llx\t0x%llx\t0x%llx\t0x%llx\n", 
				i->len+1, i->n[0], i->n[1], 0, 0, 0, 0);
#elif GRAPH_MAX_LEN==4
			fprintf(f, "%d\t0x%llx\t0x%llx\t0x%llx\t0x%llx\t0x%llx\t0x%llx\n", 
				i->len+1, i->n[0], i->n[1], i->n[2], i->n[3], 0, 0);
#else
			fprintf(f, "%d\t0x%llx\t0x%llx\t0x%llx\t0x%llx\t0x%llx\t0x%llx\n", 
				i->len+1, i->n[0], i->n[1], i->n[2], i->n[3], i->n[4], i->n[5]);
#endif
		}
		fclose(f);
	}
}

void readsym(string s, vector<bigint>& v, int target)
{
	FILE* f = fopen(s.c_str(), "r");
	if (f != 0)
	{
		fseek(f, 0, SEEK_END);
		if(ftell(f)==0)
		{
			fclose(f);
			return;
		}
		rewind(f);
		while (!feof(f))
		{
			int dd=-1;
			uint64_t n[6];
			int fields = fscanf(f, "%d\t0x%llx\t0x%llx\t0x%llx\t0x%llx\t0x%llx\t0x%llx\n", &dd, &n[0], &n[1], &n[2], &n[3], &n[4], &n[5]);
			//printf("%d %d 0x%llx\n", fields, dd, n[0]);
			if(fields==0)
				continue;
			if(dd<0)
				continue;
			if(dd>GRAPH_MAX_LEN*64)
			{
				printf("WARNING: entry exceeding GRAPH_MAX_LEN*64 in readsym()\n");
				continue;
			}
			if(target>0 && dd<target)
				continue;
			bigint y;
			y.set_len(dd - 1);
			y.n[0] = n[0];
			y.n[1] = n[1];
#if GRAPH_MAX_LEN>=4
			y.n[2] = n[2];
			y.n[3] = n[3];
#endif
#if GRAPH_MAX_LEN>=6
			y.n[4] = n[4];
			y.n[5] = n[5];
#endif
			if(y.bit(0))
			{
				//printf("%s %d 0x%llx\n", s.c_str(), dd, n[0]);
				//y.neg();
				//y.unset_from(y.len);
				//continue;
			}
			v.push_back(y);
		}
		fclose(f);
	}
}


void rm_relabelings(vector<bigint>& v, bool diagonal)
{
	if(diagonal)
	{
		for(int i=0; i<v.size(); i++)
			if(v[i].bit(0))
			{
				v[i].neg();
				v[i].unset_from(v[i].len);
			}
	}
	rem_dup(v);
#pragma omp parallel for
	for (int thr = 0; thr<nthreads; thr++)
	{
		for (int j = thr; j<v.size(); j += nthreads)
		{
			bigint x = v[j];
			vector<bigint> vnew;
			relabel(vnew, x, diagonal);
			for (int k = 0; k<vnew.size(); k++)
			{
				size_t pos = lower_bound(v.begin(), v.begin() + j, vnew[k]) - v.begin();
				if (pos == j)
					continue;
				if (!bit_equal(v[pos], vnew[k]))
					continue;
				v[j].len = 0;
				break;
			}
		}
	}
	vector<bigint> v2;
	for (int i = 0; i<v.size(); i++)
		if (v[i].len>0)
			v2.push_back(v[i]);
	v.swap(v2);
}

void cleansym(int k0, int k1, bool check_errors, const char* secondary, int sec_min_len)
{
	g_kmax[0]=k0;
	g_kmax[1]=k1;
	vector<bigint> v, vAlt;
	readsym(format_sym_name(), v, 0);
	if(secondary != 0)
		readsym(secondary, v, sec_min_len);
	//printf("%d\n", v.size());

	if(k1!=k0)
	{
		g_kmax[0]=k1;
		g_kmax[1]=k0;	
		readsym(format_sym_name(), vAlt, 0);
	}
	rem_dup(v);
	rem_dup(vAlt);

	printf("%d + %d = %d left\n", v.size(), vAlt.size(), v.size()+vAlt.size());
	neg(vAlt);
	append(v, vAlt);
	vAlt.clear();
	
	for(int i=0; i<v.size(); i++)
	{
		bigint b_inv;
		bigint b = v[i];
		invert(b_inv, b);
		if(!bit_equal(b_inv, b))
			v[i].len=0;
//		if(v[i].len+1<240)
//			v[i].len=0;
	}
	refilter(v, 2);
	
	if (check_errors)
	{
		int errors = 0;
#pragma omp parallel for
		for (int thr = 0; thr < nthreads; thr++)
		{
			graph g;
			for (int j = thr; j < v.size(); j += nthreads)
			{
				int kmax[2] = { k0, k1 };

				bigint2 mask;
				mask.first = v[j];
				mask.second = v[j];
				mask.first.neg();
				mask.first.unset_from(mask.first.len);
				g.reset();
				bool ok = build_graph(g, mask, kmax);
				if (!ok)
				{
					printf("%d 0x%llx 0x%llx 0x%llx 0x%llx\n", v[j].len+1, v[j].n[0], v[j].n[1], v[j].n[2], v[j].n[3]);
					errors++;
					v[j].len=0;
				}
			}
		}

		printf("%d errors\n", errors);
		refilter(v, 2);
	}

	bigint maxsym;
	int maxval=0, minval=INT_MAX;
	for(int j=0; j<v.size(); j++)
	{
		if(v[j].len+1>maxval)
			maxsym=v[j];
		maxval=max(maxval, v[j].len+1);
		minval=min(minval, v[j].len+1);
	}
	vector<int> counts;
	counts.resize(maxval+1);
	rm_relabelings(v, (k0==k1));
	for(int j=0; j<v.size(); j++)
	{
		counts[v[j].len+1]++;
	}

	printf("%d left\n", v.size());
	for(int j=minval; j<=maxval; j++)
		printf("%d\t%d\n", j, counts[j]);
	printf("CirculantGraph[%d, {", maxsym.len+1);
	for(int j=0; j*2<=maxsym.len; j++)
		if(maxsym.bit(j))
			printf("%d, ", j+1);
	printf("}]\n");
	printf("nx.circulant_graph(%d, [", maxsym.len+1);
	for(int j=0; j*2<=maxsym.len; j++)
		if(maxsym.bit(j))
			printf("%d, ", j+1);
	printf("])\n");
	g_kmax[0]=k0;
	g_kmax[1]=k1;
	writesym(format_sym_name(), v);

	if(k0!=k1)
	{
		v.clear();
		g_kmax[0]=k1;
		g_kmax[1]=k0;
		writesym(format_sym_name(), v);
	}
}

void do_relabeling(int d, int target, int k0, int k1, bool initial, double fraction)
{
	uint64_t t0, t1;
	t0=__rdtsc();
	vector<bigint> sym_graphs;// .clear();
	g_max_ops=0;
	vector<b128> vsig;
	if(k0!=k1)
	{
		vector<b128> vsigAlt;
		g_kmax[0] = k1;
		g_kmax[1] = k0;
		readfile(format_name(d, target), vsigAlt, d);
		neg(vsigAlt);
		append(vsig, vsigAlt);
	}

	g_kmax[0] = k0;
	g_kmax[1] = k1;
	readfile(format_name(d, target), vsig, d);
	vector<b128> vsigRef = vsig;
	//vsig.erase(vsig.begin(), vsig.end()-5000);

	g_symfn = format_sym_name();
	g_kmax[0] = k1;
	g_kmax[1] = k0;
	string altsymfn = format_sym_name();
	g_kmax[0] = k0;
	g_kmax[1] = k1;
	g_symcoverfn = format_symcover_name(d, target);

	FILE* f=fopen(g_symfn.c_str(), "r");
	vector<uint64_t> symcover;
	readsym(g_symfn, sym_graphs, target);
	if (k0 != k1)
	{
		vector<bigint> v;
		readsym(altsymfn, v, target);
		neg(v);
		append(sym_graphs, v);
	}
	printf("%d elements in sig list\n", vsig.size());

	vector<b128> vsigExcl;
	readfile(g_symcoverfn, vsigExcl, d);

	int nEx=exclude(vsig, vsigExcl);
	printf("%d previously covered\n", nEx);

	if(initial)
	{
		vector<b128> vnew;
		relabel(vnew, d, sym_graphs, k0==k1);
//		nEx = exclude(vsig, vnew);
//		printf("%d duplicate of previously processed\n", nEx);
		printf("%d left\n", vsig.size());

		if(vsig.size()>0)
		{
			random_shuffle(vsig.begin(), vsig.end());
			if(fraction<1.0)
				vsig.resize(vsig.size()*fraction);

			relabel(g_sym_graphs_ref, sym_graphs, k0==k1);
			rem_dup(g_sym_graphs_ref);
			g_max_ops=2e6;
			g_sym_target = target;
			JobArray array;
			array.submit(&vsig[0], vsig.size(), d, target, GRAPH_MAX_LEN*64, false, 0, 0, sym_hit_detect, sym_hit_miss);
		}
	}

	printf("%d symmetric\n", sym_graphs.size());
	vector<b128> vnew;
	relabel(vnew, d, sym_graphs, k0==k1);
	printf("Found %d candidate signatures\n", vnew.size());
	exclude(vnew, vsigRef);
	printf("%d not in known list\n", vnew.size());
	g_kmax[0] = k0;
	g_kmax[1] = k1;

	int new_hits=vnew.size();

	if (vnew.size()>0)
	{
		writefile(format_name(d, target), vnew, d, true);
		writefile(format_cover_name(d, target), vnew, d, true);
	}
	t1=__rdtsc();
	FILE* ff=fopen(PATH "neighbor-search-log.txt", "a");
	fprintf(ff, "%d\t%.3f\t(relabel)\n", new_hits, (t1-t0)/g_fFreq);
	fclose(ff);
}

void list_permutations(bigint x)
{
	for (int j = 1; j * 2 <= x.len; j++)
	//for (int j = 2; j <= 2; j++)
	{
		if (gcd(j, x.len+1) != 1)				
			continue;
		bigint y;
		y.clear();
		y.set_len(x.len);
		for (int k = 0; k < x.len; k++)
			if (x.bit((((k+1)*j) % (x.len+1)) - 1))
				y.set(k);
		bool bNeg=false;
		if(y.bit(0))
		{
			bNeg=true;
			y.neg();
			y.unset_from(y.len);
		}
		printf("0x%llx 0x%llx 0x%llx 0x%llx %s\n", y.n[0], y.n[1], y.n[2], y.n[3], bNeg ? "*" : "");
		bigint inv;
		invert(inv, y);
		if (!bit_equal(y, inv))
			printf("Error: not symmetric\n");
	}
}

void op_report(int d, int len)
{
	FILE* f = fopen(PATH "op_report.txt", "a");
	fprintf(f, "%d\t%d\t", d, len);
	if (!op_counts[2].empty())
	{
		printf("signature build fail: %d calls\n", op_counts[2].size());
		op_counts[2].clear();
	}
	if (!op_counts[0].empty())
	{
		sort(op_counts[0].begin(), op_counts[0].end());
		uint64_t sum = 0;
		for (int i = 0; i<op_counts[0].size(); i++)
			sum += op_counts[0][i];
		printf("misses: %d calls, median %d, 99th percentile %d, max %d, mean %d\n", op_counts[0].size(), op_counts[0][op_counts[0].size() / 2], op_counts[0][op_counts[0].size() * 99 / 100], op_counts[0].back(), (int)(sum / op_counts[0].size()));
		fprintf(f, "%d\t%d\t%d\t%d\t%d\t", op_counts[0].size(), op_counts[0][op_counts[0].size() / 2], op_counts[0][op_counts[0].size() * 99 / 100], op_counts[0].back(), (int)(sum / op_counts[0].size()));
		op_counts[0].clear();
	}
	if (!op_counts[1].empty())
	{
		sort(op_counts[1].begin(), op_counts[1].end());
		uint64_t sum = 0;
		for (int i = 0; i<op_counts[1].size(); i++)
			sum += op_counts[1][i];
		printf("hits: %d calls, median %d, 99th percentile %d, max %d, mean %d\n", op_counts[1].size(), op_counts[1][op_counts[1].size() / 2], op_counts[1][op_counts[1].size() * 99 / 100], op_counts[1].back(), (int)(sum / op_counts[1].size()));
		fprintf(f, "%d\t%d\t%d\t%d\t%d\t", op_counts[1].size(), op_counts[1][op_counts[1].size() / 2], op_counts[1][op_counts[1].size() * 99 / 100], op_counts[1].back(), (int)(sum / op_counts[1].size()));
		op_counts[1].clear();
	}
	fprintf(f, "\n");
	fclose(f);
}

void merge(int k0, int k1, int dPrev, int d, int targetPrev, int target)
{
	printf("merge(%d->%d,%d->%d)\n", dPrev, d, targetPrev, target);
	if (dPrev < d || targetPrev < target)
	{
		printf("Error: invalid argument\n");
		exit(-1);
	}
	if (dPrev == d && targetPrev == target)
		return;

	vector<b128> v;
	g_kmax[0]=k0;
	g_kmax[1]=k1;
	readfile(format_name(dPrev,targetPrev),v,dPrev);
	printf("%d -> ", v.size());
	if (d < dPrev)
	{
		for (int i = 0; i<v.size(); i++)
		{
			v[i].unset_from(d);
			v[i].set_len(d);
		}
	}
	readfile(format_name(d,target),v,d);
	printf("%d -> ", v.size());
	rem_dup(v);
	printf("%d\n", v.size());
	writefile(format_name(d,target),v,d);
//	v.clear();
//	readfile(format_name(d,targetPrev),v,d);
	readfile(format_cover_name(d,target),v,d);
	printf("%d -> ", v.size());
	rem_dup(v);
	printf("%d\n", v.size());
	writefile(format_cover_name(d,target),v,d);
}

void fix_list(int d, int k0, int k1, int target)
{
	printf("fix_list(%d,%d,%d,%d)\n", d,k0,k1,target);
	vector<bigint> v;
	g_kmax[0]=k0;
	g_kmax[1]=k1;
	readfile(format_name(d,target),v,d);
	printf("list: %d -> ", v.size());
	rem_dup(v);
	printf(" %d\n", v.size());
	writefile(format_name(d,target),v,d);
	size_t sizePrev = v.size();
	readfile(format_cover_name(d,target),v,d);
	printf("cover: %d -> ", v.size()-sizePrev);
	rem_dup(v);
	printf(" %d\n", v.size());
	writefile(format_cover_name(d,target),v,d);
}


bool test_build(graph& g, bigint b, int k0, int k1)
{
	bigint d=b;
	d.neg();
	d.unset_from(d.len);
	g.reset();

	bigint2 m;
	m.first = d;
	m.second = b;
	int kmax[2] = { k0, k1 };
	return build_graph(g, m, kmax, false);
}


bool test_build(graph& g, const graph& gRef, bigint b, int k0, int k1)
{
	if (gRef.mask.len == 0 && gRef.mask0.len == 0)
		return test_build(g, b, k0, k1);

	bigint d = b;
	d.neg();
	d.unset_from(d.len);
	g.reset();

	int kmax[2] = { k0, k1 };

	bool first = true;
	for (int i = 0; i<b.len; i++)
	{
		if (gRef.mask.bit(i) || gRef.mask0.bit(i))
		{
			if (gRef.mask.bit(i) != b.bit(i))
			{
				printf("test_build assertion failure @ bit %d\n", i);
				gRef.mask0.print();
				gRef.mask.print();
				b.print();
				exit(-1);
			}

			continue;
		}
		if (first)
		{
			if (!extend(g, gRef, b.bit(i), kmax, b.len + 1, i, 0))
				return false;
			first = false;
		}
		else
		{
			if (!extend(g, b.bit(i), kmax, b.len + 1, i))
				return false;
		}
	}
	return true;
}


int lf_len=0;
vector<b128> g_hit_list;

void general_hit_detect(uint64_t pos, b128 sig, buffer<bigint>* p, int max_len)
{
	int i;
	if(0)
	{
		for(i=0; i<p->size(); i++)
		{
			(*p)[i].print();
			bigint x;
			invert(x, (*p)[i]);
			x.len=(*p)[i].len;
			x.print();
			printf("\n");
		}
	}
	if(lf_len>0)
	{
		if(p->size()>1)
		{
			vector<b128> v;
			for(int i=0; i<p->size(); i++)
			{
				bigint x = (*p)[i];
				x.set_len(lf_len);
				x.unset_from(lf_len);
				b128 xx;
				xx = x;
				v.push_back(xx);
			}
			rem_dup(v);
			append(g_hit_list, v);
		}
		else if(p->size()==1)
		{
			bigint x = (*p)[0];
			x.set_len(lf_len);
			x.unset_from(lf_len);
			b128 xx;
			xx = x;
			g_hit_list.push_back(xx);
		}
	}
	else
	{
		g_hit_list.push_back(sig);
	}

	b128 x = sig;
	set<bigint> cover_list, dual_cover_list;
	if (p->size()>0 && do_early_enum)
	{
		int min_len = ((*p)[0]).len;
		for (int j = 1; j<p->size(); j++)
		{
			min_len = min(min_len, ((*p)[j]).len);
		}

		for (int j = max(min_len, sig.len * 2); j<sig.len * 2 + max_early_enum * 2 && j <= GRAPH_MAX_LEN * 64; j++)
		{
			int unk_bits = (j - (sig.len * 2) + 1) / 2;
			bigint y;
			y = sig;
			y.set_len(j);
			for (int k = 0; k<sig.len; k++)
				if (y.bit(k))
					y.set(j - k - 1);
			for (int opt = 0; opt<(1 << unk_bits); opt++)
			{
				for (int k = 0; k<unk_bits; k++)
				{
					int val = (opt >> k) & 1;
					if (val)
					{
						y.set(sig.len + k);
						y.set(j - (sig.len + k) - 1);
					}
					else
					{
						y.unset(sig.len + k);
						y.unset(j - (sig.len + k) - 1);
					}
				}
				if (binary_search(g_sym_graphs_ref.begin(), g_sym_graphs_ref.end(), y))
					continue;
				bigint inv;
				invert(inv, y);
				if (!bit_equal(inv, y))
				{
					printf("Error: invert test fail\n");
					exit(-1);
				}
				sym_potential_graphs.push_back(y);
				cover_list.insert(y);
			}
		}

	}
	for (i = 0; i<p->size(); i++)
	{
		bigint x = (*p)[i];
		//x.print();
		//bigint inv;
		//invert(inv, x);
		//inv.print();

		//for(int j=x.len; j<x.len+32 && j<=GRAPH_MAX_LEN*64; j++)
		for (int j = x.len; j<x.len + 32 && j <= GRAPH_MAX_LEN * 64; j++)
		{
			bigint y;
			y = x;
			y.set_len(j);
			y.unset_from(x.len);
			vector<int> errors;
			for (int k = x.len; k<j; k++)
			{
				if (y.bit(j - k - 1))
					y.set(k);
			}
			for (int k = 0; k * 2<j; k++)
			{
				if (y.bit(k) != y.bit(j - k - 1))
					errors.push_back(k);
			}
			//printf("%d\t%d\n", j, errors.size());

			if (errors.size() <= (i == 0 ? max_errors_detect_only : max_errors))
			{
				/*
				printf("%d %d\n", j, errors.size());
				y.print();
				bigint inv;
				invert(inv, y);
				inv.print();
				*/
				for (int opt = 0; opt<(1 << errors.size()); opt++)
				{
					for (int k = 0; k<errors.size(); k++)
					{
						int val = (opt >> k) & 1;
						if (val)
						{
							y.set(errors[k]);
							y.set(j - errors[k] - 1);
						}
						else
						{
							y.unset(errors[k]);
							y.unset(j - errors[k] - 1);
						}
					}
					bigint yy = y;
					if(y.bit(0) && g_kmax[0]==g_kmax[1])
					{
						yy.neg();
						yy.unset_from(yy.len);
					}

					if (!binary_search(g_sym_graphs_ref.begin(), g_sym_graphs_ref.end(), yy) && cover_list.find(y) == cover_list.end())
					{
						sym_potential_graphs.push_back(yy);
						cover_list.insert(yy);
					}
				}
			}
		}
	}
}


int lengthen(const char* fn, int k0, int k1, int len, int lenNext, int target, double fraction)
{
	g_kmax[0]=k0;
	g_kmax[1]=k1;
		
	if(target > GRAPH_MAX_LEN*64)
	{
		printf("Recompile with higher GRAPH_MAX_LEN\n");
		exit(-1);
	}

	vector<b128> vsig, vsig0, vsigExcl, vsigExclParent;
	g_hit_detect_name=format_name(lenNext, target);
	g_base_cover_name=format_cover_name(lenNext, target);
//	string parent_cover_name=format_cover_name(lenNext-5, target);
	vector<b128> vx;
	readfile(g_hit_detect_name, vx, lenNext);
	readfile(g_base_cover_name, vx, lenNext);
	readfile(fn, vsig0, len);

	vsigExcl = vsig0;
	//readfile(parent_cover_name, vsigExclParent, lenNext-5);
	uint64_t nsigs0 = vsig0.size();
	uint64_t mul = 1<<(lenNext-len);
	vsig.resize(mul*nsigs0);
	for(uint64_t i=0; i<nsigs0; i++)
		for(uint64_t j=0; j<mul; j++)
		{
			vsig[i*mul+j]=vsig0[i];
			vsig[i*mul+j].set_len(lenNext);
			for(int k=len; k<lenNext; k++)
			if(j & (One<<(k-len)))
			{
				vsig[i*mul+j].set(k);
			}
		}
	int nEx = exclude(vsig, vx);
	printf("%d previously covered\n", nEx);

	printf("%d candidates\n", vsig.size());
	shuffle(vsig);
	if(fraction < 1.0)
		vsig.resize(vsig.size()*fraction);
	int stop_size=target+1;
	JobArray array;
	g_hit_list.clear();
	sym_potential_graphs.clear();
	array.submit(&vsig[0], vsig.size(), lenNext, target, stop_size, true, 0, 0, general_hit_detect);
	writefile(g_hit_detect_name, g_hit_list, lenNext, true);
	writefile(g_base_cover_name, vsig, lenNext, true);
	printf("%d -> %d\n", vsig.size(), array.sum_populations[target]);
	vector<bigint> sym;
	finish_cyclic_search(k0, k1, g_hit_list, sym, target);
	FILE* ff=fopen(PATH "neighbor-search-log.txt", "a");
	fprintf(ff, "%d\t%d\t%.3f\t(lengthen)\n", vsig.size(), array.sum_populations[target], array.realtime/g_fFreq);
	fclose(ff);
	array.print_populations();
	array.print_cputime();
	return array.sum_populations[target];
}

int lengthen(int k0, int k1, int len, int lenNext, int targetPrev, int target, double fraction)
{
	g_kmax[0]=k0;
	g_kmax[1]=k1;
	return lengthen(format_name(len, targetPrev).c_str(), k0, k1, len, lenNext, target, fraction);
}

int lengthen_memory(const vector<b128>& v, vector<b128>& out, int k0, int k1, int len, int lenNext, int target)
{
	g_kmax[0]=k0;
	g_kmax[1]=k1;
		
	if(target > GRAPH_MAX_LEN*64)
	{
		printf("Recompile with higher GRAPH_MAX_LEN\n");
		exit(-1);
	}

	uint64_t mul = 1<<(lenNext-len);
	vector<b128> vsig;
	int nsigs0=v.size();
	vsig.resize(mul*nsigs0);
	for(uint64_t i=0; i<nsigs0; i++)
		for(uint64_t j=0; j<mul; j++)
		{
			vsig[i*mul+j]=v[i];
			vsig[i*mul+j].set_len(lenNext);
			for(int k=len; k<lenNext; k++)
			if(j & (One<<(k-len)))
			{
				vsig[i*mul+j].set(k);
			}
		}
	printf("%d candidates\n", vsig.size());
	shuffle(vsig);
	int stop_size=target+1;
	g_hit_list.clear();
	sym_potential_graphs.clear();
	JobArray array;
	array.submit(&vsig[0], vsig.size(), lenNext, target, stop_size, true, 0, 0, general_hit_detect);
	array.print_cputime();
	vector<bigint> sym;
	finish_cyclic_search(k0, k1, g_hit_list, sym, target);
	out = g_hit_list;
	return out.size();
}

int narrow(int len, int targetPrev, int target, double fraction)
{
	printf("narrow(%d, %d, %d)\n", len, targetPrev, target);
	if(target > GRAPH_MAX_LEN*64)
	{
		printf("Error: recompile with GRAPH_MAX_LEN=4 or higher\n");
		exit(-1);
	}

	string fnIn = format_name(len, targetPrev);
	string fnOut = format_name(len, target);
	string fnInCover = format_cover_name(len, targetPrev);
	//string fnInCoverParent = format_cover_name(len-5, targetPrev);

	string fnOutCover = format_cover_name(len, target);
	//string fnOutCoverParent = format_cover_name(len-5, target);

	uint64_t i=0;
	vector<b128> vsig;
	readfile(fnIn, vsig, len);
	vector<b128> vsigNextShort;
	//, vsigNextParent;
	
	readfile(fnOutCover, vsigNextShort, len);
	//readfile(fnOutCoverParent, vsigNextParent, len-5);

	vector<b128> vsigInCover;
	//, vsigInCoverParent;
	readfile(fnInCover, vsigInCover, len);
	//readfile(fnInCoverParent, vsigInCoverParent, len-5);

	int dup = rem_dup(vsig);
	printf("%d duplicate\n", dup);
	rem_dup(vsigNextShort);

	vector<b128> check;
	readfile(fnOut, check, len);
	int nEx;
	nEx = exclude(vsig, check);
	printf("%d already on the list\n", nEx);

	//append(vsigNextShort, check);
	//check = vsigNextShort;
	exclude(check, vsigNextShort);
	//exclude(check, vsigNextParent, len-5);
	if(!check.empty())
	{
		printf("WARNING: %d entries in target list are not in target cover list\n", check.size());
		append(vsigNextShort, check);
	}

	nEx = exclude(vsig, vsigNextShort);
	printf("%d excluded by target cover\n", nEx);
	//nEx = exclude(vsig, vsigNextParent, len-5);
	//printf("%d excluded by target parent cover\n", nEx);
	printf("%d remaining\n", vsig.size());


	append(vsigNextShort, vsigInCover);
	//append(vsigNextParent, vsigInCoverParent);
	rem_dup(vsigNextShort);
	//rem_dup(vsigNextParent);
	exclude(vsigNextShort, vsig);
	random_shuffle(vsig.begin(), vsig.end());
	if(fraction < 1.0)
		vsig.resize(vsig.size()*fraction);
	append(vsigNextShort, vsig);

	int retval=0;
	if(vsig.size()>0)
	{
		g_hit_detect_name = fnOut;

		JobArray array;
		bool detect_only =true;
		int stop_size=detect_only ? target+1 : GRAPH_MAX_LEN*64;
		g_hit_list.clear();
		sym_potential_graphs.clear();
		array.submit(&vsig[0], vsig.size(), len, target, stop_size, detect_only, 0, 0, general_hit_detect, 0);
		array.print_populations();
		array.print_cputime();
		printf("%d -> %d\n", vsig.size(), array.sum_populations[target]);
		FILE* ff=fopen(PATH "neighbor-search-log.txt", "a");
		fprintf(ff, "%d\t%d\t%.3f\t(oldnarrow)\n", vsig.size(), array.sum_populations[target], array.realtime/g_fFreq);
		fclose(ff);
		retval=array.sum_populations[target];
		writefile(fnOut, g_hit_list, len, true);
		vector<bigint> sym;
		finish_cyclic_search(g_kmax[0], g_kmax[1], g_hit_list, sym, target);
	}
	printf("vsigNextShort: %d elements, len %d\n", vsigNextShort.size(), len);
	writefile(fnOutCover.c_str(), vsigNextShort, len);
	//writefile(fnOutCoverParent.c_str(), vsigInCoverParent, len);

	return retval;
}


int lengthen_fast(const char* fn, int k0, int k1, int len, int lenNext, int target, double fraction, bool detect_only)
{
	g_kmax[0]=k0;
	g_kmax[1]=k1;
		
	if(target > GRAPH_MAX_LEN*64)
	{
		printf("Recompile with higher GRAPH_MAX_LEN\n");
		exit(-1);
	}

	vector<b128> vsig, vsigExcl;
	string hit_detect_name=format_name(lenNext, target);
	string hit_detect_name2=format_cover_name(lenNext, target);;
	g_base_cover_name.clear();
	vector<b128> vx;
	readfile(hit_detect_name, vx, lenNext);
	readfile(fn, vsig, len);
	int nEx = 0;//exclude(vsig, vx, len);
	printf("%d already on the list\n", nEx);
	printf("%d candidates\n", vsig.size());
	random_shuffle(vsig.begin(), vsig.end());
	if(fraction < 1.0)
		vsig.resize(vsig.size()*fraction);
	JobArray array;
	lf_len=lenNext;
	array.submit(&vsig[0], vsig.size(), len, target, target, detect_only, 0, 0, general_hit_detect);
	rem_dup(g_hit_list);
	exclude(g_hit_list, vx);
	writefile(hit_detect_name, g_hit_list, lenNext, true);
	writefile(hit_detect_name2, g_hit_list, lenNext, true);
	printf("%d -> %d\n", vsig.size(), g_hit_list.size());
	FILE* ff=fopen(PATH "neighbor-search-log.txt", "a");
	fprintf(ff, "%d\t%d\t%.3f\t(lengthen)\n", vsig.size(), array.sum_populations[target], array.realtime/g_fFreq);
	fclose(ff);
	array.print_populations();
	array.print_cputime();
	vector<bigint> sym;
	finish_cyclic_search(k0, k1, g_hit_list, sym, target);
	return array.sum_populations[target];
}

int lengthen_fast(int k0, int k1, int len, int lenNext, int targetPrev, int target, double fraction, bool detect_only)
{
	g_kmax[0]=k0;
	g_kmax[1]=k1;
	return lengthen_fast(format_name(len, targetPrev).c_str(), k0, k1, len, lenNext, target, fraction, detect_only);
}


bool g_recursive_cyclic_search = true;

void finish_cyclic_search(int k0, int k1, const vector<b128>& hits, vector<bigint>& sym, int min_save)
{
	g_kmax[0] = k0;
	g_kmax[1] = k1;
	g_symfn = format_sym_name();
	while(true)
	{
		vector<bigint> out;
		printf("%d potential symmetric\n", sym_potential_graphs.size());
		rem_dup(sym_potential_graphs);
		printf("%d unique\n", sym_potential_graphs.size());
		uint64_t t1, t2;
		t1 = __rdtsc();
	#pragma omp parallel for
		for (int thr = 0; thr<nthreads; thr++)
		{
			graph g;
			for (int i = thr; i<sym_potential_graphs.size(); i += nthreads)
			{
				g.reset();
				bigint y = sym_potential_graphs[i];
				bigint2 mask;
				mask.first = y;
				mask.second = y;
				mask.first.neg();
				mask.first.unset_from(y.len);
				bool ok = build_graph(g, mask, g_kmax);
				if (!ok)
					sym_potential_graphs[i].set_len(0);
			}
		}
		t2 = __rdtsc();

		refilter(sym_potential_graphs, 2);
		printf("%d survive build\n", sym_potential_graphs.size());

		rm_relabelings(sym_potential_graphs, (k0==k1));
		printf("%d non isomorphic\n", sym_potential_graphs.size());

		out.swap(sym_potential_graphs);
		if(out.size()==0)
			break;
		for (int64_t i = 0; i<out.size(); i++)
		{
			vector<bigint> vnew;
			relabel(vnew, out[i], (k0==k1));
			append(g_sym_graphs_ref, vnew);
		}
		rem_dup(g_sym_graphs_ref);
#if 1
		if(g_recursive_cyclic_search)
		{
			for(int i=0; i<out.size(); i++)
				for(int j=0; j*2<=out[i].len; j++)
				{
					for(int k=j; k*2<=out[i].len; k++)
					{
						bigint x = out[i];
						x.flip(j);
						if(x.len-j-1!=j)
							x.flip(x.len-j-1);
						if(k!=j)
						{
							x.flip(k);
							if(x.len-k-1!=k)
								x.flip(x.len-k-1);
						}
						sym_potential_graphs.push_back(x);
					}
				}
			rem_dup(sym_potential_graphs);
		//rm_relabelings(sym_potential_graphs, (k0==k1));
			exclude(sym_potential_graphs, g_sym_graphs_ref);
		}
#endif
		FILE* f = fopen(g_symfn.c_str(), "a");
		for (int64_t i = 0; i<out.size(); i++)
		{
			bigint y = out[i];
			int j = y.len;
			if (y.len + 1 >= min_save)
			{
				//g_sym_graphs_ref.push_back(y);
	#if GRAPH_MAX_LEN==4
				fprintf(f, "%d\t0x%llx\t0x%llx\t0x%llx\t0x%llx\t0x0\t0x0\n",
					j + 1, y.n[0], y.n[1], y.n[2], y.n[3]);
	#elif GRAPH_MAX_LEN==6
				fprintf(f, "%d\t0x%llx\t0x%llx\t0x%llx\t0x%llx\t0x%llx\t0x%llx\n",
					j + 1, y.n[0], y.n[1], y.n[2], y.n[3], y.n[4], y.n[5]);
	#else
				fprintf(f, "%d\t0x%llx\t0x%llx\t0x0\t0x0\t0x0\t0x0\n",
					j + 1, y.n[0], y.n[1]);
	#endif
			}
			sym.push_back(y);
		}
		fclose(f);

		if(!g_recursive_cyclic_search)
			break;
	}

	rem_dup(g_sym_graphs_ref);
}

void mc_memory(int len, int k0, int k1, int target, uint64_t nbases, vector<b128>& hits, vector<bigint>& sym, int min_save)
{
	if (GRAPH_MAX_LEN == 2 && target>129)
	{
		printf("Recompile with GRAPH_MAX_LEN=4\n");
		exit(-1);
	}
	g_hit_list.clear();
	sym_potential_graphs.clear();

	vector<b128> vsig;
	srand(__rdtsc());
	vsig.resize(nbases);
	g_kmax[0] = k0;
	g_kmax[1] = k1;

#pragma omp parallel for
	for (int64_t thr = 0; thr<24; thr++)
	{
		int64_t low = vsig.size()*thr / 24;
		int64_t high = vsig.size()*(thr + 1) / 24;
		if (len < 64)
		{
			for (int64_t i = low; i<high; i++)
			{
				vsig[i].clear();
				vsig[i].set_len(len);
				vsig[i].n[0] = rand64() & ((One << len) - 1);
			}
		}
		else if(len==64)
		{
			for (int64_t i = low; i<high; i++)
			{
				vsig[i].clear();
				vsig[i].set_len(len);
				vsig[i].n[0] = rand64();
			}
		}
		else
		{
			for (int64_t i = low; i<high; i++)
			{
				vsig[i].clear();
				vsig[i].set_len(len);
				vsig[i].n[0] = rand64() & (-2ll);
				vsig[i].n[1] = rand64() & ((One << (len - 64)) - 1);
			}
		}
	}
	JobArray array;
	bool detect_only = mc_detect_only;

	int stop_size = detect_only ? target + 1 : GRAPH_MAX_LEN * 64;
	array.submit(&vsig[0], vsig.size(), len, target, stop_size, detect_only, 0, 0, general_hit_detect);
	array.print_populations();
	array.print_cputime();
	printf("%d -> %d\n", vsig.size(), g_hit_list.size());
	hits = g_hit_list;
	finish_cyclic_search(k0, k1, g_hit_list, sym, min_save);
}

void mc_guided(int len, int k0, int k1, int target, uint64_t nbases, 
					  vector<b128>& known, int flips,
					  vector<b128>& hits, vector<bigint>& sym, int min_save)
{
	if (GRAPH_MAX_LEN == 2 && target>129)
	{
		printf("Recompile with GRAPH_MAX_LEN=4\n");
		exit(-1);
	}
	g_hit_list.clear();
	sym_potential_graphs.clear();

	vector<b128> vsig;
	srand(__rdtsc());
	vsig.resize(nbases);
	g_kmax[0] = k0;
	g_kmax[1] = k1;

#pragma omp parallel for
	for (int64_t thr = 0; thr<24; thr++)
	{
		int64_t low = vsig.size()*thr / 24;
		int64_t high = vsig.size()*(thr + 1) / 24;
		
		for (int64_t i = low; i<high; i++)
		{
			vsig[i]=known[rand64() % known.size()];
			b128 b;
			b.clear();
			b.set_len(len);
			while(true)
			{
				uint64_t rand_val=0;
				b.clear();
				for(int j=0; j<flips; j++)
				{
					if((j&3)==0)
						rand_val=rand64();
					b.set((rand_val >> ((j&3)*2)) % len);
				}
				if(popcnt(b)==flips)
					break;
			}
			vsig[i]=known[rand64() % known.size()]^b;
		}
	}
	exclude(vsig, known);
	shuffle(vsig);

	JobArray array;
	bool detect_only = mc_detect_only;

	int stop_size = detect_only ? target + 1 : GRAPH_MAX_LEN * 64;
	array.submit(&vsig[0], vsig.size(), len, target, stop_size, detect_only, 0, 0, general_hit_detect);
	array.print_populations();
	array.print_cputime();
	//exclude(g_hit_list, known);
	int true_hits=0;
	for(int i=0; i<g_hit_list.size(); i++)
	{
		bool good=true;
		for(int j=0; j<known.size(); j++)
		{
			if(popcnt(g_hit_list[i]^known[j])<flips)
				good=false;
		}
		if(good)
			true_hits++;
	}
	printf("%d -> %d -> %d\n", vsig.size(), g_hit_list.size(), true_hits);
	writefile(format_name(len,target), g_hit_list, len, true);
	writefile(format_cover_name(len,target), g_hit_list, len, true);
	hits = g_hit_list;
	finish_cyclic_search(k0, k1, g_hit_list, sym, min_save);
}

void list_neighbors(vector<b128>& cand, const vector<b128>& hits, size_t d, bool diagonal)
{
	size_t s0 = cand.size();
	cand.resize(s0 + d*hits.size());
	for (size_t i = 0; i < hits.size(); i++)
		for (size_t j = 0; j < d; j++)
		{
			cand[s0 + i*d + j] = hits[i];
			cand[s0 + i*d + j].flip(j);
			if (diagonal && cand[s0 + i*d + j].bit(0))
			{
				cand[s0 + i*d + j].neg();
				cand[s0 + i*d + j].unset_from(d);
			}
		}
}

void neighbor_depth_search(int len, int k0, int k1, int target, const vector<b128>& vsig0, vector<b128>& vexcl, 
	vector<b128>& hits, vector<bigint>& sym, double share, int min_save)
{
	vector<b128> vsig;
	list_neighbors(vsig, vsig0, len, k0 == k1);
	size_t nEx = rem_dup(vsig);
	exclude(vsig, vexcl);

	g_hit_list.clear();
	sym_potential_graphs.clear();

	JobArray array;
	random_shuffle(vsig.begin(), vsig.end());
	if (share < 1.0)
		vsig.resize(vsig.size()*share);
	append(vexcl, vsig);
	array.submit(&vsig[0], vsig.size(), len, target, target + 1, true, 0, 0, general_hit_detect);
	printf("%d -> %d (%.3f)\n", vsig.size(), g_hit_list.size(), array.realtime / g_fFreq);
	hits=g_hit_list;
	finish_cyclic_search(k0, k1, g_hit_list, sym, min_save);
}

void in_memory_neighbor_search_fractional(int d, int target0, double target_share, 
	vector<b128>& v, vector<b128>& cover,
	vector<bigint>& sym, int min_save)
{
	rem_dup(v);
	rem_dup(cover);
	if (v.empty())
		return;

	vector<b128> vCompl;
	vector<vector<b128> > vv;
	clusterize(vv, v);
	printf("%d clusters\n", vv.size());
	vector<int> sizes;
	for (int i = 0; i<vv.size(); i++)
		sizes.push_back(vv[i].size());
	sort(sizes.begin(), sizes.end());
	printf("Median %d, mean %.1f\n", sizes[sizes.size() / 2], double(v.size()) / vv.size());
	vector<vector<b128> > vv2;
	sort(v.begin(), v.end());
	vv2.resize(vv.size());
	for (int i = 0; i < vv.size(); i++)
	{
		for (int k = 0; k < vv[i].size(); k++)
		{
			if (is_incomplete(vv[i][k], cover, d))
				vv2[i].push_back(vv[i][k]);
			else
				vCompl.push_back(vv[i][k]);
		}
	}

	//	if (vInc.size() == 0)
	//		return;

	vv.clear();
	for (int i = 0; i < vv2.size(); i++)
		if (vv2[i].size()>0)
			vv.push_back(vv2[i]);
	printf("%d incomplete clusters\n", vv.size());
	if (vv.size() == 0)
		return;
	size_t i;
	vector<b128> cand;
	cand.resize(vv.size());
	for (i = 0; i < vv.size(); i++)
	{
		cand[i] = vv[i][rand64() % vv[i].size()];
		//	printf("%d %d 0x%llx\n", i, vv[i].size(), cand[i]);
	}
	/*
	printf("%d %d\n", cover.size(), coverParent.size());
	for (i = 0; i < coverParent.size(); i++)
	printf("0x%llx\n", coverParent[i]);
	for (i = 0; i < vv[0].size(); i++)
	printf("0x%llx\t%d\n", vv[0][i], is_incomplete(vv[0][i], cover, coverParent, d));
	*/
	vector<vector<b128> > vnewPrev = vv;
	const vector<vector<b128> > vv0 = vv;

	vector<uint64_t> cand_idx;
	for (i = 0; i < vv.size(); i++)
		cand_idx.push_back(i);
	for (int pass = 0; pass<50; pass++)
	{
		vector<b128> coverSave = cover;
		vector<b128> v;
		neighbor_depth_search(d, g_kmax[0], g_kmax[1], target0, cand, cover, v, sym, target_share, min_save);
		if (v.empty())
			break;
		vector<int> new_hits;
		new_hits.resize(vv.size());
		for (i = 0; i < new_hits.size(); i++)
			new_hits[i] = 0;

		vector<vector<b128> > vnew;
		vnew.resize(vv.size());

		vector<int> indexes;
		indexes.resize(v.size());

#pragma omp parallel for
		for (int thr = 0; thr < nthreads; thr++)
		{
			int low = v.size()*thr / nthreads;
			int high = v.size()*(thr + 1) / nthreads;
			for (int i = low; i < high; i++)
			{
				indexes[i] = -1;
				b128 neg_vi = v[i];
				neg_vi.neg();
				neg_vi.unset_from(neg_vi.len);

				for (size_t j = 0; j < cand.size(); j++)
				{
					if (dist(v[i], cand[j]) <= 1)
					{
						//		printf("0x%llx 0x%llx: %d %d\n", v[i], cand[j], j, cand_idx[j]);
						indexes[i] = cand_idx[j];
						break;
					}
					if(g_kmax[0]==g_kmax[1] && dist(neg_vi, cand[j]) <= 1)
					{
						indexes[i] = cand_idx[j];
						break;
					}
				}
			}
		}
		for (i = 0; i < v.size(); i++)
		{
			if (indexes[i] < 0)
				printf("ERROR: did not find clusters for 0x%llx\n", v[i]);
			else
			{
				//				printf("New hit: 0x%llx (idx %d, distance %d %d)\n",
				//					v[i], indexes[i], dist(v[i], vv0[indexes[i]]), dist(v[i], vv[indexes[i]]));			
				vv[indexes[i]].push_back(v[i]);
				new_hits[indexes[i]]++;
				vnew[indexes[i]].push_back(v[i]);
			}
		}
		//		cand.clear();
		//		cand.resize(vv.size());

		rem_dup(cover);
		for (i = 0; i < vv.size(); i++)
		{
			vector<b128> t;
			printf("%d: %d %d => ", i, vnew[i].size(), vnewPrev[i].size());
			for (int j = 0; j < vnew[i].size(); j++)
			{
				if (is_incomplete(vnew[i][j], cover, d))
					t.push_back(vnew[i][j]);
			}
			vnew[i].swap(t);
			t.clear();
			for (int j = 0; j < vnewPrev[i].size(); j++)
			{
				if (is_incomplete(vnewPrev[i][j], cover, d))
					t.push_back(vnewPrev[i][j]);
			}
			vnewPrev[i].swap(t);
			printf(" %d %d\n", vnew[i].size(), vnewPrev[i].size());

			if (vnew[i].size()>0)
				vnewPrev[i] = vnew[i];
			else
				vnew[i] = vnewPrev[i];
		}
		cand_idx.clear();
		cand.clear();
		int nInc = 0;
		int sum_dist = 0;
		for (i = 0; i < vv.size(); i++)
		{
			if (vnew[i].size() > 0)
			{
				int max_dist = 0;
				for (int j = 0; j<vnew[i].size(); j++)
					max_dist = max(max_dist, dist(vnew[i][j], vv0[i]));
				vector<b128> subset;
				for (int j = 0; j<vnew[i].size(); j++)
					if (dist(vnew[i][j], vv0[i]) == max_dist)
						subset.push_back(vnew[i][j]);
				int pos = rand64() % subset.size();
				cand.push_back(subset[pos]);
				cand_idx.push_back(i);//[cand.size() - 1] = i;
				sum_dist += dist(cand.back(), vv0[i]);
				nInc++;
				if (target_share>0.5 && subset.size()>1)
				{
					int pos2 = rand64() % (subset.size() - 1);
					if (pos2 >= pos)
						pos2++;
					cand.push_back(subset[pos2]);
					cand_idx.push_back(i);
					sum_dist += dist(cand.back(), vv0[i]);
				}
				//				printf("%d distance %d\n", i, dist(cand.back(), vv0[i]));
			}
		}
		target_share *= double(vv.size()) / double(v.size());
		if (target_share < 0.01)
			target_share = 0.01;
		if (target_share > 1.0)
			target_share = 1.0;
		if (nInc == 0)
			break;
		if (nInc < vv0.size() / 4)
			break;
		printf("%d -> %d (new target %f, %d incomplete, mean dist %.3f)\n", vv.size(), v.size(), target_share, nInc, double(sum_dist) / cand.size());
		vector<b128> temp = vCompl;
		for (i = 0; i < vv.size(); i++)
			append(temp, vv[i]);
		writefile(format_name(d, target0), temp, d, false);
		writefile(format_cover_name(d, target0), cover, d, false);
	}
	v = vCompl;
	for (i = 0; i < vv.size(); i++)
		append(v, vv[i]);
}

void in_memory_neighbor_search(int d, int target0, int cap, vector<b128>& v, vector<b128>& cover, vector<bigint>& sym, int min_save)
{
	rem_dup(v);
	vector<vector<b128> > vv;
	clusterize(vv, v);
	vector<int> states;
	states.resize(vv.size());
	size_t i;
	for (i = 0; i < vv.size(); i++)
		states[i] = (vv[i].size() >= cap) ? 2 : 0;
	while (true)
	{
		vector<b128> cand;
		vector<int> indexes;
		for (i = 0; i < vv.size(); i++)
		{
			if (states[i] != 0)
				continue;
			append(cand, vv[i]);
			for (int j = 0; j < vv[i].size(); j++)
				indexes.push_back(i);
		}
		if (cand.empty())
			break;
		vector<b128> v;
		neighbor_depth_search(d, g_kmax[0], g_kmax[1], target0, cand, cover, v, sym, 0.1, min_save);
		if (v.empty())
			break;
		vector<int> new_hits;
		new_hits.resize(vv.size());
		for (i = 0; i < new_hits.size(); i++)
			new_hits[i] = 0;

		indexes.resize(v.size());
#pragma omp parallel for
		for (int thr = 0; thr < nthreads; thr++)
		{
			int low = v.size()*thr / nthreads;
			int high = v.size()*(thr + 1) / nthreads;
			for (int i = low; i < high; i++)
			{
				indexes[i] = -1;
				b128 neg_vi = v[i];
				neg_vi.neg();
				neg_vi.unset_from(neg_vi.len);
				for (size_t j = 0; j < vv.size(); j++)
				{
					if (dist(v[i], vv[j]) <= 1)
					{
						indexes[i] = j;
						break;
					}
					if(g_kmax[0]==g_kmax[1] && dist(neg_vi, vv[j]) <= 1)
					{
						indexes[i] = j;
						break;
					}
				}
			}
		}
		for (i = 0; i < v.size(); i++)
		{
			if (indexes[i] < 0)
				printf("ERROR: did not find clusters for 0x%llx\n", v[i]);
			else
			{
				vv[indexes[i]].push_back(v[i]);
				new_hits[indexes[i]]++;
			}
		}

		for (i = 0; i < vv.size(); i++)
		{
			if (new_hits[i] == 0)
				states[i] = 1; // full
			if (vv[i].size() >= cap)
				states[i] = 2;
		}

		vector<b128> temp;
		for (i = 0; i < vv.size(); i++)
			append(temp, vv[i]);
		writefile(format_name(d, target0), temp, d, false);
		writefile(format_cover_name(d, target0), cover, d, false);
	}
	v.clear();
	for (i = 0; i < vv.size(); i++)
		append(v, vv[i]);
}

int narrow_memory(vector<b128>& vsig, int len, int targetPrev, int target, vector<b128>& results, vector<bigint>& sym, int min_save)
{
	if (!quiet)
		printf("narrow(%d, %d, %d)\n", len, targetPrev, target);
	int nEx = rem_dup(vsig);
	if (nEx != 0)
		printf("%d duplicate\n", nEx);
	shuffle(vsig);
	g_hit_list.clear();
	sym_potential_graphs.clear();
	JobArray array;
	int stop_size = target + 1;
	array.submit(&vsig[0], vsig.size(), len, target, stop_size, true, 0, 0, general_hit_detect, 0);
	if (!quiet)
		array.print_cputime();
	int hits = array.sum_populations[target];
	double tm = array.realtime / g_fFreq;
	if (!quiet)
		printf("%d -> %d\n", vsig.size(), hits);
	finish_cyclic_search(g_kmax[0], g_kmax[1], g_hit_list, sym, min_save);
	results.swap(g_hit_list);
	return results.size();
}

void in_memory_search(int d, int k0, int k1, int target0, int target1, vector<b128>& v, vector<b128>& cover, vector<bigint>& sym, int min_save, bool sparse, int cap)
{
	size_t size_initial = v.size();
	uint64_t t0, t1;
	t0 = __rdtsc();
	g_kmax[0] = k0;
	g_kmax[1] = k1;
	//g_max_ops = 1e4;

	if (sparse)
		in_memory_neighbor_search_fractional(d, target0, 0.1, v, cover, sym, min_save);
	else
		in_memory_neighbor_search(d, target0, cap, v, cover, sym, min_save);
	op_report(d, target0);

	if (target1 != target0)
	{
		vector<b128> v2;
		narrow_memory(v, d, target0, target1, v2, sym, min_save);
		v.swap(v2);
		op_report(d, target1);
	}
	writefile(format_name(d, target1), v, d, false);
	writefile(format_cover_name(d, target1), cover, d, false);

	t1 = __rdtsc();
	FILE* f = fopen(PATH "log_inmem.txt", "a");
	fprintf(f, "%d\t%d\t%d\t%d\t%.3f\n", target0, target1, size_initial, v.size(), (t1 - t0) / g_fFreq);
	fclose(f);
}

void in_memory_search_sparse(int k0, int k1, int d, int target)
{
	vector<b128> v, cover;
	vector<bigint> sym;
	g_kmax[0] = k0;
	g_kmax[1] = k1;
	readfile(format_name(d, target), v, d);
	readfile(format_cover_name(d, target), cover, d);
	in_memory_search(d, k0, k1, target, target, v, cover, sym, target, true);
	op_report(d, target);
}

void in_memory_search(int k0, int k1, int d, int target, int cap)
{
	vector<b128> v, cover;
	vector<bigint> sym;
	g_kmax[0] = k0;
	g_kmax[1] = k1;
	readfile(format_name(d, target), v, d);
	readfile(format_cover_name(d, target), cover, d);
	in_memory_search(d, k0, k1, target, target, v, cover, sym, target, false, cap);
	op_report(d, target);
}

void do_relabeling_new(int d, int target, int k0, int k1, double fraction)
{
	uint64_t t0, t1;
	t0=__rdtsc();
	vector<bigint> sym_graphs;// .clear();
	g_max_ops=0;
	vector<b128> vsig;
	if(k0!=k1)
	{
		vector<b128> vsigAlt;
		g_kmax[0] = k1;
		g_kmax[1] = k0;
		readfile(format_name(d, target), vsigAlt, d);
		neg(vsigAlt);
		append(vsig, vsigAlt);
	}

	g_kmax[0] = k0;
	g_kmax[1] = k1;
	readfile(format_name(d, target), vsig, d);
	vector<b128> vsigRef = vsig;

	g_symfn = format_sym_name();
	g_symcoverfn = format_symcover_name(d, target);

	FILE* f=fopen(g_symfn.c_str(), "r");
	vector<uint64_t> symcover;
	readsym(g_symfn, sym_graphs, target);
	printf("%d elements in sig list\n", vsig.size());

	vector<b128> vsigExcl;
	readfile(g_symcoverfn, vsigExcl, d);

	int nEx=exclude(vsig, vsigExcl);
	printf("%d previously covered\n", nEx);
/*
	vector<b128> vnew;
	relabel(vnew, d, sym_graphs, k0==k1);
	nEx = exclude(vsig, vnew);
	printf("%d duplicate of previously processed\n", nEx);
	printf("%d left\n", vsig.size());
*/
	if(vsig.size()==0)
		return;
	
	random_shuffle(vsig.begin(), vsig.end());
	if(fraction<1.0)
		vsig.resize(vsig.size()*fraction);

	relabel(g_sym_graphs_ref, sym_graphs, k0==k1);
	rem_dup(g_sym_graphs_ref);
	g_hit_list.clear();
	g_max_ops=1e5;
	JobArray array;
	array.submit(&vsig[0], vsig.size(), d, target, target+1, true, 0, sym_base_cover, general_hit_detect, 0);
	sym_graphs.clear();
	finish_cyclic_search(k0, k1, g_hit_list, sym_graphs, target);
	printf("%d symmetric\n", sym_graphs.size());
}

struct multilevel_state
{
	bool variable_config; // currently broken
	int max_ops;
	int default_tests_per_pass;
	multilevel_state() : diagonal(false), variable_config(true), max_ops(1e5), default_tests_per_pass(1e3), hitsPrev(0), candPrev(0) {}

	vector<b128> cand, hits, cover;
	vector<bigint> sym, sym2;
	int d;
	bool diagonal;
	int start_size;
	int max_len;
	void state_rebuild();
	void single_pass(int k0, int k1, int min_log_len, int stop_len);

	size_t hitsPrev;
	size_t candPrev;
	double timePrev;
};



void multilevel_state::state_rebuild()
{
	refilter(sym, start_size);
	hits.clear();
	relabel(hits, d, sym, diagonal);
	rem_dup(hits);
	printf("%d signatures, %d cyclic remaining\n", hits.size(), sym.size());
	cover = hits;
	cand.clear();
	list_neighbors(cand, hits, d, diagonal);
	rem_dup(cand);
	exclude(cand, cover);
	shuffle(cand);
}

uint64_t g_ml_tests[GRAPH_MAX_LEN * 64 + 1];
double g_ml_times[GRAPH_MAX_LEN * 64 + 1];


void multilevel_state::single_pass(int k0, int k1, int min_log_len, int stop_len)
{
	diagonal = (k0 == k1);
	vector<bigint> sym_new;
	g_max_ops = max_ops;
	//rem_dup(cover);
	//rem_dup(cand);
	//uint64_t t1, t2;
	//t1=__rdtsc();
	//random_shuffle(cand.begin(), cand.end());
	//shuffle(cand);
	//t2=__rdtsc();
	//printf("Shuffled %d elements in %.3f s (%e elements/s)\n", cand.size(), (t2-t1)/g_fFreq, cand.size() / ( (t2-t1)/g_fFreq ));
	//exclude(cand, cover);

	double ft0 = __rdtsc() / g_fFreq;
	printf(" %d hits, %lld candidates remaining", hits.size(), cand.size());
	if(candPrev>0)
	{
		if(cand.size() >= candPrev)
		{
			size_t repeats = cand.size() / default_tests_per_pass;
			printf(" (exponential phase, %.1f days)", repeats*(ft0-timePrev)/86400);
		}
		else
		{
			size_t delta = candPrev - cand.size();
			size_t newhits = hits.size()-hitsPrev;
			size_t repeats = cand.size() / delta;
			printf(" (project %d hits, %.1f days)", newhits*repeats, repeats*(ft0-timePrev)/86400);
		}
	}
	printf("\n");
	hitsPrev = hits.size();
	candPrev = cand.size();
	timePrev = ft0;

	if (cand.size() == 0)
		return;

	//printf("0x%llx 0x%llx %d\n", cand[0].n[0], cand[0].n[1], cand[0].len);
	//printf("0x%llx 0x%llx %d\n", cand2[0].n[0], cand2[0].n[1], cand2[0].len);
	double expect_rate = -1;
	for (int i = start_size; i >= 0; i--)
	{
		if (g_ml_tests[i] != 0)
		{
			expect_rate = g_ml_times[i] / g_ml_tests[i];
			break;
		}
	}
	int max_tests_per_pass = default_tests_per_pass;
	/*
	if (variable_config)
	{
		if (start_size > 130 && start_size < 170)
			max_tests_per_pass = 200;
		else if (expect_rate > 0)
		{
			double target_time = 30.0;
			max_tests_per_pass = target_time / expect_rate;
			if (max_tests_per_pass < 200)
				max_tests_per_pass = 200;
			if (max_tests_per_pass > 1e5)
				max_tests_per_pass = 1e5;
		}
	}
	*/
	g_hit_list.clear();
	sym_potential_graphs.clear();

	JobArray array;
	//vector<bigint> t = cand;
	vector<b128> t;
	if (cand.size() > max_tests_per_pass)
	{
//			random_shuffle(cand.begin(), cand.end());
		t.resize(max_tests_per_pass);
		memcpy(&t[0], &cand[cand.size()-max_tests_per_pass], max_tests_per_pass*sizeof(cand[0]));
		cand.resize(cand.size()-max_tests_per_pass);
	}
		else
	{
		t.swap(cand);
	}
	exclude(t, cover); // in case there are dups in cand
	append(cover, t);

	g_kmax[0] = k0;
	g_kmax[1] = k1;
	array.submit(&t[0], t.size(), d, start_size, start_size + 1, true, 0, 0, general_hit_detect);
	printf(" %d -> %d (%.3f)\n", t.size(), g_hit_list.size(), array.realtime / g_fFreq);
	if (variable_config)
	{
		FILE* f = fopen(PATH "multilevel-perf.txt", "a");
		fprintf(f, "%d\t%d\t%d\t%d\t%.3f\n", d, start_size, t.size(), g_hit_list.size(), array.realtime / g_fFreq);
		fclose(f);
	}
	g_ml_times[start_size] += array.realtime / g_fFreq;
	g_ml_tests[start_size] += t.size();
	vector<b128> v = g_hit_list;
	finish_cyclic_search(k0, k1, g_hit_list, sym_new, min_log_len);
	op_report(d, start_size);

	printf(" NS: %d new, %d cyclic, ", v.size(), sym_new.size());
	vector<b128> new_v;
	relabel(new_v, d, sym_new, (k0==k1));
	append(cover, new_v);
	append(v, new_v);
	printf(" %d relabeled\n", new_v.size());
	vector<b128> new_cand;
	list_neighbors(new_cand, v, d, (k0 == k1));
	append(hits, v);
	rem_dup(hits);
	rem_dup(cover);
	append(sym, sym_new);
	rem_dup(new_cand);
	exclude(new_cand, cover);
	if(cand.size()+new_cand.size() < 1e9) // don't bother if we're pushing limits of avail. memory
	{
		exclude(new_cand, cand);
//		shuffle(cand);
	}
//	shuffle(new_cand);
	append(cand, new_cand);
	shuffle(cand);

	for (int i = 0; i < sym_new.size(); i++)
		max_len = max(max_len, sym_new[i].len + 1);
	printf(" max_len %d\n", max_len);
	int new_target = max(start_size, max_len - 4);
	if (variable_config && new_target > start_size && start_size < stop_len)
	{
		start_size = min(new_target, stop_len);
		d = start_size / 2;
		if (start_size > stop_len - 10)
			d -= (stop_len - start_size) / 2;
		printf(" Moving to %d, %d\n", d, start_size);
//		state_rebuild();
	}
}

#include <time.h>


void unified_search(int k0, int k1, int s, int d, int tests_per_pass, int max_ops, bool self_checks)
{
	printf("unified_search(%d,%d,%d,%d,%d,%d)\n", k0, k1, s, d, tests_per_pass, max_ops);

	if(d >= GRAPH_MAX_LEN*64)
	{
		printf("Error: GRAPH_MAX_LEN=%d; recompile\n", GRAPH_MAX_LEN);
		exit(-1);
	}

	quiet = true;
	//g_sym_graphs.clear();
	vector<bigint> sym_graphs, alt_sym_graphs;
	g_sym_graphs_ref.clear();
	vector<bigint> alt_sym_graphs_ref;

	g_kmax[0] = k0;
	g_kmax[1] = k1;
	readsym(format_sym_name(), sym_graphs, d);
	if(k0!=k1)
	{
		g_kmax[0] = k1;
		g_kmax[1] = k0;
		readsym(format_sym_name(), alt_sym_graphs, d);
		neg(alt_sym_graphs);
		append(sym_graphs, alt_sym_graphs);
		alt_sym_graphs.clear();
		g_kmax[0] = k0;
		g_kmax[1] = k1;
	}
	relabel(g_sym_graphs_ref, sym_graphs, (k0==k1));
	rem_dup(g_sym_graphs_ref);

	multilevel_state state;
	state.variable_config = false;
	state.max_ops = max_ops;
	state.default_tests_per_pass = tests_per_pass;
	state.d = s;
	state.max_len = d;
	state.start_size = d;
	g_kmax[0] = k0;
	g_kmax[1] = k1;
	readfile(format_name(s, d), state.hits, s);
	readfile(format_cover_name(s, d), state.cover, s);

	if (k0 != k1)
	{
		g_kmax[0] = k1;
		g_kmax[1] = k0;
		vector<b128> hits2, cover2, sym2;
		readfile(format_name(s, d), hits2, s);
		readfile(format_cover_name(s, d), cover2, s);

		neg(hits2);
		neg(cover2);
		append(state.hits, hits2);
		append(state.cover, cover2);
	}
	relabel(state.hits, s, sym_graphs, (k0==k1));
	rem_dup(state.hits);
	append(state.cover, state.hits);
	rem_dup(state.cover);

	g_kmax[0] = k0;
	g_kmax[1] = k1;

	printf("Initial state:\n");
	printf("Known signatures: %d\n", state.hits.size());
	printf("Known cyclic graphs: %d\n", sym_graphs.size());
	//printf("0x%llx 0x%llx\n", state.hits2[0].n[0], state.hits2[0].n[1]);
	list_neighbors(state.cand, state.hits, s, (k0==k1));
	if(state.cand.size() > 2e9)
	{
		state.cand.resize(2e9);
		state.cand.reserve(2.1e9);
	}

	exclude(state.cand, state.cover);
	printf("%lld neighbors\n", state.cand.size());
	rem_dup(state.cand);
	printf("Cover size: %lld\n", state.cover.size());
	shuffle(state.cand);
	printf("Candidates left: %lld\n", state.cand.size());
	if(self_checks)
	{
		printf("Self check 1 (recorded hits)\n"); // expect near 100%
		{
			g_hit_list.clear();
			sym_potential_graphs.clear();
			JobArray array;
			vector<b128> t;
			t = state.hits;
			random_shuffle(t.begin(), t.end());
			if(t.size()>100)
				t.resize(100);
			g_max_ops = max_ops;
			g_kmax[0] = k0;
			g_kmax[1] = k1;
			array.submit(&t[0], t.size(), s, d, d + 1, true, 0, 0, general_hit_detect);
			printf(" %d -> %d (%.3f)\n", t.size(), g_hit_list.size(), array.realtime / g_fFreq);
			op_report(s,d);
		}

		printf("Self check 2 (recorded misses)\n"); // - expect near 0%
		{
			vector<b128> t;
			int n = min(10000, state.cover.size());
			t.resize(n);
			memcpy(&t[0], &state.cover[0], n*sizeof(t[0]));
			exclude(t, state.hits);

			g_hit_list.clear();
			sym_potential_graphs.clear();
			JobArray array;
			random_shuffle(t.begin(), t.end());
			if(t.size()>1000)
				t.resize(1000);
			g_max_ops = max_ops;
			g_kmax[0] = k0;
			g_kmax[1] = k1;
			array.submit(&t[0], t.size(), s, d, d + 1, true, 0, 0, general_hit_detect);
			printf(" %d -> %d (%.3f)\n", t.size(), g_hit_list.size(), array.realtime / g_fFreq);
			op_report(s,d);
		}

		printf("Self check 3 (known cyclic)\n"); // Anything other than 100% here means a bug in the code or bad entries in the cyclic file
		{
			vector<bigint> t;
			t = g_sym_graphs_ref;
			shuffle(t);
			if(t.size()>1000)
				t.resize(1000);
			vector<int> success;
			success.resize(nthreads);
	#pragma omp parallel for
			for (int thr = 0; thr < nthreads; thr++)
			{
				graph g;
				for (int pass = thr; pass < t.size(); pass += nthreads)
				{
					g.reset();
					bool ok = test_build(g, t[pass], k0, k1);
					if (ok)
						success[thr]++;
				}
			}
			int s = 0;
			for(int thr=0; thr<nthreads; thr++)
				s+=success[thr];
			printf(" %d -> %d\n", t.size(), s);
			if(t.size() != s)
			{
				printf("ERROR: self check 3 failure!\n");
			}
		}

		printf("Self check 4 (signatures of known cyclic)\n");
		// depends on g_max_ops. Too low -> increase g_max_ops. Too high - potentially inefficient.
		// Expect hits or timeouts. Any misses mean bug in the code or bad entries in the cyclic file.
		{
			g_hit_list.clear();
			sym_potential_graphs.clear();
			JobArray array;
			vector<b128> t;
			t.resize(g_sym_graphs_ref.size());
			for(int i=0; i<t.size(); i++)
			{
				t[i]=g_sym_graphs_ref[i];
				t[i].len=s;
				t[i].unset_from(s);
			}
			random_shuffle(t.begin(), t.end());
			if(t.size()>1000)
				t.resize(1000);
			g_max_ops = max_ops;
			g_kmax[0] = k0;
			g_kmax[1] = k1;
			array.submit(&t[0], t.size(), s, d, d + 1, true, 0, 0, general_hit_detect);
			printf(" %d -> %d (%.3f)\n", t.size(), g_hit_list.size(), array.realtime / g_fFreq);
			if(op_counts[0].size()>0)
			{
				printf("ERROR: self check 4 failure!\n");
			}
			op_report(s,d);
	#if 0
			if (t.size() != g_hit_list.size())
			{
				exclude(t, g_hit_list);
				for(int i=0; i<t.size(); i++)
				{
					for(int j=0; j<g_sym_graphs_ref.size(); j++)
					{
						bigint x;
						x=g_sym_graphs_ref[j];
						x.len=s;
						x.unset_from(s);
						b128 xx;
						xx = x;
						if(bit_equal(xx, t[i]))
						{
							printf("%d\t0x%llx\n", x.len+1, x.n[0]);
							break;
						}
					}
				}
				g_hit_list.clear();
				g_max_ops = max_ops*3;
				array.submit(&t[0], t.size(), s, d, d + 1, true, 0, 0, general_hit_detect);

				printf(" at max_ops=%d: %d -> %d (%.3f)\n", g_max_ops, t.size(), g_hit_list.size(), array.realtime / g_fFreq);
				op_report(s, d);
				g_max_ops = max_ops;
			}
	#endif
		}
	}

	g_kmax[0] = k0;
	g_kmax[1] = k1;
	writefile(format_name(s, d), state.hits, s);
	writefile(format_cover_name(s, d), state.cover, s);
	max_errors_detect_only=16;
	while (state.cand.size() != 0)
	{
		
		time_t t = time(0);
		struct tm* lt = localtime(&t);
		if(lt->tm_hour>=10 && lt->tm_hour<21 && lt->tm_wday>=1 && lt->tm_wday<=5)
		{
			Sleep(100000);
			continue;
		}
		if(quiet)
			printf(".");

//		printf("Current time is: %02d:%02d:%02d\n", lt->tm_hour, lt->tm_min, lt->tm_sec);

		vector<b128> coverPrev = state.cover;
		vector<b128> hitsPrev = state.hits;
		state.single_pass(k0, k1, d, d);
		vector<b128> coverNew = state.cover;
		vector<b128> hitsNew = state.hits;
		exclude(coverNew, coverPrev);
		exclude(hitsNew, hitsPrev);
		writefile(format_name(s, d), hitsNew, s, true);
		writefile(format_cover_name(s, d), coverNew, s, true);
	}
}


void signature_ladder_search(int k0, int k1, int start_size, int mc_test_count, int tests_per_pass, int min_log_len, int stop_len)
{
	for (int i = 0; i <= GRAPH_MAX_LEN * 64; i++)
	{
		g_ml_tests[i] = 0;
		g_ml_times[i] = 0;
	}

	bool diagonal = (k0 == k1);
	vector<bigint> sym_graphs;
	g_sym_graphs_ref.clear();

	g_kmax[0] = k0;
	g_kmax[1] = k1;
	readsym(format_sym_name(), sym_graphs, min_log_len);
	relabel(g_sym_graphs_ref, sym_graphs, k0==k1);
	rem_dup(g_sym_graphs_ref);

	vector<b128> hits;
	vector<bigint> sym;
	g_max_ops = 1e5;
	max_errors_detect_only = 12;
	int d = start_size / 2 - 5;
	mc_memory(d, k0, k1, start_size, mc_test_count, hits, sym, min_log_len);
	printf("%d %d\n", hits.size(), sym.size());
	if (hits.size() == 0)
		return;
	printf("%d symmetric\n", sym.size());
	if (sym.size() == 0)
		return;
	rm_relabelings(sym, diagonal);
	quiet = true;

	vector<vector<bigint> > vv;
	clusterize(vv, sym, true);

	vector<multilevel_state> vm;
	vm.resize(vv.size());
	for (int i = 0; i < vm.size(); i++)
	{
		vm[i].d = d;
		vm[i].start_size = start_size;
		vm[i].sym = vv[i];
		vm[i].variable_config = true;
		vm[i].default_tests_per_pass = tests_per_pass;
		
		relabel(vm[i].hits, d, vm[i].sym, k0==k1);
		rem_dup(vm[i].hits);
		vm[i].cover = vm[i].hits;
		list_neighbors(vm[i].cand, vm[i].hits, d, k0==k1);
		vm[i].max_len = 0;
		for (int j = 0; j < vv[i].size(); j++)
			vm[i].max_len = max(vm[i].max_len, vv[i][j].len + 1);
	}

	int nStaticCount = 0;
	int pass = 0;
	while (true)
	{
		int rem_cand = 0;
		for (int nCl = 0; nCl < vv.size(); nCl++)
		{
			if (vm[nCl].cand.size() == 0)
				continue;
			printf("=== Pass %d, cluster %d ===\n", pass, nCl);
			vm[nCl].single_pass(k0, k1, min_log_len, stop_len);
			rem_cand += vm[nCl].cand.size();
		}
		if (rem_cand == 0)
			break;

		FILE* f = fopen(PATH "multilevel2_trace.txt", "a");
		fprintf(f, "%d\t", pass);
		for (int i = 0; i < vv.size(); i++)
		{
			if (vm[i].cand.size() == 0)
				fprintf(f, "-\t");
			else
				fprintf(f, "%d\t", vm[i].max_len);
		}
		fprintf(f, "\n");
		fclose(f);
		int nMinStartSize = vm[0].start_size;
		for (int nCl = 1; nCl < vv.size(); nCl++)
			nMinStartSize = min(nMinStartSize, vm[nCl].start_size);

		if (nMinStartSize <= start_size)
			nStaticCount++;
		else
		{
			printf("New start_size %d\n", nMinStartSize);
			start_size = nMinStartSize;
			nStaticCount = 0;
		}
		if (nStaticCount >= 3 && start_size<stop_len)
		{
			start_size++;
			printf("New start_size %d\n", start_size);
			nStaticCount = 0;
		}
		refilter(g_sym_graphs_ref, start_size);
		for (int i = 0; i < vv.size(); i++)
		{
			if (vm[i].start_size < start_size)
			{
				vm[i].start_size = start_size;
				vm[i].d = start_size / 2;
				if (start_size > stop_len - 10)
					vm[i].d -= (stop_len - start_size) / 2;
				vm[i].state_rebuild();
			}
		}

		vector<bigint> ts;
		for (int i = 0; i < vv.size(); i++)
		{
			int n = vm[i].sym.size();
			rm_relabelings(vm[i].sym, diagonal);
			int n_ = vm[i].sym.size();
			printf("%d: %d symmetric, %d duplicates removed\n", i, n_, n - n_);
			append(ts, vm[i].sym);
		}
		int n = ts.size() ;
		rm_relabelings(ts, diagonal);
		int n_ = ts.size();
		printf("total: %d symmetric, %d duplicates removed\n", n_, n - n_);
		if (n - n_ != 0)
			printf("WARNING: duplicates in the sym lists\n");
		pass++;
	}
}	


void load_known_sym_graphs(int k0, int k1, int d)
{
	vector<bigint> sym_graphs;
	g_sym_graphs_ref.clear();

	g_kmax[0] = k0;
	g_kmax[1] = k1;
	readsym(format_sym_name(), sym_graphs, d);
	relabel(g_sym_graphs_ref, sym_graphs, (k0 == k1));
	rem_dup(g_sym_graphs_ref);
}

void horizontal_search(int k0, int k1, int d, int dmax, int max_tests)
{
	g_kmax[0] = k0;
	g_kmax[1] = k1;
	g_symfn = format_sym_name();
	vector<bigint> out;
	readsym(format_sym_name(), out, d);	

	g_sym_graphs_ref.clear();
	g_kmax[0] = k0;
	g_kmax[1] = k1;
	relabel(g_sym_graphs_ref, out, (k0 == k1));
	rem_dup(g_sym_graphs_ref);

	if(dmax!=INT_MAX)
	{
		vector<bigint> out2;
		for(int i=0; i<out.size(); i++)
			if(out[i].len+1<=dmax)
				out2.push_back(out[i]);
		out.swap(out2);
	}
	
	shuffle(out);
	if(max_tests!=INT_MAX)
	{
		if(out.size()>max_tests)
			out.resize(max_tests);
	}	
	while(true)
	{
		sym_potential_graphs.clear();
		for(int i=0; i<out.size(); i++)
			for(int j=0; j*2<=out[i].len; j++)
			{
				for(int k=j; k*2<=out[i].len; k++)
				{
					bigint x = out[i];
					x.flip(j);
					if(x.len-j-1!=j)
						x.flip(x.len-j-1);
					if(k!=j)
					{
						x.flip(k);
						if(x.len-k-1!=k)
							x.flip(x.len-k-1);
					}
					if(k0==k1 && x.bit(0))
					{
						x.neg();
						x.unset_from(x.len);
					}
					sym_potential_graphs.push_back(x);
				}
			}
		rem_dup(sym_potential_graphs);
		//rm_relabelings(sym_potential_graphs, (k0==k1));
		printf("%d potential symmetric\n", sym_potential_graphs.size());
		exclude(sym_potential_graphs, g_sym_graphs_ref);
		printf("%d not known\n", sym_potential_graphs.size());
		rem_dup(sym_potential_graphs);
		printf("%d unique\n", sym_potential_graphs.size());
		uint64_t t1, t2;
		t1 = __rdtsc();
	#pragma omp parallel for
		for (int thr = 0; thr<nthreads; thr++)
		{
			graph g;
			for (int i = thr; i<sym_potential_graphs.size(); i += nthreads)
			{
				g.reset();
				bigint y = sym_potential_graphs[i];
				bigint2 mask;
				mask.first = y;
				mask.second = y;
				mask.first.neg();
				mask.first.unset_from(y.len);
				bool ok = build_graph(g, mask, g_kmax);
				if (!ok)
					sym_potential_graphs[i].set_len(0);
			}
		}
		t2 = __rdtsc();

		refilter(sym_potential_graphs);
		printf("%d survive build\n", sym_potential_graphs.size());

		rm_relabelings(sym_potential_graphs, (k0==k1));
		printf("%d non isomorphic\n", sym_potential_graphs.size());

		out.swap(sym_potential_graphs);
		if(out.size()==0)
			break;

		for (int64_t i = 0; i<out.size(); i++)
		{
			vector<bigint> vnew;
			relabel(vnew, out[i], (k0==k1));
			append(g_sym_graphs_ref, vnew);
		}

		FILE* f = fopen(g_symfn.c_str(), "a");
		for (int64_t i = 0; i<out.size(); i++)
		{
			bigint y = out[i];
			int j = y.len;
			if (y.len + 1 >= d)
			{
				//g_sym_graphs_ref.push_back(y);
	#if GRAPH_MAX_LEN==4
				fprintf(f, "%d\t0x%llx\t0x%llx\t0x%llx\t0x%llx\t0x0\t0x0\n",
					j + 1, y.n[0], y.n[1], y.n[2], y.n[3]);
	#elif GRAPH_MAX_LEN==6
				fprintf(f, "%d\t0x%llx\t0x%llx\t0x%llx\t0x%llx\t0x%llx\t0x%llx\n",
					j + 1, y.n[0], y.n[1], y.n[2], y.n[3], y.n[4], y.n[5]);
	#else
				fprintf(f, "%d\t0x%llx\t0x%llx\t0x0\t0x0\t0x0\t0x0\n",
					j + 1, y.n[0], y.n[1]);
	#endif
			}
		}
		fclose(f);
	}

	rem_dup(g_sym_graphs_ref);
}

int combin(int x, int y)
{
	int n=1;
	for(int i=0; i<x; i++)
	{
		n *= y-i;
		n /= (i+1);
	}
	return n;
}

bigint reflect(bigint b, int new_len)
{
	int dx = new_len - b.len;
	int new_half_len = (new_len + 1) / 2;
	bigint y0, y0inv;
	y0.clear();
	y0inv.clear();
	y0 = b;
	if (dx>0)
		y0inv = b << dx;
	else if (dx<0)
		y0inv = b >> -dx;
	else
		y0inv = b;
	y0.unset_from(new_half_len);
	y0inv.unset_to(new_half_len);
	y0 = y0 | y0inv;
	y0.set_len(new_len);
	return y0;
}

static inline void construct_cand_ext(bigint& cand, bigint y, const bigint* shifted_ref_masks) 
{
	for(int p=0; p<GRAPH_MAX_LEN; p++)
	{
		if(y.n[p]==0)
			continue;
		do
		{
			unsigned long pos;
			_BitScanForward64(&pos, y.n[p]);
			cand &= shifted_ref_masks[pos+p*64+1];
			y.n[p] &= ~(One<<pos);
		}
		while(y.n[p]!=0);

		if(cand.zero())
		{
			p=GRAPH_MAX_LEN;
			break;
		}
	}
}


// Enables additional tracing and assertion checking inside cyclic search methods. Keep off during production runs.
bool g_accuracy_check = false;



// Generate a list of 1-incomplete cliques to be used to make pre-vet masks during cyclic search.
void do_minus_ones(graph* g0, bigint b, vector<int>& vPos, vector<pair<bigint,int> >& vMinusOnes)
{
	int half_len = (b.len + 1) / 2;
	//vector<pair<bigint,int> > vMinusOnes;
	vPos.clear();
	vMinusOnes.resize(0);
	vMinusOnes.reserve((g0->cliques[0].size()+g0->cliques[1].size())/5);
	//vPotential.clear();

	bigint allowed;
	allowed.clear();
	allowed.neg();
	allowed.unset_from(half_len);
	int nMinusOne=0;
		
	bigint shifted_ref_masks[2][GRAPH_MAX_LEN*64+1];
	for(int i=1; i<=b.len; i++)
	{
		bigint temp;
		temp.clear();
		invert(temp, g0->mask0);
		if(g0->mask0.len-i>0)
			shifted_ref_masks[0][i] = temp >> (g0->mask0.len-i);
		else if(g0->mask0.len-i<0)
			shifted_ref_masks[0][i] = temp << (i-g0->mask0.len);
		else
			shifted_ref_masks[0][i] = temp ;

		shifted_ref_masks[0][i] |= g0->mask0 << (i+1);
		temp.clear();
		invert(temp, g0->mask);
		if(g0->mask.len-i>0)
			shifted_ref_masks[1][i] = temp >> (g0->mask.len-i);
		else if(g0->mask.len-i<0)
			shifted_ref_masks[1][i] = temp << (i-g0->mask.len);
		else
			shifted_ref_masks[1][i] = temp ;
			
		shifted_ref_masks[1][i] |= g0->mask << (i+1);
	}
		
	for(int color=0; color<2; color++)
	{
		bigint base_cand;
		base_cand.clear();
		for (int k = 0; k<GRAPH_MAX_LEN; k++)
			base_cand.n[k] = -1;
		base_cand.set_len(b.len+1);
		base_cand.unset_from(b.len+1);
		base_cand &= (color ? g0->mask0 : g0->mask) << 1;
		int counts[GRAPH_MAX_LEN*64+1];
		memset(counts, 0, sizeof(counts));
		//char c2[GRAPH_MAX_LEN*64*GRAPH_MAX_LEN*64];
		//memset(c2, 0, sizeof(c2));

		int threshold = g_kmax[color]/2;
		int posArray = vMinusOnes.size();
		for(int i=0; i<g0->cliques[color].size(); i++)
		{
			const bigint& x = g0->cliques[color][i];
//				if(popcnt(x)!=x.popcnt)
//					printf("ERROR\n");
			if(x.popcnt!=g_kmax[color])
				continue;
			bigint cand=andnot(base_cand, x<<1);
			if(cand.zero())
				continue;
			construct_cand_ext(cand, x, shifted_ref_masks[color]);
			//base_cand = andnot(base_cand, cand);
			cand >>= 1;
			while(!cand.zero())
			{
				int j=cand.trailing_bit();
#ifdef DIRECT_PRE_VET_MODE
				int n=0;
				// Cheap emulation of best known logic for this operation.
				// Ideally, I want 'n' to be the number of links with indexes above half_len in this clique.
				// However, it's expensive to calculate. 
				// Instead, I do just links from root node to all nodes ('y'), and links from 'j' to all nodes with indexes below 'j' ('zz').
				// Consequently, the threshold is lower (g_kmax[color]/2 here, vs. g_kmax[color]-1 used in indirect mode)
				if(j>=half_len)
				{
					bigint y = x;
					y.set(j);
					bigint z = x;
					z.set_len(j);
					z.unset_from(j);
					bigint zz;
					invert(zz, z);
					y|=zz;

					y.unset_to(half_len);
					n = popcnt(y);
				}
				if(counts[j]==0 || (n<=threshold && j>=half_len))
#else
				//if(counts[j]<=4 || j>=x.trailing_bit())
#endif
				{
					bigint y = x;
					y.len = b.len;
					y.popcnt = color;
					vMinusOnes.push_back(make_pair(y,j));
				}
				allowed.unset(min(j,b.len-1-j));
				counts[j]++;
				//if(g0->cliques[color].size()>100)
				{
					if(j<half_len)
						base_cand.unset(j+1);
	#ifndef DIRECT_PRE_VET_MODE
					// Has no effect in most cases.
					// In cases where the number of cliques is exploding, improves performance by keeping vMinusOnes reasonably sized
					// ((e.g. see the R(7,9) call in pre_vet_accuracy_test(): this line cuts the list size from 230K to 7K elements)
					if(counts[j]*2>j)
						base_cand.unset(j+1);
	#else
					if(counts[j]>j)
						base_cand.unset(j+1);
	#endif
				}
				cand.unset(j);
			}
			if(base_cand.zero())
				break;
		}
#if 0
		if(vMinusOnes.size()-posArray < 100)
		{
			for(int i=0; i<g0->cliques[color].size(); i++)
			{
				const bigint& x = g0->cliques[color][i];
				if(x.popcnt!=g_kmax[color])
					continue;
				int v[MAX_DEPTH], vp=0;
				bigint y=x;
				while(!y.zero())
				{
					int p = y.trailing_bit();
					v[vp++]=p;
					y.unset(p);
				}
				for(int j=0; j<b.len; j++)
				{
					if(b.bit(j)==color)
						continue;
					int n=0;
					for(int k=0; k<vp; k++)
						if(b.bit(abs(j-v[k])-1)!=color)
							n++;
					//printf("%d,%d	%d\n", i, j, n);
					if(n==1)
					{
						bigint y = x;
						y.len = b.len;
						y.popcnt = color;
						vMinusOnes.push_back(make_pair(y,j));
							if(vMinusOnes.size()-posArray >= 90)
							break;
					}
				}
				if(vMinusOnes.size()-posArray >= 90)
					break;
			}
		}
#endif
	}
	int nAllowed = popcnt(allowed);
		
	vPos.resize(nAllowed);
	bigint x = allowed;
	int i=0;
	while(!x.zero())
	{
		vPos[i]=x.leading_bit();
		x.unset(vPos[i]);
		i++;
	}
}

vector<int> g_occupied;

struct CyclicSearchObject
{
public:
	static int m_free_bits;
public:
	mutable graph* g0;
	bigint* v;
	int k0, k1;
	vector<bigint>& vOut;
	vector<bool>& vOrig;
	int fwd_search_range;
	int top_search_range;
	bool diag;
	bool fwd_only;
	int min_length;
	bool random_mode;
	int nflips;
	int tests_per_call;
	int graph_pos;

public:
	bool off_by_2;
	bool do_flip_vertical;
	int fwd_flip_range;
public:
	void init()
	{
		diag=(k0==k1);
		omp_set_lock(&op_count_lock);
		//graph_pos=g_pool_usage;
		g0=0;
		for(int i=0; i<g_occupied.size(); i++)
		{
			if(g_occupied[i]==0)
			{
				g0=g_graph_pool+i;
				graph_pos=i;
				g_occupied[i]=1;
				break;
			}
		}
		if(g0==0)
		{
			g0 = g_graph_pool+g_occupied.size();
			graph_pos=g_occupied.size();
			g_occupied.push_back(1);
		}
		
	
		omp_unset_lock(&op_count_lock);
	}

	const pair<int,int>* m_pNewIdx;
	const vector<int>* m_pvPos;
	const vector<pair<bigint,int> >* m_pvPotential;
	void set_graphs(const pair<int,int>* pNewIdx, const vector<int>* pvPos, const vector<pair<bigint,int> >* pvPotential)
	{
		m_pNewIdx=pNewIdx;
		m_pvPos = pvPos;
		m_pvPotential = pvPotential;
	}

	mutable vector<pair<bigint,int> > m_vPotTemp;
	void do_minus_ones_fast(bigint b, int index) const
	{
		int half_len = (b.len+1)/2;
		m_vPos.clear();
		m_vPotential[0].clear();
		m_vPotential[1].clear();
		m_vLowHalf[0].clear();
		m_vLowHalf[1].clear();

		int nColors[2]={0,0};

#ifndef DIRECT_PRE_VET_MODE
		if(m_pNewIdx == 0)
#endif
		{
			test_build(*g0, b, k0, k1);
			m_vPotTemp.clear();
			::do_minus_ones(g0, b, m_vPos, m_vPotTemp);

			for(int i=0; i<m_vPotTemp.size(); i++)
			{
				int miss_pos=m_vPotTemp[i].second;
				bigint x = m_vPotTemp[i].first;
				int color = x.popcnt;
				nColors[color]++;
//					x.set(miss_pos);

				bigint y = x;
				y.set_len(b.len);
				y.set(miss_pos);
				int v[MAX_DEPTH];
				int vp=0;
				while(!x.zero())
				{
					int p = x.trailing_bit();
					y.set(abs(miss_pos-p)-1);
					x.unset(p);
					v[vp++]=p;
				}
			
				for(int j=0; j<vp; j++)
					for(int k=j+1; k<vp; k++)
						y.set(abs(v[j]-v[k])-1);
			
				bigint bCond = y;
				bCond.unset_to(half_len);
				int n = popcnt(bCond);
				if(n*2<=popcnt(y))
				{
					y.set_len(y.trailing_bit()+1);
					m_vPotential[color].push_back(y);
				}
				if(bCond.zero())
					m_vLowHalf[color].push_back(y);
			}
			if(g_accuracy_check && index==0)
				printf("Received %d + %d cliques\n", nColors[0], nColors[1]);

			return;
		}

		pair<int,int> pos = m_pNewIdx[index];
		int mul = -1;//pos.second;
		for(int i=1; i<=b.len; i++)
		{
			if (gcd(i, b.len + 1) != 1)
				continue;
			if (((i*pos.second) % (b.len + 1)) == 1)
			{
				mul = i;
				break;
			}
		}
		if(mul<0)
		{
			omp_set_lock(&op_count_lock);
			printf("Cant find mul\n");
			printf("pos.second %d, b.len+1 %d\n", pos.second, b.len+1);
			for(int i=1; i<=b.len; i++)
			{
				if (gcd(i, b.len + 1) != 1)
					continue;
				printf("%d => %d\n", i, ((i*pos.second) % (b.len + 1)) );
			}

			exit(-1);
		}

		const vector<pair<bigint, int> >& vPotentialRef = m_pvPotential[pos.first];

		m_vPotential[0].clear();
		m_vPotential[1].clear();
		m_vLowHalf[0].clear();
		m_vLowHalf[1].clear();
		
		const vector<int>& vPosRef = m_pvPos[pos.first];
		m_vPos.resize(vPosRef.size());

		for(int i=0; i<vPosRef.size(); i++)
		{
			m_vPos[i] = (((vPosRef[i] + 1)*mul) % (b.len + 1)) - 1;
			m_vPos[i] = min(m_vPos[i], b.len-m_vPos[i]-1);
		}
		
		vector<bigint> vShortest[4];
		vector<int> vTr[4];
		for(int j=0; j<4; j++)
		{
			vShortest[j].resize(b.len);
			vTr[j].resize(b.len);
			for(int i=0; i<b.len; i++)
			{
				vShortest[j][i].len=0;
				vTr[j][i]=b.len+1;
			}		
		}

		int mapping[GRAPH_MAX_LEN*64];
		for(int i=0; i<b.len; i++)
			mapping[i]= (((i+1)*mul) % (b.len + 1)) - 1;
		int counts[2]={0,0};
		for(int i=0; i<vPotentialRef.size(); i++)
		{
			const bigint& x__=vPotentialRef[i].first;
			int color = x__.popcnt;
			counts[color]++;
		}
		for(int i=0; i<vPotentialRef.size(); i++)
		{
			bigint x__=vPotentialRef[i].first;
			int miss_pos = mapping[vPotentialRef[i].second];//(((vPotentialRef[i].second + 1)*mul) % (b.len + 1)) - 1;
			//int color = x__.popcnt;
			//if(k0==k1 && 
			int color = 1-b.bit(miss_pos);
			const bool bAcceptAll = false;//(counts[color]<100);

			nColors[color]++;
			bigint x, x_inv;
			x.clear();
			x.set_len(b.len);
			x_inv.clear();
			x_inv.set_len(b.len);
			short v2[MAX_DEPTH];
			int vp2=0;
			int n_low=0;
			int minpos=b.len, maxpos=0;
			while(!x__.zero())
			{
				int p = x__.trailing_bit();
				int q=mapping[p];
				x.set(q);
				x_inv.set(b.len-1-q);
				x__.unset(p);
				v2[vp2++]=q;
				if(q<half_len)
					n_low++;
				maxpos=max(maxpos,q);
				minpos=min(minpos,q);
			}
			int threshold = g_kmax[color]-1;

			v2[vp2++]=miss_pos;
			int mask=3;
			if(!bAcceptAll)
			{
				if((maxpos>=vTr[color][miss_pos] && maxpos>=vTr[color+2][miss_pos])
					&& (miss_pos<half_len))// || (n_low>1 && vp2-1-n_low>1 && maxpos-minpos>=half_len)))
					mask &= ~1;
				if((b.len-1-minpos>=vTr[color][b.len-1-miss_pos] && b.len-1-minpos>=vTr[color+2][b.len-1-miss_pos])
					&& (miss_pos>=half_len))// || (n_low>1 && vp2-1-n_low>1 && maxpos-minpos>=half_len)))
					mask &= ~2;
			}
//			if(mask == 0)
//				continue;

			x.set(miss_pos);
			x_inv.set(b.len-1-miss_pos);

			bigint y1;
			y1.clear();
			y1.set_len(b.len);
#if 0
			if(b.len<255)
			{
				__m128i x0, x1, acc0, acc1;
				x0 = *(__m128i*) x.n;
				x1 = *(__m128i*) (x.n+2);
				acc0=_mm_setzero_si128();
				acc1=_mm_setzero_si128();
				for(int j=0; j<vp2; j++)
				{
					if(v2[j]+1<128)
					{
						__m128i x0_ = shift_right(x0, v2[j]+1);
						__m128i x1_ = shift_right(x1, v2[j]+1);
						__m128i x1__ = shift_left(x1, 128-(v2[j]+1));
						acc0 = _mm_or_si128(acc0, x0_);
						acc0 = _mm_or_si128(acc0, x1__);
						acc1 = _mm_or_si128(acc1, x1_);
					}
					else if(v2[j]+1==128)
					{
						acc0 = _mm_or_si128(acc0, x1);
					}
					else
					{
						__m128i x1_ = shift_right(x1, v2[j]+1-128);
						acc0 = _mm_or_si128(acc0, x1_);
					}
				}
				*(__m128i*)(y1.n) = acc0;
				*(__m128i*)(y1.n+2) = acc1;

				if(g_accuracy_check)
				{
					bigint y2;
					y2.clear();
					y2.set_len(b.len);
					for(int j=0; j<vp2; j++)
						for(int k=j+1; k<vp2; k++)
							y2.set(abs(v2[k]-v2[j])-1);
					if(!bit_equal(y2, y1))
					{
						omp_set_lock(&op_count_lock);
						printf("Bit mismatch\n");
						x.printbits();
						y1.printbits();
						y2.printbits();
						exit(-1);
					}
				}
				/*
				short buf[16];
				__m128i mv20, mv21;
				mv20=*(__m128i*)(v2+0);
				mv21=*(__m128i*)(v2+8);
				__m128i mm_1 = _mm_set1_epi16(1);
				for(int j=0; j<vp2 && j<7; j++)
				{
					__m128i s = _mm_set1_epi16(v2[j]);
					__m128i a0 = _mm_sub_epi16(mv20, s);
					__m128i a1 = _mm_sub_epi16(mv21, s);
					a0 = _mm_abs_epi16(a0);
					a1 = _mm_abs_epi16(a1);
					a0 = _mm_sub_epi16(a0, mm_1);
					a1 = _mm_sub_epi16(a1, mm_1);
					*(__m128i*)buf = a0;
					*(__m128i*)(buf+8) = a1;
					for(int k=j+1; k<vp2; k++)
						y1.set(buf[k]);
				}
				for(int j=7; j<vp2; j++)
				{
					__m128i s = _mm_set1_epi16(v2[j]);
					__m128i a1 = _mm_sub_epi16(mv21, s);
					a1 = _mm_abs_epi16(a1);
					a1 = _mm_sub_epi16(a1, mm_1);
					*(__m128i*)(buf+8) = a1;
					for(int k=j+1; k<vp2; k++)
						y1.set(buf[k]);
				}
					*/

			}
			else
#endif
			{
				/****
				This is typically THE bottleneck in cyclic search in the indirect mode. (This loop, the previous while() loop,
				and CyclicSearchObject::pre_vet() account for 50-80% of all CPU time.) 
				It might seem like a good idea to move this loop out of here and into ::do_minus_ones(),
				and feed crosslink indexes into the same relabel code above that works on node indexes,
				but it won't work (and the reason why not is far from obvious)
				*****/

				for(int j=0; j<vp2-1; j++)
					for(int k=j+1; k<vp2; k++)
						y1.set(abs(v2[k]-v2[j])-1);
			}
//			for(int j=0; j<vp2; j++)
//				y1 |= x_ >> (v2[j]+1);

			for(int dir=0; dir<2; dir++)
			{
				if(!bAcceptAll && !(mask && (1<<dir)))
					continue;
				bigint y = ((dir==0) ? x : x_inv) | y1;
				bigint bCond = y;
				bCond.unset_to(half_len);
				int n = popcnt(bCond);

				if(!bAcceptAll && n>=4 && n*2>=popcnt(y))
					continue;
				if(bCond.zero())
					m_vLowHalf[color].push_back(y);

				int m = dir ? b.len-1-miss_pos : miss_pos;
				if(bAcceptAll || (m>=half_len && n<=threshold))
				{
					y.set_len(y.trailing_bit()+1);
					m_vPotential[color].push_back(y);
					continue;
				}
				
				int tr = n;
				if(vShortest[color][m].len==0)
				{
					vShortest[color][m]=y;
					vTr[color][m]=tr;
					continue;
				}
				if(vShortest[color+2][m].len==0)
				{
					vShortest[color+2][m]=y;
					vTr[color+2][m]=tr;
					continue;
				}

				if(vTr[color][m]>tr)
				{
					vTr[color][m]=tr;
					vShortest[color][m]=y;
					continue;
				}

				if(vTr[color+2][m]>tr)
				{
					vTr[color+2][m]=tr;
					vShortest[color+2][m]=y;
					continue;
				}
			}
		}

		if(g_accuracy_check && index==0)
			printf("Received %d + %d cliques\n", nColors[0], nColors[1]);


		int nBest = (g_kmax[0]>=8 || g_kmax[1]>=8) ? 2 : 1;
		for(int i=0; i<b.len; i++)
			for(int j=0; j<2*nBest; j++)
				if(vShortest[j][i].len>0)
				{
					bigint y=vShortest[j][i];
					int color=j&1;
					y.set_len(y.trailing_bit()+1);
					m_vPotential[color].push_back(y);
				}

	}

	~CyclicSearchObject()
	{
		omp_set_lock(&op_count_lock);
		g_occupied[graph_pos]=0;
		omp_unset_lock(&op_count_lock);
	}

	// constructor for random horizontal flipping
	CyclicSearchObject(bigint* _v, int _k0, int _k1, int _nflips, int _tests_per_call, vector<bigint>& _vOut)
		: v(_v), k0(_k0), k1(_k1), nflips(_nflips), tests_per_call(_tests_per_call), vOut(_vOut),
		vOrig(vector<bool>()),
		random_mode(true), fwd_search_range(0), top_search_range(0), min_length(0),
		do_flip_vertical(false), fwd_flip_range(0), off_by_2(true)
		, m_pNewIdx(0), m_pvPos(0), m_pvPotential(0)
	{
		init();
	}

	// // constructor for fully-enumerated horizontal and vertical flipping 
	CyclicSearchObject(bigint* _v, vector<bool>& orig, int _k0, int _k1, int _fwd_search_range, int _top_search_range, bool _fwd_only, 
		int _min_length,
		vector<bigint>& _vOut) : v(_v), k0(_k0), k1(_k1), vOut(_vOut), fwd_search_range(_fwd_search_range), top_search_range(_top_search_range), fwd_only(_fwd_only), min_length(_min_length), 
		vOrig(orig),
		random_mode(false),
		nflips(0), tests_per_call(0),
		do_flip_vertical(true), fwd_flip_range(_fwd_search_range), off_by_2(true)
		, m_pNewIdx(0), m_pvPos(0), m_pvPotential(0)
	{
		init();
	}

	// constructor for random forward flipping
	CyclicSearchObject(bigint* _v, vector<bool>& orig, int _k0, int _k1, int _fwd_search_range, int _min_length, int _nflips, int _tests_per_call,
		vector<bigint>& _vOut) : v(_v), k0(_k0), k1(_k1), vOut(_vOut), fwd_search_range(_fwd_search_range), top_search_range(0), fwd_only(true), min_length(_min_length), 
		vOrig(orig),
		random_mode(true),
		nflips(_nflips),
		tests_per_call(_tests_per_call),
		do_flip_vertical(true), fwd_flip_range(_fwd_search_range), off_by_2(true)
		, m_pNewIdx(0), m_pvPos(0), m_pvPotential(0)
	{
		init();
	}
	
	// copy constructor (demanded by TBB)
	CyclicSearchObject(const CyclicSearchObject& c) : vOut(c.vOut), vOrig(c.vOrig)
	{
		random_mode=c.random_mode;
		v = c.v;
		k0=c.k0;
		k1=c.k1;
		fwd_search_range=c.fwd_search_range;
		top_search_range=c.top_search_range;
		diag=c.diag;
		fwd_only=c.fwd_only;
		min_length=c.min_length;
		nflips=c.nflips;
		tests_per_call=c.tests_per_call;
		do_flip_vertical=c.do_flip_vertical;
		fwd_flip_range=c.fwd_flip_range;
		off_by_2=c.off_by_2;
		m_pNewIdx=c.m_pNewIdx;
		m_pvPos=c.m_pvPos;
		m_pvPotential=c.m_pvPotential;
		init();
	}

    void operator()( const blocked_range<size_t>& r ) const 
	{
        for( size_t i=r.begin(); i!=r.end(); ++i ) 
		{
			bigint b = v[i];
			bigint b_inv;
			invert(b_inv, b);
			if(!bit_equal(b_inv, b))
			{
				printf("Error: graph not symmetric\n");
				continue;
			}
			do_minus_ones_fast(b, i);
			if(!random_mode)
			{
				do_bit_flip(b, vOrig[i], i);
				if(vOrig[i]==true)
					do_top_replacement(b, i);
			}
			else if(fwd_search_range>0)
			{
				do_random_flip(b, i);
			}
			else
			{
				do_random_flip_sideways(b);
			}
		}
	}

//#define ERROR_CHECK
	bool pre_vet(bigint y, int idx=-1) const 
	{
		bool bad=false;
		size_t s = m_vPotential[0].size();
		for(int j=0; j<s; j++)
		{
			if((m_vPotential[0][j] & y).zero())
			{
				if(m_vPotential[0][j].len<=y.len)
				{
#ifdef ERROR_CHECK
					bool ok = test_build(*g0, y, k0, k1);
					if(ok)
					{
						omp_set_lock(&op_count_lock);
						printf("ERROR: skipping good build idx=%d due to check clique 0/%d\n", idx, j);
						y.printbits();
						m_vPotential[0][j].printbits();
						exit(-1);
					}
#endif
					return true;
				}
			}
		}
		s  = m_vPotential[1].size();
		for(int j=0; j<s; j++)
		{
			if(andnot(m_vPotential[1][j], y).zero())
			{
				if(m_vPotential[1][j].len<=y.len)
				{
	#ifdef ERROR_CHECK
					bool ok = test_build(*g0, y, k0, k1);
					if(ok)
					{
						omp_set_lock(&op_count_lock);
						printf("ERROR: skipping good build idx=%d due to check clique 1/%d\n", idx, j);
						y.printbits();
						m_vPotential[1][j].printbits();
						exit(-1);
					}
	#endif
					return true;
				}
			}
		}
		return false;
	}

	bool pre_vet_low_half(bigint y, int idx=-1) const 
	{
		bool bad=false;
		for(int j=0; j<m_vLowHalf[0].size(); j++)
		{
			if((m_vLowHalf[0][j] & y).zero())
				return true;
		}
		for(int j=0; j<m_vLowHalf[1].size(); j++)
		{
			if(andnot(m_vLowHalf[1][j], y).zero())
				return true;
		}
		return false;
	}	

	void do_top_replacement(bigint b, int idx) const
	{
		if(fwd_only)
			return;

		int half_len = (b.len + 1) / 2;
		int base = half_len - top_search_range;
		int hit_cnt=0;
		
		bigint y0root = b;
		y0root.unset_from(base);
		y0root.set_len(base);

		bigint* temp = new bigint[(1 << top_search_range)];

		for (int n = 0; n<(1 << top_search_range); n++)
		{
			bigint y = b;
			for (int k = 0; k<top_search_range; k++)
			{
				y.unset(base + k);
				y.unset(y.len - 1 - (base + k));
				if ((n >> k) & 1)
				{
					y.set(base + k);
					y.set(y.len - 1 - (base + k));
				}
			}
			int delta = 0;
			for (int k = 0; k<top_search_range; k++)
				if (y.bit(base + k) != b.bit(base + k))
					delta++;
			//vOut.push_back(y);
			if(pre_vet(y, idx))
				continue;

			temp[hit_cnt++]=y;
		}
		size_t n = hit_cnt;
		if(n>0)
		{
			//concurrent_vector<bigint>::iterator p = vOut.grow_by(1<<top_search_range);
			omp_set_lock(&op_count_lock);
			size_t sizePrev = vOut.size();
			vOut.resize(vOut.size()+n);
			memcpy(&vOut[sizePrev], &temp[0], n*sizeof(bigint));
			omp_unset_lock(&op_count_lock);
			//memcpy((bigint*)p, temp, (1<<top_search_range)*sizeof(bigint));
		}
		
		delete[] temp;
	}

	void do_random_flip_sideways(bigint b) const
	{
		int half_len = (b.len + 1) / 2;

		vector<int>& vPos = m_vPos;
		vector<pair<bigint,int> >& vMinusOnes = m_vMinusOnes;
		int free_bits = m_free_bits;
		if(free_bits>nflips)
			free_bits=nflips;
		if(vPos.size()<nflips-free_bits)
			return;
		int nAllowed = vPos.size();

		vector<bigint> temp;
		int tests = tests_per_call;
		if(nAllowed<=10)
		{
			uint64_t nPerm = combin(nflips-free_bits, nAllowed);
			for(int i=0; i<free_bits; i++)
				nPerm *= half_len;
			if(nPerm*2<tests)
				tests=nPerm*2;
		}

		temp.resize(tests);
		int i;
		int rerun=0;
		for(i=0; i<tests; i++)
		{
			bigint y0 = b;

			bigint c;
			while(true)
			{
				c.clear();
				for(int j=0; j<nflips-free_bits; j++)
				{
					int pos = vPos[rand64() % nAllowed];
					c.set(pos);
				}
				for(int j=0; j<free_bits; j++)
					c.set(rand64() % half_len);
//				c = c & allowed;
				if(popcnt(c)==nflips)
					break;
			}

			while(!c.zero())
			{
				int pos = c.trailing_bit();
				y0.flip(pos);
				if(y0.len-1-pos!=pos)
					y0.flip(y0.len-1-pos);
				c.unset(pos);
			}
			if(pre_vet(y0))
			{
				rerun++;
				if(rerun>=1000)
					break;
				i--;
				continue;
			}
			rerun=0;
			if (diag && y0.bit(0))
			{
				y0.neg();
				y0.unset_from(y0.len);
			}
			temp[i]=y0;
			//temp.push_back(y0);
		}
//		if(i==tests)
//			printf("Got %d tests\n", i);
//		else
//			printf("Timeout at %d/%d tests\n", i, tests);
		temp.resize(i);
		rem_dup(temp);
		//printf("%d/%d ", temp.size(), tests);
		//std::copy( temp.begin(), temp.end(), vOut.grow_by(temp.size()) );
		omp_set_lock(&op_count_lock);
		size_t n = temp.size();
		size_t sizePrev = vOut.size();
		vOut.resize(vOut.size()+n);
		memcpy(&vOut[sizePrev], &temp[0], n*sizeof(bigint));
		omp_unset_lock(&op_count_lock);
	}

	mutable vector<int> m_vPos;
	mutable vector<pair<bigint,int> > m_vMinusOnes;
	mutable vector<bigint> m_vPotential[2];
	mutable vector<bigint> m_vLowHalf[2];


	void do_random_flip(bigint b, int idx) const
	{
		int half_len = (b.len + 1) / 2;
		int max_length = b.len + fwd_search_range;
		int max_half_length = max_length / 2;

		vector<int>& vPos = m_vPos;
		vector<pair<bigint,int> >& vMinusOnes = m_vMinusOnes;
		//m_free_bits=nflips;
		int min_range=max(min_length,b.len+1);
		int max_range=min(max_length, GRAPH_MAX_LEN * 64);
		int half_min_range=(min_range+1)/2;

		vector<bigint> temp;
		//for(int fl=nflips; fl<(nflips==1 ? 2 : nflips+2); fl++)
		int fl = nflips;
		{
			int free_bits = m_free_bits;
			if(free_bits>fl)
				free_bits=fl;

			if(vPos.size()<fl-free_bits)
				//continue;
					return;
			int nAllowed = vPos.size();

			int tests = tests_per_call;
			uint64_t nPerm = 1;
			if(fl<=5)
			{
				nPerm = combin(fl-free_bits, nAllowed);
				for(int i=0; i<free_bits; i++)
					nPerm *= half_len;
				tests=min(tests,nPerm*2*(max_range-min_range));
			}


			int skips=0;
			if(fl==1 && tests >= nPerm*(max_range-min_range)) 
			{
				for(int i=0; i<nPerm; i++)
				{
					bigint y0;
					y0 = b;
					y0.unset_from((min_range+1)/2);
					y0.set_len((min_range+1)/2);
					int pos = (free_bits==1) ? i : vPos[i];
					if(pos<y0.len)
						y0.flip(pos);

					if(pre_vet_low_half(y0))
					{
						skips+=max_range-min_range;
						continue;
					}

//					bool ok = test_build(*g0, y0, k0, k1);
//					if(!ok)
//						continue;
					for (int new_len = min_range; new_len<max_range; new_len++)
					{
						bigint y0 = reflect(b, new_len);

						y0.flip2(pos);
						if(pre_vet(y0))
						{
							skips++;
							continue;
						}
						if (diag && y0.bit(0))
						{
							y0.neg();
							y0.unset_from(y0.len);
						}
						bigint xinv;
						invert(xinv, y0);
						if(!bit_equal(xinv, y0))
						{
							printf("Invert 1 fail\n");
							exit(-1);
						}
						temp.push_back(y0);
					}
				}
//				if(idx==1)
//					printf("range %d .. %d, %d bits, nPerm %d, %d candidates\n", min_range, max_range, vPos.size(), nPerm, temp.size());
			}
			else
			{
				int rerun=0;
				while(temp.size()<tests && rerun<100)
				{
					bigint check;
					while(true)
					{
						check.clear();
						for(int j=0; j<fl-free_bits; j++)
						{
							int pos = vPos[rand64() % nAllowed];
							check.set(pos);
						}
						for(int j=0; j<free_bits; j++)
//							check.set(rand64() % max_half_length);
							check.set(rand64() % half_len);

						if(popcnt(check)==fl)
							break;
					}

					bigint y0;
					y0 = b;
					y0.unset_from((min_range+1)/2);
					y0.set_len((min_range+1)/2);
					bigint c = check;
					while(!c.zero())
					{
						int pos = c.trailing_bit();
						if(pos<y0.len)
							y0.flip(pos);
						c.unset(pos);
					}
					if(pre_vet_low_half(y0))
					{
						rerun++;
						continue;
					}

					size_t sizePrev=temp.size();
					for (int new_len = min_range; new_len<max_range; new_len++)
					{
						bigint xinv;
#if 0
						bigint y0 = reflect(b, new_len);
						bigint c = check;
						invert(xinv, y0);
						if(!bit_equal(xinv, y0))
						{
							printf("Invert 2 fail\n");
							exit(-1);
						}

						while(!c.zero())
						{
							int pos = c.trailing_bit();
							y0.flip2(pos);
							c.unset(pos);
						}
#else
						bigint b_ = b;
						bigint c = check;
						while(!c.zero())
						{
							int pos = c.trailing_bit();
							b_.flip2(pos);
							c.unset(pos);
						}
						bigint y0 = reflect(b_, new_len);
#endif
						if(pre_vet(y0))
						{
							skips++;
							continue;
						}
					
						if (diag && y0.bit(0))
						{
							y0.neg();
							y0.unset_from(y0.len);
						}
	//					bigint xinv;
						invert(xinv, y0);
						if(!bit_equal(xinv, y0))
						{
							printf("Invert 3 fail\n");
							exit(-1);
						}

						temp.push_back(y0);
					}
					if(temp.size()==sizePrev)
						rerun++;
					else
						rerun=0;
				}
				rem_dup(temp);
//				if(idx==1)
//					printf("range %d .. %d, %d bits, nPerm %d, %d candidates\n", min_range, max_range, vPos.size(), nPerm, temp.size());
			}
	//		printf("%d/%d  %d\n", fl, tests, temp.size());
		}
//		printf("(%d/%d/%d/%d/%d) ", nflips, m_vPos.size(), m_vPotential.size(), skips, temp.size());
		if(!temp.empty())
		{
			omp_set_lock(&op_count_lock);
			size_t n = temp.size();
			size_t sizePrev = vOut.size();
			vOut.resize(vOut.size()+n);
			memcpy(&vOut[sizePrev], &temp[0], n*sizeof(bigint));
			omp_unset_lock(&op_count_lock);
		}
	}

	void do_bit_flip(bigint b, bool isOrig, int idx) const
	{
		bool off_by_2_sideways_only = true;

		int half_len = (b.len + 1) / 2;
		int max_length = b.len + fwd_search_range;

		vector<int>& vPos = m_vPos;
		vector<pair<bigint,int> >& vMinusOnes = m_vMinusOnes;
		vector<bigint> temp;
		int skips=0;
		int max_half_length = max_length / 2;
		{
			bigint y_init = b;

			for (int pos2 = -1; pos2<max_half_length; pos2++)
				//int pos2=-1;
			{
				if (pos2 >= 0 && !off_by_2)
					break;
				if (pos2 >= 0 && fwd_only && off_by_2_sideways_only)
					break;
//				g0->reset();

				for (int pos3 = pos2; pos3<max_half_length; pos3++)
				{
					if (pos3 == pos2 && pos2 >= 0)
						continue;

					int loop_low = (fwd_only ? max(min_length,y_init.len+1) : min_length);
					int loop_high = min(max_length, GRAPH_MAX_LEN * 64);
					if(off_by_2_sideways_only && pos2 >= 0)
					{
//						if(!isOrig)
//							continue;
						loop_low = y_init.len;
						loop_high = loop_low+1;
					}
					if(pos3>=0)
					{
						if(!do_flip_vertical)
						{
							loop_low = y_init.len;
							loop_high = loop_low+1;
						}
						else
						{
							loop_high = min(loop_high, y_init.len+fwd_flip_range);
						}
					}
					for (int new_len = loop_low; new_len<loop_high; new_len++)
					{
						if (new_len == y_init.len && !isOrig)
							continue;
						int dx = new_len - y_init.len;
						int new_half_len = (new_len + 1) / 2;
						if (new_half_len <= pos3)
							continue;
						bigint y0 = reflect(y_init, new_len);
						if (pos2 >= 0)
							y0.flip2(pos2);
						if (pos3 >= 0)
							y0.flip2(pos3);

						if(pre_vet(y0, idx))
						{
							skips++;
							continue;
						}
						if (diag && y0.bit(0))
						{
							y0.neg();
							y0.unset_from(y0.len);
						}

						temp.push_back(y0);
					}
				}
			}
		}
		
		if(g_accuracy_check && idx<10)
		{
			int good_graphs=0;
			for(int i=0; i<temp.size(); i++)
			{
				bool ok = test_build(*g0, temp[i], k0, k1);
				if(ok)
					good_graphs++;
			}

			printf("%d %d allowed bits, %d pre-vet cliques, %d skipped, %d remaining, %d good\n", idx, m_vPos.size(), m_vPotential[0].size()+m_vPotential[1].size(), skips, temp.size(), good_graphs);
		}
		omp_set_lock(&op_count_lock);
		size_t n = temp.size();
		size_t sizePrev = vOut.size();
		vOut.resize(vOut.size()+n);
		memcpy(&vOut[sizePrev], &temp[0], n*sizeof(bigint));
		omp_unset_lock(&op_count_lock);
	}
};

int CyclicSearchObject::m_free_bits=0;
int g_pool_usage=0;

class CandidateCheckObject
{
	mutable graph *g0, *g1;
	vector<vector<int> >& vind;
	vector<bigint>& vCand;
	vector<bigint>& bases;
	vector<bigint>& vKnownFull;
	vector<bigint>& vMisses;
	int& nPrevKnown;
public:
	CandidateCheckObject(vector<vector<int> >& _vind,
	vector<bigint>& _vCand,
	vector<bigint>& _bases,
	vector<bigint>& _vKnownFull,
	vector<bigint>& _vMisses,
	int& _nPrevKnown)
		:
		vind(_vind), vCand(_vCand), bases(_bases), vKnownFull(_vKnownFull), vMisses(_vMisses), nPrevKnown(_nPrevKnown)
	{
		omp_set_lock(&op_count_lock);
		g0 = g_graph_pool+g_pool_usage;
		g1 = g_graph_pool+g_pool_usage+1;
		g_pool_usage+=2;
		omp_unset_lock(&op_count_lock);
	}
	CandidateCheckObject(const CandidateCheckObject& c)
		:
		vind(c.vind), vCand(c.vCand), bases(c.bases), vKnownFull(c.vKnownFull), vMisses(c.vMisses), nPrevKnown(c.nPrevKnown)
	{
		omp_set_lock(&op_count_lock);
		g0 = g_graph_pool+g_pool_usage;
		g1 = g_graph_pool+g_pool_usage+1;
		g_pool_usage+=2;
		omp_unset_lock(&op_count_lock);
	}

    void operator()( const blocked_range<size_t>& r ) const 
	{
        for( size_t i=r.begin(); i!=r.end(); ++i ) 
		{
			do_candidate_test(i);
		}
	}

	void do_candidate_test(size_t j) const
	{
		g0->reset();
		bool good_base = test_build(*g0, bases[j], g_kmax[0], g_kmax[1]);
		if (!good_base)
		{
			for (int ii = 0; ii < vind[j].size(); ii++)
			{
				int i = vind[j][ii];
				vCand[i].len = 0;
			}
			return;
		}
		int nDbHits=0;
		for (int ii = 0; ii<vind[j].size(); ii++)
		{
			int i = vind[j][ii];
			bigint x = vCand[i];
			bool db_result = true;
			size_t pos = lower_bound(vKnownFull.begin(), vKnownFull.end(), x) - vKnownFull.begin();
			if (pos == vKnownFull.size())
				db_result = false;
			if (db_result && !bit_equal(x, vKnownFull[pos]))
				db_result = false;
			bool new_hit = false;

			if (!db_result)
			{
				bigint xinv;
				invert(xinv, x);
				if(!bit_equal(xinv, x))
				{
					printf("Warning: candidate %d not symmetric\n", i);
					x.print();
					xinv.print();
					exit(-1);
					vCand[i].len = 0;
					continue;
				}
				bool test_result = good_base && test_build(*g1, *g0, x, g_kmax[0], g_kmax[1]);
				if (test_result && !db_result/* && x.len + 1 >= min_length + 1*/)
				{
					omp_set_lock(&op_count_lock);
					vMisses.push_back(x);
					omp_unset_lock(&op_count_lock);
					new_hit = true;
				}
				if (db_result && !test_result)
				{
					printf("Error: %d 0x%llx 0x%llx: in the DB but getting test_build fail\n", x.len + 1, x.n[0], x.n[1]);
					exit(-1);
				}
			}
			else
			{
				nDbHits++;
			}
			if (!new_hit)
				vCand[i].len = 0;
		}
		omp_set_lock(&op_count_lock);
		nPrevKnown+=nDbHits;
		omp_unset_lock(&op_count_lock);
	}
};

void test_candidate_cyclics(vector<bigint>& vMisses,
							vector<bigint>& vCand, int min_length, 
							vector<bigint>& vKnownFull, 
							const bigint* pSrc, int len,
							bool quiet_mode
							)
{
	if(g_accuracy_check)
		quiet_mode=false;
	if(!quiet_mode)
		printf("... ");
	uint64_t t0, t1;
	t0=__rdtsc();
	fflush(stdout);
	vector<bigint> bases = vCand;
#pragma omp parallel for
	for(int thr=0; thr<nthreads; thr++)
	{
		for (int i = thr; i<vCand.size(); i+=nthreads)
		{
			bigint& x = bases[i];
			x.unset_from(min_length / 2 - 5);
			x.set_len(min_length / 2 - 5);
		}
	}
	rem_dup(bases);
	vector<vector<int> > vind;
	vind.resize(bases.size());
	vector<int> vi;
	vi.resize(vCand.size());
#pragma omp parallel for
	for(int thr=0; thr<nthreads; thr++)
	{
		for (int i = thr; i<vCand.size(); i+=nthreads)
		{
			bigint x = vCand[i];
			x.unset_from(min_length / 2 - 5);
			x.set_len(min_length / 2 - 5);
			vi[i] = lower_bound(bases.begin(), bases.end(), x) - bases.begin();;
		}
	}
	for (int i = 0; i<vCand.size(); i++)
		vind[vi[i]].push_back(i);
	
	if(!quiet_mode)
		printf("%d candidates, %d bases @ %d, ", vCand.size(), bases.size(), min_length/2-5);

	concurrent_vector<bigint> vCandConcur;
	int nPrevKnown=0;
	g_pool_usage=0;
	size_t sizePrev=vMisses.size();
	CandidateCheckObject obj(vind, vCand, bases, vKnownFull, vMisses, nPrevKnown);
	if((g_kmax[0]==11 && g_kmax[1]==3)
		|| (g_kmax[0]==4 && g_kmax[1]==9))
	{
		parallel_for(blocked_range<size_t>(0, bases.size()), obj);
	}
	else
	{
		simple_partitioner ap;
		int grain = bases.size()/(10*nthreads);
		if(grain<1)
			grain=1;
		parallel_for(blocked_range<size_t>(0, bases.size(), grain), obj, ap);
	}
	t1=__rdtsc();
	if(!quiet_mode)
		printf("%d db hits, %d new hits in %.3f s\n", nPrevKnown, vMisses.size()-sizePrev, (t1-t0)/g_fFreq);
}

void build_prevet_lists(graph* vg, const bigint* v, vector<bigint>& vnew, int pos0, int pos1, int k0, int k1, vector<vector<int> >& vvPos, vector<vector<pair<bigint,int> > >& vvPotential, vector<pair<int,int> >& vnew_idx)
{
	vvPos.resize(pos1-pos0);
	vvPotential.resize(pos1-pos0);
#pragma omp parallel for
	for(int thr=0; thr<nthreads; thr++)
	{
		for(int i=pos0+thr; i<pos1; i+=nthreads)
		{
			bool ok = test_build(vg[i-pos0], v[i], k0, k1);
			if(!ok)
			{
				printf("Error: graph %d failed to build\n", i);
				continue;
			}
			do_minus_ones(&vg[i-pos0], v[i], vvPos[i-pos0], vvPotential[i-pos0]);
		}
	}
	
	vnew_idx.resize(vnew.size());
	for(int i=0; i<vnew.size(); i++)
		vnew_idx[i].first=-1;
#pragma omp parallel for
	for(int thr=0; thr<nthreads; thr++)
	{
		for(int i=thr; i<pos1-pos0; i+=nthreads)
		{
			vector<bigint> vr;
			vector<int> vi;
			relabel(vr, vi, v[i], k0==k1);
			for(int j=0; j<vr.size(); j++)
			{
				for(int k=0; k<vnew.size(); k++)
				{
					if(vnew[k]==vr[j])
					{
						vnew_idx[k].first=i;
						vnew_idx[k].second=vi[j];
					}
				}
			}
		}
	}
}

// Internal method to check correctness of low level cyclic search calls.
void doMultilevelPerformanceCheck(int k0, int k1, int start_length, bigint b)
{
	graph* pool=0;
	g_kmax[0]=k0;
	g_kmax[1]=k1;
	prealloc_graph_pool(pool, 1);

	vector<bigint> init_sym;
	init_sym.push_back(b);

	bool diagonal = (k0 == k1);

	int top_range = 5;
	int fwd_search_range = 20;

	vector<bigint> vKnownFull;

	uint64_t t0, t1;
	t0 = __rdtsc();
	int min_length = start_length - 1;

	vector<bigint>& v = init_sym;
	int cap = v.size();
	int cur_min_length = min_length;

	simple_partitioner ap;
	vector<bigint> vnew;
	relabel(vnew, &v[0], cap, (k0==k1));
	//shuffle(vnew);

	int grain = vnew.size()/(10*nthreads);
		if(grain<1)
			grain=1;

	vector<bool> originals;
	originals.resize(vnew.size());
	for(int i=0; i<originals.size(); i++)
		originals[i]=false;

	for(int i=0; i<originals.size(); i++)
		for(int j=0; j<cap; j++)
			if(bit_equal(vnew[i],v[j]))
			{
				originals[i]=true;
				break;
			}

	graph* vg = g_graph_pool;
#ifndef DIRECT_PRE_VET_MODE
	vector<vector<int> > vvPos;
	vector<vector<pair<bigint,int> > > vvPotential;
	vector<pair<int,int> > vnew_idx;
	build_prevet_lists(vg, &v[0], vnew, 0, cap, k0, k1, vvPos, vvPotential, vnew_idx);
#endif

	vector<bigint> vMisses;
	vector<bigint> vCand;

	g_occupied.clear();
	CyclicSearchObject obj(&vnew[0], originals, k0, k1, fwd_search_range, top_range, false, cur_min_length, vCand);
#ifndef DIRECT_PRE_VET_MODE
	obj.set_graphs(
		&vnew_idx[0],
		&vvPos[0],
		&vvPotential[0]);
#endif
	obj.do_flip_vertical=true;
	obj.fwd_flip_range=fwd_search_range;
	parallel_for(blocked_range<size_t>(0, vnew.size(), grain), obj, ap);
	size_t s = vMisses.size();
	test_candidate_cyclics(vMisses, vCand, cur_min_length, 
						vKnownFull, &v[0], cap, false);
//	rm_relabelings(vMisses, k0==k1);
}

void cyclic_ladder_search(int k0, int k1, int start_length, vector<bigint>& init_sym, int min_log_len)
{
	graph* pool=0;
	g_kmax[0]=k0;
	g_kmax[1]=k1;
	prealloc_graph_pool(pool, 1);

	bool diagonal = (k0 == k1);
	string fnFwd = format_sym_name();

	int top_range = 5;
	int fwd_search_range = 20;
	int gap_parameter = 7;

	printf("cyclic_ladder_search: %d graphs\n", init_sym.size());

	if (init_sym.size() == 0)
		return;

	vector<bigint> vKnownFull;
	relabel(vKnownFull, init_sym, diagonal);

	vector<bigint> v, vfull;
	readsym(fnFwd, v, min_log_len);
	relabel(vKnownFull, v, diagonal);

	uint64_t t0, t1;
	t0 = __rdtsc();
	int min_length = start_length - 1;

	vector<vector<bigint> > vv;
	clusterize(vv, init_sym, true);

	int pass = 0;
	int fixed_min_len = 0;
	vector<int> max_seen_len;
	vector<int> static_count;
	for (int i = 0; i < vv.size(); i++)
	{
		max_seen_len.push_back(0);
		static_count.push_back(0);
	}

	//graph vg[200];
	graph* vg = g_graph_pool;
	
	while (true)
	{
//		if(vv[0].size()>0)
//			printf("%d 0x%llx 0x%llx 0x%llx 0x%llx\n",
//				vv[0][0].len, vv[0][0].n[0], vv[0][0].n[1], vv[0][0].n[2], vv[0][0].n[3] );
		for (int z = 0; z < vv.size(); z++)
		{
			vector<bigint>& v = vv[z];
			if (v.empty())
				continue;
			int cap = v.size();
//			if (pass != 0)
				cap = min(cap, 200);
//			if(min_length==max_seen_len[z])
//				cap=min(2000, v.size());
#if 0
			if(v.size()>10000)
				cap=2000;
#endif		
			int cur_min_length = max(min_length, max_seen_len[z] - gap_parameter - 1);

			int max_ref_len=0;
			for(int i=0; i<cap; i++)
				max_ref_len=max(max_ref_len, v[i].len+1);

			int nAtMax=0, nSubMax=0;
			for(int i=0; i<cap; i++)
			{
				if(v[i].len+1==max_ref_len)
					nAtMax++;
				else
					nSubMax++;
			}
			if(cap>nAtMax)
				cap=nAtMax;


			simple_partitioner ap;
			vector<bigint> vnew;
#if 1
			relabel(vnew, &v[0], cap, (k0==k1));
#else
			for (int i = 0; i < v.size(); i++)
			{
				bigint x = v[i];
				if (diagonal && x.bit(0))
				{
					x.neg();
					x.unset_from(x.len);
				}
				vnew.push_back(x);
				for (int j = 2; j * 2 <= x.len; j++)
					//for (int j = 2; j <= 2; j++)
				{
					if (gcd(j, x.len + 1) != 1)
						continue;
					bigint y;
					y.clear();
					y.set_len(x.len);
					int p = j;
					for (int k = 0; k < x.len; k++)
					{
	//					if (x.bit((((k + 1)*j) % (x.len + 1)) - 1))
	//						y.set(k);
						if(x.bit(p-1))
							y.set(k);
						p+=j;
						if(p>=x.len+1)
							p-=x.len+1;
					}
					if (diagonal && y.bit(0))
					{
						y.neg();
						y.unset_from(y.len);
					}
					vnew.push_back(y);
				}
			}
		
#endif
			//vector<bigint> vnew_save;
			shuffle(vnew);

			int grain = vnew.size()/(10*nthreads);
				if(grain<1)
					grain=1;

			vector<bool> originals;
			originals.resize(vnew.size());
			for(int i=0; i<originals.size(); i++)
				originals[i]=false;

			for(int i=0; i<originals.size(); i++)
				for(int j=0; j<cap; j++)
					if(bit_equal(vnew[i],v[j]))
					{
						originals[i]=true;
						break;
					}

			bool test_mode = g_accuracy_check;
			bool do_flip_vertical=true;//(test_mode ? true : false);
			int fwd_flip_range=fwd_search_range;
#if 0
			if(cap >= 40)
			{
				if(static_count[z]==0)
					do_flip_vertical=false;
				else if(static_count[z]==1)
					fwd_flip_range=3;
				else
					fwd_flip_range=10;
			}
			if(cap>1000)
				do_flip_vertical=false;
#endif
		//	if(cap >= 40)
		//		do_flip_vertical=false;
#ifndef DIRECT_PRE_VET_MODE
			vector<vector<int> > vvPos;
			vector<vector<pair<bigint,int> > > vvPotential;
			vector<pair<int,int> > vnew_idx;
			build_prevet_lists(vg, &v[0], vnew, 0, cap, k0, k1, vvPos, vvPotential, vnew_idx);
#endif

			vector<bigint> vMisses;

			int random_min_length=cur_min_length;
//			if(nAtMax>23)
//				random_min_length=max(cur_min_length, max_ref_len);

			//printf("%d graphs, %d relabeled, %d at max\n", cap, vnew.size(), nAtMax);
			bool moved=false;
			int new_max_len=0;
			int random_tests_per_call = min(2000, 5000/cap);

			for(int method=0; method<3; method++)
			{
				uint64_t t2,t3;
				t2=__rdtsc();
				vector<bigint> vCand;

				if(method==0)
				{
					g_occupied.clear();
					CyclicSearchObject obj(&vnew[0], originals, k0, k1, fwd_search_range, top_range, false, cur_min_length, vCand);
//					if(!test_mode)
//						obj.off_by_2=false;
#ifndef DIRECT_PRE_VET_MODE
//					if(test_mode)
						obj.set_graphs(
							&vnew_idx[0],
							&vvPos[0],
							&vvPotential[0]);
#endif
//					obj.set_graphs(vg, vnew_idx);
					obj.do_flip_vertical=do_flip_vertical;
					obj.fwd_flip_range=fwd_flip_range;
					parallel_for(blocked_range<size_t>(0, vnew.size(), grain), obj, ap);
				}
				else
				{
					//if(cap>1000)
					//int random_tests_per_call = 100;
					//int random_tests_per_call = 500;
					g_occupied.clear();
					CyclicSearchObject obj(&vnew[0], originals, k0, k1, fwd_search_range,
						random_min_length, method+1, random_tests_per_call, vCand);
#ifndef DIRECT_PRE_VET_MODE
					obj.set_graphs(
						&vnew_idx[0],
						&vvPos[0],
						&vvPotential[0]);
#endif
					if(method==3)
						CyclicSearchObject::m_free_bits=1;
					parallel_for(blocked_range<size_t>(0, vnew.size(), grain), obj, ap);
					CyclicSearchObject::m_free_bits=0;
				}
				size_t s = vMisses.size();
				test_candidate_cyclics(vMisses, vCand, cur_min_length, 
									vKnownFull,
									&v[0], cap, test_mode ? false : true);
				t3 = __rdtsc();


				int max_len=0;
				for (int i = s; i < vMisses.size(); i++)
					max_len=max(max_len, vMisses[i].len+1);

				if(method==0 && max_len>max_ref_len)
				{
					moved=true;
					printf(" (f) ");
				}
				if(method>0 && max_len>max_ref_len)
				{
					moved=true;
					printf(" (r) ");
				}
				rm_relabelings(vMisses, k0==k1);
				int nFwd=0, nSideways=0;
				for(int i=0; i<vMisses.size(); i++)
				{
					if(vMisses[i].len+1==max_ref_len)
						nSideways++;
					else if(vMisses[i].len+1>max_ref_len)
						nFwd++;
				}
				double tPass = (t3-t2)/g_fFreq;
				/*
				if(cap==200)
				{
					if(method==0)
					{
						FILE* f=fopen("c:\\temp\\ml-fixed.txt", "a");
						fprintf(f, "%d_%d\t%d\t%d\t%d\t%.3f\n", k0+2,k1+2,max_ref_len, nFwd, nSideways, tPass);
						fclose(f);
					}
					else
					{
						FILE* f=fopen("c:\\temp\\ml-random.txt", "a");
						fprintf(f, "%d_%d\t%d\t%d\t%d\t%d\t%.3f\n", k0+2,k1+2,max_ref_len, method, random_tests_per_call, nFwd, tPass);
						fclose(f);
					}
				}
				*/
				if(moved)
				{
					new_max_len=max(new_max_len, max_len);
					break;
				}
			}

			if(new_max_len<=max_ref_len+2)
			{
				vector<bigint> vCand;

				g_occupied.clear();
				CyclicSearchObject obj0(&v[0], k0, k1, 2, 1000, vCand);
#ifndef DIRECT_PRE_VET_MODE
				vnew_idx.resize(nAtMax);
				for(int i=0; i<nAtMax; i++)
				{
					vnew_idx[i].first = i;
					vnew_idx[i].second = 1;
				}
				obj0.set_graphs(
					&vnew_idx[0],
					&vvPos[0],
					&vvPotential[0]);
#endif
				parallel_for(blocked_range<size_t>(0, nAtMax, 1), obj0, ap);

				if(cap<200)
				{
					CyclicSearchObject::m_free_bits=1;
					g_occupied.clear();
					CyclicSearchObject obj1(&v[0], k0, k1, 3, 10000, vCand);
						//(&vnew[0], originals, k0, k1, fwd_search_range, top_range, false, cur_min_length, vCand);
#ifndef DIRECT_PRE_VET_MODE
					obj1.set_graphs(
						&vnew_idx[0],
						&vvPos[0],
						&vvPotential[0]);
#endif
					parallel_for(blocked_range<size_t>(0, nAtMax, 1), obj1, ap);
				}
				if(cap<50)
				{
					g_occupied.clear();
					CyclicSearchObject obj2(&v[0], k0, k1, 4, 20000, vCand);
						//(&vnew[0], originals, k0, k1, fwd_search_range, top_range, false, cur_min_length, vCand);
#ifndef DIRECT_PRE_VET_MODE
					obj2.set_graphs(
						&vnew_idx[0],
						&vvPos[0],
						&vvPotential[0]);
#endif
					parallel_for(blocked_range<size_t>(0, nAtMax, 1), obj2, ap);
				}
				size_t s = vMisses.size();
				if(test_mode)
					printf("Horizontal (cap %d): ", cap);
				test_candidate_cyclics(vMisses, vCand, cur_min_length, 
									vKnownFull, &v[0], cap, test_mode ? false : true);
				rm_relabelings(vMisses, k0==k1);
				CyclicSearchObject::m_free_bits=0;
			}


			if(!vMisses.empty())
			{
				append(v, vMisses);
				vector<bigint> vnew, vOut;
				relabel(vnew, vMisses, k0==k1);
				vOut.resize(vnew.size()+vKnownFull.size());
				std::merge(vnew.begin(), vnew.end(), vKnownFull.begin(), vKnownFull.end(), vOut.begin());
				vKnownFull.swap(vOut);
				writesym(vMisses.begin(), vMisses.end(), fnFwd.c_str(), min_log_len);
			}

			int cur_max_len=0;
			for (int i = cap; i < v.size(); i++)
				cur_max_len=max(cur_max_len, v[i].len+1);

			vector<bigint> v2;
			if(v.size()>0)
			{
				int nMax = 0;
				for (int i = v.size()-1; i>=cap; i--)
				{
					if (v[i].len+1 == cur_max_len)
						nMax++;
				}
				int step = (nMax>1000) ? 2 : 1;

				for (int i = v.size()-1; i>=cap; i-=step)
				{
					if (v[i].len+1 == cur_max_len)
						v2.push_back(v[i]);
				}

				int nNonMax=0;
				for (int i = v.size()-1; i>=cap; i--)
				{
					if (v[i].len >= min_length && v[i].len >= cur_min_length && v[i].len+1 != cur_max_len)
						nNonMax++;
				}

				step = (nNonMax>1000) ? 2 : 1;
				for (int i = v.size()-1; i>=cap; i-=step)
				{
					if (v[i].len >= min_length && v[i].len >= cur_min_length && v[i].len+1 != cur_max_len)
						v2.push_back(v[i]);
				}
			}
			v.swap(v2);

			if(cur_max_len == max_seen_len[z])
				static_count[z]++;
			else
				static_count[z]=0;
			max_seen_len[z] = max(max_seen_len[z], cur_max_len);
			printf("Cluster %d pass %d: %d new graphs; %d graphs to go; best overall %d, best this pass %d\n", z, pass, vMisses.size(), v.size(), max_seen_len[z], cur_max_len);
		}
		double avgMaxLen=0;
		int nRemCl=0;
		for (int i = 0; i < vv.size(); i++)
		{
			if (vv[i].size()>0)
			{
				avgMaxLen+=max_seen_len[i];
				nRemCl++;
			}
		}
		if(nRemCl>0)
			avgMaxLen/=nRemCl;
		FILE* f = fopen("c:\\temp\\multilevel-trace.txt", "a");
		double tElapsed=(__rdtsc()-t0)/g_fFreq;
		if(nRemCl>=vv.size()/2)
			fprintf(f, "%.2f\t%.2f\t%.1f\t%d\t", tElapsed/vv.size(), avgMaxLen, tElapsed, pass);
		else
			fprintf(f, "%.2f\t-\t%.1f\t%d\t", tElapsed/vv.size(), tElapsed, pass);
		for (int i = 0; i < vv.size(); i++)
		{
			if (vv[i].size()>0)
				fprintf(f, "%d\t", max_seen_len[i]);
			else
				fprintf(f, "-\t");
		}
		fprintf(f, "\n");
		fclose(f);

		int m = INT_MAX;
		for (int i = 0; i < vv.size(); i++)
		{
			if (vv[i].size() == 0)
				continue;
			m = min(m, max_seen_len[i]);
		}

		int max_overall_len = 0;
		for (int i = 0; i < vv.size(); i++)
			max_overall_len = max(max_overall_len, max_seen_len[i]);

		if (m == INT_MAX)
		{
			FILE* f = fopen("c:\\temp\\multilevel.txt", "a");
			fprintf(f, "abort at pass %d, maxlen %d, 0 graphs left (%.3f s)\n", pass, max_overall_len, (__rdtsc() - t0) / g_fFreq);
			fclose(f);
			break;
		}

		int prev_len = min_length;
		min_length = max(min_length, m - gap_parameter - 1);
		if (prev_len == min_length)
			fixed_min_len++;
		else
			fixed_min_len = 0;
		if(fixed_min_len>=2 && min_length==prev_len && min_length+1<GRAPH_MAX_LEN*64)
		{
			min_length++;
			fixed_min_len=0;
		}

		printf("New min_length %d\n", min_length);
		if(min_length!=prev_len)
			refilter(vKnownFull, min_length);


		for (int z = 0; z < vv.size(); z++)
			refilter(vv[z], min_length);

		int n = 0;
		for (int i = 0; i < vv.size(); i++)
			n += vv[i].size();

		if (pass >= 100)
		{
			FILE* f = fopen("c:\\temp\\multilevel.txt", "a");
			fprintf(f, "abort at pass %d, maxlen %d, %d graphs left (%.3f s)\n", pass, max_overall_len, n, (__rdtsc() - t0) / g_fFreq);
			fclose(f);
			break;
		}

		if (n == 0)
		{
			FILE* f = fopen("c:\\temp\\multilevel.txt", "a");
			fprintf(f, "abort at pass %d, maxlen %d, %d graphs left (%.3f s)\n", pass, max_overall_len, n, (__rdtsc() - t0) / g_fFreq);
			fclose(f);
			break;
		}
		pass++;
	}
}

void cyclic_network_search(int k0, int k1, int ml, bool fwd_only, int min_log, int max_passes)
{
	graph* pool=0;
	g_kmax[0]=k0;
	g_kmax[1]=k1;
	int max_graphs_per_step=400;
	if((k0==4 && k1==7) || (k0==9 && k1==3) && nthreads>8)
		max_graphs_per_step=4000;
	prealloc_graph_pool(pool, 1, max_graphs_per_step*GRAPH_MAX_LEN*32);
	vector<bigint> vKnown;
	bool quiet_mode = false;
	int min_length = ml - 1;
	if(min_log==0)
		min_log=ml;
	if(!quiet_mode)
		printf("cyclic_network_search(%d,%d,%d)\n", k0, k1, ml);
	char fnFwd[64];
	sprintf(fnFwd, "c:\\temp\\sym%d%d.txt", k0+2, k1+2);
	char covername[64];
	sprintf(covername, "c:\\temp\\cover%d%d_%d%s.txt", k0+2, k1+2, ml, fwd_only ? "fwd" : "");

	bool diag = (k0 == k1);
	int top_range = 10;
	int fwd_search_range = (fwd_only || ml>=250) ? 50 : 30;
	bool off_by_2 = true;
	bool off_by_2_sideways_only = true;

	readsym(fnFwd, vKnown, min_length + 1);

	g_kmax[0] = k0;
	g_kmax[1] = k1;
	uint64_t ops = 0;
	sort(vKnown.begin(), vKnown.end());
	if(!quiet_mode)
		printf("%d known graphs\n", vKnown.size());
	rem_dup(vKnown);
	if(!quiet_mode)
		printf("%d unique\n", vKnown.size());

	vector<bigint> vKnownFull;
	relabel(vKnownFull, vKnown, k0==k1);

	vector<bigint> vCover;
	readsym(covername, vCover, 0);
	if(!quiet_mode)
		printf("%d previously processed\n", vCover.size());

	for (int i = 0; i<vCover.size(); i++)
	{
		vector<bigint> v;
		relabel(v, vCover[i], k0==k1);
		for (int j = 0; j < v.size(); j++)
		{
			size_t pos = lower_bound(vKnown.begin(), vKnown.end(), v[j]) - vKnown.begin();
			if (pos != vKnown.size() && bit_equal(vKnown[pos], v[j]))
				vKnown[pos].len = 0;
		}
	}
	refilter(vKnown);

	uint64_t t0 = __rdtsc();

	if(fwd_only)
	{
		for(int j=0; j<vKnown.size(); j++)
		{
			if(vKnown[j].len!=min_length)
				vKnown[j].len=0;
		}
		refilter(vKnown);
	}
	if(!quiet_mode)
		printf("%d left\n", vKnown.size());
	if(vKnown.empty())
		return;
//	graph m_graphs[nthreads], m_graphs0[nthreads], m_graphs1[nthreads];

	vector<bigint> v;
	v = vKnown;
	vKnown.clear();
	random_shuffle(v.begin(), v.end());
	int pos0 = 0, pos1 = v.size();
	if (pos1 > 40)
		pos1 = 40;

	graph* vg = g_graph_pool;

	int pass = 0;
	while (true)
	{		
		time_t t = time(0);
		struct tm* lt = localtime(&t);
		if(lt->tm_hour>=10 && lt->tm_hour<21 && lt->tm_wday>=1 && lt->tm_wday<=5)
		{
			Sleep(100000);
			continue;
		}
		if(quiet_mode)
			printf(".");
		uint64_t t1, t2;
		t1 = __rdtsc();

		vector<bigint> vCand;
		//concurrent_vector<bigint> vCandConcur;
		vector<bigint> vnew;
		relabel(vnew, &v[pos0], pos1-pos0, (k0==k1));
		shuffle(vnew);
		vector<bool> originals;
		originals.resize(vnew.size());
		for(int i=0; i<originals.size(); i++)
			originals[i]=false;

		for(int i=0; i<originals.size(); i++)
			for(int j=pos0; j<pos1; j++)
				if(bit_equal(vnew[i],v[j]))
				{
					originals[i]=true;
					break;
				}
#ifndef DIRECT_PRE_VET_MODE
		vector<vector<int> > vvPos;
		vector<vector<pair<bigint,int> > > vvPotential;
		vector<pair<int,int> > vnew_idx;
		
		build_prevet_lists(vg, &v[pos0], vnew, 0, pos1-pos0, k0, k1, vvPos, vvPotential, vnew_idx);
#endif

		g_occupied.clear();
		CyclicSearchObject obj(&vnew[0], originals, k0, k1, fwd_search_range, top_range, fwd_only, min_log-1, vCand);
		obj.fwd_flip_range=20;

#ifndef DIRECT_PRE_VET_MODE
		obj.set_graphs(
			&vnew_idx[0],
			&vvPos[0],
			&vvPotential[0]);
#endif
		if((k0==11 && k1==3)
			|| (k0==4 && k1==9))
		{
			// memory constrained scenario, avoid search object glut
			parallel_for(blocked_range<size_t>(0, vnew.size()), obj);
		}
		else
		{
			simple_partitioner ap;
			int grain = vnew.size()/(10*nthreads);
			if(grain<1)
				grain=1;
			parallel_for(blocked_range<size_t>(0, vnew.size(), grain), obj, ap);
		}
		vector<bigint> vMisses;
		test_candidate_cyclics(vMisses, vCand, min_length, 
							vKnownFull, 
							&v[pos0], pos1-pos0, quiet_mode);
		rm_relabelings(vMisses, k0==k1);
		size_t sizePrev = v.size();
		if(!fwd_only)
			append(v, vMisses);
		int nNew = v.size() - sizePrev;

		int n130=0;
		int nMaxNew=0;
		for (auto j = vMisses.begin(); j != vMisses.end(); j++)
		{
			bigint x = *j;
			if(vMisses.size()<10 && !quiet_mode)
				printf("Found a missed coloring: %d 0x%llx 0x%llx\n", x.len + 1, x.n[0], x.n[1]);
			if(x.len+1>=(k0==13 ? 130 : 180))
				n130++;
			nMaxNew=max(nMaxNew,x.len+1);
		}
		if(vMisses.size()>=10 && !quiet_mode)
			printf("Found %d missed colorings (max len %d)\n", vMisses.size(), nMaxNew);

		writesym(vMisses.begin(), vMisses.end(), fnFwd);
		writesym(v.begin()+pos0, v.begin()+pos1, covername);
		pos0 = pos1;
		pos1 = v.size();
		if (pos1>pos0 + max_graphs_per_step)
			pos1 = pos0 + max_graphs_per_step;
		if(1)
		{
			vector<bigint> vnew, vOut;
			relabel(vnew, vMisses, k0==k1);
			vOut.resize(vnew.size()+vKnownFull.size());
			std::merge(vnew.begin(), vnew.end(), vKnownFull.begin(), vKnownFull.end(), vOut.begin());
			vKnownFull.swap(vOut);
		}
		t2 = __rdtsc();
		if(fwd_only)
			nNew=n130;
		if(!quiet_mode)
			printf("Pass %d: %d / %d new graphs in %.3f s; %d to go\n", pass, vMisses.size(), nNew, (t2 - t1) / g_fFreq, v.size()-pos0);
		if(pass==max_passes)
			break;
		pass++;
		if (pos0 == pos1)
			break;
	}
	if(!quiet_mode)
		printf("Done (%lld ops), %.3f s\n", ops, (__rdtsc() - t0) / g_fFreq);
}

void random_cyclic_search(int k0, int k1, int ml, int nflips, int tests_per_call, bool forward=true)
{
	graph* pool=0;
	g_kmax[0]=k0;
	g_kmax[1]=k1;
	prealloc_graph_pool(pool, 1);

	vector<bigint> vKnown;
	bool quiet_mode = false;
	int min_length = ml - 1;
	int min_log=ml;

	char fnFwd[64];
	sprintf(fnFwd, "c:\\temp\\sym%d%d.txt", k0+2, k1+2);

	bool diag = (k0 == k1);
	int fwd_search_range = 50;

	readsym(fnFwd, vKnown, min_length + 1);

	g_kmax[0] = k0;
	g_kmax[1] = k1;
	uint64_t ops = 0;
	sort(vKnown.begin(), vKnown.end());
	if(!quiet_mode)
		printf("%d known graphs\n", vKnown.size());
	rem_dup(vKnown);
	if(!quiet_mode)
		printf("%d unique\n", vKnown.size());
	refilter(vKnown);
	vector<bigint> vKnownFull;
	relabel(vKnownFull, vKnown, k0==k1);

	uint64_t t0 = __rdtsc();

	for(int j=0; j<vKnown.size(); j++)
	{
		if(vKnown[j].len!=min_length)
			vKnown[j].len=0;
	}
	refilter(vKnown);

	if(!quiet_mode)
		printf("%d left\n", vKnown.size());
	if(vKnown.size()==0)
		return;
	graph m_graphs[nthreads], m_graphs0[nthreads], m_graphs1[nthreads];
	vector<bigint> v;
	v = vKnown;
	vKnown.clear();
	random_shuffle(v.begin(), v.end());
	int pos0 = 0, pos1 = v.size();
	if (pos1 > 40)
		pos1 = 40;

	int pass = 0;
	int nNewHits=0;
	while (true)
	{
		/**
		time_t t = time(0);
		struct tm* lt = localtime(&t);
		if(lt->tm_hour>=11 && lt->tm_hour<20)
		{
			Sleep(100000);
			continue;
		}
		**/
		uint64_t t1, t2;
		t1 = __rdtsc();



		//concurrent_vector<bigint> vCandConcur;
		vector<bigint> vCand;
		simple_partitioner ap;
		g_occupied.clear();
		if(forward)
		{
			vector<bigint> vnew;
			relabel(vnew, &v[pos0], pos1-pos0, (k0==k1));
			shuffle(vnew);
			vector<bool> originals;
			originals.resize(vnew.size());
			for(int i=0; i<originals.size(); i++)
				originals[i]=false;

			for(int i=0; i<originals.size(); i++)
				for(int j=pos0; j<pos1; j++)
					if(bit_equal(vnew[i],v[j]))
					{
						originals[i]=true;
						break;
					}
#ifndef DIRECT_PRE_VET_MODE
			vector<vector<int> > vvPos;
			vector<vector<pair<bigint,int> > > vvPotential;
			vector<pair<int,int> > vnew_idx;
		
			build_prevet_lists(g_graph_pool, &v[pos0], vnew, 0, pos1-pos0, k0, k1, vvPos, vvPotential, vnew_idx);
#endif

			CyclicSearchObject obj(&vnew[0], originals, k0, k1, fwd_search_range, 
				min_length, nflips, tests_per_call, vCand);

#ifndef DIRECT_PRE_VET_MODE
			obj.set_graphs(
				&vnew_idx[0],
				&vvPos[0],
				&vvPotential[0]);
#endif
			obj.fwd_flip_range=30;
			int grain = vnew.size()/(10*nthreads);
			if(grain<1)
				grain=1;
			parallel_for(blocked_range<size_t>(0, vnew.size(), grain), obj, ap);
		}
		else
		{
			CyclicSearchObject obj(&v[pos0], k0, k1, nflips, tests_per_call, vCand);
			int grain = (pos1-pos0)/(10*nthreads);
			if(grain<1)
				grain=1;
			parallel_for(blocked_range<size_t>(0, pos1-pos0, grain), obj, ap);
		}

		vector<bigint> vMisses;
		test_candidate_cyclics(vMisses, vCand, min_length, 
							vKnownFull, 
							&v[pos0], pos1-pos0, quiet_mode);
		rm_relabelings(vMisses, k0==k1);

		int n130=0;
		int nMaxNew=0;
		for (auto j = vMisses.begin(); j != vMisses.end(); j++)
		{
			bigint x = *j;
			if(vMisses.size()<10 && !quiet_mode)
				printf("Found a missed coloring: %d 0x%llx 0x%llx\n", x.len + 1, x.n[0], x.n[1]);
			if(x.len+1>=(k0==13 ? 130 : 180))
				n130++;
			nMaxNew=max(nMaxNew,x.len+1);
		}
		if(vMisses.size()>=10 && !quiet_mode)
			printf("Found %d missed colorings (max len %d)\n", vMisses.size(), nMaxNew);

		writesym(vMisses.begin(), vMisses.end(), fnFwd);
		nNewHits+=vMisses.size();
		if(!forward)
			append(v, vMisses);
		pos0 = pos1;
		pos1 = v.size();
		if (pos1>pos0 + 400)
			pos1 = pos0 + 400;
		if(1)
		{
			vector<bigint> vnew, vOut;
			relabel(vnew, vMisses, k0==k1);
			vOut.resize(vnew.size()+vKnownFull.size());
			std::merge(vnew.begin(), vnew.end(), vKnownFull.begin(), vKnownFull.end(), vOut.begin());
			vKnownFull.swap(vOut);
		}
		t2 = __rdtsc();
		if(!quiet_mode)
			printf("Pass %d: %d / %d new graphs in %.3f s; %d to go\n", pass, vMisses.size(), n130, (t2 - t1) / g_fFreq, v.size()-pos0);
		pass++;
		if (pos0 == pos1)
			break;
	}
	if(!quiet_mode)
		printf("Done, %d new graphs in %.3f s\n", nNewHits, (__rdtsc() - t0) / g_fFreq);
	FILE* f=fopen(forward ? "c:\\temp\\random_forward.txt" : "c:\\temp\\random_sideways.txt", "a");
	fprintf(f,"%d %d %d\t%d\t%d\t%d\t%d\t%.3f\n",
		k0, k1, ml, nflips, tests_per_call, v.size(), nNewHits, (__rdtsc() - t0) / g_fFreq);
	fclose(f);
}


void cyclic_ladder_search(int dim, int k0, int k1, int start_size, int mc_test_count, int min_log_len)
{
	bool diagonal = (k0 == k1);
	vector<bigint> sym_graphs;
	g_sym_graphs_ref.clear();

	g_kmax[0] = k0;
	g_kmax[1] = k1;
	readsym(format_sym_name(), sym_graphs, min_log_len);
	relabel(g_sym_graphs_ref, sym_graphs, k0 == k1);
	rem_dup(g_sym_graphs_ref);

	vector<b128> hits;
	vector<bigint> sym;
	g_max_ops = 1e5;
	g_recursive_cyclic_search = false;
	max_errors_detect_only = 16;
	max_errors = 16;
	//int d = start_size / 2 - 1;
	mc_memory(dim, k0, k1, start_size, mc_test_count, hits, sym, min_log_len);
	if (hits.size() == 0)
		return;
	vector<b128> cover = hits;
	while (true)
	{
		vector<vector<bigint> > vv;
		if (sym.size() > 0)
		{
			rem_dup(sym);
			rm_relabelings(sym, diagonal);
			clusterize(vv, sym, true);
		}
		printf("%d %d %d\n", hits.size(), sym.size(), vv.size());
		if (vv.size() >= 5)
			break;
		//in_memory_neighbor_search(d, start_size, INT_MAX, hits, cover, sym, min_log_len);
		vector<b128> v;
		neighbor_depth_search(dim, k0, k1, start_size, hits, cover, v, sym, 0.1, min_log_len);
		append(hits, v);
		//append(cover, v);

	}
	if (sym.size() == 0)
		return;
	quiet = true;
	vector<vector<bigint> > vv;
	clusterize(vv, sym, true);
	/*
	if(vv.size() > 10)
	{
		sym.clear();
		for(int i=0; i<10; i++)
			append(sym, vv[i]);
	}
	*/
	//printf("%d 0x%llx 0x%llx 0x%llx\n", sym[0].len, sym[0].n[0], sym[0].n[1], sym[0].n[2]);
	cyclic_ladder_search(k0, k1, start_size, sym, min_log_len);
}



void table_gen(const int* r, const int* minval, int count)
{
	vector<vector<int> > vv;
	vector<int> maxvals;
	int maxmax=0;
	int minmin=INT_MAX;
	int i;
	for(i=0; i<count; i++)
	{
		//cleansym(r[i*2+0], r[i*2+1]);
		vector<bigint> v;
		g_kmax[0]=r[i*2+0];
		g_kmax[1]=r[i*2+1];
		readsym(format_sym_name(), v, minval[i]);
		rm_relabelings(v, (g_kmax[0]==g_kmax[1]));
		int maxval=0;
		for(int j=0; j<v.size(); j++)
			maxval=max(maxval, v[j].len+1);
		vector<int> vc;
		vc.resize(maxval+1);
		for(int j=0; j<v.size(); j++)
			vc[v[j].len+1]++;
		vv.push_back(vc);
		maxvals.push_back(maxval);
		maxmax=max(maxmax, maxval);
		minmin=min(minmin, minval[i]);
	}

// \begin{longtable}[h]{|c|c|c|c|c|c|c|c|c|}\hline
// N & $(5,9)$ & $(6,8)$ & $(5,10)$ & $(7,7)$ & $(6,9)$ & $(5,11)$ & $(7,8)$ & $(6,10)$ \\\hline
	printf("N & ");
	for(i=0; i<count; i++)
	{
		int a = r[i*2+0], b=r[i*2+1];
		if(a>b)
			swap(a,b);
		printf("$(%d,%d)$ ", a+2, b+2);
		if(i<count-1)
			printf("& ");
	}
	printf("\\\n");
	printf("\\hline\n");
	for(i=minmin; i<=maxmax; i++)
	{
		printf("%d & ", i);
		for(int j=0; j<count; j++)
		{
			if(i<minval[j])
			{
			}
			else if(i>maxvals[j])
			{
			}
			else
			{
				printf("%d ", vv[j][i]);
			}
			if(j!=count-1)
				printf(" & ");
		}
		printf("\\\\ \n");
		if((i%10)==9)
			printf("\\hline\n");
	}
}

void table_gen()
{
	int r[]=
	{
		10, 3,
		11, 3,
		4, 8,
		4, 9,
		4, 10,
		5, 6,
		5, 7,
		6, 6,
	};

	int minval[]=
	{
		180,
		205,
		176,
		210,
		232,
		160,
		224,
		180
	};
	table_gen(r, minval, 8);
}

void do154mc()
{
	vector<bigint> sym;
	vector<b128> hits;
	uint64_t nbases= 1e8, npasses=1;
	g_max_ops=1e6;
	load_known_sym_graphs(13,2,125);
	vector<b128> hits_50_125, hits_50_130;
	for(int pass=0; pass<npasses; pass++)
	{
		vector<b128> v;
		mc_memory(49,13,2,125,nbases,v,sym,125);
		append(hits, v);
	}
	printf("%lld -> %d\n", nbases*npasses, hits.size());
	quiet=true;

	//writefile("c:\\temp\\144-50-110.txt", hits, 49);
	lengthen_memory(hits, hits_50_125, 13,2,49,50,125);
	lengthen_memory(hits_50_125, hits_50_130, 13,2,50,50,130);
	printf("%d at (50,125)\n", hits_50_125.size());
	printf("%d at (50,130)\n", hits_50_130.size());

	vector<b128> hits_55_130, hits_60_130, hits_65_130;
	lengthen_memory(hits_50_130, hits_55_130, 13,2,50,55,130);
	lengthen_memory(hits_55_130, hits_60_130, 13,2,55,60,130);
	lengthen_memory(hits_60_130, hits_65_130, 13,2,60,65,130);
	printf("%d at (65,130)\n", hits_65_130.size());

	vector<b128> known, cover;
	readfile("c:\\temp\\154-66-130.txt", known, 65);
	readfile("c:\\temp\\cover154-66-130.txt", cover, 65);
	exclude(cover, known);
	vector<b128> hs = hits_65_130;
	int nEx = exclude(hits_65_130, known);
	printf("%d / %d already known\n", nEx, hs.size());
	nEx = exclude(hits_65_130, cover);
	printf("%d / %d incorrectly excluded\n", nEx, hs.size());
	exclude(hs, known);
	writefile("c:\\temp\\154-66-130.txt", hs, 65, true);
	writefile("c:\\temp\\cover154-66-130.txt", hs, 65, true);
}

void ladder_searches()
{
//	cyclic_ladder_search(54,13,2,120,4e8,126); // 6 clusters, 34 s (?)
//	cyclic_ladder_search(54,9,3,115,1e6,156); // 25 clusters, 15 min / run
//	cyclic_ladder_search(64,10,3,130,1e7,180); // 18 clusters, 1 hour / run
//	cyclic_ladder_search(69,11,3,140,1e6,205); // 12 clusters, 35 min / run (?)
//	cyclic_ladder_search(54,4,7,115,5e6,160); // 20 clusters, 15 min / run
//	cyclic_ladder_search(59,4,8,120,1e5,176); // 40 clusters, 30 min / run
//	cyclic_ladder_search(69,4,9,140,2e5,200); // 10 clusters, 20 min / run
//	cyclic_ladder_search(69,4,10,140,5e4,227); // fast MC, 20 clusters, 5-6 min/cluster 
//	cyclic_ladder_search(58,5,6,125,5e6,160); // 21 clusters, 45 min / run
//	cyclic_ladder_search(67,5,7,140,5e4,220); // 15 clusters, 35 min / run
//	cyclic_ladder_search(74,5,8,150,5e4,260); 
//	cyclic_ladder_search(71,6,6,144,5e5,175); 
}

void unified_searches()
{
/**
	unified_search(4, 8, 79, 176, 1e5, 3e4); // 2 min/pass; miss rate 10/1000 @4e3, 1/1000 @3e4
	unified_search(4, 8, 79, 180, 1e5, 3e4);
	unified_search(4, 9, 100, 205, 3e5, 3e4); // 5 min/pass; miss rate 13/1000 @1e4, 10/1000 @7e4
	unified_search(4, 9, 100, 210, 3e5, 3e4); // 5 min/pass; miss rate 20/1000 @1.5e3, 10/1000 @3e4 and @9e4
	unified_search(4, 9, 100, 225, 3e5, 3e4); // miss rate 16/1000 @3e4
	unified_search(5, 6, 74, 160, 5e6, 2e4); 
	unified_search(5, 6, 74, 165, 1e6, 5e4); // 4 min/pass; miss rate 3/1000 @5e4, 13/1000 @2e3
	unified_search(5, 7, 109, 224, 5e5, 3e4); // 6 min/pass; miss rate 22/1000 @4e3, 12/1000 @3e4, 8/1000 @9e4
	unified_search(6, 6, 85, 170, 5e4, 1e5);
	unified_search(6, 6, 95, 185, 1e5, 5e4);
	unified_search(6, 6, 105, 190, 1e5, 5e4); // 1 min/pass; miss rate 5/1000 @5e4
	unified_search(6, 6, 105, 195, 1e5, 5e4);
	unified_search(9, 3 ,63, 156, 1e5, 5e4); // 1.5 min/pass; miss rate 12/1000 @5e4, 9/1000 @1.5e5
	unified_search(10, 3, 85, 181, 2e5, 3e4); // miss rate 4/1000 @3e4
	unified_search(11, 3, 99, 213, 1e5, 2e4);
	unified_search(11, 3, 99, 224, 1e4, 2e4);
	unified_search(13, 2, 65, 130, 2e5, 2e4); // 5 min/pass; cyclic miss rate 2/500 @2e4, 1/500 @3e4
	unified_search(12, 2, 54, 115, 5e5, 2e4); // 5 min/pass; cyclic miss rate <1/1000 @2e4, 10/1000 @2e3
	unified_search(12, 2, 50, 115, 2e5, 2e4); // 10 min/pass; cyclic miss rate <1/1000 @2e4, 10/1000 @8e3
**/
}

void pre_vet_accuracy_test()
{
	g_accuracy_check=true;
	bigint b;
	
	b.clear();
	b.len=171;
	b.n[0]=0x1b4be6cf977ba	;
	b.n[1]=0x6c00006cdbded9b0	;
	b.n[2]=0x2ef74f9b3e9;
	doMultilevelPerformanceCheck(6,6,171,b);
	
	b.clear();
	b.len=148;
	b.n[0]=0x14f1a23bfffbbcfc;
	b.n[1]=0xfffdc458f2868616;
	b.n[2]=0x3f3dd;
	doMultilevelPerformanceCheck(4,10,140,b);
	
	b.len=235;
	b.n[0]=0xfe45a18185a27f95;
	b.n[1]=0x7bfef7fbfcff6ba9;
	b.n[2]=0x2d13fcaeb7f9feff;
	b.n[3]=0x54ff22d0c0c;
	doMultilevelPerformanceCheck(4,10,140,b);

	b.len=226;
	b.n[0]=0x936f1f3faf94934b;
	b.n[1]=0x7f478bfac925626c;
	b.n[2]=0xf3e3db24d91a924d;
	b.n[3]=0x34b24a7d7;
	doMultilevelPerformanceCheck(5,7,220,b);

	b.len=125;
	b.n[0]=0x1c004a0882003e1c;
	b.n[1]=0x70f8008220a4007;
	b.n[2]=0x0;
	b.n[3]=0x0;
	doMultilevelPerformanceCheck(13,2,120,b);
	
/*
	bigint b;
	b.clear();
	b.len=240;
	b.n[0]=0x2efffffff75295df;
	b.n[1]=0xf18ff8c5f462fa34;
	b.n[2]=0xff742c5f462fa31f;
	b.n[3]=0xfba94aefffff;
	doMultilevelPerformanceCheck(4,12,230,b);
*/
	g_accuracy_check=false;
}

void crosscheck(int k0, int k1, int k0_, int k1_, int min_len)
{
	g_kmax[0]=k0;
	g_kmax[1]=k1;
	vector<bigint> v, vNext;
	readsym(format_sym_name(), v, min_len);
	g_kmax[0]=k0_;
	g_kmax[1]=k1_;
	string fn = format_sym_name();
	readsym(fn, vNext, 0);
	if((k0==k1) && !(k0_==k1_))
	{
		vector<bigint> v2 = v;
		neg(v2);
		append(v, v2);
	}
	if((k0_==k1_) && !(k0==k1))
	{
		for(int i=0; i<v.size(); i++)
		{
			if(v[i].bit(0))
			{
				v[i].neg();
				v[i].unset_from(v[i].len);
			}
		}
	}

	rem_dup(v);
	rem_dup(vNext);
	vector<bigint> vr;
	relabel(vr, vNext, k0_==k1_);
	rem_dup(vr);
	printf("%d (%d,%d) graphs\n", v.size(), k0+2, k1+2);
	exclude(v, vr);
	printf("%d to test\n", v.size());
#pragma omp parallel for
	for(int thr=0; thr<nthreads; thr++)
	{
		graph g;
		for(int i=thr; i<v.size(); i+=nthreads)
		{
			g.reset();
			bool ok = test_build(g, v[i], k0_, k1_);
			if(!ok)
				v[i].len=0;
		}
	}
	refilter(v);
	printf("%d graphs left\n", v.size());
	//writesym(format_sym_name(), v);
	writesym(v.begin(), v.end(), fn.c_str());
}

int cpu_init()
{
	int nPr = GetActiveProcessorCount(ALL_PROCESSOR_GROUPS);
	SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX procinfo[1024];
	DWORD len=sizeof(procinfo);
	BOOL status = GetLogicalProcessorInformationEx(RelationAll, procinfo, &len);
	vector<int> cores;
	for(int i=0; i<len; )
	{
		SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX* pRel = (SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX*) ((uint8_t*)procinfo + i);
		if(pRel->Relationship == RelationProcessorCore) // Share core
			cores.push_back(_lzcnt_u64(pRel->Processor.GroupMask[0].Mask));
		i += pRel->Size;
	}
//	printf("%d processors, %d cores\n", nPr, cores.size());
	nthreads = min(cores.size(), max_nthreads);
	if(nthreads > 10)
		nthreads--;
	if(cores.size() < nPr)
	{
		uint64_t lMask=0;
		for(int i=0; i<nPr; i++)
			lMask |= (One<<(i*2));
		SetProcessAffinityMask(GetCurrentProcess(), lMask);
	}
	int t=time(0);
	while(time(0)==t);
	uint64_t t0=__rdtsc();
	t=time(0);
	while(time(0)==t);
	uint64_t t1=__rdtsc();
	g_fFreq = t1-t0;
	printf("%d processors, %d cores, %d threads, %.3f GHz\n", nPr, cores.size(), nthreads, g_fFreq/1e9);
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
	return 0;
}


int main(int argc, char* argv[])
{
	omp_init_lock(&op_count_lock);
	init_phi();
	cpu_init();
	cyclic_ladder_search(54,4,7,115,5e6,160);
	return 0;
}
