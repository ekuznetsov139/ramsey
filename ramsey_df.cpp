// ramsey_df.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "bigint.h"
#include "graph.h"
#include <omp.h>
#ifdef WIN32
#include <Windows.h> // InterlockedIncrement
#endif

#ifndef WIN32
#define _ftelli64 ftell
#define _unlink unlink
#endif

int g_kmax[2]={5,5};
//double g_fFreq = 3.4e9;
//static const int nthreads = 4;
#ifdef KNIGHTS_CORNER
double g_fFreq = 1.238e9;
static const int nthreads = 220;
#elif defined(WIN32)
double g_fFreq = 2.7e9;
static const int nthreads = 45;
#else
double g_fFreq = 3.6e9;
static const int nthreads = 8;
#endif

#ifdef WIN32
#define PATH "c:\\temp\\"
#else
#define PATH
#endif

#if MAX_LEN==2
#define graph_pool_dim 160
#else
#define graph_pool_dim 280
#endif

#define pop_acc_dim (MAX_LEN*64+2)

graph* g_graph_pool = 0;

#define PARALLEL


template <class X>
inline void append(X& v, const X& w)
{
	if(w.size()==0)
		return;
	size_t n = v.size();
	v.resize(n+w.size());
	memcpy(&v[n], &w[0], w.size()*sizeof(w[0]));
}

bool symmetry_check(graph& g, int kmax[2], int target, bool debug);

inline bool operator<(const pair<uint64_t,int>& a, const pair<uint64_t,int>& b) { return a<b; }

template <class T>
inline size_t rem_dup(vector<T>& v)
{
	if(v.empty())
		return 0;
	sort(v.begin(), v.end());
	size_t sizePrev = v.size();
	v.resize(unique(v.begin(), v.end())-v.begin());
	return sizePrev-v.size();
}

void prealloc_graph_pool(graph*& graph_pool, int len)
{
	if(g_graph_pool!=0)
	{
		graph_pool = g_graph_pool;
		return;
	}
	g_graph_pool=new graph[nthreads*graph_pool_dim*2];
	graph_pool = g_graph_pool;
//	printf("Attempting to alloc %
	for(int i=len-1; i<graph_pool_dim; i++)
	{
		int expect_size0 = 1000*i/50;
		if(expect_size0<1000)
			expect_size0=1000;
		int expect_size1 = expect_size0;
		int def_size = expect_size0;
		if(g_kmax[0]>=7)
			expect_size0=100000;
		if(g_kmax[1]>=7)
			expect_size1=100000;
		
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


void graph_check(graph& g, int kmax[2], const char* msg)
{
	return;

		graph g_alt2;
		bool ok = build_graph3(g_alt2, make_pair(g.mask0, g.mask), kmax);
		if(!ok)
		{
			printf("graph_check(%s): build_graph3 fail\n", msg);
			g.print_short();
			g_alt2.print_short();
			exit(-1);
		}
		else
		if(memcmp(g_alt2.clique_counts, g.clique_counts, sizeof(g_alt2.clique_counts)))
		{
			printf("%s count mismatch\n", msg);
//							printf("pos %d\n", pos);
			g.print_short();
			g_alt2.print_short();
			exit(-1);
		}
}

inline bool extend_from(graph& g, const graph& gmin, bigint2 m, int kmax[2], int target)
{
	int last_bit=max(g.mask.trailing_bit(), g.mask0.trailing_bit());
	bool first_bit=true;
	int last_known_bit=max(gmin.mask.trailing_bit(), gmin.mask0.trailing_bit());
	if(last_known_bit >= last_bit)
	{
		printf("%d %d\n", last_known_bit, last_bit);
		exit(-1);
	}


	for(int i=last_known_bit+1; i<=last_bit; i++)
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
			if(!extend(g, gmin, color, kmax, target, i, 0))
				return false;
		}
		else
		{
			if(!extend(g, color, kmax, target, i))
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
	while((g.mask0.bit(first_unset) || g.mask.bit(first_unset)) && first_unset<=MAX_LEN*64)
		first_unset++;

	if(first_unset > target+3)
		return true;

	for(int i=first_unset; i<MAX_LEN*64; i++)
	{
		m.first.unset(i);
		m.second.unset(i);
	}

				
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
		g.reset();
		if(!extend_from(g, gmin[last_known], m, kmax, target+2))
			return false;
		if(g.clique_counts[0][kmax[0]]>0 || g.clique_counts[1][kmax[1]]>0)
			return false;
	}
/*	
	if(kmax[0]!=kmax[1])
	{
		graph gg;
		bool ok = build_graph3(gg, m, kmax);
		if(!ok)
			printf("*");
		if(gg.clique_counts[0][kmax[0]]>0 || gg.clique_counts[1][kmax[1]]>0)
			printf("!");
	}
*/	
	return true;
}

uint64_t g_pos_counts[2][graph_pool_dim];
uint64_t g_mean_clique_dims[2][graph_pool_dim];

vector<int> op_counts[2];
omp_lock_t op_count_lock;

int g_max_ops=0;

int depth_search(bigint sig, int kmax[2], int target, int stop_point, bool detect_only, bool quiet, int* populations, graph* graph_pool, int thr, buffer<bigint>* pHits=0, uint64_t max_ops=0)
{
	int retval=0;

	target -= 2;

	uint64_t tStart = __rdtsc();
	int tIntReportPrev = 0;
	int popPrev = populations[target];

	int ext_target = target+2;
	if(stop_point>MAX_LEN*64+1)
		stop_point=MAX_LEN*64+1;
	int bit_positions[graph_pool_dim];
	int bit_orders[graph_pool_dim];
	//graph g[128+2], gmin[128+2];
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
	bool ok = build_graph3(g[len-1], base_mask, kmax);
	if(!ok)
		return false; 

	gmin[len-1]=g[len-1];
//	build_graph3(gmin[len-1], base_mask, kmax);
	g[len-1].n=len;
	gmin[len-1].n=len;

	int pos = len, bit=0;
	int max_len=0;

	graph& gc=g[len-1];
	gc.n=len;
	gc.check_first_unset_bit();
	int x = max_extensibility(gc, kmax, ext_target, false);
	if(x < target)
		goto done;
	else
	{
		int last_bit=max(gc.mask.trailing_bit(), gc.mask0.trailing_bit());
		if(last_bit>len)
			gc.dirty=true;
	}

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
		gc.set_extensions = false;


		int next_bit=set_count;
		
		bool dup = false;
		if(delta>=6 && (delta % 6)==1)
		{
			if(g[pos-1].best_next_bit > set_count && g[pos-1].best_next_bit<target-20)
				next_bit = g[pos-1].best_next_bit;
//			else
//				next_bit = -1;
			dup=true;
		}

		int bit_val = bit;
		
		if(!dup)// && g_kmax[0]>=g_kmax[1])
		{
			int sum0=0, sum1=0;
			for(int i=0; i<set_count; i++)
			{
				sum0+=gp.mask0.bit(i);
				sum1+=gp.mask.bit(i);
			}
			if(sum1>sum0)
				bit_val=1-bit;
		}
		
//		if(delta>=6 && (delta % 6)==1)

		if(gp.mask0.bit(next_bit) || gp.mask.bit(next_bit))
			bit_positions[pos]=-1;
		else
			bit_positions[pos]=next_bit;
		if(bit_positions[pos]>=0)
			bit_orders[pos]=(bit==bit_val);

		if(gp.mask0.bit(next_bit) && bit_val==1)
			goto next;
		if(gp.mask.bit(next_bit) && bit_val==0)
			goto next;

	
		{
			uint64_t t = __rdtsc();
			int tInt = (int)floor(((t-tStart)/g_fFreq)/900.);
			if(t>=tStart && tInt>tIntReportPrev && set_count>=len+10)
			{
				double progress = 0.0;
				for(int i=len; i<len+8; i++)
					if(gp.mask.bit(i))
						progress += pow(2.0, -(i-len)-1);
				printf("...thread %d sig 0x%llx_0x%llx: %.1f%% done, max %d, %lld ops\n", thr, 
					sig.n[0], sig.n[1], progress*100.0, max_len+2, op_count);
				tIntReportPrev = tInt;
			}
		}
		op_count++;
		if(g_max_ops!=0 && op_count > g_max_ops)
		{
			return 0;
		}
		max_len = max(max_len, set_count);

		if(set_count<stop_point && extend(gc, gp, bit_val, kmax, ext_target, next_bit, 2))
		{
			int first_unset=0;
			while((gc.mask0.bit(first_unset) || gc.mask.bit(first_unset)) && first_unset<=MAX_LEN*64)
				first_unset++;
			int last_bit=max(gc.mask.trailing_bit(), gc.mask0.trailing_bit());
			gc.dirty = gp.dirty || (last_bit!=set_count);
			gc.n = gp.n+1;

			if(gc.clique_counts[0][kmax[0]]>0  || gc.clique_counts[1][kmax[1]]>0)
				goto next;

			if(delta>=6 && !(delta % 6))
			{
				bigint2 m = make_pair(gc.mask0, gc.mask);
				for(int i=set_count; i<MAX_LEN*64; i++)
				{
					m.first.unset(i);
					m.second.unset(i);
				}
				
//				build_graph3(gmin[pos], m, kmax);
				if(gmin[pos-6].n==0)
				{
					printf("Error: null gmin[%d]\n", pos-6);
					exit(-1);
				}
//				if(g[pos].parent_cliques[0].size()>=4)
//					g[pos].merge();
			
				if(gc.dirty)
				{
					bigint mask0 = gc.mask0;
					bigint mask = gc.mask;
					graph& gprev = gmin[pos-6];
					if(!extend(gmin[pos], gprev, mask.bit(set_count-5), kmax, ext_target, set_count-5, 0))
						goto next;
					for(int i=1; i<5; i++)
						if(!extend(gmin[pos], mask.bit(set_count-5+i), kmax, ext_target, set_count-5+i))
							goto next;
					if(gmin[pos].parent_cliques[0].size()>=4)
						gmin[pos].merge();
				}
				else
				{
					graph_check(gc, kmax, "save");
					gmin[pos]=gc;
				}


				int last_bit=max(gc.mask.trailing_bit(), gc.mask0.trailing_bit());
				int last_known_bit=max(gmin[pos].mask.trailing_bit(), gmin[pos].mask0.trailing_bit());
				if(last_bit == last_known_bit)
				{
					gc = gmin[pos];
					gc.dirty = false;
				}
				else
				{
					if(!extend_from(gc, gmin[pos], make_pair(gc.mask0,gc.mask), kmax, ext_target))
						goto next;
					if(last_bit+1 == popcnt(gc.mask0)+popcnt(gc.mask))
						gc.dirty = false;
					else
						gc.dirty = true;
				}
				if(gc.clique_counts[0][kmax[0]]>0  || gc.clique_counts[1][kmax[1]]>0)
					goto next;

				int min_unset_bit = 35;
				if(kmax[0]<5 || kmax[1]<5)
					min_unset_bit=10;

				if((set_count > min_unset_bit) && (set_count+20 < target))
				{
					gc.n=pos;
					gc.check_first_unset_bit();
					int x = max_extensibility(gc, kmax, ext_target, false);
					if(x < target)
					{
	//					printf("*");
						goto next;
					}
					else
					{
						int last_bit=max(gc.mask.trailing_bit(), gc.mask0.trailing_bit());
						if(last_bit>pos)
							gc.dirty=true;
					}
				}
			}
			
			if(set_count == target)
			{
				int first_unset=0;
				while((gc.mask0.bit(first_unset) || gc.mask.bit(first_unset)) && first_unset<=MAX_LEN*64)
					first_unset++;
				if(first_unset<target)
				{
					printf("thread %d, sig 0x%llx: ERROR: set_count %d, but bit %d is not set\n", thr, sig.n[0], set_count, first_unset);
					exit(-1);
				}
				if(!validate_solution(gc, kmax, target, gmin, pos))
					goto next;
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
				populations[set_count]++;

				if(set_count>=target && pHits!=0)
				{
					if(!validate_solution(gc, kmax, target, gmin, pos))
						goto next;
					int first_unset=0;
					while((gc.mask0.bit(first_unset) || gc.mask.bit(first_unset)) && first_unset<=MAX_LEN*64)
						first_unset++;
					if(first_unset<set_count)
					{
						printf("ERROR: set_count %d, but bit %d is not set\n", set_count, first_unset);
						graph gg;
						bigint2 m = make_pair(gc.mask0, gc.mask);
						for(int i=0; i<4; i++)
							printf("%llx %llx %llx\n", gc.mask0.n[i], gc.mask.n[i], gc.mask0.n[i] | gc.mask.n[i]);
						bool ok = build_graph3(gg, m, kmax);
						gg.print_short();
						exit(-1);
					}

					bigint x = gc.mask;
					x.set_len(set_count+1);
					pHits->push_back(x);
					if(detect_only)
					{
						retval=1;
						goto done;
					}
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

	if(retval==2)
	{
		sig.print();
		printf("retval==2?\n");
		return 2;
	}
	if((retval==1) || (populations[target]!=popPrev))
	{
		omp_set_lock(&op_count_lock);
		op_counts[1].push_back(op_count);
		omp_unset_lock(&op_count_lock);
		return 1;
	}
	else
	{
		omp_set_lock(&op_count_lock);
		op_counts[0].push_back(op_count);
		omp_unset_lock(&op_count_lock);
		return 0;
	}
}

int depth_search(uint64_t sig, int len, int kmax[2], int target, int stop_point, bool detect_only, bool quiet, int* populations, graph* graph_pool, int thr, buffer<bigint>* pHits=0, uint64_t max_ops=0)
{
	bigint b;
	b.clear();
	b.n[0]=sig;
	b.set_len(len);
	return depth_search(b, kmax, target, stop_point, detect_only, quiet, populations, graph_pool, thr, pHits, max_ops);
}

typedef void hit_detect_callback_fun(uint64_t pos, bigint sig, buffer<bigint>* pHits, int max_len);
typedef void base_cover_callback_fun(bigint sig, uint64_t max_len);


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
	const bigint* bases;
	uint64_t* base_masks;
	uint64_t nbases;
	int nthr;
	uint64_t t0;
	hit_detect_callback_fun* pHitDetect;
	hit_detect_callback_fun* pHitMiss;
	base_cover_callback_fun* pBaseCover;
	uint64_t max_ops;
	string base_cover_name;
};


int undivided_depth_search_job(void* arg)
{
	depth_search_state* pState =(depth_search_state*)arg;
	uint64_t t0 = __rdtsc();
	uint64_t pos = pState->thr;
	int thr = pState->thr;
	buffer<bigint> hits;
	int local_populations[pop_acc_dim];
	buffer<bigint> cover_list;
	cover_list.alloc(1.5*pState->nbases/pState->nthr);
	while(true)
	{
		bigint val = pState->bases[pos];
		hits.clear();
		memset(local_populations, 0, sizeof(local_populations));
		int ok = depth_search(val, g_kmax, pState->target, pState->stop_size, pState->detect_only, true, local_populations, pState->graph_pool+thr*graph_pool_dim*2, thr, &hits, g_max_ops);

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
			 || (ok==0 && pState->pHitMiss)
			 || pState->len<49 && (pState->pBaseCover || (pState->base_cover_name.size()>0)))
		{
			omp_set_lock(pState->pCritSect);
	//		printf("0x%llx\t%d\t%d\t<%.3f>\n", val, ok?1:0, thr, (__rdtsc()-pState->t0)/g_fFreq);
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
			//omp_set_lock(&lock);
			//pos = *(pState->pNextPos);
			//(*(pState->pNextPos))++;

			if(pState->len<49)
			{
				if(pState->pBaseCover)
				{
					pState->pBaseCover(val, 0);
				}
				else if(pState->base_cover_name.size()>0)
				{
					FILE* f=fopen(pState->base_cover_name.c_str(), "a");
					fprintf(f, "0x%llx\n", val.n[0]);
					fclose(f);
				}
			}

			omp_unset_lock(pState->pCritSect);
		}
		if(pState->len>=49)
		{
			cover_list.push_back(val);
		}

		if(pos >= pState->nbases)
			break;
/*
		if(thr==0)
		{
			for(int i=0; i<graph_pool_dim; i++)
				printf("%d\t%d\t%d\n", i, pState->graph_pool[i].cliques[0].cap, pState->graph_pool[i].cliques[1].cap);
		}
*/
		//printf("thread %d grabbing %d:%d\n", thr, pos>>10, pos%1024);
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
		else if(pState->base_cover_name.size()>0)
		{
			omp_set_lock(pState->pCritSect);
			FILE* f=fopen(pState->base_cover_name.c_str(), "a");
			if(pState->len>64)
			{
				for(int i=0; i<cover_list.size(); i++)
					fprintf(f, "0x%llx\t0x%llx\n", cover_list[i].n[0], cover_list[i].n[1]);
			}
			else
			{
				for(int i=0; i<cover_list.size(); i++)
					fprintf(f, "0x%llx\n", cover_list[i].n[0]);
			}
			fclose(f);
			omp_unset_lock(pState->pCritSect);
		}
	}
	//printf("thread %d done\n", thr);
	//fflush(stdout);
	uint64_t t1 = __rdtsc();
	pState->cputime = t1-t0;
	return 0;
}


struct JobArray
{
	omp_lock_t critSect;
	depth_search_state jobs[nthreads];
	int* populations;
	uint64_t sum_populations[pop_acc_dim];
	graph* graph_pool;
	uint64_t cputime;
	uint64_t realtime;
	string base_cover_name;

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

	void submit(const bigint* p, size_t n, int len, int target, int stop_size, bool detect_only, int fractions=1, 
		base_cover_callback_fun* pBaseCover=0, hit_detect_callback_fun* pHitDetect=0, 
		hit_detect_callback_fun* pHitFail=0, uint64_t max_ops=0);
	void print_populations();
	void print_cputime();
};

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

void JobArray::submit(const bigint* p, size_t n, int len, int target, int stop_size, bool detect_only, int fractions, base_cover_callback_fun* pBaseCover, hit_detect_callback_fun* pHitDetect, hit_detect_callback_fun* pHitFail, uint64_t max_ops)
{
	memset(populations, 0, nthreads*pop_acc_dim*sizeof(int));
	memset(sum_populations, 0, sizeof(sum_populations));
	prealloc_graph_pool(graph_pool, len);
	uint64_t* base_masks=new uint64_t[n];
	for(int i=0; i<n; i++)
		base_masks[i]=0;
	unsigned int next_pos=nthreads;
	fractions=0;
	uint64_t t0 = __rdtsc();
	for(int thr=0; thr<nthreads; thr++)
	{
		jobs[thr].cputime = 0;
		if(fractions==0 && thr>=n)
			continue;
		jobs[thr].bases = p;
		jobs[thr].nbases = n;
		jobs[thr].base_masks = base_masks;
		jobs[thr].detect_only = detect_only;
		jobs[thr].graph_pool = graph_pool;
		jobs[thr].len = len;
		jobs[thr].fractions = fractions;
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
		jobs[thr].max_ops=max_ops;
		jobs[thr].base_cover_name=base_cover_name;
		//jobs[thr].max_ops=max_ops;
//		if(fractions==0)
	}

	int nthr = min(nthreads, (int)n);

#pragma omp parallel for
	for(int thr=0; thr<nthr; thr++)
		undivided_depth_search_job((void*)(jobs+thr));

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


char hit_detect_name[128];

string hit_detect_name2;

uint64_t mc_hit_count=0;

void mc_hit_detect(uint64_t pos, bigint sig, buffer<bigint>* p, int max_len)
{
	mc_hit_count++;
	FILE* f=fopen(hit_detect_name, "a");
	if(sig.len>64)
		fprintf(f, "0x%llx\t0x%llx\n", sig.n[0], sig.n[1]);
	else
		fprintf(f, "0x%llx\n", sig.n[0]);
	fclose(f);
	if(hit_detect_name2.size()>0)
	{
		FILE* f=fopen(hit_detect_name2.c_str(), "a");
		if(sig.len>64)
			fprintf(f, "0x%llx\t0x%llx\n", sig.n[0], sig.n[1]);
		else
			fprintf(f, "0x%llx\n", sig.n[0]);
		fclose(f);
	}
}


string neighbor_hit_detect_name;
string neighbor_hit_detect_name2;
string neighbor_base_cover_name;

uint64_t neighbor_hit_count = 0;

static FILE* neighbor_hit_file=0;
string neighbor_hit_file_name="";


void neighbor_hit_detect(uint64_t pos, bigint sig, buffer<bigint>* p, int max_len)
{
	neighbor_hit_count++;
	if(neighbor_hit_file_name!=neighbor_hit_detect_name || neighbor_hit_file==0)
	{
		if(neighbor_hit_file!=0)
			fclose(neighbor_hit_file);
		neighbor_hit_file=fopen(neighbor_hit_detect_name.c_str(), "a");
		neighbor_hit_file_name=neighbor_hit_detect_name;
	}
	if(sig.len>64)
		fprintf(neighbor_hit_file, "0x%llx\t0x%llx\n", sig.n[0], sig.n[1]);
	else
		fprintf(neighbor_hit_file, "0x%llx\n", sig.n[0]);
	fflush(neighbor_hit_file);
	if(neighbor_hit_detect_name2.size()>0)
	{
		FILE* f=fopen(neighbor_hit_detect_name2.c_str(), "a");
		if(sig.len>64)
			fprintf(f, "0x%llx\t0x%llx\n", sig.n[0], sig.n[1]);
		else
			fprintf(f, "0x%llx\n", sig.n[0]);
		fclose(f);
	}
}

void neighbor_base_cover(bigint sig, uint64_t max_len)
{
	FILE* f=fopen(neighbor_base_cover_name.c_str(), "a");
	if(sig.len<=64)
		fprintf(f, "0x%llx\n", sig.n[0]);
	else
		fprintf(f, "0x%llx\t0x%llx\n", sig.n[0], sig.n[1]);
	fclose(f);
}


vector<pair<uint64_t,uint64_t> > g_mc_exclusions;

string format_name(int len, int target)
{
	char temp[128];
	sprintf(temp, PATH "%d%d-%d-%d.txt", g_kmax[0]+2, g_kmax[1]+2, len+1, target);
	return string(temp);
}

string format_cover_name(int len, int target)
{
	char temp[128];
	sprintf(temp, PATH "cover%d%d-%d-%d.txt", g_kmax[0]+2, g_kmax[1]+2, len+1, target);
	return string(temp);
}

string format_review_name(int len, int target)
{
	char temp[128];
	sprintf(temp, PATH "review%d%d-%d-%d.txt", g_kmax[0]+2, g_kmax[1]+2, len+1, target);
	return string(temp);
}


void readfile(const char* fn, vector<bigint>& sigs, int len)
{
	FILE* f=fopen(fn, "r");
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
				bigint b;
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
				bigint b;
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

void writefile(const char* fn, vector<bigint>& sigs, int len)
{
	FILE* f=fopen(fn, "w");
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

void readfile(string fn, vector<bigint>& sigs, int len)
{
	readfile(fn.c_str(), sigs, len);
}

size_t exclude(vector<bigint>& vsig, vector<bigint>& vsigExcl, int len=0)
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
	vector<bigint> vsig2[16];
#pragma omp parallel for
	for(int64_t thr=0; thr<16; thr++)
	{
		size_t low=vsig.size()*thr/16;
		size_t high=vsig.size()*(thr+1)/16;
		vsig2[thr].reserve(high-low);
		for(size_t i=low; i<high; i++)
		{
			bigint q = vsig[i];
			if(len!=0)
				q.unset_from(len);
			bool dup=binary_search(vsigExcl.begin(), vsigExcl.end(), q);
			if(!dup)
				vsig2[thr].push_back(vsig[i]);
		}
	}
//	for(int64_t thr=0; thr<16; thr++)
//		printf("%d ", vsig2[thr].size());
//	printf("\n");

	uint64_t total_size=0;
	for(int64_t thr=0; thr<16; thr++)
		total_size+=vsig2[thr].size();
	size_t retval = vsig.size()-total_size;
	vsig.resize(total_size);
	total_size=0;
	for(int64_t thr=0; thr<16; thr++)
	{
		memcpy(&vsig[total_size], &vsig2[thr][0], vsig2[thr].size()*sizeof(bigint));
		total_size+=vsig2[thr].size();
	}
	return retval;
}

bool mc_detect_only=true;

bigint random_bigint(int len)
{
	bigint x;
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

int mc(int len, int k0, int k1, int target, uint64_t nbases)
{
	if(MAX_LEN==2 && target>129)
	{
		printf("Recompile with MAX_LEN=4\n");
		exit(-1);
	}
	
	if(MAX_LEN==4 && target<129)
	{
		printf("Warning: recommended to recompile with MAX_LEN=2\n");
	}
	vector<bigint> vsig;
	srand(__rdtsc());
	vsig.resize(nbases);
	g_kmax[0]=k0;
	g_kmax[1]=k1;

	string s = format_name(len, target);
	strcpy(hit_detect_name, s.c_str());
	hit_detect_name2 = format_cover_name(len, target);

#pragma omp parallel for
	for(int64_t thr=0; thr<24; thr++)
	{
		int64_t low = vsig.size()*thr/24;
		int64_t high = vsig.size()*(thr+1)/24;
		if(len<=64)
		{
			for(int64_t i=low; i<high; i++)
			{
				vsig[i].clear();
				vsig[i].set_len(len);
				vsig[i].n[0]=rand64() & ((One<<len)-2);
			}
		}
		else
		{
			for(int64_t i=low; i<high; i++)
			{
				vsig[i].clear();
				vsig[i].set_len(len);
				vsig[i].n[0]=rand64() & (-2ll);
				vsig[i].n[1]=rand64() & ((One<<(len-64))-2);
			}
		}
	}
	/*
	uint64_t i;
	vector<bigint> vsigExcl, vsigExclParent;
	readfile(hit_detect_name, vsigExcl, len);
	size_t sizePrev = vsigExcl.size();
	rem_dup(vsigExcl);
	if(sizePrev != vsigExcl.size())
	{
		printf("WARNING: %d duplicates in the hit file\n", sizePrev-vsigExcl.size());
	}
	string exclParent = format_cover_name(len-5, target);
	readfile(exclParent, vsigExclParent, len);

	for(i=0; i<g_mc_exclusions.size(); i++)
	{
//		string s = format_name(g_mc_exclusions[i].first, g_mc_exclusions[i].second);
//		readfile1(s.c_str(), vsigExcl);
		string s = format_cover_name(g_mc_exclusions[i].first, g_mc_exclusions[i].second);
		if(g_mc_exclusions[i].first == len)
			readfile(s.c_str(), vsigExcl, g_mc_exclusions[i].first);
		else
			readfile(s.c_str(), vsigExclParent, g_mc_exclusions[i].first);
	}
	
	size_t nEx = exclude(vsig, vsigExcl);
	nEx+=exclude(vsig, vsigExclParent, len-5);
	printf("%d previously covered\n", nEx);
	nEx = rem_dup(vsig);
	printf("%d duplicate\n", nEx);
	*/
	JobArray array;
	bool detect_only=mc_detect_only;

	int stop_size=detect_only ? target+1 : MAX_LEN*64;
	mc_hit_count=0;
	array.submit(&vsig[0], vsig.size(), len, target, stop_size, detect_only, 0, 0, mc_hit_detect);
	array.print_populations();
	array.print_cputime();
	printf("%d -> %d\n", vsig.size(), mc_hit_count);
	FILE* ff=fopen(PATH "neighbor-search-log.txt", "a");
	fprintf(ff, "%d\t%d\t%.3f\t(mc)\n", vsig.size(), mc_hit_count, array.realtime/g_fFreq);
	fclose(ff);
	return mc_hit_count;
}

vector<pair<uint64_t,uint64_t> > g_ns_exclusions;
vector<vector<bigint> > g_ns_exclusion_vectors;
vector<pair<uint64_t,uint64_t> > g_ns_exclusion_vector_metadata;

vector<bigint> g_vNeighborKnownCover;
int g_iKnownCoverLen=0, g_iKnownCoverTarget=0;

void clusterize(vector<vector<bigint> >& vv, const vector<bigint>& b);

bool ns_detect_only=true;

int g_last_cluster_count=0;

vector<bigint> g_ns2_hits;
void neighbor_hit_detect_v2(uint64_t pos, bigint sig, buffer<bigint>* p, int max_len)
{
	g_ns2_hits.push_back(sig);
}

vector<bigint> neighbor_depth_search_v2(int len, int k0, int k1, int target, vector<bigint>& vsig0, vector<bigint>& vexcl)
{
	vector<bigint> vsig;
	int nsigs0 = vsig0.size();
	vsig.resize((len-1)*nsigs0);
	for(int i=0; i<nsigs0; i++)
		for(int j=0; j<len-1; j++)
		{
			vsig[i*(len-1)+j]=vsig0[i];
			if(vsig0[i].bit(j+1))
				vsig[i*(len-1)+j].unset(j+1);
			else
				vsig[i*(len-1)+j].set(j+1);

		}

	size_t nEx = rem_dup(vsig);
	exclude(vsig, vexcl);
	append(vexcl, vsig);
	g_ns2_hits.clear();
	JobArray array;
	array.base_cover_name.clear();
	array.submit(&vsig[0], vsig.size(), len, target, target+1, true, 0, 0, neighbor_hit_detect_v2);
//	printf("%d -> %d\n", vsig.size(), g_ns2_hits.size());
	return g_ns2_hits;
}

vector<bigint> narrow_v2(int len, int k0, int k1, int target, vector<bigint>& vsig)
{
	size_t nEx = rem_dup(vsig);
	g_ns2_hits.clear();
	JobArray array;
	array.base_cover_name.clear();
	array.submit(&vsig[0], vsig.size(), len, target, target+1, true, 0, 0, neighbor_hit_detect_v2);
	return g_ns2_hits;
}

bool g_offby2=false;

int neighbor_depth_search(int len, int k0, int k1, int target, int max_cluster=INT_MAX)
{
	if(MAX_LEN==2 && target>129)
	{
		printf("Recompile with MAX_LEN=4\n");
		exit(-1);
	}
	if(MAX_LEN==4 && target<129)
	{
		printf("Warning: recommended to recompile with MAX_LEN=2\n");
//		exit(-1);
	}

	g_kmax[0]=k0;
	g_kmax[1]=k1;
	vector<bigint> vsig, vsig0, vsigExcl;
	neighbor_hit_detect_name=format_name(len, target);
	neighbor_base_cover_name=format_cover_name(len, target);
	neighbor_hit_detect_name2.clear();
	readfile(neighbor_hit_detect_name.c_str(), vsig0, len);

	int nCl=-1;

	size_t sizePrev = vsig0.size();
	rem_dup(vsig0);
	if(sizePrev != vsig0.size())
	{
		printf("WARNING: %d duplicates in the hit file\n", sizePrev-vsig0.size());
	}

	if(max_cluster != INT_MAX)
	{
		vector<vector<bigint> > vCl;
		clusterize(vCl, vsig0);
		printf("%d elements in %d clusters\n", vsig0.size(), vCl.size());
		nCl = g_last_cluster_count=vCl.size();
		vsig0.clear();
		for(int i=0; i<vCl.size(); i++)
		{
			if(vCl[i].size()<=max_cluster)
			{
				for(int j=0; j<vCl[i].size(); j++)
				{
					vsig0.push_back(vCl[i][j]);
				}
			}
		}
		printf("max cluster %d: keeping %d elements\n", max_cluster, vsig0.size());
	}
	else
	{
		printf("%d elements\n", vsig0.size());
	}

	if(len!=g_iKnownCoverLen || target!=g_iKnownCoverTarget)
	{
		g_iKnownCoverLen=len;
		g_iKnownCoverTarget=target;
		g_vNeighborKnownCover.clear();
		readfile(neighbor_base_cover_name.c_str(), vsigExcl, len);
		g_vNeighborKnownCover = vsigExcl;
	}
	else
	{
		vsigExcl = g_vNeighborKnownCover;
	}



	printf("%d excluded sigs\n", vsigExcl.size());

	int nsigs0 = vsig0.size();
	if(g_offby2)
	{
		vsig.resize((len-1)*(len-1)*nsigs0);
		for(int i=0; i<nsigs0; i++)
			for(int j=0; j<len-1; j++)
				for(int k=0; k<len-1; k++)
			{
				bigint& x = vsig[i*(len-1)*(len-1)+j*(len-1)+k];
				x=vsig0[i];
				if(vsig0[i].bit(j+1))
					x.unset(j+1);
				else
					x.set(j+1);
				if(vsig0[i].bit(k+1))
					x.unset(k+1);
				else
					x.set(k+1);
			}
	}
	else
	{
		vsig.resize((len-1)*nsigs0);
		for(int i=0; i<nsigs0; i++)
			for(int j=0; j<len-1; j++)
			{
				vsig[i*(len-1)+j]=vsig0[i];
				if(vsig0[i].bit(j+1))
					vsig[i*(len-1)+j].unset(j+1);
				else
					vsig[i*(len-1)+j].set(j+1);
			}
	}
	size_t nEx = rem_dup(vsig);
	printf("%d duplicate\n", nEx);
	if(g_offby2)
	{
		for(int i=0; i<vsig.size(); i++)
		{
			uint64_t pos1 = rand64() % vsig.size();
			uint64_t pos2 = rand64() % (vsig.size()-1);
			if(pos2==pos1)
				pos2++;
			swap(vsig[pos1], vsig[pos2]);
		}
	}
	printf("%d candidates\n", vsig.size());
	nEx = exclude(vsig, vsigExcl);
	printf("%d previously covered\n", nEx);
	string exclParent = format_cover_name(len-5, target);
	vector<bigint> vsigExclParent;
	readfile(exclParent.c_str(), vsigExclParent, len-5);
	nEx = exclude(vsig, vsigExclParent, len-5);
	printf("%d previously covered at len=%d\n", nEx, len-5);

//	if(g_ns_exclusion_vectors.size()==0)
	g_ns_exclusion_vectors.resize(g_ns_exclusions.size());
	g_ns_exclusion_vector_metadata.resize(g_ns_exclusions.size());
	for(int k=0; k<g_ns_exclusions.size(); k++)
	{
		if(g_ns_exclusion_vector_metadata[k]==g_ns_exclusions[k])
			continue;
		g_ns_exclusion_vectors[k].clear();
		string fn = format_cover_name(g_ns_exclusions[k].first, g_ns_exclusions[k].second);
		readfile(fn.c_str(), g_ns_exclusion_vectors[k], g_ns_exclusions[k].first);
		g_ns_exclusion_vector_metadata[k]=g_ns_exclusions[k];
	}

	for(int k=0; k<g_ns_exclusions.size(); k++)
	{
		if(g_ns_exclusion_vectors[k].size()==0)
			continue;
		string fn = format_cover_name(g_ns_exclusions[k].first, g_ns_exclusions[k].second);
		int nEx = exclude(vsig, g_ns_exclusion_vectors[k], g_ns_exclusions[k].first);
		printf("%s: %d exclusions\n", fn.c_str(), nEx);
	}

	printf("%d remaining\n", vsig.size());

	append(g_vNeighborKnownCover, vsig);
	neighbor_hit_count=0;
	JobArray array;
	int stop_size=ns_detect_only ? target+1 : MAX_LEN*64;
	array.base_cover_name=neighbor_base_cover_name;
	if(g_offby2)
	{
		neighbor_hit_detect_name2=neighbor_base_cover_name;
		neighbor_base_cover_name.clear();
		array.base_cover_name.clear();
	}
	array.submit(&vsig[0], vsig.size(), len, target, stop_size, ns_detect_only, 0, 0, neighbor_hit_detect);
	if(neighbor_hit_file!=0)
		fclose(neighbor_hit_file);
	neighbor_hit_file=0;
	array.print_populations();
	array.print_cputime();
	printf("%d -> %d\n", vsig.size(), neighbor_hit_count);
	FILE* ff=fopen(PATH "neighbor-search-log.txt", "a");
	fprintf(ff, "%d\t%d\t%.3f\t%d\n", vsig.size(), neighbor_hit_count, array.realtime/g_fFreq, nCl);
	//fprintf(ff, "%d\t%d\t%.3f\n", vsig.size(), neighbor_hit_count, array.realtime/g_fFreq);
	fclose(ff);
	return neighbor_hit_count;
}

string review_filename="review-xl-df.txt";

void xl_hit_fail(uint64_t pos, bigint sig, buffer<bigint>* p, int max_len)
{
	printf("0x%llx: miss\n", sig.n[0]);
}

void xl_hit_detect(uint64_t pos, bigint sig, buffer<bigint>* p, int max_len)
{
	int counts[MAX_LEN*64], sym[MAX_LEN*64];
	memset(counts, 0, sizeof(counts));
	memset(sym, 0, sizeof(sym));

	int i;
	for(i=0; i<p->size(); i++)
	{
		bigint& x = (*p)[i];
		int d = x.len+1;
		counts[d]++;
	
		if(d>=130)
		{
			bigint inv;
			invert(inv, x);
			if(bit_equal(x, inv))
			{
				sym[d]++;
				if(d>=130)
				{
					FILE* f=fopen(PATH "symmetric.txt", "a");
					fprintf(f, "%d%d\t%d\t%llx\t%llx\t%llx\t%llx\n", g_kmax[0], g_kmax[1], d, x.n[0], x.n[1], x.n[2], x.n[3]);
					fclose(f);
				}
			}
		}
	}
	max_len=0;
	int max_sym=0;
	for(i=0; i<MAX_LEN*64; i++)
		if(counts[i]>0)
		{
			max_len=i;
			max_sym=sym[i];
		}

	FILE* f=fopen(review_filename.c_str(), "a");
	fprintf(f, "0x%llx\t%d\t%d\t", sig.n[0], max_len, max_sym);
	for(int i=110; i<=max_len; i++)
		fprintf(f, "%d\t", counts[i]);
	fprintf(f, "\n");
	fclose(f);
}

void narrow_hit_detect(uint64_t pos, bigint sig, buffer<bigint>* p, int size)
{
	FILE* f=fopen(review_filename.c_str(), "a");
	if(sig.len<=64)
		fprintf(f, "0x%llx\n", sig.n[0]);
	else
		fprintf(f, "0x%llx\t0x%llx\n", sig.n[0], sig.n[1]);
	fclose(f);
}

int narrow(int len, int targetPrev, int target, bool detect_only=true)
{
	printf("narrow(%d, %d, %d)\n", len, targetPrev, target);
	if(target > MAX_LEN*64)
	{
		printf("Error: recompile with MAX_LEN=4 or higher\n");
		exit(-1);
	}

	string fnIn = format_name(len, targetPrev);
	string fnOut = format_name(len, target);
	string fnInCover = format_cover_name(len, targetPrev);
	string fnInCoverParent = format_cover_name(len-5, targetPrev);

	string fnOutCover = format_cover_name(len, target);
	string fnOutCoverParent = format_cover_name(len-5, target);

	uint64_t i=0;
	vector<bigint> vsig;
	readfile(fnIn, vsig, len);
	vector<bigint> vsigNextShort, vsigNextParent;
	
	readfile(fnOutCover, vsigNextShort, len);
	readfile(fnOutCoverParent, vsigNextParent, len-5);

	vector<bigint> vsigInCover, vsigInCoverParent;
	readfile(fnInCover, vsigInCover, len);
	readfile(fnInCoverParent, vsigInCoverParent, len-5);

	vector<bigint> check;
	readfile(fnOut, vsigNextShort, len);
	exclude(check, vsigNextShort);
	exclude(check, vsigNextParent, len-5);
	if(!check.empty())
	{
		printf("WARNING: %d entries in target list are not in target cover list\n", check.size());
	}

	rem_dup(vsigNextShort);
	int dup = rem_dup(vsig);
	printf("%d duplicate\n", dup);
	int nEx = exclude(vsig, vsigNextShort);
	printf("%d excluded\n", nEx);
	nEx = exclude(vsig, vsigNextParent, len-5);
	printf("%d excluded by parent cover\n", nEx);
	nEx = exclude(vsig, check);
	printf("%d excluded by target list\n", nEx);
	printf("%d remaining\n", vsig.size());

	if(vsig.size()==0)
		return 0;

	review_filename=fnOut;
	neighbor_hit_detect_name2.clear();

	JobArray array;
	int stop_size=detect_only ? target+1 : MAX_LEN*64;
	array.submit(&vsig[0], vsig.size(), len, target, stop_size, detect_only, 0, 0, narrow_hit_detect, 0);
	array.print_cputime();
	printf("%d -> %d\n", vsig.size(), array.sum_populations[target]);
	FILE* ff=fopen(PATH "neighbor-search-log.txt", "a");
	fprintf(ff, "%d\t%d\t%.3f\t(oldnarrow)\n", vsig.size(), array.sum_populations[target], array.realtime/g_fFreq);
	fclose(ff);

	append(vsigNextShort, vsigInCover);
	append(vsigNextParent, vsigInCoverParent);
	rem_dup(vsigNextShort);
	rem_dup(vsigNextParent);
	writefile(fnOutCover.c_str(), vsigNextShort, len);
	writefile(fnOutCoverParent.c_str(), vsigInCoverParent, len);

	return array.sum_populations[target];
}

void review(int len, int k0, int k1, int target)
{
	g_kmax[0]=k0;
	g_kmax[1]=k1;
#if MAX_LEN < 4
	printf("Error: recompile with MAX_LEN=4 or higher\n");
	exit(-1);
#endif
	string fnIn = format_name(len, target);
	string fnOut = format_review_name(len, target);
	review_filename=fnOut;
	uint64_t i=0;
	vector<bigint> vsig;
	readfile(fnIn.c_str(), vsig, len);
	vector<bigint> vsigExcl;
	readfile(fnOut.c_str(),vsigExcl, len);
	exclude(vsig, vsigExcl);
	printf("%d not covered\n", vsig.size());
	if(vsig.size()==0)
		return;
	int nEx=rem_dup(vsig);
	if(nEx!=0)
		printf("%d duplicate\n", nEx);

	JobArray array;
	int stop_size=MAX_LEN*64;
	array.submit(&vsig[0], vsig.size(), len, target, stop_size, false, 0, 0, xl_hit_detect, xl_hit_fail);
	array.print_populations();
	array.print_cputime();
}

int dist(bigint a, bigint b)
{
	return __popcnt64(a.n[0]^b.n[0]) + __popcnt64(a.n[1]^b.n[1]);
}

int dist(bigint a, const vector<bigint>& vb)
{
	int min_dist=100;
	for(int i=0; i<vb.size(); i++)
	{
		if(dist(a,vb[i])==1)
			return 1;
	}
	return min_dist;
}

int dist2(bigint a, const vector<bigint>& vb)
{
	int min_dist=100;
	for(int i=0; i<vb.size(); i++)
	{
		min_dist=min(min_dist, dist(a,vb[i]));
	}
	return min_dist;
}

int dist(const vector<bigint>& va, const vector<bigint>& vb)
{
	int min_dist=100;
	for(int i=0; i<va.size(); i++)
		min_dist=min(min_dist, dist(va[i],vb));
	return min_dist;
}

void clusterize_part2(vector<vector<bigint> >& vv)
{
	{
		vector<int> pairwise_dist;
		pairwise_dist.resize(vv.size()*vv.size());
#pragma omp parallel for
		for(int thr=0; thr<4; thr++)
		{
			int low, high;
			if(thr==0)
			{
				low=0;
				high=vv.size()/8;
			}
			else if(thr==1)
			{
				low=vv.size()/8;
				high=vv.size()/4;
			}
			else if(thr==2)
			{
				low=vv.size()/4;
				high=vv.size()/2;
			}
			else
			{
				low=vv.size()/2;
				high=vv.size();
			}


			for(int i=low; i<high; i++)
			{
				for(int j=i+1; j<vv.size(); j++)
				{
					pairwise_dist[i+j*vv.size()] = dist(vv[i], vv[j]);
					pairwise_dist[j+i*vv.size()] = pairwise_dist[i+j*vv.size()] ;
				}
			}
		}
		int change=0;
		vector<vector<bigint> > vv2;
		vector<bool> covered;
		covered.resize(vv.size());
		for(int i=0; i<vv.size(); i++)
			covered[i]=false;
		for(int i=0; i<vv.size(); i++)
		{
			if(covered[i])
				continue;
			vector<uint64_t> neighbors;
			for(int j=i+1; j<vv.size(); j++)
			{
				if(pairwise_dist[i+j*vv.size()]<=1)
					neighbors.push_back(j);
			}
			while(true)
			{
				int changes=0;
				for(int j=0; j<neighbors.size(); j++)
				{
					int pos = neighbors[j];
					for(int k=i+1; k<vv.size(); k++)
					{
						if(k==pos)
							continue;
						if(pairwise_dist[pos+k*vv.size()]<=1)
						{
							bool dup=false;
							for(int m=0; m<neighbors.size(); m++)
								if(neighbors[m]==k)
									dup=true;
							if(!dup)
							{
								neighbors.push_back(k);
								changes++;
							}
						}
					}
				}
				if(changes==0)
					break;
			}
			vector<bigint> vnew = vv[i];
			for(int j=0; j<neighbors.size(); j++)
			{
				covered[neighbors[j]]=true;
				append(vnew,vv[neighbors[j]]);
			}
			vv2.push_back(vnew);
		}
		vv.swap(vv2);
	}
}

void clusterize(vector<vector<bigint> >& vv, const bigint* pb, size_t sz)
{
	if(sz==0)
		return;
	for(int i=0; i<sz; i++)
	{
		int min_dist=100;
		bool bad=false;
		int j;
		for(j=0; j<vv.size(); j++)
		{
			int d = dist(pb[i],vv[j]);
			min_dist = min(min_dist,d);
			if(min_dist<=1)
			{
				if(min_dist!=0)
					vv[j].push_back(pb[i]);
				bad=true;
				break;
			}
		}
		if(!bad)
		{
			vector<bigint> v;
			v.push_back(pb[i]);
			vv.push_back(v);
		}
	}
	//printf("%d entries, %d initial clusters\n", sz, vv.size());
	clusterize_part2(vv);
}

void clusterize(vector<vector<bigint> >& vv, const vector<bigint>& b)
{
	int max_dist=1;
	int dup=0;
	vector<bigint> b2 = b;
	int nEx = rem_dup(b2);
	if(nEx!=0)
		printf("WARNING: %d duplicates in hit file\n", nEx);
	clusterize(vv, &b[0], b.size());
}

void cover_rebuild(int len, int k0, int k1, int target, const char* secondary=0)
{
	g_kmax[0]=k0;
	g_kmax[1]=k1;
	string s=format_name(len, target);
	string sCover=format_cover_name(len, target);
	vector<bigint> v, vsig;
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


int lengthen(const char* fn, int k0, int k1, int len, int lenNext, int target)
{
	g_kmax[0]=k0;
	g_kmax[1]=k1;

	if(target > MAX_LEN*64)
	{
		printf("Recompile with higher MAX_LEN\n");
		exit(-1);
	}

	vector<bigint> vsig, vsig0, vsigExcl, vsigExclParent;
	neighbor_hit_detect_name=format_name(lenNext, target);
	neighbor_base_cover_name=format_cover_name(lenNext, target);
//	string parent_cover_name=format_cover_name(lenNext-5, target);
	vector<bigint> vx;
	readfile(neighbor_hit_detect_name, vx, lenNext);
	readfile(neighbor_base_cover_name, vx, lenNext);
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
//	nEx = exclude(vsig, vsigExclParent, lenNext-5);
//	printf("%d previously covered by parent\n", nEx);

	
	//string srcCover = format_cover_name(len, target);
	//readfile(srcCover.c_str(), vsigExcl, len);

	vector<bigint> vsigExclLong;
	/*
	mul = 1<<(lenNext-len-5);
	vsigExclLong.resize(vsigExcl.size()*mul);
	for(uint64_t i=0; i<vsigExcl.size(); i++)
		for(uint64_t j=0; j<mul; j++)
		{
			vsigExclLong[i*mul+j] = vsigExcl[i];
			vsigExclLong[i*mul+j].set_len(lenNext-5);
			for(int k=len; k<lenNext-5; k++)
				if(j & (1i64<<(k-len)))
				{
					vsigExclLong[i*mul+j].set(k);
				}
		}
	nEx = exclude(vsigExclLong, vsigExclParent, lenNext-5);
		*/
//	vsigExclLong  = vsig;
//	nEx = exclude(vsigExclLong, vx);
//	exclude(vsigExclLong, vsig);
//	vector<uint64_t> dstCover;
//	readfile1(neighbor_base_cover_name.c_str(), dstCover);
	/*
	printf("%d parent cover points already known\n", nEx);
	printf("%d new\n", vsigExclLong.size());
	FILE* f=fopen(neighbor_base_cover_name.c_str(), "a");
	for(uint64_t i=0; i<vx.size(); i++)
		fprintf(f, "0x%I64x\n", vx[i].n[0]);
	fclose(f);
	*/
	/*
	FILE* f=fopen(parent_cover_name.c_str(), "a");
	for(uint64_t i=0; i<vsigExclLong.size(); i++)
		fprintf(f, "0x%I64x\n", vsigExclLong[i].n[0]);
	fclose(f);
	*/
	printf("%d candidates\n", vsig.size());

	int stop_size=target+1;
	JobArray array;
	array.submit(&vsig[0], vsig.size(), lenNext, target, stop_size, true, 0, 0, neighbor_hit_detect);
	fclose(neighbor_hit_file);
	neighbor_hit_file=0;
	printf("%d -> %d\n", vsig.size(), array.sum_populations[target]);
	FILE* ff=fopen(PATH "neighbor-search-log.txt", "a");
	fprintf(ff, "%d\t%d\t%.3f\t(lengthen)\n", vsig.size(), array.sum_populations[target], array.realtime/g_fFreq);
	fclose(ff);
	return array.sum_populations[target];
}

void do77_part1()
{
	while(true)
	{
		g_ns_exclusions.clear();

		_unlink(PATH "77-55-110.txt");
		_unlink(PATH "77-45-110.txt");
		_unlink(PATH "77-55-115.txt");
		_unlink(PATH "77-55-120.txt");
		_unlink(PATH "77-55-125.txt");

		_unlink(PATH "cover77-55-110.txt");
		_unlink(PATH "cover77-45-110.txt");
		_unlink(PATH "cover77-55-115.txt");
		_unlink(PATH "cover77-55-120.txt");
		_unlink(PATH "cover77-55-125.txt");

		cover_rebuild(54,5,5,128);

		FILE* f;

		g_max_ops=1e5;
		uint64_t t0, t1;
		uint64_t time_zones[7];
		int hit_counts[7];
		memset(hit_counts, 0, sizeof(hit_counts));
		t0=__rdtsc();
		hit_counts[0]=0;
	//	for(int pass=0; pass<10; pass++)
		hit_counts[0]+=mc(44, 5, 5, 110, 1e7);
		lengthen(PATH "77-45-110.txt", 5, 5, 44, 54, 110);
		vector<bigint> v;
		readfile(PATH "77-45-110.txt", v, 44);
		rem_dup(v);
		vector<bigint> v2;

		uint64_t mul = 32;
		v2.resize(mul*v.size());
		for(uint64_t i=0; i<v.size(); i++)
			for(uint64_t j=0; j<mul; j++)
			{
				v2[i*mul+j]=v[i];
				v2[i*mul+j].set_len(49);
				for(int k=44; k<49; k++)
				if(j & (One<<(k-44)))
				{
					v2[i*mul+j].set(k);
				}
			}
		writefile(PATH "cover77-50-110.txt", v2, 49);
	
		g_max_ops=1e4;
		t1=__rdtsc();
		time_zones[0]=t1-t0;
		g_ns_exclusions.push_back(make_pair(54,128));

		t0=__rdtsc();
		hit_counts[1]+=neighbor_depth_search(54,5,5,110);
		hit_counts[1]+=neighbor_depth_search(54,5,5,110);
		hit_counts[1]+=neighbor_depth_search(54,5,5,110);
		hit_counts[1]+=neighbor_depth_search(54,5,5,110);
		hit_counts[2]=narrow(54, 110, 120);
		t1=__rdtsc();
		time_zones[1]=t1-t0;

		t0=__rdtsc();
		hit_counts[2]+=neighbor_depth_search(54,5,5,120);
		hit_counts[2]+=neighbor_depth_search(54,5,5,120);
		hit_counts[2]+=neighbor_depth_search(54,5,5,120);
		hit_counts[2]+=neighbor_depth_search(54,5,5,120);
		hit_counts[3]=narrow(54, 120, 125);
		t1=__rdtsc();
		time_zones[2]=t1-t0;
		t0=__rdtsc();
		for(int pass=0; pass<12; pass++)
		{
			int n = neighbor_depth_search(54,5,5,125);
			hit_counts[3]+=n;
			if(n==0)
				break;
		}
		hit_counts[4]=narrow(54, 125, 128);
		t1=__rdtsc();
		time_zones[3]=t1-t0;
		t0=__rdtsc();
		while(true)
		{
			int n = neighbor_depth_search(54,5,5,128);
			hit_counts[4]+=n;
			if(n==0)
				break;
		}
		t1=__rdtsc();
		time_zones[4]=t1-t0;
		f=fopen(PATH "log77.txt", "a");
		fprintf(f, "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%d\t%d\t%d\n",
			time_zones[0]/g_fFreq, time_zones[1]/g_fFreq, 
			time_zones[2]/g_fFreq, time_zones[3]/g_fFreq, time_zones[4]/g_fFreq, 
			hit_counts[0], hit_counts[1], hit_counts[2], hit_counts[3], hit_counts[4]);
		fclose(f);
	}
}

void test_inv(int k0, int k1, int len, int target)
{
	g_kmax[0]=k0;
	g_kmax[1]=k1;
	g_max_ops=1e5;
	while(neighbor_depth_search(len,k0,k1,target)!=0);
	return;
	vector<bigint> v;

	g_kmax[0]=k1;
	g_kmax[1]=k0;
	readfile(format_name(len,target), v, len);
	g_kmax[0]=k0;
	g_kmax[1]=k1;
	for(int i=0; i<v.size(); i++)
	{
		v[i].n[0] = ~v[i].n[0];
		v[i].n[0] &= (One<<len)-2;
	}
	neighbor_hit_detect_name=format_name(len,target);
	neighbor_hit_detect_name2=format_cover_name(len,target);
	vector<bigint> vexcl;
	readfile(neighbor_hit_detect_name2, vexcl, len);
	readfile(neighbor_hit_detect_name, vexcl, len);
	printf("%d candidates\n", v.size());
	size_t nEx = exclude(v, vexcl);
	printf("%d excluded\n", nEx);

	JobArray array;
	array.submit(&v[0], v.size(), len,target,target+1, true, 0, 0, neighbor_hit_detect);
	fclose(neighbor_hit_file);
	neighbor_hit_file=0;
	printf("%d -> %d\n", v.size(), array.sum_populations[target]);
	FILE* ff=fopen(PATH "neighbor-search-log.txt", "a");
	fprintf(ff, "%d\t%d\t%.3f\t(inv)\n", v.size(), array.sum_populations[target], array.realtime/g_fFreq);
	fclose(ff);
}


void do510_part2()
{
	g_kmax[0]=3;
	g_kmax[1]=8;
	g_max_ops=1e5;
	narrow(54,128,130);
	narrow(54,130,132);
	narrow(54,132,134);
	narrow(54,134,136);
	narrow(54,136,138);
	narrow(54,138,140);
	g_max_ops=0;
	review(54,3,8,140);
}

void do510_part1(int d1, int d2, int d3, int nRepeat1, int nRepeat2, int nRepeat3, int nMC)
{
	g_max_ops=1e4;
	g_ns_exclusions.clear();
	uint64_t t0, t1;
	uint64_t time_zones[7];
	int hit_counts[7];
	int k0=3, k1=8;
	g_kmax[0]=k0;
	g_kmax[1]=k1;

	FILE* f;

	_unlink(PATH "510-45-103.txt");
	_unlink(PATH "510-50-103.txt");
	_unlink(PATH "510-55-103.txt");
	_unlink(PATH "cover510-45-103.txt");
	_unlink(PATH "cover510-50-103.txt");
	_unlink(PATH "cover510-55-103.txt");

	_unlink(PATH "510-45-105.txt");
	_unlink(PATH "510-50-105.txt");
	_unlink(PATH "510-55-105.txt");
	_unlink(PATH "cover510-45-105.txt");
	_unlink(PATH "cover510-50-105.txt");
	_unlink(PATH "cover510-55-105.txt");

	_unlink(PATH "510-55-110.txt");
	_unlink(PATH "510-45-110.txt");
	_unlink(PATH "510-55-115.txt");
	_unlink(PATH "510-55-120.txt");
	_unlink(PATH "510-55-125.txt");
	_unlink(PATH "cover510-55-110.txt");
	_unlink(PATH "cover510-55-115.txt");
	_unlink(PATH "cover510-55-120.txt");
	_unlink(PATH "cover510-55-125.txt");
	_unlink(PATH "cover510-50-110.txt");
	_unlink(PATH "cover510-50-115.txt");
	_unlink(PATH "cover510-50-120.txt");
	_unlink(PATH "cover510-50-125.txt");
	f=fopen(PATH "cover510-55-128.txt", "rb");
	fseek(f, 0, SEEK_END);
	uint64_t pos = _ftelli64(f);
	fclose(f);
	if(pos>1e9)
		cover_rebuild(54,3,8,128);
	memset(hit_counts, 0, sizeof(hit_counts));
	t0=__rdtsc();
	hit_counts[0]=0;

	for(int pass=0; pass<5; pass++)
		hit_counts[0]+=mc(44, k0, k1, d1, nMC);

	string s = format_name(44, d1);
	lengthen(s.c_str(), k0, k1, 44, 54, d1);
	vector<bigint> v;
	readfile(s.c_str(), v, 44);
	rem_dup(v);
	vector<bigint> v2;

	uint64_t mul = 32;
	v2.resize(mul*v.size());
	for(uint64_t i=0; i<v.size(); i++)
		for(uint64_t j=0; j<mul; j++)
		{
			v2[i*mul+j]=v[i];
			v2[i*mul+j].set_len(49);
			for(int k=44; k<49; k++)
			if(j & (One<<(k-44)))
			{
				v2[i*mul+j].set(k);
			}
		}
	string s2=format_cover_name(49,d1);
	writefile(s2.c_str(), v2, 49);
	
	g_max_ops=1e4;
	t1=__rdtsc();
	time_zones[0]=t1-t0;
	g_ns_exclusions.push_back(make_pair(54,128));

	t0=__rdtsc();
	for(int pass=0; pass<nRepeat1; pass++)
		hit_counts[1]+=neighbor_depth_search(54,k0,k1,d1);
	hit_counts[2]=narrow(54, d1, d2);
	t1=__rdtsc();
	time_zones[1]=t1-t0;
	t0=__rdtsc();
	for(int pass=0; pass<nRepeat2; pass++)
		hit_counts[2]+=neighbor_depth_search(54,k0,k1,d2);
	hit_counts[3]=narrow(54, d2, 125);
	/*
	hit_counts[3]=narrow(54, d2, d3);
	for(int pass=0; pass<nRepeat3; pass++)
		neighbor_depth_search(54,k0,k1,d3);
	narrow(54, d3, 125);
	*/
	t1=__rdtsc();
	time_zones[2]=t1-t0;
	t0=__rdtsc();
	for(int pass=0; pass<16; pass++)
	{
		int n = neighbor_depth_search(54,k0,k1,125);
		hit_counts[3]+=n;
		if(n==0)
			break;
	}
	hit_counts[4]=narrow(54, 125, 128);
	t1=__rdtsc();
	time_zones[3]=t1-t0;
	t0=__rdtsc();
	g_ns_exclusions.clear();
	while(true)
	{
		int n = neighbor_depth_search(54,k0,k1,128);
		hit_counts[4]+=n;
		if(n==0)
			break;
	}
	t1=__rdtsc();
	time_zones[4]=t1-t0;
	f=fopen(PATH "log510.txt", "a");
	fprintf(f, "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%d\t%d\t%d\n",
		time_zones[0]/g_fFreq, time_zones[1]/g_fFreq, 
		time_zones[2]/g_fFreq, time_zones[3]/g_fFreq, time_zones[4]/g_fFreq, 
		hit_counts[0], hit_counts[1], hit_counts[2], hit_counts[3], hit_counts[4]);
	fclose(f);
}

int main(int argc, char* argv[])
{
	omp_set_num_threads(nthreads+1);
	omp_init_lock(&op_count_lock);
	//do77_part1();
	cover_rebuild(54, 3, 8, 128);
	do510_part1(110, 115, 120, 2, 3, 0, 1e7);
	return 0;
}
