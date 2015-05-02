#include "stdafx.h"
#include "graph.h"

#ifdef WIN32
#define align16 __declspec(align(16))
#else
#define align16 __attribute__(aligned(16))
#endif


void graph::merge()
{
	for(int color=0; color<8; color++)
	{
		for(int dim=0; dim<parent_cliques[color].size(); dim++)
		{
			size_t cur_size=cliques[color].size();
			size_t new_size=parent_cliques[color][dim]->size();
			cliques[color].resize(cur_size+new_size);
			memcpy(&cliques[color][cur_size], parent_cliques[color][dim]->ptr, new_size*sizeof(bigint));
		}
		parent_cliques[color].clear();
	}
}

void graph::check_first_unset_bit()
{
	
	first_unset_bit = 0;
	int i=0;
	while(i<MAX_LEN)
	{
		uint64_t t = ~(mask.n[i] | mask0.n[i]);
		if(t!=0)
		{
			unsigned long x;
			_BitScanForward64(&x, t);
			first_unset_bit=x+i*64+1;
			return;
		}
		i++;
	}
	first_unset_bit=MAX_LEN*64+1;
	/*
	while(first_unset_bit<MAX_LEN*64 && (mask.bit(first_unset_bit) || mask0.bit(first_unset_bit)))
		first_unset_bit++;
	first_unset_bit++;
	*/
}

void validate_graph_build(const graph& g);

void graph::update_shifted_masks(int color, int target, int n)
{
	if(search_range==0)
		return;
	/*
	bigint x;
	x.set_len(target);
	x.set(n-1);
	x.popcnt=1;
	cliques[1-color].push_back(x);
	clique_counts[1-color][0]++;
	*/
	for(int dx=0; dx<=search_range; dx++)
	{
		combo_ptrs[1-color][dx].set(target-n+dx);
	}
	for(int dx=0; dx<search_range; dx++)
	{
		if(n-1+dx+2<target)
			combo_ptrs[1-color][-dx-2].set(target-1-(n-1+dx+2));
		if(n-1<=dx)
			combo_ptrs[1-color][-dx-2].set(target-1-(dx-(n-1)));
	}
	if(color==0)
	{
		mask.set(n-1);
		if(mask.len<n)
			mask.set_len(n);
	}
	else
	{
		mask0.set(n-1);
		if(mask0.len<n)
			mask0.set_len(n);
	}
	check_first_unset_bit();
}

void graph::construct_shifted_masks(int target)
{
//	printf("graph::construct_shifted_masks(%d)\n", target);
	search_range = 0;
	check_first_unset_bit();
//	if(first_unset_bit-1 >= target-20)
//		return;

	search_range = target+1 - first_unset_bit + 1;
	if(search_range <= 0)
	{
		search_range = 0;
		return;
	}
	if(search_range >= MAX_LEN*64)
		search_range = MAX_LEN*64-1;
	combo_ptrs[0]=&combo_mask[0][MAX_LEN*64];
	combo_ptrs[1]=&combo_mask[1][MAX_LEN*64];

	invert(combo_ptrs[0][0], mask0);
	combo_ptrs[0][0].set_len(target);
	invert(combo_ptrs[1][0], mask);
	combo_ptrs[1][0].set_len(target);
	combo_ptrs[0][0] = combo_ptrs[0][0] << (target-mask0.len);
	combo_ptrs[1][0] = combo_ptrs[1][0] << (target-mask.len);
	int i;
#if MAX_LEN==2 && !defined(KNIGHTS_CORNER)
	__m128i a = _mm_loadu_si128((const __m128i*)&combo_ptrs[0][0].n[0]);
	__m128i b = _mm_loadu_si128((const __m128i*)&combo_ptrs[1][0].n[0]);
	for(i=1; i<=search_range; i++)
	{
		a = _mm_or_si128(_mm_slli_epi64(a, 1), _mm_slli_si128(_mm_srli_epi64(a, 63), 8));
		b = _mm_or_si128(_mm_slli_epi64(b, 1), _mm_slli_si128(_mm_srli_epi64(b, 63), 8));
		_mm_storeu_si128((__m128i*)(&combo_ptrs[0][i].n[0]), a);
		_mm_storeu_si128((__m128i*)(&combo_ptrs[1][i].n[0]), b);
	}
#else
	for(i=1; i<=search_range; i++)
	{
		combo_ptrs[0][i] = combo_ptrs[0][0] << i;
		combo_ptrs[1][i] = combo_ptrs[1][0] << i;
	}
#endif

	for(int t=0; t<search_range; t++)
	{
		bigint shift_mask = mask0;// << (target-mask0.len);
		shift_mask.set_len(target);
		combo_ptrs[0][-t-2].clear();
		combo_ptrs[1][-t-2].clear();
		invert(combo_ptrs[0][-t-2], shift_mask, target-(t+2));
		combo_ptrs[0][-t-2].unset(target);
		shift_mask = mask;
		shift_mask.set_len(target);
		invert(combo_ptrs[1][-t-2], shift_mask, target-(t+2));
	}
}


void graph::operator=(graph& g)
{
	n=g.n;
	mask0=g.mask0;
	mask=g.mask;
	first_unset_bit = g.first_unset_bit;
	search_range = g.search_range;
	g.merge();

	memcpy(clique_counts, g.clique_counts, sizeof(clique_counts));
	//inv_mask=g.inv_mask;
	for(int i=0; i<8; i++)
	{
		cliques[i].resize(g.cliques[i].size());
		if(g.cliques[i].size()>0)
			memcpy(&cliques[i][0], &g.cliques[i][0], g.cliques[i].size()*sizeof(bigint));
		parent_cliques[i].clear();
	}
}

void graph::print_short() const
{
	printf("dim %d\n", n);
	printf("mask:     ");
	mask0.print();
	printf("mask:     ");
	mask.print();
	int k;
	for(k=0; k<MAX_DEPTH; k++)
	{
		if(clique_counts[0][k]==0 && clique_counts[1][k]==0)
			continue;;
		printf("%d\t%d\t%d\n", k, clique_counts[0][k], clique_counts[1][k]);
	}
}

void graph::print() const
{
	printf("dim %d\n", n);
	printf("mask:     ");
	mask0.print();
	printf("mask:     ");
	mask.print();
	int k;
	for(k=0; k<MAX_DEPTH; k++)
	{
		if(clique_counts[0][k]==0 && clique_counts[1][k]==0)
			break;
		printf("%d\t%d\t%d\n", k, clique_counts[0][k], clique_counts[1][k]);
	}
	
	for(k=0; k<cliques[0].size(); k++)
	{
		printf("R ");
		cliques[0][k].print();
	}
	for(k=0; k<cliques[1].size(); k++)
	{
		printf("G ");
		cliques[1][k].print();
	}
	printf("Parent clique lists:\n");
	for(k=0; k<parent_cliques[0].size(); k++)
		printf("R %p\n", parent_cliques[0][k]);
	for(k=0; k<parent_cliques[1].size(); k++)
		printf("G %p\n", parent_cliques[1][k]);
}

void graph::reset()
{
	for(int i=0; i<8; i++)
	{
		cliques[i].clear();
		parent_cliques[i].clear();
	}
	memset(&clique_counts, 0, sizeof(clique_counts));
	n=0;
	best_next_bit=-1;
	second_best_next_bit=-1;
	save_base=-1;
}

//bool debug=false;
//bool verify_clique(const bigint& x, int target, const bigint& mask);

bool build_graph2(graph& gOut, const bigint2& mask, int kmax[2])
{
	gOut.mask0.clear();
	gOut.mask.clear();
	//gOut.inv_mask.clear();
	gOut.mask0.set_len(1);
	gOut.mask.set_len(1);
	//gOut.inv_mask.set_len(1);
	gOut.mask0.set(0);
	gOut.n=2;
	int i, j;
	for(i=0; i<2; i++)
		for(j=0; j<MAX_DEPTH; j++)
			gOut.clique_counts[i][j]=0;

	gOut.cliques[0].resize(0);
	gOut.cliques[1].resize(0);
	/*
	bigint& x = gOut.cliques[0][0];
	x.clear();
	x.set_len(mask.first.len-1);
	x.set(0);
	x.popcnt=1;
	gOut.clique_counts[0][0]=1;
	*/

	buffer<int> vLink[2];
	vLink[0].push_back(1);

	int nSetBits=0;
	if(popcnt(mask.first)!=0)
		nSetBits=mask.first.trailing_bit();
	if(popcnt(mask.second)!=0)
		nSetBits=max(nSetBits, mask.second.trailing_bit());

	int true_len = 0;//nSetBits+1;
	//mask.first.print();
	//mask.second.print();
	//while(true)
	for(int pos=1; pos<=nSetBits; pos++)
	{
		int color;
		if(mask.first.bit(gOut.n-1))
			color=0;
		else if(mask.second.bit(gOut.n-1))
			color=1;
		else
		{
			if(true_len==0)
				true_len = gOut.n;
			gOut.mask0.set_len(gOut.mask.len+1);
			gOut.mask.set_len(gOut.mask.len+1);
			gOut.n++;
			continue;
		}
		//printf("pos %d, gOut.n %d, color %d\n", pos, gOut.n, color);
		bigint target_mask;
		if(color)
			invert(target_mask, gOut.mask);
		else
			invert(target_mask, gOut.mask0);

		lite_vector& v = gOut.cliques[color];
		lite_vector& vFull = gOut.cliques[color+2];
		gOut.mask0.set_len(gOut.mask.len+1);
		gOut.mask.set_len(gOut.mask.len+1);
		if(color)
			gOut.mask.set(gOut.n-1);
		else
			gOut.mask0.set(gOut.n-1);
		gOut.n++;
		int i;

		const bigint& in_mask = color ? mask.second : mask.first;

		for(i=0; i<vLink[color].size(); i++)
		{
			int a = vLink[color][i];
			int b = gOut.n-1;
			int c = b-a;
			if(c<=0)
				printf("&");
			if(in_mask.bit(c-1))
			{
				v.resize(v.size()+1);
				bigint& x = v[v.size()-1];
				x.clear();
				x.set_len(gOut.n);
				x.set(a-1);
				x.set(b-1);
				x.popcnt=2;
				gOut.clique_counts[color][1]++;
			}
		}

		vLink[color].push_back(gOut.n-1);

		i=0;
		{

			lite_vector& vSrc = gOut.cliques[color];
			size_t nC = vSrc.size();
			for(i=0; i+16<=nC; i+=16)
			{
				int x = inv_batch_pattern_match(target_mask, &vSrc[i]);
				if(x==0)
					continue;
				for(int di=0; di<16; di++)
				//if(pattern_match(v[i], target_mask))
				if((x>>di)&1)
				{
					int pos = v.size();
					v.resize(pos+1);
					bigint& cl = v[pos];
					cl = vSrc[i+di];
					cl.set_len(gOut.n-1);
					cl.set(gOut.n-2);
					cl.popcnt++;
					gOut.clique_counts[color][cl.popcnt-1]++;
				}
			}

			for(; i<nC; i++)
			{
	//			uint64_t x = popcnt(v[i], target_mask);
	//			if(x==0)
				if(pattern_match(vSrc[i], target_mask))
				{
					bigint cl = vSrc[i];
					cl.set_len(gOut.n-1);
					cl.set(gOut.n-2);
					cl.popcnt++;
					v.push_back(cl);
					gOut.clique_counts[color][cl.popcnt-1]++;
				}
			}
		}
		if(gOut.clique_counts[color][kmax[color]]!=0)
			return false;
	}
	
	for(j=0; j<2; j++)
	{
		int i;
		lite_vector& v = gOut.cliques[j];
		size_t base = v.size();
		v.resize(base+vLink[j].size());
		for(i=0; i<vLink[j].size(); i++)
		{
			bigint& x = v[base+i];
			x.clear();
			x.set_len(gOut.n);
			x.set(vLink[j][i]-1);
			x.popcnt=1;
		}
		gOut.clique_counts[j][0]+=vLink[j].size();
	}

	for(int color=0; color<2; color++)
	{
		for(int n=0; n<gOut.cliques[color].size(); n++)
		{
			bigint& x = gOut.cliques[color][n];
			if(x.popcnt != popcnt(x))
				printf("P");
		}
	}

	if(true_len == 0)
		true_len = nSetBits+1;
	gOut.n = true_len+1;
	gOut.check_first_unset_bit();
	/*
	if(gOut.n != gOut.first_unset_bit)
	{
		printf("? %d %d\n", gOut.n, gOut.first_unset_bit);
	}
	*/
	return true;
}


bool count_c4s_on_node(const bigint& mask, const lite_vector& v, int i)
{
	int count=0;
	int j, k;
	vector<bool> observed;
//	observed.resize(v.size());
//	for(i=0; i<v.size(); i++)
//		observed[i]=false;
	{
		for(j=0; j<mask.len; j++)
		{
			if(j==i)
				continue;
			for(k=j+1; k<mask.len; k++)
			{ 
				if(k==i)
					continue;
				int c = mask.bit(i) + mask.bit(j) + mask.bit(k)
					+ mask.bit(abs(j-i)-1) + mask.bit(abs(k-j)-1) + mask.bit(abs(k-i)-1);
				if(c==6)
				{
					int hit=0;
					for(int m=0; m<v.size(); m++)
					{
						if(v[m].popcnt != 3)
							continue;
						if(v[m].bit(i) && v[m].bit(j) && v[m].bit(k))
						{
							hit++;
//								observed[m]=true;
						}
					}
				
					if(hit>1)
					{
//						printf("ERROR: duplicate clique %d, %d, %d\n", i, j, k);
					}
					if(hit==0)
					{
						printf("ERROR: missed clique %d, %d, %d\n", i, j, k);
						return false;
					}
				}
			}
		}
	}
	return true;
}

uint64_t g_set_bit_count = 0;

int g_min_depth;



struct recurse_object
{
	bare_bigint stack[MAX_DEPTH];
	bare_bigint stack2[MAX_DEPTH];
	bigint unsets[MAX_DEPTH];
	int bits[MAX_DEPTH];
	bare_bigint inverted_masks[2][MAX_LEN*64+1];
	int nMaxUnset;
	int nMinCliqueSize;
	int kmax;
	int color;
	int new_bit;
	int lastbit;

	bigint _m0, _m1;
	lite_vector* cliques;
	int* clique_counts[2];

	void prepare(int* kmax, int color, bool fastmode);
	bool unrolled_flat_search();
};

void recurse_object::prepare(int* km, int c, bool fastmode)
{
	bigint gap_mask;
	const bigint& in_mask = c ? _m1 : _m0;
	const bigint& dual_mask = c ? _m0 : _m1;

	if(fastmode)
		gap_mask = _m0 | _m1;
	else
		gap_mask = dual_mask;
	gap_mask.neg();
	gap_mask.set_len(lastbit+1);
	gap_mask.unset_from(lastbit);
	gap_mask = andnot(gap_mask, inverted_masks[c][new_bit]);
		
	bigint x = in_mask;
	x.set_len(lastbit+1);
	x.unset_from(lastbit+1);
	x.unset(new_bit);
	x &= inverted_masks[1-c][new_bit];
	stack[0] = x;
//	bits[0]=lastbit;
	unsets[0].clear();
	stack2[0]=gap_mask;
	this->kmax = km[c];
	this->color = c;
}

bool recurse_object::unrolled_flat_search()
{
	lite_vector& v = cliques[color];
	bare_bigint* inv_masks = inverted_masks[1-color];
	int depth=0;

	const bigint& in_mask = color ? _m1 : _m0;

	bigint v_curr_cl[MAX_DEPTH];
	v_curr_cl[0].clear();
	v_curr_cl[0].set_len(lastbit+1);
	v_curr_cl[0].set(new_bit);
	while(!stack[0].zero())
	{
		int b0 = stack[0].trailing_bit();
		bits[0]=b0;
		//printf("b0 %d, stack[0] %llx %llx\n", bits[0], stack[0], stack[1]);
		depth=1;
		store_and(stack[1], stack[0], inv_masks[b0]);
		//stack[1]=stack[0] & inv_masks[b0];
		bits[1]=b0-1;

		v_curr_cl[1]=v_curr_cl[0];
		v_curr_cl[1].set(b0);

		if(2>=nMinCliqueSize)
		{
			if(2>=kmax)
				return false;
			uint64_t cls = popcnt(stack[1]);
			uint64_t pos = v.size();
			v.resize(pos+cls);
			bigint y = stack[1];
			bigint current_clique;
			current_clique.clear();
			current_clique.set_len(lastbit+1);
			current_clique.set(new_bit);
			current_clique.set(b0);
			current_clique.popcnt=2;

			for(uint64_t i=0; i<cls; i++)
			{
				bigint& t = v[i+pos];
				t = current_clique;
				t.popcnt+=1;
				int b = y.trailing_bit();
				t.set(b);
				y.unset(b);
			//	printf("   ");t.printbits();
			}
			clique_counts[color][depth+1]+=cls;
				
			if(new_bit < bits[0])
			{
				v.resize(pos+cls*2);
				bigint y = stack[1];
				for(uint64_t i=0; i<cls; i++)
				{
					bigint& t = v[i+pos+cls];
					t.clear();
					bigint y = v[i+pos];
					y.set_len(bits[0]);
					invert(t, y);
					t.set(bits[0]);
					t.popcnt=depth+2;
					t.set_len(bits[0]+1);
				}
				clique_counts[color][depth+1]+=cls;
			}
		}

		int pc = popcnt(stack[1]);
		int x = 1+pc;
		if(x<nMinCliqueSize)
			goto done2;
		if(x==nMinCliqueSize)
		{
			if(pc>=2)
			{
				bigint save = stack[depth];
				bool fail=false;
				for(int d=0; d+1<nMinCliqueSize-depth; d++)
				{
					int pos = stack[depth].trailing_bit();
					stack[depth] &= inv_masks[pos];
					stack[depth].unset(pos);
					if(stack[depth].zero())
					{
						fail=true;
						break;
					}
				}
				if(!fail)
				{
					bigint current_clique=v_curr_cl[depth-1];
					current_clique.set(bits[depth-1]);
					current_clique|=save;
					current_clique.popcnt=nMinCliqueSize+1;
					clique_counts[color][nMinCliqueSize]++;
					if(nMinCliqueSize==kmax)
						return false;
					v.push_back(current_clique);
				}
			}
			goto done2;
		}

		while(true)
		{
			bits[depth]=stack[depth].trailing_bit();
			depth++;
			//stack[depth]=stack[depth-1] & inv_masks[bits[depth-1]];
			store_and(stack[depth], stack[depth-1], inv_masks[bits[depth-1]]);
			int x = depth+popcnt(stack[depth]);
			if(x<=nMinCliqueSize)
			{
				/*
				if(depth==1)
					goto done2;
				depth--;
				stack[depth].unset(bits[depth]);
				goto next_loop;
			}
			else if(x==nMinCliqueSize)
			{
			*/
				if(x==nMinCliqueSize && !stack[depth].zero())
				{
					bigint save = stack[depth];
					//bigint save = temp;
					//int rem_bits[MAX_DEPTH];
					bool fail=false;
					for(int d=0; d+1<nMinCliqueSize-depth; d++)
					{
						int pos = stack[depth].trailing_bit();
//						rem_bits[d]=pos;
						stack[depth] &= inv_masks[pos];
						stack[depth].unset(pos);
//						temp.unset(pos);
						if(stack[depth].zero())
						{
							fail=true;
							break;
						}
					}
//					if(stack[depth].zero())
//						fail=true;
					/*
					int pos = temp.trailing_bit();
					rem_bits[nMinCliqueSize-depth-1]=pos;
					
					int fail=0;
					for(int i=0; i<nMinCliqueSize-depth-1; i++)
						for(int j=i+1; j<nMinCliqueSize-depth; j++)
						{
							if(!in_mask.bit(abs(rem_bits[i]-rem_bits[j])-1))
								fail++;
						}
					
					if((fail==0) != (!stack[depth].zero()))
					{
						printf("depth=%d nMinCliqueSize=%d\n", depth, nMinCliqueSize);
						save.set_len(MAX_LEN*64);
						save.printbits();
						stack[depth].printbits();
						in_mask.printbits();
						printf("fail %d, zero() %d\n", fail);
						exit(-1);
					}
					
					//if(!stack[depth].zero())
					if(fail==0)
					*/
					if(!fail)
					{
						bigint current_clique=v_curr_cl[depth-1];
						current_clique.set(bits[depth-1]);
						current_clique|=save;
						current_clique.popcnt=nMinCliqueSize+1;
						clique_counts[color][nMinCliqueSize]++;
						if(nMinCliqueSize==kmax)
							return false;
						v.push_back(current_clique);
					}
				}
				if(depth==1)
					goto done2;
				depth--;
				stack[depth].unset(bits[depth]);
				goto next_loop;
			}
			
			bits[depth]=bits[depth-1]-1;
			v_curr_cl[depth]=v_curr_cl[depth-1];
			v_curr_cl[depth].set(bits[depth-1]);
				//stack[depth].clear();
			if(depth+1>=nMinCliqueSize)
			{
				//printf("x %llx %llx, cls %d\n", x.n[0], x.n[1], cls);
				if(!stack[depth].zero())
				{
					if(depth+1>=kmax)
						return false;
					bigint x=stack[depth];
					uint64_t cls = popcnt(x);
					uint64_t pos = v.size();
					v.resize(pos+cls);
					bigint y = x;
					bigint current_clique=v_curr_cl[depth];
					/*
					current_clique.clear();
					current_clique.set_len(lastbit+1);
					current_clique.set(new_bit);
					current_clique.set(b0);
					for(uint64_t i=1; i<depth; i++)
						current_clique.set(bits[i]);
					*/
					current_clique.popcnt=depth+1;

					for(uint64_t i=0; i<cls; i++)
					{
						bigint& t = v[i+pos];
						t = current_clique;
						t.popcnt+=1;
						int b = y.trailing_bit();
						t.set(b);
						y.unset(b);
					//	printf("   ");t.printbits();
					}
					clique_counts[color][depth+1]+=cls;
				
					if(new_bit < bits[0])
					{
						v.resize(pos+cls*2);
						bigint y = x;
						for(uint64_t i=0; i<cls; i++)
						{
							bigint& t = v[i+pos+cls];
							t.clear();
							bigint y = v[i+pos];
							y.set_len(bits[0]);
							invert(t, y);
							t.set(bits[0]);
							t.popcnt=depth+2;
							t.set_len(bits[0]+1);
						}
						clique_counts[color][depth+1]+=cls;
					}
				}
			}
		next_loop:
			while(stack[depth].zero())
			{
				if(depth==1)
					goto done2;
				depth--;
				stack[depth].unset(bits[depth]);
			}
		}
done2:
		stack[0].unset(b0);
	}
	return true;
}

bool graph::set_bit(int color, int new_bit, int kmax[2], bool exact, bigint* badmask)
{
	bigint& _m0 = mask0;
	bigint& _m1 = mask;
	if((color ? _m1 : _m0).bit(new_bit))
		return true;
	g_set_bit_count++;
	const bigint& in_mask = color ? _m1 : _m0;

//	bigint inverted_masks[MAX_LEN*64+1];
	int max_bit = (color ? _m1 : _m0).trailing_bit();
	if(color)
		_m1.set(new_bit);
	else
		_m0.set(new_bit);
	int i;
	if(new_bit > max_bit)
		max_bit = new_bit;

	recurse_object obj;
	obj.clique_counts[0] = & clique_counts[0][0];
	obj.clique_counts[1] = & clique_counts[1][0];
	obj.cliques = &cliques[0];

	for(i=0; i<=max_bit; i++)
	{
		if(!new_bit && !in_mask.bit(i))
			continue;
		obj.inverted_masks[1-color][i].clear();
		invert(obj.inverted_masks[1-color][i], in_mask, i);
	}

	check_first_unset_bit();
	int nMinCliqueSize=kmax[color]-2;
	/*
	if(kmax[color]<6 && kmax[1-color]<6)
	{
		if(first_unset_bit >= 80)
			nMinCliqueSize=kmax[color]-1;
		if(first_unset_bit >= 90)
			nMinCliqueSize=kmax[color];
	}
	*/
	obj.nMaxUnset=0;
	obj._m0 = mask0;
	obj._m1 = mask;
	obj.nMinCliqueSize = nMinCliqueSize;
	obj.new_bit = new_bit;
	obj.lastbit = max_bit;

	/*
	bigint stack[MAX_DEPTH];
	stack[0].clear();
	//int lastbit[MAX_DEPTH];
	int depth=0;
	lite_vector& v = cliques[color];
	bigint x = in_mask;
	x.set_len(max(max_bit, new_bit+1));
	x.unset(new_bit);
	lastbit[0]=max_bit-1;
	x.unset_from(max(new_bit, max_bit));
	x &= inverted_masks[new_bit];
	stack[0] = x;
	*/
	obj.prepare(kmax, color, true);
	if(!obj.unrolled_flat_search())
		return false;
		
	if(clique_counts[color][kmax[color]]!=0)
		return false;

	return true;
}

uint64_t g_build_graph_calls=0;

bool build_graph3(graph& gOut, const bigint2& mask, int kmax[2])
{
	gOut.mask0.clear();
	gOut.mask.clear();
	gOut.mask0.set_len(1);
	gOut.mask.set_len(1);
	gOut.mask0.set(0);
	gOut.n=2;
	int i, j;
	for(i=0; i<2; i++)
		for(j=0; j<MAX_DEPTH; j++)
			gOut.clique_counts[i][j]=0;

	for(i=0; i<8; i++)
	{
		gOut.cliques[i].resize(0);
		gOut.parent_cliques[i].resize(0);
	}
/*
	buffer<int> vLink[2];
	vLink[0].push_back(1);
*/
	int nSetBits=0;
	if(popcnt(mask.first)!=0)
		nSetBits=mask.first.trailing_bit();
	if(popcnt(mask.second)!=0)
		nSetBits=max(nSetBits, mask.second.trailing_bit());
	g_build_graph_calls+=popcnt(mask.first)+popcnt(mask.second);

	recurse_object obj;
	obj.clique_counts[0] = & gOut.clique_counts[0][0];
	obj.clique_counts[1] = & gOut.clique_counts[1][0];
	obj.cliques = &gOut.cliques[0];

	for(i=0; i<=nSetBits; i++)
	{
		if(mask.first.bit(i))
		{
			obj.inverted_masks[1][i].clear();
			invert(obj.inverted_masks[1][i], mask.first, i);
		}
		if(mask.second.bit(i))
		{
			obj.inverted_masks[0][i].clear();
			invert(obj.inverted_masks[0][i], mask.second, i);
		}
	}


	obj.nMaxUnset=0;

	//bigint inverted_masks[2][MAX_LEN*64];
	align16 uint64_t inverted_masks[2][MAX_LEN*64*2];
	int max_bit = max(mask.first.trailing_bit(), mask.second.trailing_bit())+1;
	if(max_bit>MAX_LEN*64)
		max_bit=MAX_LEN*64;
	for(i=0; i<max_bit; i++)
	{
		bigint m;
		invert(m, mask.first, i);
		inverted_masks[0][i*2+0]=m.n[0];
		inverted_masks[0][i*2+1]=m.n[1];
		invert(m, mask.second, i);
		inverted_masks[1][i*2+0]=m.n[0];
		inverted_masks[1][i*2+1]=m.n[1];
	}

	int true_len = 0;

	for(int pos=1; pos<=nSetBits; pos++)
	{
		int color;
		if(mask.first.bit(gOut.n-1))
			color=0;
		else if(mask.second.bit(gOut.n-1))
			color=1;
		else
		{
			if(true_len==0)
				true_len = gOut.n;
			gOut.mask0.set_len(gOut.mask.len+1);
			gOut.mask.set_len(gOut.mask.len+1);
			gOut.n++;
			continue;
		}
		gOut.mask0.set_len(gOut.mask.len+1);
		gOut.mask.set_len(gOut.mask.len+1);
		if(color)
			gOut.mask.set(gOut.n-1);
		else
			gOut.mask0.set(gOut.n-1);
		obj._m0 = gOut.mask0;
		obj._m1 = gOut.mask;

		obj.new_bit = gOut.n-1;
		obj.lastbit = gOut.n-1;
		obj.nMinCliqueSize=kmax[color]-2;
		obj.prepare(kmax, color, true);
		if(!obj.unrolled_flat_search())
			return false;
		gOut.n++;
		if(gOut.clique_counts[color][kmax[color]]!=0)
			return false;
	}

	if(true_len == 0)
		true_len = nSetBits+1;
	gOut.n = true_len+1;
	gOut.check_first_unset_bit();
	return true;
}

bool extend(graph& gOut, const graph& g, int color, int kmax[2], int target, int new_bit, int search_mode, bool debug, int force_range)
{
	if(force_range == 0)
		force_range = target;
	validate_graph_build(g);
	if(color && g.mask0.bit(new_bit))
	{
		if(debug)
			printf("Can't set %d->%d (bit set to 0)\n", new_bit, color);
		return false;
	}
	if((!color) && g.mask.bit(new_bit))
	{
		if(debug)
			printf("Can't set %d->%d (bit set to 1)\n", new_bit, color);
		return false;
	}

	//printf("extend %p -> %p\n", &g, &gOut);
	gOut.n = g.n+1;
	gOut.mask0 = g.mask0;
	gOut.mask = g.mask;
	if(gOut.mask0.len < new_bit+1)
		gOut.mask0.set_len(new_bit+1);
	if(gOut.mask.len < new_bit+1)
		gOut.mask.set_len(new_bit+1);
	for(int i=0; i<8; i++)
	{
		gOut.cliques[i].clear(); 
		gOut.parent_cliques[i] = g.parent_cliques[i];
		gOut.parent_cliques[i].push_back(&g.cliques[i]);
	}

	if(g.set_extensions)
		gOut.parent_ext = &g;
	else if(g.parent_ext != 0)
		gOut.parent_ext = g.parent_ext;
	else
		gOut.parent_ext = 0;

	memcpy(gOut.clique_counts, g.clique_counts, sizeof(g.clique_counts));
	//int new_bit = g.n-1;
	bool updated=false;
	int sizePrev = gOut.cliques[color].size();
	if(!(color ? gOut.mask.bit(new_bit) : gOut.mask0.bit(new_bit)))
	{
//		if(g.n>=52)
//		printf("extend: setting bit %d to %d\n", new_bit, color);
		if(!gOut.set_bit(color, new_bit, kmax))
		{
//			validate_graph_build(gOut);
			if(debug)
				printf("set_bit() returns error\n");

			return false;
		}
		updated=true;
	}
	if(!updated)
		return true;

	if(gOut.clique_counts[0][kmax[0]]!=0 || gOut.clique_counts[1][kmax[1]]!=0)
	{
//		validate_graph_build(gOut);
			if(debug)
				printf("kmax check\n");
		return false;
	}

	gOut.check_first_unset_bit();
	if(gOut.first_unset_bit-1 < target)
	{
		if(search_mode==1)
		{
			printf("ERROR: search_mode=1?\n");
			exit(-1);
		}
		if(search_mode==2)
		{			
			gOut.construct_shifted_masks(target);
			// really a back and forth: color 0 -> 1 -> 0 ..., could result in much cleaner/compact code
			if(gOut.search_range > 0)
			{
				buffer<int>* new_marks = gOut.new_marks;
				new_marks[0].resize(0);
				new_marks[1].resize(0);
				gOut.forced_bit_search2(kmax, target, color, new_bit, new_marks, 0, debug);//sizePrev);
				if(debug)
				{
					for(int i=0; i<2; i++)
						for(int j=0; j<new_marks[i].size(); j++)
						{
							if(debug)
								printf(" forced_bit_search2: setting %d -> %d\n", i, new_marks[i][j]-1);
						}
				}
				int changes=new_marks[0].size()+new_marks[1].size();
				if(changes!=0)
				{
					int prev_pos0=gOut.cliques[0].size();
					int prev_pos1=gOut.cliques[1].size();
					bigint badmask;
					badmask.clear();
					badmask.set_len(127);
					int x = gOut.new_clique_search(kmax, target, new_marks, debug ? &badmask : 0);
					if(x>0)
					{
						if(debug)
						{
							printf("new_clique_search return error (1)\n");
							badmask.printbits();
						}
						return false;
					}
					while(true)
					{
						if(prev_pos0==gOut.cliques[0].size()
							&& prev_pos1==gOut.cliques[1].size())
							break;
						/*
						gOut.construct_shifted_masks(target);
						if(gOut.search_range <= 0)
							break;
						*/
						for(int i=0; i<2; i++)
							for(int j=0; j<new_marks[i].size(); j++)
							{
								gOut.update_shifted_masks(1-i, target, new_marks[i][j]);
							}

						int color_mask=3;
						if(new_marks[0].size()>0)
							color_mask|=1;
						if(new_marks[1].size()>0)
							color_mask|=2;
						new_marks[0].resize(0);
						new_marks[1].resize(0);
						gOut.forced_bit_search(kmax, target, new_marks, false, color_mask, prev_pos0, prev_pos1, debug);
						int changes=new_marks[0].size()+new_marks[1].size();
						if(changes==0)
							break;
						if(debug)
						{
							for(int i=0; i<2; i++)
								for(int j=0; j<new_marks[i].size(); j++)
								{
										printf(" forced_bit_search: setting %d -> %d\n", i, new_marks[i][j]-1);
								}
						}						
						prev_pos0=gOut.cliques[0].size();
						prev_pos1=gOut.cliques[1].size();
						bigint badmask;
						badmask.clear();
						badmask.set_len(127);

						int x = gOut.new_clique_search(kmax, target, new_marks, debug ? &badmask : 0);
						if(x>0)
						{
		//					validate_graph_build(gOut);
							if(debug)
							{
								printf("new_clique_search() returns error (2)\n");
								badmask.printbits();
							}
							return false;
						}
						//break;
					}
				}
			}
		}
	}
	if(gOut.clique_counts[0][kmax[0]]!=0 || gOut.clique_counts[1][kmax[1]]!=0)
	{
		if(debug)
			printf("kmax check 2\n");

		return false;
	}

	bigint error_mask = gOut.mask & gOut.mask0;
	if(popcnt(error_mask)!=0)
	{
		if(debug)
		{
			printf("popcnt(error_mask)!=0\n");
			gOut.mask.printbits();
			gOut.mask0.printbits();
		}
		return false;
	}

	validate_graph_build(gOut);
	gOut.check_first_unset_bit();
	return true;
}

bool extend(graph& g, int color, int kmax[2], int target, int new_bit)
{
	if(color && g.mask0.bit(new_bit))
		return false;
	if((!color) && g.mask.bit(new_bit))
		return false;

	g.n++;
	if(g.mask0.len < new_bit+1)
		g.mask0.set_len(new_bit+1);
	if(g.mask.len < new_bit+1)
		g.mask.set_len(new_bit+1);
	bool updated=false;
	if(!(color ? g.mask.bit(new_bit) : g.mask0.bit(new_bit)))
	{
		if(!g.set_bit(color, new_bit, kmax))
		{
//			validate_graph_build(gOut);
			return false;
		}
		updated=true;
	}

	if(g.clique_counts[0][kmax[0]]!=0 || g.clique_counts[1][kmax[1]]!=0)
	{
//		validate_graph_build(gOut);
		return false;
	}


	bigint error_mask = g.mask & g.mask0;
	if(popcnt(error_mask)!=0)
		return false;

	g.check_first_unset_bit();
	return true;
}

int g_sym_size=0;

int g_force_sym=0;

bool symmetry_check(graph& g, int kmax[2], int target, bool debug)
{
	buffer<int>* new_marks = g.new_marks;

	if(g_sym_size==0 && g_force_sym==0)
		return true;
	
	int i;
		
	int last_set_bit=0;
	if(!g.mask0.zero())
		last_set_bit=g.mask0.trailing_bit();
	if(!g.mask.zero())
		last_set_bit=max(last_set_bit, g.mask.trailing_bit());

	if(g_force_sym!=0)
	{
		if(last_set_bit*3<=g_force_sym)
			return true;
		new_marks[0].clear();
		new_marks[1].clear();
		int sym_size = g_force_sym;
		for(i=0; i<target; i++)
		{
			if(sym_size-1-i>=target-1)
				continue;
			if(g.mask0.bit(i))
			{
				if(g.mask.bit(sym_size-1-i))
					return false;
				if(!g.mask0.bit(sym_size-1-i))
					new_marks[0].push_back(sym_size-i);
			}
			if(g.mask.bit(i))
			{
				if(g.mask0.bit(sym_size-1-i))
					return false;
				if(!g.mask.bit(sym_size-1-i))
					new_marks[1].push_back(sym_size-i);
			}
		}
		if(new_marks[0].size()>0 || new_marks[1].size()>0)
		{
			int x = g.new_clique_search(kmax, target, new_marks);
			if(x>0)
			{
				if(debug)
					printf("g.new_clique_search returns %d\n", x);
				return false; 
			}
		}
		int max_set_bit = 0;
		if(popcnt(g.mask)>0)
			max_set_bit = g.mask.trailing_bit()+1;
		if(popcnt(g.mask0)>0)
			max_set_bit = max(max_set_bit, g.mask0.trailing_bit()+1);
		g.mask.set_len(max_set_bit+1);
		g.mask0.set_len(max_set_bit+1);

		return true;
	}

	if(last_set_bit*2>g_sym_size)
	{
		int nchoices=0;
		int sym_size=0;
		for(int x=target; x<=g_sym_size; x++)
		{
			bool ok = true;
			if(g_sym_size>MAX_LEN*64)
			{
				for(i=0; i<target; i++)
				{
					if(x-1-i>=MAX_LEN*64)
						continue;
					if(g.mask0.bit(i) && g.mask.bit(x-1-i))
					{
						ok=false;
						break;
					}
					if(g.mask.bit(i) && g.mask0.bit(x-1-i))
					{
						ok=false;
						break;
					}
				}
			}
			else
			{
				bigint m = g.mask, m0=g.mask0;
				m.set_len(x);
				m0.set_len(x);
				bigint minv, m0inv;
				invert(minv, m);
				invert(m0inv, m0);
				ok = (popcnt(minv, m0)==0 && popcnt(m0inv, m)==0);
				/*
				bool ok2 = true;
				for(i=0; i<target; i++)
				{
					if(x-1-i>=MAX_LEN*64)
						continue;
					if(g.mask0.bit(i) && g.mask.bit(x-1-i))
					{
						ok2=false;
						break;
					}
					if(g.mask.bit(i) && g.mask0.bit(x-1-i))
					{
						ok2=false;
						break;
					}
				}
				if(ok != ok2)
				{
					printf("Mismatch: %d %d %d %d\n", ok, ok2, target, x);
					g.mask.print();
					g.mask0.print();
					g.mask.printbits();
					g.mask0.printbits();
					exit(-1);
				}
				*/
			}

			if(ok)
			{
				nchoices++;
				sym_size=x;
				if(nchoices>1)
					break;
			}
		}
		if(nchoices==0)
			return false;
		//int sym_size=...;
		if(nchoices==1)
		{
			new_marks[0].clear();
			new_marks[1].clear();
			for(i=0; i<target; i++)
			{
				if(sym_size-1-i>=target-1)
					continue;
				if(g.mask0.bit(i))
				{
					if(g.mask.bit(sym_size-1-i))
						return false;
					if(!g.mask0.bit(sym_size-1-i))
						new_marks[0].push_back(sym_size-i);
				}
				if(g.mask.bit(i))
				{
					if(g.mask0.bit(sym_size-1-i))
						return false;
					if(!g.mask.bit(sym_size-1-i))
						new_marks[1].push_back(sym_size-i);
				}
			}
			if(new_marks[0].size()>0 || new_marks[1].size()>0)
			{
				int x = g.new_clique_search(kmax, target, new_marks);
				if(x>0)
				{
					if(debug)
						printf("g.new_clique_search returns %d\n", x);
					return false; 
				}
			}
			int max_set_bit = 0;
			if(popcnt(g.mask)>0)
				max_set_bit = g.mask.trailing_bit()+1;
			if(popcnt(g.mask0)>0)
				max_set_bit = max(max_set_bit, g.mask0.trailing_bit()+1);
			g.mask.set_len(max_set_bit+1);
			g.mask0.set_len(max_set_bit+1);
		}
	}
	return true;
}

int max_extensibility(graph& g, int kmax[2], int target, bool debug)
{
	int i;

	g.construct_shifted_masks(target);
	int min_unset_bit = 35;
	if(kmax[0]<5 || kmax[1]<5)
		min_unset_bit=10;
	if(g.first_unset_bit-1 < min_unset_bit || g.first_unset_bit-1 >= target-20)
	{
		//printf("First unset %d, target %d\n", g.first_unset_bit, target);
		g.set_extensions = false;
		return INT_MAX;
	}

	buffer<int>* new_marks = g.new_marks;

	for(i=2; i<8; i++)
		g.cliques[i].clear();

	//int minpos0=g.cliques[0].size(), minpos1=g.cliques[1].size();
	int minpos0=0, minpos1=0;
	bool first_pass=true;
	bool first_pass_ls=true;

	buffer<pair<int,int> >& new_bits = g.new_bits;

repeat:
	while(true)
	{
		new_marks[0].resize(0);
		new_marks[1].resize(0);
		g.forced_bit_search(kmax, target, new_marks, first_pass, 3, minpos0, minpos1);
		first_pass=false;
		if(debug)
		{
			for(int i=0; i<2; i++)
				for(int j=0; j<new_marks[i].size(); j++)
					printf("forced_bit_search: found %d %d\n", i, new_marks[i][j]);
		}
		minpos0=g.cliques[0].size();
		minpos1=g.cliques[1].size();
		int changes=new_marks[0].size()+new_marks[1].size();
		if(changes==0)
			break;
		int x = g.new_clique_search(kmax, target, new_marks);
		if(x>0)
		{
			if(debug)
				printf("(1) g.new_clique_search returns %d\n", x);
			return x;
		}
		g.construct_shifted_masks(target);
	}

	for(i=0; i<target; i++)
	{
		if(g.mask0.bit(i) && g.mask.bit(i))
		{
			//printf("Bits %d set(1)\n", i);
			return i;
		}
	}

	new_bits.clear();
	if(!link_pair_search(target, kmax, g.combo_ptrs, g, new_bits, debug, first_pass_ls))
	{
		if(debug)
			printf("link_pair_search fail\n");
		return target-10;
	}
	first_pass_ls=false;

	if(new_bits.size()>0)
	{
		for(int i=0; i<new_bits.size(); i++)
		{
				if(debug)
				{
					printf("link_pair_search: found %d %d\n", 1-new_bits[i].first, new_bits[i].second);
				}

			//g.update_shifted_masks(new_bits[i].first, target, new_bits[i].second);
			new_marks[1-new_bits[i].first].push_back(new_bits[i].second);
		}

		for(int i=0; i<target; i++)
		{
			if(g.mask0.bit(i) && g.mask.bit(i))
			{
				if(debug)
					printf("Bits %d set\n", i);
				return i;
			}
		}

		int x = g.new_clique_search(kmax, target, new_marks);
		if(x>0)
		{
			if(debug)
				printf("g.new_clique_search returns %d\n", x);
			return x; 
		}
		g.construct_shifted_masks(target);
		goto repeat;
	}


	int max_set_bit = 0;
	if(popcnt(g.mask)>0)
		max_set_bit = g.mask.trailing_bit()+1;
	if(popcnt(g.mask0)>0)
		max_set_bit = max(max_set_bit, g.mask0.trailing_bit()+1);
	g.mask.set_len(max_set_bit+1);
	g.mask0.set_len(max_set_bit+1);
	return INT_MAX;
}


int bit(uint64_t a, int pos) { return (a>>pos)&1; }

int evolve(uint64_t& a, uint64_t& b, buffer<pair<uint64_t,uint64_t> > rules[2])
{
	int sum_changes=0;
	while(true)
	{
		int changes=0;

		uint64_t aprev=a, bprev=b;
		for(int j=0; j<rules[0].dim; j++)
		{
			if(a & rules[0][j].first)
				b|=rules[0][j].second;
			if(a & rules[0][j].second)
				b|=rules[0][j].first;
		}
		for(int j=0; j<rules[1].dim; j++)
		{
			if(b & rules[1][j].first)
				a|=rules[1][j].second;
			if(b & rules[1][j].second)
				a|=rules[1][j].first;
		}

		if(aprev==a && bprev==b)
			break;
		sum_changes++;
	}
	return sum_changes;
}

int evolve(uint64_t& a, uint64_t& b, buffer<pair<int,int> > rules[2])
{
	int sum_changes=0;
	while(true)
	{
		int changes=0;

		uint64_t aprev=a, bprev=b;
		
		for(int j=0; j<rules[0].dim; j++)
		{
			if(bit(a,rules[0][j].first))
				b|=One<<rules[0][j].second;
			if(bit(a,rules[0][j].second))
				b|=One<<rules[0][j].first;
		}
		for(int j=0; j<rules[1].dim; j++)
		{
			if(bit(b,rules[1][j].first))
				a|=One<<rules[1][j].second;
			if(bit(b,rules[1][j].second))
				a|=One<<rules[1][j].first;
		}
		
		if(aprev==a && bprev==b)
			break;
		sum_changes++;
	}
	return sum_changes;
}

int evolve(bigint& a, bigint& b, buffer<pair<int,int> > rules[2])
{
	int sum_changes=0;
	while(true)
	{
		int changes=0;

		bigint aprev=a, bprev=b;
		
		for(int j=0; j<rules[0].dim; j++)
		{
			if(a.bit(rules[0][j].first))
			{
//				printf("0 %d %d\n", rules[0][j].first, rules[0][j].second);
				b.set(rules[0][j].second);
			}
			if(a.bit(rules[0][j].second))
			{
				b.set(rules[0][j].first);
			}
		}
		for(int j=0; j<rules[1].dim; j++)
		{
			if(b.bit(rules[1][j].first))
				a.set(rules[1][j].second);
			if(b.bit(rules[1][j].second))
				a.set(rules[1][j].first);
		}
		

		if(aprev==a && bprev==b)
			break;
		sum_changes++;
	}
	return sum_changes;
}

int evolve_v2(bigint& a, bigint& b, buffer<int> rules_int[2], bigint rules[2][64], int base, int max_dim)
{
	int sum_changes=0;
	while(true)
	{
		int changes=0;

		bigint aprev=a, bprev=b;
		
		for(int j=0; j<rules_int[0].size(); j++)
		{
			int n = rules_int[0][j];
			if(a.bit(base+n))
				b |= rules[0][n];
		}
		for(int j=0; j<rules_int[1].size(); j++)
		{
			int n = rules_int[1][j];
			if(b.bit(base+n))
				a |= rules[1][n];
		}		

		if(aprev==a && bprev==b)
			break;
		sum_changes++;
	}
	return sum_changes;
}

void printmask(uint64_t a0, uint64_t b0, int n)
{
	for(int i=0; i<n; i++)
	{
		if(((a0>>i)&1) && ((b0>>i)&1))
			printf("#");
		else if((a0>>i)&1)
			printf("0");
		else if((b0>>i)&1)
			printf("1");
		else
			printf(".");
	}
	printf("\n");
}

void printmask(const bigint& a0, const bigint& b0, int n)
{
	for(int i=0; i<n; i++)
	{
		if(a0.bit(i) && b0.bit(i))
			printf("#");
		else if(a0.bit(i))
			printf("0");
		else if(b0.bit(i))
			printf("1");
		else
			printf(".");
	}
	printf("\n");
}

bool construct_incomplete_cliques_simple(lite_vector* part_cliques, graph& g, int target, int kmax[2], bool use_parents)
{
	int i, k;
	int base = g.first_unset_bit-1;
	int max_dim=target-base;
	g.construct_shifted_masks(target);
	bigint setmask = g.mask | g.mask0;
	
	char temp[MAX_LEN*64*MAX_LEN*64];
	char temp2[MAX_LEN*64];
	memset(temp, 0, target*target);
	memset(temp2, 0, sizeof(temp2));
	
	// retains some of the previously discovered if there was a call with use_parents=true before
	if(!use_parents)
	{
		for(i=4; i<6; i++)
		{
			int c=0;
			const lite_vector& vc = g.cliques[i];

			for(int j=0; j<vc.size(); j++)
			{
				if(!(vc[j] & setmask).zero())
					continue;
				int lb=vc[j].leading_bit();
	//			if(lb<base-1)
	//				continue;
				int tr=vc[j].trailing_bit();
				temp[tr+lb*target]|=1<<(i-4);
			}
		}
	}
	
	for(i=0; i<6; i++)
		part_cliques[i].clear();
	
	int stop_point = use_parents ? 0 : g.parent_cliques[0].size();
//	if(use_parents && g.parent_cliques[0].size()>2)
//		stop_point=g.parent_cliques[0].size()-2;
	buffer<int>& vmatch = g.vmatch;

	for(int d=g.parent_cliques[0].size(); d>=stop_point; d--)
	{
		//int d=g.parent_cliques[0].size();
		for(i=0; i<2; i++)
		{
			if(g.parent_cliques[i].size() != g.parent_cliques[0].size())
			{
				printf("?");
				exit(-1);
			}

			int c=0;
			/*
			if(g.parent_cliques[i].size()>0)
			{
				printf("%d\n", g.parent_cliques[i].size());
				exit(-1);
			}
			*/
			const lite_vector& vc = (d==g.parent_cliques[0].size()) ? g.cliques[i] : *g.parent_cliques[i][d];
			//printf("%d/%d %d %p\n", d, g.parent_cliques[0].size(), i, &vc);

			vmatch.clear();

			bigint ref_mask = (i==0) ? g.mask0 : g.mask;
			bigint dual_mask = (i==0) ? g.mask : g.mask0;

			for(int j=0; j<vc.size(); j++)
			{
				if(vc[j].popcnt != kmax[i]-1)
					continue;

				bigint x = vc[j];
				x.set_len(target);
				int last_bit = x.trailing_bit()+2;
				if(last_bit>target+1)
					last_bit=target+1;
				int array_index_base = max(base, last_bit);
				int shift = target-array_index_base+1;
				x = x << shift;
				int max_n = target;
				int min_n = base;
				vmatch.resize(0);
	
				for(int nn=min_n; nn<max_n; nn+=16)
				{
					int mask_array_index = nn-array_index_base;
					int b;
					if(mask_array_index < 0)
						b = batch_pattern_match(x, g.combo_ptrs[i]+mask_array_index);
					else
					{
						b = batch_pattern_match_2(x, g.combo_ptrs[i]+mask_array_index);
					}
					b &= ~vc[j].word(nn-1);
					b &= ~dual_mask.word(nn-1);
					if(b==0)
						continue;
					if(max_n<nn+16)
					{
						b &= (1<<(max_n-nn)) - 1;
					}
					size_t sizePrev=vmatch.size();
					vmatch.resize(sizePrev+_popcnt32(b));
					while(b!=0)
					{
						unsigned long pos;
						_BitScanForward(&pos, b);
						vmatch[sizePrev]=pos+nn;
						sizePrev++;
						b &= ~(1<<pos);
					}
				}
	
				for(k=0; k<vmatch.dim; k++)
				{
					for(int m=k+1; m<vmatch.dim; m++)
					{
						int a = vmatch[m], b=vmatch[k];
						if(a==b)
							continue;
						if(a<b)
							swap(a,b);
						int link3 = abs(a-b);
						int bit0 = ref_mask.bit(a-1);
						int bit1 = ref_mask.bit(b-1);

						int bit01 = ref_mask.bit(abs(a-b)-1);
						if(bit01==1 && b>=base-1 && a<target)
						{
							if(bit0+bit1==2)
								return false;
							if(bit0+bit1==1)
								temp2[(bit0 ? b-1 : a-1)]|=1<<i;
							else
								temp[(a-1)+(b-1)*target]|=1<<i;
						}
					}
				}
			}
		}
	}
	
	for(int d=0; d<g.parent_cliques[0].size(); d++)
	//int d=g.parent_cliques[0].size();
	for(i=4; i<6; i++)
	{
		int c=0;
		const lite_vector& vc = *g.parent_cliques[d][i];

		for(int j=0; j<vc.size(); j++)
		{
			if(!(vc[j] & setmask).zero())
				continue;
			int lb=vc[j].leading_bit();
			if(lb<base-1)
				continue;
			int tr=vc[j].trailing_bit();
			temp[tr+lb*target]|=1<<(i-4);
		}
	}
	
	for(i=base-1; i<=target; i++)
	{
		if(temp2[i]!=0)
		{
			bigint x;
			x.clear();
			x.set(i);
			if(temp2[i]&1)
				part_cliques[0].push_back(x);
			if(temp2[i]&2)
				part_cliques[1].push_back(x);
		}
	}
	for(i=base-1; i<target; i++)
	{
		for(int j=i+1; j<target; j++)
		{
			if(temp[j+i*target]!=0)
			{
				bigint x;
				x.clear();
				x.set(i);
				x.set(j);
				if(temp[j+i*target]&1)
					part_cliques[2].push_back(x);
				if(temp[j+i*target]&2)
					part_cliques[3].push_back(x);
			}
		}
	}
	return true;
}

bool link_pair_search(int target, int kmax[2], bigint* combo_ptrs[2], graph& g, buffer<pair<int,int> >& new_bits, bool debug, bool use_parents)
{
	g.check_first_unset_bit();
	if(target < g.first_unset_bit+10)
	{
		g.set_extensions = false;
		return true;
	}

//	lite_vector part_cliques[8];
	lite_vector* part_cliques = &g.cliques[2];
	int i;
	if(use_parents)
	{
		for(i=0; i<6; i++)
			part_cliques[i].clear();
	}
	if(!construct_incomplete_cliques_simple(part_cliques, g, target, kmax, use_parents))
		return false;
	int base = g.first_unset_bit-1;
	if(debug)
	{
		printf("g.first_unset_bit %d\n", base);
		g.mask0.printbits();
		g.mask.printbits();
	}
	int j, k, m;
	new_bits.resize(0);
	bigint newmask[2]={g.mask0, g.mask};

	int max_dim=target-base;

	for(i=0; i<2; i++)
	{
		bigint mask = (i==0) ? g.mask0 : g.mask;
		bigint dual_mask = (i==0) ? g.mask : g.mask0;
		const lite_vector& v = part_cliques[i];
		for(j=0; j<v.size(); j++)
		{
			//bigint x = andnot(g.cliques[2+i][j], mask);
			bigint x = v[j];
			if(popcnt(x)!=1)
				continue;
			int pos = x.trailing_bit();
			if(newmask[i].bit(pos))
				return false;
			if(!newmask[1-i].bit(pos))
			{
				if(debug)
					printf("%d => %d\n", pos, 1-i);
				newmask[1-i].set(pos);
				new_bits.push_back(make_pair(i, pos+1));
			}
		}
	}

	bigint extensions[4][MAX_LEN*64];
	/*
	extensions[0] = new bigint[max_dim];
	extensions[1] = new bigint[max_dim];
	extensions[2] = new bigint[max_dim];
	extensions[3] = new bigint[max_dim];
	*/

	for(i=0; i<4; i++)
		for(j=0; j<max_dim; j++)
		{
			extensions[i][j].clear();
			extensions[i][j].set_len(target);
		}
	for(j=0; j<max_dim; j++)
	{
		extensions[0][j].set(j+base);
		extensions[3][j].set(j+base);
	}
	for(j=0; j<new_bits.size(); j++)
	{
		int color=1-new_bits[j].first;
		int pos=new_bits[j].second-1;
		if(pos < base)
			printf("ERROR: new_bits (%d,%d), base %d\n", new_bits[j].first, new_bits[j].second, base);
		for(i=0; i<max_dim; i++)
		{
			extensions[color][i].set(pos);
			extensions[color+2][i].set(pos);
		}
	}

	while(true)
	{
		int prev_new_bits_size=new_bits.size();
		int changes=0;
		for(i=0; i<4; i++)
		{
			int c = i&1;
			bigint mask = (c==0) ? g.mask0 : g.mask;
			bigint dual_mask = (c==0) ? g.mask : g.mask0;

			//for(int d=0; d<=g.parent_cliques[0].size(); d++)		
			{
				//const lite_vector& v = (d==g.parent_cliques[0].size()) ? part_cliques[i+2] : *g.parent_cliques[i+2][d];
				const lite_vector& v = part_cliques[i+2];
				for(j=0; j<v.size(); j++)
				{
					//bigint rule = andnot(g.cliques[4+i][j], mask);
					bigint rule = v[j];
					if(rule.leading_bit() < base)
						continue;
					if(popcnt(rule)>=4)
					{
						printf("Erroneous rule\n");
						rule.printbits();
						mask.printbits();
						exit(-1);
					}
	//				rule.set_len(target);
	//				rule.printbits();

					bigint rule_save = rule;
					while(!rule.zero())
					{
						int m= rule.trailing_bit();
						rule.unset(m);
						m-=base;
						bigint rem = andnot(rule_save, extensions[c*3][m]);
						int p = popcnt(rem);
						if(p==0)
						{
							if(!newmask[1-c].bit(m+base)) // complete clique down this road
							{
								if(debug)
									printf("%d => %d\n", m+base, 1-c);
								newmask[1-c].set(m+base);
								new_bits.push_back(make_pair(c, m+base+1));
	//								changes++;
							}
						}
						if(p==1) // forced bit
						{
							int pos = rem.trailing_bit();
							if(!extensions[1-c+c*2][m].bit(pos))
							{
								extensions[1-c+c*2][m].set(pos);
								changes++;
							}
						}
					}
				}
			}
		}
		for(i=0; i<2; i++)
		{
			for(k=0; k<2; k++)
			{
				for(m=0; m<max_dim; m++)
				{
					int sum=popcnt(extensions[i*2][m])+popcnt(extensions[i*2+1][m]);
					bigint x = extensions[i*2+k][m];
					x.unset(m+base);
					while(!x.zero())
					{
						int pos = x.trailing_bit();
						if(pos-base<0 || pos-base>=target)
						{
							printf("ERROR: extensions[%d][%d] trailing bit %d , base %d, max_dim %d\n", i*2+k, m, pos, base, max_dim);
							extensions[i*2+k][m].printbits();
							//exit(-1);
						}
						extensions[i*2][m] |= extensions[k*2][pos-base];
						extensions[i*2+1][m] |= extensions[k*2+1][pos-base];
						x.unset(pos);
					}
					int sum2=popcnt(extensions[i*2][m])+popcnt(extensions[i*2+1][m]);
					if(sum!=sum2)
						changes++;
				}
			}
		}
#if 1
		for(i=0; i<2; i++)
		{
			for(m=0; m<max_dim; m++)
			{
				bigint x = extensions[i*2][m] & extensions[i*2+1][m];
				if(!x.zero())
				{
					if(!newmask[1-i].bit(m+base))
					{
						newmask[1-i].set(m+base);
						new_bits.push_back(make_pair(i, m+base+1));
					}
				}
			}
		}
#endif
		if(new_bits.size() != prev_new_bits_size)
		{
			for(i=prev_new_bits_size; i<new_bits.size(); i++)
			{
				int color = 1-new_bits[i].first;
				int pos = new_bits[i].second-1;
				for(j=0; j<max_dim; j++)
				{
					if(j+base==pos)
						continue;
					extensions[0][j] |= extensions[color*2+0][pos-base];
					extensions[1][j] |= extensions[color*2+1][pos-base];
					extensions[2][j] |= extensions[color*2+0][pos-base];
					extensions[3][j] |= extensions[color*2+1][pos-base];
				}
			}
		}
		else
		{
			if(changes==0)
				break;
		}
		if(popcnt(newmask[0], newmask[1])>0)
		{
			if(debug)
				printf("Full\n");
			return false;
		}
	}

	g.set_extensions = false;
	/*
	g.set_extensions = true;
	g.max_set_ext = base+max_dim-1;
	for(i=0; i<base; i++)
	{
		g.extensions[0][i].clear();
		g.extensions[1][i].clear();
		g.extensions[2][i].clear();
		g.extensions[3][i].clear();
	}
	for(i=base; i<base+max_dim; i++)
	{
		g.extensions[0][i] = extensions[0][i-base];
		g.extensions[1][i] = extensions[1][i-base];
		g.extensions[2][i] = extensions[2][i-base];
		g.extensions[3][i] = extensions[3][i-base];
	}
	*/
	double best_save=10, second_best_save=10;
	g.best_next_bit=-1;
	g.second_best_next_bit=-1;
	g.save_base=base;

	for(i=0; i<max_dim; i++)
	{
		if(!newmask[0].bit(base+i) && !newmask[1].bit(base+i))
		{
			int n0=0, n1=0;
			n0 = popcnt(extensions[0][i])+popcnt(extensions[1][i])-1;
			n1 = popcnt(extensions[2][i])+popcnt(extensions[3][i])-1;
			if(n0>63)
				n0=63;
			if(n1>63)
				n1=63;
			double save = 1./(One<<n0) + 1./(One<<n1);
			//g.save_rates[i] = save;
			//printf("%d  %d %d\n", i+base-2, n0, n1);
			if(save<best_save)
			{
				second_best_save = best_save;
				g.second_best_next_bit = g.best_next_bit;
				best_save = save;
				g.best_next_bit = i+base;
			}
			else if(save<second_best_save)
			{
				second_best_save = save;
				g.second_best_next_bit = i+base;
			}
		}
	}
	/*
	delete[] extensions[0];
	delete[] extensions[1];
	delete[] extensions[2]; 
	delete[] extensions[3];
	*/
	return true;
}

bool graph::set_bits(int color, buffer<int>& new_bits, int kmax[2], int target, bigint* badmask)
{
	for(int pos=0; pos<new_bits.size(); pos++)
	{
		int new_bit = new_bits[pos]-1;
//		printf("graph::set_bits(%d,%d)\n", color, new_bit);
		if(new_bit<0 || new_bit>=target+1)
			printf("?");
		if(color)
			mask.set(new_bit);
		else
			mask0.set(new_bit);
		if(color ? mask0.bit(new_bit) : mask.bit(new_bit))
			return false;
	}
	check_first_unset_bit();
	bigint ref_mask = color ? mask : mask0;

	//lite_vector& vOut = cliques[color];
	bigint target_mask[MAX_LEN*64];
	for(int pos=0; pos<new_bits.size(); pos++)
		invert(target_mask[pos], color?mask:mask0, new_bits[pos]-1);

	bigint new_mask;
	new_mask.clear();
	new_mask.set_len(MAX_LEN*64);

	for(int j=0; j<parent_cliques[color].size()+1; j++)
	{
		const lite_vector& vIn = (j==parent_cliques[color].size()) ? cliques[color] : *(parent_cliques[color][j]);
		size_t i=0;
		int dim = vIn.size();
		bool self_write = (j==parent_cliques[color].size());
			
		for(i=0; i+16<=dim; i+=16)
		{
			for(int pos=0; pos<new_bits.size(); pos++)
			{
				int x = inv_batch_pattern_match(target_mask[pos], &vIn[i]);
				if(x==0)
					continue;
				int new_bit = new_bits[pos]-1;
				for(int di=0; di<16; di++)
				//if(pattern_match(v[i], target_mask))
				if((x>>di)&1)
				{
					if(vIn[i+di].bit(new_bit))
						continue;
					lite_vector& vOut = cliques[color];
					size_t d = vOut.size();
					vOut.resize(d+1);
					bigint& clNew = vOut[d];
					clNew = vIn[i+di];
					clNew.set(new_bit);
					
//					if(popcnt(clNew, new_mask)>=2 && clNew.trailing_bit()>new_bit)
//						continue;

					if(clNew.len<new_bit+1)
						clNew.set_len(new_bit+1);
					clNew.popcnt++;
					clique_counts[color][clNew.popcnt-1]++;

					if(clNew.popcnt-1==kmax[color])
					{
						if(badmask!=0)
							*badmask=clNew;
						return false;					
					}
					int tr = vIn[i+di].trailing_bit();
					if(new_bit<tr)
					{
						bigint clExt2 = vIn[i+di];
						clExt2.unset(tr);
						clExt2.set(new_bit);
						clExt2.set_len(tr);
						bigint clExt;
						invert(clExt, clExt2);
						clExt.set(tr);
						clExt.set_len(vIn[i+di].len);
						clExt.popcnt++;
						if(!clExt.bit(new_bit))
						{
							vOut.push_back(clExt);
							clique_counts[color][clExt.popcnt-1]++;
						}						
					}
					if(self_write)
						dim = vOut.size();
				}
			}
		}
		for(; i<dim; i++)
		{
			for(int pos=0; pos<new_bits.size(); pos++)
			{
				bool match=pattern_match(vIn[i], target_mask[pos]);
				if(match)
				{
					int new_bit = new_bits[pos]-1;
					if(vIn[i].bit(new_bit))
						continue;
					lite_vector& vOut = cliques[color];
					size_t d = vOut.size();
					vOut.resize(d+1);
					bigint& clNew = vOut[d];
					clNew = vIn[i];
					clNew.set(new_bit);
//					if(popcnt(clNew, new_mask)>=2 && clNew.trailing_bit()>new_bit)
//						continue;

					if(clNew.len<new_bit+1)
						clNew.set_len(new_bit+1);
					clNew.popcnt++;
					clique_counts[color][clNew.popcnt-1]++;
					if(clNew.popcnt-1==kmax[color])
					{
						if(badmask!=0)
							*badmask=clNew;
						return false;
					}
					int tr = vIn[i].trailing_bit();
					if(new_bit<tr)
					{
						bigint clExt;
						clExt.clear();
						clExt.set_len(vIn[i].len);
						for(int k=0; k<tr; k++)
							if(clNew.bit(k))
								clExt.set(tr-k-1);
						clExt.set(tr);
						clExt.popcnt=clNew.popcnt;
						if(!clExt.bit(new_bit))
						{
							vOut.push_back(clExt);
							clique_counts[color][clExt.popcnt-1]++;
						}
					}					
					if(self_write)
						dim = vOut.size();
				}
			}
		}
	}

	for(int pos=0; pos<new_bits.size(); pos++)
	{
		int new_bit = new_bits[pos]-1;
		new_mask.set(new_bit);
		bigint x;
		x.clear();
		x.set_len(new_bit+1);
		x.set(new_bit);
		x.popcnt=1;
		cliques[color].push_back(x);
		clique_counts[color][0]++;
		if(clique_counts[0][kmax[0]]>0 || clique_counts[1][kmax[1]]>0)
			return false;
		//new_clique_search(kmax, target, 
	}

	return true;
}

int graph::new_clique_search(int kmax[2], int target, buffer<int> new_marks[2], bigint* badmask)
{
	bool fail = false;
	for(int color=0; color<2; color++)
	{
		if(new_marks[color].size()==0)
			continue;
		
		else if(new_marks[color].size()==1)
		{
			if(!set_bit(color, new_marks[color][0]-1, kmax, true, badmask))
				return target-10;
		}
		else
		{
#if MAX_LEN != 2
			if(!set_bits(color, new_marks[color], kmax, target, badmask))
					return target-10;
#else
			for(int edit=0; edit<new_marks[color].size(); edit++)
			{
				int new_bit = new_marks[color][edit]-1;
				if(new_bit<0)
					continue;
	//			if(n>=52)
				//printf("new_clique_search: setting bit %d to %d\n", new_bit, color);
				if(!set_bit(color, new_bit, kmax, true, badmask))
					return target-10;
			}
#endif
		}
	}
	if(fail)
		return target-10;
	return -1;
}

void graph::forced_bit_search2(int kmax[2], int target, int color, int new_bit, buffer<int> new_marks[2], int sizePrev, bool debug)
{
	int j;

	bigint new_bits[2];
	new_bits[0].clear();
	new_bits[1].clear();

	int i=color;
	const lite_vector& vc = cliques[i];
	bigint& ref_mask = i ? mask : mask0;
	bigint& dual_mask = i ? mask0 : mask;

	for(j=sizePrev; j<vc.size(); j++)
	{
		if(vc[j].popcnt != kmax[color])
			continue;
		if(!(vc[j].bit(new_bit)))
			continue;

		bigint x = vc[j];
		x.set_len(target);
		int last_bit = x.trailing_bit()+2;
		if(last_bit>target+1)
			last_bit=target+1;
		int array_index_base = max(first_unset_bit, last_bit);
		int shift = target-array_index_base+1;
		x = x << shift;
		int max_n = target;
		int min_n = first_unset_bit;
		for(int nn=min_n; nn<max_n; nn+=16)
		{
			int mask_array_index = nn-array_index_base;
			if(mask_array_index < -search_range || mask_array_index+15 >= MAX_LEN*64)
			{
				printf("mask_array_index out of bounds (%d), first_unset_bit %d, last_bit %d, search_range %d\n", mask_array_index, first_unset_bit, last_bit, search_range);
				exit(-1);
			}
			int b;
			if(mask_array_index < 0)
				b = batch_pattern_match(x, combo_ptrs[i]+mask_array_index);
			else
			{
				b = batch_pattern_match_2(x, combo_ptrs[i]+mask_array_index);
			}
			if(b==0)
				continue;
			for(int n=nn; n<nn+16 && n<max_n; n++)
			{
				bool match = (b>>(n-nn)) & 1;
				if(match)
				{
					if(!vc[j].bit(n-1) && !dual_mask.bit(n-1) && !new_bits[i].bit(n-1))
					{
						if(debug)
						{
							printf("%d  ", n-1);
							vc[j].printbits();
						}
						//printf(" forced_bit_search: set bit %d to %d (i %d, cl %d, j %d)\n", n, 1-i, i, cl, j);
						//verify_extension(vc[j], n, target, ref_mask);
						new_bits[i].set(n-1);
						//update_shifted_masks(i, target, n);
						new_marks[1-i].push_back(n);
					}
				}
			}
		}
	}
}

void graph::forced_bit_search(int kmax[2], int target, buffer<int> new_marks[2], bool use_parent_cliques, int color_mask, int minpos0, int minpos1, bool debug)
{
	int j;

	bigint new_bits[2];
	new_bits[0].clear();
	new_bits[1].clear();
	
	for(int i=0; i<2; i++)
	{
		if(!(color_mask & (1<<i)))
			continue;
		for(int cl=use_parent_cliques ? 0 : parent_cliques[i].size(); cl<parent_cliques[i].size()+1; cl++)
		{			
			const lite_vector& vc = (cl<parent_cliques[i].size()) ? *parent_cliques[i][cl] : cliques[i];
			bigint& ref_mask = i ? mask : mask0;
			bigint& dual_mask = i ? mask0 : mask;
			for(j=(i==0 ? minpos0 : minpos1); j<vc.size(); j++)
			{
				if(vc[j].popcnt != kmax[i])
					continue;
				bigint x = vc[j];
				x.set_len(target);
				int last_bit = x.trailing_bit()+2;
				if(last_bit>target+1)
					last_bit=target+1;
				int array_index_base = max(first_unset_bit, last_bit);
				int shift = target-array_index_base+1;
				x = x << shift;
				int max_n = target;
				int min_n = first_unset_bit;
				for(int nn=min_n; nn<max_n; nn+=16)
				{
					int mask_array_index = nn-array_index_base;
					if(mask_array_index < -search_range || mask_array_index+15 >= MAX_LEN*64)
					{
						printf("mask_array_index out of bounds (%d), first_unset_bit %d, last_bit %d, search_range %d\n", mask_array_index, first_unset_bit, last_bit, search_range);
						exit(-1);
					}
					int b;
					if(mask_array_index < 0)
						b = batch_pattern_match(x, combo_ptrs[i]+mask_array_index);
					else
					{
						b = batch_pattern_match_2(x, combo_ptrs[i]+mask_array_index);
					}
					if(b==0)
						continue;
					for(int n=nn; n<nn+16 && n<max_n; n++)
					{
						bool match = (b>>(n-nn)) & 1;
						if(match)
						{
							if(!vc[j].bit(n-1) && !dual_mask.bit(n-1) && !new_bits[i].bit(n-1))
							{
								if(debug)
								{
									printf("%d  ", n-1);
									vc[j].printbits();
								}
								//printf(" forced_bit_search: set bit %d to %d (i %d, cl %d, j %d)\n", n, 1-i, i, cl, j);
								//verify_extension(vc[j], n, target, ref_mask);
								new_bits[i].set(n-1);
								//update_shifted_masks(i, target, n);
								new_marks[1-i].push_back(n);
							}
						}
					}
				}
			}
		}
	}
}

bool verify_clique(const bigint& x, int target, const bigint& mask)
{
	int positions[MAX_DEPTH];
	int c;
	bigint_to_positions(positions, c, target+1, x);
	if(c != x.popcnt)
	{
		printf("Error: popcnt mismatch (%d / %d)\n", c, x.popcnt);
		x.print();
		for(int k=0; k<c; k++)
			printf("%d\n", positions[k]);
		return false;
	}
	int a, b;
	for(a=0; a<c; a++)
	{
		//if(positions[a]<0 || positions[a]>=
		if(!mask.bit(positions[a]-1))
		{
			printf("verify_clique error #1: link %d..%d wrong color\n", 0, positions[a]);
			mask.print();
			x.print();
			for(int k=0; k<c; k++)
				printf("%d  %d\n", positions[k], mask.bit(positions[k]-1));
		return false;
		}
	}

	for(a=0; a<c; a++)
		for(b=a+1; b<c; b++)
		{
			if(!mask.bit((positions[b]-positions[a])-1))
			{
				printf("verify_clique error #2: link %d..%d wrong color\n", positions[a], positions[b]);
				mask.print();
				x.print();
				for(int k=0; k<c; k++)
					printf("%d  %d\n", positions[k], mask.bit(positions[k]-1));
				for(int k=0; k<c; k++)
					for(int m=k+1; m<c; m++)
						printf("%d..%d  %d\n", positions[k], positions[m], mask.bit((positions[m]-positions[k])-1));
				return false;
			}
		}
	return true;
}

void verify_extension(const bigint& x, int n, int target, const bigint& mask)
{
	int positions[MAX_DEPTH];
	int c;
	bigint_to_positions(positions, c, target+1, x);
	if(c != x.popcnt)
	{
		printf("Error: popcnt mismatch (%d / %d)\n", c, x.popcnt);
		x.print();
		for(int k=0; k<c; k++)
			printf("%d\n", positions[k]);
		printf("n=%d\n", n);
		exit(-1);
	}
	positions[c++]=n;
	int a, b;

	for(a=0; a<c-1; a++)
	{
		//printf("0..%d  %d\n", positions[a], mask.bit(target-positions[a]));
		if(!mask.bit(positions[a]-1))
		{
			printf("Verify_extension error: link %d..%d wrong color\n", 0, positions[a]);
			printf("n=%d\n", n);
			printf("mask ");mask.print();
			printf("x    ");x.print();
			for(int k=0; k<c; k++)
				printf("%d  %d\n", positions[k], mask.bit(positions[k]-1));
			exit(-1);
		}
	}

	for(a=0; a<c; a++)
		for(b=a+1; b<c; b++)
		{
			//printf("%d..%d  %d\n", positions[a], positions[b], mask.bit(target-abs(positions[b]-positions[a])));
			if(!mask.bit(abs(positions[b]-positions[a])-1))
			{
				printf("Verify_extension error: link %d..%d wrong color\n", positions[a], positions[b]);
				printf("n=%d\n", n);
				mask.print();
				x.print();
				for(int k=0; k<c; k++)
					printf("%d  %d\n", positions[k], mask.bit(positions[k]-1));
				exit(-1);
			}
		}
}

void verify_extension(const bigint& x, int gap, int n, int target, const bigint& mask)
{
	int positions[MAX_DEPTH];
	int c;
	bigint_to_positions(positions, c, target+1, x);
	if(c != x.popcnt)
	{
		printf("Error: popcnt mismatch (%d / %d)\n", c, x.popcnt);
		x.print();
		for(int k=0; k<c; k++)
			printf("%d\n", positions[k]);
		printf("n=%d\n", n);
		exit(-1);
	}
	positions[c++]=n;
	int a, b;

	for(a=0; a<c; a++)
	{
		//printf("0..%d  %d\n", positions[a], mask.bit(target-positions[a]));
		if(positions[a]-1==gap)
			continue;
		if(!mask.bit(positions[a]-1))
		{
			printf("Verify_extension error: link %d..%d wrong color\n", 0, positions[a]);
			printf("n=%d\n", n);
			printf("mask ");mask.print();
			printf("x    ");x.print();
			for(int k=0; k<c; k++)
				printf("%d  %d\n", positions[k], mask.bit(positions[k]-1));
			exit(-1);
		}
	}

	for(a=0; a<c; a++)
		for(b=a+1; b<c; b++)
		{
			//printf("%d..%d  %d\n", positions[a], positions[b], mask.bit(target-abs(positions[b]-positions[a])));
			if(abs(positions[b]-positions[a])-1==gap)
				continue;
			if(!mask.bit(abs(positions[b]-positions[a])-1))
			{
				printf("Verify_extension error: link %d..%d wrong color\n", positions[a], positions[b]);
				printf("n=%d\n", n);
				printf("gap=%d\n", gap);
				mask.print();
				x.print();
				for(int k=0; k<c; k++)
					for(int m=k+1; m<c; m++)
						printf("%d..%d  %d\n", positions[k], positions[m], mask.bit(abs(positions[m]-positions[k])-1));
				exit(-1);
			}
		}
}

void verify_extension_2x(const bigint& x, int n1, int n2, int target, const bigint& mask)
{
	int positions[MAX_DEPTH];
	int c;
	bigint_to_positions(positions, c, target+1, x);
	if(c != x.popcnt)
	{
		printf("Error: popcnt mismatch (%d / %d)\n", c, x.popcnt);
		x.print();
		for(int k=0; k<c; k++)
			printf("%d\n", positions[k]);
		printf("n1=%d n2=%d\n", n1, n2);
		exit(-1);
	}
	positions[c++]=min(n1, n2);
	positions[c++]=max(n1, n2);
	int a, b;

	for(a=0; a<c-2; a++)
	{
		if(!mask.bit(positions[a]-1))
		{
			printf("Verify_extension_2x error: link %d..%d wrong color\n", 0, positions[a]);
			printf("n1=%d, n2=%d\n", n1, n2);
			printf("mask ");mask.print();
			printf("x    ");x.print();
			for(int k=0; k<c; k++)
				printf("%d  %d\n", positions[k], mask.bit(positions[k]-1));
			exit(-1);
		}
	}

	for(a=0; a<c; a++)
		for(b=a+1; b<c; b++)
		{
			if(!mask.bit(abs(positions[b]-positions[a])-1))
			{
				printf("Verify_extension_2x error: link %d..%d wrong color\n", positions[a], positions[b]);
				printf("n1=%d, n2=%d\n", n1, n2);
				mask.print();
				x.print();
				for(int k=0; k<c-2; k++)
					printf("%d  %d\n", positions[k], mask.bit(positions[k]-1));
				for(int k=0; k<c; k++)
					for(int m=k+1; m<c; m++)
						printf("%d..%d  %d\n", positions[k], positions[m], mask.bit(abs(positions[m]-positions[k])-1));
				exit(-1);
			}
		}
}

void verify_extension_2x(const bigint& x, int idx1, int idx2, int n1, int n2, int target, const bigint& mask)
{
	int positions[MAX_DEPTH];
	int c;
	bigint_to_positions(positions, c, target+1, x);
	if(c != x.popcnt)
	{
		printf("Error: popcnt mismatch (%d / %d)\n", c, x.popcnt);
		x.print();
		for(int k=0; k<c; k++)
			printf("%d\n", positions[k]);
		printf("n1=%d n2=%d\n", n1, n2);
		exit(-1);
	}
	positions[c++]=min(idx1, idx2);
	positions[c++]=max(idx1, idx2);
	int a, b;

	for(a=0; a<c; a++)
	{
		if(positions[a]==n1 || positions[a]==n2)
			continue;
		if(!mask.bit(positions[a]-1))
		{
			printf("Verify_extension_2x error: link %d..%d wrong color\n", 0, positions[a]);
			printf("n1=%d, n2=%d\n", n1, n2);
			printf("nidx1=%d, nidx2=%d\n", idx1, idx2);
			printf("mask ");mask.print();
			printf("x    ");x.print();
			for(int k=0; k<c; k++)
				printf("%d  %d\n", positions[k], mask.bit(positions[k]-1));
			exit(-1);
		}
	}

	for(a=0; a<c; a++)
		for(b=a+1; b<c; b++)
		{
			if(abs(positions[b]-positions[a])==n1 || abs(positions[b]-positions[a])==n2)
				continue;
			if(!mask.bit(abs(positions[b]-positions[a])-1))
			{
				printf("Verify_extension_2x error: link %d..%d wrong color\n", positions[a], positions[b]);
				printf("n1=%d, n2=%d\n", n1, n2);
				printf("nidx1=%d, nidx2=%d\n", idx1, idx2);
				mask.print();
				x.print();
				for(int k=0; k<c-2; k++)
					printf("%d  %d\n", positions[k], mask.bit(positions[k]-1));
				printf("%d  %d (ignore)\n", positions[c-2],mask.bit(positions[c-2]-1));
				printf("%d  %d (ignore)\n", positions[c-1],mask.bit(positions[c-1]-1));
				for(int k=0; k<c; k++)
					for(int m=k+1; m<c; m++)
						printf("%d..%d  %d\n", positions[k], positions[m], mask.bit(abs(positions[m]-positions[k])-1));
				exit(-1);
			}
		}
}


int count_triangles(const bigint& x)
{
	int count=0;
	for(int i=0; i<x.len; i++)
	{
		for(int j=i+1; j<x.len; j++)
		{
			int c = x.bit(i) + x.bit(j) + x.bit(j-i-1);
			if(c==3)
				count++;
		}
	}
	return count;
}

int count_c4s(const bigint& mask)
{
	int count=0;
	int i, j, k;
	for(i=0; i<mask.len; i++)
	{
		for(j=i+1; j<mask.len; j++)
		{
			for(k=j+1; k<mask.len; k++)
			{ 
				int c = mask.bit(i) + mask.bit(j) + mask.bit(k)
					+ mask.bit(j-i-1) + mask.bit(k-j-1) + mask.bit(k-i-1);
				if(c==6)
					count++;
			}
		}
	}
	return count;
}

int count_triangles(const bigint& mask, const lite_vector& v0, const buffer<const lite_vector*>& vv)
{
	int count=0;
	int i, j;
//	vector<bool> observed;
//	observed.resize(v.size());
//	for(i=0; i<v.size(); i++)
//		observed[i]=false;
	for(i=0; i<mask.len; i++)
	{
		for(j=i+1; j<mask.len; j++)
		{
				int c = mask.bit(i) + mask.bit(j) + mask.bit(j-i-1);
				if(c==3)
				{
					int hit=0;
					for(int n=0; n<vv.size()+1; n++)
					{
						const lite_vector& v = (n==vv.size()) ? v0 : *vv[n];
						for(int m=0; m<v.size(); m++)
						{
							if(v[m].popcnt != 2)
								continue;
							if(v[m].bit(i) && v[m].bit(j))
							{
								hit++;
//								observed[m]=true;
							}
						}
					}
					if(hit>1)
					{
						printf("ERROR: duplicate clique %d, %d\n", i, j);
					}
					if(hit==0)
					{
						printf("ERROR: missed clique %d, %d\n", i, j);
					}
//					count++;
				}
			}
		
	}
	return count;
}

int count_c4s(const bigint& mask, const lite_vector& v0, const buffer<const lite_vector*>& vv)
{
	int count=0;
	int i, j, k;
	vector<bool> observed;
//	observed.resize(v.size());
//	for(i=0; i<v.size(); i++)
//		observed[i]=false;
	for(i=0; i<mask.len; i++)
	{
		for(j=i+1; j<mask.len; j++)
		{
			for(k=j+1; k<mask.len; k++)
			{ 
				int c = mask.bit(i) + mask.bit(j) + mask.bit(k)
					+ mask.bit(j-i-1) + mask.bit(k-j-1) + mask.bit(k-i-1);
				if(c==6)
				{
					int hit=0;
					for(int n=0; n<vv.size()+1; n++)
					{
						const lite_vector& v = (n==vv.size()) ? v0 : *vv[n];
						for(int m=0; m<v.size(); m++)
						{
							if(v[m].popcnt != 3)
								continue;
							if(v[m].bit(i) && v[m].bit(j) && v[m].bit(k))
							{
								hit++;
//								observed[m]=true;
							}
						}
					}
					if(hit>1)
					{
						printf("ERROR: duplicate clique %d, %d, %d\n", i, j, k);
					}
					if(hit==0)
					{
						printf("ERROR: missed clique %d, %d, %d\n", i, j, k);
					}
//					count++;
				}
			}
		}
	}
	return count;
}

void validate_graph_build_native(const graph& g, const bigint2& mask)
{
	int count_checks[2][MAX_DEPTH];
	memset(count_checks, 0, sizeof(count_checks));
	int i, j, k;

	for(i=0; i<MAX_LEN; i++)
	{
		if(mask.first.bit(i) && mask.second.bit(i))
		{
			printf("%d %llx %llx\n", mask.first.len, mask.first.n[0], mask.first.n[1]);
			printf("%d %llx %llx\n", mask.second.len, mask.second.n[0], mask.second.n[1]);
			mask.first.print();
			mask.second.print();
			g.print_short();
			exit(-1);
		}
	}
	/*
	if(!(g.mask == mask))
	{
		printf("ERROR: graph mask does not match build mask\n");
		mask.print();
		g.print();
		exit(-1);
	}
	*/
	int n_ones = popcnt(mask.second);
	if(n_ones != g.clique_counts[1][0])
	{
		printf("ERROR: clique_counts[1][0] incorrect (%d != %d)\n", 
			n_ones, g.clique_counts[1][0]);
		printf("%d %llx %llx\n", mask.first.len, mask.first.n[0], mask.first.n[1]);
		printf("%d %llx %llx\n", mask.second.len, mask.second.n[0], mask.second.n[1]);
		mask.second.print();
		g.print();
		exit(-1);
	}
	int n_zeros = popcnt(mask.first);
	if(n_zeros != g.clique_counts[0][0])
	{
		printf("ERROR: clique_counts[0][0] incorrect (%d != %d)\n", 
			n_zeros, g.clique_counts[0][0]);
		g.print();
		exit(-1);
	}

	int n01 = count_triangles(mask.first);
	int n11 = count_triangles(mask.second);

	int n02 = count_c4s(mask.first);
	int n12 = count_c4s(mask.second);
#if 1
	if(n01!=g.clique_counts[0][1] || n11!=g.clique_counts[1][1])
	{
		printf("ERROR: clique_counts[0][1] or [1][1] incorrect (%d,%d vs %d,%d)\n",
			n01, n11, g.clique_counts[0][1], g.clique_counts[1][1]);
		count_triangles(mask.first, g.cliques[0], g.parent_cliques[0]);
		count_triangles(mask.second, g.cliques[1], g.parent_cliques[1]);
		mask.first.print();
		mask.second.print();
		g.print();
		exit(-1);
	}

	if(n02!=g.clique_counts[0][2] || n12!=g.clique_counts[1][2])
	{
		printf("ERROR: clique_counts[0][2] or [1][2] incorrect (%d,%d vs %d,%d)\n",
			n02, n12, g.clique_counts[0][2], g.clique_counts[1][2]);
		count_c4s(mask.first, g.cliques[0], g.parent_cliques[0]);
		count_c4s(mask.second, g.cliques[1], g.parent_cliques[1]);
		g.print();
		exit(-1);
	}
#endif

	for(i=0; i<2; i++)
	{
		buffer<const lite_vector*> pv;
		pv = g.parent_cliques[i];
		pv.push_back(&g.cliques[i]);
		for(k=0; k<pv.size(); k++)
		{
			const lite_vector& v = *pv[k];
			for(j=0; j<v.size(); j++)
			{
				int n = popcnt(v[j]);
				if(n!=v[j].popcnt)
				{
					printf("ERROR: %d, %d: popcnt incorrectly recorded (%d)\n", i, k, v[j].popcnt);
					v[j].print();
					exit(-1);
				}
				if(n==0 || n>=MAX_DEPTH)
				{
					printf("ERROR: invalid clique\n");
					v[j].print();
					exit(-1);
				}
				int n_match = popcnt(v[j], (i==0) ? mask.first : mask.second);
				if(n_match != n)
				{
					printf("ERROR: clique with incomplete/incorrect overlap\n");
					printf("i=%d j=%d n=%d\n", i, j, n);
					v[j].print();
					mask.first.print();
					mask.second.print();
					exit(-1);
				}
				count_checks[i][n-1]++;
			}
/*
			for(j=0; j<g.cliques[i].size(); j++)
			{
				for(k=j+1; k<g.cliques[i].size(); k++)
				{
					if(g.cliques[i][j] == g.cliques[i][k])
					{
						printf("ERROR: duplicate clique: %d = %d\n", j, k);
						g.cliques[i][j].print();
						exit(-1);
					}
				}
			}
*/
		}
	}
	if(memcmp(count_checks, g.clique_counts, sizeof(count_checks)))
	{
		printf("ERROR: count checks don't match\n");
		for(i=0; i<MAX_DEPTH; i++)
			printf("%d %d   %d %d\n", g.clique_counts[0][i], g.clique_counts[1][i], count_checks[0][i], count_checks[1][i]);
		g.print();
		exit(-1);
	}
	for(int color=0; color<2; color++)
		for(int n=0; n<g.cliques[color].size(); n++)
			if(!verify_clique(g.cliques[color][n], g.cliques[color][n].len+1, color?g.mask:g.mask0))
			{
				g.print_short();
				exit(-1);
			}
}

