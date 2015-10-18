#include "stdafx.h"
#include "graph.h"

#ifdef WIN32
#define align16 __declspec(align(16))
#else
#define align16 __attribute__(aligned(16))
#endif

bool g_high_kmax_mode = true;


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
	while(i<GRAPH_MAX_LEN)
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
	first_unset_bit=GRAPH_MAX_LEN*64+1;
	/*
	while(first_unset_bit<GRAPH_MAX_LEN*64 && (mask.bit(first_unset_bit) || mask0.bit(first_unset_bit)))
		first_unset_bit++;
	first_unset_bit++;
	*/
}

void graph::update_shifted_masks(int color, int target, int n)
{
	if(!known_shifted_ref_masks)
	{
		printf("ERROR: update() called on unset masks\n");
		exit(-1);
	}

	if(color==0 && n>=mask0.len)
	{
		printf("Error\n");
		exit(-1);
	}
	if(color==1 && n>=mask.len)
	{
		printf("Error\n");
		exit(-1);
	}

	const bigint& m = color ? mask : mask0;

	bigint temp;
	temp.clear();
	invert(temp, m);
	int shift = m.len-n-1;
	if(shift>0)
		shifted_ref_masks[n+1] = temp >> shift;
	else if(shift<0)
		shifted_ref_masks[n+1] = temp << (-shift);
	else
		shifted_ref_masks[n+1] = temp;

	shifted_ref_masks[n+1] |= m << (n+2);

//	printf("%d\n", shifted_ref_masks[1].bit(113));
	for(int i=1; i<target+2; i++)
	{
		if(m.bit(i-1))
		{
//			if(i==1)
//				printf("update_shifted_masks(%d,%d,%d): %d, %d on\n", color, target, n,  i+n+1, i-n-1);
			//shifted_ref_masks[i].clear();
			if(i+n+1<GRAPH_MAX_LEN*64)
				shifted_ref_masks[i].set(i+n+1);
			if(i-n-1>=0)
				shifted_ref_masks[i].set(i-n-1);
		}
	}
///	printf("%d\n", shifted_ref_masks[1].bit(113));
	check_first_unset_bit();
	
}

void graph::construct_shifted_masks(int target)
{
	known_shifted_ref_masks=true;
	shifted_ref_mask_len=target;
	search_range = 0;
	check_first_unset_bit();
#if 1
	for(int i=1; i<target+2; i++)
	{
		if(mask0.bit(i-1))
		{			
			shifted_ref_masks[i].clear();
			
			//if(i>0)
			{
				bigint temp;
				temp.clear();
				invert(temp, mask0);
				if(mask0.len-i>0)
					shifted_ref_masks[i] = temp >> (mask0.len-i);
				else if(mask0.len-i<0)
					shifted_ref_masks[i] = temp << (i-mask0.len);
				else
					shifted_ref_masks[i] = temp ;
			}
			
			// i=0: (b2) (b1) b0 b1 b2 
			// i=1: _ _ _ b0 b1 b2
			shifted_ref_masks[i] |= mask0 << (i+1);
		}
		if(mask.bit(i-1))
		{
			shifted_ref_masks[i] .clear();
			
			//if(i>0)
			{
				bigint temp;
				temp.clear();
				invert(temp, mask);
				if(mask.len-i>0)
					shifted_ref_masks[i] = temp >> (mask.len-i);
				else if(mask.len-i<0)
					shifted_ref_masks[i] = temp << (i-mask.len);
				else
					shifted_ref_masks[i] = temp ;
			}
			
			shifted_ref_masks[i] |= mask << (i+1);
		}
		// 2 1 0 _ 0 1 2
		// i=0: _ 0 1 2
	}

#endif
}

void graph::operator=(graph& g)
{
	n=g.n;
	mask0=g.mask0;
	mask=g.mask;
	first_unset_bit = g.first_unset_bit;
	search_range = g.search_range;
	known_shifted_ref_masks=false;

	memcpy(clique_counts, g.clique_counts, sizeof(clique_counts));
	for(int i=0; i<8; i++)
	{
		cliques[i].clear();
		parent_cliques[i] = g.parent_cliques[i];
		parent_cliques[i].push_back(&g.cliques[i]);
	}
	explored[0]=g.explored[0];
	explored[1]=g.explored[1];
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
	explored[0].clear();
	explored[1].clear();
	known_shifted_ref_masks=false;
	mask.clear();
	mask0.clear();

	implied_mask[0].clear();
	implied_mask[1].clear();
	//second_best_next_bit=-1;
	//save_base=-1;
}

uint64_t g_set_bit_count = 0;

int g_min_depth;


template <int MAX_LEN>
struct recurse_object
{
	bare_bigint<MAX_LEN> stack[MAX_DEPTH];
	int bits[MAX_DEPTH];
	bare_bigint<MAX_LEN> inverted_masks[2][MAX_LEN*64+1];
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

template <int MAX_LEN>
void recurse_object<MAX_LEN>::prepare(int* km, int c, bool fastmode)
{
	const bigint& in_mask = c ? _m1 : _m0;
	const bigint& dual_mask = c ? _m0 : _m1;
	
	ext_bigint<MAX_LEN> x;
	x = in_mask;
	x.set_len(lastbit+1);
	x.unset_from(lastbit+1);
	x.unset(new_bit);
	x &= inverted_masks[1-c][new_bit];
	stack[0] = x;
//	bits[0]=lastbit;
	this->kmax = km[c];
	this->color = c;
}

template <int MAX_LEN>
bool recurse_object<MAX_LEN>::unrolled_flat_search()
{
	lite_vector& v = cliques[color];
	size_t startpos = v.size();
	bare_bigint<MAX_LEN>* inv_masks = inverted_masks[1-color];
	int depth=0;

	const bool debug = false;//(new_bit==179 && color==1); //false;

	uint64_t ops=0;

	const bigint& in_mask = color ? _m1 : _m0;

	ext_bigint<MAX_LEN> v_curr_cl[MAX_DEPTH];
	v_curr_cl[0].clear();
	v_curr_cl[0].set_len(lastbit+1);
	v_curr_cl[0].set(new_bit);

	if(0>=nMinCliqueSize)
	{
		bigint current_clique;
		current_clique.clear();
		current_clique.set_len(lastbit+1);
		current_clique.set(new_bit);
		current_clique.popcnt=1;
		clique_counts[color][0]+=1;
		v.push_back(current_clique);
	}

	while(!stack[0].zero())
	{
		int b0 = stack[0].trailing_bit();
		int depth_cutoff = min((nMinCliqueSize+3)/2, (kmax-2+3)/2);
		int bit_cutoff = new_bit/2;

		bits[0]=b0;
		depth=1;
		if(b0 < bit_cutoff && 1 < depth_cutoff)
			break;
		store_and(stack[1], stack[0], inv_masks[b0]);
		//stack[1]=stack[0] & inv_masks[b0];
		bits[1]=b0-1;

		v_curr_cl[1]=v_curr_cl[0];
		v_curr_cl[1].set(b0);
		
		if(1>=nMinCliqueSize)
		{
			bigint current_clique;
			current_clique.clear();
			current_clique.set_len(lastbit+1);
			current_clique.set(new_bit);
			current_clique.set(b0);
			current_clique.popcnt=2;
			clique_counts[color][depth]+=1;
			v.push_back(current_clique);
		}
		
		if(2>=nMinCliqueSize)
		{
			uint64_t cls = popcnt(stack[1]);
			if(2==kmax && cls>0)
				return false;
			uint64_t pos = v.size();
			v.resize(pos+cls);
			ext_bigint<MAX_LEN> y = stack[1];
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
				ext_bigint<MAX_LEN> y = stack[1];
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
				ext_bigint<MAX_LEN> save = stack[depth];
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
					ext_bigint<MAX_LEN> current_clique=v_curr_cl[depth-1];
					current_clique.set(bits[depth-1]);
					current_clique|=save;
					current_clique.popcnt=nMinCliqueSize+1;
					clique_counts[color][nMinCliqueSize]++;
					if(nMinCliqueSize==kmax)
					{
						return false;
					}
					bigint b;
					b = current_clique;
					v.push_back(b);
				}
			}
			goto done2;
		}

		while(true)
		{
			int b;
			if(stack[depth].zero())
				goto next;
			b = stack[depth].trailing_bit();
			
			if(b < bit_cutoff && depth < depth_cutoff)
			{
				stack[depth].clear();
				goto next;
			}

			bits[depth]=b;
			depth++;

			//stack[depth]=stack[depth-1] & inv_masks[bits[depth-1]];
			//bare_bigint t = stack[depth-1], inv_masks[b]
			if(b<0 || b>GRAPH_MAX_LEN*64)
			{
//				printf("b=%d, depth=%d, new_bit=%d\n", b, depth, new_bit);
				printf("%d *", b);
			}
			store_and(stack[depth], stack[depth-1], inv_masks[b]);
			if(!stack[depth].zero())
			{
				int x = depth+popcnt(stack[depth]);
				if(x<=nMinCliqueSize)
				{
					if(x==nMinCliqueSize)
					{
						ext_bigint<MAX_LEN> save = stack[depth];
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
							ext_bigint<MAX_LEN> current_clique=v_curr_cl[depth-1];
							current_clique.set(bits[depth-1]);
							current_clique|=save;
							current_clique.popcnt=nMinCliqueSize+1;
							clique_counts[color][nMinCliqueSize]++;
							if(nMinCliqueSize==kmax)
							{
								return false;
							}
							bigint b;
							b = current_clique;
							v.push_back(b);
							
							ext_bigint<MAX_LEN> y = current_clique;
							int len = y.trailing_bit();
							ext_bigint<MAX_LEN> t;
							t.clear();
							y.set_len(len);
							invert(t, y);
							t.set(len);
							t.popcnt=y.popcnt;
							t.set_len(current_clique.len);
							b=t;
							v.push_back(b);
						}
					}
					if(depth==1)
						goto done2;
					depth--;
					stack[depth].unset(bits[depth]);
				}
				else
				{
//					bits[depth]=bits[depth-1]-1;
					v_curr_cl[depth]=v_curr_cl[depth-1];
					if(v_curr_cl[depth].bit(b))
						printf("?");
					v_curr_cl[depth].set(b);
						//stack[depth].clear();
					if(depth+1>=nMinCliqueSize)
					{
						ext_bigint<MAX_LEN> x=stack[depth];
						uint64_t cls = popcnt(x);
						if(cls!=0 && depth+1>=kmax)
							return false;

						uint64_t pos = v.size();
						v.resize(pos+cls);
						ext_bigint<MAX_LEN> y = x;
						ext_bigint<MAX_LEN> current_clique=v_curr_cl[depth];
						current_clique.popcnt=depth+1;

						for(uint64_t i=0; i<cls; i++)
						{
							bigint& t = v[i+pos];
							t = current_clique;
							t.popcnt+=1;
							int b = y.trailing_bit();
							t.set(b);
							y.unset(b);
						}
						clique_counts[color][depth+1]+=cls;
				
						if(new_bit < bits[0])
						{
							v.resize(pos+cls*2);
							ext_bigint<MAX_LEN> y = x;
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

						else if(depth+1==kmax-1)
						{
							v.resize(pos+cls*2);
							ext_bigint<MAX_LEN> y = x;
							for(uint64_t i=0; i<cls; i++)
							{
								bigint& t = v[i+pos+cls];
								t.clear();
								bigint y = v[i+pos];
								y.set_len(new_bit);
								invert(t, y);
								t.set(new_bit);
								t.popcnt=depth+2;
								t.set_len(new_bit+1);
							}
							clique_counts[color][depth+1]+=cls;
						}

					}
				}
			}
		next:
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

template<int MAX_LEN>
bool graph_set_bit_internal(graph& g, int color, int new_bit, int kmax[2], int target, bool test_only, int max_bit)
{
	recurse_object<MAX_LEN> obj;
	bigint& _m0 = g.mask0;
	bigint& _m1 = g.mask;
	const bigint& in_mask = color ? _m1 : _m0;
	obj.clique_counts[0] = & g.clique_counts[0][0];
	obj.clique_counts[1] = & g.clique_counts[1][0];
	obj.cliques = &g.cliques[0];
	int i;
	for(i=0; i<=max_bit; i++)
	{
		if(!new_bit && !in_mask.bit(i))
			continue;
		obj.inverted_masks[1-color][i].clear();
		invert(obj.inverted_masks[1-color][i], in_mask, i);
	}

	g.check_first_unset_bit();
	int nMinCliqueSize=kmax[color]-2;
	if(nMinCliqueSize<0)
		nMinCliqueSize=0;
	// Memory saving trick (triggered in cases where the potential number of cliques is so large that we're at risk of running out of memory)
	if (target>=200 && kmax[color]>=10)
	{
		nMinCliqueSize = kmax[color] - 1;
		if (nMinCliqueSize<0)
			nMinCliqueSize = 0;
	}
	if(test_only || new_bit>target)
		nMinCliqueSize=kmax[color];
	obj.nMaxUnset=0;
	obj._m0 = g.mask0;
	obj._m1 = g.mask;
	obj.nMinCliqueSize = nMinCliqueSize;
	obj.new_bit = new_bit;
	obj.lastbit = max_bit;

	obj.prepare(kmax, color, true);
	if(!obj.unrolled_flat_search())
		return false;
	return true;
}

bool graph::set_bit(int color, int new_bit, int kmax[2], int target, bool test_only)
{
	bigint& _m0 = mask0;
	bigint& _m1 = mask;
	if((color ? _m1 : _m0).bit(new_bit))
		return true;
	g_set_bit_count++;

//	bigint inverted_masks[GRAPH_MAX_LEN*64+1];
	int max_bit = (color ? _m1 : _m0).trailing_bit();
	if(color)
	{
		_m1.set(new_bit);
		if(new_bit>=_m1.len)
			_m1.set_len(new_bit+1);
	}
	else
	{
		_m0.set(new_bit);
		if(new_bit>=_m0.len)
			_m0.set_len(new_bit+1);
	}
	int i;
	if(new_bit > max_bit)
		max_bit = new_bit;

	int last_set_bit=0;
	if(!mask.zero())
		last_set_bit=mask.trailing_bit();
	if(!mask0.zero())
		last_set_bit=max(last_set_bit, mask0.trailing_bit());
	last_set_bit=max(last_set_bit, new_bit);
	last_set_bit=max(last_set_bit, target);

	bool ok;
	if(last_set_bit<128)
		ok = graph_set_bit_internal<2>(*this, color, new_bit, kmax, target, test_only, max_bit);
	else if(last_set_bit<256)
		ok = graph_set_bit_internal<4>(*this, color, new_bit, kmax, target, test_only, max_bit);
	else
		ok = graph_set_bit_internal<GRAPH_MAX_LEN>(*this, color, new_bit, kmax, target, test_only, max_bit);
	
	if(!ok)
		return false;
	if(clique_counts[color][kmax[color]]!=0)
		return false;

	update_shifted_masks(color, target, new_bit);
	return true;
}

template <int MAX_LEN>
bool build_graph_internal(graph& gOut, const bigint2& mask, int kmax[2], int nSetBits, bool do_minus2)
{
	recurse_object<MAX_LEN> obj;
	obj.clique_counts[0] = & gOut.clique_counts[0][0];
	obj.clique_counts[1] = & gOut.clique_counts[1][0];
	obj.cliques = &gOut.cliques[0];
	obj.inverted_masks[0][0].clear();
	obj.inverted_masks[1][0].clear();
	obj.nMaxUnset=0;

	if (kmax[0] > kmax[1])
	{
		for (int pos = 1; pos <= nSetBits; pos++)
		{
			int color;
			//		if(mask.first.bit(pos))
			//			color=0;
			//		else 
			if (mask.second.bit(pos))
				color = 1;
			else
			{
				continue;
			}
			obj.inverted_masks[0][pos].clear();
			ext_bigint<MAX_LEN> m;
			m = mask.second;
			invert(obj.inverted_masks[0][pos], m, pos);
			if (color)
				gOut.mask.set(pos);
			else
				gOut.mask0.set(pos);
			obj._m0 = gOut.mask0;
			obj._m1 = gOut.mask;

			obj.new_bit = pos;
			obj.lastbit = pos;
			obj.nMinCliqueSize = kmax[color] - (do_minus2 ? 2 : 1);
			if (obj.nMinCliqueSize < 0)
				obj.nMinCliqueSize = 0;

			obj.prepare(kmax, color, true);
			if (!obj.unrolled_flat_search())
				return false;
			if (gOut.clique_counts[color][kmax[color]] != 0)
				return false;
		}

		for (int pos = 1; pos <= nSetBits; pos++)
		{
			int color;
			if (mask.first.bit(pos))
				color = 0;
			//		else if (mask.second.bit(pos))
			//			color = 1;
			else
			{
				continue;
			}
			obj.inverted_masks[1][pos].clear();
			invert(obj.inverted_masks[1][pos], mask.first, pos);
			if (color)
				gOut.mask.set(pos);
			else
				gOut.mask0.set(pos);
			obj._m0 = gOut.mask0;
			obj._m1 = gOut.mask;

			obj.new_bit = pos;
			obj.lastbit = pos;
			obj.nMinCliqueSize = kmax[color] - (do_minus2 ? 2 : 1);
			if (obj.nMinCliqueSize < 1)
				obj.nMinCliqueSize = 1;
			obj.prepare(kmax, color, true);
			
			if (!obj.unrolled_flat_search())
				return false;
			if (gOut.clique_counts[color][kmax[color]] != 0)
				return false;
		}
	}
	else if(kmax[0] < kmax[1])
	{
		for (int pos = 1; pos <= nSetBits; pos++)
		{
			int color;
			if (mask.first.bit(pos))
				color = 0;
			//		else if (mask.second.bit(pos))
			//			color = 1;
			else
			{
				continue;
			}
			obj.inverted_masks[1][pos].clear();
			invert(obj.inverted_masks[1][pos], mask.first, pos);
			if (color)
				gOut.mask.set(pos);
			else
				gOut.mask0.set(pos);
			obj._m0 = gOut.mask0;
			obj._m1 = gOut.mask;

			obj.new_bit = pos;
			obj.lastbit = pos;
			obj.nMinCliqueSize = kmax[color] - (do_minus2 ? 2 : 1);
			if (obj.nMinCliqueSize < 1)
				obj.nMinCliqueSize = 1;
			obj.prepare(kmax, color, true);
			if (!obj.unrolled_flat_search())
			{
				//printf("Fail at bit %d\n", pos);
				return false;
			}
			//		gOut.n++;
			if (gOut.clique_counts[color][kmax[color]] != 0)
				return false;
		}

		for (int pos = 1; pos <= nSetBits; pos++)
		{
			int color;
			//		if(mask.first.bit(pos))
			//			color=0;
			//		else 
			if (mask.second.bit(pos))
				color = 1;
			else
			{
				continue;
			}
			obj.inverted_masks[0][pos].clear();
			invert(obj.inverted_masks[0][pos], mask.second, pos);
			if (color)
				gOut.mask.set(pos);
			else
				gOut.mask0.set(pos);
			obj._m0 = gOut.mask0;
			obj._m1 = gOut.mask;

			obj.new_bit = pos;
			obj.lastbit = pos;
			obj.nMinCliqueSize = kmax[color] - (do_minus2 ? 2 : 1);
			if (obj.nMinCliqueSize < 1)
				obj.nMinCliqueSize = 1;
			obj.prepare(kmax, color, true);
			if (!obj.unrolled_flat_search())
			{
				//printf("Fail at bit %d\n", pos);
				return false;
			}
			//		gOut.n++;
			if (gOut.clique_counts[color][kmax[color]] != 0)
				return false;
		}
	}
	else
	{
		for (int pos = 1; pos <= nSetBits; pos++)
		{
			int color;
			if (mask.first.bit(pos))
				color = 0;
			else if (mask.second.bit(pos))
				color = 1;
			else
			{
				continue;
			}
			obj.inverted_masks[1][pos].clear();
			invert(obj.inverted_masks[1][pos], mask.first, pos);
			obj.inverted_masks[0][pos].clear();
			invert(obj.inverted_masks[0][pos], mask.second, pos);
			if (color)
				gOut.mask.set(pos);
			else
				gOut.mask0.set(pos);
			obj._m0 = gOut.mask0;
			obj._m1 = gOut.mask;

			obj.new_bit = pos;
			obj.lastbit = pos;
			obj.nMinCliqueSize = kmax[color] - (do_minus2 ? 2 : 1);
			if (obj.nMinCliqueSize < 1)
				obj.nMinCliqueSize = 1;
			obj.prepare(kmax, color, true);
			if (!obj.unrolled_flat_search())
				return false;
			//		gOut.n++;
			if (gOut.clique_counts[color][kmax[color]] != 0)
				return false;
		}
	}
	return true;
}

bool build_graph(graph& gOut, const bigint2& mask, int kmax[2], bool do_minus2)
{
	gOut.mask0.clear();
	gOut.mask.clear();
	gOut.mask0.set_len(1);
	gOut.mask.set_len(1);
	if(mask.second.bit(0))
	{
		gOut.mask.set(0);
	}
	else
	{
		gOut.mask0.set(0);
	}
	gOut.explored[0].clear();
	gOut.explored[1].clear();
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

	int nSetBits=0;
	if(popcnt(mask.first)!=0)
		nSetBits=mask.first.trailing_bit();
	if(popcnt(mask.second)!=0)
		nSetBits=max(nSetBits, mask.second.trailing_bit());
	bool ok;
	if(nSetBits<128)
		ok = build_graph_internal<2>(gOut, mask, kmax, nSetBits, do_minus2);
	else if(nSetBits<256)
		ok = build_graph_internal<4>(gOut, mask, kmax, nSetBits, do_minus2);
	else
		ok = build_graph_internal<GRAPH_MAX_LEN>(gOut, mask, kmax, nSetBits, do_minus2);
	if(!ok)
		return false;
/*
	buffer<int> vLink[2];
	vLink[0].push_back(1);
*/

	int true_len = nSetBits+1;
	gOut.mask0.set_len(true_len);
	gOut.mask.set_len(true_len);
	gOut.n = true_len+1;
	gOut.check_first_unset_bit();
	return true;
}


bool extend(graph& gOut, const graph& g, int color, int kmax[2], int target, int new_bit, int search_mode)
{
	if(color && g.mask0.bit(new_bit))
	{
		return false;
	}
	if((!color) && g.mask.bit(new_bit))
	{
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

	if(g.parent_ext != 0)
		gOut.parent_ext = g.parent_ext;
	else
		gOut.parent_ext = 0;

	memcpy(gOut.clique_counts, g.clique_counts, sizeof(g.clique_counts));
	gOut.construct_shifted_masks(target);

	gOut.implied_mask[0].clear();
	gOut.implied_mask[1].clear();
	gOut.implied_mask_set = false;
	bool updated=false;
	int sizePrev = gOut.cliques[color].size();
	if(!(color ? gOut.mask.bit(new_bit) : gOut.mask0.bit(new_bit)))
	{
		if(!gOut.set_bit(color, new_bit, kmax, target, (search_mode==3)))
		{
			return false;
		}
		updated=true;
	}
	if(!updated)
		return true;

	if(gOut.clique_counts[0][kmax[0]]!=0 || gOut.clique_counts[1][kmax[1]]!=0)
	{
		return false;
	}

	if(gOut.first_unset_bit-1 < target && search_mode==2 && gOut.known_shifted_ref_masks)
	{
		bigint new_mask;
		gOut.forced_bit_search2(kmax, target, color, new_bit, new_mask, 0);
		if (gOut.implied_mask_set && !bit_equal(new_mask, gOut.implied_mask[1 - color]))
		{
			printf("extend->forced_bit_search2(%d,%d): %d forced bits\n", color, new_bit, popcnt(new_mask));
			printf("\n");
			new_mask.set_len(target);
			gOut.implied_mask[1 - color].set_len(target);
			gOut.implied_mask[color].set_len(target);
			new_mask.print();
			gOut.implied_mask[1 - color].print();
			gOut.implied_mask[color].print();
			printf("===\n");
		}
		gOut.implied_mask[1-color].clear();
		gOut.implied_mask_set = false;

		if (!new_mask.zero())
		{
			int c = 1 - color;
			while (true)
			{
				int prev_pos0 = gOut.cliques[c].size();
				bool ok = gOut.new_clique_search(kmax, c, target, new_mask);
				if (!ok)
					return false;
				gOut.forced_bit_search_oneside(kmax[c], target, c, new_mask, prev_pos0);
				gOut.implied_mask[1 - c].clear();
				gOut.implied_mask_set = false;
				if (new_mask.zero())
					break;
				c = 1 - c;
			}
		}
	}
	if(gOut.clique_counts[0][kmax[0]]!=0 || gOut.clique_counts[1][kmax[1]]!=0)
		return false;

	bigint error_mask = gOut.mask & gOut.mask0;
	if(popcnt(error_mask)!=0)
		return false;

	gOut.check_first_unset_bit();
	gOut.explored[0]=g.explored[0];
	gOut.explored[1]=g.explored[1];
	return true;
}

bool extend(graph& g, int color, int kmax[2], int target, int new_bit, bool full_reconst)
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
	if(!g.known_shifted_ref_masks || g.shifted_ref_mask_len!=target)
		g.construct_shifted_masks(target);
	if(!(color ? g.mask.bit(new_bit) : g.mask0.bit(new_bit)))
	{
		if(!g.set_bit(color, new_bit, kmax, target, !full_reconst))
		{
			return false;
		}
		updated=true;
	}

	if(g.clique_counts[0][kmax[0]]!=0 || g.clique_counts[1][kmax[1]]!=0)
	{
		return false;
	}

	bigint error_mask = g.mask & g.mask0;
	if(popcnt(error_mask)!=0)
		return false;

	g.check_first_unset_bit();
	return true;
}

bool max_extensibility(graph& g, int kmax[2], int target, bool debug)
{
	int i;

	g.construct_shifted_masks(target);

	int min_unset_bit = 35;
	if (kmax[0] < 5 || kmax[1] < 5)
		min_unset_bit = 10;
	if (g.first_unset_bit - 1 < min_unset_bit || g.first_unset_bit - 1 >= target - 20)
		return true;

	if(debug)
		printf("max_extensibility\n");

	for (i = 2; i < 8; i++)
		g.cliques[i].clear();

	int minpos0 = 0, minpos1 = 0;
	int pass_count = 0;
	int inside_pass_count = 0;
	bigint new_masks[2];

	while(true)
	{
		while (true)
		{
			new_masks[0].clear();
			new_masks[1].clear();
			int color_mask = 3;
		
			if (g_high_kmax_mode)
			{
				if (kmax[0] >= 10)
					color_mask &= ~1;
				if (kmax[1] >= 10)
					color_mask &= ~2;
			}
			g.forced_bit_search(kmax, target, new_masks, (inside_pass_count==0), color_mask, minpos0, minpos1);
			minpos0 = g.cliques[0].size();
			minpos1 = g.cliques[1].size();
			if (new_masks[0].zero() && new_masks[1].zero())
				break;
			if(!(new_masks[0] & g.mask0).zero())
			{
				printf("Error\n");
				new_masks[0].print();
				g.mask0.print();
				exit(-1);
			}
			if(!(new_masks[1] & g.mask).zero())
			{
				printf("Error\n");
				new_masks[0].print();
				g.mask0.print();
				exit(-1);
			}
			bool ok = g.new_clique_search(kmax, 0, target, new_masks[0]);
			if (!ok)
				return false;
			ok = g.new_clique_search(kmax, 1, target, new_masks[1]);
			if (!ok)
				return false;
			inside_pass_count++;
		}
		if (!(g.mask0 & g.mask).zero())
			return false;
	
		new_masks[0].clear();
		new_masks[1].clear();

		if (!link_pair_search(target, kmax, g, new_masks, debug, (pass_count==0)))
		{
			if (debug)
				printf("link_pair_search fail\n");
			return false;
		}

		if(new_masks[0].zero() && new_masks[1].zero())
			break;
		
		bool ok = g.new_clique_search(kmax, 0, target, new_masks[0]);
		if (!ok)
			return false; 
		ok = g.new_clique_search(kmax, 1, target, new_masks[1]);
		if (!ok)
			return false;

		pass_count++;
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

inline void construct_cand_ext(bigint& cand, bigint& y, const bigint* shifted_ref_masks)
{
	for(int p=0; p<GRAPH_MAX_LEN; p++)
	{
		while(y.n[p]!=0)
		{
			unsigned long pos;
			_BitScanForward64(&pos, y.n[p]);
			cand &= shifted_ref_masks[pos+p*64+1];
			if(cand.zero())
			{
				p=GRAPH_MAX_LEN;
				break;
			}
			y.n[p] &= ~(One<<pos);
		}
	}
}

void check_clique(bigint x, bigint m, int kmax)
{
#ifdef _DEBUG
	if(popcnt(x)!=kmax)
	{
		printf("Error: bad clique popcnt %d / %d\n", popcnt(x), kmax);
		exit(-1);
	}
	vector<int> v;
	while(!x.zero())
	{
		int p = x.leading_bit();
		v.push_back(p);
		x.unset(p);
	}
	if(!andnot(x,m).zero())
	{
		printf("Error: some bits of x missing in m\n");
		exit(-1);
	}
	for(int i=0; i<v.size(); i++)
		for(int j=i+1; j<v.size(); j++)
			if(!m.bit(v[j]-v[i]-1))
			{
				for(int k=0; k<v.size(); k++)
					printf("%d ", v[k]);
				printf("\n");
				m.printbits();
				printf("Bit %d missing\n", v[j]-v[i]-1);
				exit(-1);
			}
#endif
}

void construct_cand_ext_explicit(bigint& cand, bigint x, const bigint& ref_mask)
{
	vector<int> v;
	while(!x.zero())
	{
		int p = x.leading_bit();
		v.push_back(p);
		x.unset(p);
	}
	if(!andnot(x, ref_mask).zero())
	{
		printf("Error\n");
		exit(-1);
	}
	for(int i=0; i+1<GRAPH_MAX_LEN*64; i++)
	{
		if(!cand.bit(i+1))
			continue;
		if(ref_mask.bit(i))
		{
			printf("Error\n");
			exit(-1);
		}
		bool ok  = true;
		for(int j=0; j<v.size(); j++)
		{
			int p = abs(i-v[j])-1;
			if(!ref_mask.bit(p))
			{
				ok=false;
				break;
			}
		}
		if(!ok)
			cand.unset(i+1);
	}
}

int count_call=0;

bool construct_incomplete_cliques_simple(graph& g, int target, int kmax[2], bool use_parents, bool min_call=false)
{
	//printf("Call %d\n", count_call);
	//if(count_call==31)
	{
	//	printf("!");
	}
	//count_call++;
	int i, k;

	lite_vector* part_cliques = &g.cliques[2];
	int base = g.first_unset_bit-1;
	int max_dim=target-base;
	if(max_dim<0)
	{
		printf("ERROR: g.first_unset_bit %d, target %d\n", g.first_unset_bit, target);
		exit(-1);
	}
	//g.construct_shifted_masks(target);

	bigint setmask = g.mask | g.mask0;
#if 0
	bigint copy_shifted_ref_masks[GRAPH_MAX_LEN*64+5];
	memcpy(copy_shifted_ref_masks, g.shifted_ref_masks, sizeof(copy_shifted_ref_masks));
	
	g.construct_shifted_masks(target);
	int nfail=0;
	for(i=1; i<target; i++)
		if(!bit_equal(copy_shifted_ref_masks[i], g.shifted_ref_masks[i]))
		{
			printf("%d fail\n", i);
			g.shifted_ref_masks[i].set_len(255);
			copy_shifted_ref_masks[i].set_len(255);
			g.shifted_ref_masks[i].printbits();
			copy_shifted_ref_masks[i].printbits();

			nfail++;
		}
	if(nfail>0)
		exit(-1);
#endif
	
	if(popcnt(setmask)>target*9/10)
	{
		g.explored[0]=g.mask0;
		g.explored[1]=g.mask;
		for(i=0; i<6; i++)
			part_cliques[i].clear();
		return true;
	}
	
	char temp[GRAPH_MAX_LEN*64*GRAPH_MAX_LEN*64];
	char temp2[GRAPH_MAX_LEN*64];
	memset(temp, 0, target*target);
	memset(temp2, 0, sizeof(temp2));
	

	if(kmax[0]==1)
	{
		bigint b = g.mask0;
		while(!b.zero())
		{
			int p = b.trailing_bit();
			for(int j=0; j+p+1<target; j++)
			{
				if(g.mask0.bit(j))
					continue;
				if(!g.mask0.bit(j+p+1))			
				{
					//printf("^");
					//temp[(j+p+1)+j*target]|=1;
				}
				if(j>(p+1) && !g.mask0.bit(j-(p+1)))
				{
					//printf("V");
					//temp[j+(j-(p+1))*target]|=1;
				}
			}
			b.unset(p);
		}
	}


	// retains some of the previously discovered if there was a call with use_parents=true before
	//if(!use_parents)
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
				int tr=vc[j].trailing_bit();
				temp[tr+lb*target]|=1<<(i-4);
			}
		}
	}


	
	for(int d=0; d<g.parent_cliques[0].size(); d++)
	//int d=g.parent_cliques[0].size();
	for(i=4; i<6; i++)
	{
		int c=0;
		const lite_vector& vc = *g.parent_cliques[i][d];

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

	
	bool initial = g.explored[0].zero() && g.explored[1].zero();

	for(i=0; i<6; i++)
		part_cliques[i].clear();
	
	int stop_point = use_parents ? 0 : g.parent_cliques[0].size();
	if(use_parents && g.parent_cliques[0].size()>2)
		stop_point=g.parent_cliques[0].size()-2;
	buffer<int>& vmatch = g.vmatch;
	
	bigint* shifted_ref_masks = g.shifted_ref_masks;
	if(!g.known_shifted_ref_masks)
	{
		printf("ERROR: expect to have valid shifted_ref_masks\n");
		exit(-1);
	}

	bigint unexplored[2];
	unexplored[0] = andnot(g.mask0, g.explored[0]);
	unexplored[1] = andnot(g.mask, g.explored[1]);
	for(int d=g.parent_cliques[0].size(); d>=stop_point; d--)
	{
		int extensible=0;
		int local_cand=0;
		int cliques=0;
		int local_forced_miss=0;
		for(i=0; i<2; i++)
		{
			if(g.parent_cliques[i].size() != g.parent_cliques[0].size())
			{
				printf("?");
				exit(-1);
			}

			const lite_vector& vc = (d==g.parent_cliques[0].size()) ? g.cliques[i] : *g.parent_cliques[i][d];
			bool parent = (d!=g.parent_cliques[0].size());

			vmatch.clear();

			bigint ref_mask = (i==0) ? g.mask0 : g.mask;
			bigint dual_mask = (i==0) ? g.mask : g.mask0;

			bigint base_cand;
			int max_n = target;
			int min_n = base;
			base_cand.clear();
			base_cand.set_len(max_n);
			for(k=0; k<GRAPH_MAX_LEN; k++)
				base_cand.n[k]=-1;
			base_cand.unset_from(min_n);
			base_cand.neg();
			base_cand.unset_from(max_n);
			base_cand = andnot(base_cand, ref_mask<<1);
			base_cand = andnot(base_cand, dual_mask<<1);
			for(int j=0; j<vc.size(); j++)
			{
				if(vc[j].popcnt != kmax[i]-1)
					continue;
				if(parent && (vc[j]&unexplored[i]).zero())
					continue;
				cliques++;
				bigint y = vc[j];

				//bigint ref2 = ref_mask;
				//ref2.set(a-1);
				//ref2.set(b-1);
				check_clique(y, ref_mask, kmax[i]-1);
				bigint cand = base_cand;
				construct_cand_ext(cand, y, shifted_ref_masks);
#if 0			
				bigint cand2 = base_cand;
				construct_cand_ext_explicit(cand2, vc[j], ref_mask);
				if(!bit_equal(cand, cand2))
				{
					printf("construct_cand_ext failure\n");
					printf("base_cand:\n");
					base_cand.set_len(255);
					base_cand.printbits();
					printf("vc:\n");
					vc[j].printbits();
					printf("ref_mask:\n");
					ref_mask.set_len(255);
					ref_mask.printbits();

					cand.set_len(255);
					cand2.set_len(255);
					cand.printbits();
					cand2.printbits();
					exit(-1);
				}
#endif			
				int hits=popcnt(cand);
				if(hits<2)
					continue;
				vmatch.resize(hits);
				k=0;
				for(k=0; k<hits; k++)
				{
					unsigned long pos = cand.trailing_bit();
					vmatch[k]=pos;
					cand.unset(pos);
#ifdef _DEBUG
					bigint x  = vc[j];
					x.set(pos-1);
					bigint ref2 = ref_mask;
					ref2.set(pos-1);
					//ref2.set(b-1);
					check_clique(x, ref2, kmax[i]);
#endif
				}

				int local_ext=0;
				local_cand++;
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
							{
								g.explored[0]|=unexplored[0];
								g.explored[1]|=unexplored[1];
								return false;
							}
							if(bit0+bit1==1)
							{
	//							if(!(temp2[(bit0 ? b-1 : a-1)] & (1<<i)))
	//								if(!ref_mask.bit(bit0 ? b-1 : a-1))
	//									local_forced_miss++;
								temp2[(bit0 ? b-1 : a-1)]|=1<<i;
							}
							else
							{
	//							if(!(temp[(a-1)+(b-1)*target] & (1<<i)))
	//								local_ext=1;
#ifdef _DEBUG
								bigint x = vc[j];
								x.set(a-1);
								x.set(b-1);
								bigint ref2 = ref_mask;
								ref2.set(a-1);
								ref2.set(b-1);
								check_clique(x, ref2, kmax[i]+1);
#endif
								temp[(a-1)+(b-1)*target]|=1<<i;

							}
						}
					}
				}
				extensible+=local_ext;
			}
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
				{
					//printf("+");
					part_cliques[2].push_back(x);
				}
				if(temp[j+i*target]&2)
					part_cliques[3].push_back(x);
			}
		}
	}
	
	g.explored[0]|=unexplored[0];
	g.explored[1]|=unexplored[1];
	return true;
}

bool link_pair_search(int target, int kmax[2], graph& g, bigint new_masks[2], bool debug, bool use_parents)
{
//	return true;
	g.check_first_unset_bit();
	
	if(target < g.first_unset_bit+10)
	{
		g.explored[0]=g.mask0;
		g.explored[1]=g.mask;
		return true;
	}
	if (g_high_kmax_mode && (kmax[0] >= 10 || kmax[1] >= 10))
		return true;
//	lite_vector part_cliques[8];
	lite_vector* part_cliques = &g.cliques[2];
	int i;
	if(use_parents)
	{
//		for(i=0; i<6; i++)
//			part_cliques[i].clear();
	}
	if(!construct_incomplete_cliques_simple(g, target, kmax, use_parents))
		return false;
	int base = g.first_unset_bit-1;
	if(debug)
	{
		printf("g.first_unset_bit %d\n", base);
		g.mask0.printbits();
		g.mask.printbits();
	}
	int j, k, m;
	new_masks[0] = g.mask0;
	new_masks[1] = g.mask;
	new_masks[0].unset_to(base);
	new_masks[1].unset_to(base);

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
			if (pos < base)
				continue;
			if(new_masks[i].bit(pos))
				return false;
			new_masks[1 - i].set(pos);
			if(debug)
				printf("1-short: bit(%d) <= %d\n", pos, 1-i);
			if (pos < base)
			{
				printf("Error: pos<base @ line 1852 (%d,%d)\n", pos, base);
				exit(-1);
			}
		}
	}

	if(debug)
	{
		for(i=2; i<4; i++)
		{
			const lite_vector& v = part_cliques[i];
			for(j=0; j<v.size(); j++)
			{
				printf("Rule: (%d) %d %d\n", i-2, v[j].leading_bit(), v[j].trailing_bit());
			}
		}
	}
	bigint extensions[4][GRAPH_MAX_LEN*64];
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
	/*
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
	*/
	for (i = 0; i<max_dim; i++)
	{
		extensions[0][i] |= new_masks[0];
		extensions[2][i] |= new_masks[0];
		extensions[1][i] |= new_masks[1];
		extensions[3][i] |= new_masks[1];
	}

	while(true)
	{
		int prev_new_bits_size = popcnt(new_masks[0]) + popcnt(new_masks[1]);
		bigint prev_new_masks[2] = { new_masks[0], new_masks[1] };
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

					bigint rule_save = rule;
					while(!rule.zero())
					{
						int m= rule.trailing_bit();
						rule.unset(m);
						m-=base;
						if (m < 0)
						{
							printf("Assertion failure (line 1937)\n");
							exit(-1);
						}
						bigint rem = andnot(rule_save, extensions[c*3][m]);
						int p = popcnt(rem);
						if(p==0)						
							new_masks[1 - c].set(m + base);
						if(p==1) // forced bit
						{
							int pos = rem.trailing_bit();

							if (pos < base)
							{
								printf("Assertion failure (line 1950)\n");
								exit(-1);
							}

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
					new_masks[1 - i].set(m + base);
			}
		}
#endif
		//if(new_bits.size() != prev_new_bits_size)
		if(prev_new_bits_size != popcnt(new_masks[0]) + popcnt(new_masks[1]))
		{
			for (int c = 0; c < 2; c++)
			{
				bigint x = andnot(new_masks[c], prev_new_masks[c]);
				while (!x.zero())
				{
					int pos = x.leading_bit();
					extensions[0][j] |= extensions[c * 2 + 0][pos - base];
					extensions[1][j] |= extensions[c * 2 + 1][pos - base];
					extensions[2][j] |= extensions[c * 2 + 0][pos - base];
					extensions[3][j] |= extensions[c * 2 + 1][pos - base];
					x.unset(pos);
				}
			}
		}
		else
		{
			if(changes==0)
				break;
		}
		if(popcnt(new_masks[0], new_masks[1])>0)
		{
			if(debug)
				printf("Full\n");
			return false;
		}
	}

	double best_save=10, second_best_save=10;
	g.best_next_bit=-1;

	for(i=0; i<max_dim; i++)
	{
		if(!new_masks[0].bit(base+i) && !new_masks[1].bit(base+i))
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
				//g.second_best_next_bit = g.best_next_bit;
				best_save = save;
				g.best_next_bit = i+base;
			}
			else if(save<second_best_save)
			{
				second_best_save = save;
				//g.second_best_next_bit = i+base;
			}
		}
	}

	new_masks[0] = andnot(new_masks[0], g.mask0);
	new_masks[1] = andnot(new_masks[1], g.mask);
	return true;
}


bool graph::new_clique_search(int kmax[2], int color, int target, bigint new_mark)
{
	while (!new_mark.zero())
	{
		int new_bit = new_mark.leading_bit();
		if (!set_bit(color, new_bit, kmax, target))
			return false;
		new_mark.unset(new_bit);
	}
	return true;
}

void graph::forced_bit_search2(int kmax[2], int target, int color, int new_bit, bigint& new_mask, int sizePrev)
{
	new_mask.clear();

	int i=color;
	const lite_vector& vc = cliques[i];
	bigint& ref_mask = i ? mask : mask0;
	bigint& dual_mask = i ? mask0 : mask;

	int max_n = target;
	int min_n = first_unset_bit;

	bigint base_cand;
	base_cand.clear();
	base_cand.set_len(max_n);
	for(int k=0; k<GRAPH_MAX_LEN; k++)
		base_cand.n[k]=-1;
	base_cand.unset_from(min_n);
	base_cand.neg();
	base_cand.unset_from(max_n);
	int j;
	base_cand= andnot(base_cand, dual_mask<<1);

#if 1
	for(j=sizePrev; j<vc.size(); j++)
	{
		if(vc[j].popcnt != kmax[color])
			continue;
		if(!(vc[j].bit(new_bit)))
			continue;

		bigint cand = base_cand;
		bigint y = vc[j];
		construct_cand_ext(cand, y, shifted_ref_masks);
		base_cand = andnot(base_cand, cand);
		cand >>= 1;
		new_mask = new_mask | cand;
	}
#else
	buffer<int> new_marks[2];
	bool debug=false;
	bigint new_mask2;
	new_mask2.clear();
	for(j=sizePrev; j<vc.size(); j++)
	{
		if(vc[j].popcnt != kmax[color])
			continue;
		if(!(vc[j].bit(new_bit)))
			continue;

		bigint cand = base_cand;
		bigint y = vc[j];
		construct_cand_ext(cand, y, shifted_ref_masks);

		while(!cand.zero())
		{
			int n = cand.leading_bit();
			cand.unset(n);
			base_cand.unset(n);
			new_marks[1-i].push_back(n);
			//cand >>= 1;
			//new_mask2 = new_mask2 | cand;
		}
	}
	new_mask.clear();
	for(int i=0; i<new_marks[1-color].size(); i++)
	{
		new_mask.set(new_marks[1-color][i]-1);
	}
#endif
}

void graph::forced_bit_search(int kmax[2], int target, bigint new_masks[2], bool use_parent_cliques, int color_mask, int minpos0, int minpos1)
{
	int j;
	new_masks[0].clear();
	new_masks[1].clear();

	int max_n = target;
	int min_n = first_unset_bit;
	
	if(kmax[0]==1 && (color_mask & 1))
	{
		bigint b;
		b.clear();
		bigint c = mask0;
		while(!c.zero())
		{
			int p = c.trailing_bit();
			bigint d = mask0 << (p+1);
			b = b | d;
			d = mask0 >> (p+1);
			b = b | d;
			c.unset(p);
		}
		b = andnot(b, mask);
		b.unset_from(target-1);
		while(!b.zero())
		{
			int p = b.trailing_bit();
			new_masks[1].set(p);
			b.unset(p);
		}
		//color_mask &= ~1;
	}
	
	for(int i=0; i<2; i++)
	{
		if(!(color_mask & (1<<i)))
			continue;
		const bigint& ref_mask = i ? mask : mask0;
		const bigint& dual_mask = i ? mask0 : mask;

		int bad_bits=0;
		bigint bad_mask;
		bad_mask.clear();
		for(int j=0; j<target; j++)
		{
			if(!ref_mask.bit(j))
				continue;
			bigint cand;
			cand.clear();
			cand.set_len(max_n);
			for(int k=0; k<GRAPH_MAX_LEN; k++)
				cand.n[k]=-1;
			cand.unset_from(min_n);
			cand.neg();
			cand.unset_from(max_n);
//				cand = andnot(cand, dual_mask<<1);
			cand = andnot(cand, dual_mask<<1);
			cand = andnot(cand, new_masks[1-i]<<1);
			cand &= shifted_ref_masks[j+1];
			//cand &= shifted_ref_masks[k+1];
			if(cand.zero())
				//bad_bits++;
				bad_mask.set(j);
		}

		bigint base_cand;
		base_cand.clear();
		base_cand.set_len(max_n);
		for(int k=0; k<GRAPH_MAX_LEN; k++)
			base_cand.n[k]=-1;
		base_cand.unset_from(min_n);
		base_cand.neg();
		base_cand.unset_from(max_n);
		base_cand= andnot(base_cand, dual_mask<<1);
		base_cand= andnot(base_cand, new_masks[1-i]<<1);

		for(int cl=use_parent_cliques ? 0 : parent_cliques[i].size(); cl<parent_cliques[i].size()+1; cl++)
		{			
			const lite_vector& vc = (cl<parent_cliques[i].size()) ? *parent_cliques[i][cl] : cliques[i];
			for(int j=(i==0 ? minpos0 : minpos1); j<vc.size(); j++)
			{
				if(vc[j].popcnt != kmax[i])
					continue;
				if(!(vc[j]&bad_mask).zero())
					continue;
				bigint cand = base_cand;
				bigint y = vc[j];
				construct_cand_ext(cand, y, shifted_ref_masks);
				base_cand = andnot(base_cand, cand);
				cand >>= 1;
				new_masks[1 - i] |= cand;
			}
		}
	}
}




void graph::forced_bit_search_oneside(int kmax, int target, int color, bigint& new_mask, int minpos)
{
	int j;
	//bigint new_bits[2];
	//new_bits[0].clear();
	//new_bits[1].clear();
	new_mask.clear();

	int max_n = target;
	int min_n = first_unset_bit;

	const bigint& ref_mask = color ? mask : mask0;
	const bigint& dual_mask = color ? mask0 : mask;
	if (kmax == 1)
	{
		bigint b;
		b.clear();
		bigint c = ref_mask;
		while (!c.zero())
		{
			int p = c.trailing_bit();
			bigint d = ref_mask << (p + 1);
			b = b | d;
			d = ref_mask >> (p + 1);
			b = b | d;
			c.unset(p);
		}
		b = andnot(b, dual_mask);
		b.unset_from(target - 1);
		while (!b.zero())
		{
			int p = b.trailing_bit();
			new_mask.set(p);
			b.unset(p);
		}
		//color_mask &= ~1;
	}

	int bad_bits = 0;
	bigint bad_mask;
	bad_mask.clear();
	for (int j = 0; j<target; j++)
	{
		if (!ref_mask.bit(j))
			continue;
		bigint cand;
		cand.clear();
		cand.set_len(max_n);
		for (int k = 0; k<GRAPH_MAX_LEN; k++)
			cand.n[k] = -1;
		cand.unset_from(min_n);
		cand.neg();
		cand.unset_from(max_n);
		//				cand = andnot(cand, dual_mask<<1);
		cand = andnot(cand, dual_mask << 1);
		//cand = andnot(cand, new_mask << 1);
		cand &= shifted_ref_masks[j + 1];
		//cand &= shifted_ref_masks[k+1];
		if (cand.zero())
			//bad_bits++;
			bad_mask.set(j);
	}

	bigint base_cand;
	base_cand.clear();
	base_cand.set_len(max_n);
	for (int k = 0; k<GRAPH_MAX_LEN; k++)
		base_cand.n[k] = -1;
	base_cand.unset_from(min_n);
	base_cand.neg();
	base_cand.unset_from(max_n);
	base_cand = andnot(base_cand, dual_mask << 1);
	//base_cand = andnot(base_cand, new_mask << 1);

	const lite_vector& vc = cliques[color];
	for (int j = minpos; j<vc.size(); j++)
	{
		if (vc[j].popcnt != kmax)
			continue;
		if (!(vc[j] & bad_mask).zero())
			continue;
		bigint cand = base_cand;
		bigint y = vc[j];
		construct_cand_ext(cand, y, shifted_ref_masks);
		base_cand = andnot(base_cand, cand);
		cand >>= 1;
		new_mask = new_mask | cand;
	}
}
