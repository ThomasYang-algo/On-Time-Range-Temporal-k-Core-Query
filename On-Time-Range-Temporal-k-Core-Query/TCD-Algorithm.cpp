#include<iostream>
#include<fstream>
#include<cstring>
#include<string>
#include<vector>
#include<cstdio>
#include<unordered_map>
#include<unordered_set>
#include<set>
#include<algorithm>
#include<chrono>
#include<array>
using namespace std;
using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::chrono::nanoseconds;
using std::chrono::system_clock;
using std::chrono::high_resolution_clock;

/*** TEL Ultra ***/
const int VMAX = 4000000;  //Default:4000000
const int EMAX = 130000000;//Default:130000000
const int TMAX = 1665000000;//Default:1665000000

int vern = 0;
int arcn = 0;
int tmax = 0;

struct arc {
	int src, dst, t;
}arcs[EMAX];

int *hs = 0;
int *hd = 0;
int *ht = 0;
int *spre = 0, *snxt = 0;
int *dpre = 0, *dnxt = 0;
int *tpre = 0, *tnxt = 0;
int *telarc = 0;
int idx = 0;

int head = -1, tail = -1;
int *tsarc = 0;
int *tsnxt = 0, *tspre = 0;
int idt = 0;

#define _init_arr(name,size) name = new int[size]

void initmem()
{
	_init_arr(hs, VMAX);
	_init_arr(hd, VMAX);
	_init_arr(ht, TMAX);
	_init_arr(spre, EMAX), _init_arr(snxt, EMAX);
	_init_arr(dpre, EMAX), _init_arr(dnxt, EMAX);
	_init_arr(tpre, EMAX), _init_arr(tnxt, EMAX);
	_init_arr(telarc, EMAX);

	_init_arr(tsarc, EMAX);
	_init_arr(tsnxt, EMAX);
	_init_arr(tspre, EMAX);
}

void addt(int t)
{
	if (head == -1) head = idt;
	tsarc[idt] = t;
	tspre[idt] = tail, tsnxt[idt] = -1; 
	if (tail >= 0) tsnxt[tail] = idt; tail = idt;
	idt++;
}

void delt(int t)
{
	int i = (lower_bound(tsarc, tsarc + idt, t) - tsarc);
	if (head == i) head = tsnxt[i];
	else if (tsnxt[i] == -1) tsnxt[tspre[i]] = -1, tail = tspre[i];
	else {
		tsnxt[tspre[i]] = tsnxt[i];
		tspre[tsnxt[i]] = tspre[i];
	}
}

void addarc(int id, int src, int dst, int t)
{
	telarc[idx] = id;
	snxt[idx] = hs[src]; if (hs[src] >= 0) spre[hs[src]] = idx; hs[src] = idx;
	dnxt[idx] = hd[dst]; if (hd[dst] >= 0) dpre[hd[dst]] = idx; hd[dst] = idx;
	tnxt[idx] = ht[t]; if (ht[t] >= 0) tpre[ht[t]] = idx; ht[t] = idx;
	idx++;
}

void delarc_l(int *head, int *next, int* pre, int i, int idx)
{
	if (head[i] == idx) head[i] = next[idx];
	else if (next[idx] == -1) next[pre[idx]] = -1;
	else {
		next[pre[idx]] = next[idx];
		pre[next[idx]] = pre[idx];
	}
}

/*void delarc(int id, int src, int dst, int t)
{
for(int i = ht[t]; ~i; i = tnxt[i])
if(telarc[i] == id)
{
delarc_l(hs, snxt, spre, src, i);
delarc_l(hd, dnxt, dpre, dst, i);
delarc_l(ht, tnxt, tpre, t, i);
}
if(ht[t] == -1) delt(t);
}*/


/*** DataSet Part ***/
const int QMAX = 45000;
struct Q {
	int l, r, k;
}q[QMAX];
int qcnt = 0;

void loadgraph(const char* name)
{
	ifstream fin(name, ios::in);
#ifdef _DEBUG
	if (fin.is_open() == false) { printf("open graph %s fail\n", name); exit(1); }
#endif
	vector<int> v;
	int tmin = 0x7fffffff;

	string l;
	while (getline(fin, l))
	{
		int uvt[3] = { 0 };
		int p = -1, np = -1;
		for (int i = 0; i < 3; ++i)
		{
			p = np + 1, np = l.find(' ', np + 1);
			if (np == -1) np = l.size();
			uvt[i] = stoi(l.substr(p, np - p));
		}
		v.push_back(uvt[0]), v.push_back(uvt[1]);
		tmin = min(tmin, uvt[2]);
		arcs[arcn++] = { uvt[0], uvt[1], uvt[2] };
	}

	sort(v.begin(), v.end());
	v.erase(unique(v.begin(), v.end()), v.end());
	vern = v.size();

	auto get = [&](int k) {
		return (lower_bound(v.begin(), v.end(), k) - v.begin()) + 1;
	};

	for (int i = 0; i < arcn; ++i)
	{
		arcs[i].t -= tmin;
		arcs[i].src = get(arcs[i].src), arcs[i].dst = get(arcs[i].dst);
		if (arcs[i].src > arcs[i].dst) swap(arcs[i].src, arcs[i].dst);
		tmax = max(tmax, arcs[i].t);
	}

}

void loadtest(const char* name)
{
	ifstream fin(name, ios::in);
#ifdef _DEBUG
	if (fin.is_open() == false) { printf("open test %s fail\n", name); exit(1); }
#endif

	string l;
	while (getline(fin, l))
	{
		int lrk[3] = { 0 };
		int p = -1, np = -1;
		for (int i = 0; i < 3; ++i)
		{
			p = np + 1, np = l.find('\t', np + 1);
			if (np == -1) np = l.size();
			lrk[i] = stoi(l.substr(p, np - p));
		}
		q[qcnt++] = { lrk[0], lrk[1], lrk[2] };
	}

}


/*** TCD Operation Part ***/
template<typename T>
void hash_combine(size_t& seed, T const& v)
{
	seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

struct pair_hash
{
	template<typename T1, typename T2>
	size_t operator()(std::pair<T1, T2>const& p) const
	{
		size_t seed = 0;
		hash_combine(seed, p.first);
		hash_combine(seed, p.second);
		return seed;
	}
};

unordered_map<int, int> Mv;
unordered_map<pair<int, int>, int, pair_hash> Mc;
set<pair<int, int>> Hv;

bool cUpd(int src, int dst)
{
#ifdef _DEBUF
	if (Mc.count({ src, dst }) == 0) { printf("cUpd:empty update\n"); exit(1); }
#endif
	bool ret = false;
	Mc[{src, dst}] --;
	if (Mc[{src, dst}] == 0) { Mc.erase({ src, dst }); ret = true; }
	return ret;
}

void vUpd(int v)
{
#ifdef _DEBUG 
	if (Mv.count(v) == 0) { printf("vUpd:empty update\n"); exit(1); }
#endif
	int d = Mv[v];
	Hv.erase({ d, v });
	Hv.insert({ d - 1, v });
	Mv[v] --;
}

void trunc(int l, int r)
{
	int hh = head, tt = tail;
	while (hh >= 0 && tsarc[hh] < l)
	{
		int t = tsarc[hh];
		for (int i = ht[t]; ~i; i = tnxt[i])
		{
			int id = telarc[i];
			delarc_l(hs, snxt, spre, arcs[id].src, i);
			delarc_l(hd, dnxt, dpre, arcs[id].dst, i);
			delarc_l(ht, tnxt, tpre, arcs[id].t, i);
			if (ht[arcs[id].t] == -1) delt(arcs[id].t);
			if (cUpd(arcs[id].src, arcs[id].dst)) {
				vUpd(arcs[id].src);
				vUpd(arcs[id].dst);
			}
		}
		hh = tsnxt[hh];
	}
	while (tt >= 0 && tsarc[tt] > r)
	{
		int t = tsarc[tt];
		for (int i = ht[t]; ~i; i = tnxt[i])
		{
			int id = telarc[i];
			delarc_l(hs, snxt, spre, arcs[id].src, i);
			delarc_l(hd, dnxt, dpre, arcs[id].dst, i);
			delarc_l(ht, tnxt, tpre, arcs[id].t, i);
			if (ht[arcs[id].t] == -1) delt(arcs[id].t);
			if (cUpd(arcs[id].src, arcs[id].dst)) {
				vUpd(arcs[id].src);
				vUpd(arcs[id].dst);
			}
		}
		tt = tspre[tt];
	}
	head = hh, tail = tt;
}

void decomp(int k)
{
	while (Hv.size() && (Hv.begin()->first<k))
	{
		auto nv = *(Hv.begin());
		Hv.erase(Hv.begin());
		int n = nv.first, v = nv.second;
		Mv.erase(v);
		unordered_set<int> nbr;
		for (int i = hs[v]; ~i; i = snxt[i])
		{
			int id = telarc[i];
			delarc_l(hs, snxt, spre, arcs[id].src, i);
			delarc_l(hd, dnxt, dpre, arcs[id].dst, i);
			delarc_l(ht, tnxt, tpre, arcs[id].t, i);
			if (ht[arcs[id].t] == -1) delt(arcs[id].t);
			cUpd(arcs[id].src, arcs[id].dst);
			int u = arcs[id].src == v ? arcs[id].dst : arcs[id].src;
			nbr.insert(u);
		}
		for (int i = hd[v]; ~i; i = dnxt[i])
		{
			int id = telarc[i];
			delarc_l(hs, snxt, spre, arcs[id].src, i);
			delarc_l(hd, dnxt, dpre, arcs[id].dst, i);
			delarc_l(ht, tnxt, tpre, arcs[id].t, i);
			if (ht[arcs[id].t] == -1) delt(arcs[id].t);
			cUpd(arcs[id].src, arcs[id].dst);
			int u = arcs[id].src == v ? arcs[id].dst : arcs[id].src;
			nbr.insert(u);
		}
		for (auto u : nbr) vUpd(u);
	}
}

void tcdop(int l, int r, int k)
{
	trunc(l, r);
	decomp(k);
}


/*** TCD Algorithm Part ***/

int *_hs = 0;
int *_hd = 0;
int *_ht = 0;
int *_spre = 0, *_snxt = 0;
int *_dpre = 0, *_dnxt = 0;
int *_tpre = 0, *_tnxt = 0;
int *_telarc = 0;
int _idx = 0;

int _head = -1, _tail = -1;
int *_tsarc = 0;
int *_tsnxt = 0, *_tspre = 0;
int _idt = 0;

unordered_map<int, int> _Mv;
unordered_map<pair<int, int>, int, pair_hash> _Mc;
set<pair<int, int>> _Hv;

int vbytes = 0; //sizeof hd and hs
int cbytes = 0; //sizeof telarc
int tsbytes = 0;//sizeof tsarc
unsigned long long htbytes = 0ull;//sizeof ht

void _initmem()
{
	_init_arr(_hs, VMAX);
	_init_arr(_hd, VMAX);
	_init_arr(_ht, TMAX);
	_init_arr(_spre, EMAX), _init_arr(_snxt, EMAX);
	_init_arr(_dpre, EMAX), _init_arr(_dnxt, EMAX);
	_init_arr(_tpre, EMAX), _init_arr(_tnxt, EMAX);
	_init_arr(_telarc, EMAX);

	_init_arr(_tsarc, EMAX);
	_init_arr(_tsnxt, EMAX);
	_init_arr(_tspre, EMAX);
}

void bkp()
{
	memcpy(_hs, hs, vbytes);
	memcpy(_hd, hd, vbytes);
	memcpy(_ht, ht, htbytes);
	memcpy(_spre, spre, cbytes); memcpy(_snxt, snxt, cbytes);
	memcpy(_dpre, dpre, cbytes); memcpy(_dnxt, dnxt, cbytes);
	memcpy(_tpre, tpre, cbytes); memcpy(_tnxt, tnxt, cbytes);
	memcpy(_telarc, telarc, cbytes);

	_head = head, _tail = tail;
	memcpy(_tsarc, tsarc, tsbytes);
	memcpy(_tsnxt, tsnxt, tsbytes);
	memcpy(_tspre, tspre, tsbytes);

	_Mv = Mv;
	_Mc = Mc;
	_Hv = Hv;
}

void rst()
{
	memcpy(hs, _hs, vbytes);
	memcpy(hd, _hd, vbytes);
	memcpy(ht, _ht, htbytes);
	memcpy(spre, _spre, cbytes); memcpy(snxt, _snxt, cbytes);
	memcpy(dpre, _dpre, cbytes); memcpy(dnxt, _dnxt, cbytes);
	memcpy(tpre, _tpre, cbytes); memcpy(tnxt, _tnxt, cbytes);
	memcpy(telarc, _telarc, cbytes);

	head = _head, tail = _tail;
	memcpy(tsarc, _tsarc, tsbytes);
	memcpy(tsnxt, _tsnxt, tsbytes);
	memcpy(tspre, _tspre, tsbytes);

	Mv = _Mv;
	Mc = _Mc;
	Hv = _Hv;
}

const int infosum =  12;
enum INFO
{
	CELL,
	SQZ,
	RLC,
	TAG,
	PoR,
	PoU,
	PoL,
	CPoR,
	CPoU,
	CPoL,
	TOTALSIZE,
	SCCNUM,
};
long long cntinfo[infosum] = {0ll};

#ifdef CALC_PRUNE
const int MAX_SPAN = 11000;
auto prune_flag = new array<array<int, MAX_SPAN>, MAX_SPAN>; //The map for calculating prune percentage
#endif

void proc(int ts, int te);//tailored processing function

void tcd(int ql, int qr, int qk)
{
	int ts = ql;
	while (ts <= qr)
	{
		int te = qr;
		tcdop(ts, te, qk);
		if (head == -1 || tail == -1 || tsarc[head] > tsarc[tail]) break;
		bkp();
		while (te >= ts)
		{
			proc(ts, te);
			te--;
			tcdop(ts, te, qk);
			if (head == -1 || tail == -1 || tsarc[head] > tsarc[tail]) break;
		}
		ts++;
		rst();
	}
}

int Ts = 0, Te = 0; //Store query interval for indentation of when calculating prune percentage
void otcd(int ql, int qr, int qk)
{
	Ts = ql, Te = qr;
	unordered_map<int, int> endr; //end_pos for RoU
	int ts = 0, te = 0;      //current cell
	int _ts = ql, _te = qr;	 //first cell in next row
	auto rlc = [&](int __te) { 
		if (__te < te)
		{
			cntinfo[INFO::RLC]++;

#ifdef CALC_PRUNE
			cntinfo[INFO::PoR]++;
			for (int c = te - 1; c >= __te; --c)//Accumulate PoR Cell
				if ((*prune_flag)[ts - Ts][c - Ts] == 0)
			(*prune_flag)[ts - Ts][c - Ts] = 1;
#endif

			te = __te;//PoR
		}
	}; //PoR
	auto sqz = [&](int __ts, int __te) {
		if (ts < __ts || te > __te) { //Overline Trigger
			cntinfo[INFO::SQZ]++;

#ifdef CALC_PRUNE
			if (ts < __ts) {//Accumulate PoU Cell
				cntinfo[INFO::PoU]++;
				for (int r = ts + 1; r <= __ts; ++r)
					for (int c = te; c >= r; --c)
						if ((*prune_flag)[r - Ts][c - Ts] == 0)
					(*prune_flag)[r - Ts][c - Ts] = 2;
			}
			if (__te < te) {//Accumulate PoL Cell
				cntinfo[INFO::PoL]++;
				for (int r = __ts + 1; r <= __te; ++r)
					for (int c = te; c > __te; --c)
						if ((*prune_flag)[r - Ts][c - Ts] == 0)
					(*prune_flag)[r - Ts][c - Ts] = 3;
				cntinfo[INFO::PoR]++;
				for (int c = te - 1; c >= __te; --c)
					if ((*prune_flag)[ts - Ts][c - Ts] == 0)
				(*prune_flag)[ts - Ts][c - Ts] = 1;
			}
#endif

			_ts = __ts;  //PoU
			_te = __te;  //PoL
		}
	};
	auto tag = [&](int __ts)
	{
		if (ts < __ts) //Overline Triggered
		{
			cntinfo[INFO::TAG]++;
			for (int r = ts + 1; r <= __ts; ++r)
			{
				if (endr.count(r) == 0)
					endr[r] = -1;
				endr[r] = max(endr[r], te + 1); //Tag for PoU
			}

#ifdef CALC_PRUNE
			cntinfo[INFO::PoU]++;
			for (int r = ts + 1; r <= __ts; ++r)//Accumulate PoU Cell
				for (int c = te; c >= r; --c)
					if ((*prune_flag)[r - Ts][c - Ts] == 0)
				(*prune_flag)[r - Ts][c - Ts] = 2;
#endif

		}
	};
	auto pou = [&](int ts, int te) { return endr.count(ts) && endr[ts] >= te; };//PoU

	while (_ts <= _te)
	{
		ts = _ts, te = _te;
		tcdop(ts, te, qk);
		if (head == -1 || tail == -1 || tsarc[head] > tsarc[tail]) break;
		bkp();
		sqz(tsarc[head], tsarc[tail]);
		while (ts <= te)
		{
			proc(ts, te);
			rlc(tsarc[tail]);
			tag(tsarc[head]);
			if (pou(ts, te)) break;
			te--;
			tcdop(ts, te, qk);
			if (head == -1 || tail == -1 || tsarc[head] > tsarc[tail]) break;
		}
		_ts++;
		rst();
	}
}

void buildtel(int l, int r)
{
	memset(hs, -1, VMAX * sizeof(int));
	memset(hd, -1, VMAX * sizeof(int));
	memset(ht, -1, TMAX * sizeof(int));
	head = -1, tail = -1;
	idx = 0, idt = 0;
	vector<int> ts;
	for (int i = 0; i < arcn; ++i)
		if (arcs[i].t >= l && arcs[i].t <= r)
		{
			addarc(i, arcs[i].src, arcs[i].dst, arcs[i].t);
			ts.push_back(arcs[i].t);
		}

	sort(ts.begin(), ts.end());
	ts.erase(unique(ts.begin(), ts.end()), ts.end());
	for (auto t : ts) addt(t);

	vbytes = (vern + 1) * sizeof(int);
	cbytes = idx * sizeof(int);
	tsbytes = idt * sizeof(int);
	htbytes = (tmax + 1ull) * sizeof(int);
}

bool cAdd(int src, int dst)
{
	if (Mc.count({ src, dst }) == 0)
	{
		Mc[{src, dst}] = 1;
		return true;
	}
	Mc[{src, dst}] ++;
	return false;
}

void vAdd(int v)
{
	if (Mv.count(v) == 0) Mv[v] = 0;
	Hv.erase({ Mv[v], v });
	Mv[v] ++;
	Hv.insert({ Mv[v], v });
}

void initMH(int l, int r)
{
	Mv.clear();
	Mc.clear();
	Hv.clear();
	for (int i = 0; i < arcn; ++i)
	{
		int src = arcs[i].src;
		int dst = arcs[i].dst;
		int t = arcs[i].t;
		if (t < l || t > r) continue;
		if (cAdd(src, dst))
		{
			vAdd(src);
			vAdd(dst);
		}
	}
}


/*** Extension Part ***/
const int MAXTTI = 20010;
struct TTI {
	int ts, td, l;
}tti[MAXTTI];
int tticnt = 0;

bool lcmp(const TTI& x, const TTI& y)
{
	return x.l < y.l;
}	

bool tcmp(const TTI& x, const TTI& y)
{
	if (x.ts != y.ts) return x.ts < y.ts;
	return x.td < y.td;
}

const int HuangXin = 51471;
const int JianPei = 16184;
const int LinXuemin = 16185;
const int JefferyYuXu = 493399;
int jpcnt = 0;
int lxmcnt = 0;

//Disjoint Set for Retrieving Connecting Component Containing the Given User
int p[VMAX];
int find(int x)
{
	if (p[x] != x) p[x] = find(p[x]);
	return p[x];
}

void pego(int uid)
{
	set<int> vego;
	set<pair<int,int>> cego;
	static int ego_id = 1;
	if (hs[uid] != -1 || hd[uid] != -1)
	{
		vego.insert(uid);
		for (int i = hs[uid]; ~i; i = snxt[i])
		{
			int eid = telarc[i];
			vego.insert(arcs[eid].dst);
		}
		for (int i = hd[uid]; ~i; i = dnxt[i])
		{
			int eid = telarc[i];
			vego.insert(arcs[eid].src);
		}
		printf("Ego-%d of %d, size:%d\n", ego_id, uid, (int)vego.size());
		ego_id++;
		if (vego.size() > 25) return;
		for (auto v : vego) printf("%d\n", v);
		int tl = 2000000000, tr = -1;
		for (int i = head; ~i; i = tsnxt[i])
		{
			int t = tsarc[i];
			for (int j = ht[t]; ~j; j = tnxt[j])
			{
				int eid = telarc[j];
				int src = arcs[eid].src;
				int dst = arcs[eid].dst;
				int tm = arcs[eid].t;
				if (vego.count(src) && vego.count(dst))
				{
					tl = min(tl, tm);
					tr = max(tr, tm);
					cego.insert({ src, dst });
				}
			}
		}
		printf("edge timestamp range:[%d,%d]\n", tl+61, tr+61);
		for (auto e : cego) printf("%d,%d\n", e.first, e.second);

		printf("\n");
	}
}

void pcc(int uid)
{
	static int core_id = 1;
	if (hs[uid] != -1 || hd[uid] != -1)
	{
		for (int i = 1; i <= vern; ++i) p[i] = i;
		for (int i = head; ~i; i = tsnxt[i])//Traverse all edges of Core to initialize Disjoint Set
		{
			int ts = tsarc[i];
			for (int j = ht[ts]; ~j; j = tnxt[j])
			{
				int eid = telarc[j];
				int src = arcs[eid].src;
				int dst = arcs[eid].dst;
				if (find(src) != find(dst)) p[find(src)] = p[find(dst)];
			}
		}

		set<int> core_vset;
		for (int v = 1; v <= vern; ++v)
			if (find(v) == find(uid))
				core_vset.insert(v);
		printf("Size of Core-%d Connect Component Containing Jian Pei:%d\n", core_id, (int)core_vset.size());
		core_id++;
	}
}

int GetNumberOfVertex()//Get the Number of Vertex for Current Temporal K-Core
{
	unordered_set<int> v;
	for (int i = 1; i <= vern; ++i)
		if (hs[i] != -1)
			v.insert(i);
	for (int i = 1; i <= vern; ++i)
		if (hd[i] != -1)
			v.insert(i);
	return v.size();
}

void PrintCoreVertex()
{
	unordered_set<int> v;
	for (int i = 1; i <= vern; ++i)
		if (hs[i] != -1)
			v.insert(i);
	for (int i = 1; i <= vern; ++i)
		if (hd[i] != -1)
			v.insert(i);
	for (auto ver : v) printf("%d,%d\n", ver, ver);
}

void PrintCoreEdge()//Print all edges of Current Temporal K-Core
{
	for (int i = head; ~i; i = tsnxt[i])
	{
		int t = tsarc[i];
		for (int j = ht[t]; ~j; j = tnxt[j])
		{
			int eid = telarc[j];
			printf("%d,%d\n", arcs[eid].src, arcs[eid].dst);
		}
	}
}

long long GetNumberOfScc()//Compute number of Scc in Current Core
{
	for (auto vn : Mv) p[vn.first] = vn.first;
	for (int i = head; ~i; i = tsnxt[i])
	{
		int t = tsarc[i];
		for (int j = ht[t]; ~j; j = tnxt[j])
		{
			int eid = telarc[j];
			int src = arcs[eid].src;
			int dst = arcs[eid].dst;
			if (find(src) != find(dst))
			{
				p[find(src)] = find(dst);
			}
		}
	}
	long long scc_cnt = 0ll;
	for (auto vn : Mv)
		if (p[vn.first] == vn.first)
			scc_cnt++;
	return scc_cnt;
}

void proc(int ts, int te)
{
	cntinfo[INFO::CELL]++;
}

#ifdef CALC_PRUNE
void calc_prune()
{
	for(int i = 0; i < MAX_SPAN; ++ i)
		for (int j = 0; j < MAX_SPAN; ++j)
		{
			switch ((*prune_flag)[i][j])
			{
			case 1:
				cntinfo[INFO::CPoR] ++;
				break;
			case 2:
				cntinfo[INFO::CPoU] ++;
				break;
			case 3:
				cntinfo[INFO::CPoL] ++;
				break;
			}
		}
}
#endif

void pinfo(int qid, const char* graphname, int tcdcell, long long clapse, long long oclapse)
{
	printf("Graph Name:%s\n", graphname);
	printf("Platform:Omen Laptop15\n");
	printf("Query:%d\t%d\t%d\n", q[qid].l, q[qid].r, q[qid].k);
	printf("TCD Time Clapse(nanoseconds):%lld\n", clapse);
	printf("TCD Accessed Cell:%d\n", tcdcell);
	printf("OTCD Time Clapse(nanoseconds):%lld\n", oclapse);
	printf("OTCD Accessed Cell:%lld\n", cntinfo[INFO::CELL]);

/*
	printf("SQZ:%lld\n", cntinfo[INFO::SQZ]);
	printf("RLC:%lld\n", cntinfo[INFO::RLC]);
	printf("TAG:%lld\n", cntinfo[INFO::TAG]);
*/

#ifdef CALC_PRUNE
	printf("PoR:%lld\n", cntinfo[INFO::PoR]);
	printf("PoU:%lld\n", cntinfo[INFO::PoU]);
	printf("PoL:%lld\n", cntinfo[INFO::PoL]);
	printf("CPoR:%lld\n", cntinfo[INFO::CPoR]);
	printf("CPoU:%lld\n", cntinfo[INFO::CPoU]);
	printf("CPoL:%lld\n", cntinfo[INFO::CPoL]);
#endif

/*
	printf("TOTALSIZE:%lld\n", cntinfo[INFO::TOTALSIZE]);
	printf("SCCNUM:%lld\n", cntinfo[INFO::SCCNUM]);
*/

}

void proctti()
{
	sort(tti, tti + tticnt, lcmp);
	printf("min tti:%d\n", tti[0].l);
	printf("max tti:%d\n",tti[tticnt - 1].l);
	printf("tti medium:%d\n",tti[tticnt / 2].l);
	ofstream ttiset_l;
	ttiset_l.open("./youtube-tti-lorder.txt", ios::out);
	for (auto t : tti) ttiset_l << t.ts << '\t' << t.td << '\t' << t.l << '\n';

	sort(tti, tti + tticnt, tcmp);
	ofstream ttiset_t;
	ttiset_t.open("./youtube-tti-torder.txt", ios::out);
	for (auto t : tti) ttiset_t << t.ts << '\t' << t.td << '\t' << t.l << '\n';
}

/*** Main Function ***/
int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		printf("./tcd [graph] [query]\n");
		return 1;
	}
	initmem();
	_initmem();
	loadgraph(argv[1]);

#ifdef CALC_PRUNE
	for (int i = 0; i < MAX_SPAN; ++i)
		for (int j = 0; j < MAX_SPAN; ++j)
			(*prune_flag)[i][j] = 0;
#endif

	loadtest(argv[2]);
	printf("load graph complete\n");
	for (int i = 0; i < qcnt; ++i)
	{
		int l = q[i].l, r = q[i].r, k = q[i].k;
		buildtel(l, r);
		initMH(l, r);
		memset(cntinfo, 0, sizeof cntinfo);
		auto t0 = system_clock::now();
		tcd(l, r, k);
		auto t1 = system_clock::now();
		auto clapse_ns = duration_cast<nanoseconds>(t1 - t0);
		int tcdcell = cntinfo[INFO::CELL];
		buildtel(l, r);
		initMH(l, r);
		memset(cntinfo, 0, sizeof cntinfo);
		auto ot0 = system_clock::now();
		otcd(l, r, k);
		auto ot1 = system_clock::now();
		auto oclapse_ns = duration_cast<nanoseconds>(ot1 - ot0);

#ifdef CALC_PRUNE
		calc_prune();
#endif

		pinfo(i, argv[1], tcdcell, clapse_ns.count(), oclapse_ns.count());
		printf("\n");

/*
		proctti();
*/

	}

	system("pause");
	return 0;
}