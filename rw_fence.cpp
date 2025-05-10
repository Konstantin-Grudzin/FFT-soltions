/*
  _                            _
 | |                          | |
 | | _____  _ __  _ __ ___  __| |
 | |/ / _ \| '_ \| '__/ _ \/ _` |
 |   < (_) | | | | | |  __/ (_| |
 |_|\_\___/|_| |_|_|  \___|\__,_|

*/
#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
#include <numeric>
#include <unordered_map>

#include <ext/pb_ds/assoc_container.hpp> 
#include <ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds; 

template<typename T>
using ordered_set = tree<T, null_type,std::less<T>, rb_tree_tag,tree_order_statistics_node_update>;
 

// FOR HASH
#include <random>
#include <chrono>
#include <utility>
//

//#pragma comment(linker, "/STACK:400000000")

#define all(a) begin(a),end(a)
#define rall(a) rbegin(a),rend(a)
#define ACspeed ios_base::sync_with_stdio(false);cin.tie(NULL);cout.tie(NULL);
#define ll long long
#define ull  unsigned long long
#define ld  long double
#define pld pair<double,double>
#define yeah cout<<"Yes\n"
#define nah cout<<"No\n"
#define YEAH cout<<"YES\n"
#define NAH cout<<"NO\n"
#define INF 1e18+100
#define iINF 1e9+100
#define pll pair<ll,ll>
#define pii pair<int,int>
#define pull pair<ull,ull>
#define PI M_PI
#define X first
#define Y second
#define ff first
#define ss second
#define amax(x,y) x=maxn(x,y)
#define amin(x,y)  x=min(x,y)


//user vectors
template<typename T>
using v1 = std::vector<T>;

template<typename T>
using v2 = std::vector<v1<T>>;

template<typename T>
using v3 = std::vector<v2<T>>;

template<typename T>
using v4 = std::vector<v3<T>>;

using namespace std;

// --- PRINT REALISATION ---
//inline void print() { cout << endl; }
//template<typename T,typename... args>
//void print(T&& obj, args&&... Args)
//{
//	cout << obj << " ";
//	print(Args...);
//}

// --- NUMBER THEORY STUFF ---

ll is_square(ll x)
{
	ll ans = sqrt(x);
	if (ans * ans == x)
		return 1;
	else
		return 0;
}

ll mod;

template<class T>
inline T Mod(T ans, T n)
{
	return ((ans % n) + n) % n;
}

template<class T>
T gcd(T a, T b)
{
	if (b == 0)
		return a;
	return gcd(b, a % b);
}


ll phi(ll n)
{
	ll ans = n;
	for (ll i = 2; i * i <= n; ++i)
	{
		if (n % i == 0)
		{
			while (n % i == 0) { n /= i; }
			ans -= ans / i;
		}
	}
	if (n > 1)
		ans -= ans / n;
	return ans;
}

template<class T>
T pow2(T x, T n, T mod1 = 0) {
	T base = x, ans = 1;
	while (n > 0) {
		if (n & 1) {
			ans *= base;
			if (mod1)
				ans = Mod(ans, mod1);
			n--;
		}
		else {
			base *= base;
			if (mod1)
				base = Mod(base, mod1);
			n /= 2;
		}
	}
	return ans;
}

inline ll Div(ll a, ll b)
{
	return Mod(a * pow2(b, mod - 2, mod), mod);
}

template<class T>
T mulmod(T x, T y, T m)
{
	x %= m;
	T ans = 0;
	while (y > 0)
	{
		if (y & 1)
		{
			ans = (ans + x) % m;
		}
		y >>= 1;
		x = (x << 1) % m;
	}
	return ans;
}
// ---         END         ---

// ---     TEMPLATE I/0    ---

template<class T, class U>
istream& operator>>(istream& in, pair<T, U>& v) {
	in >> v.first >> v.second;
	return in;
}

template<class T, class U>
ostream& operator<<(ostream& out, pair<T, U>& v) {
	out << v.first << " " << v.second;
	return out;
}

template<class T>
istream& operator>>(istream& in, vector <T>& v) {
	for (auto& x : v) { in >> x; }
	return in;
}

template<class T>
ostream& operator<<(ostream& out, vector <T>& v) {
	for (auto& x : v) { out << x << " "; }
	return out;
}
// ---         END         ---

//use offset=0 if you construct prefix from another prefix
void pref_sum(v1<ll>& source, v1<ll>& dest, int offset = 1)
{
	dest.resize(source.size() + offset);
	partial_sum(all(source), begin(dest) + offset);
}

#define wr cout<<"-1\n"

const ll sz=1e6;
v1<ll> fac(sz);
v1<ll> rfac(sz);
ll C(ll n, ll k)
{
	if (n <( k)) return 0;
	ll x = fac[n]*rfac[k]%mod;
	x = x*rfac[n - k]%mod;
	return x;
}

inline ll bit(ll x,ll pos)
{
	return (x>>pos)&1ll;
}

#define fr(i,n,m) for(ll i=n;i<m;++i)
#define re return
//#define make_unique(a) a.resize(unique(begin(a),end(a))-begin(a))

typedef complex<double> Comp;

void fft(vector<Comp> &a) {
    int n = a.size(), l = 31 - __builtin_clz(n);
    vector<complex<long double> > R(2, {1, 0});
    vector<Comp> rt(2, {1, 0});
    for (int k = 2; k < n; k *= 2) {
        R.resize(n), rt.resize(n);
        complex<long double> x = polar(1.0L, acos(-1.0L) / k);
        for (int i = k; i < 2 * k; i++)
            rt[i] = R[i] = i & 1 ? R[i / 2] * x : R[i / 2];
    }
    vector<int> rev(n);
    for (int i = 0; i < n; i++)
        rev[i] = (rev[i / 2] | (i & 1) << l) / 2;
    for (int i = 0; i < n; i++)
        if (i < rev[i])
            swap(a[i], a[rev[i]]);
    for (int k = 1; k < n; k *= 2)
        for (int i = 0; i < n; i += 2 * k)
            for (int j = 0; j < k; j++) {
                auto *x = (double *) &rt[j + k], *y = (double *) &a[i + j + k];
                Comp z(x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0]);
                a[i + j + k] = a[i + j] - z;
                a[i + j] += z;
            }
}

vector<ll> conv(const vector<ll> &a, const vector<ll> &b) {
    if (a.empty() || b.empty()) return {};
    vector<ll> res(a.size() + b.size() - 1);
    int B = 32 - __builtin_clz(res.size()), n = 1 << B, cut = int(sqrt(mod));
    vector<Comp> L(n), R(n), outs(n), outl(n);
    for (int i = 0; i < a.size(); i++) L[i] = Comp( a[i] / cut,  a[i] % cut);
    for (int i = 0; i < b.size(); i++) R[i] = Comp( b[i] / cut,  b[i] % cut);
    fft(L), fft(R);
    for (int i = 0; i < n; i++) {
        int j = -i & (n - 1);
        outl[j] = (L[i] + conj(L[j])) * R[i] / (2.0 * n);
        outs[j] = (L[i] - conj(L[j])) * R[i] / (2.0 * n) / Comp{0, 1};
    }
    fft(outl), fft(outs);
    for (int i = 0; i < res.size(); i++) {
        ll av = llround(outl[i].real()), cv = llround(outs[i].imag());
        ll bv = llround(outl[i].imag()) + llround(outs[i].real());
        res[i] = ((av % mod * cut + bv) % mod * cut + cv) % mod;
    }
    return res;
}

vector<ll> pow2Poly(const vector<ll>& a,ll x)
{
	v1<ll> ans={1};
	v1<ll> base=a;
	while(x>0)
	{
		if(x%2==0)
		{
			base=conv(base,base);
			x>>=1;
		}
		ans=conv(base,ans);
		x--;
	}
	return ans;
}

void solve() 
{
	ll n,k;cin>>n>>k;
	map<ll,ll> cnt;
	v1<ll> a(n);cin>>a;
	for(auto el:a)
	{
		cnt[el]++;
	}

	v1<ll> b(k);cin>>b;
	ll q;cin>>q;
	v1<ll> Q(q);cin>>Q;
	v1<ll> ans(q);
	auto getAns=[&](ll x)->void
	{
		ll o=0,d=0;
		for(auto [v,k]:cnt)
		{
			if(v>=x) break;
			if(k==1)
				o++;
			else
				d++;
		}
		auto c = conv(pow2Poly({1,2},o),pow2Poly({1,2,1},d));
		for(ll i=0;i<q;++i)
		{
			ll el=Q[i];
			el-=2*x+2;
			el>>=1;
			if(el<0 || el>=c.size()) continue;
			ans[i]+=c[el];
			ans[i]%=mod;
		}
	};
	for(auto el:b)
		getAns(el);
	for(auto el:ans)
		cout<<el<<endl;
}
//#define FAC
//#define MULTEST

void before_solve()
{
	mod = 998244353;
#ifdef FAC
	fac[0]=1;
	for(ll i=1;i<=sz;++i)
	{
		fac[i] = fac[i - 1] * i % mod;
	}
#endif
}

signed main()
{
	//#define FILEIO
#ifdef FILEIO
	freopen("input.txt", "r", stdin); freopen("output.txt", "w", stdout);
#endif
	//PolyHash::base = gen_base(256, PolyHash::mod);
	ACspeed;
	ll t = 1;
#ifdef MULTEST
	cin >> t;
#endif
	cout << fixed << setprecision(50);
	before_solve();
	while (t--)
	{
		solve();
	}
}
