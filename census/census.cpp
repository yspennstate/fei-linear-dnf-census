#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <sys/resource.h>
#include <unordered_set>
#include <vector>
using namespace std;
using Family=vector<int>;
int n;
int pc(int a){return __builtin_popcount((unsigned)a);}
string encode(const Family &a){string s;for(int x:a){s.push_back(x&255);s.push_back((x>>8)&255);}return s;}
Family decode(const string &s){Family a;for(size_t j=0;j<s.size();j+=2)a.push_back((unsigned char)s[j]+256*(unsigned char)s[j+1]);return a;}
// Exact individualization/refinement. No hashes or graph shadows substitute for equality.
struct Canon {
 Family edges; string best; bool initialized=false; uint64_t aut_count=0; vector<vector<int>> twins;
 Canon(const Family &e):edges(e),twins(n,vector<int>(n,0)){
  set<int> es(e.begin(),e.end());
  for(int i=0;i<n;i++)for(int j=0;j<i;j++){
   bool ok=true;for(int a:e){int b=a;if(((a>>i)^(a>>j))&1)b^=(1<<i)|(1<<j);if(!es.count(b)){ok=false;break;}}
   twins[i][j]=twins[j][i]=ok;
  }
 }
 vector<vector<int>> refine(vector<vector<int>> cells){
  while(true){
   vector<int> color(n);for(int c=0;c<(int)cells.size();c++)for(int v:cells[c])color[v]=c;
   vector<vector<int>> out;
   for(auto cell:cells){
    map<vector<uint64_t>,vector<int>> groups;
    for(int v:cell){
     vector<uint64_t> sig;
     for(int a:edges)if(a&(1<<v)){
      uint64_t t=0;for(int u=0;u<n;u++)if(a&(1<<u))t+=uint64_t(1)<<(4*color[u]);
      sig.push_back(t);
     }
     sort(sig.begin(),sig.end());groups[sig].push_back(v);
    }
    for(auto &g:groups)out.push_back(g.second);
   }
   if(out.size()==cells.size())return out;cells=move(out);
  }
 }
 void visit(vector<vector<int>> cells,uint64_t weight=1){
  cells=refine(move(cells));int c=-1;
  for(int j=0;j<(int)cells.size();j++)if(cells[j].size()>1&&(c<0||cells[j].size()<cells[c].size()))c=j;
  if(c<0){
   int mp[9];for(int j=0;j<n;j++)mp[cells[j][0]]=j;
   Family a;for(int x:edges){int y=0;for(int j=0;j<n;j++)if(x&(1<<j))y|=1<<mp[j];a.push_back(y);}
   sort(a.begin(),a.end());string s=encode(a);if(!initialized||s<best){best=s;initialized=true;aut_count=weight;}else if(s==best)aut_count+=weight;return;
  }
  vector<int> representatives;
  for(int v:cells[c]){
   bool skip=false;for(int u:representatives)if(twins[v][u]){skip=true;break;}if(skip)continue;representatives.push_back(v);
   auto sub=cells;vector<int> rest;for(int u:cells[c])if(u!=v)rest.push_back(u);
   uint64_t multiplicity=1;for(int u:cells[c])if(u!=v&&twins[v][u])multiplicity++;
   sub[c]={v};sub.insert(sub.begin()+c+1,rest);visit(move(sub),weight*multiplicity);
  }
 }
 string run(){vector<int>a(n);iota(a.begin(),a.end(),0);visit({a});return best;}
};
string canon(const Family &a){return Canon(a).run();}
bool compatible(const Family &a,int e){for(int t:a)if(pc(t&e)>=2)return false;return true;}
void wht(vector<int64_t>&a){for(int h=1;h<(int)a.size();h*=2)for(int k=0;k<(int)a.size();k+=2*h)for(int j=0;j<h;j++){auto u=a[k+j],v=a[k+j+h];a[k+j]=u+v;a[k+j+h]=u-v;}}
long double eta(long double x){return x>0?-x*log2l(x):0;}
long double phi(long double p,long double q){return eta(p)+eta(q)-eta(p+q);}
struct Record {long double value=-1; vector<Family> families;};
void update(Record&r,long double x,const Family&a){if(x>r.value+1e-14L){r.value=x;r.families={a};}else if(fabsl(x-r.value)<=1e-14L)r.families.push_back(a);}
struct Metrics {int live,ones;long double h,inf,minr;};
Metrics metrics(const Family&a){
 int mask=0;for(int t:a)mask|=t;int live=pc(mask),N=1<<n,ones=0;
 vector<int64_t> f(N,1);for(int x=0;x<N;x++)for(int t:a)if((x&t)==t){f[x]=-1;ones++;break;}
 wht(f);long double denom=(long double)N*N,H=0,I=0,mr=1e10L;
 for(int x=0;x<N;x++){long double q=(long double)f[x]*f[x]/denom;H+=eta(q);I+=pc(x)*q;}
 for(int j=0;j<n;j++)if(mask&(1<<j)){
  long double e=0,d=0;
  for(int x=0;x<N;x++)if(!(x&(1<<j))){long double p=(long double)f[x]*f[x]/denom,q=(long double)f[x|(1<<j)]*f[x|(1<<j)]/denom;e+=phi(p,q);d+=q;}
  mr=min(mr,e/d);
 }
 return {live,ones,H,I,mr};
}
void print_families(ostream&o,const Record&r){
 o<<"{\"value\":"<<setprecision(18)<<r.value<<",\"attainers\":[";
 for(size_t k=0;k<r.families.size();k++){
  if(k)o<<",";auto &a=r.families[k];auto m=metrics(a);
  o<<"{\"terms\":[";for(size_t j=0;j<a.size();j++){if(j)o<<",";o<<a[j];}
  o<<"],\"true_count\":"<<m.ones<<",\"truth_size\":"<<(1<<n)<<",\"H\":"<<m.h<<",\"I\":"<<m.inf<<",\"min_E_Inf\":"<<m.minr<<"}";
 }
 o<<"]}";
}
int main(int argc,char**argv){
 n=argc>1?atoi(argv[1]):8;string prefix=argc>2?argv[2]:"census";if(n<1||n>9)return 2;
 int priority=getpriority(PRIO_PROCESS,0);cerr<<"priority "<<priority<<"\n";if(priority<10){cerr<<"requires inherited nice>=10\n";return 2;}
 auto start=chrono::steady_clock::now();unordered_set<string> all;all.insert("");
 ofstream log(prefix+".log");
 for(int w=n;w>=2;w--){
  vector<int> candidates;for(int e=1;e<(1<<n);e++)if(pc(e)==w)candidates.push_back(e);
  vector<string> layer(all.begin(),all.end());int depth=0;
  while(!layer.empty()){
   unordered_set<string> next;
   for(auto &s:layer){auto a=decode(s);for(int e:candidates)if(compatible(a,e)){
    a.push_back(e);auto key=canon(a);a.pop_back();if(!all.count(key))next.insert(move(key));
   }}
   for(auto&s:next)all.insert(s);layer.assign(next.begin(),next.end());depth++;
   log<<"width "<<w<<" depth "<<depth<<" new "<<layer.size()<<" total "<<all.size()<<" seconds "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<endl;
   cerr<<"width "<<w<<" depth "<<depth<<" new "<<layer.size()<<" total "<<all.size()<<"\n";
  }
 }
 ofstream data(prefix+".families");vector<long long>counts(n+1);vector<Record> rE(n+1),rH(n+1);
 for(auto&s:all){auto a=decode(s);int used=0;for(int t:a)used|=t;vector<int>isolates;for(int j=0;j<n;j++)if(!(used&(1<<j)))isolates.push_back(j);
  for(int k=0;k<=(int)isolates.size();k++){
   if(k)a.push_back(1<<isolates[k-1]);if(a.empty())continue;
   auto m=metrics(a);counts[m.live]++;update(rE[m.live],m.minr,a);update(rH[m.live],m.h/m.inf,a);
   if(m.live==n){auto b=a;sort(b.begin(),b.end());for(int t:b)data<<t<<' ';data<<'\n';}
  }
 }
 ofstream json(prefix+".json");json<<"{\"n_max\":"<<n<<",\"priority\":"<<priority<<",\"non_singleton_including_isolates\":"<<all.size()<<",\"rows\":[";
 for(int j=1;j<=n;j++){if(j>1)json<<",";json<<"{\"n\":"<<j<<",\"count\":"<<counts[j]<<",\"max_min_E_Inf\":";print_families(json,rE[j]);json<<",\"max_H_I\":";print_families(json,rH[j]);json<<"}";}
 json<<"],\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";
 return 0;
}
