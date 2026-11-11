#define main census_original_main
#include "census.cpp"
#undef main
#include <sstream>
using I128=__int128_t;
vector<int64_t> lglo,lghi;
int64_t scale;
struct Interval {long double lo,hi;};
struct Certified {Interval Hratio,Eratio;};
long double outward(I128 num,int64_t den,bool upper){
 long double direction=upper?INFINITY:-INFINITY;
 long double x=nextafterl((long double)num,direction);
 x=nextafterl(x/scale,direction);
 return nextafterl(x/den,direction);
}
Certified certified(const Family&a){
 int N=1<<n,M=N*N;
 vector<int64_t> f(N,1);for(int x=0;x<N;x++)for(int t:a)if((x&t)==t){f[x]=-1;break;}
 wht(f);vector<int> m(N);int64_t influence=0;I128 hlo=0,hhi=0;
 for(int x=0;x<N;x++){
  m[x]=int(f[x]*f[x]);influence+=pc(x)*m[x];
  hlo+=I128(m[x])*(2*n*scale-lghi[m[x]]);
  hhi+=I128(m[x])*(2*n*scale-lglo[m[x]]);
 }
 Interval er={1e20L,1e20L};
 for(int j=0;j<n;j++){
  int64_t delta=0;I128 elo=0,ehi=0;
  for(int x=0;x<N;x++)if(!(x&(1<<j))){int p=m[x],q=m[x|(1<<j)],s=p+q;delta+=q;
   elo+=I128(s)*lglo[s]-I128(p)*lghi[p]-I128(q)*lghi[q];
   ehi+=I128(s)*lghi[s]-I128(p)*lglo[p]-I128(q)*lglo[q];
  }
  if(delta){er.lo=min(er.lo,outward(elo,delta,false));er.hi=min(er.hi,outward(ehi,delta,true));}
 }
 return {{outward(hlo,influence,false),outward(hhi,influence,true)},er};
}
int main(int argc,char**argv){
 ifstream logs("log2_intervals.txt");logs>>scale;int64_t a,b;while(logs>>a>>b){lglo.push_back(a);lghi.push_back(b);}if(lglo.size()!=65537)return 2;
 n=atoi(argv[1]);string prefix=argv[2];ifstream input(prefix+".families");string line;
 vector<pair<Family,Certified>> topH,topE;long double bestHl=-1,bestEl=-1,otherHu=-1,otherEu=-1;int64_t count=0;
 // Candidates within 1e-9 of a leading lower bound are retained. Every
 // discarded candidate's certified upper endpoint is recorded separately.
 while(getline(input,line)){
  istringstream is(line);Family f;int t;while(is>>t)f.push_back(t);auto c=certified(f);count++;
  if(c.Hratio.hi+1e-9L>=bestHl){
   if(c.Hratio.lo>bestHl+1e-9L){for(auto &x:topH)otherHu=max(otherHu,x.second.Hratio.hi);topH.clear();}topH.push_back({f,c});bestHl=max(bestHl,c.Hratio.lo);
  }else otherHu=max(otherHu,c.Hratio.hi);
  if(c.Eratio.hi+1e-9L>=bestEl){
   if(c.Eratio.lo>bestEl+1e-9L){for(auto &x:topE)otherEu=max(otherEu,x.second.Eratio.hi);topE.clear();}topE.push_back({f,c});bestEl=max(bestEl,c.Eratio.lo);
  }else otherEu=max(otherEu,c.Eratio.hi);
 }
 ofstream o(prefix+".certified.json");o<<setprecision(20)<<"{\"n\":"<<n<<",\"count\":"<<count;
 auto print=[&](string name,auto &top,long double other){o<<",\""<<name<<"\":{\"other_upper\":"<<other<<",\"candidates\":[";bool first=true;for(auto &x:top){if(!first)o<<',';first=false;o<<"{\"terms\":[";for(size_t j=0;j<x.first.size();j++){if(j)o<<',';o<<x.first[j];}auto v=name=="H_I"?x.second.Hratio:x.second.Eratio;o<<"],\"lower\":"<<v.lo<<",\"upper\":"<<v.hi<<"}";}o<<"]}";};
 print("H_I",topH,otherHu);print("min_E_Inf",topE,otherEu);o<<"}\n";
}
