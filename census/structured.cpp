#define main census_original_main
#include "census.cpp"
#undef main
#include <sstream>
struct Detailed {long double H,I;int64_t ones;vector<long double> delta,erasure;};
Detailed fullstats(int dim,const Family &ts,const vector<int>& reps){
 int N=1<<dim;vector<int64_t> f(N,1);int64_t ones=0;
 for(int x=0;x<N;x++)for(int t:ts)if((x&t)==t){f[x]=-1;ones++;break;}
 wht(f);long double M=(long double)N*N,H=0,I=0;
 for(int x=0;x<N;x++){long double mu=(long double)f[x]*f[x]/M;H+=eta(mu);I+=pc(x)*mu;}
 vector<long double>d,e;
 for(int j:reps){long double ds=0,es=0;for(int x=0;x<N;x++)if(!(x&(1<<j))){long double p=(long double)f[x]*f[x]/M,q=(long double)f[x|(1<<j)]*f[x|(1<<j)]/M;ds+=q;es+=phi(p,q);}d.push_back(ds);e.push_back(es);}
 return {H,I,ones,d,e};
}
Family clean(Family a){sort(a.begin(),a.end(),[](int a,int b){return pc(a)==pc(b)?a<b:pc(a)<pc(b);});Family out;for(int t:a){bool absorbed=false;for(int s:out)if((s&t)==s){absorbed=true;break;}if(!absorbed)out.push_back(t);}return out;}
Family restrict_terms(const Family &a,int j,int b){Family out;for(int t:a){if(t&(1<<j)){if(!b)continue;t^=1<<j;}t=(t&((1<<j)-1))|((t>>(j+1))<<j);out.push_back(t);}return clean(out);}
struct Cost {long double E=0,I=0;uint64_t nodes=0;};
Cost tree(int dim,Family ts,vector<int> pending={}){
 ts=clean(ts);if(ts.empty()||ts[0]==0)return {};
 if(pending.empty()){int t=ts[0];for(int j=0;j<dim;j++)if(t&(1<<j))pending.push_back(j);}
 int j=pending[0];pending.erase(pending.begin());for(int &k:pending)k-=k>j;
 int N=1<<dim;vector<int64_t> f(N,1);for(int x=0;x<N;x++)for(int t:ts)if((x&t)==t){f[x]=-1;break;}
 wht(f);long double M=(long double)N*N,E=0,I=0;
 for(int x=0;x<N;x++)if(!(x&(1<<j))){long double p=(long double)f[x]*f[x]/M,q=(long double)f[x|(1<<j)]*f[x|(1<<j)]/M;E+=phi(p,q);I+=q;}
 vector<int64_t>().swap(f);
 auto c0=tree(dim-1,restrict_terms(ts,j,0),pending),c1=tree(dim-1,restrict_terms(ts,j,1),pending);
 return {E+(c0.E+c1.E)/2,I+(c0.I+c1.I)/2,1+c0.nodes+c1.nodes};
}
int main(int argc,char**argv){
 int width=atoi(argv[1]),q=atoi(argv[2]);int dim=width*q;if(dim>25||width>q||q>7)return 2;
 if(getpriority(PRIO_PROCESS,0)<10)return 3;
 auto start=chrono::steady_clock::now();Family ts;
 for(int a=0;a<q;a++)for(int b=0;b<q;b++){int t=0;for(int x=0;x<width;x++)t|=1<<(x*q+(a*x+b)%q);ts.push_back(t);}
 bool linear=true;for(size_t i=0;i<ts.size();i++)for(size_t j=0;j<i;j++)linear&=pc(ts[i]&ts[j])<=1;if(!linear)return 4;
 vector<int> parent(width);iota(parent.begin(),parent.end(),0);auto find=[&](int x){while(parent[x]!=x)x=parent[x];return x;};
 for(int a=1;a<q;a++)for(int b=0;b<q;b++){bool ok=true;for(int x=0;x<width;x++)if((a*x+b)%q>=width)ok=false;if(ok)for(int x=0;x<width;x++)parent[find(x)]=find((a*x+b)%q);}
 vector<int>reps;for(int x=0;x<width;x++)if(find(x)==x)reps.push_back(x*q);
 auto s=fullstats(dim,ts,reps);string file="TD_"+to_string(width)+"_"+to_string(q)+".json";ofstream o(file);o<<setprecision(18)<<"{\"width\":"<<width<<",\"q\":"<<q<<",\"n\":"<<dim<<",\"terms_count\":"<<ts.size()<<",\"degree\":"<<q<<",\"load\":"<<q*powl(2,1-width)<<",\"linear\":true,\"true_count\":"<<s.ones<<",\"truth_size\":"<<(1<<dim)<<",\"mean_sign\":"<<1-2*(long double)s.ones/(1<<dim)<<",\"H\":"<<s.H<<",\"I\":"<<s.I<<",\"H_I\":"<<s.H/s.I<<",\"coordinate_orbits\":[";
 for(size_t k=0;k<reps.size();k++){if(k)o<<",";o<<"{\"representative\":"<<reps[k]<<",\"Inf\":"<<s.delta[k]<<",\"E_Inf\":"<<s.erasure[k]/s.delta[k]<<"}";}o<<"]";
 if(dim<=25){auto c=tree(dim,ts);o<<",\"whole_term_tree_E_I\":"<<c.E/s.I<<",\"whole_term_tree_influence_sum\":"<<c.I<<",\"tree_nodes\":"<<c.nodes;if(fabsl(c.I-s.I)>1e-12L)return 5;}
 else{o<<",\"any_tree_E_I_lower\":"<<max(0.L,s.H/s.I-2/logl(2.L))<<",\"any_tree_E_I_upper\":"<<s.H/s.I;}
 struct rusage r;getrusage(RUSAGE_SELF,&r);o<<",\"peak_rss_bytes\":"<<r.ru_maxrss<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";
}
