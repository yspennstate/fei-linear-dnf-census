#define main census_original_main
#include "census.cpp"
#undef main
#include <sstream>
map<string,uint64_t> memo;
vector<int> pairs;
uint64_t pairmask(int vertexmask){uint64_t b=0;for(size_t j=0;j<pairs.size();j++)if((pairs[j]&vertexmask)==pairs[j])b|=uint64_t(1)<<j;return b;}
uint64_t labeled(uint64_t available){
 if(!available)return 1;
 Family g;for(size_t j=0;j<pairs.size();j++)if(available&(uint64_t(1)<<j))g.push_back(pairs[j]);
 string key=canon(g);auto it=memo.find(key);if(it!=memo.end())return it->second;
 int first=__builtin_ctzll(available),uv=pairs[first];uint64_t answer=2*labeled(available^(uint64_t(1)<<first));
 for(int t=1;t<(1<<n);t++)if(pc(t)>=3&&(t&uv)==uv){uint64_t b=pairmask(t);if((b&available)==b)answer+=labeled(available^b);}
 memo[key]=answer;return answer;
}
int main(int argc,char**argv){
 n=atoi(argv[1]);string prefix=argv[2];for(int t=1;t<(1<<n);t++)if(pc(t)==2)pairs.push_back(t);
 uint64_t expected=labeled((uint64_t(1)<<pairs.size())-1),observed=0,classes=0,graphs=0,factorial=1;
 for(int j=2;j<=n;j++)factorial*=j;
 ifstream input(prefix+".families");string line;
 // Removing singleton terms from each essential-n representative gives
 // precisely all hypergraphs with edges >=2, with isolated vertices allowed.
 while(getline(input,line)){
  istringstream is(line);Family a;int t;while(is>>t)if(pc(t)>=2)a.push_back(t);Canon c(a);c.run();
  if(!c.aut_count||factorial%c.aut_count)return 3;
  observed+=factorial/c.aut_count;classes++;bool graph=true;for(int t:a)graph&=pc(t)==2;graphs+=graph;
 }
 ofstream o(prefix+".labeled_check.json");o<<"{\"n\":"<<n<<",\"classes\":"<<classes<<",\"graph_classes\":"<<graphs<<",\"pair_recursion_labeled_count\":"<<expected<<",\"orbit_stabilizer_labeled_count\":"<<observed<<",\"pair_recursion_memo_states\":"<<memo.size()<<",\"match\":"<<(expected==observed?"true":"false")<<"}\n";
 return expected==observed?0:4;
}
