#include <iostream>
#include <vector>
#include <ext/hash_set>

using namespace std;
using namespace __gnu_cxx;

int main() {
  hash< vector<int> > HV;
  hash<const char*> H;
  cout << "foo -> " << H("foo shd dgshd dhdgs dgdgs dgd sgd dgs dgd") << endl;
  cout << "bar -> " << H("bar ehwe ejhe ejwhe ejhw ejhee ejhw ejhe ejhw") << endl;
  int N=10;
  vector<int> v1(N), v2;
  for (int j=0; j<N; j++) v1[j] = j*j;
  v2 = v1;
  cout << "hash de v1: " << HV(v1) << endl;
  cout << "hash de v2: " << HV(v2) << endl;
}
