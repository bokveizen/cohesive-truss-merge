#include <algorithm>
#include <bitset>
#include <cassert>
#include <chrono>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <queue>
#include <random>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "tqdm.h"

using namespace std;

// types
typedef struct {
    int u, v;
} TEdge;

bool operator<(TEdge e1, TEdge e2) {
    return e1.u < e2.u || (e1.u == e2.u && e1.v < e2.v);
}

bool operator==(TEdge e1, TEdge e2) {
    return (e1.u == e2.u) && (e1.v == e2.v);
}

struct VectorHasher {
    auto operator()(const vector<int> &V) const {
        auto hash = V.size();
        for (auto &i: V) {
            hash ^= i + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

struct TEdgeHasher {
    auto operator()(const TEdge &E) const {
        int eu = E.u;
        int ev = E.v;
        long long euv = ((long long) eu << 32) + (long long) ev;
        return hash<long long>{}(euv);
    }
};

typedef unordered_map<int, int> MII;
typedef vector<MII> VMII;
typedef vector<TEdge> VE;
typedef unordered_set<TEdge, TEdgeHasher> SE;
typedef vector<int> VI;

// variables
random_device rd;
mt19937 gen(rd());
ifstream fin;
ofstream fout;
string infile, outfile;
int n, m, nTotal, mTotal, kInput, mTk, mTkm1, mOutside, nTk, nTkm1, nOutside;
VI mapto, deg, bin, outsideNodes, Tkm1Nodes, TkNodes, Skm1Nodes, params, uniqueOutsideNodes, bestOutsideNodes;
VE binEdge, outsideEdges, Tkm1Edges, TkEdges, Skm1Edges, unstableEdges, almostUnstableEdges;
vector<VI> A, nbrsInTkm1, nbrsInTk, vTwoNewNbrs;
VMII supp, suppOrig, suppTkm1, suppTk, pos;
unordered_map<TEdge, int, TEdgeHasher> edge2trussness;
unordered_map<int, VI> v2HNNs, v2UsNbrs, v2ShNbrs;
unordered_set<int> nodesWithTkm1Nbrs;
unordered_map<TEdge, vector<TEdge>, TEdgeHasher> e2HUEs;
MII v2HUEs, v2HSEs;
chrono::time_point<chrono::system_clock> startTime, endTime;
char timeStr[100];
const int maxN = 200000;
bool debugMode_ = false;

inline bool compVertex(int i, int j) {
    return deg[i] < deg[j] || (deg[i] == deg[j] && i < j);
}

inline void orderPair(int &u, int &v) {
    if (!compVertex(u, v)) swap(u, v);
}

inline void updateSupport(int u, int v, int delta, VMII &suppInput = supp) {
    suppInput[u][v] += delta;
    suppInput[v][u] += delta;
}

inline void removeEdge(int u, int v, VMII &suppInput) {
    suppInput[u].erase(v);
    suppInput[v].erase(u);
}

void readOrderedSimpleGraph() {
    fin.open(infile.c_str());
    mTk = mTkm1 = mOutside = nTk = nTkm1 = nOutside = 0;
    fin >> n >> m;
    nTotal = n;
    mTotal = m;
    if (debugMode_) cout << nTotal << " nodes and " << mTotal << " edges in total" << endl;
    int u, v;
    deg.resize(n, 0);
    vTwoNewNbrs.resize(n);
    nbrsInTkm1.resize(n);
    nbrsInTk.resize(n);
    supp.resize(n);
    for (int i = 0; i < n; ++i) supp[i].clear();
    for (int i = 0; i < m; ++i) {
        fin >> u >> v;
        TEdge euv = {min(u, v), max(u, v)};
        edge2trussness[euv] = kInput;
        supp[u][v] = 0;
        supp[v][u] = 0;
        ++deg[u];
        ++deg[v];
    }
    mapto.resize(n);
    for (int i = 0; i < n; ++i) mapto[i] = i;
    sort(mapto.begin(), mapto.end(), compVertex);
    fin.close();
}

void intersect(const VI &a, const VI &b, VI &c) {
    c.reserve(min(a.size(), b.size()));
    unsigned j = 0;
    for (int i: a) {
        while (j < b.size() && b[j] > i) ++j;
        if (j >= b.size()) break;
        if (b[j] == i) c.emplace_back(i);
    }
}

void countTriangles(vector<MII> &suppInput = supp) {
    A.resize(n);
    for (int i = 0; i < n; ++i) A[i].clear();
    int nDeltas = 0;
    for (int vi = n - 1; vi >= 0; --vi) {
        int v = mapto[vi];
        for (auto it = suppInput[v].begin(); it != suppInput[v].end(); ++it) {
            int u = it->first;
            if (!compVertex(u, v)) continue;
            VI common;
            intersect(A[u], A[v], common);
            for (int i: common) {
                int w = mapto[i];
                ++nDeltas;
                updateSupport(u, v, 1, suppInput);
                updateSupport(v, w, 1, suppInput);
                updateSupport(w, u, 1, suppInput);
            }
            A[u].emplace_back(vi);
        }
    }
    if (debugMode_) cout << nDeltas << " triangles found." << endl;
}

void binSort(VMII &suppInput = supp, int &mRemoved = mOutside, bool updateTrussness = true) {
    bin.clear();
    bin.resize(n, 0);
    int nBins = 0;
    int mp = 0;
    for (int u = 0; u < n; ++u) {
        MII suppU = suppInput[u];
        for (auto &it: suppU) {
            int v = it.first;
            if (!compVertex(u, v)) continue;
            int supUV = it.second;
            if (supUV == 0) {
                if (updateTrussness) {
                    TEdge euv = {min(u, v), max(u, v)};
                    edge2trussness[euv] = 2;
                }
                removeEdge(u, v, suppInput);
                ++mRemoved;
                continue;
            }
            ++mp;
            ++bin[supUV];
            nBins = max(supUV, nBins);
        }
    }
    m = mp;
    ++nBins;
    int count = 0;
    for (int i = 0; i < nBins; ++i) {
        int binsize = bin[i];
        bin[i] = count;
        count += binsize;
    }
    pos.clear();
    pos.resize(n);
    for (int i = 0; i < n; ++i) pos[i].clear();
    binEdge.resize(m);
    for (int u = 0; u < n; ++u) {
        for (auto it = suppInput[u].begin(); it != suppInput[u].end(); ++it) {
            int v = it->first;
            if (!compVertex(u, v)) continue;
            int sup = it->second;
            TEdge e = {u, v};
            int &b = bin[sup];
            binEdge[b] = e;
            pos[u][v] = b++;
        }
    }
    for (int i = nBins; i > 0; --i) bin[i] = bin[i - 1];
    bin[0] = 0;
}

void updateEdge(int u, int v, int minsup, VMII &suppInput) {
    orderPair(u, v);
    int sup = suppInput[u][v];
    if (sup <= minsup) return;
    int p = pos[u][v];
    int posbin = bin[sup];
    TEdge se = binEdge[posbin];
    TEdge e = {u, v};
    if (p != posbin) {
        pos[u][v] = posbin;
        pos[se.u][se.v] = p;
        binEdge[p] = se;
        binEdge[posbin] = e;
    }
    ++bin[sup];
    updateSupport(u, v, -1, suppInput);
}

void trussDecomp() {
    bool flag = true;
    for (int s = 0; s < m; ++s) {
        int u = binEdge[s].u;
        int v = binEdge[s].v;
        orderPair(u, v);
        int supUV = supp[u][v];
        int tuv = supUV + 2;
        TEdge euv = {min(u, v), max(u, v)};
        if (tuv >= kInput) {
            suppTk = supp;
            mTk = mTotal - mTkm1 - mOutside;
            mTkm1 += mTk;
            return;
        } else if (tuv >= kInput - 1) {
            if (flag) {
                suppTkm1 = supp;
                flag = false;
            }
            edge2trussness[euv] = tuv;
            ++mTkm1;
        } else {
            edge2trussness[euv] = tuv;
            ++mOutside;
        }
        int nfound = 0;
        for (auto it = supp[u].begin(); it != supp[u].end(); ++it) {
            if (nfound == supUV)
                break;
            int w = it->first;
            if (w == v) continue;
            if (supp[v].find(w) != supp[v].end()) {
                ++nfound;
                updateEdge(u, w, supUV, supp);
                updateEdge(v, w, supUV, supp);
            }
        }
        removeEdge(u, v, supp);
    }
}

void fillInNodesEdges() {
    VI nodeTrussness;
    nodeTrussness.resize(n, 0);
    TkEdges.reserve(mTk);
    Tkm1Edges.reserve(mTkm1);
    Skm1Edges.reserve(mTkm1 - mTk);
    outsideEdges.reserve(mOutside);
    v2ShNbrs.clear();
    v2ShNbrs.reserve(nTkm1);
    for (auto &[euv, tuv]: edge2trussness) {
        int eu = euv.u;
        int ev = euv.v;
        if (tuv >= kInput) {  // edges in Tk
            TkEdges.emplace_back(euv);
            Tkm1Edges.emplace_back(euv);
            if (nodeTrussness[eu] < kInput) nodeTrussness[eu] = kInput;
            if (nodeTrussness[ev] < kInput) nodeTrussness[ev] = kInput;
        } else if (tuv == kInput - 1) {
            Tkm1Edges.emplace_back(euv);
            Skm1Edges.emplace_back(euv);
            v2ShNbrs[eu].emplace_back(ev);
            v2ShNbrs[ev].emplace_back(eu);
            if (nodeTrussness[eu] < tuv) nodeTrussness[eu] = tuv;
            if (nodeTrussness[ev] < tuv) nodeTrussness[ev] = tuv;
        } else {
            outsideEdges.emplace_back(euv);
            if (nodeTrussness[eu] < tuv) nodeTrussness[eu] = tuv;
            if (nodeTrussness[ev] < tuv) nodeTrussness[ev] = tuv;
        }
    }
    TkNodes.reserve(2 * mTk);
    Tkm1Nodes.reserve(2 * mTkm1);
    Skm1Nodes.reserve(2 * mTkm1);
    outsideNodes.reserve(nTotal);
    for (int i = 0; i < n; ++i) {
        int tvi = nodeTrussness[i];
        if (tvi >= kInput) {
            TkNodes.emplace_back(i);
            Tkm1Nodes.emplace_back(i);
        } else if (tvi == kInput - 1) {
            Tkm1Nodes.emplace_back(i);
            Skm1Nodes.emplace_back(i);
        } else {
            outsideNodes.emplace_back(i);
        }
    }
    nTk = (int) TkNodes.size();
    nTkm1 = (int) Tkm1Nodes.size();
    nOutside = (int) outsideNodes.size();
}

inline void printNodeEdgeInfo() {
    cout << "#edges outside " << mOutside << endl
         << "#edges in (k-1)-truss-shell " << mTkm1 - mTk << endl
         << "#edges in k-truss " << mTk << endl
         << "#nodes outside " << nOutside << endl
         << "#nodes in (k-1)-truss-shell " << nTkm1 - nTk << endl
         << "#nodes in k-truss " << nTk << endl;
}

inline void printNodeEdgeInfoFile() {
    fout << "#edges outside " << mOutside << endl
         << "#edges in (k-1)-truss-shell " << mTkm1 - mTk << endl
         << "#edges in k-truss " << mTk << endl
         << "#nodes outside " << nOutside << endl
         << "#nodes in (k-1)-truss-shell " << nTkm1 - nTk << endl
         << "#nodes in k-truss " << nTk << endl;
}

void computeNbrsInTkm1() {
    if (debugMode_) cout << "computing the neighbors in the (k-1)-truss" << endl;
    for (auto in_ = 0; in_ < n; ++in_) nbrsInTkm1[in_].clear();
    nodesWithTkm1Nbrs.clear();
    for (auto vp: Tkm1Nodes) {
        for (auto &it: suppOrig[vp]) {
            int w = it.first;
            nbrsInTkm1[w].emplace_back(vp);
            nodesWithTkm1Nbrs.insert(w);
        }
    }
    if (debugMode_) cout << nodesWithTkm1Nbrs.size() << " nodes have nbrs inside Tkm1" << endl;
    if (debugMode_) cout << endl;
}

void computeNbrsInTk() {
    if (debugMode_) cout << "computing the neighbors in the k-truss" << endl;
    for (auto in_ = 0; in_ < n; ++in_) nbrsInTk[in_].clear();
    for (auto vp: TkNodes) {
        for (auto &it: suppOrig[vp]) {
            int w = it.first;
            nbrsInTk[w].emplace_back(vp);
        }
    }
    if (debugMode_) cout << endl;
}

void findUnstableEdges() {
    if (debugMode_) cout << "construct the set of unstable edges" << endl;
    unstableEdges.clear();
    unstableEdges.reserve(mTkm1 - mTk);
    v2UsNbrs.clear();
    v2UsNbrs.reserve(nTkm1);
    for (auto &euv: Skm1Edges) {
        int eu = euv.u;
        int ev = euv.v;
        if (suppTkm1[eu][ev] <= kInput - 3) {
            unstableEdges.emplace_back(euv);
            v2UsNbrs[eu].emplace_back(ev);
            v2UsNbrs[ev].emplace_back(eu);
        }
    }
    if (debugMode_) cout << "there are " << unstableEdges.size() << " unstable edges" << endl;
    if (debugMode_) cout << endl;
}

void findBestOutsideNodes() {
    if (debugMode_) cout << "maximal-set-based pruning of outside-nodes" << endl;
    if (debugMode_) cout << nOutside << " outside nodes in total" << endl;
    unordered_set<VI, VectorHasher> seen;
    seen.reserve(outsideNodes.size());
    uniqueOutsideNodes.clear();
    uniqueOutsideNodes.reserve(nOutside);
    for (auto &v: outsideNodes) {
        if (!nodesWithTkm1Nbrs.contains(v)) continue;
        VI nbrV = nbrsInTkm1[v];
        if (!seen.contains(nbrV)) {
            seen.insert(nbrV);
            uniqueOutsideNodes.emplace_back(v);
        }
    }
    if (debugMode_) cout << "there are " << uniqueOutsideNodes.size() << " unique-nbr outside-nodes" << endl;

    if (uniqueOutsideNodes.size() > maxN) {
        bestOutsideNodes = uniqueOutsideNodes;
    } else {
        unordered_map<int, bitset<maxN>> setsContainingElement;
        setsContainingElement.reserve(nTkm1);
        int index = 0;
        for (auto &vo: uniqueOutsideNodes) {
            for (auto &no: nbrsInTkm1[vo]) {
                setsContainingElement[no].set(index, true);
            }
            ++index;
        }

        bestOutsideNodes.clear();
        bestOutsideNodes.reserve(uniqueOutsideNodes.size());
        for (auto &vo: uniqueOutsideNodes) {
            VI &nbrVo = nbrsInTkm1[vo];
            if (nbrVo.empty()) continue;
            bitset<maxN> result;
            bool flag = true;
            for (auto &no: nbrVo) {
                if (flag) {
                    result = setsContainingElement[no];
                    flag = false;
                } else {
                    result &= setsContainingElement[no];
                }
            }
            if (result.count() == 1) {
                bestOutsideNodes.emplace_back(vo);
            }
        }
    }
    if (debugMode_) cout << "there are " << bestOutsideNodes.size() << " best outside-nodes" << endl;
    if (debugMode_) cout << endl;
}

void computeHUEs() {
    if (debugMode_) cout << "computing HUEs" << endl;
    e2HUEs.clear();
    v2HUEs.clear();
    for (auto &euv: unstableEdges) {
        int eu = euv.u, ev = euv.v;
        VI numv, nvmu;
        unordered_set<int> nuSet;
        nuSet.reserve(suppTkm1[eu].size());
        for (auto &[n_, _]: suppTkm1[eu]) {
            if (n_ == ev) continue;
            nuSet.insert(n_);
        }
        auto &suppTkm1Ev = suppTkm1[ev];
        nvmu.reserve(min(nuSet.size(), suppTkm1Ev.size()));
        for (auto &[x_, _]: suppTkm1Ev) {
            if (x_ == eu) continue;
            if (nuSet.contains(x_)) {
                nuSet.erase(x_);
            } else {  // in Nv not in Nu
                nvmu.emplace_back(x_);
            }
        }
        numv.reserve(nuSet.size());
        for (auto &x_: nuSet) {
            numv.emplace_back(x_);
        }
        if (!numv.empty()) ++v2HUEs[ev];
        for (auto &w: numv) {
            e2HUEs[TEdge{min(ev, w), max(ev, w)}].emplace_back(euv);
            ++v2HUEs[w];
        }
        if (!nvmu.empty()) ++v2HUEs[eu];
        for (auto &w: nvmu) {
            e2HUEs[TEdge{min(eu, w), max(eu, w)}].emplace_back(euv);
            ++v2HUEs[w];
        }
    }
    if (debugMode_) cout << "there are " << e2HUEs.size() << " edges that have HUEs" << endl;
    if (debugMode_) cout << "there are " << v2HUEs.size() << " nodes that have HUEs" << endl;
    if (debugMode_) cout << endl;
}

void computeHNNs() {
    if (debugMode_) cout << "computing HNNs" << endl;
    v2HNNs.clear();
    v2HNNs.reserve(nTkm1);
    for (auto &[e_, HUEe_]: e2HUEs) {
        int x_ = e_.u;
        int y_ = e_.v;
        v2HNNs[x_].emplace_back(y_);
        v2HNNs[y_].emplace_back(x_);
    }
    if (debugMode_) cout << endl;
}

inline auto compMoreHUEs(int i_, int j_) {
    auto huesI = v2HUEs[i_];
    auto huesJ = v2HUEs[j_];
    return huesI > huesJ || (huesI == huesJ && i_ < j_);
}

inline auto compMoreNbrsInTkm1(int i_, int j_) {
    auto nbrsI = nbrsInTkm1[i_].size();
    auto nbrsJ = nbrsInTkm1[j_].size();
    return nbrsI > nbrsJ || (nbrsI == nbrsJ && i_ < j_);
}

inline auto compMoreIncidentPotential(int i_, int j_) {
    auto nbrsI = nbrsInTkm1[i_].size();
    auto nbrsJ = nbrsInTkm1[j_].size();
    auto mkI = !suppTk.empty() ? suppTk[i_].size() : 0;
    auto mkJ = !suppTk.empty() ? suppTk[j_].size() : 0;
    auto mipI = nbrsI - mkI;
    auto mipJ = nbrsJ - mkJ;
    if (mipI != mipJ) return mipI > mipJ;
    if (mkI != mkJ) return mkI > mkJ;
    return i_ < j_;
}

void updateGraph(int v1_, int v2_) {
    if (debugMode_) cout << "updateGraph: merging " << v1_ << " and " << v2_ << endl;
    mTk = mTkm1 = mOutside = nTk = nTkm1 = nOutside = 0;
    int mDiff;

    vector<bool> n1_(n, false), n2_(n, false);
    auto nSize1 = suppOrig[v1_].size();
    auto nSize2 = suppOrig[v2_].size();
    VI cn12_, n1m2, n2m1;
    cn12_.reserve(min(nSize1, nSize2));
    n1m2.reserve(nSize1);
    n2m1.reserve(nSize2);
    for (auto &it: suppOrig[v1_]) {
        auto &x_ = it.first;
        if (x_ != v2_) n1_[x_] = true;
    }
    for (auto &it: suppOrig[v2_]) {
        auto &x_ = it.first;
        if (x_ != v1_) {
            if (n1_[x_]) {
                cn12_.emplace_back(x_);
            } else {
                n2m1.emplace_back(x_);
            }
            n2_[x_] = true;
        }
    }
    for (auto &it: suppOrig[v1_]) {
        auto &x_ = it.first;
        if (!n2_[x_]) n1m2.emplace_back(x_);
    }

    if (debugMode_) cout << "updateGraph: updating m" << endl;
    if (suppOrig[v1_].contains(v2_)) {
        TEdge v1v2 = {min(v1_, v2_), max(v1_, v2_)};
        edge2trussness.erase(v1v2);
        mDiff = suppOrig[v1_][v2_] + 1;
    } else {
        mDiff = (int) cn12_.size();
    }
    mTotal -= mDiff;
    --nTotal;
    if (debugMode_)
        cout << "updateGraph: #nodes and edges after the merger = ("
             << nTotal << ", " << mTotal << ")" << endl;

    if (debugMode_) cout << "updateGraph: updating edges, degrees, supports, etc." << endl;
    deg[v1_] += deg[v2_] - mDiff;
    deg[v2_] = 0;
    for (auto v_: cn12_) {
        --deg[v_];
    }
    for (auto v_: n2m1) {
        if (v_ == v1_) continue;
        TEdge newE_ = {min(v1_, v_), max(v1_, v_)};
        edge2trussness[newE_] = kInput;
        suppOrig[v1_][v_] = 0;
        suppOrig[v_][v1_] = 0;
    }

    for (auto &it: suppOrig[v2_]) {
        auto u_ = it.first;
        suppOrig[u_].erase(v2_);
        TEdge e2u = {min(v2_, u_), max(v2_, u_)};
        edge2trussness.erase(e2u);
    }
    suppOrig[v2_].clear();

    if (debugMode_) cout << "updateGraph: for each edge among the common neighbors, support decreases by 1" << endl;
    if (cn12_.size() > 1) {
        for (auto i_ = 0; i_ < cn12_.size() - 1; ++i_) {
            int cni = cn12_[i_];
            for (auto j_ = i_ + 1; j_ < cn12_.size(); ++j_) {
                int cnj = cn12_[j_];
                if (suppOrig[cni].contains(cnj)) {
                    --suppOrig[cni][cnj];
                    --suppOrig[cnj][cni];
                }
            }
        }
    }

    if (debugMode_) cout << "updateGraph: for each (n1m2, n2m1)-type edge, support increases by 1" << endl;
    for (auto na_: n1m2) {
        for (auto nb_: n2m1) {
            if (suppOrig[na_].contains(nb_)) {
                ++suppOrig[na_][nb_];
                ++suppOrig[nb_][na_];
            }
        }
    }

    if (debugMode_) cout << "updateGraph: computing the support of incident edges" << endl;
    for (auto &it: suppOrig[v1_]) {
        auto u_ = it.first;
        int supp12_ = 0;
        for (auto &it_: suppOrig[u_]) {
            if (suppOrig[v1_].contains(it_.first)) ++supp12_;
        }
        suppOrig[v1_][u_] = supp12_;
        suppOrig[u_][v1_] = supp12_;
    }
    if (debugMode_) cout << endl;
}

void varsReInit() {
    outsideEdges.clear();
    Tkm1Edges.clear();
    TkEdges.clear();
    Skm1Edges.clear();
    outsideNodes.clear();
    Tkm1Nodes.clear();
    TkNodes.clear();
    Skm1Nodes.clear();
    for (auto &[e_, te_]: edge2trussness) {
        edge2trussness[e_] = kInput;
    }
}

int checkMergerResultIOP(int v1_, int v2_) {
    auto suppOrigBackup = suppTkm1;
    if (debugMode_) cout << "checkMergerResult: checking " << v1_ << " and " << v2_ << endl;
    // new neighbors
    if (debugMode_) cout << "checkMergerResult: computing new neighbors" << endl;
    VI n1NewNbrs;
    n1NewNbrs.reserve(suppOrig[v1_].size() + suppOrig[v2_].size());
    for (auto &it: suppOrig[v1_]) {
        auto v1n_ = it.first;
        if (v1n_ == v2_) continue;
        if (!suppOrigBackup[v1_].contains(v1n_)) {
            n1NewNbrs.emplace_back(v1n_);
            suppOrigBackup[v1_][v1n_] = 0;
            suppOrigBackup[v1n_][v1_] = 0;
        }
    }
    for (auto &it: suppOrig[v2_]) {
        auto v2n_ = it.first;
        if (v2n_ == v1_) continue;
        if (!suppOrig[v1_].contains(v2n_)) {
            n1NewNbrs.emplace_back(v2n_);
            suppOrigBackup[v1_][v2n_] = 0;
            suppOrigBackup[v2n_][v1_] = 0;
        }
    }
    int numNewNbrs = (int) n1NewNbrs.size();
    int mTk_ = mTkm1 + numNewNbrs;  // potential k-truss-edges

    // new neighbors bring new supports
    if (debugMode_) cout << "checkMergerResult: updating the supports of non-incident edges" << endl;
    if (numNewNbrs > 2) {
        for (auto i_ = 0; i_ < numNewNbrs - 1; i_++) {
            auto nn1_ = n1NewNbrs[i_];
            for (auto j_ = i_ + 1; j_ < numNewNbrs; j_++) {
                auto nn2_ = n1NewNbrs[j_];
                if (suppOrigBackup[nn1_].contains(nn2_)) {
                    updateSupport(nn1_, nn2_, 1, suppOrigBackup);
                }
            }
        }
    }
    for (auto nn_: n1NewNbrs) {
        for (auto &it: suppTkm1[v1_]) {
            auto on_ = it.first;
            if (suppOrigBackup[nn_].contains(on_)) {
                updateSupport(nn_, on_, 1, suppOrigBackup);
            }
        }
    }

    if (debugMode_) cout << "checkMergerResult: computing the support of incident edges" << endl;
    for (auto &it: suppOrigBackup[v1_]) {
        auto u_ = it.first;
        int supp12_ = 0;
        for (auto &it_: suppOrigBackup[u_]) {
            if (suppOrigBackup[v1_].contains(it_.first)) ++supp12_;
        }
        suppOrigBackup[v1_][u_] = supp12_;
        suppOrigBackup[u_][v1_] = supp12_;
    }

    if (debugMode_) cout << "checkMergerResult: bin-sorting" << endl;
    bin.clear();
    bin.resize(n, 0);
    int nBins = 0;
    int mp = 0;
    for (int u = 0; u < n; ++u) {
        auto tsupp = suppOrigBackup[u];
        for (auto &it: tsupp) {
            int v = it.first;
            if (!compVertex(u, v)) continue;
            int sup = it.second;
            if (sup == 0) {
                suppOrigBackup[u].erase(v);
                suppOrigBackup[v].erase(u);
                --mTk_;
                continue;
            }
            ++mp;
            ++bin[sup];
            nBins = max(sup, nBins);
        }
    }
    m = mp;
    ++nBins;
    int count = 0;
    for (int i = 0; i < nBins; ++i) {
        int binsize = bin[i];
        bin[i] = count;
        count += binsize;
    }
    pos.clear();
    pos.resize(n);
    for (int i = 0; i < n; ++i) pos[i].clear();
    binEdge.resize(m);
    for (int u = 0; u < n; ++u)
        for (auto it = suppOrigBackup[u].begin(); it != suppOrigBackup[u].end(); ++it) {
            int v = it->first;
            if (!compVertex(u, v)) continue;
            int sup = it->second;
            TEdge e = {u, v};
            int &b = bin[sup];
            binEdge[b] = e;
            pos[u][v] = b++;
        }
    for (int i = nBins; i > 0; --i) bin[i] = bin[i - 1];
    bin[0] = 0;

    if (debugMode_) cout << "checkMergerResult: truss decomposition" << endl;
    for (int s = 0; s < m; ++s) {
        int u = binEdge[s].u;
        int v = binEdge[s].v;
        orderPair(u, v);
        int supuv = suppOrigBackup[u][v];
        int tuv = supuv + 2;
        if (tuv >= kInput) {
            if (debugMode_) cout << endl;
            return mTk_;
        }
        --mTk_;
        int nfound = 0;
        for (auto it = suppOrigBackup[u].begin(); it != suppOrigBackup[u].end(); ++it) {
            if (nfound == supuv)
                break;
            int w = it->first;
            if (w == v) continue;
            if (suppOrigBackup[v].contains(w)) {
                ++nfound;
                if (!compVertex(u, w)) {
                    updateEdge(w, u, supuv, suppOrigBackup);
                } else {
                    updateEdge(u, w, supuv, suppOrigBackup);
                }

                if (!compVertex(v, w)) {
                    updateEdge(w, v, supuv, suppOrigBackup);
                } else {
                    updateEdge(v, w, supuv, suppOrigBackup);
                }
            }
        }
        removeEdge(u, v, suppOrigBackup);
    }
    if (debugMode_) cout << endl;
    return -1;
}

int checkMergerResultIOPInputNewNbrs(int v1_, int v2_, const VI &newNbrs_) {
    if (debugMode_) cout << "checkMergerResult: checking " << v1_ << " and " << v2_ << endl;
    auto suppOrigBackup = suppTkm1;
    // new neighbors
    if (debugMode_) cout << "checkMergerResult: computing new neighbors" << endl;
    vector<bool> n1NewNbrs(n, false);
    int numNewNbrs = 0;
    for (auto nn_: newNbrs_) {
        if (nn_ == v1_) continue;
        ++numNewNbrs;
        n1NewNbrs[nn_] = true;
        suppOrigBackup[v1_][nn_] = 0;
        suppOrigBackup[nn_][v1_] = 0;
    }
    int mTk_ = mTkm1 + numNewNbrs;  // potential k-truss-edges

    // new neighbors bring new supports
    if (debugMode_) cout << "checkMergerResult: updating the supports of non-incident edges" << endl;

    for (auto &[x_, y_]: Skm1Edges) {
        if (x_ == v1_ || x_ == v2_ || y_ == v1_ || y_ == v2_) continue;
        if ((n1NewNbrs[x_] && n1NewNbrs[y_]) ||
            (n1NewNbrs[x_] && suppTkm1[v1_].contains(y_)) ||
            (n1NewNbrs[y_] && suppTkm1[v1_].contains(x_))) {
            updateSupport(x_, y_, 1, suppOrigBackup);
        }
    }

    if (debugMode_) cout << "checkMergerResult: computing the support of incident edges" << endl;
    for (auto &it: suppOrigBackup[v1_]) {
        auto u_ = it.first;
        int supp12_ = 0;
        for (auto &it_: suppOrigBackup[u_]) {
            if (suppOrigBackup[v1_].contains(it_.first)) ++supp12_;
        }
        suppOrigBackup[v1_][u_] = supp12_;
        suppOrigBackup[u_][v1_] = supp12_;
    }

    if (debugMode_) cout << "checkMergerResult: bin-sorting" << endl;
    bin.clear();
    bin.resize(n, 0);
    int nBins = 0;
    int mp = 0;
    for (int u = 0; u < n; ++u) {
        auto tsupp = suppOrigBackup[u];
        for (auto &it: tsupp) {
            int v = it.first;
            if (!compVertex(u, v)) continue;
            int sup = it.second;
            if (sup == 0) {
                suppOrigBackup[u].erase(v);
                suppOrigBackup[v].erase(u);
                --mTk_;
                continue;
            }
            ++mp;
            ++bin[sup];
            nBins = max(sup, nBins);
        }
    }
    m = mp;
    ++nBins;
    int count = 0;
    for (int i = 0; i < nBins; ++i) {
        int binsize = bin[i];
        bin[i] = count;
        count += binsize;
    }
    pos.clear();
    pos.resize(n);
    for (int i = 0; i < n; ++i) pos[i].clear();
    binEdge.resize(m);
    for (int u = 0; u < n; ++u) {
        for (auto it = suppOrigBackup[u].begin(); it != suppOrigBackup[u].end(); ++it) {
            int v = it->first;
            if (!compVertex(u, v)) continue;
            int sup = it->second;
            TEdge e = {u, v};
            int &b = bin[sup];
            binEdge[b] = e;
            pos[u][v] = b++;
        }
    }
    for (int i = nBins; i > 0; --i) bin[i] = bin[i - 1];
    bin[0] = 0;

    if (debugMode_) cout << "checkMergerResult: truss decomposition" << endl;
    for (int s = 0; s < m; ++s) {
        int u = binEdge[s].u;
        int v = binEdge[s].v;
        orderPair(u, v);
        int supuv = suppOrigBackup[u][v];
        int tuv = supuv + 2;
        if (tuv >= kInput) {
            if (debugMode_) cout << endl;
            return mTk_;
        }
        --mTk_;
        int nfound = 0;
        for (auto it = suppOrigBackup[u].begin(); it != suppOrigBackup[u].end(); ++it) {
            if (nfound == supuv)
                break;
            int w = it->first;
            if (w == v) continue;
            if (suppOrigBackup[v].contains(w)) {
                ++nfound;
                if (!compVertex(u, w)) {
                    updateEdge(w, u, supuv, suppOrigBackup);
                } else {
                    updateEdge(u, w, supuv, suppOrigBackup);
                }
                if (!compVertex(v, w)) {
                    updateEdge(w, v, supuv, suppOrigBackup);
                } else {
                    updateEdge(v, w, supuv, suppOrigBackup);
                }
            }
        }
        removeEdge(u, v, suppOrigBackup);
    }
    if (debugMode_) cout << endl;
    return -1;
}

int checkMergerResultGeneral(int v1_, int v2_) {
    // general cases: v1_ and v2_ can be anywhere
    auto suppOrigBackup = suppTkm1;
    bool v1Outside = suppOrigBackup[v1_].empty();
    bool v2Outside = suppOrigBackup[v2_].empty();
    if (debugMode_) cout << "checkMergerResult: checking " << v1_ << " and " << v2_ << endl;

    // new neighbors
    if (debugMode_) cout << "checkMergerResult: computing new neighbors" << endl;
    VI n1NewNbrs;
    n1NewNbrs.reserve(nbrsInTkm1[v1_].size() + nbrsInTkm1[v2_].size());

    for (auto &v1n_: nbrsInTkm1[v1_]) {
        if (v1n_ == v2_) continue;
        if (v1Outside || !suppOrigBackup[v1_].contains(v1n_)) {
            n1NewNbrs.emplace_back(v1n_);
            suppOrigBackup[v1_][v1n_] = 0;
            suppOrigBackup[v1n_][v1_] = 0;
        }
    }
    for (auto &v2n_: nbrsInTkm1[v2_]) {
        if (v2n_ == v1_) continue;
        if (!suppOrigBackup[v1_].contains(v2n_)) {
            n1NewNbrs.emplace_back(v2n_);
            suppOrigBackup[v1_][v2n_] = 0;
            suppOrigBackup[v2n_][v1_] = 0;
        }
    }

    int numNewNbrs = (int) n1NewNbrs.size();

    // new neighbors bring new supports
    if (debugMode_) cout << "checkMergerResult: updating the supports of non-incident edges" << endl;
    if (numNewNbrs > 1) {
        for (auto i_ = 0; i_ < numNewNbrs - 1; i_++) {
            auto nn1_ = n1NewNbrs[i_];
            for (auto j_ = i_ + 1; j_ < numNewNbrs; j_++) {
                auto nn2_ = n1NewNbrs[j_];
                if (suppOrigBackup[nn1_].contains(nn2_)) {
                    updateSupport(nn1_, nn2_, 1, suppOrigBackup);
                }
            }
        }
    }
    for (auto nn_: n1NewNbrs) {
        for (auto &it: suppTkm1[v1_]) {
            auto on_ = it.first;
            if (suppOrigBackup[nn_].contains(on_)) {
                updateSupport(nn_, on_, 1, suppOrigBackup);
            }
        }
    }

    // if v2 is in the current (k-1)-truss, we need to deal with the loss
    int mTk_ = mTkm1 + numNewNbrs;
    if (!v2Outside) {
        VI nbr2_;
        nbr2_.reserve(suppOrigBackup[v2_].size());
        for (auto &[n2_, _]: suppTkm1[v2_]) {
            suppOrigBackup[n2_].erase(v2_);
            if (n2_ != v1_) nbr2_.emplace_back(n2_);
            --mTk_;
        }
        suppOrigBackup[v2_].clear();
        if (nbr2_.size() > 1) {
            for (auto i_ = 0; i_ < nbr2_.size() - 1; ++i_) {
                auto n2i_ = nbr2_[i_];
                for (auto j_ = i_ + 1; j_ < nbr2_.size(); ++j_) {
                    auto n2j_ = nbr2_[j_];
                    if (suppTkm1[n2i_].contains(n2j_)) {
                        updateSupport(n2i_, n2j_, -1, suppOrigBackup);
                    }
                }
            }
        }
    }

    if (debugMode_) cout << "checkMergerResult: computing the support of incident edges" << endl;
    for (auto &it: suppOrigBackup[v1_]) {
        auto u_ = it.first;
        int supp12_ = 0;
        for (auto &it_: suppOrigBackup[u_]) {
            if (suppOrigBackup[v1_].contains(it_.first)) ++supp12_;
        }
        suppOrigBackup[v1_][u_] = supp12_;
        suppOrigBackup[u_][v1_] = supp12_;
    }

    if (debugMode_) cout << "checkMergerResult: bin-sorting" << endl;
    bin.clear();
    bin.resize(n, 0);
    int nBins = 0;
    int mp = 0;
    for (int u = 0; u < n; ++u) {
        auto tsupp = suppOrigBackup[u];
        for (auto &it: tsupp) {
            int v = it.first;
            if (!compVertex(u, v)) continue;
            int sup = it.second;
            if (sup == 0) {
                suppOrigBackup[u].erase(v);
                suppOrigBackup[v].erase(u);
                --mTk_;
                continue;
            }
            ++mp;
            ++bin[sup];
            nBins = max(sup, nBins);
        }
    }
    m = mp;
    ++nBins;
    int count = 0;
    for (int i = 0; i < nBins; ++i) {
        int binsize = bin[i];
        bin[i] = count;
        count += binsize;
    }
    pos.clear();
    pos.resize(n);
    for (int i = 0; i < n; ++i) pos[i].clear();
    binEdge.resize(m);
    for (int u = 0; u < n; ++u) {
        for (auto it = suppOrigBackup[u].begin(); it != suppOrigBackup[u].end(); ++it) {
            int v = it->first;
            if (!compVertex(u, v)) continue;
            int sup = it->second;
            TEdge e = {u, v};
            int &b = bin[sup];
            binEdge[b] = e;
            pos[u][v] = b++;
        }
    }
    for (int i = nBins; i > 0; --i) bin[i] = bin[i - 1];
    bin[0] = 0;

    if (debugMode_) cout << "checkMergerResult: truss decomposition" << endl;
    for (int s = 0; s < m; ++s) {
        int u = binEdge[s].u;
        int v = binEdge[s].v;
        orderPair(u, v);
        int supuv = suppOrigBackup[u][v];
        int tuv = supuv + 2;
        if (tuv >= kInput) {
            if (debugMode_) cout << endl;
            return mTk_;
        }
        --mTk_;
        int nfound = 0;
        for (auto it = suppOrigBackup[u].begin(); it != suppOrigBackup[u].end(); ++it) {
            if (nfound == supuv)
                break;
            int w = it->first;
            if (w == v) continue;
            if (suppOrigBackup[v].contains(w)) {
                ++nfound;
                updateEdge(w, u, supuv, suppOrigBackup);
                updateEdge(w, v, supuv, suppOrigBackup);
            }
        }
        removeEdge(u, v, suppOrigBackup);
    }
    if (debugMode_) cout << endl;
    return -1;
}

int checkMergerResultIIP(int v1_, int v2_) {
    // both v1_ and v2_ are inside-nodes
    auto suppOrigBackup = suppTkm1;
    if (debugMode_) cout << "checkMergerResult: checking " << v1_ << " and " << v2_ << endl;

    // new neighbors
    if (debugMode_) cout << "checkMergerResult: computing new neighbors" << endl;
    auto &n1_ = nbrsInTkm1[v1_];
    auto &n2_ = nbrsInTkm1[v2_];
    vector<bool> n1NewNbrs(n, false);
    int numNewNbrs = 0;
    for (auto &v1n_: n1_) {
        if (v1n_ == v2_) continue;
        if (!suppOrigBackup[v1_].contains(v1n_)) {
            n1NewNbrs[v1n_] = true;
            ++numNewNbrs;
            suppOrigBackup[v1_][v1n_] = 0;
            suppOrigBackup[v1n_][v1_] = 0;
        }
    }
    for (auto &v2n_: n2_) {
        if (v2n_ == v1_) continue;
        if (!suppOrigBackup[v1_].contains(v2n_)) {
            n1NewNbrs[v2n_] = true;
            ++numNewNbrs;
            suppOrigBackup[v1_][v2n_] = 0;
            suppOrigBackup[v2n_][v1_] = 0;
        }
    }

    // new neighbors bring new supports
    if (debugMode_) cout << "checkMergerResult: updating the supports of non-incident edges" << endl;

    // deal with the loss
    int mTk_ = mTkm1 + numNewNbrs;
    vector<bool> nbr2_(n, false);
    unsigned int nnbr2_ = 0;
    for (auto &[nn2_, _]: suppTkm1[v2_]) {
        suppOrigBackup[nn2_].erase(v2_);
        if (nn2_ != v1_) {
            nbr2_[nn2_] = true;
            ++nnbr2_;
        }
        --mTk_;
    }
    suppOrigBackup[v2_].clear();

    for (auto &[x_, y_]: Tkm1Edges) {
        if (x_ == v1_ || x_ == v2_ || y_ == v1_ || y_ == v2_) continue;
        if ((n1NewNbrs[x_] && n1NewNbrs[y_]) ||
            (n1NewNbrs[x_] && suppTkm1[v1_].contains(y_)) ||
            (n1NewNbrs[y_] && suppTkm1[v1_].contains(x_))) {
            updateSupport(x_, y_, 1, suppOrigBackup);
        }
        if (nbr2_[x_] && nbr2_[y_]) updateSupport(x_, y_, -1, suppOrigBackup);
    }

    if (debugMode_) cout << "checkMergerResult: computing the support of incident edges" << endl;
    for (auto &[u_, _]: suppOrigBackup[v1_]) {
        int supp12_ = 0;
        for (auto &[w_, xx_]: suppOrigBackup[u_]) {
            if (suppOrigBackup[v1_].contains(w_)) ++supp12_;
        }
        suppOrigBackup[v1_][u_] = supp12_;
        suppOrigBackup[u_][v1_] = supp12_;
    }

    if (debugMode_) cout << "checkMergerResult: bin-sorting" << endl;
    bin.clear();
    bin.resize(n, 0);
    int nBins = 0;
    int mp = 0;
    for (int u = 0; u < n; ++u) {
        auto tsupp = suppOrigBackup[u];
        for (auto &it: tsupp) {
            int v = it.first;
            if (!compVertex(u, v)) continue;
            int sup = it.second;
            if (sup == 0) {
                suppOrigBackup[u].erase(v);
                suppOrigBackup[v].erase(u);
                --mTk_;
                continue;
            }
            ++mp;
            ++bin[sup];
            nBins = max(sup, nBins);
        }
    }
    m = mp;
    ++nBins;
    int count = 0;
    for (int i = 0; i < nBins; ++i) {
        int binsize = bin[i];
        bin[i] = count;
        count += binsize;
    }
    pos.clear();
    pos.resize(n);
    for (int i = 0; i < n; ++i) pos[i].clear();
    binEdge.resize(m);
    for (int u = 0; u < n; ++u) {
        for (auto it = suppOrigBackup[u].begin(); it != suppOrigBackup[u].end(); ++it) {
            int v = it->first;
            if (!compVertex(u, v)) continue;
            int sup = it->second;
            TEdge e = {u, v};
            int &b = bin[sup];
            binEdge[b] = e;
            pos[u][v] = b++;
        }
    }
    for (int i = nBins; i > 0; --i) bin[i] = bin[i - 1];
    bin[0] = 0;

    if (debugMode_) cout << "checkMergerResult: truss decomposition" << endl;
    for (int s = 0; s < m; ++s) {
        int u = binEdge[s].u;
        int v = binEdge[s].v;
        orderPair(u, v);
        int supuv = suppOrigBackup[u][v];
        int tuv = supuv + 2;
        if (tuv >= kInput) {
            if (debugMode_) cout << endl;
            return mTk_;
        }
        --mTk_;
        int nfound = 0;
        for (auto it = suppOrigBackup[u].begin(); it != suppOrigBackup[u].end(); ++it) {
            if (nfound == supuv)
                break;
            int w = it->first;
            if (w == v) continue;
            if (suppOrigBackup[v].contains(w)) {
                ++nfound;
                updateEdge(w, u, supuv, suppOrigBackup);
                updateEdge(w, v, supuv, suppOrigBackup);
            }
        }
        removeEdge(u, v, suppOrigBackup);
    }
    if (debugMode_) cout << endl;
    return -1;
}
