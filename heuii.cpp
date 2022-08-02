#include "truss.h"

using namespace std;

int main(int argc, char *argv[]) {
    string progName(argv[0]);
    progName = progName.substr(2);
    //    cout << progName << endl;

    assert(argc >= 3);
    string dataset(argv[1]);
    string resInfo(argv[2]);
    int argi = 4;
    kInput = (argc >= argi) ? atoi(argv[argi - 1]) : 10;

    infile = "datasets_txt/" + dataset + ".txt";
    outfile = "res_txt/" + resInfo + ".txt";
    fout.open(outfile.c_str(), ios_base::app);
    fout << progName << " ";
    fout << dataset << " ";
    fout << kInput << endl;

    readOrderedSimpleGraph();
    countTriangles(supp);
    suppOrig = supp;  // supp backup
    binSort(supp);
    trussDecomp();
    startTime = chrono::system_clock::now();
    fillInNodesEdges();
    printNodeEdgeInfoFile();

    computeNbrsInTkm1();
    computeNbrsInTk();
    findUnstableEdges();
    findBestOutsideNodes();
    computeHUEs();

    VI &v1Pruned = Tkm1Nodes;
    int maxV1 = 100;
    maxV1 = min((int)v1Pruned.size(), maxV1);
    SE bestPairsTotal;
    bestPairsTotal.reserve(maxV1 * (maxV1 - 1) * 3);

    sort(v1Pruned.begin(), v1Pruned.end(), compMoreNbrsInTkm1);
    fout << "bestV1Nbrs" << endl;
    for (auto i_ = 0; i_ < maxV1; ++i_) {
        auto vi = v1Pruned[i_];
        fout << vi << " ";
    }
    fout << endl;
    for (auto i_ = 0; i_ < maxV1 - 1; ++i_) {
        auto vi = v1Pruned[i_];
        for (auto j_ = i_ + 1; j_ < maxV1; ++j_) {
            auto vj = v1Pruned[j_];
            bestPairsTotal.insert(TEdge{min(vi, vj), max(vi, vj)});
        }
    }

    sort(v1Pruned.begin(), v1Pruned.end(), compMoreHUEs);
    fout << "bestV1HUEs" << endl;
    for (auto i_ = 0; i_ < maxV1; ++i_) {
        auto vi = v1Pruned[i_];
        fout << vi << " ";
    }
    fout << endl;
    for (auto i_ = 0; i_ < maxV1 - 1; ++i_) {
        auto vi = v1Pruned[i_];
        for (auto j_ = i_ + 1; j_ < maxV1; ++j_) {
            auto vj = v1Pruned[j_];
            bestPairsTotal.insert(TEdge{min(vi, vj), max(vi, vj)});
        }
    }

    sort(v1Pruned.begin(), v1Pruned.end(), compMoreIncidentPotential);
    fout << "bestV1IPs" << endl;
    for (auto i_ = 0; i_ < maxV1; ++i_) {
        auto vi = v1Pruned[i_];
        fout << vi << " ";
    }
    fout << endl;
    for (auto i_ = 0; i_ < maxV1 - 1; ++i_) {
        auto vi = v1Pruned[i_];
        for (auto j_ = i_ + 1; j_ < maxV1; ++j_) {
            auto vj = v1Pruned[j_];
            bestPairsTotal.insert(TEdge{min(vi, vj), max(vi, vj)});
        }
    }

    gen.seed(0);
    shuffle(v1Pruned.begin(), v1Pruned.end(), gen);
    fout << "random 0" << endl;
    for (auto i_ = 0; i_ < maxV1; ++i_) {
        auto vi = v1Pruned[i_];
        fout << vi << " ";
    }
    fout << endl;
    for (auto i_ = 0; i_ < maxV1 - 1; ++i_) {
        auto vi = v1Pruned[i_];
        for (auto j_ = i_ + 1; j_ < maxV1; ++j_) {
            auto vj = v1Pruned[j_];
            bestPairsTotal.insert(TEdge{min(vi, vj), max(vi, vj)});
        }
    }

    gen.seed(1);
    shuffle(v1Pruned.begin(), v1Pruned.end(), gen);
    fout << "random 1" << endl;
    for (auto i_ = 0; i_ < maxV1; ++i_) {
        auto vi = v1Pruned[i_];
        fout << vi << " ";
    }
    fout << endl;
    for (auto i_ = 0; i_ < maxV1 - 1; ++i_) {
        auto vi = v1Pruned[i_];
        for (auto j_ = i_ + 1; j_ < maxV1; ++j_) {
            auto vj = v1Pruned[j_];
            bestPairsTotal.insert(TEdge{min(vi, vj), max(vi, vj)});
        }
    }

    gen.seed(2);
    shuffle(v1Pruned.begin(), v1Pruned.end(), gen);
    fout << "random 2" << endl;
    for (auto i_ = 0; i_ < maxV1; ++i_) {
        auto vi = v1Pruned[i_];
        fout << vi << " ";
    }
    fout << endl;
    for (auto i_ = 0; i_ < maxV1 - 1; ++i_) {
        auto vi = v1Pruned[i_];
        for (auto j_ = i_ + 1; j_ < maxV1; ++j_) {
            auto vj = v1Pruned[j_];
            bestPairsTotal.insert(TEdge{min(vi, vj), max(vi, vj)});
        }
    }

    tqdm bar;
    int curr = 0;
    int totalPairs = (int)bestPairsTotal.size();
    unordered_map<int, VI> vOne2vTwoList;
    for (auto &[v1_, v2_] : bestPairsTotal) {
        vOne2vTwoList[min(v1_, v2_)].emplace_back(max(v1_, v2_));
    }
    for (auto &[v1_, v2List_] : vOne2vTwoList) {
        auto &n1_ = nbrsInTkm1[v1_];
        vector<bool> membershipCnt1(n, false);
        vector<bool> membershipCntK1(n, false);
        for (auto &x_ : n1_) {
            membershipCnt1[x_] = true;
        }
        for (auto &x_ : nbrsInTk[v1_]) {
            membershipCntK1[x_] = true;
        }
        for (auto v2_ : v2List_) {
            bar.progress(curr++, totalPairs);
            int scoreShell = 0, scoreTkm1 = 0, scoreHUE = 0, commonNbrTkm1 = 0, commonNbrTk = 0;
            auto &n2_ = nbrsInTkm1[v2_];
            vector<bool> membershipCnt2(n, false);
            for (auto &x_ : n2_) {
                membershipCnt2[x_] = true;
                if (x_ == v1_ || membershipCnt1[x_]) ++commonNbrTkm1;
            }
            for (auto &x_ : nbrsInTk[v2_]) {
                if (x_ == v1_ || membershipCntK1[x_]) ++commonNbrTk;
            }

            for (auto &[x_, y_] : Skm1Edges) {
                if (x_ == v1_ || x_ == v2_ || y_ == v1_ || y_ == v2_) continue;
                bool x1_ = membershipCnt1[x_], x2_ = membershipCnt2[x_];
                bool y1_ = membershipCnt1[y_], y2_ = membershipCnt2[y_];
                if ((!x1_ && !x2_) || (!y1_ && !y2_)) continue;
                if ((!x1_ && !y2_) || (!x2_ && !y1_)) {
                    ++scoreShell;
                } else if (x1_ && x2_ && y1_ && y2_) {
                    --scoreShell;
                }
            }

            for (auto &[x_, y_] : Tkm1Edges) {
                if (x_ == v1_ || x_ == v2_ || y_ == v1_ || y_ == v2_) continue;
                bool x1_ = membershipCnt1[x_], x2_ = membershipCnt2[x_];
                bool y1_ = membershipCnt1[y_], y2_ = membershipCnt2[y_];
                if ((!x1_ && !x2_) || (!y1_ && !y2_)) continue;
                if ((!x1_ && !y2_) || (!x2_ && !y1_)) {
                    ++scoreTkm1;
                } else if (x1_ && x2_ && y1_ && y2_) {
                    --scoreTkm1;
                }
            }

            for (auto &[x_, y_] : unstableEdges) {
                if (x_ == v1_ || x_ == v2_ || y_ == v1_ || y_ == v2_) continue;
                bool x1_ = membershipCnt1[x_], x2_ = membershipCnt2[x_];
                bool y1_ = membershipCnt1[y_], y2_ = membershipCnt2[y_];
                if ((!x1_ && !x2_) || (!y1_ && !y2_)) continue;
                if ((!x1_ && !y2_) || (!x2_ && !y1_)) {
                    ++scoreHUE;
                }
            }
            int resAfterMerger = checkMergerResultIIP(v1_, v2_);
            fout << v1_ << " "
                 << v2_ << " "
                 << resAfterMerger << " "
                 << scoreShell << " "
                 << scoreTkm1 << " "
                 << scoreHUE << " "
                 << commonNbrTkm1 << " "
                 << commonNbrTk << " "
                 << endl;
        }
    }
    bar.finish();
    return 0;
}
