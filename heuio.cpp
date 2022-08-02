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
    suppOrig = supp; // supp backup
    binSort(supp);
    trussDecomp();
    startTime = chrono::system_clock::now();
    fillInNodesEdges();
    printNodeEdgeInfoFile();

    computeNbrsInTkm1();
    findUnstableEdges();
    findBestOutsideNodes();
    computeHUEs();

    VI &v1Pruned = Tkm1Nodes;
    int maxV1 = 100;
    maxV1 = min((int) v1Pruned.size(), maxV1);
    unordered_set<int> bestV1Total;
    bestV1Total.reserve(maxV1 * 6);

    sort(v1Pruned.begin(), v1Pruned.end(), compMoreNbrsInTkm1);
    fout << "bestV1Nbrs" << endl;
    for (auto i_ = 0; i_ < maxV1; ++i_) {
        auto vi = v1Pruned[i_];
        bestV1Total.insert(vi);
        fout << vi << " ";
    }
    fout << endl;

    sort(v1Pruned.begin(), v1Pruned.end(), compMoreHUEs);
    fout << "bestV1HUEs" << endl;
    for (auto i_ = 0; i_ < maxV1; ++i_) {
        auto vi = v1Pruned[i_];
        bestV1Total.insert(vi);
        fout << vi << " ";
    }
    fout << endl;

    sort(v1Pruned.begin(), v1Pruned.end(), compMoreIncidentPotential);
    fout << "bestV1IPs" << endl;
    for (auto i_ = 0; i_ < maxV1; ++i_) {
        auto vi = v1Pruned[i_];
        bestV1Total.insert(vi);
        fout << vi << " ";
    }
    fout << endl;

    gen.seed(0);
    shuffle(v1Pruned.begin(), v1Pruned.end(), gen);
    fout << "random 0" << endl;
    for (auto i_ = 0; i_ < maxV1; ++i_) {
        auto vi = v1Pruned[i_];
        bestV1Total.insert(vi);
        fout << vi << " ";
    }
    fout << endl;

    gen.seed(1);
    shuffle(v1Pruned.begin(), v1Pruned.end(), gen);
    fout << "random 1" << endl;
    for (auto i_ = 0; i_ < maxV1; ++i_) {
        auto vi = v1Pruned[i_];
        bestV1Total.insert(vi);
        fout << vi << " ";
    }
    fout << endl;

    gen.seed(2);
    shuffle(v1Pruned.begin(), v1Pruned.end(), gen);
    fout << "random 2" << endl;
    for (auto i_ = 0; i_ < maxV1; ++i_) {
        auto vi = v1Pruned[i_];
        bestV1Total.insert(vi);
        fout << vi << " ";
    }
    fout << endl;

    VI &v2Pruned = bestOutsideNodes;
    sort(v2Pruned.begin(), v2Pruned.end(), compMoreNbrsInTkm1);
    int maxV2 = 50;
    maxV2 = min((int) v2Pruned.size(), maxV2);
    fout << "outside nodes" << endl;
    for (auto i_ = 0; i_ < maxV2; ++i_) {
        fout << v2Pruned[i_] << " ";
    }
    fout << endl;

    tqdm bar;
    int curr = 0;
    int totalPairs = (int) bestV1Total.size() * maxV2;
    for (auto v1_: bestV1Total) {
        vector<bool> bestN1(n, false), bestNN1(n, false);
        VI baseNewNbrs;
        baseNewNbrs.reserve(nbrsInTkm1[v1_].size());
        for (auto &[x_, _]: suppTkm1[v1_]) bestNN1[x_] = true;
        for (auto &x_: nbrsInTkm1[v1_]) {
            if (!bestNN1[x_]) baseNewNbrs.emplace_back(x_);
            bestN1[x_] = true;
        }
        unordered_map<int, VI> n2HUE, n2HSE;
        n2HUE.reserve(n);
        n2HSE.reserve(n);
        unordered_set<int> baseHUE, baseHSE;
        baseHUE.reserve(v2UsNbrs[v1_].size());
        baseHSE.reserve(v2ShNbrs[v1_].size());
        for (auto &n1_: v2UsNbrs[v1_]) {
            for (auto &[n2_, _]: suppTkm1[n1_]) {
                if (!bestN1[n2_]) {
                    n2HUE[n2_].emplace_back(n1_);
                }
            }
        }
        for (auto nn_: baseNewNbrs) {
            baseHUE.insert(n2HUE[nn_].begin(), n2HUE[nn_].end());
        }

        for (auto &n1_: v2ShNbrs[v1_]) {
            for (auto &[n2_, _]: suppTkm1[n1_]) {
                if (!bestN1[n2_]) {
                    n2HSE[n2_].emplace_back(n1_);
                }
            }
        }
        for (auto nn_: baseNewNbrs) {
            baseHSE.insert(n2HSE[nn_].begin(), n2HSE[nn_].end());
        }

        for (auto j_ = 0; j_ < maxV2; ++j_) {
            bar.progress(curr++, totalPairs);
            int v2_ = v2Pruned[j_];
            VI &n2_ = nbrsInTkm1[v2_];
            VI newNbrsV2 = baseNewNbrs;
            newNbrsV2.reserve(n2_.size() + baseNewNbrs.size());
            vector<bool> newNbrs(n, false);
            bool e12 = false;

            unordered_set<int> HUEv1v2Incident = baseHUE, HSEv1v2Incident = baseHSE;
            HUEv1v2Incident.reserve(v2UsNbrs[v1_].size());
            HSEv1v2Incident.reserve(v2ShNbrs[v1_].size());
            for (auto x_: n2_) {
                if (x_ == v1_) {
                    e12 = true;
                    continue;
                }
                if (!bestN1[x_]) {
                    newNbrs[x_] = true;
                    newNbrsV2.emplace_back(x_);
                    HUEv1v2Incident.insert(n2HUE[x_].begin(), n2HUE[x_].end());
                    HSEv1v2Incident.insert(n2HSE[x_].begin(), n2HSE[x_].end());
                }
            }
            auto numOfPHUEsV1V2Incident = HUEv1v2Incident.size();
            auto numOfPHSEsV1V2Incident = HSEv1v2Incident.size();

            auto numNewNbrsV1V2 = newNbrsV2.size();
            if (e12) --numNewNbrsV1V2;

            size_t numOfPHUEsV1V2 = 0;
            // unstable edges

            for (auto &[x_, y_]: unstableEdges) {
                if (x_ == v1_ || y_ == v1_) { // todo: more reasonable computation
//                    ++numOfPHUEsV1V2;
//                    ++numOfPHUEsV1V2Incident;
                    continue;
                }
                bool xNN = bestNN1[x_];
                bool yNN = bestNN1[y_];
                bool xN = bestN1[x_];
                bool yN = bestN1[y_];
                bool xNew = newNbrs[x_];
                bool yNew = newNbrs[y_];
                if (xNN && yNN) continue;
                if ((xN || xNew) && (yN || yNew)) ++numOfPHUEsV1V2;
            }
            numOfPHUEsV1V2 += numOfPHUEsV1V2Incident;

            // shell edges
            size_t numOfPHSEsV1V2 = 0;
            for (auto &[x_, y_]: Skm1Edges) {
                if (x_ == v1_ || y_ == v1_) { // todo: more reasonable computation
//                    ++numOfPHSEsV1V2;
//                    ++numOfPHSEsV1V2Incident;
                    continue;
                }
                bool xNN = bestNN1[x_];
                bool yNN = bestNN1[y_];
                bool xN = bestN1[x_];
                bool yN = bestN1[y_];
                bool xNew = newNbrs[x_];
                bool yNew = newNbrs[y_];
                if (xNN && yNN) continue;
                if ((xN || xNew) && (yN || yNew)) ++numOfPHSEsV1V2;
            }
            numOfPHSEsV1V2 += numOfPHSEsV1V2Incident;

            int resAfterMerger = checkMergerResultIOP(v1_, v2_);
            fout << v1_ << " " << v2_ << " "
                 << resAfterMerger << " "
                 << numNewNbrsV1V2 << " "
                 << numOfPHUEsV1V2 << " "
                 << numOfPHUEsV1V2Incident << " "
                 << numOfPHSEsV1V2 << " "
                 << numOfPHSEsV1V2Incident << " "
                 << endl;
        }
    }
    bar.finish();
    return 0;
}
