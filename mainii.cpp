#include "truss.h"

using namespace std;

int main(int argc, char *argv[]) {
    string progName(argv[0]);
    progName = progName.substr(2);
    assert(argc >= 3);
    string dataset(argv[1]);
    string resInfo(argv[2]);
    int argi = 4;
    kInput = (argc >= argi) ? atoi(argv[argi - 1]) : 10;
    ++argi;
    params.emplace_back(kInput);
    int numRounds = (argc >= argi) ? atoi(argv[argi - 1]) : 10;
    ++argi;
    params.emplace_back(numRounds);
    int numV1Check = (argc >= argi) ? atoi(argv[argi - 1]) : 100;
    ++argi;
    params.emplace_back(numV1Check);
    int numV2Check = (argc >= argi) ? atoi(argv[argi - 1]) : 50;  // acutally not used, kept for alignment
    ++argi;
    params.emplace_back(numV2Check);
    int numPairsCheck = (argc >= argi) ? atoi(argv[argi - 1]) : 100;
    ++argi;
    params.emplace_back(numPairsCheck);

    infile = "datasets_txt/" + dataset + ".txt";
    outfile = "res_txt/" + resInfo + ".txt";
    fout.open(outfile.c_str(), ios_base::app);
    // print the basic information of the experimental settings
    cout << progName << " ";
    fout << progName << " ";
    cout << dataset << " ";
    fout << dataset << " ";
    params.resize(6);  // filled with 0
    for (auto param_: params) {
        cout << param_ << " ";
        fout << param_ << " ";
    }
    cout << endl;
    fout << endl;

    auto bIIPs = numPairsCheck;  // IIPs only

    for (int i = 0; i <= numRounds; ++i) {
        cout << "round " << i << endl;
        fout << "round " << i << endl;
        auto startTimeZ = chrono::system_clock::now();  // time of truss update, etc.
        if (i == 0) {
            readOrderedSimpleGraph();
            countTriangles(supp);
            suppOrig = supp;  // supp backup
        } else {
            varsReInit();
        }
        binSort(supp);
        trussDecomp();
        if (i == 0) {
            startTime = chrono::system_clock::now();
        }
        fillInNodesEdges();
        printNodeEdgeInfo();
        printNodeEdgeInfoFile();
        if (i == numRounds) {  // last round
            cout << "final k-truss size = " << mTk << endl;
            fout << mTk << " ";
            break;
        }
        auto endTimeZ = chrono::system_clock::now();
        chrono::duration<double> elapsedSecondsZ = endTimeZ - startTimeZ;
        sprintf(timeStr, "%.6f", elapsedSecondsZ.count());
        cout << "time of truss update: " << timeStr << endl;
        fout << "time of truss update: " << timeStr << endl;

        auto startTimeA = chrono::system_clock::now();
        computeNbrsInTkm1();
        computeNbrsInTk();
        // findUnstableEdges();
        findBestOutsideNodes();
        // computeHUEs();
        // computeHNNs();

        int chosenV1, chosenV2;
        int bestRes = 0, bestResIOPs = 0, bestResIIPs = 0;
        VI &v1Pruned = Tkm1Nodes;
        int n1 = (int) v1Pruned.size();
        sort(v1Pruned.begin(), v1Pruned.end(), compMoreIncidentPotential);
        VI &v2Pruned = bestOutsideNodes;
        int n2 = (int) v2Pruned.size();
        sort(v2Pruned.begin(), v2Pruned.end(), compMoreNbrsInTkm1);

        numV1Check = (numV1Check == -1) ? n1 : min(numV1Check, n1);
        numV2Check = (numV2Check == -1) ? n2 : min(numV2Check, n2);

        auto endTimeA = chrono::system_clock::now();
        chrono::duration<double> elapsedSecondsA = endTimeA - startTimeA;
        sprintf(timeStr, "%.6f", elapsedSecondsA.count());
        cout << "time of preparation: " << timeStr << endl;
        fout << "time of preparation: " << timeStr << endl;

        // IIPs
        VE bestPairsIIPs(bIIPs, TEdge{-1, -1});
        VI bestScoresIIPs(bIIPs, -mTotal);
        auto numIIPs = numV1Check * (numV1Check - 1) / 2;
        auto startTimeB = chrono::system_clock::now();  // time of IIPs-finding pairs
        tqdm bar;
        int curr = 0;
        for (auto i_ = 0; i_ < numV1Check - 1; ++i_) {
            int v1_ = v1Pruned[i_];
            VI &n1_ = nbrsInTkm1[v1_];
            vector<bool> membershipCnt1(n, false);
            vector<bool> membershipCntK1(n, false);
            for (auto &x_: n1_) {
                membershipCnt1[x_] = true;
            }
            for (auto &x_: nbrsInTk[v1_]) {
                membershipCntK1[x_] = true;
            }
            for (auto j_ = i_ + 1; j_ < numV1Check; ++j_) {
                bar.progress(curr++, numIIPs);
                auto v2_ = v1Pruned[j_];
                VI &n2_ = nbrsInTkm1[v2_];
                vector<bool> membershipCnt2(n, false);
                int scoreV1V2 = 0;
                for (auto &x_: n2_) {
                    membershipCnt2[x_] = true;
                }
                for (auto &x_: nbrsInTk[v2_]) {
                    if (x_ == v1_ || membershipCntK1[x_]) --scoreV1V2;
                }
                for (auto &[x_, y_]: Skm1Edges) {
                    if (x_ == v1_ || x_ == v2_ || y_ == v1_ || y_ == v2_) {  // only checking non-incident edges
                        continue;
                    }
                    bool x1_ = membershipCnt1[x_], x2_ = membershipCnt2[x_];
                    bool y1_ = membershipCnt1[y_], y2_ = membershipCnt2[y_];
                    if ((!x1_ && !x2_) || (!y1_ && !y2_)) {
                        continue;
                    }
                    if ((!x1_ && !y2_) || (!x2_ && !y1_)) {
                        ++scoreV1V2;
                    } else if (x1_ && x2_ && y1_ && y2_) {
                        --scoreV1V2;
                    }
                }
                auto minScorePos = min_element(bestScoresIIPs.begin(), bestScoresIIPs.end());
                auto minScore = *minScorePos;
                auto minScoreIndex = minScorePos - bestScoresIIPs.begin();
                if (scoreV1V2 > minScore) {
                    bestPairsIIPs[minScoreIndex] = TEdge{v1_, v2_};
                    bestScoresIIPs[minScoreIndex] = scoreV1V2;
                }
            }
        }
        bar.finish();
        auto endTimeB = chrono::system_clock::now();
        chrono::duration<double> elapsedSecondsB = endTimeB - startTimeB;
        sprintf(timeStr, "%.6f", elapsedSecondsB.count());
        cout << "time of IIPs-finding pairs: " << timeStr << endl;
        fout << "time of IIPs-finding pairs: " << timeStr << endl;

        auto startTimeC = chrono::system_clock::now();  // IIPs-checking pairs
        if (numPairsCheck == 1) {
            auto &[bpu, bpv] = bestPairsIIPs[0];
            chosenV1 = bpu;
            chosenV2 = bpv;
        } else {
            for (auto &[bpu, bpv]: bestPairsIIPs) {
                if (bpu == -1 || bpv == -1)
                    continue;
                int resAfterMerger = checkMergerResultIIP(bpu, bpv);
                if (resAfterMerger > bestRes) {
                    bestRes = resAfterMerger;
                    chosenV1 = bpu;
                    chosenV2 = bpv;
                }
            }
        }
        auto endTimeC = chrono::system_clock::now();
        chrono::duration<double> elapsedSecondsC = endTimeC - startTimeC;
        sprintf(timeStr, "%.6f", elapsedSecondsC.count());
        cout << "time of IIPs-checking pairs: " << timeStr << endl;
        fout << "time of IIPs-checking pairs: " << timeStr << endl;

        auto startTimeD = chrono::system_clock::now();  // time of update the graph
        if (i == numRounds - 1) {  // all pairs have been decided
            endTime = chrono::system_clock::now();
        }
        updateGraph(chosenV1, chosenV2);
        cout << "merging " << chosenV1 << " and " << chosenV2 << endl;
        fout << "merging " << chosenV1 << " and " << chosenV2 << endl;
        supp = suppOrig;
        auto endTimeD = chrono::system_clock::now();
        chrono::duration<double> elapsedSecondsD = endTimeD - startTimeD;
        sprintf(timeStr, "%.6f", elapsedSecondsD.count());
        cout << "time of updating the graph: " << timeStr << endl;
        fout << "time of updating the graph: " << timeStr << endl;
    }
    chrono::duration<double> elapsedSeconds = endTime - startTime;
    sprintf(timeStr, "%.6f", elapsedSeconds.count());
    cout << "total runtime = " << timeStr << " seconds" << endl;
    fout << timeStr << endl;
    fout.close();
    return 0;
}
