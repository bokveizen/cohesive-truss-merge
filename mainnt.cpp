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
    int numV2Check = (argc >= argi) ? atoi(argv[argi - 1]) : 50;
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
        // findUnstableEdges();
        // findBestOutsideNodes();
        // computeHUEs();
        // computeHNNs();


        int chosenV1, chosenV2;
        int bestRes = -2, bestResIOPs = 0, bestResIIPs = 0;
        VI &v1Pruned = Tkm1Nodes;
        int n1 = (int) v1Pruned.size();
        sort(v1Pruned.begin(), v1Pruned.end(), compMoreNbrsInTkm1);
//        VI &v2Pruned = bestOutsideNodes;
        VI &v2Pruned = outsideNodes;
        int n2 = (int) v2Pruned.size();
        sort(v2Pruned.begin(), v2Pruned.end(), compMoreNbrsInTkm1);

        numV1Check = (numV1Check == -1) ? n1 : min(numV1Check, n1);
        numV2Check = (numV2Check == -1) ? n2 : min(numV2Check, n2);

        VE bestPairs(numPairsCheck, TEdge{-1, -1});
        VI bestScores(numPairsCheck, 0);

        // compute the number of possible new triangles of each inside node
        MII v2PNT;  // number of possible new triangles
        int checkPairsPNT = min(10 * numPairsCheck, (int) v1Pruned.size());
        v2PNT.reserve(checkPairsPNT);
        for (auto i_ = 0; i_ < checkPairsPNT; ++i_) {
            auto &v1_ = v1Pruned[i_];
            vector<bool> bestN1(n, false);
            for (auto &x_: nbrsInTkm1[v1_]) {
                bestN1[x_] = true;
            }
            for (auto &x_: nbrsInTkm1[v1_]) {
                for (auto &y_: nbrsInTkm1[x_]) {
                    if (!bestN1[y_]) {
                        ++v2PNT[v1_];
                    }
                }
            }
        }

        sort(v1Pruned.begin(), v1Pruned.end(),
             [&](const int &a_, const int &b_) {
                 if (v2PNT[a_] != v2PNT[b_]) {
                     return v2PNT[a_] > v2PNT[b_];
                 }
                 if (nbrsInTkm1[a_].size() != nbrsInTkm1[b_].size()) {
                     return nbrsInTkm1[a_].size() > nbrsInTkm1[b_].size();
                 }
                 return a_ < b_;
             }
        );

        auto endTimeA = chrono::system_clock::now();
        chrono::duration<double> elapsedSecondsA = endTimeA - startTimeA;
        sprintf(timeStr, "%.6f", elapsedSecondsA.count());
        cout << "time of preparation: " << timeStr << endl;
        fout << "time of preparation: " << timeStr << endl;

        auto startTimeB = chrono::system_clock::now();
        auto numIOPs = numV1Check * numV2Check;

        int curMostNTs = 0;
        for (auto i_ = 0; i_ < numV1Check; ++i_) {
            auto &v1_ = v1Pruned[i_];
            if (v2PNT[v1_] <= *min_element(bestScores.begin(), bestScores.end())) {
                break;
            }
            vector<bool> bestN1(n, false);
            for (auto &x_: nbrsInTkm1[v1_]) {
                bestN1[x_] = true;
            }
            for (auto j_ = i_ + 1; j_ < numV1Check; ++j_) {
                int scoreV1V2 = 0;
                auto &v2_ = v1Pruned[j_];
                vector<bool> bestN2(n, false);
                for (auto &n2_: nbrsInTkm1[v2_]) {
                    bestN2[n2_] = true;
                }
                for (auto &[x_, y_]: Tkm1Edges) {
                    if (x_ == v1_ || x_ == v2_ || y_ == v1_ || y_ == v2_) {  // only checking non-incident edges
                        continue;
                    }
                    bool x1_ = bestN1[x_], x2_ = bestN2[x_];
                    bool y1_ = bestN1[y_], y2_ = bestN2[y_];
                    if ((!x1_ && !x2_) || (!y1_ && !y2_)) {
                        continue;
                    }
                    if ((!x1_ && !y2_) || (!x2_ && !y1_)) {
                        ++scoreV1V2;
                    } else if (x1_ && x2_ && y1_ && y2_) {
                        --scoreV1V2;
                    }
                }

                if (scoreV1V2 > curMostNTs) {
                    chosenV1 = v1_;
                    chosenV2 = v2_;
                    curMostNTs = scoreV1V2;
                }
            }
            for (auto j_ = 0; j_ < numV2Check; ++j_) {
                int scoreV1V2 = 0;
                auto &v2_ = v2Pruned[j_];
                vector<bool> bestN2(n, false);
                for (auto &n2_: nbrsInTkm1[v2_]) {
                    bestN2[n2_] = true;
                }
                for (auto &[x_, y_]: Tkm1Edges) {
                    if (x_ == v1_ || x_ == v2_ || y_ == v1_ || y_ == v2_) {  // only checking non-incident edges
                        continue;
                    }
                    bool x1_ = bestN1[x_], x2_ = bestN2[x_];
                    bool y1_ = bestN1[y_], y2_ = bestN2[y_];
                    if ((!x1_ && !x2_) || (!y1_ && !y2_)) {
                        continue;
                    }
                    if ((!x1_ && !y2_) || (!x2_ && !y1_)) {
                        ++scoreV1V2;
                    } else if (x1_ && x2_ && y1_ && y2_) {
                        --scoreV1V2;
                    }
                }
                auto minScorePos = min_element(bestScores.begin(), bestScores.end());
                auto minScore = *minScorePos;
                auto minScoreIndex = minScorePos - bestScores.begin();
                if (scoreV1V2 > minScore) {
                    bestPairs[minScoreIndex] = TEdge{v1_, v2_};
                    bestScores[minScoreIndex] = scoreV1V2;
                }
            }
        }
        auto endTimeB = chrono::system_clock::now();
        chrono::duration<double> elapsedSecondsB = endTimeB - startTimeB;
        sprintf(timeStr, "%.6f", elapsedSecondsB.count());
        cout << "time of finding the best pair: " << timeStr << endl;
        fout << "time of finding the best pair: " << timeStr << endl;

        auto startTimeC = chrono::system_clock::now();  // time of IOPs-checking pairs
        for (auto ip = 0; ip < numPairsCheck; ++ip) {
            auto &[bestV1, bestV2] = bestPairs[ip];
            if (bestV1 == -1 || bestV2 == -1) {
                continue;
            }
            int resV1V2 = checkMergerResultGeneral(bestV1, bestV2);
            if (resV1V2 > bestRes) {
                bestRes = resV1V2;
                chosenV1 = bestV1;
                chosenV2 = bestV2;
            }
        }
        auto endTimeC = chrono::system_clock::now();
        chrono::duration<double> elapsedSecondsC = endTimeC - startTimeC;
        sprintf(timeStr, "%.6f", elapsedSecondsC.count());
        cout << "time of checking pairs: " << timeStr << endl;
        fout << "time of checking pairs: " << timeStr << endl;

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
