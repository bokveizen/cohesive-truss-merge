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

    auto bIOPs = numPairsCheck;

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

        // IOPs
        auto startTimeB = chrono::system_clock::now();  // time of IOPs-finding pairs
        auto numIOPs = numV1Check * numV2Check;
        VE bestPairsIOPs(bIOPs, TEdge{-1, -1});
        vector<tuple<size_t, size_t>> bestScoresIOPs(bIOPs, {0, 0});
        vector<VI> bestNewNbrsList(bIOPs);
        tqdm bar;
        int curr = 0;
        for (auto i_ = 0; i_ < numV1Check; ++i_) {
            int v1_ = v1Pruned[i_];
            vector<bool> bestN1(n, false), bestNN1(n, false);
            VI baseNewNbrs;  // incident prospects
            baseNewNbrs.reserve(nbrsInTkm1[v1_].size());
            for (auto &[x_, _]: suppTkm1[v1_]) {
                bestNN1[x_] = true;
            }
            for (auto &x_: nbrsInTkm1[v1_]) {
                if (!bestNN1[x_]) {
                    baseNewNbrs.emplace_back(x_);
                }
                bestN1[x_] = true;
            }
            unordered_map<int, VI> n2HSE;
            n2HSE.reserve(n);
            unordered_set<int> baseHSE;
            baseHSE.reserve(v2ShNbrs[v1_].size());
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
            for (auto j_ = 0; j_ < numV2Check; ++j_) {
                bar.progress(curr++, numIOPs);
                auto v2_ = v2Pruned[j_];
                VI &n2_ = nbrsInTkm1[v2_];
                VI newNbrsV1V2 = baseNewNbrs;
                newNbrsV1V2.reserve(n2_.size() + baseNewNbrs.size());
                vector<bool> newNbrs(n, false);
                bool e12 = false;  // edge between v1 and v2 exist?
                unordered_set<int> HSEv1v2Incident = baseHSE;
                HSEv1v2Incident.reserve(v2ShNbrs[v1_].size());
                for (auto x_: n2_) {
                    if (x_ == v1_) {
                        e12 = true;
                        continue;
                    }
                    if (!bestN1[x_]) {
                        newNbrs[x_] = true;
                        newNbrsV1V2.emplace_back(x_);
                        HSEv1v2Incident.insert(n2HSE[x_].begin(), n2HSE[x_].end());
                    }
                }
                auto numOfPHSEsV1V2Incident = HSEv1v2Incident.size();
                auto numNewNbrsV1V2 = newNbrsV1V2.size();
                if (e12) {
                    --numNewNbrsV1V2;
                }
                // shell edges
                size_t numOfPHSEsV1V2 = 0;
                for (auto &[x_, y_]: Skm1Edges) {
                    if (x_ == v1_ || y_ == v1_) {
                        continue;
                    }
                    bool xNN = bestNN1[x_];
                    bool yNN = bestNN1[y_];
                    bool xN = bestN1[x_];
                    bool yN = bestN1[y_];
                    bool xNew = newNbrs[x_];
                    bool yNew = newNbrs[y_];
                    if (xNN && yNN) {
                        continue;
                    }
                    if ((xN || xNew) && (yN || yNew)) {
                        ++numOfPHSEsV1V2;
                    }
                }
                numOfPHSEsV1V2 += numOfPHSEsV1V2Incident;

                tuple<size_t, size_t> scoreV1V2 = make_tuple(numOfPHSEsV1V2, numNewNbrsV1V2);
                auto minScorePos = min_element(bestScoresIOPs.begin(), bestScoresIOPs.end());
                auto minScore = *minScorePos;
                auto minScoreIndex = minScorePos - bestScoresIOPs.begin();
                if (scoreV1V2 > minScore) {
                    bestPairsIOPs[minScoreIndex] = TEdge{v1_, v2_};
                    bestScoresIOPs[minScoreIndex] = scoreV1V2;
                    bestNewNbrsList[minScoreIndex] = newNbrsV1V2;
                }
            }
        }
        bar.finish();
        auto endTimeB = chrono::system_clock::now();
        chrono::duration<double> elapsedSecondsB = endTimeB - startTimeB;
        sprintf(timeStr, "%.6f", elapsedSecondsB.count());
        cout << "time of IOPs-finding pairs: " << timeStr << endl;
        fout << "time of IOPs-finding pairs: " << timeStr << endl;

        auto startTimeC = chrono::system_clock::now();  // time of IOPs-checking pairs
        if (bIOPs == 1) {
            auto &[bestV1, bestV2] = bestPairsIOPs[0];
            chosenV1 = bestV1;
            chosenV2 = bestV2;
        } else {
            for (auto ip = 0; ip < bIOPs; ++ip) {
                auto &[bestV1, bestV2] = bestPairsIOPs[ip];
                if (bestV1 == -1 || bestV2 == -1) {
                    continue;
                }
                int resV1V2 = checkMergerResultIOPInputNewNbrs(bestV1, bestV2, bestNewNbrsList[ip]);
                if (resV1V2 > bestResIOPs) {
                    bestResIOPs = resV1V2;
                    chosenV1 = bestV1;
                    chosenV2 = bestV2;
                }
            }
        }
        bestRes = bestResIOPs;
        auto endTimeC = chrono::system_clock::now();
        chrono::duration<double> elapsedSecondsC = endTimeC - startTimeC;
        sprintf(timeStr, "%.6f", elapsedSecondsC.count());
        cout << "time of IOPs-checking pairs: " << timeStr << endl;
        fout << "time of IOPs-checking pairs: " << timeStr << endl;

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
