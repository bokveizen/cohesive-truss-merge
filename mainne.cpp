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
//        findBestOutsideNodes();

        int chosenV1, chosenV2;
        int bestRes = -2, bestResIOPs = -2, bestResIIPs = -2;
        VI &v1Pruned = Tkm1Nodes;
        int n1 = (int) v1Pruned.size();
        sort(v1Pruned.begin(), v1Pruned.end(), compMoreIncidentPotential);
//        VI &v2Pruned = bestOutsideNodes;
        VI &v2Pruned = outsideNodes;
        int n2 = (int) v2Pruned.size();
        sort(v2Pruned.begin(), v2Pruned.end(), compMoreNbrsInTkm1);

//        numV1Check = (numV1Check == -1) ? n1 : min(numV1Check, n1);
//        numV2Check = (numV2Check == -1) ? n2 : min(numV2Check, n2);

        auto endTimeA = chrono::system_clock::now();
        chrono::duration<double> elapsedSecondsA = endTimeA - startTimeA;
        sprintf(timeStr, "%.6f", elapsedSecondsA.count());
        cout << "time of preparation: " << timeStr << endl;
        fout << "time of preparation: " << timeStr << endl;

        // IOPs
        auto startTimeB = chrono::system_clock::now();  // time of IOPs-finding pairs
//        auto numIOPs = numV1Check * numV2Check;
//        VE bestPairsIOPs(bIOPs, TEdge{-1, -1});
//        vector<tuple<size_t, size_t>> bestScoresIOPs(bIOPs, {0, 0});
//        vector<VI> bestNewNbrsList(bIOPs);
        VE bestPairsIOPs(bIOPs, TEdge{-1, -1});
        VI bestScoresIOPs(bIOPs, -1);
//        tqdm bar;
//        int curr = 0;
        int curMostNEs = 0;
        for (auto &v2_: v2Pruned) {
//            cout << "checking " << v2_ << endl;
            auto &n2_ = nbrsInTkm1[v2_];
            int n2Size = (int) n2_.size();
            if (n2Size <= *min_element(bestScoresIOPs.begin(), bestScoresIOPs.end())) {
                break;
            }
            for (auto &v1_: v1Pruned) {
                auto &n1_ = nbrsInTkm1[v1_];
                VI newN_;
                newN_.reserve(n2_.size());
                set_difference(n2_.begin(), n2_.end(),
                               n1_.begin(), n1_.end(),
                               back_inserter(newN_));
                int ne_ = (int) newN_.size();
                if (suppOrig[v2_].contains(v1_)) {
                    --ne_;
                }
                auto minScorePos = min_element(bestScoresIOPs.begin(), bestScoresIOPs.end());
                auto minScore = *minScorePos;
                auto minScoreIndex = minScorePos - bestScoresIOPs.begin();
                if (minScore >= n2Size) {
                    break;
                }
                if (ne_ > minScore) {
                    bestPairsIOPs[minScoreIndex] = TEdge{v1_, v2_};
                    bestScoresIOPs[minScoreIndex] = ne_;
                }
            }
        }
        auto endTimeB = chrono::system_clock::now();
        chrono::duration<double> elapsedSecondsB = endTimeB - startTimeB;
        sprintf(timeStr, "%.6f", elapsedSecondsB.count());
        cout << "time of finding the best pair: " << timeStr << endl;
        fout << "time of finding the best pair: " << timeStr << endl;

        auto startTimeC = chrono::system_clock::now();  // time of IOPs-checking pairs
        for (auto &[bestV1, bestV2]: bestPairsIOPs) {
            if (bestV1 == -1 || bestV2 == -1) {
                continue;
            }
            int resV1V2 = checkMergerResultGeneral(bestV1, bestV2);
            if (resV1V2 > bestResIOPs) {
                bestResIOPs = resV1V2;
                chosenV1 = bestV1;
                chosenV2 = bestV2;
            }
        }

        auto startTimeD = chrono::system_clock::now();  // time of update the graph        
        if (i == numRounds - 1) {  // all pairs have been decided
            endTime = chrono::system_clock::now();
        }
        cout << "merging " << chosenV1 << " and " << chosenV2 << endl;
        fout << "merging " << chosenV1 << " and " << chosenV2 << endl;
        updateGraph(chosenV1, chosenV2);
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
