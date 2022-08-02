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
    int randomSeed = (argc >= argi) ? atoi(argv[argi - 1]) : 42;
    ++argi;
    params.emplace_back(randomSeed);

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
    gen.seed(randomSeed);

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
        int bestRes = -2;
        VI &v1Pruned = Tkm1Nodes;
        int n1 = (int) v1Pruned.size();
        //        sort(v1Pruned.begin(), v1Pruned.end(), compMoreIncidentPotential);
        uniform_int_distribution<> v1Rand(0, n1);
        VI &v2Pruned = bestOutsideNodes;
        int n2 = (int) v2Pruned.size();
        //        sort(v2Pruned.begin(), v2Pruned.end(), compMoreNbrsInTkm1);
        uniform_int_distribution<> v2Rand(0, n2);
        numV1Check = (numV1Check == -1) ? n1 : min(numV1Check, n1);
        numV2Check = (numV2Check == -1) ? n2 : min(numV2Check, n2);

        auto endTimeA = chrono::system_clock::now();
        chrono::duration<double> elapsedSecondsA = endTimeA - startTimeA;
        sprintf(timeStr, "%.6f", elapsedSecondsA.count());
        cout << "time of preparation: " << timeStr << endl;
        fout << "time of preparation: " << timeStr << endl;

        auto bIOPs = numPairsCheck >> 1, bIIPs = numPairsCheck >> 1;  // initially the budget is equally distributed
        auto nIOPs = n1 * n2;
        auto nIIPs = n1 * (n1 - 1) / 2;
        if (numPairsCheck == 1) {
            bernoulli_distribution bDist((float) nIOPs / (float) (nIOPs + nIIPs));
            auto usingIOPs = bDist(gen);
            if (usingIOPs) {
                ++bIOPs;
            } else {
                ++bIIPs;
            }
        }

        // IOPs
        for (auto ip_ = 0; ip_ < bIOPs; ++ip_) {
            int sampleV1 = -1, sampleV2 = -1;
            while (sampleV1 == sampleV2) {
                sampleV1 = v1Pruned[v1Rand(gen)];
                sampleV2 = v2Pruned[v2Rand(gen)];
                if (sampleV1 < 0 || sampleV1 > n || sampleV2 < 0 || sampleV2 > n) {
                    sampleV1 = -1;
                    sampleV2 = -1;
                    continue;
                }
            }
            auto resV1V2 = checkMergerResultGeneral(sampleV1, sampleV2);
            if (resV1V2 > bestRes) {
                bestRes = resV1V2;
                chosenV1 = sampleV1;
                chosenV2 = sampleV2;
            }
        }

        // IIPs
        for (auto ip_ = 0; ip_ < bIIPs; ++ip_) {
            int sampleV1 = -1, sampleV2 = -1;
            while (sampleV1 == sampleV2) {
                sampleV1 = v1Pruned[v1Rand(gen)];
                sampleV2 = v1Pruned[v1Rand(gen)];
                if (sampleV1 < 0 || sampleV1 > n || sampleV2 < 0 || sampleV2 > n) {
                    sampleV1 = -1;
                    sampleV2 = -1;
                    continue;
                }
            }
            auto resV1V2 = checkMergerResultGeneral(sampleV1, sampleV2);
            if (resV1V2 > bestRes) {
                bestRes = resV1V2;
                chosenV1 = sampleV1;
                chosenV2 = sampleV2;
            }
        }

        auto startTimeD = chrono::system_clock::now();  // time of update the graph
        if (i == numRounds - 1) {                       // all pairs have been decided
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
