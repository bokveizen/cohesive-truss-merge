#include <random>

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
    int numRounds = (argc >= argi) ? atoi(argv[argi - 1]) : 10;  // acutally not used, kept for alignment
    ++argi;
    params.emplace_back(numRounds);
    int numV1Check = (argc >= argi) ? atoi(argv[argi - 1]) : 100;  // acutally not used, kept for alignment
    ++argi;
    params.emplace_back(numV1Check);
    int numV2Check = (argc >= argi) ? atoi(argv[argi - 1]) : 50;  // acutally not used, kept for alignment
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
    readOrderedSimpleGraph();
    countTriangles(supp);
    suppOrig = supp;  // supp backup
    binSort(supp);
    trussDecomp();
    startTime = chrono::system_clock::now();
    fillInNodesEdges();
    printNodeEdgeInfo();
    printNodeEdgeInfoFile();

    computeNbrsInTkm1();
    findBestOutsideNodes();

//    VI &v1Pruned = Tkm1Nodes;
//    int n1 = (int) v1Pruned.size();
//    uniform_int_distribution<> v1Rand(0, n1);
    VI &v2Pruned = bestOutsideNodes;
    int n2 = (int) v2Pruned.size();
    uniform_int_distribution<> v2Rand(0, n2);

    int iPair = 0, bestRes = 0;
    tqdm bar;
    while (iPair < numPairsCheck) {
        bar.progress(iPair, numPairsCheck);
        int sampleV1 = -1, sampleV2 = -1;
        while (sampleV1 == sampleV2) {
            sampleV1 = v2Pruned[v2Rand(gen)];
            sampleV2 = v2Pruned[v2Rand(gen)];
        }
        auto resV1V2 = checkMergerResultGeneral(sampleV1, sampleV2);
        if (resV1V2 > bestRes) {
            bestRes = resV1V2;
        }
        fout << sampleV1 << " "
             << sampleV2 << " "
             << resV1V2 << " "
             << bestRes << endl;
        ++iPair;
    }
    bar.finish();
    return 0;
}
