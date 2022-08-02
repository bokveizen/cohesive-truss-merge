# On Improving the Cohesiveness of Graphs by Merging Nodes: Formulation, Analysis, and Algorithms

Source code for the paper **On Improving the Cohesiveness of Graphs by Merging Nodes: Formulation, Analysis, and Algorithms**, where we formulate and study the problem of improving the connectivity and robustness of graphs by merging nodes, for which we use the number of edges in the k-truss for some given k as the objective.
We prove the NP-hardness and non-submodularity of the problem.
For the problem, based on our theoretical findings regarding mergers between nodes and k-trusses, we propose **BATMAN** 
(<ins><strong>B</strong></ins>est-merger se<ins><strong>A</strong></ins>rcher for <ins><strong>T</strong></ins>russ <ins><strong>MA</strong></ins>ximizatio<ins><strong>N</strong></ins>), a fast and effective algorithm equipped with strong search-space-pruning schemes and analyze its time and space complexity.
Through extensive experiments on fourteen real-world graphs, we demonstrate the superiority of BATMAN over several baseline methods and the effectiveness of every component of it. 

Note: if a preview of the supplementary materials PDF file does not appear properly, please download the file.

## Datasets

The datasets can be found in the *datasets_txt* folder.
We make all datasets unweighted and undirected and extract the largest connected component from each dataset.

Source:

- email (EM): [https://snap.stanford.edu/data/email-Eu-core.html](https://snap.stanford.edu/data/email-Eu-core.html)
- facebook (FB): [https://snap.stanford.edu/data/ego-Facebook.html](https://snap.stanford.edu/data/ego-Facebook.html)
- enron (ER): [https://snap.stanford.edu/data/email-Enron.html](https://snap.stanford.edu/data/email-Enron.html)
- brightkite (BK): [https://snap.stanford.edu/data/loc-Brightkite.html](https://snap.stanford.edu/data/loc-Brightkite.html)
- relato (RL): [https://data.world/datasyndrome/relato-business-graph-database](https://data.world/datasyndrome/relato-business-graph-database)
- epinions (EP): [https://snap.stanford.edu/data/soc-Epinions1.html](https://snap.stanford.edu/data/soc-Epinions1.html)
- hepph (HP): [https://snap.stanford.edu/data/cit-HepPh.html](https://snap.stanford.edu/data/cit-HepPh.html)
- slashdot (SD): [https://snap.stanford.edu/data/soc-Slashdot0811.html](https://snap.stanford.edu/data/soc-Slashdot0811.html) 
- syracuse (SC): [https://networkrepository.com/socfb-Syracuse56.php](https://networkrepository.com/socfb-Syracuse56.php)
- gowalla (GW): [https://snap.stanford.edu/data/loc-Gowalla.html](https://snap.stanford.edu/data/loc-Gowalla.html)
- twitter (TT): [https://snap.stanford.edu/data/ego-Twitter.html](https://snap.stanford.edu/data/ego-Twitter.html)
- stanford (SF): [https://snap.stanford.edu/data/web-Stanford.html](https://snap.stanford.edu/data/web-Stanford.html)
- youtube (YT): [https://snap.stanford.edu/data/com-Youtube.html](https://snap.stanford.edu/data/com-Youtube.html)
- wikitalk (WT): [https://snap.stanford.edu/data/wiki-Talk.html](https://snap.stanford.edu/data/wiki-Talk.html)

## Code

The experiments consists of three parts: 
(1) the "main" experiments, (2) the "sampling" experiments, and (3) the "heuristics" experiments.   

### main experiments

Considered algorithms and the corresponding files:

- **BM (BATMAN)**: the proposed method; the code is in *mainbm.cpp*;
- **EQ (BATMAN-EQ)**: a BATMAN variant always **equally distributing** the number of pairs to check between inside-inside mergers and inside-outside mergers; the code ins in *maineq.cpp*;
- **II (BATMAN-II)**: a BATMAN variant considering **inside-inside mergers** only; the code is in *mainii.cpp*;
- **IO (BATMAN-IO)**: a BATMAN variant considering **inside-outside mergers** only; the code is in *mainio.cpp*;
- **NT (most new triangles)**: among all the inside-inside mergers and inside-outside mergers, choosing ones that increase the number of **triangles** consisting of the nodes in the current (k-1)-truss most; the code is in *mainnt.cpp*;
- **NE (most new edges)**: among all the inside-outside mergers, choosing the ones that increase the number of **edges** among the nodes in the current (k-1)-truss most; the code is in *mainne.cpp*;
- **RD (Random)**: uniform random sampling among all the IIMs and IOMs; the code is in *mainrd.cpp*

How to run the code:

	# starting at the root of the project folder; when you do "ls", you should see the cpp files
	# compile all the cpp files for the main experiments
	python compile.py mainbm.cpp maineq.cpp mainii.cpp mainio.cpp mainnt.cpp mainne.cpp mainrd.cpp;
	# create the folder for results
	mkdir -p res_txt/res_main;
	# format: ./main[algorithm:bm/eq/ii/io/nt/ne/rd] [dataset] res_main/[log_file_name] [k] [b] [n_i] [n_o] [n_c] [random_seed]
	# even when the algorithm does not need some parameters, we still keep the same inputs for consistency and alignment
	# k and b are needed for all the algorithms, for the remaining parameters:
	# BM, EQ, IO, and NT use n_i, n_o, n_c
	# II uses n_i, n_c		
	# NE uses n_c
	# RD uses n_c, random_seed
	# an example (BATMAN algorithm; dataset: email;
	# k = 10, b = 10, n_i = 100, n_o = 50, n_c = 10;
	# log stored in res_txt/res_main/email_10_10_10_50_100.txt)
	./mainbm email res_main/email_10_10_10_50_100 10 10 10 50 100 42;

In the log file, the program records the changes in each round.
Specifically, after b rounds, in the last line, the final size of the k-truss for the input k and the total running time are recorded.

### sampling experiments

Considered algorithms and the corresponding files:

- **II (inside-inside)**: random sampling among inside-inside mergers; the code is in *samplii.cpp*;
- **IO (inside-outside)**: random sampling among inside-outside mergers; the code is in *samplio.cpp*;
- **OO (outside-outside)**: random sampling among outside-outside mergers; the code is in *samploo.cpp*

How to run the code:

	# starting at the root of the project folder; when you do "ls", you should see the cpp files
	# compile all the cpp files for the sampling experiments
	python compile.py mainbm.cpp maineq.cpp mainii.cpp mainio.cpp mainnt.cpp mainne.cpp mainrd.cpp;
	# create the folder for results
	mkdir -p res_txt/res_sampling;
	# format: ./sampl[algorithm:ii/io/oo] [dataset] res_sampling/[log_file_name] [k] [b] [n_i] [n_o] [n_c] [random_seed]
	# for the sampling experiments, the algorithms only need k, n_c, and random_seed
	# but we still keep the same inputs for consistency and alignment
	# here, n_c is the number of randomly sampled mergers		
	# an example (II algorithm; dataset: email;
	# k = 10, sampling 10000 mergers;
	# log stored in res_txt/res_sampling/email_10_10000.txt)
	./samplii email res_main/email_10_10000 10 10 100 50 10000 42;

In the log file, the program records the result of each merger.
Specifically, after printing out the information of the original graph,
in each line, the information of 

- (1) first node in the merger,
- (2) second node in the merger,
- (3) performance of the merger (the size of the k-truss after the merger), and
- (4) the best performance among all the above mergers (including the current one)

is recorded.

### heuristics experiments

We have two cpp files, *heuii.cpp* and *heuio.cpp*, for inside-inside mergers and inside-outside mergers, respectively.

How to run the code:

	# starting at the root of the project folder; when you do "ls", you should see the cpp files
	# compile all the cpp files for the sampling experiments
	python compile.py heuii.cpp heuio.cpp
	# create the folder for results
	mkdir -p res_txt/res_heu;
	# format: ./heu[category:ii/io] [dataset] res_heu/[log_file_name] [k] [b] [n_i] [n_o] [n_c] [random_seed]
	# for the sampling experiments, the algorithms only need k, n_c, and random_seed
	# but we still keep the same inputs for consistency and alignment
	# here, n_c is the number of randomly sampled mergers		
	# an example (II algorithm; dataset: email;
	# k = 10, sampling 10000 mergers;
	# log stored in res_txt/res_sampling/email_10_10000.txt)
	./samplii email res_main/email_10_10000 10 10 100 50 10000 42;

In the log file, the program first records the information of the original graph,
then records the inside nodes chosen by each heuristic.
After that, for each merger that uses the chosen nodes (the union of the nodes chosen by all the heuristics),
*heuio.cpp*, the program for inside-outside mergers, records in each line the information of

- (1) the inside node in the merger,
- (2) the outside node in the merger,
- (3) performance of the merger (the size of the k-truss after the merger),
- (4) number of the "new" neighbors that the outside node brings to the inside node, 
- (5) number of the potentially helped unstable edges (shell edges with support k - 3),
- (6) number of the incident potentially helped unstable edges,
- (7) number of the potentially helped shell edges, and
- (8) number of the incident potentially helped shell edges;

and *heuii.cpp*, the program for inside-inside mergers, records in each line the information of

- (1) the first inside node in the merger,
- (2) the second inside node in the merger,
- (3) performance of the merger (the size of the k-truss after the merger),
- (4) the score using the shell edges of the merger, 
- (5) the score using all the edges in the (k-1)-truss of the merger,
- (6) the score using the unstable edges of the merger,
- (7) the number of collisions caused by the merger in the (k-1)-truss, and
- (8) the number of collisions caused by the merger in the k-truss. 
