//============================================================================
// Name        : BranchReduction.cpp
// Author      : Jasmijn Baaijens
// Version     : 0.4.0
// License     : GNU GPL v3.0
// Project     : ViralQuasispecies
// Description : Use read evidence to reduce branches in the overlap graph
//============================================================================

#include "OverlapGraph.h"

/*
IMPORTANT:
When using initial paired-end input reads as single-end sequences, the input
must be ordered singles-paired1-paired2 with corresponding integer read IDs from
0 to #SE+2*#PE. The number of SE and PE input reads must be specified when
invoking readBasedBranchReduction, as well as the minimum evidence required for
keeping a branching edge, and the initial read file(s).

APPROACH:
Identify all branches in the graph and for every branch u->(v_0,...,v_k):
1. Compute a list of all FIRST difference positions between any pair of
    branch sequences --> diff_list
2. For all i, find all common subreads between u and v_i, where PE read IDs from
    the same fragment are considered identical (i.e. modulo #PE). For all common
    subreads, check if the corresponding sequence is identical to the contig
    sequence at all positions of diff_list. If so, add to the "evidence set" of
    edge i.
3. Compute the shared evidence set, which is the intersection of evidence sets
    of all edges.
4. For all i, compute the relative evidence set by removing all initial reads
    appearing in the shared evidence set
5. If the maximal relative evidence is greater that the given threshold, remove
    all branching edges that do not have sufficient evidence.
*/

void OverlapGraph::readBasedBranchReduction(int SE_count, int PE_count, int min_evidence,
        std::string SE_file, std::string PE_file1, std::string PE_file2) {
    // adjust program settings to read from original input files
    ProgramSettings original_input = program_settings;
    original_input.singles_file = SE_file;
    original_input.paired1_file = PE_file1;
    original_input.paired2_file = PE_file2;
    // now read and store original fastq file(s)
    FastqStorage original_fastq = FastqStorage(original_input);

    // find all branches in the graph

        // find common subreads from originals dict

            // check all diff_list positions
                bool result = checkReadEvidence(contig, sequence, pos);

        // build shared evidence set

        // compute relative evidence per edge

        // remove edges with insufficient evidence 
}

void OverlapGraph::buildDiffList() {}


bool OverlapGraph::checkReadEvidence(std::string contig, std::string read, int position) {
    bool true_evidence = true;

    return true_evidence;
}
