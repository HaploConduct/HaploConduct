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
    std::vector< std::list< node_id_t > > sorted_adj_in;
    std::vector< std::list< node_id_t > > sorted_adj_out;
    sorted_adj_in = sortAdjLists(adj_in);
    sorted_adj_out = sortAdjOut(adj_out);
    std::set< node_id_t > branch_in;
    std::set< node_id_t > branch_out;
    findBranchfreeGraph(sorted_adj_in, sorted_adj_out, branch_in, branch_out);
    std::list< std::pair< node_id_t, node_id_t > > edges_to_remove;
    for (auto node : branch_in) {
        findBranchingEvidence(node, sorted_adj_in.at(node), edges_to_remove, SE_count, PE_count, min_evidence);
    }
    for (auto node : branch_out) {
        findBranchingEvidence(node, sorted_adj_out.at(node), edges_to_remove, SE_count, PE_count, min_evidence);
    }
    // now remove all selected edges from overlap graph
    for (auto node_pair : edges_to_remove) {
        Edge edge = removeEdge(node_pair.first, node_pair.second);
        branching_edges.push_back(edge);
    }
}

void OverlapGraph::findBranchingEvidence(node_id_t node1, std::list< node_id_t > neighbors,
        std::list< std::pair< node_id_t, node_id_t > > & edges_to_remove,
        int SE_count, int PE_count, int min_evidence) {
    // build list of difference positions (first difference for every pair)
    std::list< int > diff_list = buildDiffList(node1, neighbors);
    std::unordered_map< read_id_t, OriginalIndex > subreads1 = original_ID_dict.at(node1);
    for (auto node2 : neighbors) {
        // find common subreads from originals dict
        std::unordered_map< read_id_t, OriginalIndex > subreads2 = original_ID_dict.at(node2);
        for (auto subread : subreads2) {
            int subread_id = subread.first;
            std::unordered_map< read_id_t, OriginalIndex >::const_iterator common_subread;
            common_subread = subreads1.find(subread_id);
            std::unordered_map< read_id_t, OriginalIndex >::const_iterator common_PE;
            if (subread_id >= SE_count + PE_count) {
                common_PE = subreads1.find(subread_id - PE_count);
            }
            else if (subread_id >= SE_count) {
                common_PE = subreads1.find(subread_id + PE_count);
            }
            else {
                common_PE = subreads1.end();
            }
            std::string contig, sequence;
            bool result1 = checkReadEvidence(contig, sequence, diff_list);
            bool result2 = checkReadEvidence(contig, sequence, diff_list);
            std::cout << result1 << " " << result2 << std::endl;
        }

        // build shared evidence set

        // compute relative evidence per edge

        // remove edges with insufficient evidence
    }
}

std::list< int > OverlapGraph::buildDiffList(node_id_t node1, std::list< node_id_t > neighbors) {
    std::list< int > diff_list;
    // do stuff
    diff_list.push_back(1);
    return diff_list;
}


bool OverlapGraph::checkReadEvidence(std::string contig, std::string read, std::list< int > diff_list) {
    // check if subread agrees with contig on diff_list positions
    bool true_evidence = true;
    // do stuff
    return true_evidence;
}
