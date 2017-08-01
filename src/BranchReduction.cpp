//============================================================================
// Name        : BranchReduction.cpp
// Author      : Jasmijn Baaijens
// Version     : 0.4.0
// License     : GNU GPL v3.0
// Project     : ViralQuasispecies
// Description : Use read evidence to reduce branches in the overlap graph
//============================================================================

#include "BranchReduction.h"

/*
IMPORTANT:
When using initial paired-end input reads as single-end sequences, the input
must be ordered singles-paired1-paired2 with corresponding integer read IDs from
0 to #SE+2*#PE. The number of SE and PE input reads must be specified when
invoking readBasedBranchReduction, as well as the minimum evidence required for
keeping a branching edge, and the initial read file(s).
Otherwise, make sure PE_count = 0.

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

void BranchReduction::readBasedBranchReduction(int SE_count, int PE_count, int min_evidence,
        std::string SE_file, std::string PE_file1, std::string PE_file2) {
    // std::cout << "readBasedBranchReduction" << std::endl;
    std::vector< std::list< node_id_t > > sorted_adj_in;
    std::vector< std::list< node_id_t > > sorted_adj_out;
    std::vector< std::list< node_id_t > > test = overlap_graph->adj_in;
    sorted_adj_in = overlap_graph->sortAdjLists(test);
    sorted_adj_out = overlap_graph->sortAdjOut(overlap_graph->adj_out);
    std::set< node_id_t > branch_in;
    std::set< node_id_t > branch_out;
    // find all branches in the graph and process one by one
    overlap_graph->findBranchfreeGraph(sorted_adj_in, sorted_adj_out, branch_in, branch_out);
    std::list< std::pair< node_id_t, node_id_t > > edges_to_remove1;
    std::list< std::pair< node_id_t, node_id_t > > edges_to_remove2;
    for (auto node : branch_in) {
        findBranchingEvidence(node, sorted_adj_in.at(node), edges_to_remove1, SE_count, PE_count, min_evidence);
    }
    for (auto node : branch_out) {
        findBranchingEvidence(node, sorted_adj_out.at(node), edges_to_remove2, SE_count, PE_count, min_evidence);
    }
    // now remove all selected edges from overlap graph
    for (auto node_pair : edges_to_remove1) {
        Edge edge = overlap_graph->removeEdge(node_pair.second, node_pair.first);
//        branching_edges.push_back(edge);
    }
    for (auto node_pair : edges_to_remove2) {
        Edge edge = overlap_graph->removeEdge(node_pair.first, node_pair.second);
//        branching_edges.push_back(edge);
    }
}

void BranchReduction::findBranchingEvidence(node_id_t node1, std::list< node_id_t > neighbors,
        std::list< std::pair< node_id_t, node_id_t > > & edges_to_remove,
        unsigned int SE_count, unsigned int PE_count, int min_evidence) {
    // build list of difference positions (first difference for every pair)
    std::list< int > diff_list = buildDiffListOut(node1, neighbors);
    std::unordered_map< read_id_t, OriginalIndex > subreads1 = overlap_graph->original_ID_dict.at(node1);
    // build evidence list of subread IDs per neighbor
    std::vector< std::vector< read_id_t > > evidence_per_neighbor;
    for (auto node2 : neighbors) {
        // find common subreads from originals dict
        std::vector< read_id_t > evidence_list;
        Read* read = fastq_storage->get_read(node2);
        std::string contig;
        if (read->is_paired()) {
            std::cout << "TODO: resolve branches with PE reads as well" << std::endl;
            continue;
        }
        else {
            contig = read->get_seq(0);
        }
        std::unordered_map< read_id_t, OriginalIndex > subreads2 = overlap_graph->original_ID_dict.at(node2);
        for (auto subread : subreads2) {
            read_id_t subread_id = subread.first;
            std::unordered_map< read_id_t, OriginalIndex >::const_iterator common_subread;
            common_subread = subreads1.find(subread_id);
            std::unordered_map< read_id_t, OriginalIndex >::const_iterator common_PE;
            read_id_t subread_id_PE;
            if (subread_id >= SE_count + PE_count) {
                // try its /1 mate
                common_PE = subreads1.find(subread_id - PE_count);
                subread_id_PE = subread_id - PE_count;
            }
            else if (subread_id >= SE_count) {
                // try its /2 mate
                common_PE = subreads1.find(subread_id + PE_count);
                subread_id_PE = subread_id + PE_count;
            }
            else {
                // single-end subread so no mate to try
                common_PE = subreads1.end();
                subread_id_PE = 0; // dummy
            }
            bool result1 = false;
            bool result2 = false;
            if (common_subread != subreads1.end()) {
                Read* original_read = original_fastq->get_read(subread_id);
                std::string sequence = original_read->get_seq(0);
                int index = subread.second.index1;
                result1 = checkReadEvidence(contig, sequence, index, diff_list);
            }
            if (common_PE != subreads1.end()) {
                Read* original_read = original_fastq->get_read(subread_id_PE);
                std::string sequence = original_read->get_seq(0);
                int index = subread.second.index1;
                result2 = checkReadEvidence(contig, sequence, index, diff_list);
            }
            if (result1 || result2) {
                evidence_list.push_back(subread_id);
            }
        }
        std::sort(evidence_list.begin(), evidence_list.end());
        evidence_per_neighbor.push_back(evidence_list);
    }
    // compute effective evidence per neighbor: remove evidence that is shared
    // between all neighbors and count evidence load
    std::vector< int > effective_evidence_counts;
    std::vector< std::vector< read_id_t >::const_iterator > evidence_iterators;
    std::vector< int > evidence_status;
    std::cout << "initial_evidence_counts" << std::endl;
    for (auto evidence_list : evidence_per_neighbor) {
        std::cout << evidence_list.size() << std::endl;
        effective_evidence_counts.push_back(0);
        evidence_iterators.push_back(evidence_list.begin());
        evidence_status.push_back(1);
    }
    while (*std::max_element(evidence_status.begin(), evidence_status.end()) == 1) {
        // get maximum value over all evidence iterators
        read_id_t current_max = **std::max_element(evidence_iterators.begin(), evidence_iterators.end(),
            [](const std::vector< read_id_t >::const_iterator & a,
                const std::vector< read_id_t >::const_iterator & b) -> bool
            {
                return *a < *b;
            });
        unsigned int idx = 0;
        for (auto it : evidence_iterators) {
            if (*it < current_max) {
                // non-trivial evidence, i.e. not shared between all neighbors
                effective_evidence_counts.at(idx)++;
                it++;
                if (it == evidence_per_neighbor.at(idx).end()) {
                    // reached end of evidence list, mark node as finished
                    evidence_status.at(idx) = 0;
                }
            }
            idx++;
        }
    }
    std::cout << "effective_evidence_counts" << std::endl;
    for (auto count : effective_evidence_counts) {
        std::cout << count << std::endl;
    }
    // finally analyze evidence and reduce branches
    assert (neighbors.size() == effective_evidence_counts.size());
    int largest_evidence = *std::max_element(effective_evidence_counts.begin(), effective_evidence_counts.end());
    if (largest_evidence >= min_evidence) {
        // remove edges with insufficient evidence
        std::list< node_id_t >::const_iterator node2_it = neighbors.begin();
        for (unsigned int i=0; i < neighbors.size(); i++) {
            if (effective_evidence_counts.at(i) < min_evidence) {
                edges_to_remove.push_back(std::make_pair(node1, *node2_it));
            }
            node2_it++;
        }
    }
}

std::list< int > BranchReduction::buildDiffListOut(node_id_t node1, std::list< node_id_t > neighbors) {
    // build a list of all FIRST difference positions between any pair of branch sequences
    std::list< int > diff_list;
    std::vector< std::string > sequence_vec;
    std::vector< int > startpos_vec;
    for (auto node : neighbors) {
        Edge* edge = overlap_graph->getEdgeInfo(node1, node, /*reverse_allowed*/ false);
        int pos = edge->get_pos(1);
        Read* read = edge->get_read(2);
        std::string sequence = read->get_seq(0);
        sequence_vec.push_back(sequence);
        startpos_vec.push_back(pos);
    }
    // do all pairwise sequence comparisons until a difference is found;
    /* note: could be done more efficiently by evaluating all sequences in one
        for loop, but this would take extensive bookkeeping. Since the number
        of neighbors at the same branch is not expected to be big, all pairwise
        comparisons should be fine computationally. */
    std::string seq_i, seq_j;
    int pos_i, pos_j;
    int relative_pos;
    int diff_pos;
    int len;
    for (unsigned int i=0; i < neighbors.size(); i++) {
        for (unsigned int j=i+1; j < neighbors.size(); j++) {
            seq_i = sequence_vec.at(i);
            seq_j = sequence_vec.at(j);
            pos_i = startpos_vec.at(i);
            pos_j = startpos_vec.at(j);
            if (pos_i < pos_j) {
                relative_pos = pos_j - pos_i;
                len = std::min(seq_i.size()-relative_pos, seq_j.size());
                diff_pos = findDiffPos(seq_i.substr(relative_pos, len), seq_j.substr(0, len));
                diff_list.push_back(diff_pos + pos_j);
            }
            else {
                relative_pos = pos_i - pos_j;
                len = std::min(seq_j.size()-relative_pos, seq_i.size());
                diff_pos = findDiffPos(seq_i.substr(0, len), seq_j.substr(relative_pos, len));
                diff_list.push_back(diff_pos + pos_i);
            }
            assert (diff_pos < len);
        }
    }
    // remove duplicate entries
    diff_list.sort();
    diff_list.unique();
    return diff_list;
}


int BranchReduction::findDiffPos(std::string seq1, std::string seq2) {
    // given two input strings, find the first position where they disagree
    assert (seq1.size() == seq2.size());
    std::string::const_iterator it1 = seq1.begin();
    std::string::const_iterator it2 = seq2.begin();
    int diff_pos = 0;
    while (it1 != seq1.end()) {
        if (*it1 != *it2) {
            return diff_pos;
        }
        it1++;
        it2++;
        diff_pos++;
    }
    // no difference found -> return sequence length
    return diff_pos;
}


bool BranchReduction::checkReadEvidence(std::string contig, std::string read, int index, std::list< int > diff_list) {
    // check if subread agrees with contig on diff_list positions
    bool true_evidence = true;
    // do stuff
    return true_evidence;
}
