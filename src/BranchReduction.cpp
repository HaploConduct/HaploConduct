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

void BranchReduction::readBasedBranchReduction() {
    std::cout << "readBasedBranchReduction" << std::endl;
    std::vector< std::list< node_id_t > > sorted_adj_in;
    std::vector< std::list< node_id_t > > sorted_adj_out;
    sorted_adj_in = overlap_graph->sortAdjLists(overlap_graph->adj_in);
    sorted_adj_out = overlap_graph->sortAdjOut(overlap_graph->adj_out);
    std::set< node_id_t > branch_in;
    std::set< node_id_t > branch_out;
    // find all branches in the graph and process one by one
    overlap_graph->findBranchfreeGraph(sorted_adj_in, sorted_adj_out, branch_in, branch_out);
    std::list< std::pair< node_id_t, node_id_t > > edges_to_remove1;
    std::list< std::pair< node_id_t, node_id_t > > edges_to_remove2;
    for (auto node : branch_in) {
        bool outbranch = false;
        findBranchingEvidence(node, sorted_adj_in.at(node), edges_to_remove1, outbranch);
    }
    for (auto node : branch_out) {
        bool outbranch = true;
        findBranchingEvidence(node, sorted_adj_out.at(node), edges_to_remove2, outbranch);
    }
    // build a set of edges to remove to avoid duplicates
    std::set< std::pair< node_id_t, node_id_t > > edges_to_remove (edges_to_remove2.begin(), edges_to_remove2.end());
    for (auto node_pair : edges_to_remove1) {
        edges_to_remove.insert(std::make_pair(node_pair.second, node_pair.first));
    }
    // now remove all selected edges from overlap graph
    for (auto node_pair : edges_to_remove) {
        std::cout << "removing " << node_pair.first << " " << node_pair.second << std::endl;
        Edge edge = overlap_graph->removeEdge(node_pair.first, node_pair.second);
        overlap_graph->branching_edges.push_back(edge);
    }
    std::cout << edges_to_remove1.size() + edges_to_remove2.size() << " edges removed" << std::endl;
}

void BranchReduction::findBranchingEvidence(node_id_t node1, std::list< node_id_t > neighbors,
        std::list< std::pair< node_id_t, node_id_t > > & edges_to_remove, bool outbranch) {
    assert (neighbors.size() > 1);
    // build list of difference positions (first difference for every pair)
    std::vector< std::string > sequence_vec;
    std::vector< int > startpos_vec;
    std::list< int > diff_list;
    std::vector< std::pair< unsigned int, unsigned int > > missing_edges;
    std::vector< node_id_t > neighbors_vec(neighbors.begin(), neighbors.end());
    if (outbranch) {
        diff_list = buildDiffListOut(node1, neighbors_vec, sequence_vec, startpos_vec, missing_edges);
    }
    else {
        diff_list = buildDiffListIn(node1, neighbors_vec, sequence_vec, startpos_vec, missing_edges);
    }
    // build evidence list of subread IDs per neighbor
    std::unordered_map< read_id_t, OriginalIndex > subreads1 = overlap_graph->original_ID_dict.at(node1);
    std::vector< std::vector< read_id_t > > evidence_per_neighbor;
    std::vector< std::string >::const_iterator seq_it = sequence_vec.begin();
    std::vector< int >::const_iterator startpos_it = startpos_vec.begin();
    for (auto node2 : neighbors) {
        // find common subreads from originals dict
        std::vector< read_id_t > evidence_list;
        std::string contig = *seq_it;
        int startpos = *startpos_it;
        seq_it++;
        startpos_it++;
        // Read* read = fastq_storage->get_read(node2);
        // std::string contig;
        // if (read->is_paired()) {
        //     std::cout << "TODO: resolve branches with PE reads as well" << std::endl;
        //     continue;
        // }
        // else {
        //     contig = read->get_seq(0);
        // }
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
            // first check the subread itself
            bool result1 = false;
            if (common_subread != subreads1.end()) {
                Read* original_read = original_fastq->get_read(subread_id);
                std::string sequence = original_read->get_seq(0);
                int index = subread.second.index1;
                result1 = checkReadEvidence(contig, startpos, sequence, index, diff_list);
            }
            if (result1) {
                evidence_list.push_back(subread_id);
            }
            // now also try its mate
            bool result2 = false;
            if (common_PE != subreads1.end()) {
                Read* original_read = original_fastq->get_read(subread_id_PE);
                std::string sequence = original_read->get_seq(0);
                int index = subread.second.index1;
                result2 = checkReadEvidence(contig, startpos, sequence, index, diff_list);
            }
            if (result2) {
                read_id_t joint_id = program_settings.original_readcount + std::min(subread_id, subread_id_PE);
                evidence_list.push_back(joint_id);
            }
        }
        std::sort(evidence_list.begin(), evidence_list.end());
        evidence_per_neighbor.push_back(evidence_list);
    }
    // take care of 'missing edges' by removing all edges from 'missing branches'
    if (!diff_list.empty()) {
        for (auto index_pair : missing_edges) {
            // remove evidence from node2
            unsigned int idx = index_pair.second;
            evidence_per_neighbor.at(idx).clear();
            // identify newly created branches and remove corresponding edges;
            // these branches will be resolved at a next iteration of ViralQuasispecies
            node_id_t outnode = neighbors_vec.at(index_pair.first);
            node_id_t innode = neighbors_vec.at(index_pair.second);
            for (auto out_edge : overlap_graph->adj_out.at(outnode)) {
                if (outbranch) {
                    edges_to_remove.push_back(std::make_pair(outnode, out_edge.get_vertex(2)));
                }
                else { // add reverse edge because it will be reversed again later
                    edges_to_remove.push_back(std::make_pair(out_edge.get_vertex(2), outnode));
                }
            }
            for (auto in_edge : overlap_graph->adj_in.at(innode)) {
                if (outbranch) {
                    edges_to_remove.push_back(std::make_pair(in_edge, innode));
                }
                else { // add reverse edge because it will be reversed again later
                    edges_to_remove.push_back(std::make_pair(innode, in_edge));
                }
            }
        }
    }
    // compute effective evidence per neighbor: remove evidence that is shared
    // between all neighbors and count evidence load
    std::vector< int > effective_evidence_counts;
    std::vector< unsigned int > evidence_indexes;
    std::vector< int > evidence_status;
    std::cout << "initial_evidence_counts" << std::endl;
    for (unsigned int idx = 0; idx < neighbors.size(); idx++) {
        std::cout << evidence_per_neighbor.at(idx).size() << std::endl;
        effective_evidence_counts.push_back(0);
        evidence_indexes.push_back(0);
        if (evidence_per_neighbor.at(idx).empty()) {
            evidence_status.push_back(0);
        }
        else {
            evidence_status.push_back(1);
        }
    }
    // filter evidence such that we only count UNIQUE read support
    while (*std::max_element(evidence_status.begin(), evidence_status.end()) == 1) {
        // get maximum value over all evidence iterators
        std::vector< read_id_t > current_evidence;
        for (unsigned int idx = 0; idx < neighbors.size(); idx++) {
            if (evidence_status.at(idx) == 1) {
                unsigned int ev_idx = evidence_indexes.at(idx);
                current_evidence.push_back(evidence_per_neighbor.at(idx).at(ev_idx));
            }
        }
        std::sort(current_evidence.begin(), current_evidence.end());
        assert (!current_evidence.empty());
        read_id_t current_min = current_evidence.front();
        read_id_t current_max = current_evidence.back();
        std::cout << "current_max " << current_max << " current_min " << current_min << std::endl;
        assert (current_max < 2*program_settings.original_readcount);
        assert (current_min < 2*program_settings.original_readcount);
        bool unique_min;
        if (diff_list.empty()) {
            unique_min = true; // no difference positions so keep all evidence
        }
        else if (current_evidence.size() == 1) {
            unique_min = true; // only evidence remaining so definitely unique
        }
        else if (current_min < current_evidence.at(1)) {
            unique_min = true;
        }
        else {
            unique_min = false;
        }
        for (unsigned int idx = 0; idx < neighbors.size(); idx++) {
            if (evidence_status.at(idx) == 1) {
                unsigned int ev_idx = evidence_indexes.at(idx);
                if (evidence_per_neighbor.at(idx).at(ev_idx) == current_min) {
                    if (unique_min) { // keep evidence
                        effective_evidence_counts.at(idx)++;
                    }
                    std::cout << "before: " << evidence_indexes.at(idx) << std::endl;
                    evidence_indexes.at(idx)++;
                    std::cout << "after: " << evidence_indexes.at(idx) << std::endl;
                    if (evidence_indexes.at(idx) >= evidence_per_neighbor.at(idx).size()) {
                        // reached end of evidence list, mark node as finished
                        evidence_status.at(idx) = 0;
                    }
                }
            }
        }
    }
    std::cout << "effective_evidence_counts" << std::endl;
    for (auto count : effective_evidence_counts) {
        std::cout << count << std::endl;
    }
    // finally analyze evidence and reduce branches
    assert (neighbors.size() == effective_evidence_counts.size());
    //int largest_evidence = *std::max_element(effective_evidence_counts.begin(), effective_evidence_counts.end());
    //if (largest_evidence >= min_evidence) {}
    // remove edges with insufficient evidence
    std::list< node_id_t >::const_iterator node2_it = neighbors.begin();
    for (unsigned int i=0; i < neighbors.size(); i++) {
        if (effective_evidence_counts.at(i) < min_evidence) {
            edges_to_remove.push_back(std::make_pair(node1, *node2_it));
        }
        node2_it++;
    }
}

std::list< int > BranchReduction::buildDiffListOut(node_id_t node1,
        std::vector< node_id_t > neighbors,
        std::vector< std::string > & sequence_vec,
        std::vector< int > & startpos_vec,
        std::vector< std::pair< node_id_t, node_id_t > > & missing_edges) {
    // build a list of all FIRST difference positions between any pair of branch sequences
    std::cout << "buildDiffListOut" << std::endl;
    std::list< int > diff_list;
    for (auto node : neighbors) {
        Edge* edge = overlap_graph->getEdgeInfo(node1, node, /*reverse_allowed*/ false);
        int pos = edge->get_pos(1);
        Read* read = edge->get_read(2);
        assert (!read->is_paired());
        std::string sequence = read->get_seq(0);
        sequence_vec.push_back(sequence);
        startpos_vec.push_back(pos);
    }
    // do all pairwise sequence comparisons until a difference is found;
    /* note: could be done more efficiently by evaluating all sequences in one
        for loop, but this would take extensive bookkeeping. Since the number
        of neighbors at the same branch is not expected to be big, all pairwise
        comparisons should be fine computationally. */
    node_id_t node_i, node_j;
    std::string seq_i, seq_j;
    std::string subseq_i, subseq_j;
    int pos_i, pos_j;
    int relative_pos;
    std::vector< int > diff_pos;
    int len;
    int startpos;
    for (unsigned int i=0; i < neighbors.size(); i++) {
        for (unsigned int j=i+1; j < neighbors.size(); j++) {
            node_i = neighbors.at(i);
            node_j = neighbors.at(j);
            seq_i = sequence_vec.at(i);
            seq_j = sequence_vec.at(j);
            pos_i = startpos_vec.at(i);
            pos_j = startpos_vec.at(j);
            assert (node_i != node_j);
            if (pos_i < pos_j) {
                relative_pos = pos_j - pos_i;
                if (relative_pos >= int(seq_i.size())) {
                    continue;
                }
                len = std::min(seq_i.size()-relative_pos, seq_j.size());
                subseq_i = seq_i.substr(relative_pos, len);
                subseq_j = seq_j.substr(0, len);
                diff_pos = findDiffPos(subseq_i, subseq_j);
                startpos = pos_j;
            }
            else {
                relative_pos = pos_i - pos_j;
                if (relative_pos >= int(seq_j.size())) {
                    continue;
                }
                len = std::min(seq_j.size()-relative_pos, seq_i.size());
                subseq_i = seq_i.substr(0, len);
                subseq_j = seq_j.substr(relative_pos, len);
                diff_pos = findDiffPos(subseq_i, subseq_j);
                startpos = pos_i;
            }
            assert (len > 0);
            for (auto pos : diff_pos) {
                diff_list.push_back(pos + startpos);
            }
            if (diff_pos.empty()) {
                // identical overlap -> add corresponding edge
                if (pos_i == pos_j && node_i < node_j) {
                    missing_edges.push_back(std::make_pair(i, j));
                }
                else if (pos_i == pos_j && node_i > node_j) {
                    missing_edges.push_back(std::make_pair(j, i));
                }
                else if (pos_i < pos_j) {
                    missing_edges.push_back(std::make_pair(i, j));
                }
                else {
                    missing_edges.push_back(std::make_pair(j, i));
                }
            }
        }
    }
    // remove duplicate entries
    diff_list.sort();
    diff_list.unique();
    return diff_list;
}

std::list< int > BranchReduction::buildDiffListIn(node_id_t node1,
        std::vector< node_id_t > neighbors,
        std::vector< std::string > & sequence_vec,
        std::vector< int > & startpos_vec,
        std::vector< std::pair< node_id_t, node_id_t > > & missing_edges) {
    // build a list of all FIRST difference positions between any pair of branch sequences
    std::cout << "buildDiffListIn" << std::endl;
    std::list< int > diff_list;
    std::vector< int > pos_vec;
    for (auto node : neighbors) {
        Edge* edge = overlap_graph->getEdgeInfo(node, node1, /*reverse_allowed*/ false);
        int pos = edge->get_pos(1);
        Read* read = edge->get_read(1);
        assert (!read->is_paired());
        std::string sequence = read->get_seq(0);
        sequence_vec.push_back(sequence);
        pos_vec.push_back(pos);
    }
    /* since we are dealing wint an in-branch, we need to infer the startpos_vec entries
        once we have those positions, we can proceed similar to the out-branch case,
        except for the diff_pos computation where we need to reverse the sequences
        first */
    int max_pos = *std::max_element(pos_vec.begin(), pos_vec.end());
    for (auto pos : pos_vec) {
        startpos_vec.push_back(max_pos - pos);
    }
    // do all pairwise sequence comparisons until a difference is found;
    /* note: could be done more efficiently by evaluating all sequences in one
        for loop, but this would take extensive bookkeeping. Since the number
        of neighbors at the same branch is not expected to be big, all pairwise
        comparisons should be fine computationally. */
    node_id_t node_i, node_j;
    std::string seq_i, seq_j;
    std::string subseq_i, subseq_j;
    int pos_i, pos_j;
    int relative_pos;
    std::vector< int > diff_pos;
    int len;
    int startpos;
    for (unsigned int i=0; i < neighbors.size(); i++) {
        for (unsigned int j=i+1; j < neighbors.size(); j++) {
            node_i = neighbors.at(i);
            node_j = neighbors.at(j);
            seq_i = sequence_vec.at(i);
            seq_j = sequence_vec.at(j);
            pos_i = startpos_vec.at(i);
            pos_j = startpos_vec.at(j);
            if (pos_i < pos_j) {
                relative_pos = pos_j - pos_i;
                if (relative_pos >= int(seq_i.size())) {
                    continue;
                }
                len = std::min(seq_i.size()-relative_pos, seq_j.size());
                subseq_i = seq_i.substr(relative_pos, len);
                subseq_j = seq_j.substr(0, len);
                std::reverse( subseq_i.begin(), subseq_i.end() );
                std::reverse( subseq_j.begin(), subseq_j.end() );
                diff_pos = findDiffPos(subseq_i, subseq_j);
                startpos = pos_j;
            }
            else {
                relative_pos = pos_i - pos_j;
                if (relative_pos >= int(seq_j.size())) {
                    continue;
                }
                len = std::min(seq_j.size()-relative_pos, seq_i.size());
                subseq_i = seq_i.substr(0, len);
                subseq_j = seq_j.substr(relative_pos, len);
                std::reverse( subseq_i.begin(), subseq_i.end() );
                std::reverse( subseq_j.begin(), subseq_j.end() );
                diff_pos = findDiffPos(subseq_i, subseq_j);
                startpos = pos_i;
            }
            assert (len > 0);
            for (auto pos : diff_pos) {
                diff_list.push_back(len - pos + startpos);
            }
            if (diff_pos.empty()) {
                // identical overlap -> add corresponding edge
                if (pos_i == pos_j && node_i < node_j) {
                    missing_edges.push_back(std::make_pair(i, j));
                }
                else if (pos_i == pos_j && node_i > node_j) {
                    missing_edges.push_back(std::make_pair(j, i));
                }
                else if (pos_i < pos_j) {
                    missing_edges.push_back(std::make_pair(i, j));
                }
                else {
                    missing_edges.push_back(std::make_pair(j, i));
                }
            }
        }
    }
    // remove duplicate entries
    diff_list.sort();
    diff_list.unique();
    return diff_list;
}



std::vector< int > BranchReduction::findDiffPos(std::string seq1, std::string seq2) {
    // given two input strings, find the first position where they disagree
    assert (seq1.size() == seq2.size());
    std::string::const_iterator it1 = seq1.begin();
    std::string::const_iterator it2 = seq2.begin();
    std::vector< int > diff_pos;
    int pos = 0;
    while (it1 != seq1.end()) {
        if (*it1 != *it2) {
            diff_pos.push_back(pos);
            if (diff_pos.size() == 3) {
                // store at most 3 positions per sequence pair
                return diff_pos;
            }
        }
        it1++;
        it2++;
        pos++;
    }
    if (diff_pos.empty()) {
        // no difference found -> return sequence length
        std::cout << "no difference found, len = " << seq1.size() << std::endl;
    }
    return diff_pos;
}


bool BranchReduction::checkReadEvidence(std::string contig, int startpos, std::string read, int index, std::list< int > diff_list) {
    // check if subread agrees with contig on diff_list positions
    bool true_evidence = true;
    int read_start = startpos + index;
    int read_end = read_start + read.size();
    int contig_start = startpos;
    int contig_end = startpos + contig.size();
    for (auto diff_pos : diff_list) {
        if (diff_pos < read_start || diff_pos >= read_end) {
            // read does not overlap this diff_pos
            continue;
        }
        else if (diff_pos < contig_start || diff_pos >= contig_end) {
            // contig does not overlap this diff_pos
            continue;
        }
        // check if read agrees with contig base
        std::cout << diff_pos << " " << startpos << " " << index << " " << contig.size() << " " << read.size() << std::endl;
        if (read.at(diff_pos - read_start) != contig.at(diff_pos - contig_start)) {
            // disagreement --> false evidence
            true_evidence = false;
            break; // no need to continue checking
        }
    }
    return true_evidence;
}
