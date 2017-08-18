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
    std::list< Edge > missing_edges;
    std::vector< std::list< node_id_t > > final_branch_in (overlap_graph->getVertexCount(), std::list< node_id_t >());
    for (auto node : branch_in) {
        bool outbranch = false;
        std::list< node_id_t > branch = findBranchingEvidence(node, sorted_adj_in.at(node),
            missing_edges, outbranch);
        if (!branch.empty()) {
            final_branch_in.at(node) = branch;
        }
    }
    std::vector< std::list< node_id_t > > final_branch_out (overlap_graph->getVertexCount(), std::list< node_id_t >());
    for (auto node : branch_out) {
        bool outbranch = true;
        std::list< node_id_t > branch = findBranchingEvidence(node, sorted_adj_out.at(node),
            missing_edges, outbranch);
        if (!branch.empty()) {
            final_branch_out.at(node) = branch;
        }
    }
    // add missing edges to vector of removed branching edges, to ensure that
    // they will be found in the next iteration
    for (auto edge : missing_edges) {
        overlap_graph->branching_edges.push_back(edge);
    }
    // find branching components; components with false branches (i.e. missing edges)
    // are skipped automatically and all edges removed
    std::list< node_pair_t > edges_to_remove;
    findBranchingComponents(final_branch_in, final_branch_out, edges_to_remove);
    std::cout << "Final evidence_per_edge.size() = " << evidence_per_edge.size() << std::endl;
    // select UNIQUE evidence and add edges with insufficient evidence to edges_to_remove
    for (auto component : branching_components) {
        countUniqueEvidence(component, edges_to_remove);
    }
    // remove duplicates from edges_to_remove
    edges_to_remove.sort();
    edges_to_remove.unique();
    // now remove all selected edges from overlap graph
    for (auto node_pair : edges_to_remove) {
//        std::cout << "removing " << node_pair.first << " " << node_pair.second << std::endl;
        if (overlap_graph->checkEdge(node_pair.first, node_pair.second, false) < 0) {
            std::cout << "edge not found, reverse: " << overlap_graph->checkEdge(node_pair.second, node_pair.first, false) << std::endl;
            exit(1);
        }
        Edge edge = overlap_graph->removeEdge(node_pair.first, node_pair.second);
        overlap_graph->branching_edges.push_back(edge);
    }
    std::cout << edges_to_remove.size() << " edges removed" << std::endl;
}

std::list< node_id_t > BranchReduction::findBranchingEvidence(node_id_t node1, std::list< node_id_t > neighbors,
        std::list< Edge > & missing_edges, bool outbranch) {
//    std::cout << "findBranchingEvidence" << std::endl;
    assert (neighbors.size() > 1);
    std::list< node_id_t > final_branch (neighbors.begin(), neighbors.end());
    final_branch.push_front(node1);
    // build list of difference positions (first difference for every pair)
    std::vector< std::string > sequence_vec;
    std::vector< int > startpos_vec;
    std::list< int > diff_list;
    std::vector< node_pair_t > missing_inclusion_edges;
    std::vector< node_id_t > neighbors_vec(neighbors.begin(), neighbors.end());
    if (outbranch) {
        diff_list = buildDiffListOut(node1, neighbors_vec, sequence_vec, startpos_vec, missing_inclusion_edges, missing_edges);
    }
    else {
        diff_list = buildDiffListIn(node1, neighbors_vec, sequence_vec, startpos_vec, missing_edges);
    }
    // build evidence list of subread IDs per neighbor
    std::unordered_map< read_id_t, OriginalIndex > subreads1 = overlap_graph->original_ID_dict.at(node1);
    std::unordered_map< node_id_t, std::list< read_id_t > > evidence_per_neighbor;
    std::vector< std::string >::const_iterator seq_it = sequence_vec.begin();
    std::vector< int >::const_iterator startpos_it = startpos_vec.begin();
    for (auto node2 : neighbors) {
        // find common subreads from originals dict
        std::list< read_id_t > evidence_list;
        std::string contig = *seq_it;
        int startpos = *startpos_it;
        seq_it++;
        startpos_it++;
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
                std::cout << "paired evidence found!" << std::endl;
            }
        }
        evidence_list.sort();
        evidence_list.unique();
        evidence_per_neighbor.insert(std::make_pair(node2, evidence_list));
    }
    // take care of 'missing edges' due to inclusions
    for (auto node_pair : missing_inclusion_edges) {
        evidence_per_neighbor.at(node_pair.first).clear();
        if (neighbors.size() == 2) { // no actual branch
            final_branch.clear();
        }
        else { // real branch, only remove inclusion node
            final_branch.remove(node_pair.first);
        }
        // // identify newly created branches and remove corresponding edges;
        // // these branches will be resolved at a next iteration of ViralQuasispecies
        // node_id_t outnode = node_pair.first;
        // node_id_t innode = node_pair.second;
        // for (auto out_edge : overlap_graph->adj_out.at(outnode)) {
        //     edges_to_remove.push_back(std::make_pair(outnode, out_edge.get_vertex(2)));
        // }
        // for (auto in_edge : overlap_graph->adj_in.at(innode)) {
        //     edges_to_remove.push_back(std::make_pair(in_edge, innode));
        // }
    }
    // store evidence in evidence_per_edge vector (private class member)
    std::list< node_id_t >::const_iterator branch_it = final_branch.begin();
    branch_it++; // skip first entry, this contains the branching node
    for (auto neighbor : neighbors) {
        if (branch_it != final_branch.end() && neighbor == *branch_it) {
            // neighbor also in final_branch so we store its evidence
            std::unordered_map< safe_edge_count_t, std::list< read_id_t > >::iterator ev_it;
            std::list< read_id_t > current_evidence = evidence_per_neighbor.at(neighbor);
            safe_edge_count_t index;
            if (outbranch) {
                index = edgeToEvidenceIndex(node1, neighbor);
            }
            else {
                index = edgeToEvidenceIndex(neighbor, node1);
            }
            ev_it = evidence_per_edge.find(index);
            if (ev_it != evidence_per_edge.end()) {
                // intersect evidence sets
                std::list< read_id_t > existing_evidence = ev_it->second;
                std::list< read_id_t >::iterator it2 = existing_evidence.begin();
                while (it2 != existing_evidence.end()) {
                    auto find_ev = std::find(current_evidence.begin(), current_evidence.end(), *it2);
                    if (find_ev == current_evidence.end()) {
                        it2 = ev_it->second.erase(it2);
                    }
                    else {
                        it2++;
                    }
                }
            }
            else {
                // insert new evidence set
                evidence_per_edge.insert(std::make_pair(index, current_evidence));
            }
            branch_it++;
        }
    }
    assert (branch_it == final_branch.end());
    return final_branch;
}

std::list< int > BranchReduction::buildDiffListOut(node_id_t node1,
        std::vector< node_id_t > neighbors,
        std::vector< std::string > & sequence_vec,
        std::vector< int > & startpos_vec,
        std::vector< node_pair_t > & missing_inclusion_edges,
        std::list< Edge > & missing_edges) {
    // build a list of all FIRST difference positions between any pair of branch sequences
//    std::cout << "buildDiffListOut" << std::endl;
    std::list< int > diff_list;
    std::vector< Edge* > edge_vec;
    for (auto node : neighbors) {
        Edge* edge = overlap_graph->getEdgeInfo(node1, node, /*reverse_allowed*/ false);
        int pos = edge->get_pos(1);
        Read* read = edge->get_read(2);
        assert (!read->is_paired());
        std::string sequence = read->get_seq(0);
        sequence_vec.push_back(sequence);
        startpos_vec.push_back(pos);
        edge_vec.push_back(edge);
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
                if (relative_pos > int(seq_i.size() - program_settings.min_overlap_len)) {
                    missing_inclusion_edges.push_back(std::make_pair(node_i, node_j));
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
                if (relative_pos > int(seq_j.size() - program_settings.min_overlap_len)) {
                    missing_inclusion_edges.push_back(std::make_pair(node_j, node_i));
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
                // identical overlap -> add corresponding edge; NOTE: this should hardly occur!
                double score = program_settings.edge_threshold; // overlap score
                int pos1 = relative_pos;
                int pos2 = 0; // no PE-overlaps allowed
                bool ori1; // 1 if NORMAL, 0 if REVERSE (vertex 1)
                bool ori2; // similar for vertex 2
                Read* read1; // pointer to read corresponding to vertex1
                Read* read2; // pointer to read corresponding to vertex2
                std::string ord = "-";
                node_id_t vertex1; // out-vertex
                node_id_t vertex2; // in-vertex
                int overlap_perc = (int)floor(100 * len / std::min(seq_i.size(), seq_j.size())) ;

                if (pos_i < pos_j || (pos_i == pos_j && node_i < node_j)) {
                    // node_i first
                    ori1 = edge_vec.at(i)->get_ori(2);
                    ori2 = edge_vec.at(j)->get_ori(2);
                    read1 = edge_vec.at(i)->get_read(2);
                    read2 = edge_vec.at(j)->get_read(2);
                    vertex1 = node_i;
                    vertex2 = node_j;
                }
                else {
                    // node_j first
                    ori1 = edge_vec.at(j)->get_ori(2);
                    ori2 = edge_vec.at(i)->get_ori(2);
                    read1 = edge_vec.at(j)->get_read(2);
                    read2 = edge_vec.at(i)->get_read(2);
                    vertex1 = node_j;
                    vertex2 = node_i;
                }
                Edge new_edge(score, pos1, pos2, ori1, ori2, ord, read1, read2);
                new_edge.set_vertices(vertex1, vertex2);
                new_edge.set_perc(overlap_perc);
                new_edge.set_len(len, 0);
                missing_edges.push_back(new_edge);
                // mark node as false out branch
                false_out_branches.insert(node1);
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
        std::list< Edge > & missing_edges) {
    // build a list of all FIRST difference positions between any pair of branch sequences
//    std::cout << "buildDiffListIn" << std::endl;
    std::list< int > diff_list;
    std::vector< int > pos_vec;
    std::vector< Edge* > edge_vec;
    for (auto node : neighbors) {
        Edge* edge = overlap_graph->getEdgeInfo(node, node1, /*reverse_allowed*/ false);
        int pos = edge->get_pos(1);
        Read* read = edge->get_read(1);
        assert (!read->is_paired());
        std::string sequence = read->get_seq(0);
        sequence_vec.push_back(sequence);
        pos_vec.push_back(pos);
        edge_vec.push_back(edge);
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
                assert (relative_pos <= int(seq_i.size() - program_settings.min_overlap_len)); // in-branches can't be the result of inclusions
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
                assert (relative_pos <= int(seq_j.size() - program_settings.min_overlap_len)); // in-branches can't be the result of inclusions
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
                // identical overlap -> add corresponding edge; NOTE: this should hardly occur!
                double score = program_settings.edge_threshold; // overlap score
                int pos1 = relative_pos;
                int pos2 = 0; // no PE-overlaps allowed
                bool ori1; // 1 if NORMAL, 0 if REVERSE (vertex 1)
                bool ori2; // similar for vertex 2
                Read* read1; // pointer to read corresponding to vertex1
                Read* read2; // pointer to read corresponding to vertex2
                std::string ord = "-";
                node_id_t vertex1; // out-vertex
                node_id_t vertex2; // in-vertex
                int overlap_perc = (int)floor(100 * len / std::min(seq_i.size(), seq_j.size()));

                if (pos_i < pos_j || (pos_i == pos_j && node_i < node_j)) {
                    // node_i first
                    ori1 = edge_vec.at(i)->get_ori(1);
                    ori2 = edge_vec.at(j)->get_ori(1);
                    read1 = edge_vec.at(i)->get_read(1);
                    read2 = edge_vec.at(j)->get_read(1);
                    vertex1 = node_i;
                    vertex2 = node_j;
                }
                else {
                    // node_j first
                    ori1 = edge_vec.at(j)->get_ori(1);
                    ori2 = edge_vec.at(i)->get_ori(1);
                    read1 = edge_vec.at(j)->get_read(1);
                    read2 = edge_vec.at(i)->get_read(1);
                    vertex1 = node_j;
                    vertex2 = node_i;
                }
                Edge new_edge(score, pos1, pos2, ori1, ori2, ord, read1, read2);
                new_edge.set_vertices(vertex1, vertex2);
                new_edge.set_perc(overlap_perc);
                new_edge.set_len(len, 0);
                missing_edges.push_back(new_edge);
                // also mark branch node as false
                false_in_branches.insert(node1);
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
            if (diff_pos.size() == 100) {
                // store at most 100 positions per sequence pair
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
    bool true_evidence = false;
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
//        std::cout << diff_pos << " " << startpos << " " << index << " " << contig.size() << " " << read.size() << std::endl;
        if (read.at(diff_pos - read_start) != contig.at(diff_pos - contig_start)) {
            // disagreement --> false evidence
            true_evidence = false;
            break; // no need to continue checking
        }
        else {
            true_evidence = true; // at least one diff_pos covered
        }
    }
    return true_evidence;
}

void BranchReduction::findBranchingComponents(std::vector< std::list< node_id_t > > final_branch_in,
    std::vector< std::list< node_id_t > > final_branch_out,
    std::list< node_pair_t > & edges_to_remove) {
    // find all branching components by following in-branches connected to out-branches
    // and vice versa
    std::cout << "findBranchingComponents" << std::endl;
    std::unordered_map< node_id_t, bool > visited_in_branches;
    std::unordered_map< node_id_t, bool > visited_out_branches;
    std::unordered_map< node_id_t, std::list< node_id_t > > branch_in_map;
    std::unordered_map< node_id_t, std::list< node_id_t > > branch_out_map;
    // store all branches in a dict for easy access
    for (auto branch : final_branch_in) {
        if (branch.empty()) {
            continue;
        }
        node_id_t node = branch.front();
        branch.pop_front();
        visited_in_branches.insert(std::make_pair(node, false));
        branch_in_map.insert(std::make_pair(node, branch));
    }
    for (auto branch : final_branch_out) {
        if (branch.empty()) {
            continue;
        }
        node_id_t node = branch.front();
        branch.pop_front();
        visited_out_branches.insert(std::make_pair(node, false));
        branch_out_map.insert(std::make_pair(node, branch));
    }
    // now build components
    for (auto branch : branch_in_map) {
        node_id_t node = branch.first;
        if (visited_in_branches.at(node)) {
            continue;
        }
        std::list< node_id_t > neighbors = branch.second;
        std::vector< node_pair_t > new_component;
        bool has_false_branch;
        if (false_in_branches.find(node) != false_in_branches.end()) {
            has_false_branch = true;
        }
        else {
            has_false_branch = false;
        }
        for (auto in_neighbor : neighbors) {
            assert (node != in_neighbor);
            node_pair_t node_pair = std::make_pair(in_neighbor, node);
            new_component.push_back(node_pair);
        }
        visited_in_branches.at(node) = true;
        extendComponentOut(new_component, neighbors, has_false_branch, visited_in_branches,
            visited_out_branches, branch_in_map, branch_out_map);
        std::sort(new_component.begin(), new_component.end(),
            [](node_pair_t pair1, node_pair_t pair2) {
                if (pair1.first != pair2.first) {
                    return pair1.first < pair2.first;
                }
                else {
                    return pair1.second < pair2.second;
                }
            });
        auto end_it = std::unique(new_component.begin(), new_component.end());
        new_component.resize(std::distance(new_component.begin(), end_it));
        if (has_false_branch) {
            // component contains a false branch due to a missing edge;
            // remove entire component and add missing edge at next iteration
            for (auto node_pair : new_component) {
                edges_to_remove.push_back(node_pair);
            }
        }
        else {
            branching_components.push_back(new_component);
        }
    }
    // process remaining out-branches; since we already did all in-branches, the
    // remaining components are trivial
    for (auto branch : branch_out_map) {
        node_id_t node = branch.first;
        if (visited_out_branches.at(node) || visited_in_branches.find(node) != visited_in_branches.end()) {
            continue;
        }
        std::list< node_id_t > neighbors = branch.second;
        std::vector< node_pair_t > new_component;
        for (auto out_neighbor : neighbors) {
            assert (node != out_neighbor);
            node_pair_t node_pair = std::make_pair(node, out_neighbor);
            new_component.push_back(node_pair);
        }
        if (false_out_branches.find(node) != false_out_branches.end()) {
            for (auto node_pair : new_component) {
                edges_to_remove.push_back(node_pair);
            }
        }
        else {
            branching_components.push_back(new_component);
        }
        visited_out_branches.at(node) = true;
    }
    return;
}

void BranchReduction::extendComponentOut(std::vector< node_pair_t > & component,
    std::list< node_id_t > neighbors, bool & has_false_branch,
    std::unordered_map< node_id_t, bool > & visited_in_branches,
    std::unordered_map< node_id_t, bool > & visited_out_branches,
    std::unordered_map< node_id_t, std::list< node_id_t > > & branch_in_map,
    std::unordered_map< node_id_t, std::list< node_id_t > > & branch_out_map) {
    // extend component iteratively
//    std::cout << "extendComponentOut" << std::endl;
    for (auto node : neighbors) {
        auto visited = visited_out_branches.find(node);
        if (visited == visited_out_branches.end() || visited->second == true) {
            continue;
        }
        if (false_out_branches.find(node) != false_out_branches.end()) {
            has_false_branch = true;
        }
        std::list< node_id_t > branch = branch_out_map.at(node);
        for (auto out_neighbor : branch) {
            assert (node != out_neighbor);
            node_pair_t node_pair = std::make_pair(node, out_neighbor);
            component.push_back(node_pair);
        }
        visited_out_branches.at(node) = true;
        extendComponentIn(component, branch, has_false_branch, visited_in_branches,
            visited_out_branches, branch_in_map, branch_out_map);
    }
    return;
}

void BranchReduction::extendComponentIn(std::vector< node_pair_t > & component,
    std::list< node_id_t > neighbors, bool & has_false_branch,
    std::unordered_map< node_id_t, bool > & visited_in_branches,
    std::unordered_map< node_id_t, bool > & visited_out_branches,
    std::unordered_map< node_id_t, std::list< node_id_t > > & branch_in_map,
    std::unordered_map< node_id_t, std::list< node_id_t > > & branch_out_map) {
    // extend component iteratively
//    std::cout << "extendComponentIn" << std::endl;
    for (auto node : neighbors) {
        auto visited = visited_in_branches.find(node);
        if (visited == visited_in_branches.end() || visited->second == true) {
            continue;
        }
        if (false_in_branches.find(node) != false_in_branches.end()) {
            has_false_branch = true;
        }
        std::list< node_id_t > branch = branch_in_map.at(node);
        for (auto in_neighbor : branch) {
            assert (node != in_neighbor);
            node_pair_t node_pair = std::make_pair(in_neighbor, node);
            component.push_back(node_pair);
        }
        visited_in_branches.at(node) = true;
        extendComponentOut(component, branch, has_false_branch, visited_in_branches,
            visited_out_branches, branch_in_map, branch_out_map);
    }
    return;
}

void BranchReduction::countUniqueEvidence(std::vector< node_pair_t > component,
        std::list< node_pair_t > & edges_to_remove) {
    // compute effective evidence per edge: remove evidence that is shared
    // between other edges and count evidence load
//    std::cout << "countUniqueEvidence" << std::endl;
    std::unordered_map< safe_edge_count_t, std::list< read_id_t > > unique_evidence_per_edge;
    std::unordered_map< safe_edge_count_t, unsigned int > index_to_mapID;
    std::vector< bool > evidence_status;
    unsigned int idx = 0;
    for (auto node_pair : component) {
        assert (node_pair.first != node_pair.second);
        safe_edge_count_t mapID = edgeToEvidenceIndex(node_pair.first, node_pair.second);
        index_to_mapID.insert(std::make_pair(idx, mapID));
        auto evidence = evidence_per_edge.find(mapID);
        if (evidence == evidence_per_edge.end()) {
            std::cout << "mapID not found for edge " << node_pair.first << " " << node_pair.second << std::endl;
        }
        else {
            if (evidence_per_edge.at(mapID).empty()) {
                evidence_status.push_back(0);
            }
            else {
                evidence_status.push_back(1);
            }
        }
        unique_evidence_per_edge.insert(std::make_pair(mapID, std::list< read_id_t >()));
        idx++;
    }
    unsigned int component_edge_count = evidence_status.size();
    // filter evidence such that we only count UNIQUE read support
    while (*std::max_element(evidence_status.begin(), evidence_status.end()) == 1) {
        // get maximum value over all evidence iterators
        std::vector< read_id_t > current_evidence;
        for (unsigned int idx = 0; idx < component_edge_count; idx++) {
            if (evidence_status.at(idx) == 1) {
                safe_edge_count_t mapID = index_to_mapID.at(idx);
                current_evidence.push_back(evidence_per_edge.at(mapID).front());
            }
        }
        std::sort(current_evidence.begin(), current_evidence.end());
        assert (!current_evidence.empty());
        read_id_t current_min = current_evidence.front();
        read_id_t current_max = current_evidence.back();
//        std::cout << "current_max " << current_max << " current_min " << current_min << std::endl;
        assert (current_max < 2*program_settings.original_readcount);
        assert (current_min < 2*program_settings.original_readcount);
        bool unique_min;
        if (current_evidence.size() == 1) {
            unique_min = true; // only evidence remaining so definitely unique
        }
        else if (current_min < current_evidence.at(1)) {
            unique_min = true;
        }
        else {
            unique_min = false;
        }
        for (unsigned int idx = 0; idx < component_edge_count; idx++) {
            if (evidence_status.at(idx) == 1) {
                safe_edge_count_t mapID = index_to_mapID.at(idx);
                std::list< read_id_t > evidence = evidence_per_edge.at(mapID);
                if (evidence.front() == current_min) {
                    if (unique_min) { // keep evidence
                        unique_evidence_per_edge.at(mapID).push_back(current_min);
                    }
                    evidence_per_edge.at(mapID).pop_front();
                    if (evidence.size() == 1) {
                        // reached end of evidence list, mark node as finished
                        evidence_status.at(idx) = 0;
                    }
                }
            }
        }
    }
    // finally analyze evidence and reduce branches by removing edges with
    // insufficient evidence
    std::cout << "effective_evidence_counts" << std::endl;
    for (auto ev : unique_evidence_per_edge) {
        int count = ev.second.size();
        std::cout << count << " ";
        if (count < min_evidence) {
            node_pair_t node_pair = evidenceIndexToEdge(ev.first);
            if (overlap_graph->checkEdge(node_pair.first, node_pair.second, false) < 0) {
                std::cout << "edge not found, reverse: " << overlap_graph->checkEdge(node_pair.second, node_pair.first, false) << std::endl;
                exit(1);
            }
            edges_to_remove.push_back(node_pair);
        }
        node_pair_t node_pair = evidenceIndexToEdge(ev.first);
        std::cout << "evidence load for edge " << node_pair.first << "," << node_pair.second << " : ";
        for (auto id : ev.second) {
            std::cout << id << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

safe_edge_count_t BranchReduction::edgeToEvidenceIndex(node_id_t node, node_id_t neighbor) {
    // map edge (node pair) to a unique index for storage of evidence in evidence_per_edge
    assert (node < overlap_graph->getVertexCount());
    assert (neighbor < overlap_graph->getVertexCount());
    safe_edge_count_t index = node * overlap_graph->getVertexCount() + neighbor;
    assert (index < (overlap_graph->getVertexCount())*(overlap_graph->getVertexCount()));
    return index;
}

node_pair_t BranchReduction::evidenceIndexToEdge(safe_edge_count_t index) {
    // map storage index back to original node pair
    assert (index < (overlap_graph->getVertexCount())*(overlap_graph->getVertexCount()));
    node_id_t node = floor(index / overlap_graph->getVertexCount());
    node_id_t neighbor = index - node * overlap_graph->getVertexCount();
//    std::cout << "index, node, neighbor: " << index << " " << node << " " << neighbor << std::endl;
    assert (node < overlap_graph->getVertexCount());
    assert (neighbor < overlap_graph->getVertexCount());
    return std::make_pair(node, neighbor);
}

bool BranchReduction::compareNodepairs (node_pair_t pair1, node_pair_t pair2) {
    if (pair1.first != pair2.first) {
        return pair1.first < pair2.first;
    }
    else {
        return pair1.second < pair2.second;
    }
}
