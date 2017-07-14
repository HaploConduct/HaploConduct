//============================================================================
// Name        : FindNextOverlaps2.cpp
// Author      : Jasmijn Baaijens
// Version     : 0.4.0
// License     : GNU GPL v3.0
// Project     : ViralQuasispecies
// Description : An alternative algorithm for finding the new overlaps:
//               instead of reconsidering all existing edges, we check all pairs of superreads that have a subread in common.
//============================================================================

#include <list>
#include <boost/dynamic_bitset.hpp>
#include <utility>
#include <fstream>
#include <sstream>
#include <algorithm> // std::max, std::min

#include "SRBuilder.h"

void SRBuilder::findNextOverlaps2() {
    std::cout << "FindNextOverlaps2...\n";
    // build an adjacency list mapping vertices to superreads
    std::deque<Read>::iterator it;
    std::list<unsigned int>::const_iterator itv;
    for (it = single_SR_vec.begin(); it != single_SR_vec.end(); it++) {
        Read* read_ptr = & (*it);
        std::unordered_map< unsigned int, SubreadInfo > subreadMap = it->get_subreadMap();
        for (auto node_it : subreadMap) {
            unsigned int node = node_it.first;
            nodes_to_SR.at(node).push_back(read_ptr);
        }
//        std::list< unsigned int > clique = it->get_sorted_clique(0);
//        for (auto node_it : clique) {
//            nodes_to_SR.at(node_it).push_back(read_ptr);
//        }
    }
    for (it = paired_SR_vec.begin(); it != paired_SR_vec.end(); it++) {
        std::list< unsigned int > clique = it->get_sorted_clique(1);
        Read* read_ptr = & (*it);
        for (auto node_it : clique) {
            nodes_to_SR.at(node_it).push_back(read_ptr);
        }
    }
    unsigned int SR_count = single_SR_vec.size() + paired_SR_vec.size();
    unsigned int V = overlap_graph->getVertexCount();
    unsigned int c = 0;
    for (auto it : nodes_to_SR) {
        if (it.size() > c) {
            c = it.size();
        }
    }
    std::cout << "SR_count = " << SR_count << " , V = " << V <<  ", c = " << c << "\n";
//    if (V*c*c < SR_count*SR_count) {
    if (true) {
        nodeDictApproach();
    }
    else {
        bitvecApproach();
    }
}

void SRBuilder::bitvecApproach() {
    std::cout << "bitvecApproach...\n";
    // construct a bitvector for every superread, indicating its subreads by node ID
    std::vector< std::pair< Read*, boost::dynamic_bitset<> > > SR_to_bitvec;
    std::deque<Read>::iterator it;
    for (it = single_SR_vec.begin(); it != single_SR_vec.end(); it++) {
        Read* read_ptr = & (*it);
        SR_to_bitvec.push_back( std::make_pair( read_ptr, buildBitvec(read_ptr) ) );
    }
    for (it = paired_SR_vec.begin(); it != paired_SR_vec.end(); it++) {
        Read* read_ptr = & (*it);
        SR_to_bitvec.push_back( std::make_pair( read_ptr, buildBitvec(read_ptr) ) );
    }

    // build a list containing every pair of superreads having a subread in common
    std::list< OverlapCandidate > overlaps_list;
    std::pair< Read*, boost::dynamic_bitset<> > SR_bitvec1;
    std::pair< Read*, boost::dynamic_bitset<> > SR_bitvec2;
    for (unsigned int i = 0; i < SR_to_bitvec.size(); i++) {
        SR_bitvec1 = SR_to_bitvec.at(i);
        if (i % 10000 == 0 && i > 0) {
            std::cout << "i = " << i << "... ";
        }
        for (unsigned int j = i+1; j < SR_to_bitvec.size(); j++) {
            SR_bitvec2 = SR_to_bitvec.at(j);
            boost::dynamic_bitset<> product = SR_bitvec1.second & SR_bitvec2.second;
            if (!product.none()) {
                unsigned int node = product.find_first();
                OverlapCandidate SR_overlap;
                SR_overlap.SR1 = SR_bitvec1.first;
                SR_overlap.SR2 = SR_bitvec2.first;
                SR_overlap.common_node = node;
                if (program_settings.remove_trans == 1 && !SR_overlap.SR1->is_paired() && !SR_overlap.SR2->is_paired()) {
                    // make sure that the common node is first in one and second in the other superread
                    std::list< node_id_t > clique_SR1 = SR_overlap.SR1->get_sorted_clique(1);
                    std::list< node_id_t > clique_SR2 = SR_overlap.SR2->get_sorted_clique(1);
                    assert (clique_SR1.size() == 2 && clique_SR2.size() == 2);
                    if (clique_SR1.front() == clique_SR2.front() || clique_SR1.back() == clique_SR2.back()) {
                        continue;
                    }
                }
                overlaps_list.push_back(SR_overlap);
            }
        }
    }
    std::cout << "\n";
    // for every entry in this list, find the corresponding superread overlap
    std::string filename = PATH + "overlaps.txt";
    std::ofstream outfile(filename);
    int overlaps_count = 0;
    for (auto overlap_candidate : overlaps_list) {
        Overlap overlap = deduceOverlap(overlap_candidate);
        if (program_settings.no_inclusions && overlap.get_perc() == 100) {
            // ignore inclusion overlap
            continue;
        }
        if (overlap.get_len(1) > 0) {
            std::string line = overlap.get_overlap_line();
            outfile << line;
            overlaps_count++;
        }
    }
    std::cout << "Number of overlaps found: " << overlaps_count << "\n";
    outfile.close();
}

void SRBuilder::nodeDictApproach() {
    std::cout << "nodeDictApproach...\n";
    // keep track of superread edges already found
    overlaps_found = std::vector< std::set< read_id_t >> (new_read_count, std::set< read_id_t >());
    // build a list containing every pair of superreads having a subread in common
    std::list< OverlapCandidate > overlaps_list;
    unsigned int node_count = nodes_to_SR.size();
    int status = 0;
    std::cout << "Building list of overlapping superreads: \n";
    for (unsigned int node=0; node < node_count; node++) {
        std::vector< Read* > SR_list = nodes_to_SR.at(node);
        int perc = static_cast<int>((node*100)/node_count);
        if (perc % 10 == 0 && perc > status) {
            status = perc;
            std::cout << status << "%.. \n";
        }
        unsigned int count = SR_list.size();
        for (unsigned int i=0; i < count; i++) {
            Read* SR1 = SR_list.at(i);
            read_id_t id1 = SR1->get_read_id();
            for (unsigned int j=i+1; j < count; j++) {
                Read* SR2 = SR_list.at(j);
                read_id_t id2 = SR2->get_read_id();
                // check if overlap was already found
                unsigned int smallest = std::min(id1, id2);
                unsigned int largest = std::max(id1, id2);
                if (overlaps_found.at(smallest).count(largest) != 0) { // overlap was found before
                    continue;
                }
                overlaps_found.at(smallest).insert(largest);
                OverlapCandidate SR_overlap;
                SR_overlap.SR1 = SR1;
                SR_overlap.SR2 = SR2;
                SR_overlap.common_node = node;
                if (program_settings.remove_trans == 1 && !SR1->is_paired() && !SR2->is_paired()) {
                    // make sure that the common node is first in one and second in the other superread
                    std::list< node_id_t > clique_SR1 = SR_overlap.SR1->get_sorted_clique(0);
                    std::list< node_id_t > clique_SR2 = SR_overlap.SR2->get_sorted_clique(0);
                    assert (clique_SR1.size() == 2 && clique_SR2.size() == 2);
                    if (clique_SR1.front() == clique_SR2.front() || clique_SR1.back() == clique_SR2.back()) {
                        continue;
                    }
                }
                overlaps_list.push_back(SR_overlap);
            }
        }
    }
    std::cout << "100%\n";
    // for every entry in this list, find the corresponding superread overlap
    std::string filename = PATH + "overlaps.txt";
    std::ofstream outfile(filename);
    int overlaps_count = 0;
    std::cout << "Deducing and writing overlaps: \n";
    unsigned int total_count = overlaps_list.size();
    status = 0;
    for (auto overlap_candidate : overlaps_list) {
        int perc = static_cast<int>((overlaps_count*100)/total_count);
        if (perc % 10 == 0 && perc > status) {
            status = perc;
            std::cout << status << "%.. \n";
        }
        Overlap overlap = deduceOverlap(overlap_candidate);
        if (program_settings.no_inclusions && overlap.get_perc()) {
            // ignore inclusion overlap
            continue;
        }
        if (overlap.get_len(1) > 0) {
            std::string line = overlap.get_overlap_line();
            outfile << line;
            overlaps_count++;
        }
    }
    std::cout << "100% \n";
    std::cout << "Number of overlaps found: " << overlaps_count << "\n";
    next_overlaps_count += overlaps_count;
    outfile.close();
}

boost::dynamic_bitset<> SRBuilder::buildBitvec(Read* superread) {
    unsigned int V = overlap_graph->getVertexCount();
    boost::dynamic_bitset<> bitvec(V);
    std::list< unsigned int > clique;
    if (superread->is_paired()) {
        clique = superread->get_sorted_clique(1);
    }
    else {
        clique = superread->get_sorted_clique(0);
    }
    for (auto node_it : clique) {
        assert (node_it < V);
        bitvec[node_it] = 1;
    }
    return bitvec;
}

Overlap SRBuilder::deduceOverlap(OverlapCandidate candidate) {
    // overlap parameters to determine:
    read_id_t id1;
    read_id_t id2;
    unsigned int pos1;
    unsigned int pos2;
    std::string ord;
    std::string ori1;
    std::string ori2;
    unsigned int perc1;
    unsigned int perc2;
    unsigned int len1;
    unsigned int len2;
    std::string type1;
    std::string type2;

    Read* SR1 = candidate.SR1;
    Read* SR2 = candidate.SR2;
    unsigned int node = candidate.common_node;
    if (!(SR1->is_paired()) && !(SR2->is_paired())) { // S-S overlap
        int idx1 = findCliqueIndex(node, SR1, true, false);
        int idx2 = findCliqueIndex(node, SR2, true, false);
        unsigned int lenA = (SR1->get_seq(0)).length();
        unsigned int lenB = (SR2->get_seq(0)).length();
        if (idx1-idx2 >= 0) {
            id1 = SR1->get_read_id();
            id2 = SR2->get_read_id();
            pos1 = idx1 - idx2;
            if (pos1 > lenA) { // no overlap (due to opposite subread ends being removed)
                assert (program_settings.error_correction); // this has to happen due to error correction
                return Overlap(0, 0, 0, 0, "-", "-", "-", 0, 0, 0, 0, "s", "s"); // this overlap will be ignored
            }
            len1 = std::min(lenA - pos1, lenB);
        }
        else {
            id1 = SR2->get_read_id();
            id2 = SR1->get_read_id();
            pos1 = idx2 - idx1;
            if (pos1 > lenB) { // no overlap (due to opposite subread ends being removed)
                assert (program_settings.error_correction); // this has to happen due to error correction
                return Overlap(0, 0, 0, 0, "-", "-", "-", 0, 0, 0, 0, "s", "s"); // this overlap will be ignored
            }
            len1 = std::min(lenA, lenB - pos1);
        }
        perc1 = (int)floor(std::max(len1/float(lenA), len1/float(lenB))*100);
        if (perc1 > 100) {
            std::cout << "len1=" << len1 << ", lenA=" << lenA << ", lenB=" << lenB << ", idx1=" << idx1 << ", idx2=" << idx2 << ", pos=" << pos1 << "\n";
        }
        // the remaining parameters are redundant here
        pos2 = 0;
        ord = "-";
        ori1 = "+";
        ori2 = "+";
        perc2 = 0;
        len2 = 0;
        type1 = "s";
        type2 = "s";
    }
    else if ((SR1->is_paired()) && !(SR2->is_paired())) { // P-S overlap
        int idx1l = findCliqueIndex(node, SR1, true, false);
        int idx1r = findCliqueIndex(node, SR1, false, false);
        int idx2l = findCliqueIndex(node, SR2, true, false);
        int idx2r = findCliqueIndex(node, SR2, true, true);
        if (idx1l-idx2l >= 0) {
            id1 = SR1->get_read_id();
            id2 = SR2->get_read_id();
            pos1 = idx1l - idx2l;
            len1 = (SR1->get_seq(1)).length() - pos1;
            if (len1 <= 0) {
                return Overlap(0, 0, 0, 0, "-", "-", "-", 0, 0, 0, 0, "p", "s"); // this overlap will be ignored
            }
            type1 = "p";
            type2 = "s";
        }
        else {
            id1 = SR2->get_read_id();
            id2 = SR1->get_read_id();
            pos1 = idx2l - idx1l;
            len1 = (SR1->get_seq(1)).length();
            if (pos1 >= len1) {
                return Overlap(0, 0, 0, 0, "-", "-", "-", 0, 0, 0, 0, "p", "s"); // this overlap will be ignored
            }
            type1 = "s";
            type2 = "p";
        }
        perc1 = (int)floor(len1/float((SR1->get_seq(1)).length())*100);
        pos2 = idx2r - idx1r;
        len2 = std::min((SR1->get_seq(2)).length(), (SR2->get_seq(0)).length()-pos2);
        if (len2 <= 0 || pos2 < 0) {
            return Overlap(0, 0, 0, 0, "-", "-", "-", 0, 0, 0, 0, "p", "s"); // this overlap will be ignored
        }
        perc2 = (int)floor(len2/float((SR1->get_seq(2)).length())*100);
        ord = "-";
        ori1 = "+";
        ori2 = "+";
    }
    else if (!(SR1->is_paired()) && (SR2->is_paired())) { // S-P overlap
        int idx1l = findCliqueIndex(node, SR1, true, false);
        int idx1r = findCliqueIndex(node, SR1, true, true);
        int idx2l = findCliqueIndex(node, SR2, true, false);
        int idx2r = findCliqueIndex(node, SR2, false, false);
        if (idx1l-idx2l >= 0) {
            id1 = SR1->get_read_id();
            id2 = SR2->get_read_id();
            pos1 = idx1l - idx2l;
            len1 = (SR2->get_seq(1)).length();
            if (pos1 >= len1) {
                return Overlap(0, 0, 0, 0, "-", "-", "-", 0, 0, 0, 0, "s", "p"); // this overlap will be ignored
            }
            type1 = "s";
            type2 = "p";
        }
        else {
            id1 = SR2->get_read_id();
            id2 = SR1->get_read_id();
            pos1 = idx2l - idx1l;
            len1 = (SR2->get_seq(1)).length() - pos1;
            if (len1 <= 0) {
                return Overlap(0, 0, 0, 0, "-", "-", "-", 0, 0, 0, 0, "s", "p"); // this overlap will be ignored
            }
            type1 = "p";
            type2 = "s";
        }
        perc1 = (int)floor(len1/float((SR2->get_seq(1)).length())*100);
        pos2 = idx1r - idx2r;
        len2 = std::min((SR2->get_seq(2)).length(), (SR1->get_seq(0)).length()-pos2);
        if (len2 <= 0 || pos2 < 0) {
            return Overlap(0, 0, 0, 0, "-", "-", "-", 0, 0, 0, 0, "s", "p"); // this overlap will be ignored
        }
        perc2 = (int)floor(len2/float((SR2->get_seq(2)).length())*100);
        ord = "-";
        ori1 = "+";
        ori2 = "+";
    }
    else { // P-P overlap
        int idx1l = findCliqueIndex(node, SR1, true, false);
        int idx1r = findCliqueIndex(node, SR1, false, false);
        int idx2l = findCliqueIndex(node, SR2, true, false);
        int idx2r = findCliqueIndex(node, SR2, false, false);
        unsigned int lenA = (SR1->get_seq(1)).length();
        unsigned int lenB = (SR2->get_seq(1)).length();
        unsigned int lenC = (SR1->get_seq(2)).length();
        unsigned int lenD = (SR2->get_seq(2)).length();
        bool front_ord, back_ord;
        if (idx1l-idx2l >= 0) {
            id1 = SR1->get_read_id();
            id2 = SR2->get_read_id();
            pos1 = idx1l - idx2l;
            len1 = std::min(lenA - pos1, lenB);
            front_ord = 1;
        }
        else {
            id1 = SR2->get_read_id();
            id2 = SR1->get_read_id();
            pos1 = idx2l - idx1l;
            len1 = std::min(lenA, lenB - pos1);
            front_ord = 0;
        }
        if (idx1r-idx2r >= 0) {
            pos2 = idx1r - idx2r;
            len2 = std::min(lenC - pos2, lenD);
            back_ord = 1;
        }
        else {
            pos2 = idx2r - idx1r;
            len2 = std::min(lenC, lenD - pos2);
            back_ord = 1;
        }
        if (len1 <= 0 || len2 <= 0) {
            return Overlap(0, 0, 0, 0, "-", "-", "-", 0, 0, 0, 0, "p", "p"); // this overlap will be ignored
        }
        perc1 = (int)floor(std::max(len1/float(lenA), len1/float(lenB))*100);
        perc2 = (int)floor(std::max(len2/float(lenC), len2/float(lenD))*100);
        if ((front_ord && back_ord) || (!front_ord && !back_ord)) {
            ord = "1";
        }
        else {
            ord = "2";
        }
        ori1 = "+";
        ori2 = "+";
        type1 = "p";
        type2 = "p";
    }
    // write data to overlap
    Overlap overlap(id1, id2, pos1, pos2, ord, ori1, ori2, perc1, perc2, len1, len2, type1, type2);
    return overlap;
}
