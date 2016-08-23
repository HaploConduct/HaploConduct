//============================================================================
// Name        : Edge.h
// Author      : Jasmijn Baaijens
// Version     : 0.01 Beta
// License     : GNU GPL v3.0
// Project     : ViralQuasispecies
// Description : Edge class for overlap graph
//============================================================================

#ifndef EDGE_H_
#define EDGE_H_

#include <assert.h>
#include <algorithm> // std::swap

#include "Read.h"

class Edge
{
private:
    double score; // overlap score
    int pos1; // first overlap position
    int pos2; // second overlap position (for PAIRED overlap)
    int pos3; // first overlap position of reversed overlap w.r.t. read1
    int pos4; // second overlap position of reversed overlap w.r.t. read1
    bool ori1; // 1 if NORMAL, 0 if REVERSE (vertex 1)
    bool ori2; // similar for vertex 2
    Read* read1; // pointer to read corresponding to vertex1
    Read* read2; // pointer to read corresponding to vertex2
    char ord; // in case of PAIRED-PAIRED overlap denotes the read that comes first in overlap2
    unsigned int vertex1; // out-vertex
    unsigned int vertex2; // in-vertex
    int overlap_perc; // overlap percentage
    int overlap_len; // overlap length
    double mismatch_rate; // mismatch rate in overlap
    
public: 
    Edge(double s, int p1, int p2, bool orientation1, bool orientation2, std::string o, Read* r1, Read* r2) : 
        score(s), pos1(p1), pos2(p2), ori1(orientation1), ori2(orientation2), read1(r1), read2(r2) {
        ord = *(o.c_str());
        if (score < 0 && score != -1) {
            std::cout << score << std::endl;
        }
        assert (score == 0 || score == -1 || score > 0);
        overlap_len = -1;
        overlap_perc = -1;
        mismatch_rate = -1;
//        assert (orientation1 != 0 || orientation2 != 0);
//        vertex1 = read1->get_vertex_id(orientation1); // orientation determines whether we take normal or reverse vertex
//        vertex2 = read2->get_vertex_id(orientation2); // idem
    }
    
    unsigned int get_nonoverlap_len() {
        unsigned int len1 = read1->get_len();
        unsigned int len2 = read2->get_len();
        unsigned int nonoverlap_len = len1 + len2 - 2*overlap_len;
        return nonoverlap_len;
    }
    
    void set_mismatch(double mm_rate) {
        mismatch_rate = mm_rate;
    }
    
    double get_mismatch_rate() {
        assert ((mismatch_rate >= 0 && mismatch_rate <= 1) || mismatch_rate == -1);
        return mismatch_rate;
    }
    
    void swap_reads() {
        assert (pos1 == 0); // only swap reads when the order is undetermined, i.e. the overlap starts at position 0,
        assert (vertex1 > vertex2); // and such that the second vertex becomes largest
        std::swap(read1, read2);
        std::swap(vertex1, vertex2);
        std::swap(ori1, ori2);
        if (ord == '1') {
            ord = '2';
        }
        else if (ord == '2') {
            ord = '1';
        }
        pos3 = -pos3;
        pos4 = -pos4;
    }
    
    bool switch_edge_orientation() { // switch respectively ++, --, +- or -+ to --, ++, -+ or +-.
        bool ori_changed = false;
        std::swap(pos1, pos3);
        std::swap(pos2, pos4);
        ori1 = !ori1;
        ori2 = !ori2;
        assert (ord == '-' || ord == '1' || ord == '2');
        if (pos1 < 0) {
            std::swap(read1, read2);
            std::swap(vertex1, vertex2);
            std::swap(ori1, ori2);
            pos1 = -pos1;
            if (pos2 < 0) {
                ord = '1';
                pos2 = -pos2;
            }
            else if (ord != '-') {
                ord = '2';
            }
            ori_changed = true;
        }
        else {
            if (pos2 < 0) {
                pos2 = -pos2;
                ord = '2';
            }
            else if (ord != '-') {
                ord = '1';
            }
        }
        return ori_changed;
    }
    
    double get_score() const { return score; }
    
    void set_vertices(unsigned int v1, unsigned int v2) { 
        vertex1 = v1;
        vertex2 = v2;
    }
        
    unsigned int get_vertex(int i) const {
        assert (i == 1 || i == 2);
        if (i == 1) {
            return vertex1;
        }
        else {
            return vertex2;
        }
    }
        
    int get_pos(int i) const {
        assert (i == 1 || i == 2);
        if (i == 1) {
            return pos1;
        }
        else {
            return pos2;
        }
    }
        
    bool get_ori(int i) const {
        assert (i == 1 || i == 2);
        if (i == 1) {
            return ori1;
        }
        else {
            return ori2;
        }
    }
        
    char get_ord() const { return ord; }
        
    Read* get_read(int i) const {
        assert (i == 1 || i == 2);
        if (i == 1) {
            return read1;
        }
        else {
            return read2;
        }
    }
    
    void set_extra_pos(int p3, int p4=0) {
        pos3 = p3;
        pos4 = p4;
    }
    
    int get_extra_pos(int i) {
        assert (i == 1 || i == 2);
        if (i == 1) {
            return pos3;
        }
        else {
            return pos4;
        }
    }
    
    int get_perc() const {
        assert (overlap_perc >= 0);
        return overlap_perc;
    }
    
    void set_perc(int perc) {
        overlap_perc = perc;
    }
    
    int get_len() const {
        assert (overlap_len >= 0);
        return overlap_len;
    }
    
    void set_len(int len) {
        overlap_len = len;
    }
};

#endif /* EDGE_H_ */
