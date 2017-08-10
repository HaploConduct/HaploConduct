//============================================================================
// Name        : BranchReduction.h
// Author      : Jasmijn Baaijens
// Version     : 0.4.0
// License     : GNU GPL v3.0
// Project     : ViralQuasispecies
// Description : Use read evidence to reduce branches in the overlap graph
//============================================================================

#ifndef BRANCHREDUCTION_H_
#define BRANCHREDUCTION_H_

#include <list>
#include <stack>
#include <vector>
#include <set>
#include <iostream>

#include "Overlap.h"
#include "Types.h"
#include "Edge.h"
#include "FastqStorage.h"
#include "OverlapGraph.h"


class BranchReduction {
private:
    std::string PATH;
    std::shared_ptr<FastqStorage> fastq_storage;
    std::shared_ptr<FastqStorage> original_fastq;
    std::shared_ptr<OverlapGraph> overlap_graph;
    ProgramSettings program_settings;
    unsigned int SE_count;
    unsigned int PE_count;
    int min_evidence;
//    std::unordered_map< std::pair< node_id_t, node_id_t >, std::vector< read_id_t > > evidence_per_edge;

public:
    BranchReduction(std::shared_ptr<FastqStorage> fastq,
            std::shared_ptr<FastqStorage> originals,
            ProgramSettings ps,
            std::shared_ptr<OverlapGraph> graph) {
        PATH = ps.output_dir;
        fastq_storage = fastq;
        original_fastq = originals;
        overlap_graph = graph;
        program_settings = ps;
        SE_count = ps.branch_SE_c;
        PE_count = ps.branch_PE_c;
        min_evidence = ps.branch_min_ev;
    }
    ~BranchReduction() {}

    // BranchReduction.cpp
    void readBasedBranchReduction();
    std::vector< node_id_t > findBranchingEvidence(node_id_t node1, std::list< node_id_t > neighbors,
            std::list< std::pair< node_id_t, node_id_t > > & edges_to_remove,
            bool outbranch);
    std::list< int > buildDiffListOut(node_id_t node1,
            std::vector< node_id_t > neighbors,
            std::vector< std::string > & sequence_vec,
            std::vector< int > & startpos_vec,
            std::vector< std::pair< unsigned int, unsigned int > > & missing_edges);
    std::list< int > buildDiffListIn(node_id_t node1,
            std::vector< node_id_t > neighbors,
            std::vector< std::string > & sequence_vec,
            std::vector< int > & startpos_vec,
            std::vector< std::pair< unsigned int, unsigned int > > & missing_edges);
    std::vector< int > findDiffPos(std::string seq1, std::string seq2);
    bool checkReadEvidence(std::string contig, int startpos, std::string read, int index, std::list< int > diff_list);
};


#endif /* BRANCHREDUCTION_H_ */
