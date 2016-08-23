//============================================================================
// Name        : ViralQuasispecies.cpp
// Author      : Jasmijn Baaijens
// Version     : 0.01 Beta
// License     : GNU GPL v3.0
// Project     : ViralQuasispecies
// Description : Additional graph algorithms to extend OverlapGraph.cpp
//============================================================================

#include <fstream>
#include <boost/dynamic_bitset.hpp>
#include <boost/timer.hpp>
#include <set>
#include <algorithm> // std::random_shuffle
#include <iterator> // std::next, std::back_inserter
#include <cstdlib> // std::rand, std::srand

#include "OverlapGraph.h"

std::vector< std::vector< node_id_t > > OverlapGraph::getEdgesForMerging() {
    boost::dynamic_bitset<> bitvec(vertex_count);
    std::vector< std::vector< node_id_t > > node_vec;
    node_id_t node = 0;
    for (auto adj_list : adj_out) {
        if (!bitvec[node] && adj_list.size() > 0) {
            std::vector< std::pair<node_id_t, int> > pairs;
            node_id_t outneighbor;
            for (auto edge_it : adj_list) {
                outneighbor = edge_it.get_vertex(2);
                assert (node == edge_it.get_vertex(1));
                int perc = edge_it.get_perc();
                pairs.push_back(std::make_pair(outneighbor, perc));
            }
//            // sort outgoing edges by decreasing overlap percentages
//            std::sort(pairs.begin(), pairs.end(), [=](const std::pair<node_id_t, int>& a, const std::pair<node_id_t, int>& b)
//            {
//                return a.second > b.second;
//            }
//            );
            for (auto pair_it : pairs) {
                outneighbor = pair_it.first;
                if (!bitvec[outneighbor]) {
                    std::vector< node_id_t > edge_nodes = {node, outneighbor};
                    node_vec.push_back(edge_nodes);
                    bitvec[node] = 1;
                    bitvec[outneighbor] = 1;
                    break;
                }
            }
        }
        node++;
    }
    return node_vec;
}

std::vector< node_id_t > OverlapGraph::sortVerticesByIndegree() {
    std::vector< std::pair<node_id_t, unsigned int> > pairs;
    node_id_t v = 0;
    for (auto it : adj_in) {
        pairs.push_back(std::make_pair(v, it.size()));
        v++;
    }
    std::sort(pairs.begin(), pairs.end(), [=](const std::pair<node_id_t, unsigned int>& a, const std::pair<node_id_t, unsigned int>& b)
    {
        return a.second < b.second;
    }
    );
    assert (pairs.size() == vertex_count);
    std::vector< node_id_t > sorted_vertices;
    for (auto pair_it : pairs) {
        node_id_t node = pair_it.first;
        assert (node >= 0 && node < vertex_count);
        sorted_vertices.push_back(node);
    }
    assert (sorted_vertices.size() == vertex_count);
    return sorted_vertices;
}

void OverlapGraph::vertexLabellingHeuristic(unsigned int & conflict_count) {
    std::cout << "Applying vertex labelling heuristic...\n";
    std::list< Edge > min_edges_to_be_moved;
    std::list< Edge > min_edges_to_be_deleted;
    boost::dynamic_bitset<> opt_orientations;
    opt_orientations.resize(vertex_count, true); // all-ones vector: initially all labels forward
    if (program_settings.add_duplicates) { // second half of vertices are duplicates, hence reversed
        assert (!program_settings.resolve_orientations);
        unsigned int readcount = fastq_storage->get_readcount();
        for (node_id_t i=readcount; i<vertex_count; i++) { // change labels for duplicates
            opt_orientations[i] = 0;
        }
    }
    else if (program_settings.resolve_orientations) {
        std::cout << "resolving vertex orientations by BFS\n";
        labelVertices(min_edges_to_be_moved, min_edges_to_be_deleted, opt_orientations);
        // try k possible labellings and choose the best one
        std::list< Edge > edges_to_be_moved;
        std::list< Edge > edges_to_be_deleted;
        unsigned int delete_count = min_edges_to_be_deleted.size();
        boost::dynamic_bitset<> orientations;
        int count = 1;
        while (count < 100 && delete_count > 0) { // k = 100 tries
            count++;
            labelVertices(edges_to_be_moved, edges_to_be_deleted, orientations);
            if (edges_to_be_deleted.size() < delete_count) {
                min_edges_to_be_deleted = edges_to_be_deleted;
                min_edges_to_be_moved = edges_to_be_moved;
                delete_count = min_edges_to_be_deleted.size();
                opt_orientations = orientations;
            }
            edges_to_be_moved.clear();
            edges_to_be_deleted.clear();
        }   
        // move edges
        std::cout << "moving edges where necessary..\n";  
        for (auto edge_it : min_edges_to_be_moved) {
            node_id_t u = edge_it.get_vertex(1);
            node_id_t v = edge_it.get_vertex(2);
            removeEdge(v,u);
            addEdge(edge_it);
        }
        if (delete_count > 0) {
            // remove edges that were conflicting
            std::cout << "deleting " << delete_count << " conflicting edges..\n";  
            conflict_count = delete_count;
            for (auto edge_it : min_edges_to_be_deleted) {
                node_id_t u = edge_it.get_vertex(1);
                node_id_t v = edge_it.get_vertex(2);
                removeEdge(u,v);
            }
        }
        else {
            std::cout << "A perfect labelling was found after " << count << " tries, no conflicting edges." << std::endl;
        }
    }
    vertex_orientations = opt_orientations;
}

void OverlapGraph::labelVertices(std::list< Edge > & edges_to_be_moved, std::list< Edge > & edges_to_be_deleted, boost::dynamic_bitset<> & orientations) {
    // label vertices by BFS -> O(V+E)
    orientations.resize(vertex_count, true); // all-ones vector: initially all labels forward
    boost::dynamic_bitset<> visited(vertex_count);
    std::list<node_id_t> bfs;
    // sort vertices by increasing indegree
    std::vector< node_id_t > sorted_vertices = sortVerticesByIndegree();
    // process vertices in this order
    for (node_id_t i = 0; i < vertex_count; i++) {
        node_id_t start_node = sorted_vertices.at(i);
        if (!visited[start_node]) { 
            bfs.push_back(start_node);
            visited[start_node] = 1; // orientation = 1 (i.e. forward) by default
        }
        while (!bfs.empty()) {
            node_id_t node = bfs.front();
            bfs.pop_front();
            // copy all neighbors (both in and out) to vector and order randomly
            std::vector< node_id_t > adj_vec(adj_in.at(node).begin(), adj_in.at(node).end());
            for (auto edge_it : adj_out.at(node)) {
                node_id_t neighbor = edge_it.get_vertex(2);
                adj_vec.push_back(neighbor);
            }
            std::srand( unsigned( std::time(0) ) ); // set the random seed
            std::random_shuffle( adj_vec.begin(), adj_vec.end() );
            // recursively check all neighbors
            for (auto neighbor_it : adj_vec) {
                if (!visited[neighbor_it]) {
                    bfs.push_back(neighbor_it);
                    visited[neighbor_it] = 1;
                    Edge* edge = getEdgeInfo(node, neighbor_it);
                    if (edge->get_ori(1) == edge->get_ori(2)) {
                        orientations[neighbor_it] = orientations[node];
                    }
                    else {
                        orientations[neighbor_it] = !orientations[node];
                    }                   
                }
            }
        }
    }
    for (node_id_t i=0; i < vertex_count; i++) {
        assert (visited[i]);
    }
    // check all edges -> O(E)
    std::vector< std::list<Edge> >::iterator it1;
    std::list<Edge>::iterator it2;
    int switch_count = 0;
    int overall_count = 0;
    int move_count = 0;
    node_id_t i = 0;
    assert (edges_to_be_deleted.empty());
    assert (edges_to_be_moved.empty());
    for (it1 = adj_out.begin(); it1 != adj_out.end(); it1++) {
        std::list< Edge >::iterator it2 = it1->begin();
        unsigned int size = it1->size();
        for (unsigned int j=0; j < size; j++) {
            overall_count++;
            node_id_t u = it2->get_vertex(1);
            assert (u == i);
            node_id_t v = it2->get_vertex(2);
            bool ori1 = it2->get_ori(1);
            bool ori2 = it2->get_ori(2);
            bool type_v1 = orientations[u];
            bool type_v2 = orientations[v];
            if (ori1 == type_v1 && ori2 == type_v2) { // edge ok
                it2++;
            }
            else if ((ori1 == ori2 && type_v1 != type_v2) || (ori1 != ori2 && type_v1 == type_v2)) { // contradiction
//                    std::cout << "unsolvable vertex orientation due to conflicting edges... removing.\n"; 
                edges_to_be_deleted.push_back(*it2);
                it2++;
            }
            else { // edge agrees with labelling but the read orientations need to be flipped
                Edge edge = *it2;
                bool move = edge.switch_edge_orientation();
                assert (it2->get_vertex(1) == u);
                switch_count++;
                if (move) { // edge orientation was changed, so move edge to other adjacency list 
                    move_count++;
                    edges_to_be_moved.push_back(edge);
//                        edges_to_be_inserted.push_back(*it2); // store edge that has to be moved
                    assert (edge.get_ori(1) == type_v2 && edge.get_ori(2) == type_v1);
                    assert (edge.get_vertex(1) == v);
                    it2++;
//                        it2 = adj_out.at(u).erase(it2); 
//                        // also update adj_in because this is important for cycle check
//                        adj_in.at(v).remove(u);
//                        adj_in.at(u).push_back(v);          
                }
                else {
                    it2->switch_edge_orientation();
                    assert (it2->get_ori(1) == type_v1 && it2->get_ori(2) == type_v2);
                    it2++;
                }
            }
        }
        i++;
    }
}



void OverlapGraph::dfs_helper(node_id_t parent, node_id_t node, boost::dynamic_bitset<> &marked, boost::dynamic_bitset<> &visited, std::vector<node_id_t>& path, bool remove, bool threecycles, std::vector< std::vector<node_id_t> > &twocycle_table, std::set< std::vector<node_id_t> > &threecycle_set) {
//    std::cout << "dfs_helper..\n";
    if (marked[node]) {
//        std::cout << "marked\n";
        reportCycle(parent, node, remove);
        backedge_count++;
//        std::cout << "backedge found: " << parent << " to " << node << "\n";
//        std::cout << "backedge_count: " << backedge_count << "\n";
        if (threecycles) {
            // construct 2 adjacency bitvectors:
            boost::dynamic_bitset<> adj1(vertex_count); // going into parent
            for (auto it : adj_in.at(parent)) {
                adj1[it] = 1;
            }
//            std::cout << "adj1 done \n";
            boost::dynamic_bitset<> adj2(vertex_count); // going out of node
            for (auto edge_it : adj_out.at(node)) {
                node_id_t v = edge_it.get_vertex(2);
                adj2[v] = 1;
            }
//            std::cout << "adj2 done \n";
            // check if this backedges gives rise to a 2-cycle
            if (adj2[parent] == 1) {
                twocycle_table.at(node).push_back(parent);
                twocycle_table.at(parent).push_back(node);
                twocycle_count++;
            }
            // find all 3-cycles induced by this backedge
            boost::dynamic_bitset<> threecycle_nodes(vertex_count);
//            std::cout << "computing bitwise AND \n";
            threecycle_nodes = adj1 & adj2;
//            std::cout << "searching for 3-cycles... ";
            for (node_id_t i = 0; i < vertex_count; i++) {
                if (threecycle_nodes[i]) {
                    std::vector<node_id_t> cycle;
                    if (std::min({parent, node, i}) == i) {
                        cycle.push_back(i);
                        cycle.push_back(parent);
                        cycle.push_back(node);
                    
                    }
                    else if (std::min({parent, node, i}) == node) {
                        cycle.push_back(node);
                        cycle.push_back(i);
                        cycle.push_back(parent);
                    }
                    else {
                        cycle.push_back(parent);
                        cycle.push_back(node);
                        cycle.push_back(i);
                    }
                    threecycle_set.insert(cycle);
//                    threecycle_count++;
                }
            }
//            std::cout << "done \n";
        }
        int len_path = 0;
//        std::cout << "cycle path: " << node << " ";
        for (auto it = path.rbegin(); it != path.rend(); it++) {
            len_path++;
//            std::cout << *it << " ";
            if (*it == node) {                
//                std::cout << "cycle of length " << len_path << "\n";
                break;
            }
        }
//        std::cout << "\n"; 
    }
    else if (!visited[node]) { 
//        std::cout << "not visited\n";
        marked[node] = 1;
        path.push_back(node);
        if (path.size() > graph_depth) {
            graph_depth = path.size(); // NOTE: this only the length of the longest DFS path that was traversed, the graph depth could be larger
        }
        // sort neighbors by increasing value of pos1
//        std::cout << "sort neighbors\n";
        std::vector< std::pair<node_id_t, int> > pairs;
        for (auto it : adj_out.at(node)) {
            pairs.push_back(std::make_pair(it.get_vertex(2), it.get_pos(1)));
        }
        std::sort(pairs.begin(), pairs.end(), [=](const std::pair<node_id_t, int>& a, const std::pair<node_id_t, int>& b)
        {
            return a.second < b.second;
        }
        );
        std::vector< node_id_t > sorted_neighbors;
        for (auto pair_it : pairs) {
            node_id_t v = pair_it.first;
            assert (v >= 0 && v < vertex_count);
            sorted_neighbors.push_back(v);
        }
//        std::cout << "sorting done, now processing " << sorted_neighbors.size() << " neighbors\n";
        // process neighbors according to sorted list
        for (auto it_v : sorted_neighbors) {
            dfs_helper(node, it_v, marked, visited, path, remove, threecycles, twocycle_table, threecycle_set);
        }
//        std::cout << "processing done\n";
        marked[node] = 0;
        visited[node] = 1;
        assert (path.back() == node);
        path.pop_back();
    }
}

void OverlapGraph::removeCycles() {
    std::string filename = PATH + "cycles.txt";
    remove(filename.c_str());
    std::cout << "removeCycles.. current edge count " << edge_count << std::endl;
    // sort vertices by increasing order of indegree and process nodes in this order
    std::vector< node_id_t > sorted_vertices = sortVerticesByIndegree();
    std::cout << "vertices sorted by indegree..\n";
    // DFS to find cycles -> O(V+E)
    boost::dynamic_bitset<> visited(vertex_count);
    boost::dynamic_bitset<> marked(vertex_count);
    std::vector< std::vector<node_id_t> > empty_table;
    std::set< std::vector<node_id_t> > empty_set;
    std::vector<node_id_t> path;
    for (auto i : sorted_vertices) {
        if (!visited[i]) { 
            dfs_helper(vertex_count, i, marked, visited, path, true, false, empty_table, empty_set);
        }
    }
    if (backedge_count == 0) {
        std::cout << "Overlap graph is cycle-free :)\n";
    }
    else {
        std::cout << "\nTHE GRAPH CONTAINS CYCLES!!! Number of back-edges: " << backedge_count << "\n";
        std::cout << "New edge count " << edge_count << std::endl;
    }
}


void OverlapGraph::findCycles() {
    std::string filename = PATH + "cycles.txt";
    remove(filename.c_str());
    std::cout << "findCycles..\n";
    // sort vertices by increasing order of indegree and process nodes in this order
    std::vector< node_id_t > sorted_vertices = sortVerticesByIndegree();
    std::cout << "vertices sorted by indegree..\n";
    // DFS to find cycles -> O(V+E)
    boost::dynamic_bitset<> visited(vertex_count);
    boost::dynamic_bitset<> marked(vertex_count);
    std::vector< std::vector<node_id_t> > empty_table;
    std::set< std::vector<node_id_t> > empty_set;
    std::vector<node_id_t> path;
    for (auto i : sorted_vertices) {
        if (!visited[i]) { 
            dfs_helper(vertex_count, i, marked, visited, path, false, false, empty_table, empty_set);
        }
    }
    if (backedge_count == 0) {
        std::cout << "Overlap graph is cycle-free :)\n";
    }
    else {
        std::cout << "\nTHE GRAPH CONTAINS CYCLES!!! Number of back-edges: " << backedge_count << "\n";
    }
}

void OverlapGraph::findInduced3Cycles() {
    std::string filename = PATH + "cycles.txt";
    remove(filename.c_str());
    std::cout << "findInduced3Cycles..\n";
    // sort vertices by increasing order of indegree and process nodes in this order
    std::vector< node_id_t > sorted_vertices = sortVerticesByIndegree();
    // DFS to find cycles -> O(V+E)
    boost::dynamic_bitset<> visited(vertex_count);
    boost::dynamic_bitset<> marked(vertex_count);
    std::vector< node_id_t > empty_vec = {};
    std::vector< std::vector<node_id_t> > twocycle_vec(vertex_count, empty_vec); // vector storing 2-cycles in adj.list format
    std::set< std::vector<node_id_t> > threecycle_set; // vector storing 3-cycles
    std::vector<node_id_t> path;
    for (auto i : sorted_vertices) {
//    for (node_id_t i = 0; i < vertex_count; i++) {
        if (!visited[i]) { 
//            std::cout << "visiting " << i << "\n";
            dfs_helper(vertex_count, i, marked, visited, path, false, true, twocycle_vec, threecycle_set);
        }
    }
    threecycle_count = threecycle_set.size();
    std::cout << "Longest DFS path: " << graph_depth << "\n";
    std::cout << "Number of back-edges: " << backedge_count << "\n";
    // process 3-cycles: check for induced 3-cycle subgraphs
    std::cout << "Checking subgraphs induced by 3-cycles...\n";
    std::cout << "Number of 3-cycles: " << threecycle_count << "\n";
    std::cout << "Number of 2-cycles: " << twocycle_count << "\n";
    std::vector< std::vector<node_id_t> > induced_3cycles;
    unsigned int count1 = 0;
    unsigned int count2 = 0;
    for (auto it : twocycle_vec) {
        count1 += it.size();
        count2 += (std::set<node_id_t>(it.begin(), it.end()).size());
    }
    assert (count1 == count2);
    assert (count1 == 2*twocycle_count);
    node_id_t node1;
    node_id_t node2;
    node_id_t node3;
    for (auto cycle_it : threecycle_set) {
        assert (cycle_it.size() == 3);
        node1 = cycle_it.at(0);
        node2 = cycle_it.at(1);
        node3 = cycle_it.at(2);
//        std::cout << node1 << " " << node2 << " " << node3 << "\n";
        // check for two-cycles
        bool twocycle_found = 0;
        for (auto node_it : twocycle_vec.at(node1)) { 
            if (node_it == node2) { // node1-node2
                twocycle_found = 1;
                break;
            }
            if (node_it == node3) { // node1-node3
                twocycle_found = 1;
                break;
            }
        }
        for (auto node_it : twocycle_vec.at(node2)) {
            if (node_it == node3) { // node2-node3
                twocycle_found = 1;
                break;
            }
        }
        if (!twocycle_found) {
            induced_3cycles.push_back(cycle_it);
        }
    }
    std::cout << "Number of induced 3-cycle subgraphs: " << induced_3cycles.size() << "\n";
}


void OverlapGraph::findBranches() {
    std::cout << "findBranches..." << std::endl;
    if (program_settings.min_overlap_perc > 0) {
        std::cout << "NOTE: min_overlap_perc > 0" << std::endl;
//        return;
    }
    unsigned int removal_count_in = 0;
    unsigned int removal_count_out = 0;
    unsigned int node_count_in = 0;
    unsigned int node_count_out = 0;
    // find all outgoing branches
    for (node_id_t i = 0; i < vertex_count; i++) {
        bool cont = true;
        std::list< Edge > adj_list = adj_out.at(i);
        for (auto edge1 = adj_list.begin(); edge1 != adj_list.end() && cont; ++edge1) {
            node_id_t v1 = edge1->get_vertex(2);
            for (auto edge2 = edge1; ++edge2 != adj_list.end() && cont; /**/) {
                node_id_t v2 = edge2->get_vertex(2);
                if (checkEdge(v1, v2) == -1) {
                    removal_count_out += adj_list.size();
                    node_count_out++;
                    cont = false;
                }
            }
        }
    }
    // find all incoming branches
    for (node_id_t i = 0; i < vertex_count; i++) {
        bool cont = true;
        std::list< node_id_t > adj_list = adj_in.at(i);
        for (auto v1 = adj_list.begin(); v1 != adj_list.end() && cont; ++v1) {
            for (auto v2 = v1; ++v2 != adj_list.end() && cont; /**/) {
                if (checkEdge(*v1, *v2) == -1) {
                    removal_count_in += adj_list.size();
                    node_count_in++;
                    cont = false;
                }
            }
        }
    }
    std::cout << "Number of edges to be removed in order to eliminate all out-branches: " << removal_count_out << std::endl;
    std::cout << "Number of edges to be removed in order to eliminate all in-branches: " << removal_count_in << std::endl;
    std::cout << "Number of nodes out-disconnected: " << node_count_out << std::endl;
    std::cout << "Number of nodes in-disconnected: " << node_count_in << std::endl;
//    unsigned int remaining_out = edge_count - removal_count;
//    std::cout << "Number of edges remaining: " << remaining << std::endl;
}

void OverlapGraph::findBranchfreeGraph(std::vector< std::list< node_id_t > > & cur_adj_in, std::vector< std::list< node_id_t > > & cur_adj_out, std::set< node_id_t > & remove_in, std::set< node_id_t > & remove_out) {
    std::cout << "findBranchfreeGraph..." << std::endl;
    // assumes that transitive edges have been removed first
    if (program_settings.min_overlap_perc > 0) {
        std::cout << "NOTE: min_overlap_perc > 0" << std::endl;
    }
    // find all outgoing branches
    for (node_id_t i = 0; i < vertex_count; i++) {
        std::list< node_id_t > adj_list = cur_adj_out.at(i);
        if (adj_list.size() > 1) {
            remove_out.insert(i);
            remove_in.insert(adj_list.begin(), adj_list.end());
        }
    }
    // find all incoming branches
    for (node_id_t i = 0; i < vertex_count; i++) {
        std::list< node_id_t > adj_list = cur_adj_in.at(i);
        if (adj_list.size() > 1) {
            remove_in.insert(i);
            remove_out.insert(adj_list.begin(), adj_list.end());
        }
    }
    std::cout << "Number of nodes out-disconnected: " << remove_out.size() << std::endl;
    std::cout << "Number of nodes in-disconnected: " << remove_in.size() << std::endl;
}


void OverlapGraph::findTransEdges(std::vector< std::list< node_id_t > > & cur_adj_in, std::vector< std::list< node_id_t > > & cur_adj_out, std::vector< std::list< node_id_t > > & new_adj_in, std::vector< std::list< node_id_t > > & new_adj_out, bool removeTrans) {
    std::cout << "Find transitive edges..." << std::endl;
    // assumes adjacency lists are sorted
    unsigned int total_edges_kept = 0;
    unsigned int total_edges_removed = 0;
    node_id_t node1 = 0;
    node_id_t node2;
    std::list< node_id_t > list1;
    std::list< node_id_t > list2;
    for (auto adj_per_node : cur_adj_out) {
        list1 = adj_per_node;
        for (auto node_it : adj_per_node) {
            node2 = node_it;
            list2 = cur_adj_in.at(node2);
            bool transitive = nonemptyIntersect(list1, list2);
            if (transitive != removeTrans) {
                new_adj_out.at(node1).push_back(node2);
                new_adj_in.at(node2).push_back(node1);
                total_edges_kept++;
            }
            else {
                total_edges_removed++;
            }
        }
        node1++;
    }
    std::cout << "Total edges kept: " << total_edges_kept << std::endl;
    std::cout << "Total edges removed: " << total_edges_removed << std::endl;
}


bool OverlapGraph::nonemptyIntersect(std::list< node_id_t > & list1, std::list< node_id_t > & list2) {
    // assumes list1 and list2 are sorted
    auto it1 = list1.begin();
    auto it2 = list2.begin();
    while (it1 != list1.end() && it2 != list2.end()) {
        if (*it1 == *it2) { // nonempty intersection
            return true;
        }
        else if (*it1 < *it2) {
            it1++;
        }
        else {
            it2++;
        }
    }
    return false;
}

std::vector< std::list< node_id_t > > OverlapGraph::sortAdjLists(std::vector< std::list< node_id_t > > input_lists) {
    std::vector< std::list< node_id_t > > output_lists;
    for (auto list_it : input_lists) {
        list_it.sort();
        output_lists.push_back(list_it);
    }
    return output_lists;
}

std::vector< std::list< node_id_t > > OverlapGraph::sortAdjOut() {
    // sort adj_out by increasing outneighbor ID
    std::vector< std::list< node_id_t > > sorted_adj_out;
    for (auto adj_list : adj_out) {
        std::list< node_id_t > neighbors;
        node_id_t outneighbor;
        for (auto edge_it : adj_list) {
            outneighbor = edge_it.get_vertex(2);
            neighbors.push_back(outneighbor);
        }
        neighbors.sort();
        sorted_adj_out.push_back(neighbors);
    }
    return sorted_adj_out;
}

void OverlapGraph::removeBranches() {
    std::cout << "removeBranches..." << std::endl;
    // create sorted adjacency lists
    std::vector< std::list< node_id_t > > sorted_adj_in;
    std::vector< std::list< node_id_t > > sorted_adj_out;
    sorted_adj_in = sortAdjLists(adj_in);
    sorted_adj_out = sortAdjOut();
    // obtain graph without transitive edges
    std::vector< std::list< node_id_t > > new_adj_in;
    std::vector< std::list< node_id_t > > new_adj_out;
    new_adj_in = std::vector< std::list< node_id_t > > (vertex_count, std::list< node_id_t >{}); // vector of empty adjacency lists
    new_adj_out = std::vector< std::list< node_id_t > > (vertex_count, std::list< node_id_t >{}); // vector of empty adjacency lists
    findTransEdges(sorted_adj_in, sorted_adj_out, new_adj_in, new_adj_out, true);
    // find all branches
    std::set< node_id_t > remove_in;
    std::set< node_id_t > remove_out;
    findBranchfreeGraph(new_adj_in, new_adj_out, remove_in, remove_out);
    // remove branches from new graph
    for (auto node : remove_in) {
        new_adj_in.at(node).clear();
    }
    for (auto node: remove_out) {
        new_adj_out.at(node).clear();
    }
//    std::cout << "branches removed" << std::endl;
    // assign vertices to connected components in branch-free graph
    std::map< node_id_t, unsigned int > component_map;
    boost::dynamic_bitset<> visited(vertex_count);
    unsigned int current_component = 0;
    for (node_id_t i = 0; i < vertex_count; i++) {
        if (!visited[i]) {
            std::list< node_id_t > stack;
            stack.push_back(i);
            visited[i] = 1;
            while (!stack.empty()) {
                node_id_t node = stack.front();
                stack.pop_front();
                component_map.insert(std::make_pair(node, current_component));
                for (auto out_nb : new_adj_out.at(node)) {
                    if (!visited[out_nb]) {
                        stack.push_back(out_nb);
                        visited[out_nb] = 1;
                    }
                }
                for (auto in_nb : new_adj_in.at(node)) {
                    if (!visited[in_nb]) {
                        stack.push_back(in_nb);
                        visited[in_nb] = 1;
                    }
                }
            }
            current_component++;
        }
    }
    std::cout << "Total number of components " << current_component << std::endl;
    // 2. remove all edges of overlap graph between different components of branch-free graph
    std::list< std::pair< node_id_t, node_id_t > > edges_to_remove;
    for (node_id_t i = 0; i < vertex_count; i++) {
        for (auto edge : adj_out.at(i)) {
            assert (edge.get_vertex(1) == i);
            node_id_t j = edge.get_vertex(2);
            if (component_map.at(i) != component_map.at(j)) {
                edges_to_remove.push_back(std::make_pair(i, j));
            }
        }
    }
    for (auto node_pair : edges_to_remove) {
        removeEdge(node_pair.first, node_pair.second);
    }
    std::cout << "Number of edges remaining: " << edge_count << std::endl;
}

void OverlapGraph::removeTransitiveEdges() {
    if (program_settings.remove_trans == 0) {
        return;
    }
    std::cout << "removeTransitiveEdges..." << std::endl;
    // create sorted adjacency lists
    std::vector< std::list< node_id_t > > sorted_adj_in;
    std::vector< std::list< node_id_t > > sorted_adj_out;
    sorted_adj_in = sortAdjLists(adj_in);
    sorted_adj_out = sortAdjOut();
    // obtain graph of all transitive edges
    std::vector< std::list< node_id_t > > new_adj_in;
    std::vector< std::list< node_id_t > > new_adj_out;
    new_adj_in = std::vector< std::list< node_id_t > > (vertex_count, std::list< node_id_t >{}); // vector of empty adjacency lists
    new_adj_out = std::vector< std::list< node_id_t > > (vertex_count, std::list< node_id_t >{}); // vector of empty adjacency lists
    findTransEdges(sorted_adj_in, sorted_adj_out, new_adj_in, new_adj_out, false);
    if (program_settings.remove_trans == 1) { // remove double transitive edges
        node_id_t node1 = 0;
        for (auto neighbors : new_adj_out) {
            for (auto node2 : neighbors) {
                removeEdge(node1, node2);
            }
            node1++;
        }
    }
    else {
        // find all double transitive edges
        std::vector< std::list< node_id_t > > new_adj_in2;
        std::vector< std::list< node_id_t > > new_adj_out2;
        new_adj_in2 = std::vector< std::list< node_id_t > > (vertex_count, std::list< node_id_t >{}); // vector of empty adjacency lists
        new_adj_out2 = std::vector< std::list< node_id_t > > (vertex_count, std::list< node_id_t >{}); // vector of empty adjacency lists
        findTransEdges(new_adj_in, new_adj_out, new_adj_in2, new_adj_out2, false);
        if (program_settings.remove_trans == 2) { // remove double transitive edges
            node_id_t node1 = 0;
            for (auto neighbors : new_adj_out2) {
                for (auto node2 : neighbors) {
                    removeEdge(node1, node2);
                }
                node1++;
            }
        }
        else {
            // find all triple transitive edges
            std::vector< std::list< node_id_t > > new_adj_in3;
            std::vector< std::list< node_id_t > > new_adj_out3;
            new_adj_in3 = std::vector< std::list< node_id_t > > (vertex_count, std::list< node_id_t >{}); // vector of empty adjacency lists
            new_adj_out3 = std::vector< std::list< node_id_t > > (vertex_count, std::list< node_id_t >{}); // vector of empty adjacency lists
            findTransEdges(new_adj_in2, new_adj_out2, new_adj_in3, new_adj_out3, false);
            // remove triple transitive edges
            node_id_t node1 = 0;
            for (auto neighbors : new_adj_out3) {
                for (auto node2 : neighbors) {
                    removeEdge(node1, node2);
                }
                node1++;
            }
        }
    }
    std::cout << "Number of edges remaining: " << edge_count << std::endl;
}


///**************************************************
//                 Test functions
//**************************************************/                 

//void OverlapGraph::testCycles() {
//    std::vector< std::list< unsigned int> > test_adj_list;
//    for (unsigned int i = 0; i < 5; i++) {
//        std::list<unsigned int> list;
//        list.push_back((i + 1) % 5);
//        test_adj_list.push_back(list);
//    }
//    removeCycles(test_adj_list, 5);
//}

//void OverlapGraph::removeCycles(std::vector< std::list< node_id_t >> &tmp_adj_out, unsigned int V) {
//    std::cout << "removeCycles.. testing\n";
//    // DFS to find cycles -> O(V+E)
//    boost::dynamic_bitset<> visited(V);
//    boost::dynamic_bitset<> marked(V);
//    for (unsigned int i = 0; i < V; i++) {
//        if (!visited[i]) { 
//            dfs_helper(V, i, marked, visited, tmp_adj_out);
//        }
//    }
//    std::cout << "Number of back-edges: " << backedge_count << "\n";
//}


//void OverlapGraph::dfs_helper(node_id_t parent, node_id_t node, boost::dynamic_bitset<> &marked, boost::dynamic_bitset<> &visited, std::vector< std::list< node_id_t >> &tmp_adj_out) {
////    std::cout << "dfs_helper..\n";
//    if (marked[node]) {
//        reportCycle(parent, node, true);
//        backedge_count++;
//    }
//    else if (!visited[node]) {
//        marked[node] = 1;
//        std::list< node_id_t >::const_iterator it;
//        for (it = tmp_adj_out[node].begin(); it != tmp_adj_out[node].end(); it++) {
//            node_id_t neighbor = *it;
//            dfs_helper(node, neighbor, marked, visited, tmp_adj_out);
//        }
//        visited[node] = 1;
//        marked[node] = 0;
//    }
//}
