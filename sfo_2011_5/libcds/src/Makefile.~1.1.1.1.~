CPP=g++

#CPPFLAGS=-g3 -Wall -O0
CPPFLAGS=-O9 -Wall -DNDEBUG 
#CPPFLAGS=-g

INCL=-I../includes/

CODERS_DIR=coders
CODERS_OBJECTS=$(CODERS_DIR)/huff.o $(CODERS_DIR)/huffman_codes.o

STATIC_PERMUTATION_DIR=static_permutation
STATIC_PERMUTATION_OBJECTS=$(STATIC_PERMUTATION_DIR)/perm.o $(STATIC_PERMUTATION_DIR)/static_permutation.o $(STATIC_PERMUTATION_DIR)/static_permutation_mrrr.o $(STATIC_PERMUTATION_DIR)/static_permutation_builder_mrrr.o

STATIC_BITSEQUENCE_DIR=static_bitsequence
STATIC_BITSEQUENCE_OBJECTS=$(STATIC_BITSEQUENCE_DIR)/static_bitsequence.o $(STATIC_BITSEQUENCE_DIR)/static_bitsequence_naive.o $(STATIC_BITSEQUENCE_DIR)/table_offset.o $(STATIC_BITSEQUENCE_DIR)/static_bitsequence_rrr02.o $(STATIC_BITSEQUENCE_DIR)/static_bitsequence_brw32.o $(STATIC_BITSEQUENCE_DIR)/static_bitsequence_builder_rrr02.o $(STATIC_BITSEQUENCE_DIR)/static_bitsequence_builder_brw32.o $(STATIC_BITSEQUENCE_DIR)/static_bitsequence_rrr02_light.o $(STATIC_BITSEQUENCE_DIR)/static_bitsequence_builder_rrr02_light.o $(STATIC_BITSEQUENCE_DIR)/static_bitsequence_sdarray.o $(STATIC_BITSEQUENCE_DIR)/sdarray.o $(STATIC_BITSEQUENCE_DIR)/static_bitsequence_builder_sdarray.o

STATIC_SEQUENCE_DIR=static_sequence
STATIC_SEQUENCE_OBJECTS=$(STATIC_SEQUENCE_DIR)/static_sequence.o $(STATIC_SEQUENCE_DIR)/static_sequence_wvtree.o $(STATIC_SEQUENCE_DIR)/wt_coder_binary.o $(STATIC_SEQUENCE_DIR)/wt_coder_huff.o $(STATIC_SEQUENCE_DIR)/wt_node_internal.o $(STATIC_SEQUENCE_DIR)/wt_node_leaf.o $(STATIC_SEQUENCE_DIR)/wt_coder.o $(STATIC_SEQUENCE_DIR)/wt_node.o $(STATIC_SEQUENCE_DIR)/static_sequence_gmr_chunk.o  $(STATIC_SEQUENCE_DIR)/static_sequence_builder_gmr_chunk.o  $(STATIC_SEQUENCE_DIR)/static_sequence_gmr.o $(STATIC_SEQUENCE_DIR)/static_sequence_builder_wvtree.o $(STATIC_SEQUENCE_DIR)/static_sequence_builder_gmr.o $(STATIC_SEQUENCE_DIR)/static_sequence_wvtree_noptrs.o $(STATIC_SEQUENCE_DIR)/static_sequence_builder_wvtree_noptrs.o $(STATIC_SEQUENCE_DIR)/static_sequence_bs.o

UTILS_DIR=utils
UTILS_OBJECTS=$(UTILS_DIR)/alphabet_mapper_none.o $(UTILS_DIR)/alphabet_mapper.o $(UTILS_DIR)/alphabet_mapper_cont.o

%.o: %.cpp
	@echo " [C++] Compiling $<"
	@$(CPP) $(CPPFLAGS) $(INCL) -c $< -o $@

all: lib 

clean:
	@echo " [CLN] Removing object files"
	@rm -f $(CODERS_OBJECTS) $(STATIC_BITSEQUENCE_OBJECTS) $(STATIC_SEQUENCE_OBJECTS) $(UTILS_OBJECTS) $(STATIC_PERMUTATION_OBJECTS)

lib: pre $(CODERS_OBJECTS) $(STATIC_BITSEQUENCE_OBJECTS) $(STATIC_SEQUENCE_OBJECTS) $(UTILS_OBJECTS) $(STATIC_PERMUTATION_OBJECTS)
	@echo " [LIB] Packing the object files"
	@ar rcs ../lib/libcds.a $(CODERS_OBJECTS) $(STATIC_BITSEQUENCE_OBJECTS) $(STATIC_SEQUENCE_OBJECTS) $(UTILS_OBJECTS) $(STATIC_PERMUTATION_OBJECTS)

pre:
	@echo " [HDR] Populating the includes folder"
	@cp basics.h ../includes/
	@cp */*.h ../includes/

