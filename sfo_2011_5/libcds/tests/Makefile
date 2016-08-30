CPP=g++

#CPPFLAGS=-g3 -Wall -I../includes/
CPPFLAGS=-O9 -Wall -DNDEBUG -I../includes/

OBJECTS=make_bitmap.o static_bitsequence_tester.o static_sequence_tester.o static_sequence_wvtree_test.o static_sequence_gmr_test.o static_sequence_gmr_chunk_test.o static_sequence_wvtree_noptrs_test.o static_bitsequence_test.o static_sequence_bs_test.o text_to_int.o
BIN=make_bitmap static_sequence_wvtree_test static_sequence_gmr_test static_sequence_gmr_chunk_test static_sequence_wvtree_noptrs_test static_bitsequence_test text_to_int static_sequence_bs_test

LIB=../lib/libcds.a

%.o: %.cpp
	@echo " [C++] Compiling $<"
	@$(CPP) $(CPPFLAGS) -c $< -o $@

all: $(OBJECTS) $(BIN)

static_bitsequence_test:
	@echo " [C++] Building static_bitsequence_test"
	@$(CPP) $(CPPFLAGS) -o static_bitsequence_test static_bitsequence_test.o static_bitsequence_tester.o $(LIB)

make_bitmap:
	@echo " [C++] Building make_bitmap"
	@$(CPP) $(CPPFLAGS) -o make_bitmap make_bitmap.o $(LIB)

text_to_int:
	@echo " [C++] Building text_to_int"
	@$(CPP) $(CPPFLAGS) -o text_to_int text_to_int.o $(LIB)

static_sequence_wvtree_test:
	@echo " [C++] Building static_sequence_wvtree_test"
	@$(CPP) $(CPPFLAGS) -o static_sequence_wvtree_test static_sequence_wvtree_test.o static_sequence_tester.o $(LIB) 

static_sequence_gmr_test:
	@echo " [C++] Building static_sequence_gmr_test"
	@$(CPP) $(CPPFLAGS) -o static_sequence_gmr_test static_sequence_gmr_test.o static_sequence_tester.o $(LIB) 

static_sequence_wvtree_noptrs_test:
	@echo " [C++] Building static_sequence_wvtree_noptrs_test"
	@$(CPP) $(CPPFLAGS) -o static_sequence_wvtree_noptrs_test static_sequence_wvtree_noptrs_test.o static_sequence_tester.o $(LIB) 

static_sequence_gmr_chunk_test:
	@echo " [C++] Building static_sequence_gmr_chunk_test"
	@$(CPP) $(CPPFLAGS) -o static_sequence_gmr_chunk_test static_sequence_gmr_chunk_test.o static_sequence_tester.o $(LIB) 

static_sequence_bs_test:
	@echo " [C++] Building static_sequence_bs_test"
	@$(CPP) $(CPPFLAGS) -o static_sequence_bs_test static_sequence_bs_test.o static_sequence_tester.o $(LIB) 

clean:
	@echo " [CLN] Cleaning object files"
	@rm -f $(OBJECTS) $(BIN)
