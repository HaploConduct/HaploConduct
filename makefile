SRC_DIR    = src
BIN_DIR    = bin
SFO_DIR    = sfo_2011_5
QC_DIR     = quick-cliques

CC = g++
CPPFLAGS = -Wall -fopenmp -std=c++11 -g -O0
LIBS = -l boost_timer -l boost_system -l boost_program_options

EXEC_NAMES = ViralQuasispecies

EXECS = $(BIN_DIR)/ViralQuasispecies

VPATH = src

.PHONY : all

all: $(EXECS)
	@$(MAKE) -C $(QC_DIR)
	@$(MAKE) -C $(SFO_DIR)

.PHONY : clean

clean: 
	rm -rf $(EXECS) $(BIN_DIR)
	@$(MAKE) -C $(QC_DIR) clean
	@$(MAKE) -C $(SFO_DIR) clean

$(BIN_DIR)/ViralQuasispecies: $(SRC_DIR)/*.cpp $(SRC_DIR)/*.h ${BIN_DIR}
	$(CC) $(CPPFLAGS) $(SRC_DIR)/*.cpp -o $@ ${LIBS}

$(BIN_DIR):
	mkdir -p $(BIN_DIR)
