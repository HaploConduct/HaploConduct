SRC_DIR    = src
BIN_DIR    = bin
QC_DIR     = quick-cliques

CC = g++
CPPFLAGS = -Wall -fopenmp -std=c++11 -g -O2
LIBS = -l boost_timer -l boost_system -l boost_program_options

EXEC_NAMES = ViralQuasispecies

EXECS = $(BIN_DIR)/ViralQuasispecies

SYM_LINKS = haploconduct

VPATH = src

.PHONY : all

all: $(EXECS) $(SYM_LINKS)
	@$(MAKE) -C $(QC_DIR)

.PHONY : clean

clean:
	rm -rf $(EXECS) $(BIN_DIR)
	rm -f $(SYM_LINKS)
	@$(MAKE) -C $(QC_DIR) clean

$(BIN_DIR)/ViralQuasispecies: $(SRC_DIR)/*.cpp $(SRC_DIR)/*.h ${BIN_DIR}
	$(CC) $(CPPFLAGS) $(SRC_DIR)/*.cpp -o $@ ${LIBS}

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(SYM_LINKS):
	ln -sf haploconduct.py haploconduct
	chmod +x haploconduct
	chmod +x savage.py polyte.py polyte-split.py
