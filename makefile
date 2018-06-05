SRC_DIR    = src
BIN_DIR    = bin
QC_DIR     = quick-cliques

SYMLINK1   = savage
SYMLINK2   = polyte

CC = g++
CPPFLAGS = -Wall -fopenmp -std=c++11 -g -O2
LIBS = -l boost_timer -l boost_system -l boost_program_options

EXEC_NAMES = ViralQuasispecies

EXECS = $(BIN_DIR)/ViralQuasispecies

# SYM_LINKS = savage polyte

VPATH = src

.PHONY : all

all: $(EXECS) $(SYM_LINKS)
	@$(MAKE) -C $(QC_DIR)

.PHONY : clean

clean:
	rm -rf $(EXECS) $(BIN_DIR)
	# rm -f $(SYM_LINKS)
	@$(MAKE) -C $(QC_DIR) clean

$(BIN_DIR)/ViralQuasispecies: $(SRC_DIR)/*.cpp $(SRC_DIR)/*.h ${BIN_DIR}
	$(CC) $(CPPFLAGS) $(SRC_DIR)/*.cpp -o $@ ${LIBS}

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# $(SYM_LINKS):
# 	ln -sf savage.py savage
# 	ln -sf freq_est.py freq_est
