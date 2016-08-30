
all: clean libcompact tests


doc:
	@echo " [DOC] Generating documentation"
	@doxygen

libcompact: 
	@echo " [MSG] Entering directory src"
	@make --no-print-directory -C src
	
tests: libcompact
	@echo " [MSG] Entering directory tests"
	@make --no-print-directory -C tests

clean:
	@echo " [MSG] Entering directory src"
	@make --no-print-directory -C src clean
	@echo " [MSG] Entering directory tests"
	@make --no-print-directory -C tests clean
	@echo " [CLN] Cleaning lib folder"
	@rm -f lib/*.a
	@touch lib/delete_me
	@echo " [CLN] Cleaning includes folder"
	@rm -f includes/*.h
	@touch includes/delete_me


