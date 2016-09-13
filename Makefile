OBJ_EXT=src/ExtractRNA.o
OBJ_FUNC=src/assem.o src/gene_predict.o src/parse_seed.o
CC=g++
SYSLIB=-fopenmp
HASHLIB=-Wno-deprecated
BUILDFLG=-ffunction-sections -fdata-sections -fmodulo-sched -msse
EXE_TAX=bin/parallel-meta
EXE_RNA=bin/extract-rna
EXE_CLT=bin/class-tax
EXE_CLF=bin/class-func
EXE_TAL=bin/taxa-sel
EXE_FUL=bin/func-sel
EXE_CMP=bin/comp-sam
EXE_CMF=bin/comp-sam-func
EXE_FUNC=bin/parallel-meta-func
EXE_CMC=bin/comp-corr
EXE_PIP=bin/pipeline
EXE_SMG=bin/split-matrix
EXE_STG=bin/split-table
EXE_SSQ=bin/split-seq
EXE_FSQ=bin/format-seq
EXE_RCV=bin/rare-curv
EXE_CFN=bin/class-func-nsti
EXE_RAR=bin/rand-rare
EXE_UTX=bin/update-taxa

tax:$(OBJ_TAX) src/frame.cpp
	$(CC) -c -o $(OBJ_EXT) src/ExtractRNA.cpp $(HASHLIB)
	$(CC) -o $(EXE_TAX) src/frame.cpp $(OBJ_MAP) $(OBJ_EXT) $(SYSLIB) $(HASHLIB)
	$(CC) -o $(EXE_CLT) src/class_tax.cpp $(HASHLIB) $(BUILDFLG) $(SYSLIB)
	$(CC) -o $(EXE_CLF) src/class_func.cpp $(HASHLIB) $(BUILDFLG) $(SYSLIB)
	$(CC) -o $(EXE_RNA) src/ExtractRNA_plus.cpp $(OBJ_EXT) $(HASHLIB)	
	$(CC) -o $(EXE_TAL) src/taxa_sel.cpp $(HASHLIB) $(BUILDFLG) $(SYSLIB)
	$(CC) -o $(EXE_FUL) src/func_sel.cpp $(HASHLIB) $(BUILDFLG)
	$(CC) -o $(EXE_CMP) src/comp_sam.cpp $(HASHLIB) $(BUILDFLG) $(SYSLIB)
	$(CC) -o $(EXE_CMF) src/comp_sam_func.cpp $(HASHLIB) $(BUILDFLG) $(SYSLIB)
	$(CC) -o $(EXE_CMC) src/comp_corr.cpp $(HASHLIB) $(BUILDFLG) $(SYSLIB)
	$(CC) -o $(EXE_PIP) src/pipeline.cpp $(HASHLIB) $(BUILDFLG) $(SYSLIB)
	$(CC) -o $(EXE_SMG) src/split_matrix.cpp $(HASHLIB) $(BUILDFLG)
	$(CC) -o $(EXE_STG) src/split_table.cpp $(HASHLIB) $(BUILDFLG)
	$(CC) -o $(EXE_SSQ) src/split_seq.cpp $(HASHLIB) $(BUILDFLG)
	$(CC) -o $(EXE_FSQ) src/format_seq.cpp $(HASHLIB) $(BUILDFLG)
	$(CC) -o $(EXE_RCV) src/rare_curv.cpp $(HASHLIB) $(BUILDFLG) $(SYSLIB)
	$(CC) -o $(EXE_CFN) src/class_func_nsti.cpp $(HASHLIB) $(BUILDFLG) $(SYSLIB)
	$(CC) -o $(EXE_RAR) src/rand_rare.cpp $(HASHLIB) $(BUILDFLG)
	$(CC) -o $(EXE_UTX) src/update_taxa.cpp $(HASHLIB) $(BUILDFLG)

clean:
	rm -rf bin/* src/*.o
