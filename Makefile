##########################################################################################
# This pipeline will take a fasta file in ./fasta and align its first two sequences
# according to the toycoati model.
#
# Example: make aln/example-001.fasta
#   This aligns fasta/example-001.fasta and outputs the result to aln/

RSCRIPT=Rscript --vanilla

default: all

.PHONY: default all coati

all: fst/toycoati.fst

.DELETE_ON_ERROR:

##########################################################################################
# Construct fst/toycoati.fst

# Construct an FST for a MG94 codon model
fst/mutation.fst: scripts/mutation.R fst/nuc_syms.txt
	$(RSCRIPT) $< | fstcompile --isymbols=fst/nuc_syms.txt --osymbols=fst/nuc_syms.txt - \
		| fstrmepsilon | fstarcsort --sort_type=olabel > $@

# Construct an FST for a geometric indel model
fst/indel.fst: scripts/indel.R fst/nuc_syms.txt
	$(RSCRIPT) $< | fstcompile --arc_type=standard --isymbols=fst/nuc_syms.txt --osymbols=fst/nuc_syms.txt - \
		| fstrmepsilon | fstarcsort --sort_type=ilabel > $@

# Construct an FST for a geometric indel model
fst/indel_raw.fst: scripts/indel.R fst/nuc_syms.txt
	$(RSCRIPT) $< | fstcompile --arc_type=standard --isymbols=fst/nuc_syms.txt --osymbols=fst/nuc_syms.txt - \
		| fstarcsort --sort_type=ilabel > $@

# Compose mutations and indels together
fst/toycoati_temp.fst: fst/mutation_coati.fst fst/indel_coati.fst
	fstcompose $^ > $@

# Optimize the toy-COATi
fst/toycoati.fst: fst/toycoati_temp.opt.fst
	fstrmepsilon $< $@

# Optimize mutation
fst/mutation_coati.fst: fst/mutation.opt.fst
	fstarcsort --sort_type=olabel $< $@

# optimize indels
fst/indel_coati.fst: fst/indel.opt.fst
	fstarcsort --sort_type=ilabel $< $@

# encode an FST
%.enc.fst %.codex: %.fst
	fstencode --encode_labels $< $*.codex $*.enc.fst

# reduce an FST to a more efficient form
%.min.fst: %.enc.fst
	fstrmepsilon $< | fstdeterminize | fstminimize > $@

# decode an FST
%.opt.fst: %.min.fst %.codex
	fstencode --decode $< $*.codex $@

##########################################################################################
# Construct pairwise alignment

# Extract input acceptor
work/in_tape/%.fst: fasta/%
	$(RSCRIPT) scripts/acceptor.R $< 1 \
		| fstcompile --isymbols=fst/nuc_syms.txt --osymbols=fst/nuc_syms.txt - \
		| fstarcsort --sort_type=olabel > $@

# Extract output acceptor
work/out_tape/%.fst: fasta/%
	$(RSCRIPT) scripts/acceptor.R $< 2 \
		| fstcompile --isymbols=fst/nuc_syms.txt --osymbols=fst/nuc_syms.txt - \
		| fstarcsort --sort_type=ilabel > $@

# Find alignment graph
work/graph/%.fst: work/in_tape/%.fst work/out_tape/%.fst fst/toycoati.fst
	fstcompose work/in_tape/$*.fst fst/toycoati.fst \
		| fstarcsort --sort_type=olabel \
		| fstcompose - work/out_tape/$*.fst \
		> $@

# Find shortest path through graph
work/path/%.fst: work/graph/%.fst
	fstshortestpath $< | fsttopsort > $@

# Covert shortest path into an alignment
aln/%: work/path/%.fst scripts/fasta.R
	fstprint --isymbols=fst/nuc_syms.txt --osymbols=fst/nuc_syms.txt $< \
		| $(RSCRIPT) scripts/fasta.R - > $@

##########################################################################################
# MISC

# create a graph of an FST
%.dot: %.fst
	fstdraw --isymbols=fst/nuc_syms.txt --osymbols=fst/nuc_syms.txt $< > $@

# print an FST
%.pdf: %.dot
	dot -Tpdf -o$@ $<

# remove intermediate files
clean:
	rm -f fst/toycoati* fst/mutation* fst/indel*

.PHONY: clean
