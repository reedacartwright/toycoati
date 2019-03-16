RSCRIPT=Rscript --vanilla

default: all

.PHONY: default all coati

all: fst/toycoati.fst

.DELETE_ON_ERROR:

##########################################################################################

# Construct an FST for a MG94 codon model
fst/mutation.fst: scripts/mutation.R fst/nuc_syms.txt
	$(RSCRIPT) $< | fstcompile --isymbols=fst/nuc_syms.txt --osymbols=fst/nuc_syms.txt - \
		| fstrmepsilon | fstarcsort --sort_type=olabel > $@

# Construct an FST for a geometric indel model
fst/indel.fst: scripts/indel.R fst/nuc_syms.txt
	$(RSCRIPT) $< | fstcompile --arc_type=standard --isymbols=fst/nuc_syms.txt --osymbols=fst/nuc_syms.txt - \
		| fstrmepsilon | fstarcsort --sort_type=ilabel > $@

# Compose mutations and indels together
fst/toycoati_.fst: fst/mutation_coati.fst fst/indel_coati.fst
	fstcompose $^ > $@

# Optimize the toy-COATi
fst/toycoati.fst: fst/toycoati_.opt_fst
	fstrmepsilon $< $@

# Optimize mutation
fst/mutation_coati.fst: fst/mutation.opt_fst
	fstarcsort --sort_type=olabel $< $@

# optimize indels
fst/indel_coati.fst: fst/indel.opt_fst
	fstarcsort --sort_type=ilabel $< $@

# encode an FST
fst/%.enc_fst fst/%.codex: fst/%.fst
	fstencode --encode_labels $< fst/$*.codex fst/$*.enc_fst

# reduce an FST to a more efficient form
fst/%.min_fst: fst/%.enc_fst
	fstrmepsilon $< | fstdeterminize | fstminimize > $@

# decode an FST
fst/%.opt_fst: fst/%.min_fst fst/%.codex
	fstencode --decode $< fst/$*.codex $@

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
