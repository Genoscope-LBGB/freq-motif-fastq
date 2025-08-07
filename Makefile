all:

install:
	module load rust && cargo install --root . --path freq-motif-fastq
	$(MAKE) install-bin

install-bin:
	[ -d $(PREFIX)/bin ] || mkdir $(PREFIX)/bin
	for bin in $(wildcard ./scripts/*); do \
		bin_basename=$$(basename $$bin) ; \
		target_file=$(PREFIX)/bin/$${bin_basename} ; \
		if [ -f $$target_file ] ; then \
			echo $${bin} is not unique or already exists ; \
			exit 1 ; \
		fi ; \
		cp $$bin $$target_file ; \
		chmod +x $$target_file ; \
	done
	for bin in $(wildcard ./bin/* ); do \
		bin_basename=$$(basename $$bin) ; \
		target_file=$(PREFIX)/bin/$${bin_basename} ; \
		if [ -f $$target_file ] ; then \
			echo $${bin} is not unique or already exists ; \
			exit 1 ; \
		fi ; \
		cp $$bin $$target_file ; \
		chmod +x $$target_file ; \
	done