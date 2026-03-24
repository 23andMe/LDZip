
install: deps build
	R CMD INSTALL LDZipMatrix_1.0.tar.gz

build: 
	Rscript -e "Rcpp::compileAttributes()"
	R CMD build ./

deps:
	Rscript -e "p<-c('Rcpp','RSQLite','DBI'); p<-setdiff(p, rownames(installed.packages())); if(length(p)) install.packages(p, repos='https://cloud.r-project.org')"

doc: 
	Rscript -e "roxygen2::roxygenise()"
	R CMD Rd2pdf . -o man/LDZipMatrix.pdf --force

test: test-data test-index-query test-rsid-query test-get-neighbor test-pairwise test-variant-parsing test-doc test-order test-fail test-na

test-index-query:
	Rscript -e "tinytest::run_test_file('inst/tinytest/test_index_query.R')"

test-doc:
	Rscript -e "tinytest::run_test_file('inst/tinytest/test_doc.R')"

test-order:
	Rscript -e "tinytest::run_test_file('inst/tinytest/test_order.R')"

test-fail:
	Rscript -e "tinytest::run_test_file('inst/tinytest/test_fail.R')"

test-rsid-query:
	Rscript -e "tinytest::run_test_file('inst/tinytest/test_rsid_query.R')"

test-get-neighbor:
	Rscript -e "tinytest::run_test_file('inst/tinytest/test_get_neighbor.R')"

test-pairwise:
	Rscript -e "tinytest::run_test_file('inst/tinytest/test_pairwise.R')"

test-variant-parsing:
	Rscript -e "tinytest::run_test_file('inst/tinytest/test_variant_parsing.R')"

test-na:
	Rscript -e "tinytest::run_test_file('inst/tinytest/test_na.R')"

BITS_LIST := 8 16 32 99
MIN_LIST  := 0 0.0001 0.001 0.01 0.1

test-data: clean-test-data test-data-matrix test-data-tabular test-data-na

clean-test-data:
	@rm -rf tests/data
	@mkdir -p tests/data

test-data-matrix:
	@Rscript ../scripts/generate_random_matrix.R 1000 tests/data/test_index
	@for bits in $(BITS_LIST); do \
	  for min in $(MIN_LIST); do \
	    out=tests/data/compressed_bits_index_$${bits}_min_$${min}; \
	    echo "Generating data with bits=$${bits}, min=$${min} -> $$out"; \
	    ../cpp/bin/ldzip compress plinkSquare \
	      --ld_file tests/data/test_index.bin \
	      --snp_file tests/data/test_index.bin.vars \
	      --output_prefix $$out \
	      --min $$min \
	      --bits $$bits > /dev/null; \
	  done; \
	done

test-data-tabular:
	@Rscript ../scripts/generate_random_tabular.R 1000 tests/data/test_tabular
	@for bits in $(BITS_LIST); do \
	  for min in $(MIN_LIST); do \
	    out=tests/data/compressed_bits_tabular_$${bits}_min_$${min}; \
	    echo "Generating data with bits=$${bits}, min=$${min} -> $$out"; \
	    ../cpp/bin/ldzip compress plinkTabular \
	      --ld_file tests/data/test_tabular.vcor \
	      --snp_file tests/data/test_tabular.vars \
	      --output_prefix $$out \
	      --min $$min \
	      --bits $$bits \
		  --min_col UNPHASED_R> /dev/null; \
	  done; \
	done
	
test-data-na:
	@Rscript ../scripts/generate_random_tabular.R 1000 tests/data/test_tabular_na D
	out=tests/data/compressed_na; \
	echo "Generating data with D with NA -> $$out"; \
	../cpp/bin/ldzip compress plinkTabular \
		--ld_file tests/data/test_tabular_na.vcor \
		--snp_file tests/data/test_tabular_na.vars \
		--output_prefix $$out \
		--min_col D > /dev/null


quick-test: clean-test-data test-data-tabular test-rsid-query
