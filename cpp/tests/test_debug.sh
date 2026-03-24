cd output/
rm -rf *
Rscript ../../../scripts/generate_debug_test.R 



# debug compress tabular
../../bin/ldzip compress plinkTabular -l test.vcor -s test.vars -o test1 -m 0.0 -b 8 --min_col UNPHASED_R
../../bin/ldzip decompress -i test1 -o back -t tabular

# debug concat


../../bin/ldzip compress plinkSquare -l test1.bin -s test1.bin.vars -o test1 -m 0.0 -b 99
../../bin/ldzip compress plinkSquare -l test2.bin -s test2.bin.vars -o test2 -m 0.0 -b 99
../../bin/ldzip compress plinkSquare -l test3.bin -s test3.bin.vars -o test3 -m 0.0 -b 99
../../bin/ldzip compress plinkSquare -l test4.bin -s test4.bin.vars -o test4 -m 0.0 -b 99
../../bin/ldzip concat -i test1 test2 test3 test4 -o concat 
../../bin/ldzip decompress -i concat -o concat_full

