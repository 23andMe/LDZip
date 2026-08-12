
test:
	rm -rf .nf-test .nf-test.log
	CONSTRAIN_MEMORY=1 nf-test test tests/stub.nf.test --verbose
	CONSTRAIN_MEMORY=1 nf-test test tests/unit.nf.test --verbose
