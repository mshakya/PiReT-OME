#!/usr/bin/env bash
set -e
rootdir=$( cd $(dirname $0) ; pwd -P )

if [ -z ${EDGE_HOME+x} ]; then
	EDGE_HOME="$rootdir/../../"
fi

test_result(){
	MainErrLog=$rootdir/TestOutput/error.log
	TestLog=$rootdir/TestOutput/log.txt
	Test=$rootdir/pipeline_test_both/differential_gene/prokaryote/test_prok/EdgeR/EdgeR.liver__spleen__sig.txt
	Expect=$rootdir/EdgeR.liver__spleen__sig.txt
	testName="PiReT RNA seq test";
	if cmp -s "$Test" "$Expect"
	then
		echo "$testName passed!"
		touch "$rootdir/TestOutput/test.success"
	else
		echo "$testName failed!"
		if [ -f "$TestLog" ]
		then
			cat $TestLog >> $MainErrLog
		fi
		touch "$rootdir/TestOutput/test.fail"
	fi
}

cd $rootdir
echo "Working Dir: $rootdir";
echo "EDGE HOME Dir: $EDGE_HOME";

if [ ! -f "$rootdir/TestOutput/test.success" ]
then
	rm -rf $rootdir/TestOutput
fi
#TODO: change runPipeline-piret to runPipeline after integrating runPipeline-piret RNAseq subroutine in runPipeline
perl $EDGE_HOME/runPipeline-piret -c $rootdir/config.txt -o $rootdir/TestOutput -cpu 4 -noColorLog || true

cp $rootdir/pipeline_test_both/differential_gene/prokaryote/test_prok/EdgeR/EdgeR.liver__spleen__sig.txt .

test_result;
