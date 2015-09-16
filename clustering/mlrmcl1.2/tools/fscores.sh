if [[ $# < 3 ]]
then
	echo "Usage: $0 <predictedClusterFile> <groundTruthClusterFile> <groundTruthClusterSizesFile>"
	echo "See the comments inside $0 for further details on the format of these files"
	exit
fi

#Argument 1 is predicted clustering file. can be
#either single column or double column.
#If 2 columns are present, it is assumed that first is the element
#id and the second is the cluster id; in the case of 1 column, the
#line no is assumed to be the element id.
predictFile=$1

#Argument 2 is file with ground truth clusters. It must be double
#column, with similar semantics as above. The file should be
#sorted by the first column. Each element can belong
#to multiple clusters.
gtFile=$2

#Argument 3 is groundTruthSizes. This 
#is double column, with first column being size
#of the cluster followed by the cluster id.
gtSizesFile=$3

scriptDir=`dirname $0`
evalScript=${scriptDir}/evaluateClusters.py

#total no. of vertices in the graph
nvtxs=`wc -l $predictFile | cut -f1 -d" "`

#predicted clusters with size less than min will not be
#considered for evaluation
min=2

#predicted clusters with size greater than max will not be
#considered for evaluation
max=$nvtxs

for f in $predictFile
do
	outFile=${f}.fscores
	python $evalScript $f $gtFile $gtSizesFile $min $nvtxs \
	$max 1 > $outFile
	echo "Output in ${outFile}"
done

