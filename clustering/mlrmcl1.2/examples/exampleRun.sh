../mlrmcl -c 500 synthetic.graph > synthetic.graph.c500.i2.0.b0.5.out
cut -f2 -d" " synthetic.groundtruth | sort -n | uniq -c > synthetic.groundtruth.sizes
../tools/fscores.sh synthetic.graph.c500.i2.0.b0.5 synthetic.groundtruth synthetic.groundtruth.sizes
python ../tools/normMI.py synthetic.groundtruth synthetic.graph.c500.i2.0.b0.5 > synthetic.graph.c500.i2.0.b0.5.nmi
