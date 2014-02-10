

default:
#	(cd Programs ; make ; cd .. )
	(cd Lib ; make ) 
	(cd Programs ; make ) 
#	(cd Test ; make ) 

clean:
	(cd Lib ; make clean )
	(cd Programs ; make clean )
#	(cd Test ; make clean )

realclean:
	(cd Lib ; make realclean )
	(cd Programs ; make realclean )
#	(cd Test ; make realclean )



soft_rmcl_spec:
	./srmcl -b 1.0 -w -o ./output/WIPHI_softMLR-MCL ./graph/WIPHI_graph
	cp ./output/* ../7.evaluation_PPI/clustering_result/



toy_example:
	./srmcl -o ./output/toy1soft ./graphs/toy_example1
	./srmcl -t 50 -o ./output/toy2soft ./graphs/toy_example2
	./srmcl -o ./output/toy3soft ./graphs/toy_example3

	





