#while true; do
	k=$1
	n=$2
	
	./FACT -g $k $n 0 > trees.nex
	factres=`./FACT++ trees.nex`
	
	python gen_tnt_input.py $n $k trees.nex $factres > tnt_input.tre

	./tnt64/tnt proc tnt_input.tre


#	rm results.tree
#	rm trees.nex
#done

