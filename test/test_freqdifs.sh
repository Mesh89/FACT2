k=$1
n=$2
while true; do
	./FACT -g $k $n 0 > trees.nex
	factres=`./FACT++ trees.nex 1`
	
	python gen_tnt_input.py $n $k trees.nex $factres > tnt_input.tre

	./tnt64/tnt proc tnt_input.tre

	res=`cat results.tre | grep '(' | tail -2`
	tntres=`echo $res | cut -d" " -f 1`
	factres=`echo $res | cut -d" " -f 2`
	
	if [ "$tntres" != "$factres" ]; then	
		echo "TNT: $tntres"
		echo "FACT: $factres"
		break
	fi
	
	rm results.tre
	rm trees.nex
done

