k=$1
n=$2
./FACT -g $k $n 0 > trees.nex

python gen_tnt_input.py $n $k trees.nex "(0 1)" > tnt_input.tre

