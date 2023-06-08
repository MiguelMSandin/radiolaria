#!/bin/bash

LIKELIHOODS="RAxML.likelihoods_100BS"

echo "Searching best tree"
for INFO in $(ls | grep "RAxML_info"); do
	like=$(grep "Final ML Optimization Likelihood: " $INFO | sed -e 's/.* //g')
	echo "$INFO: $like" >> $LIKELIHOODS
	echo "$INFO: $like"
done

MIN=$(cut -d ":" -f 2 $LIKELIHOODS | sort -n | tail -1)
BEST_TREE=$(grep "$MIN" $LIKELIHOODS | sed -e 's/: .*//g')
BEST_TREE=${BEST_TREE/.RAxML_info.txt/}

echo "Best tree found in: $BEST_TREE"

echo "" >> $LIKELIHOODS
echo "Best tree found in: $BEST_TREE" >> $LIKELIHOODS
