{ head -n 1 test.tsv; tail -n +2 test.tsv | sort -t $'	' -k8,8nr; } | column -ts $'	' | less -S # descending sort pval; a large pval and a tiny effect size - subsample = population
