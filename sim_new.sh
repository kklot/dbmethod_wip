#!/bin/bash
#SBATCH --job-name=dbmt
#SBATCH --partition=fuchs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=3000
#SBATCH --mail-type=ALL
#SBATCH --time=32:00:00

NSV=(2 3 4 5 6 7)
SAMPLE=(1000) 
BIAS=(none logis)
TREND=(none increase decrease)

for sps in ${SAMPLE[@]} 
do
	for td in ${TREND[@]}
	do
		for nsv in ${NSV[@]} 
		do
			for ba in ${BIAS[@]} 
			do
				srun Rscript -e "rmarkdown::render(\"sim_new.Rmd\",\"html_document\", \"o${nsv}${sps}${ba}${td}\", \"html\", params=list(nsv=${nsv},sample_size=${sps},bias=\"${ba}\",trend=\"${td}\"))"
			done
		done
	done
done
exit 0
