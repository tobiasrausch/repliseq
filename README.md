repliseq installation
---------------------

`git clone --recursive https://github.com/tobiasrausch/repliseq.git`

`cd repliseq`

`make all`


Running repliseq
----------------

The order of bam files on the command line must follow the cell-cycle.

`./src/repliseq -r <ref.fa> -o outprefix <g1b.bam> <s1.bam> <s2.bam> <s3.bam> <s4.bam> <g2.bam>`

Plotting tag densities for each cell-cycle fraction.

`Rscript R/reppattern.R -f outprefix.profile.tsv -r chr12:24000000-26000000`

Plotting the replication time.

`Rscript R/reptime.R -f outprefix.reptime.tsv -r chr12`
