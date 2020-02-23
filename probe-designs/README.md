Probe sets
==========

This directory contains some of the probe sets designed with CATCH:
* `V-All.250k.201810.fasta.gz`: 249,999 probes that target genomes of all viruses known to infect humans, using sequence data through 2018-10 (588 species)
* `V-All.350k.201810.fasta.gz`: 349,999 probes that target genomes of all viruses known to infect humans, using sequence data through 2018-10 (588 species)
* `V-All.700k.201810.fasta.gz`: 699,999 probes that target genomes of all viruses known to infect humans, using sequence data through 2018-10 (588 species)
* `V-All.201606.fasta.gz`: 349,998 probes that target genomes of all viruses known to infect humans, using sequence data through 2016-06 (356 species); `V-All.201606.trimmed.fasta.gz` is this same probe set with bases trimmed on the 3' end to meet synthesis cycle limits
* `V-WAfr.201506.fasta.gz`: 44,995 probes that target genomes of viruses commonly circulating in West Africa, using sequence data through 2015-06 (23 species); this includes reverse complements as well, for 89,990 probes in total
* `V-Flu.45k.201810.fasta.gz`: 44,848 probes that target genomes of influenza A, B, and C viruses, using sequence data through 2018-10; this includes reverse complements as well, for 89,696 probes in total
* `V-Flu.6k.201810.fasta.gz`: 5,999 probes that target genomes of influenza A, B, and C viruses, using sequence data through 2018-10; this includes reverse complements as well, for 11,998 probes in total
* `V-Respiratory.100k.202001.fasta.gz`: 99,809 probes that target genomes of respiratory-related viruses, including SARS-CoV-2, using sequence data through 2020-02-18 (namely: enteroviruses, HPIV, HRSV, human-infecting coronaviruses, human mastadenoviruses, HMPV, Rhinovirus A/B/C, influenza A/B/C). Note that this targets all sequence data from the SARS-related coronavirus species, including bat and pangolin SARS-like CoV sequences.
* `V-MM.201603.fasta.gz`: 6,219 probes that target genomes of measles and mumps viruses, using sequence data through 2016-03; this includes reverse complements as well, for 12,438 probes in total
* `V-ZC.201602.fasta.gz`: 6,171 probes that target genomes of chikungunya and Zika viruses, using sequence data through 2016-02; this includes reverse complements as well, for 12,342 probes in total

Note that there are 20 nt PCR adapters on both ends of probes in `V-WAfr.201506.fasta.gz`, `V-Flu.45k.201810.fasta.gz`, `V-Flu.6k.201810.fasta.gz`, `V-MM.201603.fasta.gz`, and `V-ZC.201602.fasta.gz`.
