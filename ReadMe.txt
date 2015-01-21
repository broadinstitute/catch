Steps for making hybridization probe sequences (These are the Lassa steps but it would be the same for any):

All consensus sequence were concatenated into two single files (one for Nigeria, one for Sierra Leone).  The two regions were treated completely separately.

In some cases, duplicate sequences were found in the file.  This likely occurred because of problems building a consensus sequence for one of the segments (i.e. where the bad sequence was replaced by a general consensus sequence).  Duplicate sequences were removed.  Sequenced labeled as “Consensus” sequence were removed as these contained a lot of degenerate bases.

While we originally used the Picard Bait designer, we found that enough customization was needed to warrant making our own designer.
 
We designed the baits to tile across the Lassa genome creating a probe every 50 bases.  There would be two sets of adapters for each bait set (Nigeria, Sierra Leone, and the Spike Ins).  The adapter sets alternated for each probe:
 
Adapter A  ____________________|____________________
Adapter B            ____________________|______________
 
Where each line is a 100 bp sequence so that the probes with Adapter A do not overlap with other probes using Adapter A.
 
The adapters are as follows:
 
1) Sierra Leone
 
Adapter A:
Default universal
(for reference, ATCGCACCAGCGTGT and CACTGCGGCTCCTCA)
 
Adapter B:
---Adapter------SLLASVB------Adapter---
5’-TCGCGCCCATAACTCNNNNNNNNNNTGGTCGTAGAGCGCA-3’
 
2) Nigerian Lassa:
 
Adapter A:
---Adapter------NGLASVA------Adapter---
5’-TCAGGCTATGCGGCTNNNNNNNNNNGTCCTGGCGACGATG-3’
 
Adapter B:
---Adapter------NGLASVB------Adapter---
5’-AGCCGGTCAGTCGATAACACNNNNNNNNNNGTGCCTCGCTGGAAGTAGAC-3’
 
3) Spike In
 
Adapter A:
---Adapter------SPIKEINA------Adapter---
5’-GGTTCCGCGTCACAGNNNNNNNNNNCAACCGGACGGATCA-3’
 
Adapter B:
---Adapter------SPIKEINB------Adapter---
5’-ATCAGGCCCACAGATGGTAGNNNNNNNNNNGGCAATCCAGGGACAAAGTA-3’
%}
 
 
We found that the reverse complement of the probes produced by the BaitDesigner also reverse complemented the primers.  We checked this several times.
 
Duplicate probes were removed using Tim Fennell's RemoveSimilarSequences program.  Because we had a lot of room on our tiles and wanted to maximize the variety in our baits, the program removed all duplicates with a duplicate defined at 0 mismatches.   

We found that many of the probes had runs of of Ns in them.  This was caused by difficulties in the original consensus sequence.  While we can build probes with Ns in them Chris was worried that strings of N longer than 2 bases might cause problems making the probes and very long strings might pull random sequences.  

The Matlab script does not write probe sequences that have stretches of Ns in them.  At the end it goes back to the original file and adds back in a probe sequence that goes right up to the left and right of the stretch of Ns and adds a probe for these sequences.  In this case, both sequences denoted with X as well as the reverse complement would be added in:  

XXXXXXXXXXXXXXXXXNNNNNNNXXXXXXXXXXXXXXXX


This, of course, only occurs if there are at least 100 bases on either side, and those 100 bases do not also include stretches of Ns.

The Matlab script used is called makeViralProbesAlternatingAdapters.m
which is documented in the text itself

If you want to use the offset function in your duplicate removal you need to do some additional coding to make the probes and then add the primers later.  It considers the primers part of the probe so RemoveSimilarSequences doesn't work for that, so you have to modify the code or send it blank adapters and then put the adapters on later.  You can probably do this with a simple unix script.

It might be easier to rewrite the thing in whatever language you know best.
        
example of remove similar sequences:
        java -jar RemoveSimilarSequences.jar OFFSET_ALLOWED=20 MISMATCHES_ALLOWED=4 I=C:\Temp\ViralGenomeProbes_withOverlap.fasta O=Probes_unique4.fasta
   

Tarjei doesn't want the header sequences.  He just wants the probes so you can remove the > lines with a grep
        e.g. grep -v '>' Probes.fasta > justSeqs.txt

