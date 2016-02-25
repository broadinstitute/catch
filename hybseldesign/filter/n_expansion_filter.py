"""Returns probes where each 'N' base is replaced by real bases.

The 'N' base in a probe indicates an unknown -- i.e., the base
can either 'A', 'T', 'C', or 'G'.

This acts as a filter on the probes by returning, for each probe p:
 - if p does not contain an 'N' base, then p itself.
 - if p does contain one or more 'N' bases, then 4 or more probes
   in which the sequence is maintained except 'N' is expanded to
   be 'A', 'T', 'C', and 'G'.

For example, if a probe is 'ANA', then the following probes are
returned in place of 'ANA': 'AAA', 'ATA', 'ACA', 'AGA'. If a
probe is 'ANN', then the following probes are returned in place
of 'ANN': 'AAA', 'AAT', 'AAC', 'AAG', 'ATA', 'ATT', 'ATC', 'ATG',
'ACA', 'ACT', 'ACC', 'ACG', 'AGA', 'AGT', 'AGC', 'AGG'. If a
probe contains n 'N' bases, then 4^n probes are returned in place
of the probe.

The number of output probes is dependent on the number of 'N'
bases within each probe. The order of the input probes is
conserved -- i.e., the new, "expanded" probes are returned among
the other probes.
"""

from hybseldesign.filter.base_filter import BaseFilter
from hybseldesign import probe

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class NExpansionFilter(BaseFilter):
    """Filter that expands 'N' bases within probes.
    """

    def _filter(self, input):
        """Return input probes where 'N' bases are replaced with real bases.
        """
        real_bases = ['A', 'T', 'C', 'G']

        output = []
        for p in input:
            if 'N' not in p.seq_str:
                # p has no 'N' bases, so there is nothing to expand
                output += [p]
                continue

            expanded_probe_seqs = [p.seq_str]
            # Keep iterating (expanding) while there are still 'N'
            # bases left
            while [s for s in expanded_probe_seqs if 'N' in s]:
                expanded_probe_seqs_updated = []
                for s in expanded_probe_seqs:
                    N_pos = s.index('N')
                    if N_pos == -1:
                        # There is no need to expand s because there is no 'N'
                        expanded_probe_seqs_updated += [s]
                        continue
                    # Expand the first 'N' in s (at position N_pos)
                    s_list = list(s)
                    for b in real_bases:
                        s_list[N_pos] = b
                        expanded_probe_seqs_updated += [''.join(s_list)]
                expanded_probe_seqs = expanded_probe_seqs_updated

            for seq in expanded_probe_seqs:
                output += [probe.Probe.from_str(seq)]
        return output
