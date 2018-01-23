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

Since the number of output probes grows exponentially in the
number of 'N' bases within a probe, we can limit this with
the 'limit_n_expansion_randomly' parameter. When set to a nonnegative
integer, only limit_n_expansion_randomly 'N' bases are expanded;
these ones are chosen randomly. The other 'N' bases, if there
exist others, are randomly replaced with an unambiguous base.
When set to None, all 'N' bases are expanded.
"""

import random

from catch.filter.base_filter import BaseFilter
from catch import probe

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class NExpansionFilter(BaseFilter):
    """Filter that expands 'N' bases within probes.
    """

    def __init__(self, limit_n_expansion_randomly=3):
        """
        Args:
            limit_n_expansion_randomly: when set to a nonnegative integer,
                only this number of 'N' bases are expanded, and they
                are randomly chosen; the rest of are replaced with
                random unambiguous bases. When None, all 'N' bases
                are expanded
        """
        self.limit_n_expansion_randomly = limit_n_expansion_randomly

    def _filter(self, input):
        """Return input probes where 'N' bases are replaced with real bases.
        """
        real_bases = ['A', 'T', 'C', 'G']

        output = []
        for p in input:
            num_n = p.seq_str.count('N')
            if num_n == 0:
                # p has no 'N' bases, so there is nothing to expand
                output += [p]
                continue

            p_seq_init = p.seq_str
            if (self.limit_n_expansion_randomly is not None and
                    num_n > self.limit_n_expansion_randomly):
                # Randomly replace (num_n - self.limit_n_expansion_randomly)
                # 'N' bases with random unambiguous bases
                occurrences = [i for i, base in enumerate(p_seq_init)
                               if base == 'N']
                p_seq_init_list = list(p_seq_init)
                while len(occurrences) > self.limit_n_expansion_randomly:
                    occ_to_replace = random.choice(occurrences)
                    replacement = random.choice(real_bases)
                    p_seq_init_list[occ_to_replace] = replacement
                    occurrences.remove(occ_to_replace)
                p_seq_init = ''.join(p_seq_init_list)

            expanded_probe_seqs = [p_seq_init]
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
