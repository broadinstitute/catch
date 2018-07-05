"""Functions to use in testing custom_cover_range_fn."""

def covers_abc(probe_seq, sequence, kmer_start, kmer_end,
                       full_probe_len, full_sequence_len):
    # probe_seq is an array; sequence is a str
    probe_seq = ''.join(probe_seq)

    # If 'ABC' is in the probe and sequence, this says that the
    # probe covers the 'ABC' (and that substring only)
    if 'ABC' in probe_seq and 'ABC' in sequence:
        i = sequence.index('ABC')
        return (i, i + len('ABC'))
    else:
        return None
