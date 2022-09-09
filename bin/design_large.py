#!/usr/bin/env python3
"""Design probes for genome capture, with options and parameters that
optimize resource usage for large, highly diverse input.

The downside of these options is usually a small increase in the number of
designed probes. As with design.py, this program still supports full
customization of the argument values.

This wraps design.py and offers a way to run design.py without requiring
deep familiarity with CATCH's options. That is, it takes into account
recommendations that often work well in practice.
"""

import design

__author__ = 'Hayden Metsky <hayden@broadinstitute.org>'


if __name__ == "__main__":
    args = design.init_and_parse_args(args_type='large')
    design.main(args)
