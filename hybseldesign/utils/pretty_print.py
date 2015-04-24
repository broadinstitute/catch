"""Functions for producing pretty strings to print.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'


def table(data, col_justify, header_underline=True):
    """Return a string with data formatted into a table.

    Args:
        data: 2D array in which data[i][j] is a string with the
            contents of the i'th row and j'th column
        col_justify: list, with the same number of columns as data,
            in which the j'th entry gives the justification of
            the j'th column of the table ("left", "right", or "center")
        header_underline: when True, includes a row of dashes immediately
            under the first row (the header)

    Returns:
        string in which data is formatted as a table
    """
    if len(data) == 0:
        return ''

    num_rows = len(data)
    num_cols = len(data[0])

    # Check that there is a consistent number of columns
    for row in data:
        if len(row) != num_cols:
            raise ValueError("data has inconsistent number of columns")
    if len(col_justify) != num_cols:
        raise ValueError("col_justify has incorrect number of entries")

    # Determine each column's width (maximum length of an entry)
    col_widths = [0 for j in xrange(num_cols)]
    for row in data:
        for j, col in enumerate(row):
            entry = str(col).rstrip()
            # The length of the entry is the maximum length over
            # all of its lines
            entry_len = max(len(line) for line in entry.split('\n'))
            col_widths[j] = max(col_widths[j], entry_len)

    # Determine each row's height (maximum lines of an entry)
    row_heights = [0 for i in xrange(num_rows)]
    for i, row in enumerate(data):
        for j, col in enumerate(row):
            entry = str(col).rstrip()
            num_lines = 1 + entry.count('\n')
            row_heights[i] = max(row_heights[i], num_lines)

    table = ''
    for i, row in enumerate(data):
        row_height = row_heights[i]
        for h in xrange(row_height):
            row_str = ''
            for j, col in enumerate(row):
                if j > 0:
                    # Add a space between columns
                    row_str += ' '
                entry = str(data[i][j]).rstrip()
                # Extract line h of entry
                if h > entry.count('\n'):
                    # Make this line empty
                    val = ""
                else:
                    val = entry.split('\n')[h]
                # Justify the value in the column
                if col_justify[j] == "left":
                    row_str += val.ljust(col_widths[j])
                elif col_justify[j] == "right":
                    row_str += val.rjust(col_widths[j])
                elif col_justify[j] == "center":
                    row_str += val.center(col_widths[j])
                else:
                    raise ValueError("Unknown column justification at %d" % j)
            table += row_str + '\n'
        if i == 0 and header_underline:
            # Add a row of dashes
            row_str = ''
            for j, width in enumerate(col_widths):
                if j > 0:
                    # Add a space between columns
                    row_str += ' '
                row_str += '-' * width
            table += row_str + '\n'
    return table
