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
    col_widths = [0 for i in xrange(num_rows)]
    for row in data:
        for j, col in enumerate(row):
            col_widths[j] = max(col_widths[j], len(col))

    table = ''
    for i, row in enumerate(data):
        row_str = ''
        for j, col in enumerate(row):
            if j > 0:
                # Add a space between columns
                row_str += ' '
            if col_justify[j] == "left":
                row_str += str(data[i][j]).ljust(col_widths[j])
            elif col_justify[j] == "right":
                row_str += str(data[i][j]).rjust(col_widths[j])
            elif col_justify[j] == "center":
                row_str += str(data[i][j]).center(col_widths[j])
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
