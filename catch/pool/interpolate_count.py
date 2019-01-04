"""Functions for interpolating probe count from parameter values.

We only computed probe counts at particular combinations of parameter
values (e.g., along a grid). When optimizing parameters, the search
will explore a parameter space that includes values for which the
probe count was not computed. This will interpolate probe counts
at those values.
"""

from collections import defaultdict
import logging
import math

import numpy as np
from scipy import interpolate

__author__ = 'Hayden Metsky <hayden@mit.edu>'

logger = logging.getLogger(__name__)


def _round_up(x, b):
    """Round float x up to the nearest multiple of int b
    """
    return int(math.ceil(float(x) / b)) * b

def _round_down(x, b):
    """Round float x down to the nearest multiple of int b
    """
    return int(math.floor(float(x) / b)) * b


def _make_interp_probe_count_for_dataset_standard_fn(probe_counts,
        cover_extension_scale=1.0/10):
    """Generate and return a function that interpolates probe count for a dataset.

    This operates only on the mismatches and cover_extension parameters.

    Args:
        probe_counts: dict giving number of probes for each dataset and
            choice of parameters
        cover_extension_scale: scale the cover_extension parameter by this
            amount relative to the mismatches parameter when calculating
            the area of a bounding box

    Returns:
        function whose input is a dataset, value for the mismatches parameter,
        and value for the cover extension parameter. The function linearly
        interpolates the number of probes required in that dataset for
        those parameter values, based on the values (which were explicitly
        calculated) in probe_counts.
    """
    memoized_bounding_boxes = {dataset: {} for dataset in probe_counts.keys()}
    def immediate_bounding_box(mismatches, cover_extension):
        """Return the box around the values rounded down/up to the nearest int
        to use when memoizing bounding boxes.
        """
        return (_round_down(mismatches, 1),
                _round_up(mismatches, 1),
                _round_down(cover_extension, 1),
                _round_up(cover_extension, 1))

    def find_bounding_box_around_point(dataset, mismatches, cover_extension):
        """Return a rectangular bounding box around given parameters.

        Based on probe_counts[dataset], choose parameter values that form the
        smallest box around (mismatches, cover_extension) so that we can
        use actual computed probe counts (based on the box's parameter values)
        to interpolate the probe count for (mismatches, cover_extension).
        """
        # 'mismatches' may be a float between the min/max mismatches for
        # dataset. 'cover_extension' may be a float between the min/max
        # cover_extension for dataset. Find the mismatches and cover_extension
        # parameters for dataset that form the smallest rectangle encompassing
        # 'mismatches' and 'cover_extension', so that we can use the points of
        # this rectangle to perform bilinear interpolation. We brute force over
        # all possible rectangles encompassing 'mismatches' and
        # 'cover_extension' (with some heuristics to not be completely
        # exhaustive), which takes O(n^3) time; there are faster methods, but
        # this should suffice here.
        # Consider mismatches on the x-axis (left to right) and cover_extension
        # on the y-axis (bottom to top)
        points = set(probe_counts[dataset].keys())
        points_topleft = set()
        points_topright = set()
        points_bottomleft = set()
        points_bottomright = set()
        for p in points:
            m, ce = p
            if m == mismatches:
                if ce == cover_extension:
                    points_topleft.add(p)
                    points_topright.add(p)
                    points_bottomleft.add(p)
                    points_bottomright.add(p)
                elif ce > cover_extension:
                    points_topleft.add(p)
                    points_topright.add(p)
                else:
                    points_bottomleft.add(p)
                    points_bottomright.add(p)
            elif m > mismatches:
                if ce == cover_extension:
                    points_topright.add(p)
                    points_bottomright.add(p)
                elif ce > cover_extension:
                    points_topright.add(p)
                else:
                    points_bottomright.add(p)
            else:
                if ce == cover_extension:
                    points_topleft.add(p)
                    points_bottomleft.add(p)
                elif ce > cover_extension:
                    points_topleft.add(p)
                else:
                    points_bottomleft.add(p)

        points_topright_by_y = defaultdict(set)
        for p in points_topright:
            m, ce = p
            points_topright_by_y[ce].add(p)
        points_bottomleft_by_x = defaultdict(set)
        for p in points_bottomleft:
            m, ce = p
            points_bottomleft_by_x[m].add(p)

        min_rectangle, min_area = None, float('inf')
        for p_topleft in points_topleft:
            p_topleft_m, p_topleft_ce = p_topleft
            # find all points in the top-right with the same cover_extension
            # (y-value) as p_topleft
            for p_topright in points_topright_by_y[p_topleft_ce]:
                p_topright_m, p_topright_ce = p_topright
                # find all points in the bottom-left with the same mismatches
                # (x-value) as p_topleft
                for p_bottomleft in points_bottomleft_by_x[p_topleft_m]:
                    p_bottomleft_m, p_bottomleft_ce = p_bottomleft
                    # to form a rectangle, we need the point (p_topright_m,
                    # p_bottomleft_ce); this point should be in the
                    # bottom-right, so check if it exists
                    p_bottomright = (p_topright_m, p_bottomleft_ce)
                    if p_bottomright in points_bottomright:
                        # we found a valid rectangle; now compute its 'area'
                        width = p_topright_m - p_topleft_m
                        height = ((p_topright_ce - p_bottomleft_ce) *
                            cover_extension_scale)
                        # add pseudocounts to width and height because if
                        # width or height is 0, we still want the other
                        # dimension to be accounted for
                        area = (width + 0.001) * (height + 0.001)
                        if area < min_area:
                            min_rectangle = (p_topleft, p_bottomright)
                            min_area = area
        return min_rectangle

    def interp_probe_count_for_dataset(dataset, param_vals):
        """
        Using the given probe counts at particular parameter values, interpolate
        the number of probes for 'dataset' and given mismatches (param_vals[0])
        and cover_extension (param_vals[1]), where each of these may be floats
        """
        mismatches, cover_extension = param_vals

        immediate_bb = immediate_bounding_box(mismatches, cover_extension)
        if immediate_bb in memoized_bounding_boxes[dataset]:
            # The bounding box for (mismatches, cover_extension) has been
            # memoized
            min_rectangle = memoized_bounding_boxes[dataset][immediate_bb]
        else:
            # Compute and memoize the bounding box for (mismatches,
            # cover_extension)
            min_rectangle = find_bounding_box_around_point(dataset,
                                                           mismatches,
                                                           cover_extension)
            if min_rectangle is None:
                raise Exception(("Unable to find rectangular bounding box around "
                                 "(mismatches, cover_extension)=(%f, %f) for "
                                 "dataset %s") % (mismatches, cover_extension,
                                 dataset))
            memoized_bounding_boxes[dataset][immediate_bb] = min_rectangle

        rect_topleft, rect_bottomright = min_rectangle
        mismatches_floor, cover_extension_ceil = rect_topleft
        mismatches_ceil, cover_extension_floor = rect_bottomright

        # At cover_extension_floor and at cover_extension_ceil, interpolate the
        # number of probes at 'mismatches' mismatches (interpolate linearly)
        for ce in [cover_extension_floor, cover_extension_ceil]:
            count_left = probe_counts[dataset][(mismatches_floor, ce)]
            count_right = probe_counts[dataset][(mismatches_ceil, ce)]
            mismatches_diff = mismatches_ceil - mismatches_floor
            if mismatches_diff == 0:
                # count_left should equal count_right
                assert count_left == count_right
                count = count_left
            elif count_left <= count_right:
                count_diff = count_right - count_left
                f = float(mismatches - mismatches_floor) / mismatches_diff
                count = f * count_diff + count_left
            else:
                count_diff = count_left - count_right
                f = float(mismatches - mismatches_floor) / mismatches_diff
                count = count_left - f * count_diff
            if ce == cover_extension_floor:
                count_floor = count
            if ce == cover_extension_ceil:
                count_ceil = count

        # Interpolate between cover_extension_floor and cover_extension_ceil
        # using count_floor and count_ceil (interpolate linearly)
        cover_extension_diff = cover_extension_ceil - cover_extension_floor
        if cover_extension_diff == 0:
            # count_floor should equal count_ceil
            assert count_floor == count_ceil
            final_interp = count_floor
        elif count_floor <= count_ceil:
            count_diff = count_ceil - count_floor
            f = float(cover_extension - cover_extension_floor) / cover_extension_diff
            final_interp = f * count_diff + count_floor
        else:
            count_diff = count_floor - count_ceil
            f = float(cover_extension - cover_extension_floor) / cover_extension_diff
            final_interp = count_floor - f * count_diff

        return final_interp

    return interp_probe_count_for_dataset


_interp_nd_fn_memoized = {}
def _make_interp_probe_count_for_dataset_nd_fn(probe_counts):
    """Generate and return a function that interpolates probe count for a dataset.

    This uses a function from scipy's interpolate package to operate on
    an arbitrary number of parameters.

    Args:
        probe_counts: dict giving number of probes for each dataset and
            choice of parameters

    Returns:
        function whose input is a dataset and values for arbitrary
        parameters. The function linearly interpolates the number of
        probes required in that dataset for those parameter values, based
        on the values (which were explicitly calculated) in probe_counts.
    """
    # Reset the memoized dict for a new call (this is useful for unit tests,
    # which may call this function multiple times with different inputs --
    # i.e., values in probe_counts)
    _interp_nd_fn_memoized = {}

    def interp_probe_count_for_dataset(dataset, param_vals):
        """
        Using the given probe counts at particular parameter values, interpolate
        the number of probes for 'dataset' and the parameter values given
        in 'param_vals', where each of these may be floats.
        """
        if dataset in _interp_nd_fn_memoized:
            nd_fn = _interp_nd_fn_memoized[dataset]
        else:
            points = []
            values = []
            for p in probe_counts[dataset].keys():
                points += [p]
                values += [probe_counts[dataset][p]]
            points = np.array(points)
            values = np.array(values)

            nd_fn = interpolate.LinearNDInterpolator(points, values,
                rescale=True)
            _interp_nd_fn_memoized[dataset] = nd_fn

        try:
            return nd_fn(np.array(param_vals))[0]
        except ValueError:
            raise ValueError(param_vals, dataset, probe_counts[dataset])

    return interp_probe_count_for_dataset


def _make_total_probe_count_across_datasets_fn(probe_counts,
        interp_fn_type='standard'):
    """Generate and return a function that interpolates probe count.

    Args:
        probe_counts: dict giving number of probes for each dataset and
            choice of parameters
        interp_fn_type: 'standard' (only perform interpolation on mismatches
            and cover_extension parameters) or 'nd' (use scipy's interpolate
            package to interpolate over n-dimensions)

    Returns:
        function whose input is choices of parameter values for all datasets.
        The function interpolates the number of probes required for each
        dataset, and sums across all datasets to determine a total number
        of probes.
    """
    assert interp_fn_type in ['standard', 'nd']
    if interp_fn_type == 'standard':
        interp_fn = _make_interp_probe_count_for_dataset_standard_fn
    elif interp_fn_type == 'nd':
        interp_fn = _make_interp_probe_count_for_dataset_nd_fn
    interp_probe_count_for_dataset = interp_fn(probe_counts)

    def total_probe_count_across_datasets(x):
        """
        Sum the (interpolated) probe counts across datasets.

        Let the number of datasets (len(probe_counts)) be N.
        x is a list giving all the parameter values across datasets,
        such that x_i is the (i % N)'th parameter of the (i/N)'th dataset,
        for i=0,1,2,...
        """
        num_datasets = len(probe_counts)
        # The number of parameter values must be a multiple of the number
        # of datasets
        assert len(x) % num_datasets == 0

        num_params = int(len(x) / num_datasets)

        s = 0
        for i, dataset in enumerate(sorted(probe_counts.keys())):
            param_vals = [x[num_params * i + j] for j in range(num_params)]
            s += interp_probe_count_for_dataset(dataset, param_vals)
        return s

    return total_probe_count_across_datasets
