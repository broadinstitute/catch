"""Tests for interval module.
"""

import unittest

from catch.utils import interval
from catch.utils.interval import IntervalSet

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class TestIntervalSet(unittest.TestCase):
    """Tests the IntervalSet class and its methods.
    """

    def test_eq(self):
        # assertEquals should call the __eq__ method
        self.assertEqual(IntervalSet([(1, 3), (3, 5)]), IntervalSet([(1, 5)]))
        self.assertNotEqual(IntervalSet([(1, 3)]), IntervalSet([(5, 7)]))

    def test_num_elements(self):
        self.assertEqual(len(IntervalSet([])), 0)
        self.assertEqual(len(IntervalSet([(1, 5)])), 4)
        self.assertEqual(len(IntervalSet([(1, 3), (3, 5)])), 4)
        self.assertEqual(len(IntervalSet([(1, 3), (5, 9)])), 6)
        self.assertEqual(len(IntervalSet([(0, 10), (3, 6)])), 10)
        self.assertEqual(len(IntervalSet([(1, 5), (3, 7), (10, 13),
                                          (20, 30)])), 19)

    def compare_intersect(self, input_a, input_b, desired_output):
        a = IntervalSet(input_a)
        b = IntervalSet(input_b)
        o = IntervalSet(desired_output)
        # Intersect is commutative, so check both orders
        self.assertEqual(a.intersection(b), o)
        self.assertEqual(b.intersection(a), o)

    def test_intersect(self):
        self.compare_intersect([], [], [])
        self.compare_intersect([], [(1, 3)], [])
        self.compare_intersect([(1, 3)], [], [])
        self.compare_intersect([(1, 3), (10, 15)], [(5, 7), (17, 20)], [])
        self.compare_intersect([(1, 3), (5, 7)], [(1, 3), (5, 7)], [(1, 3),
                                                                    (5, 7)])
        self.compare_intersect([(1, 10)], [(3, 7)], [(3, 7)])
        self.compare_intersect([(2, 100)], [(0, 50)], [(2, 50)])
        self.compare_intersect([(2, 100), (101, 150)], [(0, 50)], [(2, 50)])
        self.compare_intersect([(0, 7)], [(4, 10)], [(4, 7)])
        self.compare_intersect([(1, 5), (10, 15)], [(1, 5),
                                                    (15, 20)], [(1, 5)])
        self.compare_intersect([(1, 5), (10, 15)], [(3, 12)], [(3, 5),
                                                               (10, 12)])
        self.compare_intersect([(1, 5), (10, 15)], [(12, 20), (25, 30)],
                               [(12, 15)])
        self.compare_intersect([(100, 200), (250, 300), (350, 450),
                                (475, 600)], [(175, 275), (250, 350),
                                              (400, 500)],
                               [(175, 200), (250, 300), (400, 450),
                                (475, 500)])

    def compare_union(self, input_a, input_b, desired_output):
        a = IntervalSet(input_a)
        b = IntervalSet(input_b)
        o = IntervalSet(desired_output)

        # Union is commutative, so check both orders
        self.assertEqual(a.union(b), o)
        self.assertEqual(b.union(a), o)

        # Making a new IntervalSet should merge overlapping
        # intervals, effectively taking the union
        ab = IntervalSet(input_a + input_b)
        self.assertEqual(ab, o)

    def test_union(self):
        self.compare_union([], [], [])
        self.compare_union([], [(1, 3)], [(1, 3)])
        self.compare_union([(1, 3)], [], [(1, 3)])
        self.compare_union([(1, 3), (10, 15)], [(5, 7), (17, 20)],
                           [(1, 3), (5, 7), (10, 15), (17, 20)])
        self.compare_union([(1, 3), (5, 7)], [(1, 3), (5, 7)], [(1, 3),
                                                                (5, 7)])
        self.compare_union([(1, 10)], [(3, 7)], [(1, 10)])
        self.compare_union([(2, 100)], [(0, 50)], [(0, 100)])
        self.compare_union([(0, 7)], [(4, 10)], [(0, 10)])
        self.compare_union([(4, 10)], [(0, 7)], [(0, 10)])
        self.compare_union([(1, 5), (10, 15)], [(1, 5), (15, 20)], [(1, 5),
                                                                    (10, 20)])
        self.compare_union([(1, 5), (10, 15)], [(3, 12)], [(1, 15)])
        self.compare_union([(1, 5), (10, 15)], [(12, 20), (25, 30)],
                           [(1, 5), (10, 20), (25, 30)])
        self.compare_union([(100, 200), (250, 300), (350, 450), (475, 600)],
                           [(175, 275), (250, 350), (400, 500)], [(100, 600)])

    def compare_difference(self, input_a, input_b, desired_output):
        a = IntervalSet(input_a)
        b = IntervalSet(input_b)
        o = IntervalSet(desired_output)
        self.assertEqual(a.difference(b), o)

    def test_difference(self):
        self.compare_difference([], [], [])
        self.compare_difference([], [(1, 3)], [])
        self.compare_difference([(1, 3)], [], [(1, 3)])
        self.compare_difference([(1, 3)], [(5, 10)], [(1, 3)])
        self.compare_difference([(5, 10)], [(1, 3)], [(5, 10)])
        self.compare_difference([(1, 10)], [(3, 7)], [(1, 3), (7, 10)])
        self.compare_difference([(1, 10)], [(3, 7), (10, 12)], [(1, 3),
                                                                (7, 10)])
        self.compare_difference([(1, 5), (10, 15)], [(3, 7)], [(1, 3),
                                                               (10, 15)])
        self.compare_difference([(3, 7)], [(1, 5), (10, 15)], [(5, 7)])
        self.compare_difference([(1, 5), (10, 15)], [(3, 12)], [(1, 3),
                                                                (12, 15)])
        self.compare_difference([(3, 12)], [(1, 5), (10, 15)], [(5, 10)])
        self.compare_difference([(1, 5), (10, 15), (20, 30)], [(3, 12)],
                                [(1, 3), (12, 15), (20, 30)])
        self.compare_difference([(100, 200), (250, 300), (350, 450),
                                 (475, 600)], [(175, 275), (250, 350),
                                               (400, 500)],
                                [(100, 175), (350, 400), (500, 600)])

    def compare_overlaps_interval(self, interval_set_intervals,
                                  interval_to_check, desired_output):
        interval_set = IntervalSet(interval_set_intervals)
        start, end = interval_to_check
        self.assertEqual(interval_set.overlaps_interval(start, end),
                         desired_output)

    def test_overlaps_interval(self):
        self.compare_overlaps_interval([(1, 5), (10, 14)], (-2, 6), True)
        self.compare_overlaps_interval([(1, 5), (10, 14)], (3, 5), True)
        self.compare_overlaps_interval([(1, 5), (10, 14)], (4, 8), True)
        self.compare_overlaps_interval([(1, 5), (10, 14)], (-2, 8), True)
        self.compare_overlaps_interval([(1, 5), (10, 14)], (9, 12), True)
        self.compare_overlaps_interval([(1, 5), (10, 14)], (10, 12), True)
        self.compare_overlaps_interval([(1, 5), (10, 14)], (8, 16), True)
        self.compare_overlaps_interval([(1, 5), (10, 14)], (12, 16), True)
        self.compare_overlaps_interval([(1, 5), (10, 14)], (3, 12), True)
        self.compare_overlaps_interval([(1, 5), (10, 14)], (3, 16), True)
        self.compare_overlaps_interval([(1, 5), (10, 14)], (-2, 16), True)
        self.compare_overlaps_interval([(1, 5), (10, 14)], (-2, 1), False)
        self.compare_overlaps_interval([(1, 5), (10, 14)], (5, 8), False)
        self.compare_overlaps_interval([(1, 5), (10, 14)], (6, 9), False)
        self.compare_overlaps_interval([(1, 5), (10, 14)], (14, 20), False)
        self.compare_overlaps_interval([(1, 5), (10, 14)], (16, 20), False)
        self.compare_overlaps_interval([(1, 5), (10, 14)], (5, 10), False)


class TestMergeOverlapping(unittest.TestCase):
    """Tests the merge_overlapping function.
    """

    def compare(self, input, desired_output):
        self.assertEqual(interval.merge_overlapping(input), desired_output)

    def test_single(self):
        self.compare([(1, 2)], [(1, 2)])

    def test_two_nonoverlapping(self):
        self.compare([(1, 2), (4, 5)], [(1, 2), (4, 5)])

    def test_two_overlapping(self):
        self.compare([(1, 3), (2, 5)], [(1, 5)])

    def test_two_touching(self):
        self.compare([(1, 3), (3, 5)], [(1, 5)])

    def test_out_of_order(self):
        self.compare([(4, 5), (1, 2)], [(1, 2), (4, 5)])
        self.compare([(3, 5), (1, 3)], [(1, 5)])

    def test_three(self):
        self.compare([(1, 5), (3, 7), (9, 12)], [(1, 7), (9, 12)])

    def test_within(self):
        self.compare([(1, 5), (3, 4), (6, 8)], [(1, 5), (6, 8)])

    def test_loner_in_middle(self):
        self.compare([(1, 5), (3, 5), (5, 8), (10, 12), (15, 18), (17, 20)],
                     [(1, 8), (10, 12), (15, 20)])


class TestSchedule(unittest.TestCase):
    """Tests the schedule function.
    """

    def compare(self, input, desired_output):
        self.assertEqual(interval.schedule(input), desired_output)

    def test_single(self):
        self.compare([((1, 2), 1)], [1])

    def test_two_nonoverlapping(self):
        self.compare([((1, 2), 1), ((4, 5), 2)], [1, 2])
        self.compare([((4, 5), 2), ((1, 2), 1)], [1, 2])

    def test_two_overlapping(self):
        self.compare([((1, 3), 1), ((2, 5), 2)], [1])
        self.compare([((2, 5), 2), ((1, 3), 1)], [1])

    def test_three(self):
        self.compare([((1, 5), 1), ((3, 7), 2), ((9, 12), 3)], [1, 3])

    def test_within(self):
        self.compare([((1, 5), 1), ((3, 4), 2), ((6, 8), 3)], [2, 3])

    def test_loner_in_middle(self):
        self.compare([((1, 5), 1), ((3, 4), 2), ((5, 8), 3), ((10, 12), 4),
                      ((15, 18), 5), ((17, 20), 6)], [2, 3, 4, 5])
