"""AF3 execution coverage lives at the local runner boundary.

The former tests for parsing ``sacct`` output and polling array elements were
removed because the controller no longer schedules or monitors a child array.
"""
