# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
__author__ = 'walzer'

from operator import attrgetter
from Base import fred2_attrgetter


def uniquify_list(mylist, key=None):
    seen = set()  # similar to F8 http://www.peterbe.com/plog/uniqifiers-benchmark (but mine's nicer: key=attrgetter !)
    if not isinstance(key, attrgetter) and not isinstance(key, fred2_attrgetter):
        return [seen.add(id(obj)) or obj for obj in mylist if id(obj) not in seen]
    return [seen.add(key(obj)) or obj for obj in mylist if key(obj) not in seen]


def lengthrestrict_list(mylist, length):
    assert(isinstance(length, int))
    return [obj for obj in mylist if len(obj) == length]

