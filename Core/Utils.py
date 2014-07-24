# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
__author__ = 'walzer'

def uniquify_list(mylist, key=None):
    seen = set()  # similar to F8 http://www.peterbe.com/plog/uniqifiers-benchmark (but mine's nicer: key=attrgetter !)
    if not key:
        key = id
    return [seen.add(key(obj)) or obj for obj in mylist if key(obj) not in seen]

