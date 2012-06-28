

import _cow
import os
print os.getcwd()
d = _cow.cow_domain_new()
f = _cow.cow_dfield_new(d, "primitive")

_cow.cow_domain_commit(d)
_cow.cow_dfield_commit(f)

print _cow.cow_dfield_getname(f)

_cow.cow_dfield_del(f)
_cow.cow_domain_del(d)
