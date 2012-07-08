

from _cow import *
import os, sys
import atexit

atexit.register(cow_finalize)

cow_init()
domain = cow_domain_new()
dfield = cow_dfield_new(domain, "prim")
cow_dfield_addmember(dfield, "vx")
cow_dfield_addmember(dfield, "vy")
cow_dfield_addmember(dfield, "vz")

cow_domain_setndim(domain, 2);
cow_domain_setguard(domain, 3);
cow_domain_setsize(domain, 0, 10);
cow_domain_setsize(domain, 1, 12);
cow_domain_commit(domain)
cow_dfield_commit(dfield)

print test_trans, TEST_TRANS
print cow_dfield_getname(dfield)
print cow_dfield_getdata(dfield)

A = getarray3(dfield)
print A.shape
print A.flags

#test_trans(None, None, None, None)
#cow_dfield_loop(dfield, TEST_TRANS, None)
cow_dfield_del(dfield)
cow_domain_del(domain)
