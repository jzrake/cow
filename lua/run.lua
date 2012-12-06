
local cow = require "cow"
local domain = cow.cow_domain_new()
local dfield = cow.cow_dfield_new()

local modes = bit32.bor(cow.COW_DISABLE_MPI, cow.COW_NOREOPEN_STDOUT)

cow.cow_init(0, nil, modes)
cow.cow_domain_setndim(domain, 3)
cow.cow_domain_setsize(domain, 0, 10)
cow.cow_domain_setsize(domain, 1, 10)
cow.cow_domain_setsize(domain, 2, 10)
cow.cow_domain_commit(domain)

assert(cow.cow_domain_getndim(domain) == 3)

cow.cow_dfield_addmember(dfield, "vx")
cow.cow_dfield_addmember(dfield, "vy")
cow.cow_dfield_addmember(dfield, "vz")
cow.cow_dfield_setname(dfield, "velocity")
cow.cow_dfield_setdomain(dfield, domain)
cow.cow_dfield_commit(dfield)

assert(cow.cow_dfield_getname(dfield) == "velocity")

cow.cow_dfield_write(dfield, "outfile.h5")
cow.cow_dfield_read(dfield, "outfile.h5")
cow.cow_dfield_del(dfield)
cow.cow_domain_del(domain)
cow.cow_finalize()
