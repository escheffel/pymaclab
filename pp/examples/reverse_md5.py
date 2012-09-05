#!/usr/bin/env python
# File: reverse_md5.py
# Author: Vitalii Vanovschi
# Desc: This program demonstrates parallel computations with pp module
# It tries to reverse an md5 hash in parallel
# Parallel Python Software: http://www.parallelpython.com

import math
import sys
import md5
import pp


def md5test(hash, start, end):
    """Calculates md5 of the integerss between 'start' and 'end' and
       compares it with 'hash'"""
    for x in xrange(start, end):
        if md5.new(str(x)).hexdigest() == hash:
            return x


print """Usage: python reverse_md5.py [ncpus]
    [ncpus] - the number of workers to run in parallel,
    if omitted it will be set to the number of processors in the system
"""

# tuple of all parallel python servers to connect with
#ppservers = ("*",) # auto-discover
#ppservers = ("10.0.0.1","10.0.0.2") # list of static IPs
ppservers = ()

if len(sys.argv) > 1:
    ncpus = int(sys.argv[1])
    # Creates jobserver with ncpus workers
    job_server = pp.Server(ncpus, ppservers=ppservers)
else:
    # Creates jobserver with automatically detected number of workers
    job_server = pp.Server(ppservers=ppservers)

print "Starting pp with", job_server.get_ncpus(), "workers"

#Calculates md5 hash from the given number
hash = md5.new("1829182").hexdigest()
print "hash =", hash
#Now we will try to find the number with this hash value

start = 1
end = 2000000

# Since jobs are not equal in the execution time, division of the problem
# into a 128 of small subproblems leads to a better load balancing
parts = 128

step = (end - start) / parts + 1
jobs = []

for index in xrange(parts):
    starti = start+index*step
    endi = min(start+(index+1)*step, end)
    # Submit a job which will test if a number in the range starti-endi
    # has given md5 hash
    # md5test - the function
    # (hash, starti, endi) - tuple with arguments for md5test
    # () - tuple with functions on which function md5test depends
    # ("md5",) - tuple with module names which must be imported before
    # md5test execution
    # jobs.append(job_server.submit(md5test, (hash, starti, endi),
    # globals=globals()))
    jobs.append(job_server.submit(md5test, (hash, starti, endi),
            (), ("md5", )))

# Retrieve results of all submited jobs
for job in jobs:
    result = job()
    if result:
        break

# Print the results
if result:
    print "Reverse md5 for", hash, "is", result
else:
    print "Reverse md5 for", hash, "has not been found"

job_server.print_stats()

# Properly finalize all tasks (not necessary)
job_server.wait()

# Parallel Python Software: http://www.parallelpython.com
