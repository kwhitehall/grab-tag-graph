# Status

### Profiling
* I first profiled the entire code to get a feel for how the profile and code would run on my system. Despite an egregious completion time, it was otherwise fine.
* Pursuant to Kim's instructions, I profiled the find_MCC function. This function is called once in the course of the entire code. The results for this are in the file findMCCProfileResultsv1.
* Looking at the profile results, it became clear that two functions were the majority of the time used by find_MCC, checked_nodes_MCC and check_criteria.
* As would be expected, the sum time for each function was slightly smaller than the sum for the time of the funciton that calls it. However, the time difference was extremely small.
* Once I noticed these two check functions, I wrote a new function into find_MCC called bothCriteriaCheck. This function serves as a wrapper for the checking section of the code in find_MCC, to aid as a profiler tool to confirm my suspicions about checked_nodes_MCC and check_criteria.
* The profiler was run again on find_MCC with my new function incorporated. As expected bothCriteriaCheck had a similar sum time to find_MCC, checked_nodes_MCC and check_criteria.
* Of note is that the way the profiler runs in the code right now does not allow for the code to run to completion. The profile stops after running on find_MCC for some reason. This needs to be fixed, obviously.
* Also of note is that the implementation of bothCriteriaCheck needs to be looked at and potentially modified. After the implementation of that function, check_criteria was called seven fewer times than in the original find_MCC. Depending upon the cause of this, it may save time or be hurting the validity of the data. In the latter's case, we are not to worry as this function serves only to point me in the direction of what is using the most time for find_MCC. In the former's case, this may indicate a useful way to cut down on runtime.

### check_criteria
* This function appears to be where the majority of time for find_MCC is spent.
* This function will be picked apart and thoroughly examined to determine what can be optimized here.
* Some of these for loops look like they might hide inefficiencies that can be optimized and cut time.