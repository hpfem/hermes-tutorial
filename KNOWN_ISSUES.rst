KNOWN ISSUES
============

P04-adaptivity/03-system
~~~~~~~~~~~~~~~~~~~~~~~~ 

Adaptivity works but probably the exact solution of the first component is wrong. 
Everything works, but the exact error of the first component stays over 90%.
Noteworthy is the fact, that it did not work (i.e. was exhibiting the same behavior) in 
these commits:

edc33284dc01fe98e9c9daec105a5a3139c9ab28  2011-03-26 19:58:49
1bc3eff209a6cf71a3f1f98053d612e0970e8b0d  2011-05-07 20:11:47
cf24c53579aa9bda1f6181602b9b52a2a371582e  2011-06-04 12:52:28
a22cc8a8efb5750c093423c381a4cdb11fac9683  (the current HEAD)

in hermes-legacy.