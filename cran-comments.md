## New patch submission
This is a submission for a new patch release. In this version I have:

* Addressed valgrind errors on CRAN check results.

* Fixed the memory leaks (debugged with valgrind 3.23.0).

## Previous CRAN check: valgrind

```
==1768996== LEAK SUMMARY:
==1768996==    definitely lost: 240,000 bytes in 3 blocks
==1768996==    indirectly lost: 0 bytes in 0 blocks
==1768996==      possibly lost: 80,000 bytes in 1 blocks
==1768996==    still reachable: 346,246,427 bytes in 74,931 blocks
==1768996==         suppressed: 0 bytes in 0 blocks
==1768996== Reachable blocks (those to which a pointer was found) are not shown.
==1768996== To see them, rerun with: --leak-check=full --show-leak-kinds=all
==1768996== 
==1768996== For lists of detected and suppressed errors, rerun with: -s
==1768996== ERROR SUMMARY: 402 errors from 3 contexts (suppressed: 0 from 0)
```

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
