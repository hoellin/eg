# Average rank method over BENGI for GM12878

## Prerequisites

Detailed procedure is explained [here](/notes_BENGI/avg_rank_method/avg_rank_method_with_code).

## Not correcting the error in `Run-Average-Rank.sh` gives same results as in the paper

If we keep the small mistake found in `Run-Average-Rank.sh`, we find (hopefully!) the very same Precision - Recall curves as authors:

![Image: Precision-Recall curves and AUPR with small mistake in Run-Average-Rank.sh](same_AUPR_as_authors.png)

## Correct results

If we correct the small error found in `Run-Average-Rank.sh`, we obtain the following results.

![Image: Precision-Recall curves and AUPR](correct_results_with_aupr.png)

![Image: Precision-Recall curves](correct_results.png)

