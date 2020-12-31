#define STRICT_R_HEADERS
#include <Rcpp.h>

using namespace Rcpp;

static const double kSmallEpsilon = 0.00000000000005684341886080801486968994140625;
static const double kExactTestBias = 0.00000000000000000000000010339757656912845935892608650874535669572651386260986328125;
static const double kExactTestEpsilon2 = 0.0000000000009094947017729282379150390625;

/*

Call R function in C++ code

*/

// [[Rcpp::export]]
double SNPHWE2(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, uint32_t midp) {
  // This function implements an exact SNP test of Hardy-Weinberg
  // Equilibrium as described in Wigginton, JE, Cutler, DJ, and
  // Abecasis, GR (2005) A Note on Exact Tests of Hardy-Weinberg
  // Equilibrium. American Journal of Human Genetics. 76: 887 - 893.
  //
  // The original version was written by Jan Wigginton.
  //
  // This version was written by Christopher Chang.  It contains the following
  // improvements over the original SNPHWE():
  // - Proper handling of >64k genotypes.  Previously, there was a potential
  //   integer overflow.
  // - Detection and efficient handling of floating point overflow and
  //   underflow.  E.g. instead of summing a tail all the way down, the loop
  //   stops once the latest increment underflows the partial sum's 53-bit
  //   precision; this results in a large speedup when max heterozygote count
  //   >1k.
  // - No malloc() call.  It's only necessary to keep track of a few partial
  //   sums.
  // - Support for the mid-p variant of this test.  See Graffelman J, Moreno V
  //   (2013) The mid p-value in exact tests for Hardy-Weinberg equilibrium.
  //
  // Note that the SNPHWE_t() function below is a lot more efficient for
  // testing against a p-value inclusion threshold.  SNPHWE2() should only be
  // used if you need the actual p-value.
  intptr_t obs_homc;
  intptr_t obs_homr;
  if (obs_hom1 < obs_hom2) {
    obs_homc = obs_hom2;
    obs_homr = obs_hom1;
  } else {
    obs_homc = obs_hom1;
    obs_homr = obs_hom2;
  }
  const int64_t rare_copies = 2LL * obs_homr + obs_hets;
  const int64_t genotypes2 = (obs_hets + obs_homc + obs_homr) * 2LL;
  if (!genotypes2) {
    if (midp) {
      return 0.5;
    }
    return 1;
  }
  int32_t tie_ct = 1;
  double curr_hets_t2 = obs_hets;
  double curr_homr_t2 = obs_homr;
  double curr_homc_t2 = obs_homc;
  double tailp = (1 - kSmallEpsilon) * kExactTestBias;
  double centerp = 0;
  double lastp2 = tailp;
  double lastp1 = tailp;

  if (obs_hets * genotypes2 > rare_copies * (genotypes2 - rare_copies)) {
    // tail 1 = upper
    while (curr_hets_t2 > 1.5) {
      // het_probs[curr_hets] = 1
      // het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0)
      curr_homr_t2 += 1;
      curr_homc_t2 += 1;
      lastp2 *= (curr_hets_t2 * (curr_hets_t2 - 1)) / (4 * curr_homr_t2 * curr_homc_t2);
      curr_hets_t2 -= 2;
      if (lastp2 < kExactTestBias) {
	tie_ct += (lastp2 > (1 - 2 * kSmallEpsilon) * kExactTestBias);
	tailp += lastp2;
	break;
      }
      centerp += lastp2;
      // doesn't seem to make a difference, but seems best to minimize use of
      // INFINITY
      if (centerp > std::numeric_limits<double>::max()) {
	return 0;
      }
    }
    if ((centerp == 0) && (!midp)) {
      return 1;
    }
    while (curr_hets_t2 > 1.5) {
      curr_homr_t2 += 1;
      curr_homc_t2 += 1;
      lastp2 *= (curr_hets_t2 * (curr_hets_t2 - 1)) / (4 * curr_homr_t2 * curr_homc_t2);
      curr_hets_t2 -= 2;
      const double preaddp = tailp;
      tailp += lastp2;
      if (tailp <= preaddp) {
	break;
      }
    }
    double curr_hets_t1 = obs_hets + 2;
    double curr_homr_t1 = obs_homr;
    double curr_homc_t1 = obs_homc;
    while (curr_homr_t1 > 0.5) {
      // het_probs[curr_hets + 2] = het_probs[curr_hets] * 4 * curr_homr * curr_homc / ((curr_hets + 2) * (curr_hets + 1))
      lastp1 *= (4 * curr_homr_t1 * curr_homc_t1) / (curr_hets_t1 * (curr_hets_t1 - 1));
      const double preaddp = tailp;
      tailp += lastp1;
      if (tailp <= preaddp) {
	break;
      }
      curr_hets_t1 += 2;
      curr_homr_t1 -= 1;
      curr_homc_t1 -= 1;
    }
  } else {
    // tail 1 = lower
    while (curr_homr_t2 > 0.5) {
      curr_hets_t2 += 2;
      lastp2 *= (4 * curr_homr_t2 * curr_homc_t2) / (curr_hets_t2 * (curr_hets_t2 - 1));
      curr_homr_t2 -= 1;
      curr_homc_t2 -= 1;
      if (lastp2 < kExactTestBias) {
	tie_ct += (lastp2 > (1 - 2 * kSmallEpsilon) * kExactTestBias);
	tailp += lastp2;
	break;
      }
      centerp += lastp2;
      if (centerp > std::numeric_limits<double>::max()) {
	return 0;
      }
    }
    if ((centerp == 0) && (!midp)) {
      return 1;
    }
    while (curr_homr_t2 > 0.5) {
      curr_hets_t2 += 2;
      lastp2 *= (4 * curr_homr_t2 * curr_homc_t2) / (curr_hets_t2 * (curr_hets_t2 - 1));
      curr_homr_t2 -= 1;
      curr_homc_t2 -= 1;
      const double preaddp = tailp;
      tailp += lastp2;
      if (tailp <= preaddp) {
	break;
      }
    }
    double curr_hets_t1 = obs_hets;
    double curr_homr_t1 = obs_homr;
    double curr_homc_t1 = obs_homc;
    while (curr_hets_t1 > 1.5) {
      curr_homr_t1 += 1;
      curr_homc_t1 += 1;
      lastp1 *= (curr_hets_t1 * (curr_hets_t1 - 1)) / (4 * curr_homr_t1 * curr_homc_t1);
      const double preaddp = tailp;
      tailp += lastp1;
      if (tailp <= preaddp) {
	break;
      }
      curr_hets_t1 -= 2;
    }
  }
  if (!midp) {
    return tailp / (tailp + centerp);
  }
  return (tailp - ((1 - kSmallEpsilon) * kExactTestBias * 0.5) * tie_ct) / (tailp + centerp);
}

// [[Rcpp::export]]
int32_t SNPHWEX_tailsum(uint32_t high_het_side, double* base_probp, double* saved_hetsp, double* saved_hom1p, double* saved_hom2p, uint32_t* tie_ctp, double *totalp) {
  // similar to fisher23_tailsum()
  double total = 0;
  double cur_prob = *base_probp;
  double tmp_hets = *saved_hetsp;
  double tmp_hom1 = *saved_hom1p;
  double tmp_hom2 = *saved_hom2p;
  double tmps_hets;
  double tmps_hom1;
  double tmps_hom2;
  // identify beginning of tail
  if (high_het_side) {
    if (cur_prob > kExactTestBias) {
      double prev_prob = tmp_hom1 * tmp_hom2;
      while (prev_prob > 0.5) {
	tmp_hets += 2;
	cur_prob *= (4 * prev_prob) / (tmp_hets * (tmp_hets - 1));
	tmp_hom1 -= 1;
	tmp_hom2 -= 1;
	if (cur_prob <= kExactTestBias) {
	  break;
	}
	prev_prob = tmp_hom1 * tmp_hom2;
      }
      *base_probp = cur_prob;
      tmps_hets = tmp_hets;
      tmps_hom1 = tmp_hom1;
      tmps_hom2 = tmp_hom2;
    } else {
      tmps_hets = tmp_hets;
      tmps_hom1 = tmp_hom1;
      tmps_hom2 = tmp_hom2;
      while (1) {
	const double prev_prob = cur_prob;
	tmp_hom1 += 1;
	tmp_hom2 += 1;
	cur_prob *= (tmp_hets * (tmp_hets - 1)) / (4 * tmp_hom1 * tmp_hom2);
	if (cur_prob < prev_prob) {
	  // this should never happen, but better to play it safe re: rounding
	  // error
	  return 1;
	}
	tmp_hets -= 2;
	if (cur_prob > (1 - 2 * kExactTestEpsilon2) * kExactTestBias) {
	  // throw in extra (1 - kSmallEpsilon) multiplier to prevent rounding
	  // errors from causing this to keep going when the left-side test
	  // stopped
	  if (cur_prob > (1 - kSmallEpsilon) * kExactTestBias) {
	    break;
	  }
          *tie_ctp += 1;
	}
	total += cur_prob;
      }
      const double prev_prob = cur_prob;
      cur_prob = *base_probp;
      *base_probp = prev_prob;
    }
  } else {
    if (cur_prob > kExactTestBias) {
      while (tmp_hets > 1.5) {
	tmp_hom1 += 1;
	tmp_hom2 += 1;
	cur_prob *= (tmp_hets * (tmp_hets - 1)) / (4 * tmp_hom1 * tmp_hom2);
	tmp_hets -= 2;
	if (cur_prob <= kExactTestBias) {
	  break;
	}
      }
      *base_probp = cur_prob;
      tmps_hets = tmp_hets;
      tmps_hom1 = tmp_hom1;
      tmps_hom2 = tmp_hom2;
    } else {
      tmps_hets = tmp_hets;
      tmps_hom1 = tmp_hom1;
      tmps_hom2 = tmp_hom2;
      while (1) {
	const double prev_prob = cur_prob;
	tmp_hets += 2;
	cur_prob *= (4 * tmp_hom1 * tmp_hom2) / (tmp_hets * (tmp_hets - 1));
	if (cur_prob < prev_prob) {
	  return 1;
	}
	tmp_hom1 -= 1;
	tmp_hom2 -= 1;
	if (cur_prob > (1 - 2 * kExactTestEpsilon2) * kExactTestBias) {
	  if (cur_prob > kExactTestBias) {
	    break;
	  }
          *tie_ctp += 1;
	}
	total += cur_prob;
      }
      const double prev_prob = cur_prob;
      cur_prob = *base_probp;
      *base_probp = prev_prob;
    }
  }
  *saved_hetsp = tmp_hets;
  *saved_hom1p = tmp_hom1;
  *saved_hom2p = tmp_hom2;
  if (cur_prob > (1 - 2 * kExactTestEpsilon2) * kExactTestBias) {
    if (cur_prob > kExactTestBias) {
      // even most extreme table on this side is too probable
      *totalp = 0;
      return 0;
    }
    *tie_ctp += 1;
  }
  // sum tail to floating point precision limit
  if (high_het_side) {
    while (1) {
      const double prev_tot = total;
      total += cur_prob;
      if (total <= prev_tot) {
	break;
      }
      tmps_hets += 2;
      cur_prob *= (4 * tmps_hom1 * tmps_hom2) / (tmps_hets * (tmps_hets - 1));
      tmps_hom1 -= 1;
      tmps_hom2 -= 1;
    }
  } else {
    while (1) {
      const double prev_tot = total;
      total += cur_prob;
      if (total <= prev_tot) {
	break;
      }
      tmps_hom1 += 1;
      tmps_hom2 += 1;
      cur_prob *= (tmps_hets * (tmps_hets - 1)) / (4 * tmps_hom1 * tmps_hom2);
      tmps_hets -= 2;
    }
  }
  *totalp = total;
  return 0;
}

// [[Rcpp::export]]
double SNPHWEX(int32_t female_hets, int32_t female_hom1, int32_t female_hom2, int32_t male1, int32_t male2, uint32_t midp) {
  // See Graffelman J, Weir BS (2016) Testing for Hardy-Weinberg equilibrium at
  // biallelic genetic markers on the X chromosome.
  // Evaluation strategy is similar to fisher23().
  if ((!male1) && (!male2)) {
    return SNPHWE2(female_hets, female_hom1, female_hom2, midp);
  }
  double cur_prob = (1 - kExactTestEpsilon2) * kExactTestBias;
  double tailp = cur_prob;
  double centerp = 0;
  uint32_t tie_ct = 1;
  // 1. Determine relative tail vs. center masses for the male1/male2-unchanged
  //    slice.
  double cur_female_hetd = (double)female_hets;
  double cur_female_hom1d = (double)female_hom1;
  double cur_female_hom2d = (double)female_hom2;
  double n1 = cur_female_hetd + 2 * cur_female_hom1d;
  double n2 = cur_female_hetd + 2 * cur_female_hom2d;
  double tmp_hets = cur_female_hetd;
  // "left" = low hets side, "right" = high hets side
  double orig_base_probl;
  double orig_base_probr;
  double orig_saved_lhets;
  double orig_saved_lhom1;
  double orig_saved_lhom2;
  double orig_saved_rhets;
  double orig_saved_rhom1;
  double orig_saved_rhom2;
  if (cur_female_hetd * (n1 + n2) > n1 * n2) {
    // current het count is greater than expected 2f(1-f), so we're on the
    // "right" side
    orig_base_probr = cur_prob;
    orig_saved_rhets = cur_female_hetd;
    orig_saved_rhom1 = cur_female_hom1d;
    orig_saved_rhom2 = cur_female_hom2d;

    // scan leftwards
    double tmp_hom1 = cur_female_hom1d;
    double tmp_hom2 = cur_female_hom2d;
    while (tmp_hets > 1.5) {
      tmp_hom1 += 1;
      tmp_hom2 += 1;
      cur_prob *= (tmp_hets * (tmp_hets - 1)) / (4 * tmp_hom1 * tmp_hom2);
      tmp_hets -= 2;
      if (cur_prob < kExactTestBias) {
	tie_ct += (cur_prob > (1 - 2 * kExactTestEpsilon2) * kExactTestBias);
	tailp += cur_prob;
	break;
      }
      centerp += cur_prob;
      if (centerp > std::numeric_limits<double>::max()) {
	return 0;
      }
    }
    orig_saved_lhets = tmp_hets;
    orig_saved_lhom1 = tmp_hom1;
    orig_saved_lhom2 = tmp_hom2;
    orig_base_probl = cur_prob;
    while (tmp_hets > 1.5) {
      tmp_hom1 += 1;
      tmp_hom2 += 1;
      cur_prob *= (tmp_hets * (tmp_hets - 1)) / (4 * tmp_hom1 * tmp_hom2);
      tmp_hets -= 2;
      const double preaddp = tailp;
      tailp += cur_prob;
      if (tailp <= preaddp) {
	break;
      }
    }
    tmp_hets = cur_female_hetd;
    tmp_hom1 = cur_female_hom1d;
    tmp_hom2 = cur_female_hom2d;
    cur_prob = orig_base_probr;
    while (1) {
      tmp_hets += 2;
      cur_prob *= (4 * tmp_hom1 * tmp_hom2) / (tmp_hets * (tmp_hets - 1));
      const double preaddp = tailp;
      tailp += cur_prob;
      if (tailp <= preaddp) {
	break;
      }
      tmp_hom1 -= 1;
      tmp_hom2 -= 1;
    }
  } else {
    // on the "left" side
    orig_base_probl = cur_prob;
    orig_saved_lhets = cur_female_hetd;
    orig_saved_lhom1 = cur_female_hom1d;
    orig_saved_lhom2 = cur_female_hom2d;

    // scan rightwards
    double tmp_hom1 = cur_female_hom1d;
    double tmp_hom2 = cur_female_hom2d;
    double quarter_numer;
    while (1) {
      quarter_numer = tmp_hom1 * tmp_hom2;
      if (quarter_numer <= 0.5) {
	break;
      }
      tmp_hets += 2;
      cur_prob *= (4 * quarter_numer) / (tmp_hets * (tmp_hets - 1));
      tmp_hom1 -= 1;
      tmp_hom2 -= 1;
      if (cur_prob < kExactTestBias) {
	tie_ct += (cur_prob > (1 - 2 * kExactTestEpsilon2) * kExactTestBias);
	tailp += cur_prob;
	quarter_numer = tmp_hom1 * tmp_hom2;
	break;
      }
      centerp += cur_prob;
      if (centerp > std::numeric_limits<double>::max()) {
	return 0;
      }
    }
    orig_saved_rhets = tmp_hets;
    orig_saved_rhom1 = tmp_hom1;
    orig_saved_rhom2 = tmp_hom2;
    orig_base_probr = cur_prob;
    while (quarter_numer > 0.5) {
      tmp_hets += 2;
      cur_prob *= (4 * quarter_numer) / (tmp_hets * (tmp_hets - 1));
      tmp_hom1 -= 1;
      tmp_hom2 -= 1;
      const double preaddp = tailp;
      tailp += cur_prob;
      if (tailp <= preaddp) {
	break;
      }
      quarter_numer = tmp_hom1 * tmp_hom2;
    }
    tmp_hets = cur_female_hetd;
    tmp_hom1 = cur_female_hom1d;
    tmp_hom2 = cur_female_hom2d;
    cur_prob = orig_base_probl;
    while (tmp_hets > 1.5) {
      tmp_hom1 += 1;
      tmp_hom2 += 1;
      cur_prob *= (tmp_hets * (tmp_hets - 1)) / (4 * tmp_hom1 * tmp_hom2);
      const double preaddp = tailp;
      tailp += cur_prob;
      if (tailp <= preaddp) {
	break;
      }
      tmp_hets -= 2;
    }
  }
  // a "row" holds male1/male2 constant.
  const double orig_row_prob = tailp + centerp;
  n1 += male1;
  n2 += male2;
  for (uint32_t male1_decreasing = 0; male1_decreasing < 2; ++male1_decreasing) {
    double cur_male1 = male1;
    double cur_male2 = male2;
    double row_prob = orig_row_prob;
    double cur_lhets = orig_saved_lhets;
    double cur_lhom1 = orig_saved_lhom1;
    double cur_lhom2 = orig_saved_lhom2;
    double cur_rhets = orig_saved_rhets;
    double cur_rhom1 = orig_saved_rhom1;
    double cur_rhom2 = orig_saved_rhom2;
    double base_probl = orig_base_probl;
    double base_probr = orig_base_probr;
    uint32_t iter_ct;
    if (male1_decreasing) {
      iter_ct = 2 * female_hom2 + female_hets;
      if (iter_ct > ((uint32_t)male1)) {
	iter_ct = male1;
      }
    } else {
      iter_ct = 2 * female_hom1 + female_hets;
      if (iter_ct > ((uint32_t)male2)) {
	iter_ct = male2;
      }
    }
    for (uint32_t iter_idx = 0; iter_idx < iter_ct; ++iter_idx) {
      if (male1_decreasing) {
	const double old_male1 = cur_male1;
	const double old_female2 = n2 - cur_male2;
	cur_male2 += 1;
	cur_male1 -= 1;
	// row likelihood is ((n1 choose male1) * (n2 choose male2)) /
	//   ((n1 + n2) choose (male1 + male2))
	row_prob *= (old_male1 * old_female2) / (cur_male2 * (n1 - cur_male1));
	// bugfix (19 Apr 2017): We cannot move to the right of the mode here.
	// Otherwise, if the mode itself is more probable than our initial
	// table, but the table to the immediate right of the mode is not,
	// we'll fail to count the mode.
	// ("right" = high het count, "left" = low het count.)
	if (cur_lhets) {
	  cur_lhom1 += 1;
	  base_probl *= (old_male1 * cur_lhets) / (2 * cur_male2 * cur_lhom1);
	  cur_lhets -= 1;
	} else {
	  cur_lhets += 1;
	  base_probl *= (2 * old_male1 * cur_lhom2) / (cur_male2 * cur_lhets);
	  cur_lhom2 -= 1;
	}
      } else {
	const double old_male2 = cur_male2;
	const double old_female1 = n1 - cur_male1;
	cur_male1 += 1;
	cur_male2 -= 1;
	row_prob *= (old_male2 * old_female1) / (cur_male1 * (n2 - cur_male2));
	if (cur_lhets) {
	  cur_lhom2 += 1;
	  base_probl *= (old_male2 * cur_lhets) / (2 * cur_male1 * cur_lhom2);
	  cur_lhets -= 1;
	} else {
	  cur_lhets += 1;
	  base_probl *= (2 * old_male2 * cur_lhom1) / (cur_male1 * cur_lhets);
	  cur_lhom1 -= 1;
	}
      }
      double tail_incr1;
      if (SNPHWEX_tailsum(0, &base_probl, &cur_lhets, &cur_lhom1, &cur_lhom2, &tie_ct, &tail_incr1)) {
	// all tables in this row, and all subsequent rows, are less probable
	// than the initial table.
	double cur_female1 = n1 - cur_male1;
	double cur_female2 = n2 - cur_male2;
	if (male1_decreasing) {
	  while (1) {
	    const double preaddp = tailp;
	    tailp += row_prob;
	    if (tailp == preaddp) {
	      break;
	    }
	    cur_male2 += 1;
	    cur_female1 += 1;
	    row_prob *= (cur_male1 * cur_female2) / (cur_male2 * cur_female1);
	    cur_male1 -= 1;
	    cur_female2 -= 1;
	  }
	} else {
	  while (1) {
	    const double preaddp = tailp;
	    tailp += row_prob;
	    if (tailp == preaddp) {
	      break;
	    }
	    cur_male1 += 1;
	    cur_female2 += 1;
	    row_prob *= (cur_male2 * cur_female1) / (cur_male1 * cur_female2);
	    cur_male2 -= 1;
	    cur_female1 -= 1;
	  }
	}
	break;
      }
      tailp += tail_incr1;
      if (male1_decreasing) {
	const double old_male1 = cur_male1 + 1;
	if (cur_rhom2) {
	  cur_rhets += 1;
	  base_probr *= (2 * old_male1 * cur_rhom2) / (cur_male2 * cur_rhets);
	  cur_rhom2 -= 1;
	} else {
	  cur_rhom1 += 1;
	  base_probr *= (old_male1 * cur_rhets) / (2 * cur_male2 * cur_rhom1);
	  cur_rhets -= 1;
	}
      } else {
	const double old_male2 = cur_male2 + 1;
	if (cur_rhom1) {
	  cur_rhets += 1;
	  base_probr *= (2 * old_male2 * cur_rhom1) / (cur_male1 * cur_rhets);
	  cur_rhom1 -= 1;
	} else {
	  cur_rhom2 += 1;
	  base_probr *= (old_male2 * cur_rhets) / (2 * cur_male1 * cur_rhom2);
	  cur_rhets -= 1;
	}
      }
      double tail_incr2 = 0.0; // maybe-uninitialized warning
      SNPHWEX_tailsum(1, &base_probr, &cur_rhets, &cur_rhom1, &cur_rhom2, &tie_ct, &tail_incr2);
      tailp += tail_incr2;
      centerp += row_prob - tail_incr1 - tail_incr2;
      if (centerp > std::numeric_limits<double>::max()) {
	return 0;
      }
    }
  }
  if (!midp) {
    return tailp / (tailp + centerp);
  }
  return (tailp - ((1 - kExactTestEpsilon2) * kExactTestBias * 0.5) * ((int32_t)tie_ct)) / (tailp + centerp);
}


