#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void calc_segmentation_magnitude_impl(NumericVector segmentation_magnitude_overall,NumericVector index_keypresses,NumericVector gauss_values,int gauss_n_indexes_per_side,int indexes_gauss_offset)
{
  int n_segmentation_magnitude_overall = segmentation_magnitude_overall.size();
  NumericVector segmentation_magnitude_particpant(n_segmentation_magnitude_overall);

  /** - Center a gaussian around each keypress of a participant
   *  - If two gaussians overlap then the bigger value is used (not summed up!)
   *    --> restrict the influence of a single participant on the overall segmentation magnitude
   *        (if participant constantly presses a key at one point in time this will not
   *         cause an excessive increase to the overall segmentation magnitude)
   */

  for (int i = 0; i < index_keypresses.size(); i++)
  {
    int index_cur_keypress = index_keypresses[i];
    
    for (int gauss_i = 0; gauss_i < gauss_n_indexes_per_side*2+1; gauss_i++)
    {
      // Maximum of gaussian is located at gauss_n_indexes_per_side
      // --> maximum centered on index_cur_keypress (may vary because of indexes_gauss_offset)
      int index_current = index_cur_keypress - gauss_n_indexes_per_side + gauss_i + indexes_gauss_offset;
      
      if (index_current >= 0 && index_current < n_segmentation_magnitude_overall)
      {
        if (segmentation_magnitude_particpant[index_current] < gauss_values[gauss_i])
        {
          segmentation_magnitude_particpant[index_current] = gauss_values[gauss_i];
        }
      }
    }
  }
  
  /** Add segmentation magnitude of current participant to overall segmentation magnitude
   *  of all participants
   *  += : add values in-place --> R object is directly manipulated without copying to new
   *       memory area (speed benefit)
   */ 
  
  for (int i = 0; i < n_segmentation_magnitude_overall; i++)
  {
    segmentation_magnitude_overall[i] += segmentation_magnitude_particpant[i];
  }
}
