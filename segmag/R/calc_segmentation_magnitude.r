calc_segmentation_magnitude <- function(segmag)
{
  # Baut ein Array in das die Segmentierungsstaerke geschrieben wird
  # Um jeden Tastendruck wird Gauss gelegt
  # Gauss vorberechnet mit cutoff
  # Eine VPn jeweils Huellfunktion, damit maximaler Beitrag beschraenkt ist
  # Statt in Zeit wird in Indizes des Arrays gerechnet, da schneller geht.. Dazu alle mittels time_steps umgerechnet, also 1 / time_steps Arrayfelder je Sekunde

  if (! is.segmag(segmag)) stop("segmag must be an object of class segmag")
  
  # Vektor mit Segmentierungsstaerke ueber Zeit in time_steps als Indizes
  segmentation_magnitude_overall <- numeric(segmag$index_time_max+1)
  
  for (id in levels(segmag$ids))
  {    
    index_keypresses <- segmag$index_keypresses[segmag$ids == id]
    
    calc_segmentation_magnitude_impl(segmentation_magnitude_overall,index_keypresses,segmag$gauss_values,segmag$gauss_n_indexes_per_side,segmag$indexes_gauss_offset)
  }
  
  # as.numeric(as.character()) and plyr::round_any: Fix floating point issue causing problems in addressing specific time points (round_any and as.numeric(as.character()) fix different occurances of the issue)
  # Example:
  # tmp <- segmag(factor(c(1)),c(0),gauss_sd = 0.8)
  # tmp$data$segmentation_magnitude[tmp$data$time==0.00]
  # Before Fix returns: numeric(0)
  # After Fix returns: [1] 0.4986779
  return( data.frame(
    time=as.numeric(as.character(
      plyr::round_any(
        seq(
          segmag$time_min,
          segmag$time_min + (segmag$index_time_max*segmag$time_steps),
          segmag$time_steps
        ),
        segmag$time_steps
      )
    )),
    segmentation_magnitude=segmentation_magnitude_overall)
  )
}