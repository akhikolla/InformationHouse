#include "image_helper.h"

bool IsOverThresholdFrame(unsigned short *piData, unsigned short threshold)
{
  for(int i = 0; i < IMAGE_SIZE; i++)
    if (piData[i] > threshold)
      return true;

  return false;
}
