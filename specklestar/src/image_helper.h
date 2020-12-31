#ifndef _SPECKLE_IMAGE_HELPER_H_
#define _SPECKLE_IMAGE_HELPER_H_

#ifndef IMAGE_SIZE
#define IMAGE_SIZE (512 * 512)

#endif // IMAGE_SIZE

bool IsOverThresholdFrame(unsigned short *piData, unsigned short threshold = 50000);

#endif // _SPECKLE_IMAGE_HELPER_H_
