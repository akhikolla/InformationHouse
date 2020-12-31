#ifndef __FAMSIZE_H__
#define __FAMSIZE_H__
static void famSize(int famsize[], int fam_group[], int fam_group_size);
#endif


static void famSize(int famsize[], int fam_group[], int fam_group_size)
{
  int i_fgrp = 0;
  int j_fgrp = 0;
  int j_fgrp_temp = -1;
  while(j_fgrp < fam_group_size-1)
  {
    if(fam_group[j_fgrp]!= fam_group[j_fgrp+1])
    {
      famsize[i_fgrp] = j_fgrp - j_fgrp_temp;
      j_fgrp_temp = j_fgrp;
      i_fgrp ++;
    }
    j_fgrp++;
  }
  famsize[i_fgrp] = j_fgrp - j_fgrp_temp;
  return;
}

