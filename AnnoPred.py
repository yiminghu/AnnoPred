#!/usr/bin/env python

from annopred import pred_main, LD_PyWrapper, prior_generating, coord_trimmed

if __name__ == '__main__':
  print(pred_main.main)
  print(LD_PyWrapper.callLDSC)
  print(prior_generating.generate_h2_pT)
  print(prior_generating.generate_h2_from_user)
  print(coord_trimmed.main)
