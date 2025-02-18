import numpy as np
from dvsg import calculate_DVSG

sv = np.ones([2, 2])

gv = np.ones([3, 3])

dvsg, dvsg_map = calculate_DVSG(sv, gv)