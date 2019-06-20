import numpy as np


WINDOW_OPEN=70  # ns before peak
#window_close=300  #ns after peak
WINDOW_CLOSE=1000 # ns (for data processing)


BG_MIN=0    # ns for avg calculation
BG_MAX=400  # ns for avg calculation




# for rolling average- not to interfere with pulse

N_AVG=40
WINDOW_PRE=50 #ns before window
