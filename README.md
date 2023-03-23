# RIFL

The code are organized into 4 four folders. 

"helper" contains helper functions for the RIFL method as well as other methods for comparison. These codes are imported when running the 3 examples.

"example1_low_dim" contains code for the low-dimensional prediction problem in Appendix C. "mnb_parametric.R" implements the m-out-of-n bootstrap in a separate R script. "RIFL_parametric.R" implements all other methods including RIFL, and also contains the code to generate the figures in Appendix C. 

In "example2_high_dim", "RIFL_high.R" contains code to implement all methods for the high-dimensional prediction example in Section 6.1. The folder “R” contains all necessary code to compute de-biased global dissimilarity measure, which are imported by "RIFL_high.R".

"example3_tate" contains code for the TATE example in Section 6.2. "mnb_tate.R" implements the m-out-of-n bootstrap, while "RIFL_TATE.R" implements all the other methods. "tate.R" implements our doubly robust TATE estimator, which is imported by the other two R scripts.
