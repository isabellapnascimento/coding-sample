#   This script is adapted from an econometrics course where I served as a TA, 
#   rewritten in a pedagogical format for easy replication and learning of basic
#   econometric models. :)

#   It covers three standard tasks:
#     (Section 1) Binary outcome: Linear Probability Model (LPM) vs Probit
#     (Section 2) Difference-in-Differences (DiD)
#     (Section 3) Two-way Fixed Effects (firm & year)
#
# How to use?
#   1) Adjust the working directory below.
#   2) Ensure data files are present in your computer:
#        - smoking.dta  
#        - jtrain.dta  
#      The DiD section uses organ_donations from the 'causaldata' package.
#   3) Run in RStudio. Output is printed to the console.
