##################################################################################################################
# Paper: Property Prediction of Crystalline Solids from Composition and Crystal Structure                        #
# Authors: Bruno A. Calfa (bruno.calfa@wisc.edu), John R. Kitchin (jkitchin@andrew.cmu.edu)                      #
# Supplementary Material: npregkernel.R                                                                          #
##################################################################################################################

npregkernel <- function(txdat,tydat,bw,
                        leave.one.out=FALSE,ukertype="liracine",ckertype="epanechnikov",okertype="liracine",...) {
  npksum(txdat,
         tydat,
         leave.one.out=leave.one.out,
         bws=bw,
         ukertype=ukertype,
         ckertype=ckertype,
         okertype=okertype,...)$ksum/
    npksum(txdat,
           leave.one.out=leave.one.out,
           bws=bw,
           ukertype=ukertype,
           ckertype=ckertype,
           okertype=okertype,...)$ksum
}
