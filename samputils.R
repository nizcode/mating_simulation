# utilities from sampmate script
# gengt(): function to generate a new gt assignment given two gt's.
gengt <- function(gt1, gt2)
{
    if((gt1 == 0) && (gt2 ==0)) { 
        return(0) # two ref alleles? always homozygous for ref allele
    } else if(((gt1 == 0) && (gt2 ==2)) || ((gt1 == 2) && (gt2 ==0))) {
        return(1) # homozyg for ref and homozyg for alt? Then always heterozyg
    } else if((gt1 == 2) && (gt2 ==2)) {
        return(2) # two alt alleles? always homozygous for alt allele
    } else if(((gt1 == 0) && (gt2 ==1)) || ((gt1 == 1) && (gt2 ==0))) {
        return(round(runif(1))) # homozyg for ref and heterozyg? 50% chance of homozyg ref, 50% chance heterozyg
    } else if(((gt1 == 2) && (gt2 ==1)) || ((gt1 == 1) && (gt2 ==2))) {
        return(round(1+runif(1))) # homozyg for alt and heterozyg? 50% chance of homozyg alt, 50% chance heterozyg
    } else if((gt1 == 1) && (gt2 ==1)) {
        # both heterozyg's? slightly more complicated: 25% chance homozyg ref, 50% heterozyg, 25 homozyg alt
        r <- runif(1)
        if(r < 0.25) {
            return(0)
        } else if(r < 0.75) {
            return(1)
        } else {
            return(2)
        }
    }
}

givesx <- function(numsexes)
{
    return(round(runif(numsexes))) # entirely random, could be male or female either way.
}
