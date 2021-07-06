recode_nonwords <- function(x, part_of_speech) {
    x[x == ""] <- NA
    x[x == "xx"] <- NA
    x[part_of_speech == 'bab'] <- '__babble__'
    x[part_of_speech == 'chi'] <- '__childinvent__'
    x[part_of_speech == 'fam'] <- '__familyinvent__'
    x[part_of_speech == 'neo'] <- '__neologism__'
    x[part_of_speech == 'n:prop'] <- '__name__'
    x[part_of_speech == 'sing'] <- '__sing__'
    x[part_of_speech == 'wplay'] <- '__babble__'
    return(x)
}

