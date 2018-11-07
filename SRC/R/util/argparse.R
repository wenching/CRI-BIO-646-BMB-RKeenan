# [Parsing boolean values with argparse](https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse)

str2bool <- function(s) {
  if(any(tolower(s) %in% c('yes', 'true', 't', 'y', '1')))
    return(TRUE)
  else if(any(tolower(s) %in% c('no', 'false', 'f', 'n', '0')))
    return(FALSE)
  else
    stop(paste("ArgumentTypeError('Boolean value expected.'), but", s))
}

# parser.add_argument("--nice", type=str2bool, nargs='?',
#                     const=True, default=NICE,
#                     help="Activate nice mode.")