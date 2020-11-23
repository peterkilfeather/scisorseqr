## Original by HT

#if line starts with ">" 
#if not on the first row of the file, 
#use the newline separator to delimit
#print the name of the fasta sequence without ">"
#if on line 1, do the same
#otherwise if does not start with ">" print whole line (FASTA sequence)
function FastaToTbl()
{
        awk '{
                if (substr($1,1,1)==">")
                        if (NR>1)
                                printf "\n%s ", substr($1,2,length($1)-1)
                        else
                                printf "%s ", substr($1,2,length($1)-1)
                else
                        printf "%s", $0
        }END{printf "\n"}'  "$@"
}


