POP2=JIGA
for POP in MAQU MBNS
  do
            echo $POP
                    realSFS Results/$POP.saf.idx Results/$POP2.saf.idx > Results/$POP.$POP2.sfs
                    done

# we also need the comparison between MAQU and MBNS 
realSFS Results/MAQU.saf.idx Results/MBNS.saf.idx > Results/MAQU.MBNS.sfs

https://github.com/nt246/physalia-lcwgs/blob/main/day_4/markdowns/02_fst.md
realSFS fst index Results/MAQU.saf.idx Results/MBNS.saf.idx Results/JIGA.saf.idx -sfs Results/MAQU.MBNS.sfs -sfs Results/MAQU.JIGA.sfs -sfs Results/MBNS.JIGA.sfs -fstout Results/JIGA.pbs -whichFst 1
