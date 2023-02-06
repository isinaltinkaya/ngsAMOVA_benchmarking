#total contig size: 100e6 == 1e8
mkdir -pv winranges

#1e5
paste <(seq 0 1e5 999e5) <(seq 1e5 1e5 100e6) > winranges/winranges-contig_100-win_5.bed

#1e6
paste <(seq 0 1e6 99e6) <(seq 1e6 1e6 100e6) > winranges/winranges-contig_100-win_6.bed

#1e7
paste <(seq 0 1e7 9e7) <(seq 1e7 1e7 100e6) > winranges/winranges-contig_100-win_7.bed
