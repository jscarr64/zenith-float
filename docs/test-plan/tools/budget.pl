sub p{my $l=shift; chomp $l; my @f; my $fld=''; my $st=0; my @c=split //,$l; for(my $i=0;$i<@c;$i++){ my $c=$c[$i];
  if($st==0){ if($c eq '"'){$st=1} elsif($c eq ','){push @f,$fld;$fld=''} else {$fld.=$c} }
  else { if($c eq '"'){ if($i+1<@c && $c[$i+1] eq '"'){$fld.='"';$i++} else {$st=0} } else {$fld.=$c} } } push @f,$fld; @f }
my %thr=(x86_64=>210,i686=>210); # PM-0 corrected survey (1703-pm0-survey.log)
my $W=$ENV{PM_WORKERS}//8; # RC worker count (Spence pmrev3 7)
 my $AX=96;   # survey threshold cells (zf only), heavy axis 96 pts both targets
my %cost=(zf=>4.1,zb=>12); my %fac=(x86_64=>1,i686=>($ENV{PM0_FAC}//2.70)); # PM-0 i686 factor: median 1.91, weighted 2.70 (default), max 5.39
open F,$ARGV[0]; my %c; my $hd=0;
while(<F>){ next if /^#/; if(!$hd){$hd=1;next} my @f=p($_); next if $f[15] eq 'NA_PROOF'; my ($cr,$cls,$modes)=@f[0,5,8]; my $m7 = $modes=~/7 modes/ ? 7:1;
  next if $cls eq 'doc-probable-prime';
  my ($fast,$heavy,$cut);
  if($cls=~/^rm-/){ $fast=1; $heavy=$AX; $cut=6 }
  elsif($cls=~/refuse/){ $fast=1; $heavy=1; $cut=1 }
  elsif($cls=~/^(?:prec-MAX-|prec-2\^(?:30-1|30|60-1)$|m30-)/){ $c{"$cr MAX x86_64"}++; $c{"$cr MAX i686"}++; next } # pmrev8 P1: one row runs on both targets (execution config, not an outcome)
  elsif($cls eq 'prec-1-2-3'){ $fast=3*$m7; $heavy=3*$m7; $cut=3*$m7 }
  elsif($cls=~/^trig-huge/){ $fast=1; $heavy=$AX*$m7; $cut=6*$m7 }
  elsif($cls eq 'prec-axis'){ $fast=6; $heavy=$AX*$m7; $cut=6*$m7 }
  elsif($modes=~/^fast: all 7/){ $fast=7; $heavy=$AX*7; $cut=6*7 }
  elsif($modes=~/^ToEven at p in/){ $fast=6; $heavy=$AX*$m7; $cut=6*$m7 }
  else { $fast=1; $heavy=$AX; $cut=6 }
  $c{"$cr fast"}+=$fast; $c{"$cr heavy"}+=$heavy; $c{"$cr cut"}+=$cut; }
# thresholds (zf): T-1,T,T+1 cells x normal; fast ToEven, cut and heavy x 7 modes
for my $t (keys %thr){ $c{"zf thr $t"}=$thr{$t} }
for my $cr (qw(zf zb)){ for my $t (qw(x86_64 i686)){ my $th= $cr eq 'zf' ? $thr{$t} : 0;
  my $f=($c{"$cr fast"}+$th)*$cost{$cr}*$fac{$t}/1000; my $h=($c{"$cr heavy"}+$th*7)*$cost{$cr}*$fac{$t}/3600000; my $u=($c{"$cr cut"}+$th*7)*$cost{$cr}*$fac{$t}/3600000;
  printf "%s %-6s cells fast=%d cut=%d full=%d | fast %.1f min | cut-heavy %.2f h | full-heavy %.2f h (i686 factor %.2f)\n",$cr,$t,$c{"$cr fast"}+$th,$c{"$cr cut"}+$th*7,$c{"$cr heavy"}+$th*7,$f/60,$u,$h, $fac{$t} } }
my $tot=0; for my $p ([qw(zf x86_64)],[qw(zf i686)],[qw(zb x86_64)],[qw(zb i686)]){ my ($cr,$t)=@$p; my $th=$cr eq "zf"?$thr{$t}:0; $tot+=($c{"$cr heavy"}+$th*7)*$cost{$cr}*$fac{$t}/3600000 }
printf "full heavy total %.2f slot-h -> %.2f h on %d workers\n",$tot,$tot/$W,$W;
my $mx=0; for my $k (sort grep {/MAX/} keys %c){ $mx+=$c{$k}; printf "%s rows %d -> worst case %.1f slot-h at 10 min/row\n",$k,$c{$k},$c{$k}*10/60 }
my %fam; $fam{$_=~s/ MAX .*//r}+=0 for grep {/MAX/} keys %c;
printf "MAX worst case: no skip %.1f slot-h -> %.2f h on %d workers; with MAX-timeout skip (i686 row SKIPPED_AFTER_TIMEOUT when x86_64 times out) %.1f slot-h -> %.2f h\n",$mx*10/60,$mx*10/60/$W,$W,$mx/2*10/60,$mx/2*10/60/$W;
printf "RC wall on %d workers (full heavy + MAX worst case with skip): %.2f h; with no MAX timeouts (MAX rows at 4.1-12 ms): %.2f h\n",$W,($tot+$mx/2*10/60)/$W,$tot/$W;
