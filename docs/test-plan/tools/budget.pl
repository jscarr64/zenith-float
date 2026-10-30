sub p{my $l=shift; chomp $l; my @f; my $fld=''; my $st=0; my @c=split //,$l; for(my $i=0;$i<@c;$i++){ my $c=$c[$i];
  if($st==0){ if($c eq '"'){$st=1} elsif($c eq ','){push @f,$fld;$fld=''} else {$fld.=$c} }
  else { if($c eq '"'){ if($i+1<@c && $c[$i+1] eq '"'){$fld.='"';$i++} else {$st=0} } else {$fld.=$c} } } push @f,$fld; @f }
my %thr=(x86_64=>140,i686=>123); my $AX=96;   # survey threshold cells (zf only), heavy axis 96 pts both targets
my %cost=(zf=>4.1,zb=>12); my %fac=(x86_64=>1,i686=>4);
open F,$ARGV[0]; my %c; my $hd=0;
while(<F>){ next if /^#/; if(!$hd){$hd=1;next} my @f=p($_); next if $f[15] eq 'NA_PROOF'; my ($cr,$cls,$modes)=@f[0,5,8]; my $m7 = $modes=~/7 modes/ ? 7:1;
  next if $cls eq 'doc-probable-prime';
  my ($fast,$heavy,$cut);
  if($cls=~/^rm-/){ $fast=1; $heavy=$AX; $cut=6 }
  elsif($cls=~/refuse/){ $fast=1; $heavy=1; $cut=1 }
  elsif($cls=~/^prec-MAX-(\S+)/){ $c{"$cr MAX $1"}++; next }
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
  printf "%s %-6s cells fast=%d cut=%d full=%d | fast %.1f min | cut-heavy %.2f h | full-heavy %.2f h (i686 at 1.9x: fast %.1f min cut %.2f h full %.2f h)\n",$cr,$t,$c{"$cr fast"}+$th,$c{"$cr cut"}+$th*7,$c{"$cr heavy"}+$th*7,$f/60,$u,$h, ($t eq 'i686'?($f/60*1.9/4,$u*1.9/4,$h*1.9/4):(0,0,0)) } }
printf "full heavy with 4 popslot workers (Paul): %s\n", join("; ",map { my ($cr,$t)=@$_; my $th=$cr eq "zf"?$thr{$t}:0; sprintf("%s %s %.2f h (1.9x: %.2f h)",$cr,$t,($c{"$cr heavy"}+$th*7)*$cost{$cr}*$fac{$t}/3600000/4,$t eq "i686"?($c{"$cr heavy"}+$th*7)*$cost{$cr}*$fac{$t}/3600000/4*1.9/4:0) } ([qw(zf x86_64)],[qw(zf i686)],[qw(zb x86_64)],[qw(zb i686)]));
for my $k (sort grep {/MAX/} keys %c){ printf "%s rows %d -> worst case %.1f h at 10 min/row (single worker), %.1f h on 4 workers\n",$k,$c{$k},$c{$k}*10/60,$c{$k}*10/60/4 }
