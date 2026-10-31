#!/usr/bin/perl
# Harper budget v4 (pmrev12): measured factors for every active target in inputs/targets.csv (pmrev12 sample log pop-logs/harper/2026-10-30/1931-2036-targets.log);
# prec-PM_MAX_P+1 rows are refused up front by the precision check (zf num.rs:55-57, zb precision.rs:17) before any allocation, so they are costed as one
# ordinary row (4.1 ms zf, 12 ms zb, an upper bound: the check runs before the work the per-row cost measures), not at the 10-minute MAX budget.
# MAX rows (prec-PM_MAX_P, prec-2^60-1, m30-*) are shown both ways: worst case 10 min per row per active target, and no-timeout at the per-row cost.
use strict; use warnings;
sub p{my $l=shift; chomp $l; my @f; my $fld=''; my $st=0; my @c=split //,$l; for(my $i=0;$i<@c;$i++){ my $c=$c[$i];
  if($st==0){ if($c eq '"'){$st=1} elsif($c eq ','){push @f,$fld;$fld=''} else {$fld.=$c} }
  else { if($c eq '"'){ if($i+1<@c && $c[$i+1] eq '"'){$fld.='"';$i++} else {$st=0} } else {$fld.=$c} } } push @f,$fld; @f }
my $THR=210; # threshold cells per target, PM-0 corrected survey (1703-pm0-survey.log)
my $W=$ENV{PM_WORKERS}//8; my $SL=2; my $AX=49; my %cost=(zf=>4.1,zb=>12); my (%fac,@ACT);
{ my $tf=$ENV{PM_TARGETS}//"$ENV{HOME}/harper-work/zf-precision-matrix/docs/test-plan/inputs/targets.csv"; open my $th,'<',$tf or die "$tf"; <$th>;
  while(<$th>){ chomp; next if /^\s*$/; my ($tg,$run,$st,$in,$fa,$fl)=split /,/,$_,-1; $fa//=''; $fl//='';
    if($fa ne ''){ die "budget FAIL: $tg factor $fa has no measurement log ($fl)\n" unless $fl ne '' && -f "$ENV{HOME}/$fl" }
    if($st eq 'active'){ die "budget FAIL: active $tg has no measured factor\n" if $fa eq ''; $fac{$tg}=$fa; push @ACT,$tg; print "factor $tg: $fa (measured, $fl)\n" }
    else { print "factor $tg: none ($st, NOT-RUN, not budgeted)\n" } } }
open F,$ARGV[0] or die; my %c; my $hd=0;
while(<F>){ next if /^#/; if(!$hd){$hd=1;next} my @f=p($_); next if $f[15] eq 'NA_PROOF'; my ($cr,$cls,$modes)=@f[0,5,8]; my $m7 = $modes=~/7 modes/ ? 7:1;
  next if $cls=~/^doc-/;
  my ($fast,$heavy,$cut);
  if($cls=~/^rm-/){ $fast=1; $heavy=$AX; $cut=6 }
  elsif($cls=~/^m32-bytes/ || $cls=~/refuse/ || $cls eq 'prec-PM_MAX_P+1'){ $fast=1; $heavy=1; $cut=1; $c{"$cr refuse-upfront"}++ if $cls eq 'prec-PM_MAX_P+1' }
  elsif($cls=~/^(?:prec-MAX-|prec-PM_MAX_P$|prec-2\^60-1$|m30-)/){ $c{"$cr MAX"}++; next }
  elsif($cls eq 'prec-1-2-3'){ $fast=3*$m7; $heavy=3*$m7; $cut=3*$m7 }
  elsif($cls=~/^trig-huge/){ $fast=1; $heavy=$AX*$m7; $cut=6*$m7 }
  elsif($cls eq 'prec-axis'){ $fast=6; $heavy=$AX*$m7; $cut=6*$m7 }
  elsif($modes=~/^fast: all 7/){ $fast=7; $heavy=$AX*7; $cut=6*7 }
  elsif($modes=~/^ToEven at p in/){ $fast=6; $heavy=$AX*$m7; $cut=6*$m7 }
  else { $fast=1; $heavy=$AX; $cut=6 }
  $c{"$cr fast"}+=$fast; $c{"$cr heavy"}+=$heavy; $c{"$cr cut"}+=$cut; }
printf "prec-PM_MAX_P+1 rows costed as up-front refusals: zf %d, zb %d\n",$c{'zf refuse-upfront'}//0,$c{'zb refuse-upfront'}//0;
my $mxrows=($c{'zf MAX'}//0)+($c{'zb MAX'}//0); printf "MAX rows (prec-PM_MAX_P, prec-2^60-1, m30-*): zf %d, zb %d\n",$c{'zf MAX'}//0,$c{'zb MAX'}//0;
my ($T,$TW,$TN)=(0,0,0);
for my $t (@ACT){ my $ht=0; my $nt=0;
  for my $cr (qw(zf zb)){ my $th=$cr eq 'zf'?$THR:0; my $h=($c{"$cr heavy"}+$th*7)*$cost{$cr}*$fac{$t}/3600000; my $f=($c{"$cr fast"}+$th)*$cost{$cr}*$fac{$t}/60000;
    my $mn=($c{"$cr MAX"}//0)*$cost{$cr}*$fac{$t}/3600000; $ht+=$h; $nt+=$mn;
    printf "  %s %s: fast %.1f min, full heavy %.2f slot-h (cells %d x %.1f ms x %.2f)\n",$t,$cr,$f,$h,$c{"$cr heavy"}+$th*7,$cost{$cr},$fac{$t} }
  my $wc=$mxrows*10/60; printf "%s: full heavy %.2f slot-h; MAX worst %.1f slot-h (%d rows x 10 min); MAX no-timeout %.4f slot-h\n",$t,$ht,$wc,$mxrows,$nt;
  $T+=$ht; $TW+=$wc; $TN+=$nt }
printf "TOTAL over %d active targets: full heavy %.2f slot-h; MAX worst %.1f slot-h; MAX no-timeout %.4f slot-h\n",scalar(@ACT),$T,$TW,$TN;
for my $n ($SL,$W){ printf "RC wall on %d %s: worst case (full heavy + every MAX row times out) %.2f h; no MAX timeouts %.2f h\n",$n,$n==$SL?'popslot slots (Harper)':'workers (overnight)',($T+$TW)/$n,($T+$TN)/$n }
