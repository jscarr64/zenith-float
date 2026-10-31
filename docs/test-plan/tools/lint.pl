#!/usr/bin/perl
# Harper lint v7 (pmrev8, pop-logs/PRINCIPLES.md P1: lint rejects any per-target gold, limit, refusal variant or budget, and FAILs any row whose expected outcome could differ by target, via pm_target_dependent in tools/pm_target.pl). Was v6 (pmrev7, Jeff ruling 2026-10-30 6:00 PM CT, ZF-M27: any target-specific gold, pin or class is BAD through the shared rule in tools/pm_target.pl). Was v5 (pmrev6, Spence review 5 item 1: every hedge form rejected through the shared rule in tools/pm_hedge.pl; exit 1 on any BAD). Was v4: UNCLEAR + either/or. Was v3: structure + error variants checked against enums parsed from source (item 10).
use strict; use warnings; use FindBin; require "$FindBin::Bin/pm_maxp.pl"; use FindBin; require "$FindBin::Bin/pm_hedge.pl"; require "$FindBin::Bin/pm_target.pl";
my $LM="design-stage"; if(@ARGV && $ARGV[0] eq "--mode"){ shift; $LM=shift; die "mode must be rc or design-stage\n" unless $LM=~/^(?:rc|design-stage)$/ }
our @MAXROWS; my ($csv,$T)=@ARGV; $T//="$ENV{HOME}/harper-work/zf-precision-matrix";
sub en { my ($f,$name)=@_; open my $h,'<',"$T/$f" or die "$f"; local $/; my $s=<$h>; my ($b)=$s=~/pub enum \Q$name\E\b[^{]*\{(.*?)\n\}/s or die "enum $name not found in $f"; $b=~s{//[^\n]*}{}g; my %v; my $d=0; for my $l (split /\n/,$b){ $v{$1}=1 if $d==0 && $l=~/^\s*([A-Z]\w*)\b/; $d+=()=$l=~/[{(]/g; $d-=()=$l=~/[})]/g } %v }
my %ZF=en('zenith-float-num/src/defs.rs','Error'); my %ZB=en('zenith-bound/src/error.rs','BoundError');
print "zf Error (zenith-float-num/src/defs.rs): ",join(' ',sort keys %ZF),"\n";
print "zb BoundError (zenith-bound/src/error.rs): ",join(' ',sort keys %ZB),"\n";
sub p{my $l=shift; chomp $l; my @f; my $fld=''; my $st=0; my @c=split //,$l; for(my $i=0;$i<@c;$i++){ my $c=$c[$i];
  if($st==0){ if($c eq '"'){$st=1} elsif($c eq ','){push @f,$fld;$fld=''} else {$fld.=$c} }
  else { if($c eq '"'){ if($i+1<@c && $c[$i+1] eq '"'){$fld.='"';$i++} else {$st=0} } else {$fld.=$c} } } push @f,$fld; @f }
open F,$csv or die; my ($n,$ln,%bad,%cls,%vs,%ex,$ung)=(0,0); my $hdr=0;
my %allv=(%ZF,%ZB);
while(<F>){ $ln++; if(/^#/){ next } if(!$hdr){ $hdr=1; $bad{"P6 32-bit reference in header '$_'"}++ for pm_p6($_); $bad{"per-target column in header '$1' (targets live in inputs/targets.csv)"}++ while /\b(\w*(?:x86_64|aarch64|i686|windows|darwin|wasm32|apple|linux)\w*)\b/gi; next } $n++; my @f=p($_);
 my $B=sub { my $k=shift; $bad{$k}++; $ex{$k}//="line $ln: $f[1] $f[5]" };
 $B->("P6 32-bit reference '$_' (tools/pm_target.pl pm_p6; 64-bit minimum, ZF-M31)") for pm_p6($_);
 $B->("fields!=17") if @f!=17; my $t=join " ",@f[10,11];
 $B->("R4 MISSING") if $t=~/R4 MISSING/; $B->("UNCLEAR") if /UNCLEAR/;
  $B->("hedge $_ (tools/pm_hedge.pl rule, cols 10-11)") for pm_hedges($f[10],$f[11]); $B->("'or documented'") if /or documented/; $B->("target-specific $_ (tools/pm_target.pl rule, ZF-M27: one gold per row for both targets)") for pm_target_specific(@f);
 my $own = $f[0] eq 'zf' ? \%ZF : \%ZB;
 if($f[0] eq 'zf'){ $B->("zf row names BoundError::") if $t=~/BoundError::/ }
 else { $B->("zb row names zf Error::") if $t=~/(?<!Bound)\bError::/ }
 while($t=~/(?<![A-Za-z])(?:BoundError|Error)::([A-Za-z]\w*)/g){ $B->("unknown variant $1 for $f[0]") unless $own->{$1} }
 # bare CamelCase tokens that are a variant of the OTHER crate only
 while($t=~/\b([A-Z][a-z]+(?:[A-Z][a-z]+)+|Domain|Resource|Singular|Representation)\b/g){ my $w=$1; next unless $allv{$w}; $B->("variant $w not in $f[0] enum") unless $own->{$w} }
 $B->("ruling-needed") if /RULING NEEDED|ruling-needed/i; $B->("waiver/escalat") if /waiv|escalat/i;
 $B->("P1 target-dependent $_ (tools/pm_target.pl pm_target_dependent)") for pm_target_dependent(@f); $B->("NA without proof") if $f[15] eq 'NA_PROOF' && $f[16] eq '';
 $B->("empty oracle") if $f[15] ne 'NA_PROOF' && $f[12] eq '';
 $ung++ if $f[15] ne 'NA_PROOF' && "$f[12] $f[13]"=~/UNGENERATED/;
 push @MAXROWS,[$ln,"$f[6] $f[10] $f[11]"] if $f[5]=~/^prec-(?:PM_MAX_P|over-MAX)/; $B->("literal MAX precision in a MAX-class row (Paul M32 ruling 2026-10-30 7:18 PM CT: use the named pin PM_MAX_P)") if $f[5]=~/^prec-(?:PM_MAX_P|over-MAX)/ && "$f[6] $f[10] $f[11]"=~/2\^32(?!\d)/; $B->("old MAX class $f[5] (renamed prec-PM_MAX_P, prec-PM_MAX_P+1)") if $f[5]=~/^prec-2\^32/; $cls{$f[5]=~s/ .*//r}++; $vs{"$f[0] $f[15]"}++ }
{ my $tf=$ENV{PM_TARGETS}//"$T/docs/test-plan/inputs/targets.csv"; open my $th,'<',$tf or die "$tf"; my $h=<$th>; chomp $h; $bad{"targets.csv header '$h'"}++ unless $h eq 'target,runner,status,install,factor,factor_log'; my $tn=0;
  while(<$th>){ chomp; next if /^\s*$/; $tn++; my @c=split /,/,$_,-1; my ($tg,$run,$st,$fac,$fl)=@c[0,1,2,4,5]; $fl//='';
   $bad{"targets.csv $tg: status '$st' not active|pending-install|pending-CI"}++ unless $st=~/^(?:active|pending-install|pending-CI)$/;
   $bad{"targets.csv $tg: active target has no runner"}++ if $st eq 'active' && ($run eq '' || $run eq 'none');
   $bad{"targets.csv $tg: active target has no measured factor"}++ if $st eq 'active' && ($fac//'') eq '';
   $bad{"targets.csv $tg: factor '$fac' not from a measurement (factor_log '$fl' missing)"}++ if ($fac//'') ne '' && ($fl eq '' || ! -f "$ENV{HOME}/$fl");
   $bad{"targets.csv $tg: pending-CI target names a Pop runner '$run'"}++ if $st eq 'pending-CI' && $run ne 'none' }
  print "targets.csv: $tn targets\n" }
{ my $cf=$ENV{PM_CONSTANTS}//"$T/docs/test-plan/inputs/pm-constants.csv"; my $sc; { open my $hh,'<',$csv or die; while(<$hh>){ last unless /^#/; ($sc)=/^# source commit ([0-9a-f]{40})/ if !$sc } } die "no source commit in $csv\n" unless $sc;
  my $X=eval { pm_maxp($T,$sc,$cf) }; if(!$X){ my $m=$@; chomp $m; $bad{"PM_MAX_P does not resolve: $m"}++ } else {
   print "PM_MAX_P = zf MAX_PRECISION = $$X{value} ($$X{form}) from $$X{file}:$$X{line} at $sc sha256:$$X{sha}; expected $$X{expected} ($$X{expected_commit} $$X{expected_at}); M32 in source: ".($$X{m32}?"yes":"no")."\n";
   if(!($$X{m32} && $$X{match})){ if($LM eq "rc"){ $bad{"RC: PM_MAX_P resolves from a source without M32 ($$X{expected_commit} not an ancestor of $sc, or value differs from the expected $$X{expected})"}++ } else { print "OPEN PIN PM_MAX_P: the pinned source predates M32; RC lint FAILs until it includes $$X{expected_commit}\n" } }
   my @lits=pm_lits($$X{value},$$X{expected});
   for my $r (@MAXROWS){ my ($ln,$txt)=@$r; for my $q (@lits){ if(index($txt,$q)>=0){ $bad{"literal MAX value '$q' in a MAX-class row (PM_MAX_P resolves from zf's constant, never a literal)"}++; $ex{"literal MAX value '$q' in a MAX-class row (PM_MAX_P resolves from zf's constant, never a literal)"}//="line $ln"; last } } }
   my $TD=$ENV{PM_TOOLS}//"$T/docs/test-plan/tools"; my @tl=pm_lits($$X{expected}, $$X{m32}?($$X{value}):());
   for my $tf (sort glob("$TD/*.pl")){ open my $th,'<',$tf or die; my $tn=0; while(my $tline=<$th>){ $tn++; for my $q (@tl){ if(index($tline,$q)>=0){ $bad{"literal MAX value '$q' in tool $tf:$tn"}++; } } } } } }
print "lint mode $LM\n"; print "rows $n\n"; print "BAD $_ $bad{$_} (first: ".($ex{$_}//"file-level check, no row").")\n" for sort keys %bad; print "lint: ",(%bad?"FAIL":"clean (no BAD lines)"),"\n";
print "rows with an UNGENERATED oracle pin (pm-coverage FAIL, blocks RC until pinned): ",$ung//0,"\n";
print join(" ",map {"$_=$cls{$_}"} sort keys %cls),"\n"; print "$_ $vs{$_}\n" for sort keys %vs;
exit(%bad?1:0);
