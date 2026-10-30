#!/usr/bin/perl
# Harper lint v5 (pmrev6, Spence review 5 item 1: every hedge form rejected through the shared rule in tools/pm_hedge.pl; exit 1 on any BAD). Was v4: UNCLEAR + either/or. Was v3: structure + error variants checked against enums parsed from source (item 10).
use strict; use warnings; use FindBin; require "$FindBin::Bin/pm_hedge.pl";
my ($csv,$T)=@ARGV; $T//="$ENV{HOME}/harper-work/zf-precision-matrix";
sub en { my ($f,$name)=@_; open my $h,'<',"$T/$f" or die "$f"; local $/; my $s=<$h>; my ($b)=$s=~/pub enum \Q$name\E\b[^{]*\{(.*?)\n\}/s or die "enum $name not found in $f"; $b=~s{//[^\n]*}{}g; my %v; my $d=0; for my $l (split /\n/,$b){ $v{$1}=1 if $d==0 && $l=~/^\s*([A-Z]\w*)\b/; $d+=()=$l=~/[{(]/g; $d-=()=$l=~/[})]/g } %v }
my %ZF=en('zenith-float-num/src/defs.rs','Error'); my %ZB=en('zenith-bound/src/error.rs','BoundError');
print "zf Error (zenith-float-num/src/defs.rs): ",join(' ',sort keys %ZF),"\n";
print "zb BoundError (zenith-bound/src/error.rs): ",join(' ',sort keys %ZB),"\n";
sub p{my $l=shift; chomp $l; my @f; my $fld=''; my $st=0; my @c=split //,$l; for(my $i=0;$i<@c;$i++){ my $c=$c[$i];
  if($st==0){ if($c eq '"'){$st=1} elsif($c eq ','){push @f,$fld;$fld=''} else {$fld.=$c} }
  else { if($c eq '"'){ if($i+1<@c && $c[$i+1] eq '"'){$fld.='"';$i++} else {$st=0} } else {$fld.=$c} } } push @f,$fld; @f }
open F,$csv or die; my ($n,$ln,%bad,%cls,%vs,%ex,$ung)=(0,0); my $hdr=0;
my %allv=(%ZF,%ZB);
while(<F>){ $ln++; if(/^#/){ next } if(!$hdr){ $hdr=1; next } $n++; my @f=p($_);
 my $B=sub { my $k=shift; $bad{$k}++; $ex{$k}//="line $ln: $f[1] $f[5]" };
 $B->("fields!=17") if @f!=17; my $t=join " ",@f[10,11];
 $B->("R4 MISSING") if $t=~/R4 MISSING/; $B->("UNCLEAR") if /UNCLEAR/;
  $B->("hedge $_ (tools/pm_hedge.pl rule, cols 10-11)") for pm_hedges($f[10],$f[11]); $B->("'or documented'") if /or documented/;
 my $own = $f[0] eq 'zf' ? \%ZF : \%ZB;
 if($f[0] eq 'zf'){ $B->("zf row names BoundError::") if $t=~/BoundError::/ }
 else { $B->("zb row names zf Error::") if $t=~/(?<!Bound)\bError::/ }
 while($t=~/(?<![A-Za-z])(?:BoundError|Error)::([A-Za-z]\w*)/g){ $B->("unknown variant $1 for $f[0]") unless $own->{$1} }
 # bare CamelCase tokens that are a variant of the OTHER crate only
 while($t=~/\b([A-Z][a-z]+(?:[A-Z][a-z]+)+|Domain|Resource|Singular|Representation)\b/g){ my $w=$1; next unless $allv{$w}; $B->("variant $w not in $f[0] enum") unless $own->{$w} }
 $B->("ruling-needed") if /RULING NEEDED|ruling-needed/i; $B->("waiver/escalat") if /waiv|escalat/i;
 $B->("NA without proof") if $f[15] eq 'NA_PROOF' && $f[16] eq '';
 $B->("empty oracle") if $f[15] ne 'NA_PROOF' && $f[12] eq '';
 $ung++ if $f[15] ne 'NA_PROOF' && "$f[12] $f[13]"=~/UNGENERATED/;
 $cls{$f[5]=~s/ .*//r}++; $vs{"$f[0] $f[15]"}++ }
print "rows $n\n"; print "BAD $_ $bad{$_} (first: $ex{$_})\n" for sort keys %bad; print "lint: ",(%bad?"FAIL":"clean (no BAD lines)"),"\n";
print "rows with an UNGENERATED oracle pin (pm-coverage FAIL, blocks RC until pinned): ",$ung//0,"\n";
print join(" ",map {"$_=$cls{$_}"} sort keys %cls),"\n"; print "$_ $vs{$_}\n" for sort keys %vs;
exit(%bad?1:0);
