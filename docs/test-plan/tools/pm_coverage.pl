#!/usr/bin/perl
# Harper pm-coverage v2 (pmrev5: rule f parses full records, --manifest). Was v1 (Spence pmrev3 item 1). Separate tool; lint.pl stays the structural/variant lint.
# usage: perl pm_coverage.pl [--mode rc|design-stage] [--results verdicts.tsv --manifest run.manifest] <edges.csv> <shards.csv> [tree]
#  rc (default, strict): every rule below is enforced; any count > 0 is FAIL and exit 1.
#  design-stage: only the structural rules (S1-S3, e) are enforced; rules a-d are listed as GAP lines with
#    their full counts (never hidden, never PASS); output is marked NOT AN RC VERDICT. Exit 1 on a structural FAIL.
# rules: a unsigned proof | b cited proof file missing | c UNGENERATED oracle pin on an executed row
#        d unpinned PM-1 row ([PM1-UNPINNED]) | e shard per-family row count below matrix count
#        f results records (design s17.1): row_id TAB verdict TAB zf_sha TAB zb_sha TAB env_name TAB env_hash TAB gen_sha256 [TAB target TAB secs].
#          FAIL per: duplicate row_id (never overwritten), row_id not in matrix, malformed record, any of the 5 provenance fields
#          != the --manifest value (key=value lines: zf_sha zb_sha env_name env_hash gen_sha256), executed row missing, verdict != PASS
#          (SKIPPED_AFTER_TIMEOUT / FAIL count as FAIL).
#        S1 every prec-MAX-<target> (timeout-budget) row's verdict space includes SKIPPED_AFTER_TIMEOUT | S2 every NA_PROOF row cites a
#        pm-proofs path | S3 every executed row's verdict space is PASS|FAIL[|SKIPPED_AFTER_TIMEOUT]
use strict; use warnings;
my $mode='rc'; my $res; my $man; my @a;
while(@ARGV){ my $x=shift; if($x eq '--mode'){ $mode=shift } elsif($x eq '--results'){ $res=shift } elsif($x eq '--manifest'){ $man=shift } else { push @a,$x } }
die "mode must be rc or design-stage\n" unless $mode=~/^(rc|design-stage)$/;
my ($csv,$shards,$T)=@a; die "usage\n" unless $shards; $T//="$ENV{HOME}/harper-work/zf-precision-matrix";
sub p{my $l=shift; chomp $l; my @f; my $fld=''; my $st=0; my @c=split //,$l; for(my $i=0;$i<@c;$i++){ my $c=$c[$i];
  if($st==0){ if($c eq '"'){$st=1} elsif($c eq ','){push @f,$fld;$fld=''} else {$fld.=$c} }
  else { if($c eq '"'){ if($i+1<@c && $c[$i+1] eq '"'){$fld.='"';$i++} else {$st=0} } else {$fld.=$c} } } push @f,$fld; @f }
my (%cnt,%ex,%files,%exp,$rows,$na,$exe); my $hdr=0; my $ln=0; my @ids;
my $G=sub { my ($r,$m)=@_; $cnt{$r}++; push @{$ex{$r}},$m if @{$ex{$r}//[]}<3 };
open F,$csv or die "$csv"; while(<F>){ $ln++; next if /^#/; if(!$hdr){$hdr=1;next} $rows++; my @f=p($_); my $id="$f[0]|$f[1]|$f[5]";
 my $shard = $f[15] eq 'NA_PROOF' ? 'proof' : $f[9]=~/pm-max-i686/ ? 'pm-max-i686' : $f[9]=~/pm-max-x86_64/ ? 'pm-max-x86_64' : $f[9]=~/fast/ ? 'fast' : 'heavy';
 $exp{"$shard|$f[0]|$f[1]"}++;
 if($f[15] eq 'NA_PROOF'){ $na++; my ($pf)=$f[16]=~/^(docs\/test-plan\/pm-proofs\/\S+?\.md)#/;
   if(!$pf){ $G->('S2',"line $ln $f[1]"); next }
   $files{$pf}{n}++; if(! -f "$T/$pf"){ $G->('b',"line $ln $f[1] -> $pf missing") }
   else { open my $h,'<',"$T/$pf"; local $/; my $s=<$h>; $G->('a',"line $ln $f[1] -> $pf UNSIGNED") unless $s=~/^Second-agent sign-off: (?!UNSIGNED)\S+ \S+ [0-9a-f]{7,}/m }
   next }
 $exe++; push @ids,$id;
 $G->('S3',"line $ln $f[1] vs=$f[15]") unless $f[15]=~/^PASS\|FAIL(\|SKIPPED_AFTER_TIMEOUT)?$/;
 $G->('S1',"line $ln $f[1] $f[5]") if $f[5]=~/^prec-MAX-/ && $f[15] !~ /SKIPPED_AFTER_TIMEOUT/;
 $G->('c',"line $ln $f[1] $f[5]") if "$f[12] $f[13]"=~/UNGENERATED/;
 $G->('d',"line $ln $f[1] $f[5]") if $f[10]=~/\[PM1-UNPINNED\]/; }
close F;
my %have; open S,$shards or die "$shards"; <S>; while(<S>){ my @f=p($_); $have{"$f[0]|$f[1]|$f[2]"}=$f[3] } close S;
for my $k (sort keys %exp){ my $h=$have{$k}//0; $G->('e',"shard $k has $h rows, matrix has $exp{$k}") if $h<$exp{$k} }
my %fe = map { $_=>(-f "$T/$_"?1:0) } keys %files; my $fx=grep {$fe{$_}} keys %fe; my $fm=(keys %fe)-$fx;
my @MK=qw(zf_sha zb_sha env_name env_hash gen_sha256);
if(defined $res){ die "rule f: --results requires --manifest (expected provenance comes from the manifest, never from the results file)\n" unless defined $man;
   my %M; open my $mh,'<',$man or die "$man"; while(<$mh>){ chomp; next if /^\s*(#|$)/; my ($k,$v)=/^\s*(\w+)\s*=\s*(\S+)\s*$/ or die "manifest: bad line $.\n"; $M{$k}=$v } close $mh;
   for (@MK){ die "manifest: missing $_\n" unless defined $M{$_} }
   my %known=map {$_=>1} @ids; my (%v,%seen); open R,$res or die "$res"; my $rl=0;
   while(<R>){ $rl++; chomp; next if /^#/ || $_ eq ''; my @r=split /\t/,$_,-1;
     if(@r<7){ $G->('f',"results line $rl malformed (".scalar(@r)." fields, need 7)"); next }
     my $i=$r[0];
     if($seen{$i}++){ $G->('f',"results line $rl duplicate row_id $i (first record kept, verdict $v{$i}; this one $r[1])"); next }
     $v{$i}=$r[1];
     $G->('f',"results line $rl row_id $i not in matrix") unless $known{$i};
     for my $k (0..4){ my $got=$r[2+$k]; $G->('f',"results line $rl $i $MK[$k]=$got != manifest $M{$MK[$k]}") unless $got eq $M{$MK[$k]} } }
   close R;
   for my $i (@ids){ my $v=$v{$i}//'MISSING'; $G->('f',"$i verdict $v (counts as FAIL)") unless $v eq 'PASS' } }
print "pm-coverage mode=$mode rows=$rows NA_PROOF=$na executed=$exe; cited proof files: ".(keys %fe)." distinct, $fx exist, $fm missing\n";
my %name=(a=>'unsigned proof',b=>'cited proof file missing',c=>'UNGENERATED oracle pin',d=>'unpinned PM-1 row',e=>'shard per-family count below matrix',f=>'verdict not PASS (SKIPPED_AFTER_TIMEOUT/FAIL/missing = FAIL)',S1=>'MAX row lacks SKIPPED_AFTER_TIMEOUT',S2=>'NA_PROOF row without pm-proofs path',S3=>'bad executed verdict space');
my $fail=0; my %struct=map {$_=>1} qw(S1 S2 S3 e);
for my $r (qw(a b c d e f S1 S2 S3)){ my $n=$cnt{$r}//0;
  if($r eq 'f' && !defined $res){ print "rule f ($name{f}): not evaluated, no --results file (RC requires one)\n"; $fail++ if $mode eq 'rc'; next }
  my $enf = $mode eq 'rc' || $struct{$r}; my $tag = $n==0 ? 'PASS' : $enf ? 'FAIL' : 'GAP';
  $fail++ if $n && $enf; printf "rule %s (%s): %s count=%d%s\n",$r,$name{$r},$tag,$n,($n?" first: ".join(' ; ',@{$ex{$r}}):'') }
print $mode eq 'design-stage' ? "DESIGN-STAGE (structural rules only; GAP lines are open gaps, NOT AN RC VERDICT): ".($fail?"FAIL":"structural PASS")."\n" : "pm-coverage RC: ".($fail?"FAIL":"PASS")."\n";
exit($fail?1:0);
